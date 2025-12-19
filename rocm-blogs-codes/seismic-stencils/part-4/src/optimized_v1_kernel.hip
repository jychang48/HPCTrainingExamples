#include "helper.hpp"
#include "optimized_v1_kernel.hpp"

// Constant memory for finite difference coefficients (local to this TU)
template <int R> __constant__ float opt_v1_dx[2 * R + 1];
template <int R> __constant__ float opt_v1_dy[2 * R + 1];
template <int R> __constant__ float opt_v1_dz[2 * R + 1];

// Initialize constant memory arrays for optimized_v1 kernel
template <int R>
void init_fd_optimized_v1(const float *d) {
    HIP_CHECK(hipMemcpyToSymbol(opt_v1_dx<R>, d, (2 * R + 1) * sizeof(float)));
    HIP_CHECK(hipMemcpyToSymbol(opt_v1_dy<R>, d, (2 * R + 1) * sizeof(float)));
    HIP_CHECK(hipMemcpyToSymbol(opt_v1_dz<R>, d, (2 * R + 1) * sizeof(float)));
}

// Explicit instantiation
template void init_fd_optimized_v1<RADIUS>(const float*);

// Thread block dimensions
#if RADIUS == 1
#define BLOCK_DIM_X 256
#define BLOCK_DIM_Y 1
#elif RADIUS == 2
#define BLOCK_DIM_X 128
#define BLOCK_DIM_Y 2
#elif RADIUS == 3
#define BLOCK_DIM_X 64
#define BLOCK_DIM_Y 3
#elif RADIUS == 4
#define BLOCK_DIM_X 64
#define BLOCK_DIM_Y 4
#elif RADIUS == 6
#define BLOCK_DIM_X 32
#define BLOCK_DIM_Y 4
#elif RADIUS == 8
#define BLOCK_DIM_X 32
#define BLOCK_DIM_Y 4
#else
#error "Unsupported RADIUS value"
#endif

/**
 * Optimized v1 kernel with ILP + explicit FMA
 *
 * Optimizations:
 * 1. Compute z-stencil BEFORE barrier (ILP: overlaps with barrier wait)
 * 2. Use __builtin_fmaf for explicit FMA operations
 * 3. Preload coefficients into registers
 */
template <int R>
__launch_bounds__((BLOCK_DIM_X) * (BLOCK_DIM_Y))
__global__ void optimized_v1_kernel_impl(float *p_out, const float *p_in, int line, int
        slice, int x0, int x1, int y0, int y1, int z0, int z1, int dimz) {

    const size_t i = (x0 + VEC_LEN * (threadIdx.x + blockIdx.x * blockDim.x));
    const size_t j = y0 + threadIdx.y + blockIdx.y * blockDim.y;

    if (i >= x1) return;

    const int kbegin = z0 + blockIdx.z * dimz;
    const int kend = kbegin + dimz > z1 ? z1 : kbegin + dimz;

    size_t pos = i + line * j + slice * kbegin;
    size_t slice_vec = slice >> VEC_EXP;
    ptrdiff_t line_vec = line >> VEC_EXP;

    p_in += pos - R * slice;
    p_out += pos;
    const vec *p_in_vec = reinterpret_cast<const vec*>(p_in);
    vec *p_out_vec = reinterpret_cast<vec*>(p_out);

    // LDS for y-direction stencil
    const int lds_y = BLOCK_DIM_Y + 2*R;
#if BLOCK_DIM_Y == RADIUS
    const int sj = y0 + threadIdx.y;
    size_t spos = threadIdx.x + sj * BLOCK_DIM_X;
#else
    size_t spos = threadIdx.x + (R + threadIdx.y) * BLOCK_DIM_X;
#endif
    __shared__ vec smem[BLOCK_DIM_X * lds_y];

    // Preload FD coefficients from constant memory into registers
    float coeff_x[2 * R + 1];
    float coeff_y[2 * R + 1];
    float coeff_z[2 * R + 1];
    for (int r = 0; r <= 2 * R; ++r) {
        coeff_x[r] = opt_v1_dx<R>[r];
        coeff_y[r] = opt_v1_dy<R>[r];
        coeff_z[r] = opt_v1_dz<R>[r];
    }
    #define COEFF_X(r) coeff_x[r]
    #define COEFF_Y(r) coeff_y[r]
    #define COEFF_Z(r) coeff_z[r]

    // Z-direction sliding window
    vec w[2 * R + 1];

    // Output registers
    vec out[R+1];

    // X-direction stencil window
    float x_win[2 * XWIN + VEC_LEN];
    vec *x_win_vec = reinterpret_cast<vec*>(x_win);

    // Prime the z sliding window
    for (int r = 0; r < R; ++r) {
        w[r] = p_in_vec[0];
        p_in_vec += slice_vec;
    }

    // Prime the xy pipeline
    for (int r = R; r < 2 * R; ++r) {
        for (int r2 = 0; r2 < 2*XWIN_VEC + 1; ++r2)
            x_win_vec[r2] = p_in_vec[0 - XWIN_VEC + r2];

        __syncthreads();
        {
#if BLOCK_DIM_Y == RADIUS
            smem[spos - (BLOCK_DIM_X * R)] = p_in_vec[0 - R * line_vec];
            smem[spos] = x_win_vec[XWIN_VEC];
            smem[spos + (BLOCK_DIM_X * BLOCK_DIM_Y)] = p_in_vec[0 + line_vec * BLOCK_DIM_Y];
#else
            smem[spos] = x_win_vec[XWIN_VEC];
            // Loop form in initial load (optimal for R=8)
            for (int h = threadIdx.y; h < R; h += BLOCK_DIM_Y) {
                smem[threadIdx.x + h * BLOCK_DIM_X] = p_in_vec[(h - int(threadIdx.y) - R) * line_vec];
                smem[threadIdx.x + (R + BLOCK_DIM_Y + h) * BLOCK_DIM_X] = p_in_vec[(BLOCK_DIM_Y + h - int(threadIdx.y)) * line_vec];
            }
#endif
        }
        __syncthreads();

        // Compute xy stencils with explicit FMA (inlined coefficients)
        out[r-R] = {0.0f};
        for (int r2 = 0; r2 <= 2 * R; ++r2) {
            for (int ii = 0; ii < VEC_LEN; ++ii) {
                out[r-R][ii] = __builtin_fmaf(smem[spos + (r2 - R) * BLOCK_DIM_X][ii], COEFF_Y(r2), out[r-R][ii]);
                out[r-R][ii] = __builtin_fmaf(x_win[XWIN_OFF + r2 + ii], COEFF_X(r2), out[r-R][ii]);
            }
        }

        w[r] = x_win_vec[XWIN_VEC];
        p_in_vec += slice_vec;
    }

    // Main z-loop with ILP + explicit FMA
    for (int k = kbegin; k < kend; ++k) {
        // Step 1: Load x into registers
        for (int r2 = 0; r2 < 2*XWIN_VEC+1; ++r2)
            x_win_vec[r2] = p_in_vec[0 - XWIN_VEC + r2];

        // Step 2: Update z sliding window
        w[2*R] = x_win_vec[XWIN_VEC];

        // Step 3: COMPUTE Z-STENCIL EARLY with explicit FMA (ILP: overlaps with barrier)
        vec lapz = {0.0f};
        for (int r = 0; r <= 2 * R; ++r) {
            for (int ii = 0; ii < VEC_LEN; ++ii) {
                lapz[ii] = __builtin_fmaf(w[r][ii], COEFF_Z(r), lapz[ii]);
            }
        }

        // Step 4: Load y into LDS (barrier here)
        __syncthreads();
        {
#if BLOCK_DIM_Y == RADIUS
            smem[spos - (BLOCK_DIM_X * R)] = p_in_vec[0 - R * line_vec];
            smem[spos] = x_win_vec[XWIN_VEC];
            smem[spos + (BLOCK_DIM_X * BLOCK_DIM_Y)] = p_in_vec[0 + line_vec * BLOCK_DIM_Y];
#else
            smem[spos] = x_win_vec[XWIN_VEC];
#if RADIUS == 8
            // Explicit unrolling for R=8 in main z-loop
            smem[threadIdx.x + threadIdx.y * BLOCK_DIM_X] = p_in_vec[-R * line_vec];
            smem[threadIdx.x + (R + BLOCK_DIM_Y + threadIdx.y) * BLOCK_DIM_X] = p_in_vec[BLOCK_DIM_Y * line_vec];
            smem[threadIdx.x + (threadIdx.y + BLOCK_DIM_Y) * BLOCK_DIM_X] = p_in_vec[(BLOCK_DIM_Y - R) * line_vec];
            smem[threadIdx.x + (R + 2*BLOCK_DIM_Y + threadIdx.y) * BLOCK_DIM_X] = p_in_vec[2*BLOCK_DIM_Y * line_vec];
#else
            for (int h = threadIdx.y; h < R; h += BLOCK_DIM_Y) {
                smem[threadIdx.x + h * BLOCK_DIM_X] = p_in_vec[(h - int(threadIdx.y) - R) * line_vec];
                smem[threadIdx.x + (R + BLOCK_DIM_Y + h) * BLOCK_DIM_X] = p_in_vec[(BLOCK_DIM_Y + h - int(threadIdx.y)) * line_vec];
            }
#endif
#endif
        }
        __syncthreads();

        // Step 5: Compute XY stencils with explicit FMA (needs LDS)
        vec lapxy = {0.0f};
        for (int r = 0; r <= 2 * R; ++r) {
            for (int ii = 0; ii < VEC_LEN; ++ii) {
                lapxy[ii] = __builtin_fmaf(smem[spos + (r - R) * BLOCK_DIM_X][ii], COEFF_Y(r), lapxy[ii]);
                lapxy[ii] = __builtin_fmaf(x_win[XWIN_OFF + r + ii], COEFF_X(r), lapxy[ii]);
            }
        }

        // Step 6: Combine z + xy
        for (int ii = 0; ii < VEC_LEN; ++ii) {
            out[0][ii] += lapz[ii];
        }
        out[R] = lapxy;

        // Step 7: Write result
        if (j < y1)
            ntstore(p_out_vec[0], out[0]);

        // Step 8: Shift sliding windows
        for (int r = 0; r < R; ++r)
            out[r] = out[r+1];
        for (int r = 0; r < 2*R; ++r)
            w[r] = w[r+1];

        p_in_vec += slice_vec;
        p_out_vec += slice_vec;
    }

#undef COEFF_X
#undef COEFF_Y
#undef COEFF_Z
}

template <int R>
void optimized_v1(float *p_out, const float *p_in, const float *d, int line, int
        slice, int x0, int x1, int y0, int y1, int z0, int z1, int dimwin) {

    dimwin = dimwin == -1 ? z1 - z0 : dimwin;

    dim3 block (BLOCK_DIM_X, BLOCK_DIM_Y);
    dim3 grid;
    grid.x = ceil(x1 - x0, VEC_LEN * block.x);
    grid.y = ceil(y1 - y0, block.y);
    grid.z = ceil(z1 - z0, dimwin);

    optimized_v1_kernel_impl<R><<<grid, block>>>(p_out, p_in, line, slice, x0, x1,
            y0, y1, z0, z1, dimwin);
    HIP_CHECK(hipGetLastError());
}

template void optimized_v1<RADIUS>(float*, const float*, const float*, int, int, int, int, int, int, int, int, int);

#undef BLOCK_DIM_X
#undef BLOCK_DIM_Y
