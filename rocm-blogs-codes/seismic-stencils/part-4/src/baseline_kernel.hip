#include "helper.hpp"
#include "baseline_kernel.hpp"

// Constant memory for finite difference coefficients
template <int R> __constant__ float d_dx[2 * R + 1];
template <int R> __constant__ float d_dy[2 * R + 1];
template <int R> __constant__ float d_dz[2 * R + 1];

// Initialize constant memory arrays
template <int R>
void init_fd_xy_gpu(const float *dx, const float *dy) {
    HIP_CHECK(hipMemcpyToSymbol(d_dx<R>, dx, (2 * R + 1) * sizeof(float)));
    HIP_CHECK(hipMemcpyToSymbol(d_dy<R>, dy, (2 * R + 1) * sizeof(float)));
}

template <int R>
void init_fd_z_gpu(const float *d) {
    HIP_CHECK(hipMemcpyToSymbol(d_dz<R>, d, (2 * R + 1) * sizeof(float)));
}

// Explicit instantiations for the init functions
template void init_fd_xy_gpu<RADIUS>(const float*, const float*);
template void init_fd_z_gpu<RADIUS>(const float*);

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
#define BLOCK_DIM_X 64
#define BLOCK_DIM_Y 4
#else
#error "Unsupported RADIUS value"
#endif

template <int R>
__launch_bounds__((BLOCK_DIM_X) * (BLOCK_DIM_Y))
__global__ void baseline_kernel_impl(float *p_out, const float *p_in, int line, int
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

    for (int r = R; r < 2 * R; ++r) {
        // Load x into registers
        for (int r2 = 0; r2 < 2*XWIN_VEC + 1; ++r2)
            x_win_vec[r2] = p_in_vec[0 - XWIN_VEC + r2];

        // Load y into LDS
        __syncthreads();
        {
#if BLOCK_DIM_Y == RADIUS
            smem[spos - (BLOCK_DIM_X * R)] = p_in_vec[0 - R * line_vec];
            smem[spos] = x_win_vec[XWIN_VEC];
            smem[spos + (BLOCK_DIM_X * BLOCK_DIM_Y)] = p_in_vec[0 + line_vec * BLOCK_DIM_Y];
#else
            smem[spos] = x_win_vec[XWIN_VEC];
            // Load y-direction halos
            for (int h = threadIdx.y; h < R; h += BLOCK_DIM_Y) {
                int y_off_lo = (h - int(threadIdx.y) - R);
                int y_off_hi = (BLOCK_DIM_Y + h - int(threadIdx.y));
                smem[threadIdx.x + h * BLOCK_DIM_X] = p_in_vec[y_off_lo * line_vec];
                smem[threadIdx.x + (R + BLOCK_DIM_Y + h) * BLOCK_DIM_X] = p_in_vec[y_off_hi * line_vec];
            }
#endif
        }
        __syncthreads();

        // Compute xy stencils
        out[r-R] = {0.0f};
        for (int r2 = 0; r2 <= 2 * R; ++r2) {
            out[r-R] += smem[spos + (r2 - R) * BLOCK_DIM_X] * d_dy<R>[r2];
            for (int ii = 0; ii < VEC_LEN; ++ii)
                out[r-R][ii] += x_win[XWIN_OFF + r2 + ii] * d_dx<R>[r2];
        }

        w[r] = x_win_vec[XWIN_VEC];
        p_in_vec += slice_vec;
    }

    // Main z-loop
    for (int k = kbegin; k < kend; ++k) {
        // Load x into registers
        for (int r2 = 0; r2 < 2*XWIN_VEC+1; ++r2)
            x_win_vec[r2] = p_in_vec[0 - XWIN_VEC + r2];

        // Load y into LDS
        __syncthreads();
        {
#if BLOCK_DIM_Y == RADIUS
            smem[spos - (BLOCK_DIM_X * R)] = p_in_vec[0 - R * line_vec];
            smem[spos] = x_win_vec[XWIN_VEC];
            smem[spos + (BLOCK_DIM_X * BLOCK_DIM_Y)] = p_in_vec[0 + line_vec * BLOCK_DIM_Y];
#else
            smem[spos] = x_win_vec[XWIN_VEC];
            // Load y-direction halos
            for (int h = threadIdx.y; h < R; h += BLOCK_DIM_Y) {
                int y_off_lo = (h - int(threadIdx.y) - R);
                int y_off_hi = (BLOCK_DIM_Y + h - int(threadIdx.y));
                smem[threadIdx.x + h * BLOCK_DIM_X] = p_in_vec[y_off_lo * line_vec];
                smem[threadIdx.x + (R + BLOCK_DIM_Y + h) * BLOCK_DIM_X] = p_in_vec[y_off_hi * line_vec];
            }
#endif
        }
        __syncthreads();

        // Compute xyz stencils
        w[2*R] = x_win_vec[XWIN_VEC];
        out[R] = {0.0f};
        for (int r = 0; r <= 2 * R; ++r) {
            out[0] += w[r] * d_dz<R>[r];
            out[R] += smem[spos + (r - R) * BLOCK_DIM_X] * d_dy<R>[r];
            for (int ii = 0; ii < VEC_LEN; ++ii)
                out[R][ii] += x_win[XWIN_OFF + r + ii] * d_dx<R>[r];
        }

        // Write result
        if (j < y1)
            ntstore(p_out_vec[0], out[0]);

        // Shift sliding windows
        for (int r = 0; r < R; ++r)
            out[r] = out[r+1];
        for (int r = 0; r < 2*R; ++r)
            w[r] = w[r+1];

        p_in_vec += slice_vec;
        p_out_vec += slice_vec;
    }
}

template <int R>
void baseline(float *p_out, const float *p_in, const float *d, int line, int
        slice, int x0, int x1, int y0, int y1, int z0, int z1, int dimwin) {

    dimwin = dimwin == -1 ? z1 - z0 : dimwin;

    dim3 block (BLOCK_DIM_X, BLOCK_DIM_Y);
    dim3 grid;
    grid.x = ceil(x1 - x0, VEC_LEN * block.x);
    grid.y = ceil(y1 - y0, block.y);
    grid.z = ceil(z1 - z0, dimwin);

    baseline_kernel_impl<R><<<grid, block>>>(p_out, p_in, line, slice, x0, x1,
            y0, y1, z0, z1, dimwin);
    HIP_CHECK(hipGetLastError());
}

template void baseline<RADIUS>(float*, const float*, const float*, int, int, int, int, int, int, int, int, int);

#undef BLOCK_DIM_X
#undef BLOCK_DIM_Y
