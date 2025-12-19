#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string.h>
#include <assert.h>
#include <hip/hip_runtime.h>

#ifndef RADIUS
#define RADIUS 4
#endif

#include "helper.hpp"
#include "compare.hpp"
#include "baseline_kernel.hpp"
#include "fd_coefficients.hpp"
#include "initialize.hpp"

void print_help(const char* program_name) {
    printf("Seismic Stencils - 3D Finite Difference Stencil Computation\n");
    printf("=============================================================\n\n");
    printf("Usage: %s [OPTIONS]\n\n", program_name);
    printf("Options:\n");
    printf("  -nx <val>     Grid size in x-direction          [default: 512]\n");
    printf("  -ny <val>     Grid size in y-direction          [default: 512]\n");
    printf("  -nz <val>     Grid size in z-direction          [default: 512]\n");
    printf("  -nt <val>     Number of iterations              [default: 100]\n");
    printf("  -nw <val>     Z-window size for sliding window  [default: radius-dependent]\n");
    printf("  -align <val>  Alignment factor for leading dim  [default: radius-dependent]\n");
    printf("  -h, --help    Show this help message\n\n");
    printf("Compile-time settings (current build):\n");
    printf("  RADIUS  = %d (stencil radius, %s order accuracy)\n", RADIUS,
           RADIUS == 1 ? "2nd" : RADIUS == 2 ? "4th" : RADIUS == 3 ? "6th" :
           RADIUS == 4 ? "8th" : RADIUS == 6 ? "12th" : RADIUS == 8 ? "16th" : "unknown");
    printf("  VEC_LEN = %d (vector length, 2^VEC_EXP)\n\n", VEC_LEN);
    printf("Examples:\n");
    printf("  %s                                    # 512^3 grid, 100 iterations\n", program_name);
    printf("  %s -nx 256 -ny 256 -nz 256            # 256^3 grid\n", program_name);
    printf("  %s -nx 1024 -ny 1024 -nz 1024 -nt 50 # 1024^3 grid, 50 iterations\n", program_name);
    printf("  %s -nx 512 -ny 512 -nz 512 -nw 128 -align 64\n\n", program_name);
    printf("Note: nx should be divisible by VEC_LEN (%d) for best performance.\n", VEC_LEN);
}

int main(int argc, char **argv) {

    // Default grid size and number of iterations
    int nx = 512;
    int ny = 512;
    int nz = 512;
    int nt = 100;
    // Default nw (z-window size) and align (leading dimension alignment)
    // align defaults to BLOCK_DIM_X for optimal memory access
    int nw = (RADIUS <= 4) ? 128 : (RADIUS <= 6) ? 80 : 160;
    int align = (RADIUS == 1) ? 256 : (RADIUS == 2) ? 128 : 64;

    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-nx") == 0 && i + 1 < argc) {
            nx = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-ny") == 0 && i + 1 < argc) {
            ny = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-nz") == 0 && i + 1 < argc) {
            nz = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-nt") == 0 && i + 1 < argc) {
            nt = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-nw") == 0 && i + 1 < argc) {
            nw = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-align") == 0 && i + 1 < argc) {
            align = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_help(argv[0]);
            return 0;
        } else {
            printf("Unknown option: %s\n", argv[i]);
            print_help(argv[0]);
            return 1;
        }
    }
    
    printf("Settings: nx = %d ny = %d nz = %d nt = %d nw = %d align = %d\n",
            nx, ny, nz, nt, nw, align);

    // Finite difference coefficients
    const int R = RADIUS;
    float h = 1.0;
    float d[2 * R + 1];

    // Finite difference approximation to use
    switch (R) {
        case 1:
            // 2nd order
            fd_coefficients_d2_order_2(d, R, h);
            break;
        case 2:
            // 4th order
            fd_coefficients_d2_order_4(d, R, h);
            break;
        case 3:
            // 6th order
            fd_coefficients_d2_order_6(d, R, h);
            break;
        case 4:
            // 8th order
            fd_coefficients_d2_order_8(d, R, h);
            break;
        case 6:
            // 12th order
            fd_coefficients_d2_order_12(d, R, h);
            break;
        case 8:
            // 16th order
            fd_coefficients_d2_order_16(d, R, h);
            break;
        default:
            printf("ERROR: unsupported radius (-DRADIUS=%d), exiting program...\n", R);
            exit(0);
    }

    // Padded grid size
    int mx = nx + 2 * XWIN;
    // Make the leading dimension a multiple of align
    // Align is typically chosen as a multiple of the cache line size (64 B)
    mx = ceil(mx, align) * align;
    const int my = ny + 2 * R;
    const int mz = nz + 2 * R;

    const int line = mx;
    const int slice = mx * my;

    // Total number of grid points
    const size_t m = (size_t)mx * (size_t)my * (size_t)mz;

    // Align the offset, if any
    int offset = align ? align - XWIN : 0;
    offset = offset < 0 ? 0 : offset;

    float *d_p_in; alloc(&d_p_in, m, offset);
    float *d_p_out; alloc(&d_p_out, m, offset);
    float *d_p_ref; alloc(&d_p_ref, m, offset);
    
    // Initialize input, output, and reference arrays
    // initialize a polynomial of the form: a.x * x^s.x + a.y * y^s.y + a.z * z^s.z
    // i.e., x + y + z
    float3 a = {1.0f, 1.0f, 1.0f};
    float3 s = {1.0f, 1.0f, 1.0f};
    initialize_polynomial(d_p_in, a, s, 0, mx, 0, my, 0, mz, line, slice);
    HIP_CHECK(hipMemset(d_p_out, 0, m * sizeof(float)));
    HIP_CHECK(hipMemset(d_p_ref, 0, m * sizeof(float)));

    HIP_CHECK(hipDeviceSynchronize());

    init_fd_xy_gpu<R>(d, d);
    init_fd_z_gpu<R>(d);
    
    // Performance metrics
    float total_elapsed, elapsed, maxval;
    hipEvent_t start, stop;
    double theoretical_fetch_xyz = nx * ny * nz + 2 * R * (ny * nz + nx * nz + nx * ny);
    double theoretical_write = nx * ny * nz;
    double total_bytes_xyz = (theoretical_fetch_xyz + theoretical_write) * sizeof(float); // GB 
    HIP_CHECK( hipEventCreate(&start) );
    HIP_CHECK( hipEventCreate(&stop)  );
    
    printf("\nApplying baseline stencil kernel...\n");
    if (nx % VEC_LEN) {
        printf("\nWarning: nx = %d not divisible by vec length %d, skipping\n", nx, VEC_LEN);
    } else {
        HIP_CHECK(hipMemset(d_p_out, 0, m * sizeof(float)));
        total_elapsed = 0;
        for (int i = 0; i < nt; ++i) {
            HIP_CHECK( hipDeviceSynchronize()                     );
            HIP_CHECK( hipEventRecord(start)                      );
            baseline<R>(d_p_out, d_p_in, d, line, slice, XWIN, nx + XWIN, R, ny + R, R, R + nz, nw);
            HIP_CHECK( hipEventRecord(stop)                       );
            HIP_CHECK( hipEventSynchronize(stop)                  );
            HIP_CHECK( hipEventElapsedTime(&elapsed, start, stop) );
            total_elapsed += elapsed;
        }
        maxval = max_absval_gpu(d_p_out, d_p_ref, line, slice, XWIN, nx + XWIN, R, ny + R, R, nz + R);
        printf("\tMaximum absolute pointwise difference: %g\n", maxval);
        printf("\tAverage kernel time: %g ms\n", total_elapsed / nt);
        printf("\tEffective memory bandwidth %g GB/s\n", total_bytes_xyz * nt / total_elapsed / 1e6);
    }

    dealloc(d_p_in, offset);
    dealloc(d_p_out, offset);
    dealloc(d_p_ref, offset);

    return 0;
}

