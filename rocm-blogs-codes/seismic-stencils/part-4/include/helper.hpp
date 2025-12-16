#pragma once
#include <math.h>
#include <iostream>
#include <string.h>
#include <assert.h>
#include <hip/hip_runtime.h>


#define POS(i,j,k) ( (size_t)(i) + ((size_t)line * (size_t)(j)) + ((size_t)slice * (size_t)(k)) )

// HIP error check
#define HIP_CHECK(stat)                                           \
{                                                                 \
    if(stat != hipSuccess)                                        \
    {                                                             \
        std::cerr << "HIP error: " << hipGetErrorString(stat) <<  \
        " in file " << __FILE__ << ":" << __LINE__ << std::endl;  \
        exit(-1);                                                 \
    }                                                             \
}

// Initialization functions for finite difference coefficients
// Declared here, defined in baseline_kernel.hip
template <int R>
void init_fd_xy_gpu(const float *dx, const float *dy);

template <int R>
void init_fd_z_gpu(const float *d);

// Vectorized floats
#ifndef VEC_EXP
#define VEC_EXP 0
#endif
#define VEC_LEN (1 << (VEC_EXP))
using vec = __attribute__((__vector_size__(VEC_LEN * sizeof(float)))) float;

// x window stencil
// XWIN must be >= RADIUS and a multiple of VEC_LEN for correct stencil computation
// XWIN_VEC = XWIN / VEC_LEN (number of vector loads needed)
#if VEC_LEN == 1
#define XWIN RADIUS
#define XWIN_VEC RADIUS
#elif VEC_LEN == 2
#if RADIUS <= 2
#define XWIN 2
#define XWIN_VEC 1
#elif RADIUS <= 4
#define XWIN 4
#define XWIN_VEC 2
#elif RADIUS <= 6
#define XWIN 6
#define XWIN_VEC 3
#else
#define XWIN 8
#define XWIN_VEC 4
#endif
#elif VEC_LEN == 4
#if RADIUS <= 4
#define XWIN 4
#define XWIN_VEC 1
#elif RADIUS <= 8
#define XWIN 8
#define XWIN_VEC 2
#else
#define XWIN 12
#define XWIN_VEC 3
#endif
#else
// VEC_LEN >= 8
#if RADIUS <= VEC_LEN
#define XWIN VEC_LEN
#define XWIN_VEC 1
#else
#define XWIN (2 * VEC_LEN)
#define XWIN_VEC 2
#endif
#endif
#define XWIN_OFF (XWIN-RADIUS)

// Nontemporal stores
#define ntload(ref) __builtin_nontemporal_load(&ref)
#define ntstore(ref, val) __builtin_nontemporal_store(val, &ref)

// Integer division: round up
__inline__ int ceil(int x, int y) {
    assert(y > 0);
    return (x - 1) / y + 1;
}

template <typename T>
void alloc(T **ptr, size_t num_elements, size_t offset) {
    HIP_CHECK(hipMalloc(ptr, (num_elements + offset) * sizeof(T)));
    *ptr += offset;
}

template <typename T>
void dealloc(T *ptr, size_t offset) {
    ptr -= offset;
    HIP_CHECK(hipFree(ptr));
}

template <typename T>
void print_array(T *arr, int line, int slice, int i0, int i1, int j0, int j1, int k0, int k1) {

    for (int k = k0; k < k1; ++k) {
        printf("k = %d \n", k);
        for (int j = j0; j < j1; ++j) {
            for (int i = i0; i < i1; ++i)
                printf("%2.2f ", arr[i + line * j + slice * k]);
            printf("\n");
        }
    }

}

