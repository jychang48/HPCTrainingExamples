#pragma once

/**
 * @brief Initialize finite difference coefficients for optimized_v2 kernel
 *
 * @tparam R Stencil radius
 * @param d Array of length 2*R+1 containing FD coefficients
 */
template <int R>
void init_fd_optimized_v2(const float *d);

/**
 * @brief Optimized v2 3D finite difference stencil kernel
 *
 * Optimizations over baseline:
 * - Improved ILP via loop unrolling (process 2 z-levels per iteration)
 * - Software prefetching for next iteration's data
 * - Better register utilization for reused data
 * - Reduced memory traffic through register blocking
 *
 * @tparam R Stencil radius (determines order of accuracy: 2*R order)
 * @param p_out Array to write the result to
 * @param p_in Array to apply the approximation to
 * @param d Array of length 2*R+1 containing FD coefficients (scaled by grid spacing)
 * @param line Line stride (leading dimension)
 * @param slice Slice stride (line * ny)
 * @param x0, x1 X-direction bounds [x0, x1)
 * @param y0, y1 Y-direction bounds [y0, y1)
 * @param z0, z1 Z-direction bounds [z0, z1)
 * @param dimwin Number of grid points to cover in a sliding window (-1 for full z-extent)
 */
template <int R>
void optimized_v2(float *p_out, const float *p_in, const float *d, int line, int
        slice, int x0, int x1, int y0, int y1, int z0, int z1, int dimwin=-1);

