#pragma once

/**
 * @brief Baseline 3D finite difference stencil kernel
 * 
 * Computes a central high order finite difference (FD) approximation using
 * vectorized LDS (Local Data Share) with a sliding window approach.
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
void baseline(float *p_out, const float *p_in, const float *d, int line, int
        slice, int x0, int x1, int y0, int y1, int z0, int z1, int dimwin=-1);

