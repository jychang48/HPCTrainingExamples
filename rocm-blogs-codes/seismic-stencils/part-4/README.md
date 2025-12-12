# Seismic Stencils - Part 4

High-performance 3D finite difference stencil computation on AMD GPUs using HIP. This implementation uses vectorized LDS (Local Data Share) with a sliding window approach for optimal memory bandwidth utilization.

## Features

- Vectorized computation using HIP vector types
- LDS-based data sharing for y-direction stencil
- Sliding window approach for z-direction
- Support for multiple stencil radii (2nd to 8th order accuracy)
- Auto-detection of GPU architecture

## Requirements

- ROCm with HIP support
- AMD GPU (MI250X, MI300X, RX 7000 series, etc.)

## Building

### Basic Build (Auto-detect GPU)

```bash
make
```

The Makefile will automatically detect your GPU architecture using `rocminfo`.

### Specify GPU Architecture

```bash
make ARCH=gfx90a      # MI250X
make ARCH=gfx942      # MI300X
make ARCH=gfx1100     # RX 7900 XT
make ARCH=gfx1201     # Strix Point
```

### Compile-Time Options

| Option       | Description                          | Default | Valid Values |
|--------------|--------------------------------------|---------|--------------|
| `ARCH`       | GPU architecture                     | auto    | gfx90a, gfx942, gfx1100, etc. |
| `radius`     | Stencil radius (determines order)    | 4       | 1, 2, 3, 4   |
| `vec`        | Vector exponent (VEC_LEN = 2^vec)    | 2       | 0, 1, 2, ... |
| `SAVE_TEMPS` | Save intermediate files (assembly, IR) | off   | 1 (to enable) |

**Stencil Radius to Order Mapping:**
- `radius=1`: 2nd order accuracy
- `radius=2`: 4th order accuracy
- `radius=3`: 6th order accuracy
- `radius=4`: 8th order accuracy

**Vector Length:**
- `vec=0`: VEC_LEN = 1 (scalar)
- `vec=1`: VEC_LEN = 2
- `vec=2`: VEC_LEN = 4 (default)
- `vec=3`: VEC_LEN = 8

### Build Examples

```bash
# 8th order stencil with vec length 4 (default)
make

# 4th order stencil with vec length 2
make radius=2 vec=1

# 6th order stencil for MI250X
make ARCH=gfx90a radius=3

# Build with assembly/IR output for debugging
make SAVE_TEMPS=1

# Clean build artifacts
make clean
```

## Running

### Basic Usage

The binary name includes the configuration: `baseline_R<radius>_vec<vec>_<arch>.x`

```bash
./baseline_R4_vec2_gfx1201.x [nx] [ny] [nz] [nt] [nw] [align]
```

### Runtime Arguments

| Argument | Position | Description                    | Default |
|----------|----------|--------------------------------|---------|
| `nx`     | 1        | Grid size in x-direction       | 512     |
| `ny`     | 2        | Grid size in y-direction       | 512     |
| `nz`     | 3        | Grid size in z-direction       | 512     |
| `nt`     | 4        | Number of iterations           | 100     |
| `nw`     | 5        | Z-window size for sliding window | 1     |
| `align`  | 6        | Alignment factor for leading dimension | 1 |

### Runtime Examples

```bash
# Show help (example with gfx1201 architecture)
./baseline_R4_vec2_gfx1201.x --help
./baseline_R4_vec2_gfx1201.x -h

# Default 512x512x512 grid, 100 iterations
./baseline_R4_vec2_gfx1201.x

# Custom grid size
./baseline_R4_vec2_gfx1201.x 256 256 256

# 1024^3 grid with 50 iterations
./baseline_R4_vec2_gfx1201.x 1024 1024 1024 50

# Full customization
./baseline_R4_vec2_gfx1201.x 512 512 512 100 1 64

# Example with MI250X binary
./baseline_R4_vec2_gfx90a.x 512 512 512
```

## Output

The program outputs:
- Configuration settings (grid size, iterations, alignment)
- Maximum absolute pointwise difference (for correctness verification)
- Average kernel execution time (ms)
- Effective memory bandwidth (GB/s)

## Performance Tips

1. **Grid Size**: Use grid sizes divisible by VEC_LEN for best performance
2. **Alignment**: Set `align=64` or `align=128` for cache-line aligned access
3. **Vector Length**: Higher vector lengths may improve throughput but require nx divisible by VEC_LEN
4. **Window Size (nw)**: Tune based on GPU memory and register pressure

## Project Structure

```
part-4/
├── build/                          # Build artifacts
│   └── R4_vec2_gfx1201/            # Subdirectory per configuration
│       ├── main.o                  # Object files
│       ├── baseline_kernel.o
│       ├── *-hip-*.s               # GPU assembly (temp files)
│       └── *-host-*.s              # Host assembly (temp files)
├── include/                        # Header files
│   ├── baseline_kernel.hpp         # Kernel function declarations
│   ├── compare.hpp                 # Comparison utilities
│   ├── fd_coefficients.hpp         # Finite difference coefficients
│   ├── helper.hpp                  # Utility macros and functions
│   └── initialize.hpp              # Grid initialization functions
├── src/                            # Source files
│   ├── baseline_kernel.hip         # Baseline kernel implementation
│   └── main.cpp                    # Main driver program
├── baseline_R4_vec2_gfx1201.x      # Output binary (config in name)
├── Makefile                        # Build configuration
└── README.md                       # This file
```

Build subdirectories are organized by configuration (`R<radius>_vec<vec>_<arch>/`), allowing multiple configurations to coexist.

## License

See the main repository for license information.

