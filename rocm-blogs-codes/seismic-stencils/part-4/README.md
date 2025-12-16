# Seismic Stencils - Part 4

High-performance 3D finite difference stencil computation on AMD GPUs using HIP. This implementation uses vectorized LDS (Local Data Share) with a sliding window approach for optimal memory bandwidth utilization.

## Features

- Vectorized computation using HIP vector types
- LDS-based data sharing for y-direction stencil
- Sliding window approach for z-direction
- Support for multiple stencil radii (2nd to 16th order accuracy)
- Auto-detection of GPU architecture

## Requirements

- ROCm with HIP support
- AMD GPU (MI250X, MI300X, RX 7000 series, etc.)

## Building

### Basic Build (Auto-detect GPU)

```bash
make
```

By default, this builds all 18 configurations (6 radii × 3 vec values). Binaries are placed in the `bin/` directory.

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
| `rad`        | Stencil radius (determines order)    | all     | 1, 2, 3, 4, 6, 8 |
| `vec`        | Vector exponent (VEC_LEN = 2^vec)    | all     | 0, 1, 2 |
| `SAVE_TEMPS` | Save intermediate files (assembly, IR) | off   | 1 (to enable) |

**Stencil Radius to Order Mapping:**
- `rad=1`: 2nd order accuracy
- `rad=2`: 4th order accuracy
- `rad=3`: 6th order accuracy
- `rad=4`: 8th order accuracy
- `rad=6`: 12th order accuracy
- `rad=8`: 16th order accuracy

**Vector Length:**
- `vec=0`: VEC_LEN = 1 (scalar)
- `vec=1`: VEC_LEN = 2
- `vec=2`: VEC_LEN = 4

### Build Examples

```bash
# Build all 18 configurations (6 radii × 3 vecs)
make

# Build all radii with vec=2 (6 binaries)
make vec=2

# Build all vecs with rad=4 (3 binaries)
make rad=4

# Build specific configuration
make rad=4 vec=2

# Build for specific GPU architecture
make ARCH=gfx90a rad=6 vec=2

# Build with assembly/IR output for debugging
make SAVE_TEMPS=1

# Clean everything
make clean

# Clean specific configurations
make clean rad=4 vec=2    # Clean only R4 vec2
make clean rad=4          # Clean all vecs for R4
make clean vec=2          # Clean all rads for vec2
```

## Running

### Basic Usage

Binaries are located in the `bin/` directory. The binary name includes the configuration: `baseline_R<rad>_vec<vec>_<arch>.x`

```bash
bin/baseline_R4_vec2_gfx90a.x [OPTIONS]
```

### Runtime Arguments

All arguments are optional and use flag-based syntax:

| Flag      | Description                    | Default |
|-----------|--------------------------------|---------|
| `-nx`     | Grid size in x-direction       | 512     |
| `-ny`     | Grid size in y-direction       | 512     |
| `-nz`     | Grid size in z-direction       | 512     |
| `-nt`     | Number of iterations           | 100     |
| `-nw`     | Z-window size for sliding window | radius-dependent |
| `-align`  | Alignment factor for leading dimension | radius-dependent |
| `-h, --help` | Show help message            | -       |

### Runtime Examples

```bash
# Show help
bin/baseline_R4_vec2_gfx90a.x -h
bin/baseline_R4_vec2_gfx90a.x --help

# Default 512x512x512 grid, 100 iterations
bin/baseline_R4_vec2_gfx90a.x

# Custom grid size
bin/baseline_R4_vec2_gfx90a.x -nx 256 -ny 256 -nz 256

# 1024^3 grid with 50 iterations
bin/baseline_R4_vec2_gfx90a.x -nx 1024 -ny 1024 -nz 1024 -nt 50

# Full customization
bin/baseline_R4_vec2_gfx90a.x -nx 512 -ny 512 -nz 512 -nt 100 -nw 128 -align 64

# Flags can be specified in any order
bin/baseline_R4_vec2_gfx90a.x -nt 50 -nx 256 -nz 256 -ny 256
```

## Output

The program outputs:
- Configuration settings (grid size, iterations, alignment)
- Maximum absolute pointwise difference (for correctness verification)
- Average kernel execution time (ms)
- Effective memory bandwidth (GB/s)

## Performance Tips

1. **Grid Size**: Use grid sizes divisible by VEC_LEN for best performance
2. **Alignment**: Default alignment values are optimized per radius. For R=4, use `-align 64`; for R=6, use `-align 1`; for R=8, use `-align 64`
3. **Vector Length**: Higher vector lengths may improve throughput but require nx divisible by VEC_LEN
4. **Window Size (nw)**: Default values are optimized per radius (R=4: 128, R=6: 80, R=8: 160). Tune based on GPU memory and register pressure
5. **Build System**: Use filtered builds (`make rad=4 vec=2`) to build only the configurations you need, reducing build time

## Project Structure

```
part-4/
├── bin/                            # All compiled binaries
│   ├── baseline_R1_vec0_gfx90a.x
│   ├── baseline_R1_vec1_gfx90a.x
│   ├── baseline_R1_vec2_gfx90a.x
│   └── ...                         # 18 total configurations
├── build/                          # Build artifacts
│   ├── rad1/                       # Organized by radius
│   │   ├── vec0/
│   │   │   └── gfx90a/             # Architecture-specific objects
│   │   │       ├── main.o
│   │   │       └── baseline_kernel.o
│   │   ├── vec1/
│   │   └── vec2/
│   ├── rad2/
│   └── ...                         # rad3, rad4, rad6, rad8
├── include/                        # Header files
│   ├── baseline_kernel.hpp         # Kernel function declarations
│   ├── compare.hpp                 # Comparison utilities
│   ├── fd_coefficients.hpp         # Finite difference coefficients
│   ├── helper.hpp                  # Utility macros and functions
│   └── initialize.hpp              # Grid initialization functions
├── src/                            # Source files
│   ├── baseline_kernel.hip         # Baseline kernel implementation
│   └── main.cpp                    # Main driver program
├── Makefile                        # Build configuration
└── README.md                       # This file
```

Build directories are organized hierarchically as `build/rad{number}/vec{number}/{arch}/`, allowing multiple configurations to coexist. All binaries are placed in the `bin/` directory.

## License

See the main repository for license information.

