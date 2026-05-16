# simplicialdepth

An R package for computing simplicial depth in R² and R³.

## Features

- **2D Simplicial Depth**: Computing depth of points in the plane (all-points and query-based)
- **3D Simplicial Depth**: Fast O(n²) algorithm for 3D depth computation
- **Angular Simplicial Depth**: Directional depth computation
- **Cross-platform Support**: macOS, Linux, and Windows
- **Optimized Performance**: Compiled with -O3 optimization on all platforms

## Installation

```r
# Install from GitHub
devtools::install_github("lecitin/simplicial-depth")

# Or build from source
R CMD INSTALL .
```

## Quick Start

```r
library(simplicialdepth)

# 2D example: compute depths at fixed query points
X <- matrix(rnorm(100), ncol=2)
query_points <- matrix(c(0, 0, 1, 1), ncol=2, byrow=TRUE)
depths_2d <- simplicialdepth(X, query_points)

# 3D example
X3 <- matrix(rnorm(300), ncol=3)
q3 <- c(0, 0, 0)
depth_3d <- simplicialdepth(X3, q3)
```

## Build System

### Compilation Flags

The package is compiled with optimizations for performance and code quality:

- **`-O3`**: Maximum compiler optimization for fastest runtime execution
- **`-Wall`**: Compiler warnings to catch potential bugs (clean compilation, no warnings)

#### Optional: Link-Time Optimization (LTO)

For 5-15% additional performance improvement (at the cost of 2-5x longer build time), enable Link-Time Optimization:

```bash
# Set environment variable before installation
export PKG_CXXFLAGS="-O3 -Wall -flto"
R CMD INSTALL .

# Or edit src/Makevars.in and uncomment the -flto section
```

**Requirements for LTO:**
- GCC 4.8+ or Clang 3.5+ (Linux/macOS)
- MSVC 2015+ (Windows)

### Cross-Platform Support

The build system automatically detects and configures:

- **macOS**: Clang compiler with proper dylib runtime paths
- **Linux**: GCC/Clang with ELF runtime paths
- **Windows**: MSVC or MinGW with RcppParallel integration

**Build scripts:**
- `configure` - Unix/Linux/macOS configuration
- `configure.win` - Windows configuration (with graceful fallback)
- `src/Makevars.in` - Unix/Linux/macOS template
- `src/Makevars.win.in` - Windows template

### Continuous Integration

GitHub Actions automatically test the package on:

- Windows Server 2022 (latest R)
- macOS 12+ (latest R)
- Ubuntu 22.04 LTS (latest R and R-devel)

Each test run performs:
- Full R CMD check (CRAN standards)
- Complete test suite (16 tests)
- Code coverage measurement (optional)

View CI status and logs: https://github.com/lecitin/simplicial-depth/actions

## Testing

Run the full test suite:

```r
devtools::test()
```

Or in the shell:

```bash
R CMD check .
```

All 16 tests should pass:
- 2 angularSD 2D tests
- 3 angularSD 3D tests  
- 10 simplicialdepth tests
- 1 spherical ASD consistency test

## Performance Notes

### Compilation

- Default build (with -O3): ~30-60 seconds
- With LTO enabled (-flto): ~2-5 minutes
- Build time varies by hardware and compiler version

### Runtime

- 2D queries: O(n² log n) for all-points mode
- 3D queries: O(n²) for query-based depth (major improvement over brute force)
- Parallel execution: Automatic via RcppParallel for multi-core systems

## Compiler Warnings

The package is compiled with `-Wall` (all common warnings enabled) and achieves **zero warnings** on:
- Apple Clang 17.0+
- GCC 11+
- Clang 14+

This ensures code quality and early detection of potential issues.

## Dependencies

- **Rcpp**: C++ to R bindings
- **RcppParallel**: Multi-core parallel computing with TBB
- **R >= 3.5.0**: For package infrastructure

## Development

### Building Locally

```bash
# Configure (Unix-like systems)
./configure

# Clean build artifacts
./cleanup

# Install with custom flags
PKG_CXXFLAGS="-O3 -Wall" R CMD INSTALL .
```

### Modifying Build Flags

Edit `src/Makevars.in` (Unix) or `src/Makevars.win.in` (Windows) to adjust compiler flags.

Available options:
```makefile
# Optimization levels
PKG_CXXFLAGS = -O0 -O1 -O2 -O3  # Choose one

# Additional optimizations
PKG_CXXFLAGS = -O3 -march=native  # Optimize for local machine (not portable!)
PKG_CXXFLAGS = -O3 -flto           # Link-time optimization (slower build)

# Strict warnings
PKG_CXXFLAGS = -O3 -Wall -Werror   # Fail on any warnings
```

## Algorithm Overview

### 2D Simplicial Depth

Counts the number of triangles formed by sample points that contain each query point.

- All-points mode: O(n³) via brute force
- Query-based: O(n² log n) via topological sweep

### 3D Simplicial Depth

O(n²) algorithm counting containing tetrahedra:

1. Project sample points onto unit sphere centered at query point
2. Gnomonic projection to plane z=1
3. Color points blue (upper hemisphere) or red (lower hemisphere)
4. Count:
   - N1: Red points inside blue triangles
   - N2: Blue points inside red triangles
   - N3: Red-blue segment crossings (using topological sweep)
5. Depth = N1 + N2 + N3

**Reference**: See `src/sd_3d.md` for mathematical details.

## Citation

If you use this package in published research, please cite:

```bibtex
@software{simplicialdepth2024,
  title={simplicialdepth: Simplicial Depth in R² and R³},
  author={Author Name},
  year={2024},
  url={https://github.com/lecitin/simplicial-depth}
}
```

## License

[License information to be added]

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new features
4. Ensure all tests pass with clean compilation
5. Submit a pull request

The CI/CD pipeline will automatically test your changes on Windows, macOS, and Linux.

## Troubleshooting

### Build fails on Windows

If the configure.win script fails, the build will fall back to a pre-generated configuration. 
Check that RcppParallel is installed:

```r
install.packages("RcppParallel")
```

### Compiler warnings during build

The package is tested to compile with zero warnings. If you see warnings:

1. Ensure you have a recent compiler (GCC 11+, Clang 14+, MSVC 2015+)
2. Check the build output for details
3. Report the issue on GitHub

### Tests fail

Ensure dependencies are installed:

```r
install.packages(c("Rcpp", "RcppParallel", "testthat", "devtools"))
```

Then run tests with verbose output:

```r
devtools::test(reporter = "progress")
```

## Related Publications

- [Add references to research papers using simplicial depth]

## Contact

[Maintainer contact information]
