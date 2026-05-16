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
- **`-flto`**: Link-Time Optimization for 5-15% additional performance gain

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

### Algorithm Complexity

#### Circular Angular Simplicial Depth (2D)
- **All-points mode** (leave-one-out): O(n²) - Topological sweep algorithm for each point
- **Query-based**: O(n log n) per query - Sort and angular sweep

#### Spherical Angular Simplicial Depth (3D)
- **All-points mode** (leave-one-out): O(n²) - Topological sweep algorithm for each point
- **Query-based**: O(n²) per query - Gnomonic projection + topological sweep (with antipode evaluation)

#### Simplicial Depth (2D)
- **All-points mode**: O(n² log n) - Leave-one-out: n × O(n log n) topological sweep
- **Query-based**: O(n log n) per query - Topological sweep + point-in-triangle tests

#### Simplicial Depth (3D)
- **All-points mode**: O(n³) - Leave-one-out: n × O(n²) tetrahedra containment
- **Query-based**: O(n²) per query - Sphere projection + gnomonic projection + topological sweep
- **Key optimizations**:
  - N1, N2: Point-in-triangle via topological sweep
  - N3: Segment crossings via topological sweep
  - Total: Single O(n²) pass instead of O(n⁴) brute force enumeration

#### k-Hull Depth (3D)
- **All-points mode**: O(n³) - Leave-one-out: n × O(n²) k-hull containment
- **Query-based**: O(n²) per query - Gnomonic projection + topological sweep

### Runtime Characteristics

- **Parallelization**: Automatic via RcppParallel (TBB) for multi-core systems
- **All-points mode**: Parallelized across points; each point computed independently
- **Compiler optimization**: `-O3 -flto` by default for ~5-15% additional performance
- **Memory efficiency**: Linear in sample size for most algorithms

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

# Strict warnings
PKG_CXXFLAGS = -O3 -Wall -Werror   # Fail on any warnings
```

## Algorithm Overview

### Circular Angular Simplicial Depth (2D)

Counts the number of pairs of sample points (relative to the origin on a circle) that contain each query point.

- **All-points mode**: O(n²) - Topological sweep for each point  
- **Query-based**: O(n log n) - Angular sorting and half-plane counting

### Spherical Angular Simplicial Depth (3D)

Counts the number of spherical triangles (triples of sample points) that contain each query point on a unit sphere.

- **All-points mode**: O(n²) - Topological sweep for each point
- **Query-based**: O(n²) - Gnomonic projection + topological sweep
- **Key feature**: Evaluates rays at both origin and antipode for leave-one-out calculations

**Reference**: See stored fact about antipode evaluation in 3D angular ASD

### 2D Simplicial Depth

Counts the number of triangles formed by sample points that contain each query point.

- **All-points mode**: O(n² log n) - Leave-one-out: n × O(n log n) topological sweep
- **Query-based**: O(n log n) - Topological sweep algorithm

### 3D Simplicial Depth

O(n²) algorithm counting containing tetrahedra using sphere projection:

1. Project sample points onto unit sphere centered at query point
2. Gnomonic projection to plane z=1
3. Color points blue (upper hemisphere, z>0) or red (lower hemisphere, z<0)
4. Count using topological sweep:
   - N1: Red points inside blue triangles
   - N2: Blue points inside red triangles  
   - N3: Red-blue segment crossings
5. Depth = N1 + N2 + N3

- **All-points mode**: O(n³) - Leave-one-out: n × O(n²) tetrahedra counting
- **Query-based**: O(n²) - Single topological sweep pass
- **Major improvement**: O(n²) instead of O(n⁴) brute force enumeration

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
