# CLAUDE.md

This file provides guidance for Claude Code when working with the triasymq repository.

## Project Overview

**triasymq** is a Fortran 77 library for high-precision numerical quadrature (integration) rules on geometric shapes: triangles, squares, tetrahedra, and cubes. It retrieves pre-computed optimal quadrature nodes and weights for numerical integration of smooth functions.

- **Authors:** Zydrunas Gimbutas (NYU), Hong Xiao (UC Davis)
- **Version:** 1.11
- **License:** Modified 3-clause BSD (see COPYING)

## Build Commands

```bash
# Build default (triasymq)
make

# Build specific modules
make -f makefile.triasymq      # Triangle D_3 symmetric quadrature
make -f makefile.triarotq      # Triangle rotational quadrature
make -f makefile.triaarbq      # Triangle arbitrary quadrature
make -f makefile.triaintq      # Triangle interpolation nodes
make -f makefile.squaresymq    # Square D_4 symmetric quadrature
make -f makefile.squaresymvq   # Square variable symmetric quadrature
make -f makefile.squarearbq    # Square arbitrary quadrature
make -f makefile.squareintq    # Square interpolation nodes
make -f makefile.tetraarbq     # Tetrahedron quadrature
make -f makefile.tetraintq     # Tetrahedron interpolation nodes
make -f makefile.cubearbq      # Cube quadrature

# Clean build artifacts
make clean        # Remove .o files
make distclean    # Remove .o files and fort.* output files
```

Set the HOST environment variable before building:
- `linux-gfortran` - Linux with gfortran
- `linux-gfortran-openmp` - Linux with OpenMP support
- `macos-gfortran` - macOS with gfortran
- `macos-gfortran-openmp` - macOS with OpenMP support

Example: `HOST=linux-gfortran make -f makefile.triasymq`

## Code Architecture

### File Naming Conventions

| Pattern | Purpose |
|---------|---------|
| `*.f` | Core Fortran implementation |
| `*_dr.f` | Driver/test programs |
| `*_table.txt` | Pre-computed quadrature data tables |
| `makefile.*` | Module-specific makefiles |

### Main Modules

| Module | Symmetry | Shape | Orders |
|--------|----------|-------|--------|
| `triasymq.f` | D_3 (full) | Triangle | 1-50 |
| `triarotq.f` | C_3 (rotational) | Triangle | 1-50 |
| `triaarbq.f` | Arbitrary | Triangle | 1-50 |
| `triaintq.f` | D_3 (interp) | Triangle | 0-20 |
| `squaresymq.f` | D_4 (full) | Square | 1-21 |
| `squaresymvq.f` | Variable | Square | 1-29 |
| `squarearbq.f` | Arbitrary | Square | 1-30 |
| `squareintq.f` | Interp | Square | 0-20 |
| `tetraarbq.f` | Arbitrary | Tetrahedron | 1-15 |
| `tetraintq.f` | Interp | Tetrahedron | 0-10 |
| `cubearbq.f` | Arbitrary | Cube | 1-15 |

### Supporting Files

- `ortho2eva.f` - Orthogonal polynomials on triangles
- `ortho3eva.f` - Orthogonal polynomials on tetrahedra
- `lege2eva.f` - Orthogonal polynomials on unit square
- `lege3eva.f` - Orthogonal polynomials on unit cube
- `legeexps.f` - Legendre expansion routines
- `prini.f` - Printing/formatting utilities
- `tetragauc.f` - Tensor product Gaussian quadrature for tetrahedra
- `gaussq.f` - Gaussian quadratures (from netlib.org)

### Data Directory

`tables/` contains raw quadrature data for import into MATLAB/Octave/Fortran.

## Fortran Conventions

- **Standard:** Fortran 77 (legacy) - compiled with `-std=legacy` flag
- **Fixed-form source:** 6-character indentation, 72-character line limit
- **Implicit typing:** Uses `implicit real *8` convention
- **Precision:** All tables computed in quad precision for accuracy
- **Comments:** Lines starting with `c` or `C` in column 1

## Running Tests

Each module has a driver program (`*_dr.f`) that tests the quadrature:

```bash
# Build and run triasymq test
HOST=linux-gfortran make -f makefile.triasymq
./int2

# Output appears in fort.* files (Fortran unit numbers)
```

## Key Subroutines

Most modules follow a similar API pattern:

```fortran
c     Retrieve quadrature for order n on triangle with vertices (vert)
      call triasymq(n, vert1, vert2, vert3, rnodes, weights, npts)
```

- `n` - Quadrature order (polynomial degree of exactness)
- `vert*` - Vertex coordinates defining the shape
- `rnodes` - Output: quadrature node coordinates
- `weights` - Output: quadrature weights
- `npts` - Output: number of quadrature points

## Reference

H. Xiao, Z. Gimbutas, "A numerical algorithm for the construction of efficient quadrature rules in two and higher dimensions," Computers and Mathematics with Applications, 59 (2009), pp. 663-676
