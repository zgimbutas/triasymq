# CLAUDE.md

Fortran 77 library of precomputed high-order quadrature and interpolation
rules on triangles, squares, tetrahedra and cubes (Xiao and Gimbutas).
README.md is the authoritative reference: per-module descriptions,
makefile invocations, calling sequences, available orders and the
`tables/` data.

## Build and run

- One `makefile.<module>` per module (`makefile.triasymq`, `makefile.squarearbq`,
  ...); bare `make` runs `makefile.triasymq`. Each builds `int2` and runs it.
- `HOST` selects the compiler flags from `make.inc`: `linux-gfortran`,
  `linux-gfortran-openmp`, `macos-gfortran`, `macos-gfortran-openmp`. Every
  `makefile.*` sets `HOST=linux-gfortran-openmp` itself, so override it on
  the make command line, not in the environment:
  `make -f makefile.triasymq HOST=macos-gfortran`.
- Driver output lands in `fort.*` files (`prini(6,13)` → stdout and
  `fort.13`; `triasymq_dr` also writes nodes and gnuplot files to
  `fort.11`, `fort.12`, `fort.14`). `make distclean` removes them.

## Fortran conventions

- Compiled with `-std=legacy`: fixed-form source, 72-column lines,
  `implicit real *8`.
- All tables were computed in quad precision.

## Reference

H. Xiao, Z. Gimbutas, "A numerical algorithm for the construction of
efficient quadrature rules in two and higher dimensions," Computers and
Mathematics with Applications, 59 (2009), pp. 663-676.
