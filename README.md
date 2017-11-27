# Quadratures for triangles, squares, cubes and tetrahedra

Copyright (C) 2009-2012: Zydrunas Gimbutas and Hong Xiao

Contact: Zydrunas Gimbutas <gimbutas@cims.nyu.edu>
         Hong Xiao <hxiao@ucdavis.edu>

Date: November 26, 2017

Version 1.8


### Contents

```
triasymq.f - fully symmetric (D_3) quadrature for triangle
triarotq.f - rotationally symmetric quadrature for triangle
triaarbq.f - arbitrary symmetric (or asymmetric) quadrature for triangle
triaintq.f - interpolation nodes/quadrature for smooth functions on triangle
```
```
squaresymq.f - fully symmetric (D_4) quadrature for a unit square
squaresymvq.f - centro/rotationally symmetric quadrature for a unit square
squarearbq.f - arbitrary symmetric (or asymmetric) quadrature for a unit square
```

```
tetraarbq.f - arbitrary symmetric (or asymmetric) quadrature for tetrahedron
tetraintq.f - interpolation nodes/quadrature for smooth functions on tetrahedron
```
```
cubearbq.f - arbitrary symmetric (or asymmetric) quadrature for a unit cube
```


### Quadratures for triangles


```
makefile - run a simple test (triasymq)
```

```
triasymq.f - construct (or rather retrieve) a fully symmetric (D_3) 
quadrature formula on the user-defined triangle in the plane.

triasymq_dr.f - driver routines for triasymq.f

makefile.triasymq - makefile for triasymq.f
                       ( make -f makefile.triasymq )
```

```
c
c   n       1     2     3     4     5     6     7     8     9    10
c  -----------------------------------------------------------------
c nodes     1     3     6     6     7    12    15    16    19    25
c
c
c   n      11    12    13    14    15    16    17    18    19    20
c  -----------------------------------------------------------------
c nodes    28    33    37    42    49    55    60    67    73    79
c
c
c   n      21    22    23    24    25    26    27    28    29    30
c  -----------------------------------------------------------------
c nodes    87    96   103   112   120   130   141   150   159   171
c
c
c   n      31    32    33    34    35    36    37    38    39    40
c  -----------------------------------------------------------------
c nodes   181   193   204   214   228   243   252   267   282   295
c
c
c   n      41    42    43    44    45    46    47    48    49    50
c  -----------------------------------------------------------------
c nodes   309   324   339   354   370   385   399   423   435   453
```

```
triarotq.f - construct (or rather retrieve) a rotationally symmetric
quadrature formula on the user-defined triangle in the plane.

triarotq_dr.f - driver routines for triarotq.f

makefile.triarotq - makefile for triarotq.f
                       ( make -f makefile.triarotq )
```

```
c   n       1     2     3     4     5     6     7     8     9    10
c  -----------------------------------------------------------------
c nodes     1     3     6     6     7    12    12    16    19    24
c
c
c   n      11    12    13    14    15    16    17    18    19    20
c  -----------------------------------------------------------------
c nodes    27    33    36    42    48    54    60    66    70    78
c
c
c   n      21    22    23    24    25    26    27    28    29    30
c  -----------------------------------------------------------------
c nodes    85    93   102   111   117   129   138   147   156   168
c
c
c   n      31    32    33    34    35    36    37    38    39    40
c  -----------------------------------------------------------------
c nodes   177   189   201   213   225   237   249   262   277   291
c
c
c   n      41    42    43    44    45    46    47    48    49    50
c  -----------------------------------------------------------------
c nodes   303   318   333   348   363   378   393   414   429   444
```

```
triaarbq.f - construct (or rather retrieve) an asymmetric
quadrature formula on the user-defined triangle in the plane.

triaarbq_dr.f - driver routines for triaarbq.f

makefile.triaarbq - makefile for triaarbq.f
                       ( make -f makefile.triaarbq )
```

```
c   n       1     2     3     4     5     6     7     8     9    10
c  -----------------------------------------------------------------
c nodes     1     3     4     6     7    11    12    16    19    24
c
c
c   n      11    12    13    14    15    16    17    18    19    20
c  -----------------------------------------------------------------
c nodes    27    32    36    41    46    53    58    65    71    78
c
c
c   n      21    22    23    24    25    26    27    28    29    30
c  -----------------------------------------------------------------
c nodes    87    94   101   109   118   128   136   148   156   170
c
c
c   n      31    32    33    34    35    36    37    38    39    40
c  -----------------------------------------------------------------
c nodes   178   192   204   212   225   236   249   265   274   290
c
c
c   n      41    42    43    44    45    46    47    48
c  -----------------------------------------------------------------
c nodes   302   318   332   348   363   380   393   412 
```

```
triaintq.f - construct (or rather retrieve) fully symmetric (D_3) 
interpolation nodes on the user-defined triangle in the plane.

triaintq_dr.f - driver routines for triaintq.f

makefile.triaintq - makefile for triaintq.f
                       ( make -f makefile.triaintq )
```

```
c interp    0     1     2     3     4     5     6     7     8     9    10
c  -----------------------------------------------------------------------
c quadr     1     2     4     5     7     8    10    12    14    15    17
c  -----------------------------------------------------------------------
c nodes     1     3     6    10    15    21    28    36    45    55    66
c  -----------------------------------------------------------------------
c cond #   1.0   1.0   1.4   1.9   2.1   3.4   4.3   4.8   4.8   6.5   8.1  
c
c
c interp   11    12    13    14    15    16    17    18    19    20
c  -----------------------------------------------------------------------
c quadr    19    20    22    24    25    27    28    30    32    33
c  -----------------------------------------------------------------------
c nodes    78    91   105   120   136   153   171   190   210   231
c  -----------------------------------------------------------------------
c cond #  15.7  19.2  21.4  38.6  31.4  44.3  75.3  117   153   194
```

### Quadratures for squares

```
squaresymq.f - construct (or rather retrieve) a fully symmetric (D_4)
quadrature formula for the unit square [-1,1]^2.

squaresymq_dr.f - driver routines for squaresymq.f

makefile.squaresymq - makefile for squaresymq.f 
                       ( make -f makefile.squaresymq )
```

```
c   n       1     2     3     4     5     6     7     8     9    10
c  -----------------------------------------------------------------
c nodes     1     4     4     8     8    12    12    20    20    28
c
c
c   n      11    12    13    14    15    16    17    18    19    20
c  -----------------------------------------------------------------
c nodes    28    37    37    48    48    57    57    72    72    85
c
c
c   n      21
c  -----------------------------------------------------------------
c nodes    85
```

```
squaresymvq.f - construct (or rather retrieve) a variably symmetric
quadrature formula for the unit square [-1,1]^2.

squaresymvq_dr.f - driver routines for squaresymvq.f

makefile.squaresymvq - makefile for squaresymvq.f 
                       ( make -f makefile.squaresymvq )
```

```
c   n       1     2     3     4     5     6     7     8     9    10
c  -----------------------------------------------------------------
c nodes     1f    4f    4f    7c    7c   12f   12f   17r   17r*  24r
c
c
c   n      11    12    13    14    15    16    17    18    19    20
c  -----------------------------------------------------------------
c nodes    24r*  33r   33r*  43c   43c*  54c   54c*  67c   67c*  81r
c
c
c   n      21    22    23    24    25    26    27    28    29   
c  -----------------------------------------------------------------
c nodes    81r*  96r   96r* 115c  115c* 132r  132r* 152r  152r* 
c
c
c       f  - fully symmetric (D_4 symmetry)
c       r  - rotational symmetry only (rotate by 90 degrees)
c       c  - center symmetry only (rotate by 180 degrees)
c
c       *  - number of quadrature nodes is less than (n+1)*(n+2)/2/3
```

```
squarearbq.f - construct (or rather retrieve) a quadrature formula 
for the unit square [-1,1]^2.

squarearbq_dr.f - driver routines for squarearbq.f

makefile.squarearbq - makefile for squarearbq.f ( make -f makefile.squarearbq )
```

```
c   n       1     2     3     4     5     6     7     8     9    10
c  -----------------------------------------------------------------
c nodes     1f    3    4f     6    7c    10    12f   16    17r*  22
c
c
c   n      11    12    13    14    15    16    17    18    19    20
c  -----------------------------------------------------------------
c nodes    24r*  31    33r*  40    43c*  52   54c*   64   67c*   78
c
c
c   n      21    22    23    24    25    26    27    28    29    30
c  -----------------------------------------------------------------
c nodes    81r*  93    96r* 109   115c* 127   132r* 146   152r* 167
c
c
c       f  - fully symmetric (D_4 symmetry)
c       r  - rotational symmetry only (rotate by 90 degrees)
c       c  - center symmetry only (rotate by 180 degrees)
c
c       *  - number of quadrature nodes is less than (n+1)*(n+2)/2/3
```

### Quadratures for tetrahedra

```
tetraarbq.f - construct (or rather retrieve) a quadrature formula for
the standard tetrahedron, defined by vertices
(-1,-1/Sqrt[3],-1/Sqrt[6]), (0,2/Sqrt[3],-1/Sqrt[6]), 
(1,-1/Sqrt[3],-1/Sqrt[6]), and (0,0,3/Sqrt[6]).
Only orders 1 and 5 are fully symmetric.

tetraarbq_dr.f - driver routines for tetraarbq.f

makefile.tetraarbq - makefile for tetraarbq.f ( make -f makefile.tetraarbq )
```

```
c   n       1     2     3     4     5     6     7     8     9    10
c  -----------------------------------------------------------------
c nodes     1f    4     6    11    14f   23    31    44    56    74
c
c   n      11    12    13    14    15
c  -----------------------------------------------------------------
c nodes    95   122   146   177   214
c
c
c       f  - fully symmetric
```

```
tetraintq.f - construct (or rather retrieve) fully symmetric
interpolation nodes for
the standard tetrahedron, defined by vertices
(-1,-1/Sqrt[3],-1/Sqrt[6]), (0,2/Sqrt[3],-1/Sqrt[6]), 
(1,-1/Sqrt[3],-1/Sqrt[6]), and (0,0,3/Sqrt[6]).

tetraintq_dr.f - driver routines for tetraintq.f

makefile.tetraintq - makefile for tetraintq.f ( make -f makefile.tetraintq )
```

```
c interp    0     1     2     3     4     5     6     7     8     9    10
c  -----------------------------------------------------------------------
c quadr     1     2     3     5     6     7     9    10    11    13    15
c  -----------------------------------------------------------------------
c nodes     1     4    10    20    35    56    84   120   165   220   286
c  -----------------------------------------------------------------------
c cond #   1.0   1.0   1.6   1.9   3.4   3.5   4.6  10.7  12.7  16.8  37.0
```

### Quadratures for cubes

```
cubearbq.f - construct (or rather retrieve) a quadrature formula for
the unit cube [-1,1]^3. These are non-symmetric quadrature nodes in
general, but come in pairs for orders 3, 5, 7, 9, 11, 13, and 15. 

cubearbq_dr.f - driver routines for cubearbq.f

makefile.cubearbq - makefile for cubearbq.f ( make -f makefile.cubearbq )

```

```
c   n       1     2     3     4     5     6     7     8     9    10
c  -----------------------------------------------------------------
c nodes     1c    4     6c   10    13c*  22    26c*  42    48c*  73
c
c
c   n      11    12    13    14    15   
c  -----------------------------------------------------------------
c nodes    82c* 115   129c* 173   188c*
c
c       
c       c  - centro-symmetric
c       *  - number of quadrature nodes is less than (n+1)*(n+2)*(n+3)/6/4
```


### Testing and debugging routines

```
lege2eva.f - orthogonal polynomials on a unit square
lege3eva.f - orthogonal polynomials on a unit cube
ortho2eva.f - orthogonal polynomials on a triangle
ortho3eva.f - orthogonal polynomials on a tetrahedron
```

```
tetragauc.f - tensor product rule for smooth functions on a tetrahedron
gaussq.f - Gaussian type quadratures (from netlib.org)
```

```
legeexps.f - Legendre expansion routines
prini.f - printing routines
```

```
gen_table - helper script to generate *_table.txt files
```

```
tables - raw data files for import to matlab, octave, and fortran
```

### Notes

1. All tables have bee recomputed in quad precision, moderate
improvements have been made for square*q and cube*q routines.

2. triaarbq quadratures are not-optimized, see triarotq.

3. In this release, 14th order quadrature in cubearbq table has 173
nodes (in the reference paper, the corresponding quadrature has 172
nodes but one negative weight).

4. The library is released under a modified 3-clause BSD license.


### References

H. Xiao, Z. Gimbutas, "A numerical algorithm for the construction
of efficient quadrature rules in two and higher dimensions,"
Computers and Mathematics with Applications, 59 (2009), pp. 663-676







 
