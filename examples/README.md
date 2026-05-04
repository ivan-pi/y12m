# Y12M Examples

This folder contains example driver programs that demonstrate how to use the
Y12M sparse linear-algebra library.  Drivers are grouped into three categories:

* **Synthetic matrix tests** – parametric matrices from the original Zlatev
  benchmark suite (classes D, E and F2), plus timing drivers.
* **Discretised PDE examples** – physically-motivated problems arising from
  finite-difference and finite-element discretisations of PDEs.
* **Legacy drivers** – the original fixed-form Fortran 77 programs distributed
  with the Y12M package.

A shared utility module `y12m_example_util.f90` provides the three matrix
generators (`matrd1`, `matre1`, `matrf2`) and a wall-clock timing helper
(`time`) used by several drivers.

---

## Contents

- [API usage table](#api-usage-table)
- [Synthetic matrix tests](#synthetic-matrix-tests)
  - [`timings_de.f90`](#timings_def90)
  - [`matf_bench.f90`](#matf_benchf90)
- [Discretised PDE examples](#discretised-pde-examples)
  - [`poisson_9pt.F90`](#poisson_9ptf90)
  - [`biharmonic_13pt.f90`](#biharmonic_13ptf90)
  - [`fem_anisotropic.f90`](#fem_anisotropicf90)
  - [`darcy_flow.f90`](#darcy_flowf90)
- [Legacy drivers](#legacy-drivers)

---

## API usage table

The table below lists every modern (free-form Fortran 90+) driver together
with the Y12M API routines it exercises.  Columns will be extended as
coverage grows.

| Driver | `Y12MA` |
|--------|:-------:|
| `timings_de` | ✓ |
| `matf_bench` | ✓ |
| `poisson_9pt` | ✓ |
| `biharmonic_13pt` | ✓ |
| `fem_anisotropic` | ✓ |
| `darcy_flow` | ✓ |

---

## Synthetic matrix tests

These drivers use the parametric sparse matrices introduced by Zlatev et al.
to benchmark the solver across a wide range of problem sizes and fill levels.
The right-hand side is constructed as the row sum of the matrix so that the
exact solution is the all-ones vector, giving a built-in correctness check.

> **Note:** Annotated sparsity-pattern images for classes D, E and F2 will be
> added in a future update.

### `timings_de.f90`

Timing driver for **class D** and **class E** matrices, reproducing the
benchmark tables from Zlatev et al. in double precision.

* **Class D(n, c):** square, unsymmetric; NNZ = 4·n + 55; requires n > 22,
  1 ≤ c ≤ n − 13.
* **Class E(n, c):** symmetric positive-definite five-point pattern;
  NNZ = 5·n − 2·c − 2; requires n ≥ 3, 2 ≤ c ≤ n − 1.

The driver has an object-oriented design: an abstract `coo_matrix` base type
provides the storage and the solver call (`y12ma`); concrete subtypes
`matrix_d` and `matrix_e` generate their respective sparsity patterns.

```
Usage: timings_de [--help] [-M {d|e|de}] [-n RANGE] [-c RANGE]
       RANGE: single integer or begin:end:step
```

Defaults reproduce the Zlatev benchmark table (`-M d -n 650:1000:50 -c 4:204:40`).

### `matf_bench.f90`

Timing driver for **class F2** matrices.  F2 generalises class D by adding a
lower-left 10×10 triangular block and a cyclic off-diagonal band of width r − 1
located at a column distance c from the main diagonal.

* **Class F2(n, c, r, α):** square; NNZ = r·n + 110; requires n ≥ 22,
  11 ≤ c ≤ n − 11, 2 ≤ r ≤ n − c − 9, α ≥ 1.  The parameter α controls
  the condition number.

```
Usage: matf_bench [--help] [-n RANGE] [-c RANGE] [-r RANGE] [--alpha VALUE]
       RANGE: single integer or begin:end:step
```

---

## Discretised PDE examples

These drivers build sparse linear systems arising from the discretisation of
partial differential equations on structured grids, solve them with `y12ma`,
and verify the solution against an analytically known exact solution using
the Method of Manufactured Solutions (MMS) or a Fourier-series reference.

Convergence shell scripts (`*_convergence.sh`) perform a grid-refinement
study and save the solution data for plotting with gnuplot (`plot_poisson.gp`).

### `poisson_9pt.F90`

Solves the **2-D Laplace equation** (steady-state heat diffusion) on the unit
square with a parabolic Dirichlet boundary condition on the top edge:

```
  Δu = 0   on (0,1)×(0,1)
  u(x,1) = 4x(1−x),   u = 0 on all other edges
```

The discretisation uses the **9-point isotropic Laplacian stencil** (weight
matrix `[1 4 1; 4 −20 4; 1 4 1] / (6h²)`).  The reference solution is a
Fourier sine series (odd modes only).  Expected convergence order: 4th order
in the mesh spacing h.

```
Usage: poisson_9pt [--help] [N] [output_file]
```

### `biharmonic_13pt.f90`

Solves the **2-D Biharmonic equation** (Kirchhoff plate bending) on the unit
square with simply-supported boundary conditions:

```
  ∇⁴u = f(x,y)   on (0,1)×(0,1),   u = 0,  ∇²u = 0  on ∂Ω
```

Uses a **13-point finite-difference stencil** for ∇⁴.  A command-line switch
also enables a 21-point alternative stencil derived from the convolution of the
5-point and 9-point Laplacians.  The right-hand side and exact solution are set
via a manufactured solution (two difficulty levels available).

```
Usage: biharmonic_13pt [--help] [N] [--mms-level {1|2}] [--stencil {13|21}]
```

### `fem_anisotropic.f90`

Solves a **2-D anisotropic diffusion problem** via the finite-element method:

```
  −∇·(K ∇u) = f(x,y)   on (0,1)×(0,1),   u = 0 on ∂Ω
```

where K is a full 2×2 diffusion tensor.  Features:

* **Bilinear Q1 quadrilateral elements** on a structured NX×NY grid.
* 2×2 Gauss quadrature for exact stiffness integration.
* 3×3 Gauss quadrature for highly accurate force integration.
* Independent NX and NY grid dimensions (rectangular elements).
* Exact solution: `u(x,y) = exp(x) sin(πy)`.

```
Usage: fem_anisotropic [--help] [NX] [NY]
```

### `darcy_flow.f90`

Solves a **mixed Darcy flow** problem near a stagnation point on the unit
square:

```
  q + ∇p = 0
  ∇·q = 0
```

using a **staggered Marker-and-Cell (MAC) grid** (equivalent to an RT0-P0
mixed finite-element discretisation).  The resulting indefinite saddle-point
system has velocity and pressure degrees of freedom interleaved cell-by-cell
to minimise LU fill-in.  The pressure nullspace is removed by pinning one cell
to its exact value.  The exact solution is a potential-flow stagnation-point
field.

Note: MAC/RT0-P0 is not the optimal solver architecture for this class of
mixed problem (a Schur-complement or block-preconditioned approach would be
more efficient), but this driver is used as a proxy application to test Y12M
on an indefinite saddle-point system with a challenging sparsity pattern.

```
Usage: darcy_flow [--help] [N]
```

---

## Legacy drivers

The files below are the original **fixed-form Fortran 77** drivers distributed
with the Y12M package.  They read their parameters from standard input (using
`READ`) and write results to standard output.  They are retained for
historical reference and regression testing; the modern replacements
(`timings_de` and `matf_bench`) supersede them.

| File | Matrix class | API routine | Description |
|------|-------------|-------------|-------------|
| `maind.f` | D | `Y12MAE` | Timing/accuracy test for class D and/or E matrices; reads 7 integers from stdin (`indexp nstart nincr nend cstart cincr cend`). |
| `maine.f` | D / E | `Y12MFE` | Variant of `maind.f` with iterative-refinement solver (`y12mfe`); reads the same 7 integers plus AFLAG/IFLAG control parameters. |
| `mainf.f` | F2 | `Y12MBE` / `Y12MCE` / `Y12MDE` / `Y12MFE` | Full test for class F2 matrices; exercises both direct (`y12mbe`+`y12mce`+`y12mde`) and iterative-refinement (`y12mfe`) paths; reads matrix parameters and AFLAG/IFLAG from stdin. |
| `mainf2.f` | F2 | `Y12MBE` / `Y12MCE` / `Y12MDE` | Debugging variant of `mainf.f` with extra `WRITE` trace statements; useful for step-by-step inspection of the solver. |
| `test1.f` | small dense | `Y12MAE` | Minimal worked example from the Y12M manual; solves a small system read entirely from stdin (n, nnz, then COO entries, then RHS). |
| `test2.f` | small dense | `Y12MFE` | Same as `test1.f` but calls the iterative-refinement driver `y12mfe`; compares solution in `x` with the all-ones vector. |

The data file `maind.d` contains a sample input for `maind` (class D,
n = 600, five values of c).
