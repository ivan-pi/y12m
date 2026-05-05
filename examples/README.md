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
  - [`poisson_2d.F90`](#poisson_2df90)
  - [`poisson_3d.f90`](#poisson_3df90)
  - [`biharmonic_13pt.f90`](#biharmonic_13ptf90)
  - [`fem_anisotropic.f90`](#fem_anisotropicf90)
  - [`darcy_flow.f90`](#darcy_flowf90)
  - [`heat_implicit.f90`](#heat_implicitf90)
  - [`multiple_loads.f90`](#multiple_loadsf90)
  - [`newton_bratu.f90`](#newton_bratuf90)
- [Legacy drivers](#legacy-drivers)

---

## API usage table

The table below lists every modern (free-form Fortran 90+) driver together
with the Y12M API routines it exercises.  Columns will be extended as
coverage grows.

| Driver | `Y12MA` | `Y12MB` | `Y12MC` | `Y12MD` |
|--------|:-------:|:-------:|:-------:|:-------:|
| `timings_de` | ✓ | | | |
| `matf_bench` | ✓ | | | |
| `poisson_2d` | ✓ | | | |
| `poisson_3d` | ✓ | | | |
| `biharmonic_13pt` | ✓ | | | |
| `fem_anisotropic` | ✓ | | | |
| `darcy_flow` | ✓ | | | |
| `heat_implicit` | | ✓ | ✓ | ✓ |
| `multiple_loads` | | ✓ | ✓ | ✓ |
| `newton_bratu` | | ✓ | ✓ | ✓ |

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

### `poisson_2d.F90`

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
Usage: poisson_2d [--help] [N] [output_file]
```

### `poisson_3d.f90`

Solves the **3-D Poisson equation** on the unit cube with homogeneous Dirichlet
boundary conditions using a **7-point finite-difference stencil** and a
manufactured sine-mode solution.  Expected convergence order: 2nd order in h.

```
Usage: poisson_3d [--help] [N] [kx] [ky] [kz]
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

### `heat_implicit.f90`

Solves the **2-D heat equation** (parabolic PDE) on the unit square using
implicit (backward Euler) time integration:

```
  u_t = κ (u_xx + u_yy)   on (0,1)×(0,1),  t ∈ (0, T]
  u = 0 on the boundary,   u(x,y,0) = sin(πx) sin(πy)
```

The exact solution (separation of variables) is
`u(x,y,t) = exp(−2π²κt) sin(πx) sin(πy)`.

At each time step the implicit Euler scheme requires solving the linear
system `A u^{n+1} = u^n` where `A = I + r K` with `r = dt κ / h²` and
K the 5-point stiffness matrix.  Because A is **constant across all time
steps**, the driver demonstrates **Case (ii) – same matrix, many
right-hand sides**:

* `y12mb` + `y12mc` (with `IFLAG(5) = 2`) are called **once** to factorize A.
* `y12md` is called at **every time step** using the stored LU factors,
  switching from `IFLAG(5) = 2` (first step) to `IFLAG(5) = 3`
  (all subsequent steps).

This amortises the factorization cost over all time steps, which is
particularly advantageous for long integration runs.

```
Usage: heat_implicit [--help] [N] [nsteps] [T] [kappa] [output_file]
       N           total grid size (>= 3, default 42)
       nsteps      number of time steps (default 200)
       T           final time (default 0.1)
       kappa       thermal diffusivity > 0 (default 1.0)
       output_file output file (default: heat_implicit.dat)
```

### `multiple_loads.f90`

Solves a family of **2-D Poisson problems** on the unit square with
homogeneous Dirichlet boundary conditions and a set of modal load cases:

```
  −∇²u_k = f_k(x,y) = sin(kπx) sin(πy)   on (0,1)×(0,1),  k = 1, ..., nrhs
```

Exact solutions: `u_k(x,y) = sin(kπx) sin(πy) / ((k²+1) π²)`.

The stiffness matrix K (5-point Laplacian) is **identical for all load
cases**.  The driver demonstrates **Case (ii) – same matrix, many
right-hand sides**:

* `y12mb` + `y12mc` (with `IFLAG(5) = 2`) factorize K **once**.
* `y12md` is called once per load case (`IFLAG(5) = 2` for the first
  RHS pre-substituted by `y12mc`, then `IFLAG(5) = 3` for the rest).

A naive implementation would refactorize K for each RHS at cost
`nrhs × O(nnz_LU)`; the LU-reuse strategy reduces this to
`1 × O(nnz_LU)` plus `nrhs` cheap substitution passes.

```
Usage: multiple_loads [--help] [N] [nrhs] [output_file]
       N           total grid size (>= 3, default 42)
       nrhs        number of load cases (>= 1, default 8)
       output_file output file (default: multiple_loads.dat)
```

### `newton_bratu.f90`

Solves the **nonlinear Bratu (Liouville–Gelfand) problem** on the unit
square using Newton's method:

```
  −∇²u = λ exp(u)   on (0,1)×(0,1),   u = 0 on ∂Ω
```

The problem has two solutions for λ < λ_c ≈ 6.808 and no solution for
λ > λ_c.  The Newton linearisation at step k reads:

```
  J(u^k) δu = −F(u^k)
  u^{k+1}   = u^k + δu
```

where the Jacobian J = K − h² λ diag(exp(u)) has the **same non-zero
pattern** as the 5-point Laplacian K at every iteration (only diagonal
values change).  The driver demonstrates **Case (iii) – same sparsity,
changing values (structural reuse)**:

* First Newton step: `IFLAG(4) = 1`.  `y12mb` computes the Markowitz
  column ordering and stores it in HA; `y12mc` + `y12md` complete the
  first solve.
* Subsequent steps: `IFLAG(4) = 2`.  `y12mb` reuses the stored ordering
  and skips the expensive Markowitz search, only inserting updated
  numerical values at the already-known positions.

`IFLAG` is initialised **once outside the Newton loop** so that the
fill-in counts stored in `IFLAG(9)`/`IFLAG(10)` by `y12mc` on the first
call are preserved for subsequent iterations.

```
Usage: newton_bratu [--help] [N] [lambda] [max_iter] [output_file]
       N           total grid size (>= 3, default 22)
       lambda      reaction parameter in (0, 6.808) (default 1.0)
       max_iter    max Newton iterations (default 20)
       output_file output file (default: newton_bratu.dat)
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
