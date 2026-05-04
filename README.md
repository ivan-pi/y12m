# y12m

Sparse direct solver for systems of linear equations **Ax = b**.  
Sequential, in-core LU factorization with partial pivoting and Markowitz
reordering.

> **Disclaimer:** This is a modified version of the Y12M library originally
> distributed through Netlib. The original source and documentation are
> available at <http://www.netlib.org/y12m/>.  All algorithmic credit belongs
> to the original authors.

**Book:**

> Zlatev, Z., Wasniewski, J., & Schaumburg, K. (1981). *Y12M: Solution of
> Large and Sparse Systems of Linear Algebraic Equations*, Lecture Notes in
> Computer Science, Vol. 121. Springer, Berlin.
> <https://doi.org/10.1007/3-540-10874-2>

**Home page of original author Zahari Zlatev:** <https://www.dmu.dk/atmosphericenvironment/staff/zlatev.htm>

---

## Contents

- [When to use Y12M](#when-to-use-y12m)
  - [Advantages](#advantages)
  - [Limitations](#limitations)
  - [Alternatives worth considering](#alternatives-worth-considering)
- [Building with CMake](#building-with-cmake)
- [API overview](#api-overview)
- [Citation](#citation)

---

## When to use Y12M

### Advantages

- **Simple and self-contained.** Pure Fortran, no external dependencies, easy
  to embed in existing projects.
- **Robust pivoting.** Combines Markowitz minimum-degree reordering with
  threshold pivoting; reliable on the difficult matrices that arise in stiff
  ODE integration and atmospheric chemistry.
- **Developed for production use.** Designed in the early 1980s by the
  original authors specifically for large-scale stiff ODE integrators and
  atmospheric chemistry models.
- **Multiple-RHS and structural reuse.** The lower-level API (`Y12MB` /
  `Y12MC` / `Y12MD`) lets you reuse an LU factorization or a sparsity
  ordering across multiple solves (see [docs/multiple_rhs.md](docs/multiple_rhs.md)).
- **Iterative refinement.** `Y12MF` provides built-in iterative refinement for
  extra accuracy.

### Limitations

- **32-bit integer indexing.** Internal integer arrays use default `INTEGER`
  (32-bit), so the practical limit is matrices of a few tens of thousands of
  rows/columns.  Very large problems (n ≳ 46 000) may trigger integer overflow
  in some internal Markowitz counters.
- **Sequential only.** There is no shared-memory or distributed-memory
  parallelism; a single thread is used throughout.
- **Fortran 77–style API.** The calling convention requires pre-allocated
  workspace arrays and integer flag vectors; it is more verbose than modern
  solver interfaces.
- **No fill-in reduction ordering.** `Y12MB` only builds the internal data
  structures (compressed row/column lists and the Markowitz linked list); it
  does **not** compute any global, pre-factorization fill-in reducing ordering
  such as AMD, COLAMD, or RCM.  During factorization `Y12MC` applies a greedy
  *local* Markowitz strategy that selects, at each step, the pivot minimising
  the product `(nnz_in_row − 1) × (nnz_in_col − 1)`, but this is no substitute
  for a global graph-based ordering.

  **Workaround — permute yourself.**  Because Y12M accepts its input in
  coordinate (COO) format (`SNR`/`RNR` index arrays), it is straightforward to
  apply an external permutation before calling `Y12MB`:

  1. Compute a fill-reducing permutation `perm(1:n)` with an external library
     (e.g. [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview),
     [AMD](https://people.engr.tamu.edu/davis/suitesparse.html), or the
     Reverse Cuthill–McKee routine in
     [Sparsekit](https://people.sc.fsu.edu/~jburkardt/f_src/sparsekit/sparsekit.html)).
  2. Renumber the stored indices:
     `SNR(i) = perm(SNR(i))` and `RNR(i) = perm(RNR(i))` for `i = 1 … Z`.
  3. Permute the right-hand side before the solve: `b_perm(perm(i)) = b(i)`.
  4. After `Y12MD` returns, recover the original ordering:
     `x(i) = x_perm(perm(i))`.

### Alternatives worth considering

For a broader survey see:
[*Direct Solvers for Sparse Matrices*](https://portal.nersc.gov/project/sparse/superlu/SparseDirectSurvey.pdf),
X. Li, April 2022.

| Library | Notes |
|---------|-------|
| [UMFPACK](https://people.engr.tamu.edu/davis/suitesparse.html) | Multifrontal LU; part of SuiteSparse. Widely used general-purpose choice. |
| [SuperLU](https://portal.nersc.gov/project/sparse/superlu/) | Supernodal LU; sequential and distributed (SuperLU_DIST) variants available. |
| [Eigen `SparseLU`](https://eigen.tuxfamily.org) | Header-only C++ with a modern API; good for moderate-size problems. |
| [MUMPS](https://mumps-solver.org) | Multifrontal, supports MPI; excellent for very large distributed problems. |
| [MKL PARDISO / DSS](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/current/pardiso.html) | High-performance threaded solver bundled with Intel oneAPI MKL. |
| [Panua PARDISO](https://panua.ch/pardiso/) | Commercial successor to the original Basel PARDISO; supports shared-memory parallelism. |
| [Apple Accelerate Sparse](https://developer.apple.com/documentation/accelerate/sparse_solvers) | Native sparse direct solvers on macOS/iOS via the Accelerate framework. |

---

## Building with CMake

```sh
git clone https://github.com/ivan-pi/y12m.git && cd y12m
cmake -B build -DCMAKE_Fortran_COMPILER=gfortran
cmake --build build
ctest --test-dir build          # run the test suite
```

The build produces a static library `liby12m_legacy` and, optionally, the
example programs.  The `y12m` Fortran module (`.mod` file) is placed in
`build/include/`.

---

## API overview

The `y12m` Fortran module exposes **generic interfaces** (`y12ma`, `y12mb`, …)
that dispatch automatically to the single-precision (`E`) or double-precision
(`F`) variant based on the type of the actual arguments.  The
precision-specific external procedures (`y12mbe` / `y12mbf`, etc.) can also be
called directly without the module.

### High-level drivers

| Subroutine | Purpose |
|------------|---------|
| `Y12MA` | Black-box driver: single matrix, single RHS. Calls `Y12MB` + `Y12MC` + `Y12MD` internally. |
| `Y12MF` | Factorize and solve with iterative refinement in one call. |

### Lower-level subroutines

| Subroutine | Purpose |
|------------|---------|
| `Y12MB` | Reorder and prepare the matrix for factorization. |
| `Y12MC` | Compute the LU factorization. |
| `Y12MD` | Solve using an existing LU factorization. |
| `Y12MG` | Estimate the reciprocal condition number. *(optional)* |
| `Y12MH` | Compute the one-norm of **A**. *(optional, must precede `Y12MC`)* |

### Calling order

<p align="center">
   <img src="callseq.png" width="30%" alt="Y12M calling sequence diagram">
</p>

For detailed usage scenarios (multiple RHS, structure reuse, condition
estimation) see [docs/multiple_rhs.md](docs/multiple_rhs.md) and the full
[API reference](docs/API.md).

---

## Citation

If you use Y12M in published work, please cite the accompanying book:

```bibtex
@book{zlatev1981y12m,
  author    = {Zlatev, Zahari and Wasniewski, Jerzy and Schaumburg, Kjeld},
  title     = {{Y12M}: Solution of Large and Sparse Systems of Linear
               Algebraic Equations},
  series    = {Lecture Notes in Computer Science},
  volume    = {121},
  publisher = {Springer},
  address   = {Berlin, Heidelberg, New York},
  year      = {1981},
  doi       = {10.1007/3-540-10874-2}
}
```

