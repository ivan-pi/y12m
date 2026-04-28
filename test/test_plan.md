# Y12M Test Plan

This document analyses the current test coverage of the y12m sparse-solver
package and lists tests worth adding in the future.

---

## 1. Summary of current tests

| Test file | Routine(s) exercised | Precision | Matrix types |
|---|---|---|---|
| `test_y12ma_sp.f90` | y12ma | SP | diagonal, tridiagonal, arrow |
| `test_y12ma_dp.f90` | y12ma | DP | diagonal, tridiagonal, arrow |
| `test_y12ma_5x5.F90` | y12ma | SP + DP | UMFPACK 5×5 |
| `test_y12ma_rand.F90` | y12ma | SP + DP | random 9×9 (LAPACK ref, optional) |
| `test_y12mb_mc_md_sp.f90` | y12mb + y12mc + y12md | SP | diagonal, full 2×2, tridiagonal, arrow |
| `test_y12mb_mc_md_dp.f90` | y12mb + y12mc + y12md | DP | same |
| `test_y12mb_errors.F90` | y12mb | SP + DP | error-path coverage |
| `test_y12mc_errors.F90` | y12mc | SP + DP | error-path coverage |
| `test_y12mc_z_intent.F90` | y12mc | SP + DP | regression: z intent violation |
| `test_y12md_errors.F90` | y12md | SP + DP | error-path coverage |
| `test_y12mf_sp.f90` | y12mf | SP | tridiagonal n=3,5,7,10 |
| `test_y12mf_errors.f90` | y12mf | SP | 6×6 reference; error/output invariants; LU-reuse |
| `test_y12mg_mh.f90` | y12mg + y12mh | SP + DP | tridiagonal |
| `y12m_solve.f90` | y12ma | SP + DP | Matrix Market files mat0–mat5 |
| `y12m_mm.f90` | y12ma | SP + DP | arbitrary Matrix Market (optional) |
| `test_y12ma_y12mf_bench_sp.f90` | y12ma + y12mf | SP | UMFPACK 5×5, HB 4×4, Pardiso 8×8, Templates 6×6 |

---

## 2. Coverage gaps

### 2.1 Missing precision variant for y12mf

`y12mfe` (single precision) is the only iterative-refinement routine in
`src/legacy/`.  No double-precision `y12mff` exists in the current source
tree, so there is no DP variant to test.  If a double-precision IR routine
is added in the future, tests analogous to `test_y12mf_sp.f90` should be
created.

### 2.2 C and C++ interfaces not tested

`src/y12m_c.f90` wraps the Fortran routines with `bind(C)` linkage, and the
public API is exposed through `include/y12m.h` (C) and `include/y12m.hpp`
(C++ templates).  No tests call these interfaces at all.

### 2.3 Large-n Markowitz overflow bug

The 32-bit overflow fix in `y12mce`/`y12mcf` (variable `nr = n*n` for large
`n` replaced by `nr = huge(nr)`) is not covered by any current test.

### 2.4 Ill-conditioned and near-singular matrices

All current test matrices are well-conditioned.  There are no tests that
verify error flags and iterative-refinement benefit on ill-conditioned
problems, nor tests that exercise the near-singularity detection path
(IFAIL = 7 or 8).

### 2.5 LU structure reuse (multiple same-structure systems)

The `y12mb + y12mc + y12md` reuse path (`IFLAG(4)=2`, `IFLAG(5)=3`) is
tested only inside `test_y12mf_errors.f90` (for y12mf).  The analogous
three-step API reuse path is not tested directly.

### 2.6 Limited variety of benchmark matrices

Most tests use tridiagonal or arrow matrices.  Only one external benchmark
matrix (UMFPACK 5×5) is tested with the direct API before the new benchmark
driver is added.  Broader structural variety (block-diagonal, banded, general
unsymmetric) is missing.

### 2.7 y12mh/y12mg with non-tridiagonal inputs

`test_y12mg_mh.f90` tests the 1-norm and condition-number estimate only on
tridiagonal matrices of sizes 5 and 10.  Non-trivial sparsity patterns and
poorly conditioned matrices are not covered.

### 2.8 Mixed-precision comparison

There is no test that solves the same problem in single and double precision
and verifies that the DP solution is significantly more accurate (e.g., using
an ill-conditioned matrix where SP gives a large error but DP does not).

---

## 3. Suggested future tests

### SP-1  C-interface smoke test (`test_y12mae_c.c`)
Call `y12mae_c` from a C translation unit on the UMFPACK 5×5 matrix.
Verifies correct argument passing through the `bind(C)` wrappers and that
the `ifail` return value is zero.

### SP-2  C++ template smoke test (`test_y12m_hpp.cpp`)
Instantiate the C++ wrapper templates from `include/y12m.hpp` for `float`
and `double` on the same 5×5 matrix.  Checks that both specialisations
compile and produce the correct solution.

### SP-3  Large-n Markowitz overflow regression (`test_y12mc_large_n_dp.f90`)
Construct a tridiagonal system with n = 50 000 and solve it with
`y12mb + y12mc + y12md` (double precision).  The test passes if IFAIL = 0
and the solution error is below a tight tolerance; it fails on unpatched
builds due to 32-bit overflow in the Markowitz initialisation (`nr = n*n`).

### SP-4  Ill-conditioned system: SP vs DP accuracy comparison
Use a Hilbert matrix (size 6 or 8), which is notoriously ill-conditioned.
Solve in both SP and DP with y12ma.  Check that:
- The DP solution satisfies a tight residual tolerance.
- The SP solution is demonstrably less accurate (larger relative residual).
This documents the expected precision difference explicitly.

### SP-5  Iterative-refinement improvement on an ill-conditioned matrix
On a matrix where the direct SP solve has a large backward error, run y12mf
and verify that `AFLAG(9)/AFLAG(11)` (relative-error estimate) decreases
monotonically across refinement iterations and reaches a value below a
strict tolerance that the direct solve cannot meet.

### SP-6  LU structure-reuse with y12mb + y12mc + y12md
Set up two right-hand sides for the same sparsity structure.  Solve the
first with `IFLAG(4) = 0`, `IFLAG(5) = 1` to factorize.  Solve the second
with `IFLAG(4) = 2`, `IFLAG(5) = 3` to reuse the factorization.  Verify
that both solutions are correct and that the second call does not modify
the LU factors.

### SP-7  Near-singular matrix: IFAIL = 7 or 8 detection
Construct a nearly singular matrix (one row/column approximately linearly
dependent on another) and call y12mc with tight stability settings.  Verify
that IFAIL is non-zero (7 or 8) and that no undefined-behaviour occurs.

### SP-8  y12mh and y12mg with general unsymmetric matrices
Call y12mh and y12mg on the UMFPACK 5×5 and the Pardiso 8×8 matrices.
Verify that 1-norm values match hand-computed column sums, and that RCOND
is positive and consistent with a well-conditioned matrix.

### SP-9  Matrix Market regression suite
Read all `.mtx` files from the `test/data/` directory with `y12m_solve`,
compare the computed residuals against stored tolerances, and add at least
one non-square or structurally symmetric matrix to the suite.

### SP-10  FPM build test
Run `fpm test` to verify that the `fpm.toml` package description is
up-to-date and that all example programs compile and execute correctly
under the FPM build system.
