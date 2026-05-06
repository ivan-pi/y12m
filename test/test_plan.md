# Y12M Test Plan

This document analyses the current test coverage of the y12m sparse-solver
package and lists tests worth adding in the future.

---

## 1. Summary of current tests

Files with a `.F90` extension are compiled twice via a `TEST_SINGLE_PRECISION`
preprocessor macro, producing separate `_sp` and `_dp` CTest executables.
Files with a `.f90` extension produce a single executable.

| Test file | CTest name(s) | Routine(s) exercised | Precision | Matrix types / notes |
|---|---|---|---|---|
| `test_y12ma_sp.f90` | `test_y12ma_sp` | y12ma | SP | n=1 diagonal (IFAIL=12 error path), n=2–4 tridiagonal, n=5 arrow |
| `test_y12ma_dp.f90` | `test_y12ma_dp` | y12ma | DP | n=6–7,9–10 tridiagonal, n=8 arrow |
| `test_y12ma_5x5.F90` | `test_y12ma_5x5_sp`, `test_y12ma_5x5_dp` | y12ma | SP + DP | UMFPACK 5×5; residual check (2-norm and ∞-norm) |
| `test_y12ma_rand.F90` | `test_y12ma_rand_sp`, `test_y12ma_rand_dp` | y12ma | SP + DP | random 9×9 diag.-dominant, 3 RHS; cross-checked vs LAPACK sgesv/dgesv (optional, requires LAPACK) |
| `test_y12ma_hilbert.f90` | `test_y12ma_hilbert` | y12ma | SP + DP | Hilbert matrices: compares SP vs DP forward/backward errors at n=8 (no fixed cutoffs), and checks SP singularity detection on larger n |
| `test_y12mb_mc_md_sp.f90` | `test_y12mb_mc_md_sp` | y12mb + y12mc + y12md | SP | n=1 diagonal (IFAIL=12 error path), n=2 full 2×2, n=3–4 tridiagonal, n=5 arrow |
| `test_y12mb_mc_md_dp.f90` | `test_y12mb_mc_md_dp` | y12mb + y12mc + y12md | DP | n=6,8,9 tridiagonal, n=7 arrow, n=10 pentadiagonal (bandwidth 2) |
| `test_y12mb_errors.F90` | `test_y12mb_errors_sp`, `test_y12mb_errors_dp` | y12mb | SP + DP | error-path coverage (IFAIL 5,6,11–18,24,25); valid-call postconditions (IFLAG(1)=−1, AFLAG(6)=max\|A\|) |
| `test_y12mc_errors.F90` | `test_y12mc_errors_sp`, `test_y12mc_errors_dp` | y12mc | SP + DP | error-path coverage (IFAIL 2–4,7–9,19–22); success block verifies IFLAG(1)=−2, AFLAG(6,8), pivot sequence, and residual |
| `test_y12mc_z_intent.F90` | `test_y12mc_z_intent_sp`, `test_y12mc_z_intent_dp` | y12mc | SP + DP | regression: `z` must be unchanged on success path, on early-error exit, and when passed as an integer literal |
| `test_y12md_errors.F90` | `test_y12md_errors_sp`, `test_y12md_errors_dp` | y12mb + y12mc + y12md | SP + DP | error-path coverage (IFAIL=1); success with IFLAG(5)=1 and IFLAG(5)=2; multiple-RHS reuse (IFLAG(5)=3) on 6×6 reference matrix |
| `test_y12mf_sp.f90` | `test_y12mf_sp` | y12mf | SP | tridiagonal n=3,5,7,10 |
| `test_y12mf_errors.f90` | `test_y12mf_errors` | y12mf | SP | 6×6 reference; error-path (IFAIL=10,23); output-parameter invariants (IFLAG(12), AFLAG(9–11), B, B1); LU-reuse (IFLAG(5)=3) |
| `test_y12mg_mh.f90` | `test_y12mg_mh` | y12mg + y12mh | SP + DP | 1-norm of n=5 (SP) and n=10 (DP) tridiagonal; condition estimate after full y12mb→y12mc→y12md→y12mg sequence |
| `test_y12ma_y12mf_bench_sp.f90` | `test_y12ma_y12mf_bench_sp` | y12ma + y12mf | SP | UMFPACK 5×5, Harwell-Boeing 4×4, Pardiso 8×8, Templates 6×6 |
| `y12m_solve.f90` | `y12m_solve_mat0` … `y12m_solve_mat5` (no mat4) | y12ma (DP) | DP | five `.mtx` data files (mat0–mat3, mat5); file-based regression driver with `--verbose` output |
| `y12m_mm.f90` | *(excluded when `WITH_NIST_MMIO=OFF`)* | y12ma | SP + DP | arbitrary Matrix Market format; requires NIST mmio library |

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

`test_y12ma_hilbert` now improves mixed-precision coverage by comparing SP and
DP forward/backward errors without hardcoded absolute cutoffs and by checking
that single precision eventually triggers singularity detection for larger
Hilbert sizes.  Iterative-refinement benefit on ill-conditioned matrices
(SP-5 from section 3) remains future work.

### 2.5 Limited variety of benchmark matrices

Most tests use tridiagonal or arrow matrices.  `test_y12ma_y12mf_bench_sp`
adds four literature matrices (UMFPACK 5×5, Harwell-Boeing 4×4, Pardiso 8×8,
Templates 6×6) but only in single precision and only for y12ma/y12mf.  Broader
structural variety (block-diagonal, banded, general unsymmetric) at double
precision and through the y12mb/y12mc/y12md API is still missing.

### 2.6 y12md input-validation checks not implemented

`test_y12md_errors.F90` documents (as commented-out blocks) three input checks
that `y12md` does not currently perform: `NN < 2*Z` (IFLAG(5)=1), `NN < 3*Z`
(IFLAG(5)=2), and `IHA < N`.  These tests cannot be activated until the
corresponding validation logic is added to the `y12md` source.

### 2.7 y12mh/y12mg with non-tridiagonal inputs

`test_y12mg_mh.f90` tests the 1-norm and condition-number estimate only on
tridiagonal matrices of sizes 5 and 10.  Non-trivial sparsity patterns and
poorly conditioned matrices are not covered.

### 2.8 Mixed-precision comparison

There is no test that solves the same problem in single and double precision
and verifies that the DP solution is significantly more accurate (e.g., using
an ill-conditioned matrix where SP gives a large error but DP does not).

### 2.9 Disconnected examples and test suite

The `examples/` directory contains drivers (`maind.f`, `maine.f`, `mainf.f`,
`mainf2.f`) and a unit test (`test_matrix_generators.f90`) for the
class D, E and F sparse matrix generators.  These generators can produce
parameterised matrices of varying size and density.  They are registered as
CTest tests but are separate from the main `test/` suite.  Bringing the
class-D/E/F generators into the test suite would let the automated tests
cover a much wider range of matrix structures and sizes without duplicating
test infrastructure.

### 2.10 y12m_solve missing mat4 data file

The `y12m_solve` driver is registered for mat0, mat1, mat2, mat3, and mat5,
but `test/data/mat4.mtx` does not exist.  Adding this data file (or explaining
why it is intentionally absent) would complete the mat0–mat5 series.

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

### SP-4  ~~Ill-conditioned system: SP vs DP accuracy comparison~~  **DONE**
~~Use a Hilbert matrix (size 6 or 8), which is notoriously ill-conditioned.
Solve in both SP and DP with y12ma.  Check that:
- The DP solution satisfies a tight residual tolerance.
- The SP solution is demonstrably less accurate (larger relative residual).
This documents the expected precision difference explicitly.~~

Implemented in `test_y12ma_hilbert.f90`.  It compares SP and DP forward and
backward errors at n=8 without hardcoded absolute thresholds, and also scans
larger Hilbert sizes to ensure SP singularity detection is exercised.

### SP-5  Iterative-refinement improvement on an ill-conditioned matrix
On a matrix where the direct SP solve has a large backward error, run y12mf
and verify that `AFLAG(9)/AFLAG(11)` (relative-error estimate) decreases
monotonically across refinement iterations and reaches a value below a
strict tolerance that the direct solve cannot meet.

### SP-6  Near-singular matrix: IFAIL = 7 or 8 detection
Construct a nearly singular matrix (one row/column approximately linearly
dependent on another) and call y12mc with tight stability settings.  Verify
that IFAIL is non-zero (7 or 8) and that no undefined-behaviour occurs.

### SP-7  y12mh and y12mg with general unsymmetric matrices
Call y12mh and y12mg on the UMFPACK 5×5 and the Pardiso 8×8 matrices.
Verify that 1-norm values match hand-computed column sums, and that RCOND
is positive and consistent with a well-conditioned matrix.

### SP-8  `y12m_solve` regression suite with stored tolerances
The `test/data/` directory contains five sparse systems in `.mtx` format
(mat0–mat3, mat5).  Extend the `y12m_solve` driver to record expected
residual norms alongside each data file and fail if the computed residual
exceeds those stored tolerances.  Add `mat4.mtx` and at least one additional
problem to broaden structural coverage.

### SP-9  Integrate class D/E/F matrix generators into the test suite
Move the class D, E and F matrix generators from `examples/` into the shared
test infrastructure so that `y12ma`, `y12mb + y12mc + y12md`, and `y12mf`
are all automatically exercised on parameterised matrices of varying size and
density.  This would significantly increase structural diversity without
adding hand-crafted matrices.

### SP-10  y12md input-validation checks
Add `NN` and `IHA` precondition checks to the `y12md` implementation and
activate the commented-out test blocks in `test_y12md_errors.F90`
(see gap 2.6 above).  The expected error codes are IFAIL=5 for `NN` violations
and IFAIL=15 for `IHA < N`.
