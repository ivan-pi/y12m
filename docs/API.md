# Y12M — API Reference

Y12M solves large, sparse systems of linear algebraic equations **Ax = b** by
Gaussian elimination with sparse matrix techniques and partial pivoting.

> **Disclaimer:** This document was reworked from the original Netlib `doc` file
> using AI-assisted tools.  The information may differ from, or contain errors
> relative to, the original text.  We have **not** cross-checked against the
> authoritative Springer publication:
>
> > Zlatev, Z., Wasniewski, J., & Schaumburg, K. (1981). *Y12M: Solution of
> > Large and Sparse Systems of Linear Algebraic Equations*. Lecture Notes in
> > Computer Science, Vol. 121. Springer.
> > https://doi.org/10.1007/3-540-10874-2
>
> Known discrepancies between this document and the `doc` file (or the source
> code) are listed in the [Discrepancies](#discrepancies) section.

---

**Contents:**

- [Overview](#overview)
  - [Common API points](#common-api-points)
  - [Subroutine summary](#subroutine-summary)
- [Procedure Reference](#procedure-reference)
  - [Input storage format](#input-storage-format)
  - [Y12MA](#y12ma) — single-system black-box driver
  - [Y12MB](#y12mb) — reorder matrix and prepare work arrays
  - [Y12MC](#y12mc) — sparse LU factorization
  - [Y12MD](#y12md) — triangular back-substitution
  - [Y12MF](#y12mf) — black-box driver with iterative refinement
  - [Y12MH](#y12mh) — compute matrix one-norm
  - [Y12MG](#y12mg) — reciprocal condition-number estimate
- [Flag arrays](#flag-arrays)
  - [AFLAG](#aflag)
  - [IFLAG](#iflag)
- [Error codes](#error-codes)
- [Discrepancies](#discrepancies)
- [References](#references)

---

## Overview

### Common API points

- **No internal allocation.** All working storage (`A`, `SNR`, `RNR`, `HA`,
  `PIVOT`, …) must be provided by the caller.  Arrays must be large enough to
  accommodate fill-in during factorization.

- **Settings via `AFLAG` and `IFLAG`.** Numerical tolerances, pivoting strategy,
  and reuse flags are passed through `AFLAG(8)` (real) and `IFLAG(10)` (integer).
  The black-box drivers Y12MA and Y12MF set their own defaults; when calling
  Y12MB / Y12MC / Y12MD directly the caller must initialize the relevant entries
  before each call.

- **Errors via `IFAIL`.** `IFAIL` is always the **last** argument of any
  subroutine that performs error checking.  On exit `IFAIL = 0` means success;
  a positive value means an error was detected and computation was stopped
  immediately.  See [Error codes](#error-codes) for the complete list.

- **Precision variants.** Each subroutine exists in a single-precision variant
  (suffix `E`, e.g. `Y12MAE`) and a double-precision variant (suffix `F`,
  e.g. `Y12MAF`).  The Fortran module `y12m` provides generic interfaces
  (`Y12MA`, `Y12MB`, …) that dispatch based on the actual-argument type.

- **Calling order.** The normal sequence for a single solve is
  `Y12MB` → `Y12MC` → `Y12MD`.  Y12MA and Y12MF wrap this sequence internally.
  The optional Y12MH (one-norm) must be called **before** Y12MC; the optional
  Y12MG (condition estimate) must be called **after** Y12MC.

- **State tracking via `IFLAG(1)`.** Set `IFLAG(1) ≥ 0` before the first call
  to the package; do not modify it between calls.  Y12MB sets it to `-1` on
  exit; Y12MC sets it to `-2`.  Y12MD requires it to be `-2` on entry.

- **Multiple right-hand sides / matrix reuse.** Use `IFLAG(4)` and `IFLAG(5)` to
  skip redundant factorization work when solving a sequence of systems.  See
  [docs/multiple_rhs.md](multiple_rhs.md) for worked examples.

### Subroutine summary

| Subroutine | Purpose |
|------------|---------|
| [`Y12MA`](#y12ma) | Black-box driver: prepare, factorize, and solve in one call |
| [`Y12MB`](#y12mb) | Prepare and reorder the sparse matrix |
| [`Y12MC`](#y12mc) | LU factorization |
| [`Y12MD`](#y12md) | Triangular solve (back-substitution) |
| [`Y12MF`](#y12mf) | Black-box driver with iterative refinement |
| [`Y12MH`](#y12mh) | One-norm of the original matrix (input to Y12MG) |
| [`Y12MG`](#y12mg) | Reciprocal condition-number estimate |

---

## Procedure Reference

### Input storage format

The non-zero elements of the coefficient matrix **A** are stored in
coordinate (COO) format on entry to Y12MB / Y12MA / Y12MF:

- `A(i)`, `i = 1…Z` — the non-zero values in *arbitrary* order.
- `SNR(i)` — column index of `A(i)`.
- `RNR(i)` — row index of `A(i)`.

Arrays `A` and `SNR` have length `NN ≥ 2*Z`; array `RNR` has length `NN1 ≥ Z`.
The extra space beyond `Z` is working storage for fill-in during factorization.

**Example** — `n = 5`, `Z = 12`:

```
    ┌                ┐
    │  5  0  0  3  0 │
    │  2  4  0  0  1 │
A = │  0  1  3  0  2 │
    │  0  0  0  2  3 │
    │  0  0  0  2  1 │
    └                ┘
```

The 12 non-zeros may be listed in any order.  A row-major ordering gives:

```
Index  :  1   2   3   4   5   6   7   8   9  10  11  12
A      :  5   3   2   4   1   1   3   2   2   3   2   1
SNR    :  1   4   1   2   5   2   3   5   4   5   4   5
RNR    :  1   1   2   2   2   3   3   3   4   4   5   5
```

Arrays `A`, `SNR`, and `RNR` must be allocated longer than `Z` to accommodate
fill-in.  With `NN = NN1 = 2*Z = 24` (the minimum), the layout looks like:

```
A  (1:NN=24) :  5  3  2  4  1  1  3  2  2  3  2  1  _  _  _  _  _  _  _  _  _  _  _  _
               |<──────── non-zeros (1:Z=12) ────────>|<──── workspace (Z+1:NN) ─────────>|

SNR(1:NN=24) :  1  4  1  2  5  2  3  5  4  5  4  5  _  _  _  _  _  _  _  _  _  _  _  _

RNR(1:NN1=24):  1  1  2  2  2  3  3  3  4  4  5  5  _  _  _  _  _  _  _  _  _  _  _  _
               |<──────── non-zeros (1:Z=12) ────────>|<──── workspace (Z+1:NN1) ─────────>|
```

The recommended allocation is `2*Z ≤ NN ≤ 3*Z` and `2*Z ≤ NN1 ≤ 3*Z`; for
matrices with high fill-in a larger `NN` may be needed.  Increase `NN` (or
`NN1`) when `iflag(6)` (or `iflag(7)`) is large after a solve.

---

### Y12MA

Black-box driver for a **single** sparse system with a **single** right-hand
side.  Internally calls Y12MB, Y12MC, and Y12MD after setting recommended
defaults for `AFLAG(1–4)` and `IFLAG(2–5)`.

```fortran
subroutine y12ma[e|f](n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, aflag, iflag, b, ifail)
  integer,             intent(in)    :: n, z, nn, nn1, iha
  real(kind=[sp,dp]),  intent(inout) :: a(nn), aflag(8)
  integer,             intent(inout) :: snr(nn), rnr(nn1), ha(iha,11), iflag(10)
  real(kind=[sp,dp]),  intent(out)   :: pivot(n)
  real(kind=[sp,dp]),  intent(inout) :: b(n)
  integer,             intent(out)   :: ifail
```

**Arguments:**

- `n`: Number of equations.
- `z`: Number of non-zero elements in **A**.
- `a(nn)`: On entry, the first `z` locations hold the non-zero values of **A**
  in arbitrary order.  On exit, holds the non-zeros of the upper triangular
  matrix **U** (diagonal elements in `pivot`).
- `snr(nn)`: Column indices of `a`.  On exit, updated to reflect the structure
  of **U**.
- `nn`: Length of `a` and `snr`.  Restriction: `nn ≥ 2*z`.
  Recommended: `2*z ≤ nn ≤ 3*z`.
- `rnr(nn1)`: Row indices of `a`.  On exit, normally all zero.
- `nn1`: Length of `rnr`.  Restriction: `nn1 ≥ z`.
  Recommended: `2*z ≤ nn1 ≤ 3*z`.
- `pivot(n)`: On exit, the diagonal elements of **U**.  Small values indicate
  numerical singularity (see also `aflag(8)`).
- `ha(iha,11)`: Integer work array.
- `iha`: First dimension of `ha`.  Restriction: `iha ≥ n`.
- `aflag(8)`: Real algorithm flags.  Y12MA **overwrites** `aflag(1–4)` with its
  own defaults.  See [AFLAG](#aflag).
- `iflag(10)`: Integer algorithm flags.  Y12MA **overwrites** `iflag(2–5)`.
  `iflag(1)` must be ≥ 0 before the first call.  See [IFLAG](#iflag).
- `b(n)`: On entry, the right-hand side **b**.  On exit, the solution **x**.
- `ifail`: Error indicator.  0 = success.  See [Error codes](#error-codes).

**Defaults set by Y12MA:**

These values reflect the recommended settings from the original documentation
and have been found satisfactory for a wide range of sparse matrices with
elements of order 1.  Y12MA **overwrites** these entries on every call, so
user-provided values are ignored.  When implementing a custom calling sequence
with Y12MB / Y12MC / Y12MD directly, these same values are a useful starting
point — only adjust them if you have specific knowledge of your matrix.

| Flag | Value | Rationale |
|------|-------|-----------|
| `aflag(1)` | 16.0 | Stability factor; limits pivot to elements ≥ max_row / 16 |
| `aflag(2)` | 1.0 × 10⁻¹² | Drop-tolerance; suited for well-scaled matrices (elements ≈ 1) |
| `aflag(3)` | 1.0 × 10¹⁶ | Growth-factor limit; high threshold minimises false stops |
| `aflag(4)` | 1.0 × 10⁻¹² | Singularity threshold relative to `aflag(6)` |
| `iflag(2)` | 2 | Markowitz search width; scan 2 rows with fewest non-zeros |
| `iflag(3)` | 1 | General Markowitz pivoting strategy |
| `iflag(4)` | 0 | No structure reuse (single-system solve) |
| `iflag(5)` | 1 | Discard **L** after factorization (saves memory for single-RHS solves) |

**Calls:** Y12MB, Y12MC, Y12MD.

---

### Y12MB

Reorders the non-zero elements of **A** by rows and prepares the auxiliary
array `ha` needed by Y12MC.

```fortran
subroutine y12mb[e|f](n, z, a, snr, nn, rnr, nn1, ha, iha, aflag, iflag, ifail)
  integer,            intent(in)    :: n, z, nn, nn1, iha
  real(kind=[sp,dp]), intent(inout) :: a(nn), aflag(8)
  integer,            intent(inout) :: snr(nn), rnr(nn1), ha(iha,11), iflag(10)
  integer,            intent(out)   :: ifail
```

**Arguments:**

- `n`: Number of equations.
- `z`: Number of non-zeros.
- `a(nn)`: On entry, the first `z` elements are the non-zeros in arbitrary
  order.  On exit, the first `z` elements are ordered **by row**.
- `snr(nn)`: Column indices of `a`.  Updated to match the row-ordered layout.
- `nn`: Length of `a` and `snr`.  Restriction: `nn ≥ 2*z`.
- `rnr(nn1)`: Row indices of `a`.  On exit, holds row numbers ordered
  **by column** (the column-linked list used by Y12MC).
- `nn1`: Length of `rnr`.  Restriction: `nn1 ≥ z`.
- `ha(iha,11)`: Integer work array.  On exit, columns 1, 3, 4, 6, 7, 8, and 11
  hold structural information required by Y12MC.
- `iha`: First dimension of `ha`.  Restriction: `iha ≥ n`.
- `aflag(8)`: Only `aflag(6)` is written (set to max|a(i,j)|); all other entries
  are ignored.
- `iflag(10)`: `iflag(1)` must be ≥ 0 on the first call; `iflag(4)` must be 0,
  1, or 2.  On exit `iflag(1) = -1`.
- `ifail`: Error indicator.  0 = success.  See [Error codes](#error-codes).

**`ha` columns on exit from Y12MB:**

| Column | Contents |
|--------|----------|
| `ha(i,1)` | Position in `a` of the first non-zero in row *i* |
| `ha(i,3)` | Position in `a` of the last non-zero in row *i* |
| `ha(i,4)` | Position in `rnr` of the first row-number for column *i* |
| `ha(i,6)` | Position in `rnr` of the last row-number for column *i* |
| 7, 8, 11 | Pivotal-search metadata (used by Y12MC) |

> Do **not** alter `n`, `z`, `a`, `snr`, `nn`, `rnr`, `nn1`, columns 1, 3, 4,
> 6, 7, 8, 11 of `ha`, `aflag(6)`, `iflag(1)`, `iflag(4)`, or `ifail` between
> calls to Y12MB and Y12MC.

---

### Y12MC

LU factorization of a sparse matrix **A**: computes **L** and **U** such that
**LU = PAQ** (P, Q permutation matrices).  Simultaneously transforms **b** to
**c = L⁻¹Pb**.

```fortran
subroutine y12mc[e|f](n, z, a, snr, nn, rnr, nn1, pivot, b, ha, iha, aflag, iflag, ifail)
  integer,            intent(in)    :: n, z, nn, nn1, iha
  real(kind=[sp,dp]), intent(inout) :: a(nn), b(n), aflag(8)
  integer,            intent(inout) :: snr(nn), rnr(nn1), ha(iha,11), iflag(10)
  real(kind=[sp,dp]), intent(out)   :: pivot(n)
  integer,            intent(out)   :: ifail
```

**Arguments:**

- `n`: Number of equations.
- `z`: Number of non-zeros (as on exit from Y12MB).
- `a(nn)`: On entry, non-zeros ordered by row (from Y12MB).  On exit, holds the
  off-diagonal non-zeros of **U** and—when `iflag(5) = 2`—also the off-diagonal
  non-zeros of **L**.
- `snr(nn)`: Column indices of `a`.  Updated to reflect **U** (and optionally
  **L**).
- `nn`: Length of `a` and `snr`.  Recommended: `2*z ≤ nn ≤ 3*z` when
  `iflag(5) = 1`; `3*z ≤ nn ≤ 5*z` when `iflag(5) = 2`.
- `rnr(nn1)`: Column-linked row-number list from Y12MB.  Normally all zero on
  exit.
- `nn1`: Length of `rnr`.
- `pivot(n)`: On exit, the diagonal elements of **U**.
- `b(n)`: On entry, the right-hand side **b**.  On exit, the transformed vector
  **c = L⁻¹Pb**.
- `ha(iha,11)`: Work array prepared by Y12MB; further updated by Y12MC.
- `iha`: First dimension of `ha`.  Restriction: `iha ≥ n`.
- `aflag(8)`: Numerical tolerances and output diagnostics.  See [AFLAG](#aflag).
- `iflag(10)`: Algorithm flags.  `iflag(1)` must be `-1` on entry (set by
  Y12MB); set to `-2` on successful exit.  See [IFLAG](#iflag).
- `ifail`: Error indicator.  0 = success.  See [Error codes](#error-codes).

> Do **not** alter `n`, `a`, `snr`, `nn`, `b`, `pivot`, columns 1, 2, 3, 7, 8
> of `ha`, `iha`, `aflag`, `iflag(1)`, `iflag(3)`, `iflag(4)`, or `ifail`
> between calls to Y12MC and Y12MD.

---

### Y12MD

Solves **Ax = b** using the LU factorization produced by Y12MC.

```fortran
subroutine y12md[e|f](n, a, nn, b, pivot, snr, ha, iha, iflag, ifail)
  integer,            intent(in)    :: n, nn, iha
  real(kind=[sp,dp]), intent(in)    :: a(nn), pivot(n)
  integer,            intent(in)    :: snr(nn), ha(iha,11), iflag(10)
  real(kind=[sp,dp]), intent(inout) :: b(n)
  integer,            intent(out)   :: ifail
```

**Arguments:**

- `n`: Number of equations.
- `a(nn)`: Non-zeros of **U** (and optionally **L**) as left by Y12MC.  Not
  modified.
- `nn`: Length of `a` and `snr`.
- `b(n)`: On entry (when `iflag(5) ≠ 3`): the vector **c = L⁻¹Pb** from Y12MC.
  On entry (when `iflag(5) = 3`): a new right-hand side **b** (reusing the
  existing LU factorization).  On exit: the solution **x**.
- `pivot(n)`: Diagonal of **U** from Y12MC.  Not modified.
- `snr(nn)`: Column indices from Y12MC.  Not modified.
- `ha(iha,11)`: Work array from Y12MC.  Not modified.
- `iha`: First dimension of `ha`.  Restriction: `iha ≥ n`.
- `iflag(10)`: `iflag(1)` must be `-2`.  Use the same values of `iflag(3)` and
  `iflag(4)` as in the preceding Y12MC call.  `iflag(5)`: 1 or 2 if Y12MC was
  called for this system; 3 to reuse an existing factorization (supply only a
  new `b`).  Not modified on exit.
- `ifail`: Error indicator.  0 = success.  See [Error codes](#error-codes).

---

### Y12MF

Black-box driver with **iterative refinement**.  Solves **Ax = b** and
improves the solution by computing successive corrections
**d(k) = QU⁻¹L⁻¹Pr(k-1)** where **r(k-1) = b − Ax(k-1)**.

> Residuals must be accumulated in higher than working precision for the
> refinement to be meaningful.  `Y12MFE` (single precision) accumulates in
> double precision and rounds back to single.  `Y12MFF` (double precision)
> accumulates in **quad precision** (`selected_real_kind(33)`) when the
> compiler provides it, and otherwise falls back to a **double-double**
> compensated accumulation (TwoProd + TwoSum; Dekker's splitting by
> default, or the Fortran 2018 `IEEE_FMA` intrinsic with
> `-DY12M_USE_IEEE_FMA=ON`).  The strategy is fixed at compile time; the
> module constants `y12mff_uses_quad` (logical) and
> `y12mff_accumulator_kind` (integer kind value) report which one was
> built in.  The CMake option `-DY12M_FORCE_DOUBLE_DOUBLE=ON` forces the
> double-double path even when quad precision is available (useful for
> reproducibility across compilers and for testing the fallback); the
> configure step prints which accumulator was selected.

```fortran
subroutine y12mf[e|f](n, a, snr, nn, rnr, nn1, a1, sn, nz, ha, iha, b, b1, x, y, aflag, iflag, ifail)
  integer,            intent(in)    :: n, nn, nn1, nz, iha
  real(kind=[sp,dp]), intent(inout) :: a(nn), aflag(11)
  integer,            intent(inout) :: snr(nn), rnr(nn1), ha(iha,13), sn(nz), iflag(12)
  real(kind=[sp,dp]), intent(inout) :: a1(nz), b(n)
  real(kind=[sp,dp]), intent(out)   :: b1(n), x(n), y(n)
  integer,            intent(out)   :: ifail
```

**Arguments:**

- `n`: Number of equations.
- `a(nn)`: When factorizing (`iflag(5) = 2`): first `nz` entries hold the
  non-zeros of **A** in arbitrary order; modified by Y12MF.  When reusing LU
  (`iflag(5) = 3`): holds the non-zeros of **U** and **L**; unchanged.
- `snr(nn)`: Column indices of `a`.  Modified when factorizing.
- `nn`: Length of `a` and `snr`.  Restriction: `nn ≥ 2*nz`.
- `rnr(nn1)`: Row indices of `a`.  Ignored when reusing LU.
- `nn1`: Length of `rnr`.  Restriction: `nn1 ≥ nz`.
  Recommended: `1.5*nz ≤ nn1 ≤ 2*nz`.
- `a1(nz)`: Copy of the original non-zeros ordered by row (for residual
  computation).  Set by Y12MF on first factorization.  Do not alter between
  successive calls.
- `sn(nz)`: Column indices of `a1`.  Set by Y12MF on first factorization.
- `nz`: Number of non-zeros in **A**.
- `ha(iha,13)`: Work array (columns 1–11 as in Y12MC; columns 12–13 hold
  row-start/end positions of the original matrix in `a1`/`sn`).
- `iha`: First dimension of `ha`.  Restriction: `iha ≥ n`.
- `b(n)`: On entry: right-hand side **b**.  On exit: last correction vector
  **d(p-1)**.
- `b1(n)`: On exit: the original right-hand side **b** (saved internally).
- `x(n)`: On exit: the corrected solution **x(p)**.
- `y(n)`: On exit: the pivot elements (diagonal of **U**).  Small values
  indicate numerical singularity (see also `aflag(8)`).
- `aflag(11)`: Real algorithm flags and output diagnostics.  See [AFLAG](#aflag).
- `iflag(12)`: Integer algorithm flags.  See [IFLAG](#iflag).
- `ifail`: Error indicator.  0 = success.  See [Error codes](#error-codes).

**Calls:** Y12MB, Y12MC, Y12MD.

---

### Y12MH

Computes the **one-norm** of sparse matrix **A** for use as `anorm` in Y12MG.
Must be called **before** Y12MC, which overwrites array `a`.

```fortran
subroutine y12mh[e|f](n, nz, a, snr, work, anorm)
  integer,            intent(in)  :: n, nz
  real(kind=[sp,dp]), intent(in)  :: a(nz)
  integer,            intent(in)  :: snr(nz)
  real(kind=[sp,dp]), intent(out) :: work(n), anorm
```

**Arguments:**

- `n`: Matrix dimension.
- `nz`: Number of non-zeros.
- `a(nz)`: Non-zero values of **A**.  Not modified.
- `snr(nz)`: Column indices of `a`.  Not modified.
- `work(n)`: Work array; contents are not meaningful on exit.
- `anorm`: On exit, the one-norm (maximum column-sum of absolute values) of **A**.

> Y12MH has **no** `ifail` parameter and performs no error checking.

---

### Y12MG

Estimates the **reciprocal condition number** of **A** using the LINPACK
algorithm.  Must be called **after** Y12MC, while the LU factorization is
intact.  Requires `iflag(5) ≥ 2` (matrix **L** must have been retained).

```fortran
subroutine y12mg[e|f](n, nn, a, snr, w, pivot, anorm, rcond, iha, ha, iflag, ifail)
  integer,            intent(in)    :: n, nn, iha
  real(kind=[sp,dp]), intent(in)    :: a(nn), pivot(n), anorm
  integer,            intent(in)    :: snr(nn), ha(iha,3), iflag(5)
  real(kind=[sp,dp]), intent(out)   :: w(n), rcond
  integer,            intent(inout) :: ifail
```

**Arguments:**

- `n`: Matrix dimension.
- `nn`: Length of `a` and `snr`.
- `a(nn)`: Non-zeros of **U** and **L** as left by Y12MC.
- `snr(nn)`: Column indices from Y12MC.
- `w(n)`: Work array; overwritten during computation.
- `pivot(n)`: Diagonal of **U** from Y12MC.
- `anorm`: One-norm of the original **A** (from Y12MH).
- `rcond`: On exit, an estimate of the reciprocal condition number.  The
  condition number is approximately `1/rcond`.  Set to −1 on error.
- `iha`: First dimension of `ha`.
- `ha(iha,3)`: Uses only the first three columns of the `ha` array from Y12MC.
- `iflag(5)`: Only `iflag(5)` is examined: if `iflag(5) = 1` (L was discarded),
  `ifail` is set to 26.
- `ifail`: On entry, checked; if already non-zero, returns immediately.  On
  exit: 0 = success; 26 = L was discarded.  See [Error codes](#error-codes).

> **Note:** The argument order places `iha` **before** `ha`, which differs from
> the convention used in Y12MB, Y12MC, and Y12MD.  Y12MG only reads the first
> three columns of `ha` and the first five entries of `iflag`.  However, when
> Y12MG is called as part of the Y12MB → Y12MC → Y12MG workflow, `ha` and
> `iflag` must still be allocated at their full sizes (`(iha,11)` and `(10)`
> respectively) because those larger extents are required by Y12MB and Y12MC.
> Declaring them smaller solely for Y12MG would be incorrect.

---

## Flag arrays

The `AFLAG` and `IFLAG` arrays are the primary mechanism for communicating
numerical tolerances, pivoting strategy, storage options, and diagnostic
results between the caller and the solver.  The caller initializes the
relevant input entries before each call; on return the output entries contain
diagnostics that can be used to assess solution quality and tune future calls.

### AFLAG

`AFLAG` is a real array of length 8 (or 11 for Y12MF).  Entries 1–4 are
**input** tolerances; entries 5–8 are **output** diagnostics.

| Index | In/Out | Description |
|-------|--------|-------------|
| `aflag(1)` | in | **Stability factor.** An element is accepted as a pivot only if its absolute value is at least `max_row / aflag(1)`. Must be > 1.0 (else overridden to 1.0005). Recommended: 4.0–16.0. |
| `aflag(2)` | in | **Drop-tolerance.** Elements whose absolute value falls below this threshold during elimination are discarded. Recommended: 1.0×10⁻¹² for well-scaled matrices. |
| `aflag(3)` | in | **Growth-factor threshold.** Elimination stops with `ifail = 4` when `aflag(5) > aflag(3)`. Must be ≥ 1.0×10⁵ (else overridden). Recommended: 1.0×10⁶ for Y12MC (direct use); 1.0×10¹⁶ for Y12MF. |
| `aflag(4)` | in | **Singularity threshold.** Elimination stops with `ifail = 3` when a pivot `\|a(i,i)\| < aflag(4)*aflag(6)`. Recommended: 1.0×10⁻¹². |
| `aflag(5)` | out | **Growth factor** = `aflag(7)/aflag(6)` after each elimination step. Large values indicate potentially inaccurate results. |
| `aflag(6)` | out | **Maximum element** of the original matrix **A** (set by Y12MB). |
| `aflag(7)` | out | **Running maximum** element encountered during any elimination step. |
| `aflag(8)` | out | **Minimum pivot** (absolute value). Small values indicate numerical singularity. |

Additional entries used only by Y12MF:

| Index | In/Out | Description |
|-------|--------|-------------|
| `aflag(9)` | out | Max-norm of the last correction vector **d(p-1)**. Should be small for an accurate solution. |
| `aflag(10)` | out | Max-norm of the last residual vector **r(p-1)**. Should be small for an accurate solution. |
| `aflag(11)` | out | Max-norm of the corrected solution **x(p)**. The ratio `aflag(9)/aflag(11)` gives a relative-error estimate. |

### IFLAG

`IFLAG` is an integer array of length 10 (or 12 for Y12MF).

| Index | In/Out | Description |
|-------|--------|-------------|
| `iflag(1)` | inout | **State flag.** Must be ≥ 0 before the first call; do not modify between calls.  Set to `-1` by Y12MB and to `-2` by Y12MC. |
| `iflag(2)` | in | **Markowitz search width.** Pivotal search is carried out in the `iflag(2)` rows with fewest non-zeros. Ignored when `iflag(3) = 0`. Recommended: ≤ 3. |
| `iflag(3)` | in | **Pivoting strategy.** 0 = no pivoting; 1 = general Markowitz (recommended); 2 = diagonal elements only. |
| `iflag(4)` | in | **Structure reuse.** 0 = no reuse; 1 = first system of a same-structure sequence; 2 = subsequent system in the same-structure sequence. |
| `iflag(5)` | in | **L storage.** 1 = discard **L** after factorization; 2 = keep **L**; 3 = reuse existing LU (Y12MD only — do not call Y12MB/Y12MC). |
| `iflag(6)` | out | **Row garbage collections.** If large, increase `nn` for faster performance. |
| `iflag(7)` | out | **Column garbage collections.** If large, increase `nn1` for faster performance. |
| `iflag(8)` | out | **Peak non-zero count** in `a` at any elimination step.  If much less than `nn`, a smaller `nn` could be used. |
| `iflag(9)` | out | Minimum `nn1` that avoids column garbage collections (written when `iflag(4) = 1`). |
| `iflag(10)` | out | Minimum `nn` that avoids row garbage collections (written when `iflag(4) = 1`). |

Additional entries used only by Y12MF:

| Index | In/Out | Description |
|-------|--------|-------------|
| `iflag(11)` | in | **Maximum refinement iterations.** Restriction: `iflag(11) ≥ 2`. Recommended: < 33. |
| `iflag(12)` | out | **Actual iterations performed.** |

---

## Error codes

`IFAIL = 0` on exit means success.  A positive value indicates an error; the
computation was stopped at the point of detection.

| `ifail` | Cause and remedy |
|---------|-----------------|
| 1 | `Y12MD` was called without a preceding `Y12MC` call for the current system (including the first system in a same-coefficient-matrix sequence Ax₁=b₁, Ax₂=b₂, …).  Set `iflag(1) ≥ 0` before the first call to any routine in the package to ensure reliable detection. |
| 2 | `Y12MC` was called without a preceding `Y12MB` call.  As with `ifail = 1`, set `iflag(1) ≥ 0` before the first package call. |
| 3 | A pivot `\|a(i,i)\| < aflag(4)*aflag(6)` was selected.  With a sufficiently small `aflag(4)` this indicates numerical singularity of the coefficient matrix. |
| 4 | The growth factor `aflag(5)` exceeded the threshold `aflag(3)`.  With a sufficiently large `aflag(3)` this indicates that matrix elements grow so rapidly during factorization that continuing is not justified.  A smaller stability factor `aflag(1)` may give better results. |
| 5 | The length `nn` of arrays `a` and `snr` is too small.  Increase `nn` (and possibly `nn1`). |
| 6 | The length `nn1` of array `rnr` is too small.  Increase `nn1` (and possibly `nn`). |
| 7 | A row with only zero elements in its active part was found during factorization.  With a sufficiently small drop-tolerance `aflag(2)` this indicates numerical singularity.  With a large `aflag(2)` singularity is uncertain; re-run with a smaller `aflag(2)` and/or check `aflag(8)` and `aflag(5)`. |
| 8 | A column with only zero elements in its active part was found during factorization.  The same interpretation and remedies apply as for `ifail = 7`. |
| 9 | A pivot element is missing.  This can occur when `aflag(2) > 0` and `iflag(4) = 2` (a positive drop-tolerance used for a subsequent system in a same-structure sequence); decrease the drop-tolerance and refactorize.  The error can also occur when a special pivoting strategy (`iflag(3) = 0` or `iflag(3) = 2`) is used with an unsuitable matrix. |
| 10 | `Y12MF` was called with `iflag(5) = 1`, which requests removal of the lower triangular factor **L**.  Use `iflag(5) = 2` to retain **L**. |
| 11 | The coefficient matrix **A** has at least two entries at the same position (i, j).  Check the input data. |
| 12 | The number of equations satisfies `n < 2`.  Check the value of `n`. |
| 13 | The non-zero count is non-positive (`z ≤ 0`, or `nz ≤ 0` in Y12MF).  Check the parameter `z` (renamed `nz` in Y12MF). |
| 14 | The non-zero count is less than the number of equations (`z < n`).  If `z` (or `nz` in Y12MF) is correctly set, the coefficient matrix is structurally singular. |
| 15 | The leading dimension `iha` of array `ha` is less than `n`.  Set `iha ≥ n`. |
| 16 | `iflag(4)` has an invalid value; it must be 0, 1, or 2.  See the description of `iflag(4)` for details. |
| 17 | A row with all zero elements was found in **A** before Gaussian elimination begins — the matrix is structurally singular. |
| 18 | A column with all zero elements was found in **A** before Gaussian elimination begins — the matrix is structurally singular. |
| 19 | `iflag(2) < 1`.  Set `iflag(2)` to a positive integer (`iflag(2) = 3` is recommended). |
| 20 | `iflag(3)` is out of range; it must be 0, 1, or 2. |
| 21 | `iflag(5)` is out of range; it must be 1, 2, or 3.  When `iflag(5) = 3`, call only `Y12MD` — do not call `Y12MB` or `Y12MC` (see also `ifail = 22`). |
| 22 | `Y12MC` was called with `iflag(5) = 3`.  Call `Y12MC` with `iflag(5) = 1` or `2`. |
| 23 | The maximum number of refinement iterations `iflag(11)` (Y12MF only) is less than 2.  Set `iflag(11) ≥ 2`. |
| 24 | At least one element has a column index outside [1, n]. |
| 25 | At least one element has a row index outside [1, n]. |
| 26 | `Y12MG` requires the lower triangular factor **L** to be available, but **L** was discarded during factorization (`iflag(5) = 1` was used).  Retain **L** by setting `iflag(5) = 2`. |

---

## Discrepancies

The following issues were found by cross-checking the original `doc` file
against the source files in `src/legacy/`.

1. **`aflag(3)` default in Y12MA — doc says 1.0×10⁶, code sets 1.0×10¹⁶.**
   Both `y12mae` and `y12maf` (see `y12ma.f`) assign `aflag(3) = 1.0e+16` / `1.0d+16`.
   This document uses the value found in the source.

2. **`iflag(2)` default in Y12MA — doc says 3, code sets 2.**
   Both `y12mae` and `y12maf` (see `y12ma.f`) assign `iflag(2) = 2`.  This document uses
   the value found in the source.

3. **`ifail = 26` absent from the original error-diagnostics section.**
   Y12MG can return `ifail = 26` when **L** has been discarded, but this code
   does not appear in the `doc` file.  It is documented here from `y12mg.f`.

4. **Y12MG and Y12MH not covered by the `doc` file.**
   The original `doc` file documents only Y12MA, Y12MB, Y12MC, Y12MD, and Y12MF.
   Y12MG and Y12MH are present in `src/legacy/` but undocumented in `doc`.

5. **Y12MG: argument order — `iha` precedes `ha`.**
   All other subroutines declare `…, ha, iha, …`; Y12MG reverses this to
   `…, anorm, rcond, iha, ha, …`.

6. **Y12MG uses smaller array extents than the other subroutines.**
   The source declares `ha(iha,3)` (only 3 columns) and `iflag(5)` (only
   5 elements), while every other subroutine uses `ha(iha,11)` and
   `iflag(10)`.

7. **Typo in Y12MF's `iflag(2)` description in `doc`.**
   The original text reads "If `IFLAG((2) .ge. 0` then the pivotal search…"
   (mismatched parenthesis).  The intended condition is consistent with
   Y12MC's description: the search is performed when `iflag(3) ≥ 1`.

8. **`IFAIL = 22` — Y12MB does not check `iflag(5) = 3`.**
   The original documentation (and the problem description) states that both
   Y12MB and Y12MC return `IFAIL = 22` when called with `iflag(5) = 3`.
   However, in the source code only Y12MC (`y12mce`/`y12mcf` in `y12mc.f`)
   performs this check; Y12MB (`y12mbe`/`y12mbf` in `y12mb.f`) contains no
   such guard.
   Consequently, calling Y12MB with `iflag(5) = 3` will silently proceed
   rather than return `IFAIL = 22`.  The description of `IFAIL = 22` in this
   document reflects what the code actually does.

---

## References

See [references.md](references.md) for the full reference list from the original
Netlib documentation (40 entries).
