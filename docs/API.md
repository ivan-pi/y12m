# Y12M API Reference

Y12M solves large, sparse systems of linear algebraic equations **Ax = b** by
Gaussian elimination with sparse matrix techniques.

## Common API points

- **No internal allocation.** No arrays are allocated internally. The caller
  must supply all working storage (`A`, `SNR`, `RNR`, `HA`, `PIVOT`, …) and
  must make them large enough to accommodate fill-in during factorization.
- **Settings via `AFLAG` and `IFLAG`.** Algorithm parameters (stability factor,
  drop-tolerance, growth-factor threshold, pivoting strategy, …) are passed and
  returned through the `AFLAG(8)` (real) and `IFLAG(10)` (integer) arrays.
  Y12MA and Y12MF set their own defaults; when calling Y12MB/Y12MC/Y12MD
  directly the caller must initialize the relevant entries.
- **Errors via `IFAIL`.** Every subroutine that can detect an error has `IFAIL`
  as its **last** argument. On entry `IFAIL` is ignored; on exit `IFAIL = 0`
  means success, `IFAIL > 0` indicates an error and the computations were
  terminated immediately. See the [Error Codes](#error-codes) section for the
  full list.
- **Precision variants.** Each subroutine comes in a single-precision variant
  (suffix **E**, e.g. `Y12MAE`) and a double-precision variant (suffix **F**,
  e.g. `Y12MAF`). The Fortran module `y12m` exposes generic interfaces (`Y12MA`,
  `Y12MB`, …) that dispatch based on the actual-argument type.
- **Calling order.** The normal sequence for one solve is:
  `Y12MB` → `Y12MC` → `Y12MD`.  `Y12MA` and `Y12MF` wrap this sequence
  internally.  The optional `Y12MH` (one-norm) must precede `Y12MC`; the
  optional `Y12MG` (condition estimate) must follow `Y12MC`.
- **State tracking via `IFLAG(1)`.** The package uses `IFLAG(1)` as an internal
  state flag.  Initialize `IFLAG(1) ≥ 0` before the first call; do not modify
  it between calls.  Y12MB sets it to `-1` on exit; Y12MC sets it to `-2`.
  Y12MD requires it to be `-2` on entry.
- **Multiple right-hand sides / matrix reuse.** Use `IFLAG(4)` and `IFLAG(5)`
  to avoid redundant work when solving several systems; see
  [docs/multiple_rhs.md](multiple_rhs.md).

---

## Subroutine index

| Subroutine | Purpose |
|------------|---------|
| [`Y12MA`](#y12ma) | Black-box driver: orders, factorizes, and solves in one call |
| [`Y12MB`](#y12mb) | Reorders and prepares the sparse matrix |
| [`Y12MC`](#y12mc) | LU factorization |
| [`Y12MD`](#y12md) | Triangular solves (back-substitution) |
| [`Y12MF`](#y12mf) | Black-box driver with iterative refinement |
| [`Y12MH`](#y12mh) | One-norm of the original matrix (for condition estimation) |
| [`Y12MG`](#y12mg) | Reciprocal condition-number estimate |

---

## Y12MA

### Purpose

Black-box driver for a **single** sparse system **Ax = b** with a **single**
right-hand side.  Internally calls Y12MB, Y12MC, and Y12MD after setting
recommended defaults for `AFLAG(1–4)` and `IFLAG(2–5)`.

### Signatures

```fortran
! Single precision
SUBROUTINE Y12MAE(N, Z, A, SNR, NN, RNR, NN1, PIVOT, HA, IHA, AFLAG, IFLAG, B, IFAIL)
REAL             :: A(NN), PIVOT(N), B(N), AFLAG(8)
INTEGER          :: SNR(NN), RNR(NN1), HA(IHA,11), IFLAG(10)
INTEGER          :: N, Z, NN, NN1, IHA, IFAIL

! Double precision
SUBROUTINE Y12MAF(N, Z, A, SNR, NN, RNR, NN1, PIVOT, HA, IHA, AFLAG, IFLAG, B, IFAIL)
DOUBLE PRECISION :: A(NN), PIVOT(N), B(N), AFLAG(8)
INTEGER          :: SNR(NN), RNR(NN1), HA(IHA,11), IFLAG(10)
INTEGER          :: N, Z, NN, NN1, IHA, IFAIL
```

### Parameters

| Name | Type | Intent | Description |
|------|------|--------|-------------|
| `N` | INTEGER | in | Number of equations. |
| `Z` | INTEGER | in | Number of non-zero elements in **A**. |
| `A` | REAL/DP array(`NN`) | inout | On entry: first `Z` locations hold the non-zero elements of the coefficient matrix in arbitrary order. On exit: non-zero elements of the upper triangular matrix **U** (diagonal in `PIVOT`). |
| `SNR` | INTEGER array(`NN`) | inout | On entry: `SNR(j)`, j=1…Z, is the column number of the non-zero element in `A(j)`. On exit: column numbers of the non-zero elements of **U**. |
| `NN` | INTEGER | in | Length of `A` and `SNR`. Restriction: `NN ≥ 2*Z`. Recommended: `2*Z ≤ NN ≤ 3*Z`. |
| `RNR` | INTEGER array(`NN1`) | inout | On entry: `RNR(i)`, i=1…Z, is the row number of the non-zero element in `A(i)`. On exit: normally all zero. |
| `NN1` | INTEGER | in | Length of `RNR`. Restriction: `NN1 ≥ Z`. Recommended: `2*Z ≤ NN1 ≤ 3*Z`. |
| `PIVOT` | REAL/DP array(`N`) | out | Diagonal elements of **U** on exit. Small values indicate numerical singularity (see also `AFLAG(8)`). |
| `HA` | INTEGER array(`IHA`,11) | work | Integer work array. |
| `IHA` | INTEGER | in | First dimension of `HA`. Restriction: `IHA ≥ N`. |
| `AFLAG` | REAL/DP array(8) | inout | Algorithm flags; see [AFLAG table](#aflag-y12ma). Y12MA **overwrites** `AFLAG(1–4)` with its own defaults. |
| `IFLAG` | INTEGER array(10) | inout | Algorithm flags; see [IFLAG table](#iflag-y12ma). Y12MA **overwrites** `IFLAG(2–5)` with its own defaults. `IFLAG(1)` must be ≥ 0 before the first call. |
| `B` | REAL/DP array(`N`) | inout | On entry: right-hand side vector **b**. On exit: computed solution vector **x**. |
| `IFAIL` | INTEGER | out | Error indicator. 0 = success. See [Error Codes](#error-codes). |

#### AFLAG — Y12MA {#aflag-y12ma}

| Index | Description | Default set by Y12MA |
|-------|-------------|----------------------|
| `AFLAG(1)` | **Stability factor.** A pivot is accepted only if `|a(i,j)| ≥ max_row / AFLAG(1)`. | 16.0 |
| `AFLAG(2)` | **Drop-tolerance.** Elements smaller than this (in absolute value) during elimination are discarded. | 1.0 × 10⁻¹² |
| `AFLAG(3)` | **Growth-factor threshold.** Elimination stops with `IFAIL=4` when `AFLAG(5) > AFLAG(3)`. | 1.0 × 10¹⁶ |
| `AFLAG(4)` | **Singularity threshold.** Elimination stops with `IFAIL=3` when the pivot `|a(i,i)| < AFLAG(4)*AFLAG(6)`. | 1.0 × 10⁻¹² |
| `AFLAG(5)` | **Growth factor** (output). Set to `AFLAG(7)/AFLAG(6)` after each elimination step. Large values indicate possible inaccuracy. | — |
| `AFLAG(6)` | **Maximum element** (output). Set by Y12MB to max\|a(i,j)\|. | — |
| `AFLAG(7)` | **Running maximum** (output). Largest element found during any elimination step. | — |
| `AFLAG(8)` | **Minimum pivot** (output). Smallest \|pivot\| encountered. Small values indicate numerical singularity. | — |

#### IFLAG — Y12MA {#iflag-y12ma}

| Index | Description | Default set by Y12MA |
|-------|-------------|----------------------|
| `IFLAG(1)` | Internal state flag. Must be ≥ 0 before first call; do not alter between calls. | (work space) |
| `IFLAG(2)` | **Markowitz search width.** Pivotal search is carried out in the `IFLAG(2)` rows with fewest non-zeros. | 2 |
| `IFLAG(3)` | **Pivoting strategy.** 0 = no pivoting; 1 = general (recommended); 2 = diagonal only. | 1 |
| `IFLAG(4)` | **Structure reuse flag.** 0 = no reuse; 1 = first system of a same-structure sequence; 2 = subsequent system. | 0 |
| `IFLAG(5)` | **L storage flag.** 1 = discard L after factorization; 2 = keep L. | 1 |
| `IFLAG(6)` | **Row garbage collections** (output). Large value → increase `NN`. | — |
| `IFLAG(7)` | **Column garbage collections** (output). Large value → increase `NN1`. | — |
| `IFLAG(8)` | **Peak non-zero count** (output). If much less than `NN`, `NN` can be reduced. | — |
| `IFLAG(9)` | Minimum `NN1` avoiding column garbage collections (set when `IFLAG(4)=1`). | — |
| `IFLAG(10)` | Minimum `NN` avoiding row garbage collections (set when `IFLAG(4)=1`). | — |

### Auxiliary subroutines

Calls Y12MB, Y12MC, and Y12MD.

---

## Y12MB

### Purpose

Reorders the non-zero elements of matrix **A** by rows and prepares the
auxiliary array `HA` for subsequent factorization by Y12MC.

### Signatures

```fortran
! Single precision
SUBROUTINE Y12MBE(N, Z, A, SNR, NN, RNR, NN1, HA, IHA, AFLAG, IFLAG, IFAIL)
REAL    :: A(NN), AFLAG(8)
INTEGER :: SNR(NN), RNR(NN1), HA(IHA,11), IFLAG(10)
INTEGER :: N, Z, NN, NN1, IHA, IFAIL

! Double precision
SUBROUTINE Y12MBF(N, Z, A, SNR, NN, RNR, NN1, HA, IHA, AFLAG, IFLAG, IFAIL)
DOUBLE PRECISION :: A(NN), AFLAG(8)
INTEGER          :: SNR(NN), RNR(NN1), HA(IHA,11), IFLAG(10)
INTEGER          :: N, Z, NN, NN1, IHA, IFAIL
```

### Parameters

| Name | Type | Intent | Description |
|------|------|--------|-------------|
| `N` | INTEGER | in | Number of equations. Unchanged on exit. |
| `Z` | INTEGER | in | Number of non-zeros. Unchanged on exit. |
| `A` | REAL/DP array(`NN`) | inout | On entry: first `Z` elements are the non-zeros in arbitrary order. On exit: first `Z` elements are the non-zeros ordered **by row**. |
| `SNR` | INTEGER array(`NN`) | inout | On entry: column numbers of non-zeros in `A`. On exit: column numbers of the row-ordered non-zeros. |
| `NN` | INTEGER | in | Length of `A` and `SNR`. Restriction: `NN ≥ 2*Z`. |
| `RNR` | INTEGER array(`NN1`) | inout | On entry: row numbers of non-zeros. On exit: row numbers ordered **by column** (column-linked list). |
| `NN1` | INTEGER | in | Length of `RNR`. Restriction: `NN1 ≥ Z`. |
| `HA` | INTEGER array(`IHA`,11) | inout | Work array. On exit columns 1, 3, 4, 6, 7, 8, 11 hold structural information needed by Y12MC. |
| `IHA` | INTEGER | in | First dimension of `HA`. Restriction: `IHA ≥ N`. |
| `AFLAG` | REAL/DP array(8) | inout | Only `AFLAG(6)` is set (to max\|a(i,j)\|); other entries are ignored. |
| `IFLAG` | INTEGER array(10) | inout | `IFLAG(1)` must be ≥ 0 on first call. `IFLAG(4)` must be 0, 1, or 2 (see Y12MC). On exit `IFLAG(1) = -1`. Other entries used as work space. |
| `IFAIL` | INTEGER | out | Error indicator. 0 = success. See [Error Codes](#error-codes). |

#### HA contents on exit from Y12MB

| Column | Contents |
|--------|----------|
| `HA(i,1)` | Position in `A` of the first non-zero element of row *i*. |
| `HA(i,3)` | Position in `A` of the last non-zero element of row *i*. |
| `HA(i,4)` | Position in `RNR` of the first row-number of column *i*. |
| `HA(i,6)` | Position in `RNR` of the last row-number of column *i*. |
| cols 7, 8, 11 | Pivotal-search metadata (for use by Y12MC). |

### Notes

- Do **not** alter `N`, `Z`, `A`, `SNR`, `NN`, `RNR`, `NN1`, columns 1, 3, 4,
  6, 7, 8, and 11 of `HA`, `AFLAG(6)`, `IFLAG(1)`, `IFLAG(4)`, or `IFAIL`
  between calls to Y12MB and Y12MC.
- If `IFAIL > 0` on exit there is no point calling Y12MC or Y12MD.

---

## Y12MC

### Purpose

LU factorization of a sparse matrix **A**: computes **L** and **U** such that
**LU = PAQ** (P, Q are permutation matrices).  The right-hand side vector **b**
is simultaneously transformed to **c = L⁻¹ P b**.

### Signatures

```fortran
! Single precision
SUBROUTINE Y12MCE(N, Z, A, SNR, NN, RNR, NN1, PIVOT, B, HA, IHA, AFLAG, IFLAG, IFAIL)
REAL    :: A(NN), PIVOT(N), B(N), AFLAG(8)
INTEGER :: SNR(NN), RNR(NN1), HA(IHA,11), IFLAG(10)
INTEGER :: N, Z, NN, NN1, IHA, IFAIL

! Double precision
SUBROUTINE Y12MCF(N, Z, A, SNR, NN, RNR, NN1, PIVOT, B, HA, IHA, AFLAG, IFLAG, IFAIL)
DOUBLE PRECISION :: A(NN), PIVOT(N), B(N), AFLAG(8)
INTEGER          :: SNR(NN), RNR(NN1), HA(IHA,11), IFLAG(10)
INTEGER          :: N, Z, NN, NN1, IHA, IFAIL
```

### Parameters

| Name | Type | Intent | Description |
|------|------|--------|-------------|
| `N` | INTEGER | in | Number of equations. Unchanged on exit. |
| `Z` | INTEGER | in | Number of non-zeros. Unchanged on exit. |
| `A` | REAL/DP array(`NN`) | inout | On entry: non-zeros ordered by row (output of Y12MB). On exit: non-zeros of **U** (off-diagonal), and—when `IFLAG(5)=2`—also non-zeros of **L** (off-diagonal). |
| `SNR` | INTEGER array(`NN`) | inout | Column numbers corresponding to `A`. Updated to reflect **U** (and optionally **L**) structure. |
| `NN` | INTEGER | in | Length of `A` and `SNR`. Restriction: `NN ≥ 2*Z`. Recommended: `2*Z ≤ NN ≤ 3*Z` when `IFLAG(5)=1`; `3*Z ≤ NN ≤ 5*Z` when `IFLAG(5)=2`. |
| `RNR` | INTEGER array(`NN1`) | inout | Column-ordered row-number list (from Y12MB). On exit: normally all zero. |
| `NN1` | INTEGER | in | Length of `RNR`. Restriction: `NN1 ≥ Z`. |
| `PIVOT` | REAL/DP array(`N`) | out | Diagonal elements of **U**. |
| `B` | REAL/DP array(`N`) | inout | On entry: right-hand side **b**. On exit: transformed vector **c = L⁻¹ P b**. |
| `HA` | INTEGER array(`IHA`,11) | inout | Work array prepared by Y12MB. Updated by Y12MC. |
| `IHA` | INTEGER | in | First dimension of `HA`. Restriction: `IHA ≥ N`. |
| `AFLAG` | REAL/DP array(8) | inout | See [AFLAG table](#aflag-y12mc). |
| `IFLAG` | INTEGER array(10) | inout | See [IFLAG table](#iflag-y12mc). On exit `IFLAG(1) = -2`. |
| `IFAIL` | INTEGER | out | Error indicator. 0 = success. See [Error Codes](#error-codes). |

#### AFLAG — Y12MC {#aflag-y12mc}

| Index | In/Out | Description |
|-------|--------|-------------|
| `AFLAG(1)` | in | Stability factor. Must be > 1.0; if not, set to 1.0005. Recommended: 4.0–16.0. Unchanged on exit. |
| `AFLAG(2)` | in | Drop-tolerance. Positive small number or zero. Recommended: 1.0×10⁻¹² for well-scaled matrices. Unchanged on exit. |
| `AFLAG(3)` | in | Growth-factor threshold. Must be ≥ 1.0×10⁵; smaller values are overridden. Recommended: 1.0×10⁶. Unchanged on exit. |
| `AFLAG(4)` | in | Singularity threshold multiplier. Non-negative; negative values replaced by their absolute value. Recommended: 1.0×10⁻¹². Unchanged on exit. |
| `AFLAG(5)` | out | Growth factor = `AFLAG(7)/AFLAG(6)` after each step. Large values indicate potential inaccuracy. |
| `AFLAG(6)` | in | Maximum element of original **A** (set by Y12MB). Unchanged on exit. |
| `AFLAG(7)` | out | Maximum element encountered during elimination. |
| `AFLAG(8)` | out | Minimum pivot (absolute value). Small values indicate numerical singularity. |

#### IFLAG — Y12MC {#iflag-y12mc}

| Index | In/Out | Description |
|-------|--------|-------------|
| `IFLAG(1)` | inout | State flag. Must be `-1` on entry (set by Y12MB). Set to `-2` on successful exit. |
| `IFLAG(2)` | in | Markowitz search width. Positive integer < N. Recommended: ≤ 3. Ignored when `IFLAG(3)=0`. Unchanged on exit. |
| `IFLAG(3)` | in | Pivoting strategy: 0 = none, 1 = general, 2 = diagonal only. Unchanged on exit. |
| `IFLAG(4)` | in | Structure reuse: 0 = no reuse; 1 = first of same-structure sequence; 2 = subsequent. Unchanged on exit. |
| `IFLAG(5)` | in | L storage: 1 = discard L; 2 = keep L. (Value 3 is invalid here; triggers `IFAIL=22`.) Unchanged on exit. |
| `IFLAG(6)` | out | Row garbage collections. If large, increase `NN`. |
| `IFLAG(7)` | out | Column garbage collections. If large, increase `NN1`. |
| `IFLAG(8)` | out | Peak non-zero count in `A` at any elimination step. |
| `IFLAG(9)` | out | Minimum `NN1` for no column garbage (set when `IFLAG(4)=1`). |
| `IFLAG(10)` | out | Minimum `NN` for no row garbage (set when `IFLAG(4)=1`). |

### Notes

- Do **not** alter `N`, `A`, `SNR`, `NN`, `B`, `PIVOT`, columns 1, 2, 3, 7, 8
  of `HA`, `IHA`, `AFLAG`, `IFLAG(1)`, `IFLAG(3)`, `IFLAG(4)`, or `IFAIL`
  between calls to Y12MC and Y12MD.
- If `IFAIL > 0` on exit do not call Y12MD.

---

## Y12MD

### Purpose

Solves **Ax = b** using the LU factorization computed by Y12MC.

### Signatures

```fortran
! Single precision
SUBROUTINE Y12MDE(N, A, NN, B, PIVOT, SNR, HA, IHA, IFLAG, IFAIL)
REAL    :: A(NN), PIVOT(N), B(N)
INTEGER :: SNR(NN), HA(IHA,11), IFLAG(10)
INTEGER :: N, NN, IHA, IFAIL

! Double precision
SUBROUTINE Y12MDF(N, A, NN, B, PIVOT, SNR, HA, IHA, IFLAG, IFAIL)
DOUBLE PRECISION :: A(NN), PIVOT(N), B(N)
INTEGER          :: SNR(NN), HA(IHA,11), IFLAG(10)
INTEGER          :: N, NN, IHA, IFAIL
```

### Parameters

| Name | Type | Intent | Description |
|------|------|--------|-------------|
| `N` | INTEGER | in | Number of equations. Unchanged on exit. |
| `A` | REAL/DP array(`NN`) | in | Non-zeros of **U** (and optionally **L**) as prepared by Y12MC. Not modified. |
| `NN` | INTEGER | in | Length of `A` and `SNR`. Not modified. |
| `B` | REAL/DP array(`N`) | inout | On entry (when `IFLAG(5) ≠ 3`): vector **c = L⁻¹ P b** from Y12MC. On entry (when `IFLAG(5) = 3`): new right-hand side **b** (LU already available). On exit: solution **x**. |
| `PIVOT` | REAL/DP array(`N`) | in | Diagonal of **U** (from Y12MC). Not modified. |
| `SNR` | INTEGER array(`NN`) | in | Column-number array (from Y12MC). Not modified. |
| `HA` | INTEGER array(`IHA`,11) | in | Work array (from Y12MC). Not modified. |
| `IHA` | INTEGER | in | First dimension of `HA`. Restriction: `IHA ≥ N`. |
| `IFLAG` | INTEGER array(10) | in | `IFLAG(1)` must be `-2`. Use same values of `IFLAG(3)` and `IFLAG(4)` as in Y12MC. `IFLAG(5)`: 1 or 2 if Y12MC was called for this system; 3 if reusing an existing LU (then only supply new **b**). Not modified. |
| `IFAIL` | INTEGER | out | Error indicator. 0 = success. See [Error Codes](#error-codes). |

### Notes

- `IFLAG(1)` is `-2` both on entry and on successful exit (unchanged).
- Do **not** call Y12MB or Y12MC between calls to Y12MD for successive
  right-hand sides when reusing the same LU (`IFLAG(5) = 3`).

---

## Y12MF

### Purpose

Black-box driver with **iterative refinement**.  Factorizes and solves
**Ax = b**, then iteratively improves the solution by computing successive
corrections **d(k) = Q U⁻¹ L⁻¹ P r(k-1)** where **r(k-1) = b − A x(k-1)**.
Only a single-precision version (`Y12MFE`) is available.

### Signature

```fortran
SUBROUTINE Y12MFE(N, A, SNR, NN, RNR, NN1, A1, SN, NZ, HA, IHA, B, B1, X, Y, AFLAG, IFLAG, IFAIL)
REAL    :: A(NN), B(N), B1(N), X(N), Y(N), A1(NZ), AFLAG(11)
INTEGER :: SNR(NN), RNR(NN1), HA(IHA,13), SN(NZ), IFLAG(12)
INTEGER :: N, NN, NN1, NZ, IHA, IFAIL
```

> **Note:** The inner products in the iterative-refinement loop are accumulated
> in double precision and rounded to single precision.  If double-precision
> arithmetic is significantly more expensive than single, consider using the
> plain Y12MA (single) driver.

### Parameters

| Name | Type | Intent | Description |
|------|------|--------|-------------|
| `N` | INTEGER | in | Number of equations. Unchanged on exit. |
| `A` | REAL array(`NN`) | inout | When factorizing (`IFLAG(5)=2`): non-zeros of **A** in arbitrary order in first `NZ` locations; modified by Y12MF. When reusing LU (`IFLAG(5)=3`): non-zeros of **U** and **L**; unchanged on exit. |
| `SNR` | INTEGER array(`NN`) | inout | Column numbers corresponding to `A`. Modified when factorizing; unchanged when reusing LU. |
| `NN` | INTEGER | in | Length of `A` and `SNR`. Restriction: `NN ≥ 2*NZ`. |
| `RNR` | INTEGER array(`NN1`) | inout | Row numbers of non-zeros (ignored when reusing LU). |
| `NN1` | INTEGER | in | Length of `RNR`. Restriction: `NN1 ≥ NZ`. Recommended: `1.5*NZ ≤ NN1 ≤ 2*NZ`. |
| `A1` | REAL array(`NZ`) | inout | Copy of the original non-zeros (ordered by row) for residual computation. Set by Y12MF on first factorization; unchanged when reusing LU. Do not alter between successive calls. |
| `SN` | INTEGER array(`NZ`) | inout | Column numbers of original non-zeros (ordered by row). Set by Y12MF on first factorization; unchanged when reusing LU. |
| `NZ` | INTEGER | in | Number of non-zeros in **A**. Unchanged on exit. |
| `HA` | INTEGER array(`IHA`,13) | inout | Work array (columns 1–11 as in Y12MC; columns 12–13 hold row-start/end positions in `A1`/`SN`). |
| `IHA` | INTEGER | in | First dimension of `HA`. Restriction: `IHA ≥ N`. |
| `B` | REAL array(`N`) | inout | On entry: right-hand side **b**. On exit: last correction vector **d(p-1)**. |
| `B1` | REAL array(`N`) | out | On exit: original right-hand side **b** (saved internally). |
| `X` | REAL array(`N`) | out | On exit: corrected solution vector **x(p)**. |
| `Y` | REAL array(`N`) | out | On exit: diagonal elements of **U** (pivots). Small values indicate singularity. |
| `AFLAG` | REAL array(11) | inout | See [AFLAG table](#aflag-y12mf). |
| `IFLAG` | INTEGER array(12) | inout | See [IFLAG table](#iflag-y12mf). |
| `IFAIL` | INTEGER | out | Error indicator. 0 = success. See [Error Codes](#error-codes). |

#### AFLAG — Y12MF {#aflag-y12mf}

| Index | In/Out | Description |
|-------|--------|-------------|
| `AFLAG(1)` | in | Stability factor (> 1.0; else set to 1.0005). Recommended: 4.0–16.0. |
| `AFLAG(2)` | in | Drop-tolerance. Recommended: in interval `(a×10⁻⁴, a×10⁻¹)` where `a` is the order of magnitude of the elements. If negative, Y12MF computes the row-minimum `a` and uses `|AFLAG(2)|×a` as the tolerance. |
| `AFLAG(3)` | in | Growth-factor threshold (must be ≥ 1.0×10⁵). Recommended: 1.0×10¹⁶. |
| `AFLAG(4)` | in | Singularity threshold multiplier (non-negative). Recommended: 1.0×10⁻¹². |
| `AFLAG(5)` | out | Growth factor. |
| `AFLAG(6)` | out | Maximum element of original **A**. |
| `AFLAG(7)` | out | Maximum element found during elimination. |
| `AFLAG(8)` | out | Minimum pivot. Small values indicate singularity. |
| `AFLAG(9)` | out | Max-norm of last correction vector **d(p-1)**. Should be small for accurate solutions. |
| `AFLAG(10)` | out | Max-norm of last residual vector **r(p-1)**. Should be small for accurate solutions. |
| `AFLAG(11)` | out | Max-norm of corrected solution **x(p)**. The ratio `AFLAG(9)/AFLAG(11)` gives an estimate of the relative error. |

#### IFLAG — Y12MF {#aflag-y12mf}

| Index | In/Out | Description |
|-------|--------|-------------|
| `IFLAG(1)` | work | Used as work space by Y12MF. |
| `IFLAG(2)` | in | Markowitz search width. Positive integer. Recommended: ≤ 3. Ignored when `IFLAG(3)=0`. Unchanged on exit. |
| `IFLAG(3)` | in | Pivoting strategy: 0 = none, 1 = general, 2 = diagonal only. Unchanged on exit. |
| `IFLAG(4)` | in | Structure reuse: 0, 1, or 2 (same semantics as Y12MC). Unchanged on exit. |
| `IFLAG(5)` | in | 2 = factorize; 3 = reuse existing LU. (Value 1 triggers `IFAIL=10`.) Unchanged on exit. |
| `IFLAG(6)` | out | Row garbage collections. If large, increase `NN`. |
| `IFLAG(7)` | out | Column garbage collections. If large, increase `NN1`. |
| `IFLAG(8)` | out | Peak non-zero count. |
| `IFLAG(9)` | out | Minimum `NN1` for no column garbage (set when `IFLAG(4)=1`). |
| `IFLAG(10)` | out | Minimum `NN` for no row garbage (set when `IFLAG(4)=1`). |
| `IFLAG(11)` | in | Maximum number of refinement iterations. Restriction: `IFLAG(11) ≥ 2`. Recommended: < 33. Unchanged on exit. |
| `IFLAG(12)` | out | Number of iterations actually performed. |

### Auxiliary subroutines

Calls Y12MB, Y12MC, and Y12MD.

---

## Y12MH

### Purpose

Computes the **one-norm** of a sparse matrix **A** for use as `ANORM` in Y12MG.
Must be called **before** Y12MC (which overwrites `A`).

### Signatures

```fortran
! Single precision
SUBROUTINE Y12MHE(N, NZ, A, SNR, WORK, ANORM)
REAL    :: A(NZ), WORK(N), ANORM
INTEGER :: SNR(NZ), N, NZ

! Double precision
SUBROUTINE Y12MHF(N, NZ, A, SNR, WORK, ANORM)
DOUBLE PRECISION :: A(NZ), WORK(N), ANORM
INTEGER          :: SNR(NZ), N, NZ
```

### Parameters

| Name | Type | Intent | Description |
|------|------|--------|-------------|
| `N` | INTEGER | in | Number of equations (matrix dimension). |
| `NZ` | INTEGER | in | Number of non-zeros in **A**. |
| `A` | REAL/DP array(`NZ`) | in | Non-zero values of **A**. Not modified. |
| `SNR` | INTEGER array(`NZ`) | in | Column indices of `A`. Not modified. |
| `WORK` | REAL/DP array(`N`) | work | Work space; contents are not meaningful on exit. |
| `ANORM` | REAL/DP | out | One-norm (maximum column-sum) of **A**. |

> **Note:** Y12MH has no `IFAIL` parameter and performs no error checking.

---

## Y12MG

### Purpose

Estimates the **reciprocal condition number** of **A** using the LINPACK
algorithm (Cline, Moler, Stewart, Wilkinson, 1979).  Must be called **after**
Y12MC, while the LU factorization is still intact.

### Signatures

```fortran
! Single precision
SUBROUTINE Y12MGE(N, NN, A, SNR, W, PIVOT, ANORM, RCOND, IHA, HA, IFLAG, IFAIL)
REAL    :: A(NN), W(N), PIVOT(N), RCOND, ANORM
INTEGER :: SNR(NN), HA(IHA,3), IFLAG(5)
INTEGER :: N, NN, IHA, IFAIL

! Double precision
SUBROUTINE Y12MGF(N, NN, A, SNR, W, PIVOT, ANORM, RCOND, IHA, HA, IFLAG, IFAIL)
DOUBLE PRECISION :: A(NN), W(N), PIVOT(N), RCOND, ANORM
INTEGER          :: SNR(NN), HA(IHA,3), IFLAG(5)
INTEGER          :: N, NN, IHA, IFAIL
```

### Parameters

| Name | Type | Intent | Description |
|------|------|--------|-------------|
| `N` | INTEGER | in | Number of equations. |
| `NN` | INTEGER | in | Length of `A` and `SNR`. |
| `A` | REAL/DP array(`NN`) | in | Non-zeros of **U** and **L** as left by Y12MC. |
| `SNR` | INTEGER array(`NN`) | in | Column indices (from Y12MC). |
| `W` | REAL/DP array(`N`) | work | Work array; overwritten during computation. |
| `PIVOT` | REAL/DP array(`N`) | in | Diagonal of **U** (from Y12MC). |
| `ANORM` | REAL/DP | in | One-norm of original **A** (from Y12MH). |
| `RCOND` | REAL/DP | out | Reciprocal condition number estimate. The actual condition number is approximately `1/RCOND`. On error, set to −1. |
| `IHA` | INTEGER | in | First dimension of `HA`. |
| `HA` | INTEGER array(`IHA`,3) | in | Work array (uses first 3 columns from Y12MC output). |
| `IFLAG` | INTEGER array(5) | in | Only `IFLAG(5)` is checked: if `IFLAG(5) = 1` (L was discarded), sets `IFAIL = 26`. |
| `IFAIL` | INTEGER | inout | On entry: checked; if already non-zero, returns immediately. On exit: 0 = success; 26 = L was discarded (call with `IFLAG(5) ≥ 2`). |

---

## Error Codes

`IFAIL` is the error-diagnostics parameter.  On exit from each subroutine,
`IFAIL = 0` indicates success.  A positive value means an error was detected
and computation was halted.

| `IFAIL` | Cause |
|---------|-------|
| 1 | Matrix was not factorized: Y12MD called without a preceding Y12MC for the current system. Ensure `IFLAG(1) ≥ 0` before the first package call. |
| 2 | Matrix was not ordered: Y12MC called without a preceding Y12MB. |
| 3 | A pivotal element with `|a(i,i)| < AFLAG(4)*AFLAG(6)` was selected. Likely numerical singularity when `AFLAG(4)` is small. |
| 4 | Growth factor `AFLAG(5) > AFLAG(3)`. Computation stopped to prevent overflow. A smaller stability factor `AFLAG(1)` may help. |
| 5 | `NN` (length of `A` and `SNR`) is too small. Increase `NN` (and possibly `NN1`). |
| 6 | `NN1` (length of `RNR`) is too small. Increase `NN1` (and possibly `NN`). |
| 7 | A row without non-zeros in its active part was found during factorization. Likely numerical singularity (especially when `AFLAG(2)` is small). |
| 8 | A column without non-zeros in its active part was found during factorization. Same interpretation as `IFAIL=7`. |
| 9 | A pivotal element is missing. Occurs with `AFLAG(2) > 0` and `IFLAG(4) = 2`, or with special pivoting strategies (`IFLAG(3) = 0` or `2`). Decrease drop-tolerance or refactorize. |
| 10 | Y12MF called with `IFLAG(5) = 1`. Use `IFLAG(5) = 2` instead. |
| 11 | At least two non-zeros occupy the same position (i, j). Check input data. |
| 12 | `N < 2`. Check the value of `N`. |
| 13 | `Z ≤ 0` (renamed `NZ` in Y12MF). Check the number of non-zeros. |
| 14 | `Z < N`. If `Z` is correct, the matrix is structurally singular. |
| 15 | `IHA < N`. Increase the first dimension of `HA`. |
| 16 | `IFLAG(4)` is not 0, 1, or 2. |
| 17 | A row of the coefficient matrix contains no non-zeros. Structurally singular. |
| 18 | A column of the coefficient matrix contains no non-zeros. Structurally singular. |
| 19 | `IFLAG(2) < 1`. Set `IFLAG(2)` to a positive integer (recommended: 3). |
| 20 | `IFLAG(3)` is not 0, 1, or 2. |
| 21 | `IFLAG(5)` is not 1, 2, or 3. |
| 22 | Y12MB or Y12MC called with `IFLAG(5) = 3`. These subroutines require `IFLAG(5) = 1` or `2`. |
| 23 | `IFLAG(11) < 2` (Y12MF). Increase the maximum iteration count. |
| 24 | A column index is outside the range [1, N]. |
| 25 | A row index is outside the range [1, N]. |
| 26 | Y12MG called when `IFLAG(5) = 1`, i.e., matrix **L** was discarded. Call with `IFLAG(5) ≥ 2` to retain **L**. |

---

## Discrepancies between this documentation and the implementation

The following issues were found by cross-checking the original `doc` file
against the source files in `src/legacy/`.  They are flagged here for human
review.

1. **`AFLAG(3)` default in Y12MA (doc vs. code).**
   The `doc` file states *"Y12MA sets `AFLAG(3)` = 1.0E6"*, but both
   `y12mae.f` and `y12maf.f` set `aflag(3) = 1.0e+16` / `1.0d+16`.
   This documentation uses the value found in the code (1.0 × 10¹⁶).

2. **`IFLAG(2)` default in Y12MA (doc vs. code).**
   The `doc` file states *"Subroutine Y12MA sets `IFLAG(2)` = 3"*, but the
   code in `y12mae.f` / `y12maf.f` sets `iflag(2) = 2`.
   This documentation uses the value found in the code (2).

3. **`IFAIL = 26` not in the original error-diagnostics section.**
   Y12MG can return `IFAIL = 26` (L was discarded), but this code does not
   appear in the error-diagnostics section of the original `doc` file.
   It is documented here based on the source in `y12mg.f`.

4. **Y12MG and Y12MH not in the original `doc` file.**
   The `doc` file documents only Y12MA, Y12MB, Y12MC, Y12MD, and Y12MF.
   Y12MG and Y12MH are present in `src/legacy/y12mg.f` and
   `src/legacy/y12mh.f` but have no corresponding section in `doc`.
   They are documented here from the source-code comments.

5. **Y12MG signature quirk: `IHA` precedes `HA`.**
   In Y12MG the argument order is `…, ANORM, RCOND, IHA, HA, IFLAG, IFAIL`,
   which places `IHA` before `HA`.  This is the reverse of the ordering used
   in Y12MB, Y12MC, and Y12MD.

6. **Y12MG uses only first 3 columns of `HA` and only 5 elements of `IFLAG`.**
   The array dimensions in `y12mg.f` are `HA(IHA,3)` and `IFLAG(5)`, whereas
   all other subroutines declare `HA(IHA,11)` and `IFLAG(10)`.  Callers must
   ensure the arrays are compatible.

7. **`IFLAG(2)` description in Y12MF (typo in original doc).**
   The original `doc` text for Y12MF's `IFLAG(2)` contains a copy-paste error:
   it says "If `IFLAG((2) .ge. 0` then the pivotal search…" (parenthesis
   mismatch), and the intended condition should likely be `IFLAG(3) ≥ 0` or
   `IFLAG(2) ≥ 1`, consistent with Y12MC's description.

---

## References

1. Zlatev, Z., Wasniewski, J., & Schaumburg, K. (1981). *Y12M: solution of large and sparse systems of linear algebraic equations*. Lecture Notes in Computer Science, Vol. 121. Springer. https://doi.org/10.1007/3-540-10874-2
2. Cline, A.K., Moler, C.B., Stewart, G.W., and Wilkinson, J.H. (1979). "An estimate for the condition number of a matrix." *SIAM J. Numer. Anal.* 16, 368–375.
3. Dongarra, J.J., Bunch, J.R., Moler, C.B., and Stewart, G.W. (1979). *LINPACK User's Guide*. SIAM, Philadelphia.
4. Gustavson, F.G. (1972). "Some basic techniques for solving sparse systems of linear equations." In: *Sparse Matrices and Their Applications* (Rose & Willoughby, eds.), pp. 41–52. Plenum Press, New York.
