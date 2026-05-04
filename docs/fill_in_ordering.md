# Fill-in reduction ordering in Y12M

## Background

Y12M applies a **greedy local Markowitz strategy** during factorization
(`Y12MC`): at each elimination step it selects the pivot that minimizes the
product `(nnz_in_row − 1) × (nnz_in_col − 1)`.  This limits fill-in locally
but is not a substitute for a global, pre-factorization fill-in reducing
ordering such as AMD, COLAMD, or RCM.

`Y12MB` builds the internal CSR/CSC data structures and the Markowitz linked
list required by `Y12MC`.  It does **not** compute any global ordering
permutation.

## Applying an external permutation yourself

Because Y12M accepts input in coordinate (COO) format — value array `A`, column
index array `SNR`, and row index array `RNR` — it is straightforward to apply a
fill-reducing permutation before calling `Y12MB`.

### Step 1 — Compute a permutation

Use an external library to compute a symmetric fill-reducing permutation
`perm(1:n)` of the row/column indices.  Common choices:

- [AMD](https://people.engr.tamu.edu/davis/suitesparse.html) — Approximate
  Minimum Degree (part of SuiteSparse)
- [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) — multilevel
  graph partitioning and reordering
- Reverse Cuthill–McKee (RCM) — available in
  [Sparsekit](https://people.sc.fsu.edu/~jburkardt/f_src/sparsekit/sparsekit.html)
  and many other sparse toolkits

### Step 2 — Renumber the index arrays

Renumber the stored row and column indices.  The matrix **values** (`A`) do not
change — only the indices need updating.  This can be done in-place:

```fortran
do i = 1, z
    snr(i) = perm(snr(i))   ! column indices
    rnr(i) = perm(rnr(i))   ! row indices
end do
```

Each iteration reads `perm` at an indirect (gathered) address and writes to a
distinct element of `snr`/`rnr`.  There are no loop-carried dependencies, so
the loop is safe to vectorize with a gather instruction.

### Step 3 — Permute the right-hand side

The right-hand side vector must be permuted to form **b' = P b**.  Because
`perm` is a permutation (all values are distinct), the scatter writes go to
distinct output locations and the loop is safe to vectorize with a scatter
instruction.  A temporary is still required because `b_perm` and `b` may
alias.  The scatter form is:

```fortran
real :: b_perm(n)
do i = 1, n
    b_perm(perm(i)) = b(i)   ! scatter: P b
end do
```

An equivalent gather form, which is sometimes preferable for vectorization, is:

```fortran
real :: b_perm(n)
do i = 1, n
    b_perm(i) = b(iperm(i))  ! gather: P b using inverse permutation
end do
```

### Step 4 — Solve

Call `Y12MB`, `Y12MC`, and `Y12MD` (or `Y12MA`) using `b_perm` in place of the
original right-hand side.

### Step 5 — Recover the solution

After `Y12MD` returns, map the solution back to the original ordering:

```fortran
do i = 1, n
    x(i) = b_perm(perm(i))
end do
```

## Notes

### Permutation arrays and permutation matrices

Let **P** be the *n × n* permutation matrix corresponding to the reordering.
In Fortran the permutation is stored as an integer array `perm(1:n)` where
`perm(i) = j` means row/column `i` of the original matrix maps to position `j`
in the permuted matrix; equivalently, `P(j, i) = 1`.

The inverse permutation `iperm` satisfies `iperm(perm(i)) = i` for all `i`,
and corresponds to the matrix **Pᵀ = P⁻¹**:

```fortran
integer :: iperm(n)
do i = 1, n
    iperm(perm(i)) = i
end do
```

The overall transformation is the symmetric (congruence) permutation
**A' = P A Pᵀ**.  The right-hand side is permuted as **b' = P b**
(Steps 2–3 above), and the solution recovered via **x = Pᵀ x'** (Step 5).

### Numerical stability

Y12M's factorization (`Y12MC`) selects pivots using a **local** Markowitz
criterion: at each elimination step it chooses the entry that minimizes the
product `(nnz_in_row − 1) × (nnz_in_col − 1)`, subject to a numerical
threshold controlled by `AFLAG(1)` (the stability parameter *u*).  A
fill-reducing pre-ordering changes the sparsity pattern that `Y12MC` sees,
which in turn changes which candidates pass the threshold test and are
available for selection.

In the worst case, a permutation that strongly concentrates non-zeros can
reduce the pool of numerically acceptable pivots and degrade stability.
In practice, for well-scaled matrices the effect is small, because `Y12MC`
already incorporates a threshold criterion.  Users who observe poor backward
error after applying a fill-reducing ordering should consider relaxing the
threshold (lowering `AFLAG(1)`) or using a stability-aware ordering such as
MC64 (from HSL) that combines sparsity and numerical considerations.

