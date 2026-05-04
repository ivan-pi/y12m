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

### Step 3 — Permute the right-hand side

The right-hand side vector must be permuted before the solve.  This is a
scatter operation and **requires a temporary** array — an in-place scatter can
overwrite values that are still needed:

```fortran
real :: b_perm(n)
do i = 1, n
    b_perm(perm(i)) = b(i)
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

- `Y12MC` and `Y12MD` operate entirely in the permuted index space and require
  no changes.
- If the permutation is symmetric (i.e., the same permutation `perm` is applied
  to both rows and columns, as in `P A Pᵀ`), the solution vector permutation is
  simply the inverse permutation.  For a permutation stored as an index array
  `perm`, the inverse is `iperm` where `iperm(perm(i)) = i`.
