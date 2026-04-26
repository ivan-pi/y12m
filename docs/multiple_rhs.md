# Solving Systems with Multiple Right-Hand Sides

Y12M supports several reuse patterns that avoid redundant work when solving
more than one system.

## Case 1 – Multiple right-hand sides, same matrix

When the coefficient matrix **A** does not change between solves, the LU
factorization only needs to be computed once.  Subsequent right-hand sides are
solved by calling `Y12MD` alone.

**Required setup:** call `Y12MB` and `Y12MC` with `IFLAG(5) = 2` so that the
non-zero elements of the lower triangular factor **L** are retained after
factorization.

**Subsequent solves:** set `IFLAG(5) = 3`, load the new right-hand side into
array `B`, and call `Y12MD`.  `IFLAG(1)` must not be changed between calls;
it remains `-2` after a successful factorization.

```
! First solve
IFLAG(5) = 2                         ! keep L after factorization
call Y12MB(...)
call Y12MC(...)                      ! computes LU, modifies B -> L^{-1}Pb
call Y12MD(...)                      ! back-solves, x in B

! Subsequent solves (new B only)
IFLAG(5) = 3                         ! signal: LU already available
B = <new right-hand side>
call Y12MD(...)                      ! back-solves, x in B
```

> **Note:** Do *not* call `Y12MB` or `Y12MC` between right-hand sides; doing
> so would destroy the factorization.

## Case 2 – Multiple matrices with the same sparsity structure

When successive matrices share the same set of non-zero positions (same
`SNR`/`RNR` sparsity pattern) but different numerical values, the reordering
step (`Y12MB`) was already optimised for that structure.  Each new matrix still
requires a full factorization, but pivoting can reuse the column ordering found
for the first matrix.

**First system:** call with `IFLAG(4) = 1`.

**Subsequent systems (same structure):** reload the new values into array `A`,
set `IFLAG(4) = 2`, and call `Y12MB` + `Y12MC` + `Y12MD` as usual.

```
! First matrix
IFLAG(4) = 1
IFLAG(5) = 1                         ! or 2 if multiple RHS are expected
call Y12MB(...)
call Y12MC(...)
call Y12MD(...)

! Second matrix – same sparsity structure, new values
A = <new non-zero values>            ! same positions as before
IFLAG(4) = 2                         ! reuse structural information
call Y12MB(...)
call Y12MC(...)
call Y12MD(...)
```

## Combining both cases

Cases 1 and 2 can be combined: use `IFLAG(4) = 2` together with `IFLAG(5) = 2`
to retain **L** for multiple right-hand sides while still processing successive
same-structure matrices efficiently.

## Condition number estimation

To estimate the condition number at any factorization:

1. Call `Y12MH` **before `Y12MC`** to store the one-norm of the original matrix
   in `ANORM` (Y12MC overwrites array `A`).
2. Call `Y12MG` **after `Y12MC`**, passing `ANORM`, to obtain the reciprocal
   condition number `RCOND`.  The actual condition number is `1.0 / RCOND`.

Both subroutines are optional and independent of the right-hand-side reuse
patterns described above.
