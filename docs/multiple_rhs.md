# Usage Scenarios

Y12M is designed for five distinct solving scenarios, each requiring a
different combination of `IFLAG(4)` and `IFLAG(5)` settings.

## Case (i) — One system with a single right-hand side

Use `Y12MA` (the black-box driver) or call `Y12MB` + `Y12MC` + `Y12MD` with
`IFLAG(4) = 0` and `IFLAG(5) = 1`.  After the solve, the matrix factors are
discarded and no state is kept.

```
IFLAG(4) = 0   ! no structural reuse
IFLAG(5) = 1   ! discard L after factorization
call Y12MB(...)
call Y12MC(...)
call Y12MD(...)
```

## Case (ii) — Several systems with the same coefficient matrix

When the coefficient matrix **A** does not change between solves, the LU
factorization only needs to be computed once.  Set `IFLAG(5) = 2` to retain
the **L** factor, then reuse it by calling only `Y12MD` with `IFLAG(5) = 3`
for every subsequent right-hand side.

`IFLAG(1)` must not be changed between calls; it remains `-2` after a
successful factorization.

```
! First solve
IFLAG(4) = 0   ! or 1 if the structure may be reused later
IFLAG(5) = 2   ! keep L after factorization
call Y12MB(...)
call Y12MC(...)
call Y12MD(...)

! Subsequent solves — new right-hand side only
IFLAG(5) = 3   ! signal: LU already available
B = <new right-hand side>
call Y12MD(...)
```

> **Note:** Do *not* call `Y12MB` or `Y12MC` between right-hand sides; doing
> so would destroy the factorization.

## Case (iii) — Several systems with the same sparsity structure

When successive matrices share the same set of non-zero positions (same
`SNR`/`RNR` pattern) but carry different numerical values, the ordering
information computed on the first call to `Y12MB` can be reused.

Set `IFLAG(4) = 1` for the **first** system in the sequence, and `IFLAG(4) = 2`
for every **subsequent** system with the same structure.

```
! First matrix
IFLAG(4) = 1   ! compute and store structural ordering
IFLAG(5) = 1   ! or 2 if multiple RHS are expected
call Y12MB(...)
call Y12MC(...)
call Y12MD(...)

! Second matrix — same non-zero positions, new numerical values
A   = <new non-zero values>   ! positions unchanged
IFLAG(4) = 2                  ! reuse structural ordering
call Y12MB(...)
call Y12MC(...)
call Y12MD(...)
```

## Case (iv) — Combined: same structure, same coefficient matrix appears successively

This is the most general reuse scenario: a sequence of matrices sharing the
same sparsity structure, some of which are numerically identical.  Use
`IFLAG(4) = 1` for the first matrix, `IFLAG(4) = 2` for subsequent ones, and
`IFLAG(5) = 2` / `IFLAG(5) = 3` to avoid re-factorizing when the numerical
values have not changed.

```
! First matrix
IFLAG(4) = 1
IFLAG(5) = 2   ! keep L for potential RHS reuse
call Y12MB(...)
call Y12MC(...)
call Y12MD(...)

! Same matrix, new right-hand side
IFLAG(5) = 3
B = <new right-hand side>
call Y12MD(...)

! New matrix, same sparsity structure
A = <new non-zero values>
IFLAG(4) = 2
IFLAG(5) = 2   ! keep L again
call Y12MB(...)
call Y12MC(...)
call Y12MD(...)
```

## Case (v) — Several systems with different coefficient matrices of different structure

Each system must be solved independently; no reuse is possible.  Call `Y12MA`
(or `Y12MB` + `Y12MC` + `Y12MD` with `IFLAG(4) = 0`) for every system in the
sequence.

```
do i = 1, num_systems
    ! load A, SNR, RNR, B for system i
    IFLAG(4) = 0
    IFLAG(5) = 1
    call Y12MB(...)
    call Y12MC(...)
    call Y12MD(...)
end do
```

## Condition number estimation (any case)

To estimate the condition number during any factorization:

1. Call `Y12MH` **before `Y12MC`** to compute the one-norm of the original
   matrix into `ANORM` (Y12MC overwrites array `A`).
2. Call `Y12MG` **after `Y12MC`**, passing `ANORM`, to obtain the reciprocal
   condition number `RCOND`.  The actual condition number is `1.0 / RCOND`.

Both subroutines are optional and independent of the reuse patterns described
above.
