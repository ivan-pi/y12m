# y12m

Solution of Large and Sparse Systems of Linear Algebraic Equations

**Original version at Netlib:** http://www.netlib.org/y12m/

**Book:**

> Zlatev, Z., Wasniewski, J., & Schaumburg, K. (1981). Y12M: solution of large and sparse systems of linear algebraic equations (Vol. 121). Berlin, Heidelberg, New York: Springer. https://doi.org/10.1007/3-540-10874-2

**Home page of author Zahari Zlatev:** https://www.dmu.dk/atmosphericenvironment/staff/zlatev.htm

## Calling Sequence

The Y12M package provides subroutines at two levels. Every subroutine is available in a single-precision variant (suffix `E`) and a double-precision variant (suffix `F`), for example `Y12MBE` / `Y12MBF`.

### High-level drivers

| Subroutine | Purpose |
|------------|---------|
| `Y12MA` | Black-box driver for a single system with a single right-hand side. Calls `Y12MB`, `Y12MC`, and `Y12MD` internally. |
| `Y12MF` | Factorizes and solves a system in one call with iterative refinement to improve accuracy. |

### Lower-level subroutines

For finer control—solving the same system for multiple right-hand sides, reusing an existing LU factorization, or processing a sequence of matrices that share the same sparsity structure—the lower-level subroutines should be called directly:

| Subroutine | Purpose |
|------------|---------|
| `Y12MH` | Computes the one-norm of matrix A. *(optional)* |
| `Y12MB` | Prepares and reorders the matrix for factorization. |
| `Y12MC` | Computes the LU factorization of the matrix. |
| `Y12MG` | Computes the reciprocal of the condition number. *(optional)* |
| `Y12MD` | Solves the system using the LU factorization. |

### Calling order

```mermaid
flowchart TD
    start(["Start"])

    ymh["Y12MH
    Compute one-norm of A
    (optional)"]

    ymb["Y12MB
    Prepare and reorder matrix"]

    ymc["Y12MC
    LU factorization"]

    ymg["Y12MG
    Reciprocal condition number
    (optional)"]

    ymd["Y12MD
    Solve system Ax = b"]

    done(["Done"])

    start --> ymh --> ymb --> ymc --> ymg --> ymd --> done
    ymd -->|"new right-hand side,
    same matrix (IFLAG(5)=3)"| ymd
    ymc -->|"new matrix,
    same sparsity structure"| ymh

    style ymh fill:#ffffcc,stroke:#999
    style ymg fill:#ffffcc,stroke:#999
```

The two optional subroutines have positional constraints:

- **`Y12MH`** must be called **before `Y12MC`**, because the LU factorization overwrites the matrix values stored in array `A`. It is recommended to call `Y12MH` before `Y12MB`.
- **`Y12MG`** must be called **after `Y12MC`**, while the LU factorization is still intact. It takes the one-norm computed by `Y12MH` as an input argument.

