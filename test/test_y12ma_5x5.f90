!> @file test_y12ma_5x5.f90
!> @brief Regression test for the Y12M high-level driver `y12ma` on a 5x5 sparse system.
!>
!> This test solves a fixed 5x5 sparse linear system using the high-level routine `y12ma`,
!> then checks the solution by computing the residual r = b - A*x using the ORIGINAL input
!> triplets (A values and their row/column indices), because `y12ma` overwrites its input
!> arrays during factorization/solve.
!>
!> PASS/FAIL is determined by a relative residual check:
!>   relerr = ||r||_2 / max(||b||_2, sqrt(tiny(1.0)))
!>
!> Matrix provenance:
!>   The matrix and RHS in this test were taken from the “UMFPACK” Fortran examples page
!>   (section “C original form / UMFPACK”):
!>   https://geo.mff.cuni.cz/~lh/Fortran/UMFPACK/#c-original-form-umfpack
program test_y12ma_5x5
  use y12m, only: y12ma
  implicit none

  integer, parameter :: n   = 5
  integer, parameter :: nz  = 12
  integer, parameter :: nn  = 50
  integer, parameter :: nn1 = 50
  integer, parameter :: iha = n

  real    :: a(nn), b(n), pivot(n)
  integer :: snr(nn), rnr(nn1)
  integer :: ha(iha, 11), iflag(10)
  real    :: aflag(8)

  real    :: a0(nz), b0(n)
  integer :: snr0(nz), rnr0(nz)

  real    :: resid(n)
  real    :: normrinf, normbinf, denominf, relinf
  real    :: rtol2, rtoli
  integer :: i, ifail

  ! Workspace constraints from doc
  if (nn  < 2*nz) error stop "FAIL: invalid workspace (NN must be >= 2*NZ)"
  if (nn1 < nz)   error stop "FAIL: invalid workspace (NN1 must be >= NZ)"
  if (iha < n)    error stop "FAIL: invalid workspace (IHA must be >= N)"

  ! --- Define A in triplet form ---
  a0   = [ 2.0,  3.0,  3.0,  4.0,  6.0, -1.0, -3.0,  2.0,  1.0,  4.0,  2.0,  1.0]

  ! Per doc:
  !   SNR(j) = column number of A(j)
  !   RNR(j) = row number of A(j)
  rnr0 = [ 1, 1, 2, 2, 2, 3, 3, 3, 4, 5, 5, 5 ]   ! rows
  snr0 = [ 1, 2, 1, 3, 5, 2, 3, 4, 3, 2, 3, 5 ]   ! cols

  b0 = [ 8.0, 45.0, -3.0, 3.0, 19.0 ]

  ! Copy into working arrays (y12ma overwrites these)
  a(1:nz)   = a0
  snr(1:nz) = snr0
  rnr(1:nz) = rnr0
  b         = b0

  ! --- Solve (factorize + solve) ---
  call y12ma(n, nz, a, snr, nn, rnr, nn1, pivot, ha, iha, aflag, iflag, b, ifail)
  if (ifail /= 0) then
    print *, "FAIL: y12ma returned IFAIL=", ifail
    error stop 1
  end if

  ! --- Residual r = b0 - A0*x, where x is now stored in b ---
  resid = b0
  do i = 1, nz
    resid(rnr0(i)) = resid(rnr0(i)) - a0(i) * b(snr0(i))
  end do

  ! --- Print solution and residual ---
  print *, "solution x (stored in b):"
  do i = 1, n
    print "(1x,g0.7)", b(i)
  end do

  print *, "residual r = b - A*x:"
  do i = 1, n
    print "(1x,g0.7)", resid(i)
  end do

  ! --- Relative tolerance tests (2-norm and infinity norm) ---
  normr2 = norm2(resid)
  normb2 = norm2(b0)

  denom  = max(normb2, sqrt(tiny(1.0)))
  relerr = normr2 / denom

  normrinf = maxval(abs(resid))
  normbinf = maxval(abs(b0))

  denominf = max(normbinf, tiny(1.0))
  relinf   = normrinf / denominf

  rtol2  = 1.0e-5
  rtoli  = 1.0e-5

  if (relerr <= rtol2 .and. relinf <= rtoli) then
    print "(a,1x,g0.7,1x,a,1x,g0.7)", "PASS: rel2=", relerr, "<= rtol2=", rtol2
    print "(a,1x,g0.7,1x,a,1x,g0.7)", "      relinf=", relinf, "<= rtoli=", rtoli
  else
    print "(a,1x,g0.7,1x,a,1x,g0.7)", "FAIL: rel2=", relerr, ">  rtol2=", rtol2
    print "(a,1x,g0.7,1x,a,1x,g0.7)", "      relinf=", relinf, ">  rtoli=", rtoli
    error stop 2
  end if
end program test_y12ma_5x5