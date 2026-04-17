! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4.5
!
! test_y12mc_z_intent.f90 - Regression test for y12mce/y12mcf z intent violation.
!
! The y12m documentation states that z is "Unchanged on exit" for y12mc, but
! the legacy implementation (y12mce, y12mcf) internally modifies z and only
! attempts to restore it at the end.  The restore is broken on the early-exit
! path: zz=z is saved AFTER the goto that handles initial error conditions, so
! when ifail>0 is set before that save point, z is overwritten with an
! uninitialised local variable.
!
! This test verifies three invariants:
!
!  1. z is unchanged after a successful y12mb->y12mc->y12md sequence.
!  2. z is unchanged after y12mc returns early with an error (ifail=22,
!     triggered by iflag(5)=3 before the zz=z save point is reached).
!  3. y12mc can be called with z as an integer literal constant.  On an
!     unpatched build the implementation writes to the argument storage,
!     which causes a segmentation fault on any system where constant data
!     is mapped read-only.
!
! Before the fix test 2 fails because z is assigned from the uninitialised
! local zz, which is typically 0 on first entry.  Test 3 causes a segfault.
!
program test_y12mc_z_intent
  use y12m
  implicit none

  integer, parameter :: NMAX  = 10
  integer, parameter :: NNP   = 200
  integer, parameter :: NN1P  = 200

  real    :: a(NNP), pivot(NMAX), b(NMAX), aflag(8)
  integer :: snr(NNP), rnr(NN1P), ha(NMAX,11), iflag(10)
  integer :: n, z, z_saved, ifail, nfail
  real    :: err

  nfail = 0

  ! -----------------------------------------------------------------------
  ! Test 1: z must be unchanged after a successful factorization (normal path).
  ! -----------------------------------------------------------------------
  n = 4
  call build_tridiag_sp(n, NNP, NN1P, z, a, snr, rnr, b)
  z_saved = z

  iflag(1) = 0
  iflag(2) = 3
  iflag(3) = 1
  iflag(4) = 0
  iflag(5) = 1
  aflag(1) = 16.0
  aflag(2) = 1.0e-12
  aflag(3) = 1.0e+16
  aflag(4) = 1.0e-12

  call y12mb(n, z, a, snr, NNP, rnr, NN1P, ha, NMAX, aflag, iflag, ifail)
  if (ifail /= 0) then
    write(*,'(a,i0)') 'FAIL test1: y12mb failed with ifail=', ifail
    nfail = nfail + 1
    go to 10
  end if

  call y12mc(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, NMAX, &
      aflag, iflag, ifail)
  if (ifail /= 0) then
    write(*,'(a,i0)') 'FAIL test1: y12mc failed with ifail=', ifail
    nfail = nfail + 1
    go to 10
  end if

  call y12md(n, a, NNP, b, pivot, snr, ha, NMAX, iflag, ifail)

  if (z /= z_saved) then
    write(*,'(a,i0,a,i0)') &
        'FAIL test1: z modified on success path: before=', z_saved, &
        ' after=', z
    nfail = nfail + 1
  else
    write(*,'(a)') 'PASS test1: z unchanged after successful factorization'
  end if

10 continue

  ! -----------------------------------------------------------------------
  ! Test 2: z must be unchanged after y12mc returns with an error before
  ! it has a chance to save the original value of z.
  !
  ! Setting iflag(5)=3 causes ifail=22 to be set before the zz=z save in
  ! y12mce/y12mcf, so the restore "z=zz" at the exit label uses an
  ! uninitialised local variable.  Without the fix z is corrupted (typically
  ! set to 0); with the fix z is never written to.
  ! -----------------------------------------------------------------------
  n = 4
  call build_tridiag_sp(n, NNP, NN1P, z, a, snr, rnr, b)

  ! Prepare the matrix first so iflag(1) is set to -1.
  iflag(1) = 0
  iflag(2) = 3
  iflag(3) = 1
  iflag(4) = 0
  iflag(5) = 1
  aflag(1) = 16.0
  aflag(2) = 1.0e-12
  aflag(3) = 1.0e+16
  aflag(4) = 1.0e-12

  call y12mb(n, z, a, snr, NNP, rnr, NN1P, ha, NMAX, aflag, iflag, ifail)
  if (ifail /= 0) then
    write(*,'(a,i0)') 'FAIL test2: y12mb failed with ifail=', ifail
    nfail = nfail + 1
    go to 20
  end if

  z_saved = z

  ! Trigger the early-exit path: iflag(5)=3 sets ifail=22 before zz=z.
  iflag(5) = 3

  call y12mc(n, z, a, snr, NNP, rnr, NN1P, pivot, b, ha, NMAX, &
      aflag, iflag, ifail)

  if (ifail /= 22) then
    write(*,'(a,i0)') &
        'FAIL test2: expected ifail=22 for iflag(5)=3, got ', ifail
    nfail = nfail + 1
  end if

  if (z /= z_saved) then
    write(*,'(a,i0,a,i0)') &
        'FAIL test2: z modified on error-exit path: before=', z_saved, &
        ' after=', z
    nfail = nfail + 1
  else
    write(*,'(a)') 'PASS test2: z unchanged after early-error exit'
  end if

20 continue

  ! -----------------------------------------------------------------------
  ! Test 3: call y12mb and y12mc with z supplied as an integer literal
  ! constant (5).  The 3-by-3 matrix used here has exactly 5 nonzeros:
  !
  !   1.0  0.0  3.0
  !   0.0 -1.0  4.0
  !   2.0  0.0  0.0
  !
  ! A*[1,1,1]^T = [4, 3, 2], so the expected solution is x = [1,1,1].
  !
  ! y12mbe (called by y12mb) only reads z and never writes to it, so
  ! passing a literal is harmless there.  y12mce (called by y12mc)
  ! modifies z internally; on an unpatched build that write hits the
  ! read-only constant storage and causes a segmentation fault.
  ! -----------------------------------------------------------------------
  n = 3
  a(1) = 1.0  ; snr(1) = 1 ; rnr(1) = 1
  a(2) = 3.0  ; snr(2) = 3 ; rnr(2) = 1
  a(3) = -1.0 ; snr(3) = 2 ; rnr(3) = 2
  a(4) = 4.0  ; snr(4) = 3 ; rnr(4) = 2
  a(5) = 2.0  ; snr(5) = 1 ; rnr(5) = 3
  b(1) = 4.0 ; b(2) = 3.0 ; b(3) = 2.0

  iflag(1) = 0
  iflag(2) = 3
  iflag(3) = 1
  iflag(4) = 0
  iflag(5) = 1
  aflag(1) = 16.0
  aflag(2) = 1.0e-12
  aflag(3) = 1.0e+16
  aflag(4) = 1.0e-12

  call y12mb(n, 5, a, snr, NNP, rnr, NN1P, ha, NMAX, aflag, iflag, ifail)
  if (ifail /= 0) then
    write(*,'(a,i0)') 'FAIL test3: y12mb failed with ifail=', ifail
    nfail = nfail + 1
    go to 30
  end if

  call y12mc(n, 5, a, snr, NNP, rnr, NN1P, pivot, b, ha, NMAX, &
      aflag, iflag, ifail)
  if (ifail /= 0) then
    write(*,'(a,i0)') 'FAIL test3: y12mc failed with ifail=', ifail
    nfail = nfail + 1
    go to 30
  end if

  call y12md(n, a, NNP, b, pivot, snr, ha, NMAX, iflag, ifail)

  err = maxval(abs(b(1:n) - 1.0))
  if (err > 1.0e-4) then
    write(*,'(a,es10.3)') 'FAIL test3: z as literal, max_err=', err
    nfail = nfail + 1
  else
    write(*,'(a,es10.3)') 'PASS test3: z as literal, max_err=', err
  end if

30 continue

  ! -----------------------------------------------------------------------
  ! Summary
  ! -----------------------------------------------------------------------
  if (nfail /= 0) then
    write(*,'(i0,a)') nfail, ' test(s) FAILED'
    stop 1
  end if
  write(*,'(a)') 'All test_y12mc_z_intent tests PASSED'

contains

  ! n-by-n tridiagonal (diag=3, off-diag=-1), solution x=[1,...,1].
  subroutine build_tridiag_sp(n, nnmax, nn1max, z, a, snr, rnr, b)
    integer, intent(in)  :: n, nnmax, nn1max
    integer, intent(out) :: z
    real,    intent(out) :: a(nnmax), b(n)
    integer, intent(out) :: snr(nnmax), rnr(nn1max)
    integer :: i
    z = 0
    do i = 1, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i   ; a(z) =  3.0
    end do
    do i = 2, n
      z = z + 1 ; rnr(z) = i ; snr(z) = i-1 ; a(z) = -1.0
    end do
    do i = 1, n-1
      z = z + 1 ; rnr(z) = i ; snr(z) = i+1 ; a(z) = -1.0
    end do
    b(1:n) = 0.0
    do i = 1, z
      b(rnr(i)) = b(rnr(i)) + a(i)
    end do
  end subroutine build_tridiag_sp

end program test_y12mc_z_intent
