! SPDX-License-Identifier: GPL-2.0-only
!
! test_matrix_generators.f90 — Unit tests for the three sparse-matrix
! generators in y12m_example_util.
!
! Each generator is exercised as follows:
!
!  matrd1 — class D(n,c):
!    * Verify NNZ formula: z == 4*n + 55.
!    * Verify all row/column indices lie in [1, n].
!    * Verify the matrix solves A x = b (x = all-ones) with y12mae.
!
!  matre1 — class E(n,c):
!    * Verify NNZ formula: z == 5*n - 2*c - 2.
!    * Verify all row/column indices lie in [1, n].
!    * Verify the matrix solves A x = b (x = all-ones) with y12mae.
!
!  matrf2 — class F2(m,n,c,index,alpha):
!    * Verify error codes for each out-of-range parameter.
!    * Verify all row/column indices lie in [1, m] / [1, n] for valid input.
!    * Verify the matrix solves A x = b (x = all-ones, square case) with y12mae.

program test_matrix_generators
   use y12m_example_util, only: matrd1, matre1, matrf2
   implicit none

   integer :: nfail

   nfail = 0

   call test_matrd1(nfail)
   call test_matre1(nfail)
   call test_matrf2_errors(nfail)
   call test_matrf2_valid(nfail)

   if (nfail /= 0) then
      write(*, '(i0, a)') nfail, ' test(s) FAILED'
      stop 1
   end if
   write(*, '(a)') 'All test_matrix_generators tests PASSED'

contains

   ! ------------------------------------------------------------------ !
   !  Helper: check NNZ formula, index bounds, and A x = b (x = 1)     !
   ! ------------------------------------------------------------------ !
   subroutine check_solve(label, n, z, a, snr, rnr, nn, nn1, nfail, tol)
      character(len=*), intent(in)    :: label
      integer,          intent(in)    :: n, z, nn, nn1
      real,             intent(inout) :: a(nn)
      integer,          intent(inout) :: snr(nn), rnr(nn1)
      integer,          intent(inout) :: nfail
      real,             intent(in)    :: tol

      real    :: b(n), pivot(n), aflag(8)
      integer :: ha(n, 11), iflag(10)
      integer :: i, ifail
      real    :: err

      ! Build b = A * [1, 1, ..., 1] (row sums)
      b = 0.0
      do i = 1, z
         b(rnr(i)) = b(rnr(i)) + a(i)
      end do

      call y12mae(n, z, a, snr, nn, rnr, nn1, pivot, ha, n, &
                  aflag, iflag, b, ifail)

      if (ifail /= 0) then
         write(*, '(3a, i0)') 'FAIL ', label, ' y12mae ifail=', ifail
         nfail = nfail + 1
         return
      end if

      err = maxval(abs(b(1:n) - 1.0))
      if (err > tol) then
         write(*, '(3a, es10.3, a, es10.3)') 'FAIL ', label, &
            ' max_err=', err, ' > tol=', tol
         nfail = nfail + 1
      else
         write(*, '(3a, es10.3)') 'PASS ', label, ' max_err=', err
      end if
   end subroutine check_solve


   ! ------------------------------------------------------------------ !
   !  matrd1 tests                                                        !
   ! ------------------------------------------------------------------ !
   subroutine test_matrd1(nfail)
      integer, intent(inout) :: nfail

      integer, parameter :: N = 50, C = 5
      integer, parameter :: NN = 25*N, NN1 = 20*N

      real    :: a(NN)
      integer :: snr(NN), rnr(NN1)
      integer :: z, i

      call matrd1(N, z, C, NN, NN1, a, snr, rnr)

      ! NNZ formula
      if (z /= 4*N + 55) then
         write(*, '(a, i0, a, i0)') &
            'FAIL matrd1 NNZ: expected ', 4*N + 55, ', got ', z
         nfail = nfail + 1
      else
         write(*, '(a, i0)') 'PASS matrd1 NNZ z=', z
      end if

      ! Index bounds
      do i = 1, z
         if (snr(i) < 1 .or. snr(i) > N .or. rnr(i) < 1 .or. rnr(i) > N) then
            write(*, '(a, i0, a, i0, a, i0)') &
               'FAIL matrd1 index bounds at i=', i, &
               ' snr=', snr(i), ' rnr=', rnr(i)
            nfail = nfail + 1
            return
         end if
      end do
      write(*, '(a)') 'PASS matrd1 index bounds'

      ! Solve A x = b, x should be all-ones
      call matrd1(N, z, C, NN, NN1, a, snr, rnr)
      call check_solve('matrd1 solve', N, z, a, snr, rnr, NN, NN1, nfail, 1.0e-2)
   end subroutine test_matrd1


   ! ------------------------------------------------------------------ !
   !  matre1 tests                                                        !
   ! ------------------------------------------------------------------ !
   subroutine test_matre1(nfail)
      integer, intent(inout) :: nfail

      integer, parameter :: N = 50, C = 5
      integer, parameter :: NN = 10*N, NN1 = 10*N

      real    :: a(NN)
      integer :: snr(NN), rnr(NN1)
      integer :: z, i

      call matre1(N, z, C, NN, NN1, a, snr, rnr)

      ! NNZ formula
      if (z /= 5*N - 2*C - 2) then
         write(*, '(a, i0, a, i0)') &
            'FAIL matre1 NNZ: expected ', 5*N - 2*C - 2, ', got ', z
         nfail = nfail + 1
      else
         write(*, '(a, i0)') 'PASS matre1 NNZ z=', z
      end if

      ! Index bounds
      do i = 1, z
         if (snr(i) < 1 .or. snr(i) > N .or. rnr(i) < 1 .or. rnr(i) > N) then
            write(*, '(a, i0, a, i0, a, i0)') &
               'FAIL matre1 index bounds at i=', i, &
               ' snr=', snr(i), ' rnr=', rnr(i)
            nfail = nfail + 1
            return
         end if
      end do
      write(*, '(a)') 'PASS matre1 index bounds'

      ! Solve A x = b, x should be all-ones
      call matre1(N, z, C, NN, NN1, a, snr, rnr)
      call check_solve('matre1 solve', N, z, a, snr, rnr, NN, NN1, nfail, 1.0e-4)
   end subroutine test_matre1


   ! ------------------------------------------------------------------ !
   !  matrf2 error-code tests                                             !
   ! ------------------------------------------------------------------ !
   subroutine test_matrf2_errors(nfail)
      integer, intent(inout) :: nfail

      ! Valid base parameters (square, n=50, c=15, index=3, alpha=1)
      integer, parameter :: M = 50, N = 50, C = 15, IDX = 3
      real,    parameter :: ALPHA = 1.0
      integer, parameter :: NN = 100*M, NN1 = 100*M

      real    :: a(NN)
      integer :: snr(NN), rnr(NN1)
      integer :: nz, ifejlm

      ! n < 22 → ifejlm=1
      call matrf2(22, 21, C, IDX, ALPHA, NN, NN1, nz, a, snr, rnr, ifejlm)
      call check_err('matrf2 n<22', ifejlm, 1, nfail)

      ! m < n → ifejlm=2
      call matrf2(N - 1, N, C, IDX, ALPHA, NN, NN1, nz, a, snr, rnr, ifejlm)
      call check_err('matrf2 m<n', ifejlm, 2, nfail)

      ! c < 11 → ifejlm=3
      call matrf2(M, N, 5, IDX, ALPHA, NN, NN1, nz, a, snr, rnr, ifejlm)
      call check_err('matrf2 c<11', ifejlm, 3, nfail)

      ! index < 2 → ifejlm=4 (note: generator still runs, error is non-fatal)
      call matrf2(M, N, C, 1, ALPHA, NN, NN1, nz, a, snr, rnr, ifejlm)
      call check_err('matrf2 index<2', ifejlm, 4, nfail)

      ! alpha <= 0 → ifejlm=7
      call matrf2(M, N, C, IDX, -1.0, NN, NN1, nz, a, snr, rnr, ifejlm)
      call check_err('matrf2 alpha<=0', ifejlm, 7, nfail)

   end subroutine test_matrf2_errors


   subroutine check_err(label, got, expected, nfail)
      character(len=*), intent(in)    :: label
      integer,          intent(in)    :: got, expected
      integer,          intent(inout) :: nfail
      if (got == expected) then
         write(*, '(3a, i0)') 'PASS ', label, ' ifejlm=', got
      else
         write(*, '(3a, i0, a, i0)') 'FAIL ', label, &
            ' expected ifejlm=', expected, ', got ', got
         nfail = nfail + 1
      end if
   end subroutine check_err


   ! ------------------------------------------------------------------ !
   !  matrf2 valid-input test (square case, x = all-ones)                !
   ! ------------------------------------------------------------------ !
   subroutine test_matrf2_valid(nfail)
      integer, intent(inout) :: nfail

      integer, parameter :: N = 50, C = 15, IDX = 3
      real,    parameter :: ALPHA = 1.0
      integer, parameter :: NN = 100*N, NN1 = 100*N

      real    :: a(NN)
      integer :: snr(NN), rnr(NN1)
      integer :: nz, ifejlm, i

      call matrf2(N, N, C, IDX, ALPHA, NN, NN1, nz, a, snr, rnr, ifejlm)

      if (ifejlm /= 0) then
         write(*, '(a, i0)') 'FAIL matrf2 valid ifejlm=', ifejlm
         nfail = nfail + 1
         return
      end if
      write(*, '(a, i0)') 'PASS matrf2 valid NNZ nz=', nz

      ! Index bounds: rows in [1,N], columns in [1,N]
      do i = 1, nz
         if (snr(i) < 1 .or. snr(i) > N .or. rnr(i) < 1 .or. rnr(i) > N) then
            write(*, '(a, i0, a, i0, a, i0)') &
               'FAIL matrf2 valid index bounds at i=', i, &
               ' snr=', snr(i), ' rnr=', rnr(i)
            nfail = nfail + 1
            return
         end if
      end do
      write(*, '(a)') 'PASS matrf2 valid index bounds'

      ! Solve A x = b
      call matrf2(N, N, C, IDX, ALPHA, NN, NN1, nz, a, snr, rnr, ifejlm)
      call check_solve('matrf2 valid solve', N, nz, a, snr, rnr, NN, NN1, nfail, 1.0e-4)
   end subroutine test_matrf2_valid

end program test_matrix_generators
