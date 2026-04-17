! SPDX-License-Identifier: GPL-2.0-only
!
! y12m_example_util.f90 — Utility module for the y12m example programs.
!
! Provides three sparse-matrix generators (matrd1, matre1, matrf2) and a
! timing helper (time).  These routines were previously compiled as separate
! fixed-form Fortran files; collecting them here gives all example drivers a
! single, well-defined interface.
!
! Matrix classes
! --------------
! matrd1  — class D(n,c): square, fixed NNZ = 4*n + 55.
! matre1  — class E(n,c): symmetric five-point stencil, NNZ = 5*n - 2*c - 2.
! matrf2  — class F2(m,n,c,index,alpha): general rectangular/square matrix,
!            parametrised sparsity pattern and condition number.
!
! Timing
! ------
! time(iusec)  — fills iusec with the current wall-clock reading in
!                microseconds (integer).

module y12m_example_util
   implicit none
   private

   integer, parameter :: sp = kind(1.0e0)

   public :: matrd1
   public :: matre1
   public :: matrf2
   public :: time

contains

   ! --------------------------------------------------------------------- !
   !  matrd1 — class D sparse matrix generator                             !
   ! --------------------------------------------------------------------- !
   !
   ! Generates a square n-by-n sparse matrix of class D(n,c) in COO format.
   ! The matrix depends on two parameters:
   !   n  — dimension (must be > 22 for use with the y12m benchmark driver)
   !   c  — structure parameter (positive integer < n-2)
   !
   ! The non-zero elements of the matrix (in any order) are stored in the
   ! first z positions of array a.  Their column and row numbers are stored
   ! in the corresponding positions of arrays snr and rnr respectively.
   !
   ! On exit:
   !   z       — number of non-zero elements; always z = 4*n + 55.
   !   a(1:z)  — non-zero values.
   !   snr(1:z)— column indices of the stored values.
   !   rnr(1:z)— row indices of the stored values.
   !
   subroutine matrd1(n, z, c, nn, nn1, a, snr, rnr)
      integer, intent(in)  :: n, c, nn, nn1
      integer, intent(out) :: z
      real,    intent(out) :: a(nn)
      integer, intent(out) :: snr(nn), rnr(nn1)

      integer :: i, r, s, l, k, rr1, rr2, rr3

      ! Diagonal band: a(i) = 1, column i, row i
      do i = 1, n
         a(i)   = 1.0
         snr(i) = i
         rnr(i) = i
      end do

      ! First off-diagonal band: a = i+1, shifted column c+i (mod n)
      do i = 1, n
         r    = n + i
         s    = c + i
         a(r) = real(i + 1)
         if (s <= n) then
            snr(r) = s
         else
            snr(r) = s - n
         end if
         rnr(r) = i
      end do

      ! Second off-diagonal band: a = -i, shifted column c+i+1 (mod n)
      l = 2*n
      do i = 1, n
         r    = l + i
         s    = c + i + 1
         a(r) = real(-i)
         if (s <= n) then
            snr(r) = s
         else
            snr(r) = s - n
         end if
         rnr(r) = i
      end do

      ! Third off-diagonal band: a = 16, shifted column c+i+2 (mod n)
      k = 3*n
      do i = 1, n
         r    = k + i
         s    = c + i + 2
         a(r) = 16.0
         if (s <= n) then
            snr(r) = s
         else
            snr(r) = s - n
         end if
         rnr(r) = i
      end do

      ! Dense triangular block in the upper-left corner (55 entries total)
      rr1 = 10
      rr2 = 4*n
      rr3 = 1
      do
         do i = 1, rr1
            a(rr2 + i)   = 100.0 * real(i)
            snr(rr2 + i) = n - rr1 + i
            rnr(rr2 + i) = rr3
         end do
         if (rr1 == 1) exit
         rr2 = rr2 + rr1
         rr1 = rr1 - 1
         rr3 = rr3 + 1
      end do

      z = 4*n + 55
   end subroutine matrd1


   ! --------------------------------------------------------------------- !
   !  matre1 — class E sparse matrix generator                             !
   ! --------------------------------------------------------------------- !
   !
   ! Generates a symmetric n-by-n sparse matrix of class E(n,c) in COO
   ! format.  The pattern resembles a 1-D finite-difference stencil plus
   ! two additional off-diagonal bands at distance c.
   !
   ! On exit:
   !   z = 5*n - 2*c - 2   non-zero entries.
   !
   subroutine matre1(n, z, c, nn, nn1, a, snr, rnr)
      integer, intent(in)  :: n, c, nn, nn1
      integer, intent(out) :: z
      real,    intent(out) :: a(nn)
      integer, intent(out) :: snr(nn), rnr(nn1)

      integer :: i, r, r1, r2

      ! Diagonal: a = 4
      do i = 1, n
         a(i)   = 4.0
         snr(i) = i
         rnr(i) = i
      end do

      ! Sub-diagonal: a = -1, column i+1, row i
      r = n - 1
      do i = 1, r
         r1       = n + i
         a(r1)    = -1.0
         snr(r1)  = i + 1
         rnr(r1)  = i
      end do

      ! Super-diagonal: a = -1, column i, row i+1
      r2 = 2*n - 1
      do i = 1, r
         r1       = r2 + i
         a(r1)    = -1.0
         snr(r1)  = i
         rnr(r1)  = i + 1
      end do

      ! Lower band at distance c: column i+c, row i
      r  = n - c
      r2 = 3*n - 2
      do i = 1, r
         r1       = r2 + i
         a(r1)    = -1.0
         snr(r1)  = i + c
         rnr(r1)  = i
      end do

      ! Upper band at distance c: column i, row i+c
      r2 = 4*n - 2 - c
      do i = 1, r
         r1       = r2 + i
         a(r1)    = -1.0
         snr(r1)  = i
         rnr(r1)  = i + c
      end do

      z = 5*n - 2*c - 2
   end subroutine matre1


   ! --------------------------------------------------------------------- !
   !  matrf2 — class F2 sparse matrix generator                            !
   ! --------------------------------------------------------------------- !
   !
   ! Generates sparse (rectangular or square) matrices of class F2 in COO
   ! format.  The number of rows and columns, the average number of
   ! non-zero elements per row, the sparsity pattern, and the condition
   ! number of the matrix can all be specified by the caller.
   !
   ! The non-zero elements of the matrix are accumulated (in arbitrary
   ! order) in the first nz positions of array a.  The column and row
   ! numbers of the non-zero element stored in a(i), i = 1 to nz, are
   ! found in snr(i) and rnr(i) respectively.
   !
   ! Input parameters
   ! ----------------
   !   m      — number of rows.  n <= m < HUGE(m) must hold.
   !   n      — number of columns.  21 < n < HUGE(n) must hold.
   !   c      — sparsity-pattern parameter.  10 < c < n-10 must hold.
   !   index  — average non-zeros per row.  1 < index, n-c-index >= 9
   !            must hold.
   !   alpha  — condition-number parameter.  alpha > 0.0 must hold.
   !            alpha ≈ 1.0 produces a well-conditioned matrix; large
   !            values of alpha typically produce ill-conditioned matrices.
   !            Note: setting alpha = 2**i (for any integer i whose result
   !            is representable on the machine) avoids all round-off
   !            errors during matrix generation.
   !   nn     — allocated length of a and snr.
   !            index*m + 109 < nn < HUGE(nn) must hold.
   !   nn1    — allocated length of rnr.
   !            index*m + 109 < nn1 < HUGE(nn1) must hold.
   !
   ! Output parameters
   ! -----------------
   !   nz       — number of non-zero entries stored.
   !   a(1:nz)  — non-zero values.
   !   snr(1:nz)— column numbers: the column number of a(i) is in snr(i).
   !   rnr(1:nz)— row numbers: the row number of a(i) is in rnr(i).
   !   ifejlm   — 0 on success; positive error code otherwise:
   !                1 — n out of range
   !                2 — m out of range
   !                3 — c out of range
   !                4 — index out of range (non-fatal: generation continues)
   !                5 — nn out of range
   !                6 — nn1 out of range
   !                7 — alpha out of range
   !
   subroutine matrf2(m, n, c, index, alpha, nn, nn1, nz, a, snr, rnr, ifejlm)
      integer, intent(in)  :: m, n, c, index, nn, nn1
      real,    intent(in)  :: alpha
      integer, intent(out) :: nz, ifejlm
      real,    intent(out) :: a(nn)
      integer, intent(out) :: snr(nn), rnr(nn1)

      integer, parameter :: maxint = huge(0)
      integer :: m1, nz1, k, m2, n2, index1
      integer :: i, j, j1, rr1, rr2, rr3, l
      real    :: alpha_initial, alpha_cur

      m1      = m
      ifejlm  = 0
      nz1     = index*m + 110
      k       = 1
      alpha_cur     = alpha
      alpha_initial = alpha
      index1  = index - 1

      ! ---- parameter validation ----
      if (n < 22 .or. n > maxint) then
         ifejlm = 1
      else if (m < n .or. m > maxint) then
         ifejlm = 2
      else if (c < 11 .or. n - c < 11) then
         ifejlm = 3
      else if (index < 2 .or. n - c - index < 9) then
         ! Error 4 is non-fatal: the original source did not return on this
         ! condition, so generation continues with ifejlm=4 set.
         ifejlm = 4
      else if (nn < nz1 .or. nn > maxint) then
         ifejlm = 5
      else if (nn1 < nz1 .or. nn1 > maxint) then
         ifejlm = 6
      else if (alpha <= 0.0) then
         ifejlm = 7
      end if
      if (ifejlm > 0 .and. ifejlm /= 4) return

      ! ---- generate the first n rows (the square part) ----
      ! Diagonal
      do i = 1, n
         a(i)   = 1.0
         snr(i) = i
         rnr(i) = i
      end do
      nz = n

      ! Off-diagonal bands in the square block
      j1 = 1
      do j = 1, index1
         j1 = -j1
         do i = 1, n
            a(nz + i) = real(j1 * j * i)
            if (i + c + j - 1 <= n) then
               snr(nz + i) = i + c + j - 1
            else
               snr(nz + i) = c + i + j - 1 - n
            end if
            rnr(nz + i) = i
         end do
         nz = nz + n
      end do

      ! Dense triangular corner block (55 entries)
      rr1 = 10
      rr2 = nz
      rr3 = 1
      do
         do i = 1, rr1
            a(rr2 + i)   = alpha_cur * real(i)
            snr(rr2 + i) = n - rr1 + i
            rnr(rr2 + i) = rr3
         end do
         if (rr1 == 1) exit
         rr2 = rr2 + rr1
         rr1 = rr1 - 1
         rr3 = rr3 + 1
      end do
      nz = nz + 55

      ! ---- generate the remaining m-n rows (if m > n) ----
      alpha_cur = 1.0 / alpha_cur
      m1 = m1 - n
      do while (m1 > 0)
         n2 = k * n
         if (m1 >= n) then
            m2 = n
         else
            m2 = m1
         end if
         l = 0
         if (m2 <= 10) l = 11 - m2
         do i = 1, m2
            a(nz + i)   = alpha_cur * real(k + 1)
            snr(nz + i) = i + l
            rnr(nz + i) = n2 + i
         end do
         nz = nz + m2
         j1 = 1
         do j = 1, index1
            j1 = -j1
            do i = 1, m2
               a(nz + i) = alpha_cur * real(j) * real(j1) * real((k + 1)*i + 1)
               if (i + c + j - 1 <= n) then
                  snr(nz + i) = i + c + j - 1
               else
                  snr(nz + i) = c + i + j - 1 - n
               end if
               rnr(nz + i) = n2 + i
            end do
            nz = nz + m2
         end do
         k  = k + 1
         alpha_cur = 1.0 / alpha_cur
         m1 = m1 - n
      end do

      ! Dense triangular corner block in the last rows (55 entries)
      alpha_cur = 1.0 / alpha_initial
      rr1 = 1
      rr2 = nz
      do
         do i = 1, rr1
            a(rr2 + i)   = alpha_cur * real(rr1 + 1 - i)
            snr(rr2 + i) = i
            rnr(rr2 + i) = m - 10 + rr1
         end do
         if (rr1 == 10) exit
         rr2 = rr2 + rr1
         rr1 = rr1 + 1
      end do
      nz = nz + 55
   end subroutine matrf2


   ! --------------------------------------------------------------------- !
   !  time — wall-clock timer                                               !
   ! --------------------------------------------------------------------- !
   !
   ! Returns the current wall-clock reading as an integer number of
   ! microseconds.  The result should be differenced between two calls to
   ! obtain an elapsed time; dividing the difference by 1e6 then gives
   ! seconds (as done in mainf/mainf2).
   !
   subroutine time(iusec)
      use, intrinsic :: iso_fortran_env, only: int64
      integer, intent(out) :: iusec
      integer(int64) :: count, rate
      call system_clock(count, rate)
      iusec = int(real(count, kind=sp) / real(rate, kind=sp) * 1.0e6_sp)
   end subroutine time

end module y12m_example_util
