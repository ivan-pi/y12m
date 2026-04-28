! SPDX-License-Identifier: GPL-2.0-only
!
! matf_bench.f90 — Timing driver for class F2 sparse matrices.
!
! The F2 matrices are square (n-by-n) and generalise class D by adding a
! lower-left 10-by-10 triangular block of elements.  (r-1) is the width of
! a cyclic off-diagonal band located at a distance c from the main diagonal.
!
! Usage
! -----
!   matf_bench [--help] [-n RANGE] [-c RANGE] [-r RANGE] [--alpha VALUE]
!
!   RANGE  a single integer (e.g. 100) or begin:end:step (e.g. 100:400:100)
!
! Defaults:
!   -n 100:400:100  -c 25:75:25  -r 3  --alpha 1.0
!
! Class F2(n,c,r,alpha): square, NNZ = r*n + 110
!   n >= 22,  11 <= c <= n-11,  2 <= r <= n-c-9,  alpha >= 1
!
! Output
! ------
!   For every (n, c, r): elapsed time, max error, and the same aflag /
!   iflag diagnostics as the original example drivers (mainf.f, mainf2.f).
!   At the end: total elapsed wall-clock time.
!
! ==========================================================================
! Module f_matrix_mod
! ==========================================================================

module f_matrix_mod
   use, intrinsic :: iso_fortran_env, only: int64
   implicit none
   private

   integer, parameter, public :: dp = kind(1.0d0)

   !> F2(n,c,r,alpha) sparse matrix in COO format (square case, m = n).
   !>
   !> Generalises class D by adding a lower-left 10-by-10 triangular block.
   !> (r-1) is the width of the cyclic off-diagonal band at distance c from
   !> the main diagonal.
   type, public :: f2_matrix
      integer  :: n     = 0         !< Matrix dimension (n-by-n)
      integer  :: r     = 0         !< Band-width parameter (index in matrf2)
      real(dp) :: alpha = 1.0_dp    !< Condition-number parameter
      integer  :: nnz   = 0         !< Non-zeros in the original matrix
      integer  :: nn    = 0         !< Allocated length of a / snr
      integer  :: nn1   = 0         !< Allocated length of rnr
      real(dp), allocatable :: a(:)
      integer,  allocatable :: snr(:), rnr(:)
   contains
      procedure :: fill  => fill_f2
      procedure :: solve => f2_solve
   end type f2_matrix

contains

   ! ------------------------------------------------------------------
   !> Populate COO arrays for the square F2(n,c,r,alpha) matrix.
   !>
   !> Pattern (square, m = n):
   !>   a(i,i)                               = 1      (diagonal)
   !>   a(i, mod(c+i+j-2, n)+1) = (-1)^j*j*i          (j = 1..r-1)
   !>   a(1..10, n-9..n)        = alpha * col          (upper-right triangle)
   !>   a(n-9..n, 1..10)        = (1/alpha) * values   (lower-left triangle)
   !>
   !> NNZ = r*n + 110.
   subroutine fill_f2(self, n, c, r, alpha)
      class(f2_matrix), intent(inout) :: self
      integer,  intent(in) :: n, c, r
      real(dp), intent(in) :: alpha

      integer  :: i, j, j1, rr1, rr2, rr3, nz, index1
      real(dp) :: alpha_inv

      self%n     = n
      self%r     = r
      self%alpha = alpha
      self%nnz   = r * n + 110
      self%nn    = max(25 * n, 5 * self%nnz)
      self%nn1   = max(20 * n, 4 * self%nnz)

      if (allocated(self%a))   deallocate(self%a)
      if (allocated(self%snr)) deallocate(self%snr)
      if (allocated(self%rnr)) deallocate(self%rnr)

      allocate(self%a(self%nn))
      allocate(self%snr(self%nn))
      allocate(self%rnr(self%nn1))

      index1 = r - 1

      ! Diagonal: a(i,i) = 1
      do i = 1, n
         self%a(i)   = 1.0_dp
         self%snr(i) = i
         self%rnr(i) = i
      end do
      nz = n

      ! Off-diagonal cyclic bands (r-1 bands of n entries each)
      j1 = 1
      do j = 1, index1
         j1 = -j1
         do i = 1, n
            self%a(nz + i) = real(j1 * j * i, dp)
            if (i + c + j - 1 <= n) then
               self%snr(nz + i) = i + c + j - 1
            else
               self%snr(nz + i) = c + i + j - 1 - n
            end if
            self%rnr(nz + i) = i
         end do
         nz = nz + n
      end do

      ! Dense triangular block in the upper-right corner (rows 1..10, 55 entries)
      rr1 = 10
      rr2 = nz
      rr3 = 1
      do
         do i = 1, rr1
            self%a(rr2 + i)   = alpha * real(i, dp)
            self%snr(rr2 + i) = n - rr1 + i
            self%rnr(rr2 + i) = rr3
         end do
         if (rr1 == 1) exit
         rr2 = rr2 + rr1
         rr1 = rr1 - 1
         rr3 = rr3 + 1
      end do
      nz = nz + 55

      ! Dense triangular block in the lower-left corner (rows n-9..n, 55 entries)
      rr1 = 1
      rr2 = nz
      alpha_inv = 1.0_dp / alpha
      do
         do i = 1, rr1
            self%a(rr2 + i)   = alpha_inv * real(rr1 + 1 - i, dp)
            self%snr(rr2 + i) = i
            self%rnr(rr2 + i) = n - 10 + rr1
         end do
         if (rr1 == 10) exit
         rr2 = rr2 + rr1
         rr1 = rr1 + 1
      end do
      nz = nz + 55

      self%nnz = nz
   end subroutine fill_f2

   ! ------------------------------------------------------------------
   !> Solve A x = b where b = A * 1 (row sums; exact solution is all-ones).
   !>
   !> Uses the y12ma generic interface (double-precision path).
   !>
   !> Returns:
   !>   aflag(8)  — diagnostic reals as set by y12ma
   !>   iflag(10) — diagnostic integers as set by y12ma
   !>   ifail     — error diagnostic parameter (0 = success)
   !>   elapsed   — wall-clock time in seconds for the y12ma call
   !>   max_err   — max(|x - 1|); or -1 if ifail /= 0
   subroutine f2_solve(self, aflag, iflag, ifail, elapsed, max_err)
      use y12m, only: y12ma
      class(f2_matrix), intent(inout) :: self
      real(dp), intent(out) :: aflag(8)
      integer,  intent(out) :: iflag(10)
      integer,  intent(out) :: ifail
      real(dp), intent(out) :: elapsed
      real(dp), intent(out) :: max_err

      real(dp), allocatable :: b(:), pivot(:)
      integer,  allocatable :: ha(:,:)
      integer(int64) :: t1, t2, rate
      integer :: i, z_local

      allocate(b(self%n), pivot(self%n), ha(self%n, 11))

      ! Build b = A * [1, 1, ..., 1] (row sums; exact solution is all-ones)
      b = 0.0_dp
      do i = 1, self%nnz
         b(self%rnr(i)) = b(self%rnr(i)) + self%a(i)
      end do

      z_local = self%nnz
      ifail   = 0

      call system_clock(t1, rate)
      call y12ma(self%n, z_local, self%a, self%snr, self%nn, &
                 self%rnr, self%nn1, pivot, ha, self%n, &
                 aflag, iflag, b, ifail)
      call system_clock(t2)

      elapsed = real(t2 - t1, dp) / real(rate, dp)

      if (ifail == 0) then
         max_err = maxval(abs(b(1:self%n) - 1.0_dp))
      else
         max_err = -1.0_dp
      end if
   end subroutine f2_solve

end module f_matrix_mod


! ==========================================================================
! Program matf_bench
! ==========================================================================

program matf_bench
   use f_matrix_mod
   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit, int64
   implicit none

   ! ---- Local type for command-line arguments (includes defaults) --------
   type :: cmdargs
      integer  :: nstart =  100
      integer  :: nstep  =  100
      integer  :: nend   =  400
      integer  :: cstart =   25
      integer  :: cstep  =   25
      integer  :: cend   =   75
      integer  :: rstart =    3
      integer  :: rstep  =    1
      integer  :: rend   =    3
      real(dp) :: alpha  = 1.0_dp
   end type cmdargs

   ! ---- Working variables -----------------------------------------------
   type(cmdargs) :: args
   type(f2_matrix) :: mat
   real(dp) :: aflag(8)
   integer  :: iflag(10), ifail
   real(dp) :: elapsed, max_err, total_elapsed
   integer(int64) :: t_prog_start, t_prog_end, wall_rate
   integer :: n, c, r

   args = parse_args()

   call system_clock(t_prog_start, wall_rate)
   total_elapsed = 0.0_dp

   ! ---- Main loop ---------------------------------------------------------
   ! F2(n,c,r,alpha): n >= 22, 11 <= c <= n-11, 2 <= r <= n-c-9, alpha >= 1
   do n = args%nstart, args%nend, args%nstep
      do c = args%cstart, args%cend, args%cstep
         do r = args%rstart, args%rend, args%rstep

            if (n < 22) then
               write(output_unit, '(a,i0,a)') &
                  '  Skipping F2(n=', n, '): n must be >= 22'
               cycle
            end if
            if (c < 11 .or. c > n - 11) then
               write(output_unit, '(a,i0,a,i0,a,i0,a,i0,a)') &
                  '  Skipping F2(n=', n, ', c=', c, &
                  '): c must be in [11, n-11=', n - 11, ']'
               cycle
            end if
            if (r < 2 .or. r > n - c - 9) then
               write(output_unit, '(a,i0,a,i0,a,i0,a,i0,a)') &
                  '  Skipping F2(n=', n, ', c=', c, ', r=', r, &
                  '): r must be in [2, n-c-9=', n - c - 9, ']'
               cycle
            end if

            write(output_unit, '(/,a,i0,3(a,i0),a,1pe10.2)') &
               'Matrix of class F2:  n=', n, '  c=', c, '  r=', r, &
               '  NNZ=', r*n + 110, '  alpha=', args%alpha
            call mat%fill(n, c, r, args%alpha)
            write(output_unit, '(a,i0)') &
               '  Non-zeros in original matrix: ', mat%nnz
            call mat%solve(aflag, iflag, ifail, elapsed, max_err)
            total_elapsed = total_elapsed + elapsed
            call print_solve_stats(aflag, iflag, ifail, elapsed, max_err)

         end do
      end do
   end do

   ! ---- Total elapsed wall-clock time ------------------------------------
   call system_clock(t_prog_end)
   write(output_unit, '(/,a,G12.3,a)') &
      'Total elapsed wall-clock time: ', &
      real(t_prog_end - t_prog_start, dp) / real(wall_rate, dp), ' s'

contains

   ! -----------------------------------------------------------------------
   !> Print the per-solve diagnostics in the same format as the original
   !> example drivers (mainf.f, mainf2.f).
   subroutine print_solve_stats(aflag, iflag, ifail, elapsed, max_err)
      real(dp), intent(in) :: aflag(8)
      integer,  intent(in) :: iflag(10)
      integer,  intent(in) :: ifail
      real(dp), intent(in) :: elapsed
      real(dp), intent(in) :: max_err

      write(output_unit, '(a,f12.3,a)') &
         '  elapsed time:                    ', elapsed, ' s'
      if (max_err >= 0.0_dp) &
         write(output_unit, '(a,1pe12.3)') &
            '  largest error:                   ', max_err
      write(output_unit, '(a,i8)') &
         '  largest no. of elements in A:   ', iflag(8)
      write(output_unit, '(a,i4)') &
         '  collections in row list:         ', iflag(6)
      write(output_unit, '(a,i4)') &
         '  collections in column list:      ', iflag(7)
      write(output_unit, '(a,1pe12.3)') &
         '  largest element in orig. matrix: ', aflag(6)
      write(output_unit, '(a,1pe12.3)') &
         '  largest element in LU matrix:    ', aflag(7)
      write(output_unit, '(a,1pe12.3)') &
         '  growth factor:                   ', aflag(5)
      write(output_unit, '(a,1pe12.3)') &
         '  minimal pivotal element:         ', aflag(8)
      write(output_unit, '(a,1pe12.3,a,1pe12.3)') &
         '  drop tolerance: ', aflag(2), &
         '   stability factor: ', aflag(1)
      write(output_unit, '(a,i0)') &
         '  error diagnostic parameter:      ', ifail
   end subroutine print_solve_stats

   ! -----------------------------------------------------------------------
   !> Parse command-line arguments and return a cmdargs instance.
   !> Unknown flags cause a help printout followed by a non-zero exit.
   function parse_args() result(args)
      type(cmdargs) :: args

      integer :: i, nargs, ios
      character(len=64) :: arg, val

      nargs = command_argument_count()
      i = 1
      do while (i <= nargs)
         call get_command_argument(i, arg)
         select case (trim(arg))
         case ('--help', '-h')
            call print_help()
            stop 0

         case ('-n')
            i = i + 1
            if (i > nargs) then
               write(error_unit, '(a)') 'Error: -n requires an argument'
               stop 1
            end if
            call get_command_argument(i, val)
            call parse_range(val, args%nstart, args%nend, args%nstep)

         case ('-c')
            i = i + 1
            if (i > nargs) then
               write(error_unit, '(a)') 'Error: -c requires an argument'
               stop 1
            end if
            call get_command_argument(i, val)
            call parse_range(val, args%cstart, args%cend, args%cstep)

         case ('-r')
            i = i + 1
            if (i > nargs) then
               write(error_unit, '(a)') 'Error: -r requires an argument'
               stop 1
            end if
            call get_command_argument(i, val)
            call parse_range(val, args%rstart, args%rend, args%rstep)

         case ('--alpha')
            i = i + 1
            if (i > nargs) then
               write(error_unit, '(a)') 'Error: --alpha requires an argument'
               stop 1
            end if
            call get_command_argument(i, val)
            read(val, *, iostat=ios) args%alpha
            if (ios /= 0) then
               write(error_unit, '(2a)') &
                  'Error: cannot parse alpha value: ', trim(val)
               stop 1
            end if

         case default
            write(error_unit, '(2a)') 'Error: unknown argument: ', trim(arg)
            call print_help()
            stop 1
         end select
         i = i + 1
      end do
   end function parse_args

   ! -----------------------------------------------------------------------
   !> Parse a range string of the form:
   !>   VALUE            -> vstart=VALUE, vend=VALUE, vstep=1
   !>   BEGIN:END        -> vstart=BEGIN, vend=END,   vstep=1
   !>   BEGIN:END:STEP   -> vstart=BEGIN, vend=END,   vstep=STEP
   subroutine parse_range(str, vstart, vend, vstep)
      character(len=*), intent(in)  :: str
      integer, intent(out) :: vstart, vend, vstep

      integer :: pos1, pos2, ios
      character(len=64) :: part1, part2, part3

      pos1 = index(str, ':')
      if (pos1 == 0) then
         ! Single value
         read(str, *, iostat=ios) vstart
         if (ios /= 0) then
            write(error_unit, '(2a)') 'Error: cannot parse value: ', trim(str)
            stop 1
         end if
         vend  = vstart
         vstep = 1
         return
      end if

      ! Look for the second colon in the substring after pos1
      pos2 = index(str(pos1 + 1:), ':')

      if (pos2 == 0) then
         ! Two parts: BEGIN:END
         part1 = str(1:pos1 - 1)
         part2 = str(pos1 + 1:)
         read(part1, *, iostat=ios) vstart
         if (ios /= 0) then
            write(error_unit, '(2a)') &
               'Error: cannot parse range start: ', trim(part1)
            stop 1
         end if
         read(part2, *, iostat=ios) vend
         if (ios /= 0) then
            write(error_unit, '(2a)') &
               'Error: cannot parse range end: ', trim(part2)
            stop 1
         end if
         vstep = 1
      else
         ! Three parts: BEGIN:END:STEP (pos2 is relative, make it absolute)
         pos2  = pos2 + pos1
         part1 = str(1:pos1 - 1)
         part2 = str(pos1 + 1:pos2 - 1)
         part3 = str(pos2 + 1:)
         read(part1, *, iostat=ios) vstart
         if (ios /= 0) then
            write(error_unit, '(2a)') &
               'Error: cannot parse range start: ', trim(part1)
            stop 1
         end if
         read(part2, *, iostat=ios) vend
         if (ios /= 0) then
            write(error_unit, '(2a)') &
               'Error: cannot parse range end: ', trim(part2)
            stop 1
         end if
         read(part3, *, iostat=ios) vstep
         if (ios /= 0) then
            write(error_unit, '(2a)') &
               'Error: cannot parse range step: ', trim(part3)
            stop 1
         end if
      end if
   end subroutine parse_range

   ! -----------------------------------------------------------------------
   !> Print the help message and return.
   subroutine print_help()
      write(output_unit, '(a)') &
         'Usage: matf_bench [--help] [-n RANGE] [-c RANGE]' // &
         ' [-r RANGE] [--alpha VALUE]'
      write(output_unit, '(a)') ''
      write(output_unit, '(a)') &
         'Timing driver for class F2 sparse matrices (double precision).'
      write(output_unit, '(a)') &
         'Reproduces the F-matrix benchmark from Zlatev et al.'
      write(output_unit, '(a)') ''
      write(output_unit, '(a)') &
         'Class F2(n,c,r,alpha): square, NNZ = r*n + 110'
      write(output_unit, '(a)') &
         '  n >= 22,  11 <= c <= n-11,  2 <= r <= n-c-9,  alpha >= 1'
      write(output_unit, '(a)') ''
      write(output_unit, '(a)') 'Options:'
      write(output_unit, '(a)') &
         '  --help, -h     Show this help message and exit.'
      write(output_unit, '(a)') &
         '  -n RANGE       Matrix dimension (default: 100:400:100).'
      write(output_unit, '(a)') &
         '  -c RANGE       Band-offset parameter (default: 25:75:25).'
      write(output_unit, '(a)') &
         '  -r RANGE       Band-width parameter (default: 3).'
      write(output_unit, '(a)') &
         '  --alpha VALUE  Condition-number parameter (default: 1.0).'
      write(output_unit, '(a)') ''
      write(output_unit, '(a)') 'RANGE syntax:'
      write(output_unit, '(a)') &
         '  VALUE            single integer, e.g. 100'
      write(output_unit, '(a)') &
         '  BEGIN:END        step defaults to 1'
      write(output_unit, '(a)') &
         '  BEGIN:END:STEP   Fortran-slice style, e.g. 100:400:100'
      write(output_unit, '(a)') ''
      write(output_unit, '(a)') 'Defaults:'
      write(output_unit, '(a)') &
         '  -n 100:400:100  -c 25:75:25  -r 3  --alpha 1.0'
   end subroutine print_help

end program matf_bench
