! SPDX-License-Identifier: GPL-2.0-only
!
! timings_de.f90 — Timing driver for class D and E sparse matrices.
!
! Reproduces the timing benchmark from Zlatev et al. using double precision
! arithmetic and dynamic memory allocation.  The design is object-oriented:
! an abstract COO matrix parent type provides storage and the solver call;
! the two concrete children matrix_d and matrix_e generate their respective
! patterns.
!
! Usage
! -----
!   timings_de [--help] [-class {d|e|de}] [-n RANGE] [-c RANGE]
!
!   RANGE  a single integer (e.g. 650) or begin:end:step (e.g. 650:1000:50)
!
! Defaults (Zlatev benchmark table):
!   -class de  -n 650:1000:50  -c 4:204:40
!
! Class D(n,c): square,          NNZ = 4*n + 55
! Class E(n,c): symmetric SPD,   NNZ = 5*n - 2*c - 2
!
! Output
! ------
!   For every (class, n, c): elapsed time, max error, and the same aflag /
!   iflag diagnostics as the original example drivers (maind.f, maine.f).
!   At the end: total elapsed wall-clock time.
!
! ==========================================================================
! Module coo_matrix_mod
! ==========================================================================

module coo_matrix_mod
   use, intrinsic :: iso_fortran_env, only: real64, int64
   implicit none
   private

   integer, parameter, public :: dp = real64

   ! ------------------------------------------------------------------
   ! Abstract base type
   ! ------------------------------------------------------------------

   !> Abstract COO sparse matrix.
   !>
   !> Concrete subclasses override `fill` to populate the COO arrays.
   !> The parent provides `prepare` (allocates storage, then calls fill)
   !> and `solve` (invokes y12maf and reports diagnostics).
   type, public, abstract :: coo_matrix
      integer :: n   = 0  !< Matrix dimension
      integer :: c   = 0  !< Structure parameter
      integer :: nnz = 0  !< Non-zeros in the original matrix
      integer :: nn  = 0  !< Allocated length of a / snr
      integer :: nn1 = 0  !< Allocated length of rnr
      real(dp), allocatable :: a(:)
      integer,  allocatable :: snr(:), rnr(:)
   contains
      procedure(fill_iface), deferred :: fill
      procedure :: prepare => coo_prepare
      procedure :: solve   => coo_solve
   end type coo_matrix

   abstract interface
      !> Populate the COO arrays.  Called by coo_prepare after allocation.
      !> On entry self%n, self%c, self%nn, self%nn1 are set and
      !> self%a(1:nn), self%snr(1:nn), self%rnr(1:nn1) are allocated.
      !> On exit self%nnz holds the non-zero count.
      subroutine fill_iface(self)
         import :: coo_matrix
         class(coo_matrix), intent(inout) :: self
      end subroutine fill_iface
   end interface

   ! ------------------------------------------------------------------
   ! Concrete subclasses
   ! ------------------------------------------------------------------

   !> Class D(n,c) — square, NNZ = 4*n + 55.
   type, public, extends(coo_matrix) :: matrix_d
   contains
      procedure :: fill => fill_d
   end type matrix_d

   !> Class E(n,c) — symmetric SPD five-point stencil, NNZ = 5*n - 2*c - 2.
   type, public, extends(coo_matrix) :: matrix_e
   contains
      procedure :: fill => fill_e
   end type matrix_e

contains

   ! ------------------------------------------------------------------
   !> Allocate storage sized for LU fill-in and populate the matrix.
   subroutine coo_prepare(self, n, c)
      class(coo_matrix), intent(inout) :: self
      integer, intent(in) :: n, c

      self%n   = n
      self%c   = c
      self%nn  = 25 * n
      self%nn1 = 20 * n

      if (allocated(self%a))   deallocate(self%a)
      if (allocated(self%snr)) deallocate(self%snr)
      if (allocated(self%rnr)) deallocate(self%rnr)

      allocate(self%a(self%nn))
      allocate(self%snr(self%nn))
      allocate(self%rnr(self%nn1))

      call self%fill()
   end subroutine coo_prepare

   ! ------------------------------------------------------------------
   !> Solve A x = b where b = A * 1 (row sums; exact solution is all-ones).
   !>
   !> Uses y12maf (double precision) directly to avoid the incorrect
   !> intent declarations in the y12m generic wrapper module.
   !>
   !> Returns:
   !>   aflag(8)  — diagnostic reals as set by y12maf
   !>   iflag(10) — diagnostic integers as set by y12maf
   !>   ifail     — error diagnostic parameter (0 = success)
   !>   elapsed   — wall-clock time in seconds for the y12maf call
   !>   max_err   — max(|x - 1|); or -1 if ifail /= 0
   subroutine coo_solve(self, aflag, iflag, ifail, elapsed, max_err)
      class(coo_matrix), intent(inout) :: self
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

      ! Build b = A * 1 (row sums; exact solution is the all-ones vector)
      b = 0.0_dp
      do i = 1, self%nnz
         b(self%rnr(i)) = b(self%rnr(i)) + self%a(i)
      end do

      z_local = self%nnz
      ifail   = 0

      call system_clock(t1, rate)
      call y12maf(self%n, z_local, self%a, self%snr, self%nn, &
                  self%rnr, self%nn1, pivot, ha, self%n, &
                  aflag, iflag, b, ifail)
      call system_clock(t2)

      elapsed = real(t2 - t1, dp) / real(rate, dp)

      if (ifail == 0) then
         max_err = maxval(abs(b(1:self%n) - 1.0_dp))
      else
         max_err = -1.0_dp
      end if

      deallocate(b, pivot, ha)
   end subroutine coo_solve

   ! ------------------------------------------------------------------
   ! fill_d — populate COO arrays for class D(n,c)
   !
   ! Pattern:
   !   a(i,i)                       = 1    (diagonal)
   !   a(i, c+i wrapping around n)  = i+1  (band 1, cyclic column c+i)
   !   a(i, c+i+1 wrapping)         = -i   (band 2, cyclic column c+i+1)
   !   a(i, c+i+2 wrapping)         = 16   (band 3, cyclic column c+i+2)
   !   a(i, n-10+j) = 100*j, i=1..11-j, j=1..10  (dense 10×10 corner)
   ! NNZ = 4*n + 55
   ! ------------------------------------------------------------------
   subroutine fill_d(self)
      class(matrix_d), intent(inout) :: self

      integer :: n, c, i, r, s, l, k, rr1, rr2, rr3

      n = self%n
      c = self%c

      ! Diagonal: a(i,i) = 1
      do i = 1, n
         self%a(i)   = 1.0_dp
         self%snr(i) = i
         self%rnr(i) = i
      end do

      ! Band 1: a(i, c+i) = i+1, column wraps cyclically around n
      do i = 1, n
         r = n + i
         s = c + i
         self%a(r) = real(i + 1, dp)
         if (s <= n) then
            self%snr(r) = s
         else
            self%snr(r) = s - n
         end if
         self%rnr(r) = i
      end do

      ! Band 2: a(i, c+i+1) = -i, column wraps cyclically around n
      l = 2 * n
      do i = 1, n
         r = l + i
         s = c + i + 1
         self%a(r) = real(-i, dp)
         if (s <= n) then
            self%snr(r) = s
         else
            self%snr(r) = s - n
         end if
         self%rnr(r) = i
      end do

      ! Band 3: a(i, c+i+2) = 16, column wraps cyclically around n
      k = 3 * n
      do i = 1, n
         r = k + i
         s = c + i + 2
         self%a(r) = 16.0_dp
         if (s <= n) then
            self%snr(r) = s
         else
            self%snr(r) = s - n
         end if
         self%rnr(r) = i
      end do

      ! Dense triangular block in the upper-right corner (55 entries)
      rr1 = 10
      rr2 = 4 * n
      rr3 = 1
      do
         do i = 1, rr1
            self%a(rr2 + i)   = 100.0_dp * real(i, dp)
            self%snr(rr2 + i) = n - rr1 + i
            self%rnr(rr2 + i) = rr3
         end do
         if (rr1 == 1) exit
         rr2 = rr2 + rr1
         rr1 = rr1 - 1
         rr3 = rr3 + 1
      end do

      self%nnz = 4 * n + 55
   end subroutine fill_d

   ! ------------------------------------------------------------------
   ! fill_e — populate COO arrays for class E(n,c)
   !
   ! Pattern (symmetric matrix):
   !   a(i,i)     = 4       (diagonal)
   !   a(i,i+1)   = -1      (upper side-diagonal, i=1..n-1)
   !   a(i+1,i)   = -1      (lower side-diagonal, i=1..n-1)
   !   a(i,i+c)   = -1      (upper band at distance c, i=1..n-c)
   !   a(i+c,i)   = -1      (lower band at distance c, i=1..n-c)
   ! NNZ = 5*n - 2*c - 2
   ! ------------------------------------------------------------------
   subroutine fill_e(self)
      class(matrix_e), intent(inout) :: self

      integer :: n, c, i, r, r1, r2

      n = self%n
      c = self%c

      ! Diagonal: a(i,i) = 4
      do i = 1, n
         self%a(i)   = 4.0_dp
         self%snr(i) = i
         self%rnr(i) = i
      end do

      ! Upper side-diagonal: a(i, i+1) = -1 for i = 1..n-1
      r = n - 1
      do i = 1, r
         r1           = n + i
         self%a(r1)   = -1.0_dp
         self%snr(r1) = i + 1
         self%rnr(r1) = i
      end do

      ! Lower side-diagonal: a(i+1, i) = -1 for i = 1..n-1
      r2 = 2 * n - 1
      do i = 1, r
         r1           = r2 + i
         self%a(r1)   = -1.0_dp
         self%snr(r1) = i
         self%rnr(r1) = i + 1
      end do

      ! Upper band at distance c: a(i, i+c) = -1 for i = 1..n-c
      r  = n - c
      r2 = 3 * n - 2
      do i = 1, r
         r1           = r2 + i
         self%a(r1)   = -1.0_dp
         self%snr(r1) = i + c
         self%rnr(r1) = i
      end do

      ! Lower band at distance c: a(i+c, i) = -1 for i = 1..n-c
      r2 = 4 * n - 2 - c
      do i = 1, r
         r1           = r2 + i
         self%a(r1)   = -1.0_dp
         self%snr(r1) = i
         self%rnr(r1) = i + c
      end do

      self%nnz = 5 * n - 2 * c - 2
   end subroutine fill_e

end module coo_matrix_mod


! ==========================================================================
! Program timings_de
! ==========================================================================

program timings_de
   use coo_matrix_mod
   use, intrinsic :: iso_fortran_env, only: output_unit, error_unit, int64
   implicit none

   ! ---- Parameter ranges (default = Zlatev benchmark) ------------------
   integer :: nstart, nstep, nend
   integer :: cstart, cstep, cend
   logical :: run_d, run_e

   ! ---- Working variables -----------------------------------------------
   type(matrix_d) :: mat_d
   type(matrix_e) :: mat_e
   real(dp) :: aflag(8)
   integer  :: iflag(10), ifail
   real(dp) :: elapsed, max_err, total_elapsed
   integer(int64) :: t_prog_start, t_prog_end, wall_rate
   integer :: n, c

   ! ---- Defaults (Zlatev benchmark: n=650(50)1000, c=4(40)204) ---------
   nstart = 650
   nstep  = 50
   nend   = 1000
   cstart = 4
   cstep  = 40
   cend   = 204
   run_d  = .true.
   run_e  = .true.

   call parse_args(nstart, nstep, nend, cstart, cstep, cend, run_d, run_e)

   call system_clock(t_prog_start, wall_rate)
   total_elapsed = 0.0_dp

   ! ---- Class D -----------------------------------------------------------
   ! Valid parameter range: n > 22, 1 <= c <= n-13
   if (run_d) then
      do n = nstart, nend, nstep
         do c = cstart, cend, cstep
            if (n <= 22) then
               write(output_unit, '(a,i0,a)') &
                  '  Skipping D(n=', n, '): n must be > 22'
               cycle
            end if
            if (c < 1 .or. c > n - 13) then
               write(output_unit, '(a,i0,a,i0,a,i0,a)') &
                  '  Skipping D(n=', n, ', c=', c, &
                  '): c must be in [1, n-13=', n - 13, ']'
               cycle
            end if
            write(output_unit, '(/,a,i0,a,i0)') &
               'Matrix of class D(n,c):  n=', n, '  c=', c
            call mat_d%prepare(n, c)
            write(output_unit, '(a,i0)') &
               '  Non-zeros in original matrix: ', mat_d%nnz
            call mat_d%solve(aflag, iflag, ifail, elapsed, max_err)
            total_elapsed = total_elapsed + elapsed
            call print_solve_stats(aflag, iflag, ifail, elapsed, max_err)
         end do
      end do
   end if

   ! ---- Class E -----------------------------------------------------------
   ! Valid parameter range: n >= 3, 2 <= c <= n-1
   if (run_e) then
      do n = nstart, nend, nstep
         do c = cstart, cend, cstep
            if (n < 3) then
               write(output_unit, '(a,i0,a)') &
                  '  Skipping E(n=', n, '): n must be >= 3'
               cycle
            end if
            if (c < 2 .or. c > n - 1) then
               write(output_unit, '(a,i0,a,i0,a,i0,a)') &
                  '  Skipping E(n=', n, ', c=', c, &
                  '): c must be in [2, n-1=', n - 1, ']'
               cycle
            end if
            write(output_unit, '(/,a,i0,a,i0)') &
               'Matrix of class E(n,c):  n=', n, '  c=', c
            call mat_e%prepare(n, c)
            write(output_unit, '(a,i0)') &
               '  Non-zeros in original matrix: ', mat_e%nnz
            call mat_e%solve(aflag, iflag, ifail, elapsed, max_err)
            total_elapsed = total_elapsed + elapsed
            call print_solve_stats(aflag, iflag, ifail, elapsed, max_err)
         end do
      end do
   end if

   ! ---- Total elapsed wall-clock time ------------------------------------
   call system_clock(t_prog_end)
   write(output_unit, '(/,a,f12.3,a)') &
      'Total elapsed wall-clock time: ', &
      real(t_prog_end - t_prog_start, dp) / real(wall_rate, dp), ' s'

contains

   ! -----------------------------------------------------------------------
   !> Print the per-solve diagnostics in the same format as the original
   !> example drivers (maind.f, maine.f).
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
   !> Parse command-line arguments.  Unknown flags cause a help printout
   !> followed by a non-zero exit.
   subroutine parse_args(nstart, nstep, nend, cstart, cstep, cend, &
                         run_d, run_e)
      integer, intent(inout) :: nstart, nstep, nend
      integer, intent(inout) :: cstart, cstep, cend
      logical, intent(inout) :: run_d, run_e

      integer :: i, nargs
      character(len=64) :: arg, val

      nargs = command_argument_count()
      i = 1
      do while (i <= nargs)
         call get_command_argument(i, arg)
         select case (trim(arg))
         case ('--help', '-help', '-h')
            call print_help()
            stop

         case ('-n')
            i = i + 1
            if (i > nargs) then
               write(error_unit, '(a)') 'Error: -n requires an argument'
               stop 1
            end if
            call get_command_argument(i, val)
            call parse_range(val, nstart, nend, nstep)

         case ('-c')
            i = i + 1
            if (i > nargs) then
               write(error_unit, '(a)') 'Error: -c requires an argument'
               stop 1
            end if
            call get_command_argument(i, val)
            call parse_range(val, cstart, cend, cstep)

         case ('-class')
            i = i + 1
            if (i > nargs) then
               write(error_unit, '(a)') 'Error: -class requires an argument'
               stop 1
            end if
            call get_command_argument(i, val)
            select case (trim(val))
            case ('d')
               run_d = .true.
               run_e = .false.
            case ('e')
               run_d = .false.
               run_e = .true.
            case ('de', 'ed')
               run_d = .true.
               run_e = .true.
            case default
               write(error_unit, '(2a)') 'Error: unknown -class value: ', trim(val)
               stop 1
            end select

         case default
            write(error_unit, '(2a)') 'Error: unknown argument: ', trim(arg)
            call print_help()
            stop 1
         end select
         i = i + 1
      end do
   end subroutine parse_args

   ! -----------------------------------------------------------------------
   !> Parse a range string of the form:
   !>   VALUE            → vstart=VALUE, vend=VALUE, vstep=1
   !>   BEGIN:END        → vstart=BEGIN, vend=END,   vstep=1
   !>   BEGIN:END:STEP   → vstart=BEGIN, vend=END,   vstep=STEP
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
            write(error_unit, '(2a)') 'Error: cannot parse range start: ', trim(part1)
            stop 1
         end if
         read(part2, *, iostat=ios) vend
         if (ios /= 0) then
            write(error_unit, '(2a)') 'Error: cannot parse range end: ', trim(part2)
            stop 1
         end if
         vstep = 1
      else
         ! Three parts: BEGIN:END:STEP (pos2 is relative, make it absolute)
         pos2 = pos2 + pos1
         part1 = str(1:pos1 - 1)
         part2 = str(pos1 + 1:pos2 - 1)
         part3 = str(pos2 + 1:)
         read(part1, *, iostat=ios) vstart
         if (ios /= 0) then
            write(error_unit, '(2a)') 'Error: cannot parse range start: ', trim(part1)
            stop 1
         end if
         read(part2, *, iostat=ios) vend
         if (ios /= 0) then
            write(error_unit, '(2a)') 'Error: cannot parse range end: ', trim(part2)
            stop 1
         end if
         read(part3, *, iostat=ios) vstep
         if (ios /= 0) then
            write(error_unit, '(2a)') 'Error: cannot parse range step: ', trim(part3)
            stop 1
         end if
      end if
   end subroutine parse_range

   ! -----------------------------------------------------------------------
   !> Print the help message and return.
   subroutine print_help()
      write(output_unit, '(a)') &
         'Usage: timings_de [--help] [-class {d|e|de}] [-n RANGE] [-c RANGE]'
      write(output_unit, '(a)') ''
      write(output_unit, '(a)') &
         'Timing driver for sparse class D and E matrices (double precision).'
      write(output_unit, '(a)') &
         'Reproduces the benchmark table from Zlatev et al.'
      write(output_unit, '(a)') ''
      write(output_unit, '(a)') 'Options:'
      write(output_unit, '(a)') &
         '  --help            Show this help message and exit.'
      write(output_unit, '(a)') &
         '  -class {d|e|de}   Matrix classes to benchmark (default: de).'
      write(output_unit, '(a)') &
         '  -n RANGE          Matrix dimension range (default: 650:1000:50).'
      write(output_unit, '(a)') &
         '  -c RANGE          Structure parameter range (default: 4:204:40).'
      write(output_unit, '(a)') ''
      write(output_unit, '(a)') 'RANGE syntax:'
      write(output_unit, '(a)') &
         '  VALUE            single integer, e.g. 650'
      write(output_unit, '(a)') &
         '  BEGIN:END        step defaults to 1'
      write(output_unit, '(a)') &
         '  BEGIN:END:STEP   Fortran-slice style, e.g. 650:1000:50'
      write(output_unit, '(a)') ''
      write(output_unit, '(a)') 'Defaults (Zlatev benchmark):'
      write(output_unit, '(a)') &
         '  -class de  -n 650:1000:50  -c 4:204:40'
      write(output_unit, '(a)') ''
      write(output_unit, '(a)') 'Matrix descriptions:'
      write(output_unit, '(a)') &
         '  D(n,c): square, NNZ = 4*n + 55.'
      write(output_unit, '(a)') &
         '          Diagonal ones, three cyclic bands above diagonal,'
      write(output_unit, '(a)') &
         '          and a 10x10 dense triangle in the upper-right corner.'
      write(output_unit, '(a)') &
         '  E(n,c): symmetric SPD, NNZ = 5*n - 2*c - 2.'
      write(output_unit, '(a)') &
         '          Diagonal 4, side-diagonals -1, two off-diagonal'
      write(output_unit, '(a)') &
         '          bands at distance c with value -1.'
      write(output_unit, '(a)') &
         '          Similar to the five-point Laplacian on a 1-D grid.'
   end subroutine print_help

end program timings_de
