! SPDX-License-Identifier: GPL-2.0-only
!
! y12m_example_helpers.f90 — Shared helper routines for the y12m example programs.
!
! Provides:
!   parse_range        — CLI range parser (VALUE | BEGIN:END | BEGIN:END:STEP)
!   rms_diff           — root-mean-square of element-wise differences (1D and 2D)
!   write_field_output — write a 2D numerical + exact solution to a gnuplot file
!   write_three_column_output — write 1D x, numerical, exact columns

module y12m_example_helpers
   use, intrinsic :: iso_fortran_env, only: error_unit
   implicit none
   private

   integer, parameter :: dp = kind(1.0d0)

   public :: parse_range
   public :: rms_diff
   public :: write_field_output
   public :: write_three_column_output
   public :: compress_matrix

   !> Root-mean-square of element-wise differences.
   !> Overloaded for 1-D and 2-D assumed-shape real(dp) arrays.
   interface rms_diff
      module procedure rms_diff_1d
      module procedure rms_diff_2d
   end interface rms_diff

contains

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
         pos2  = pos2 + pos1
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
   !> Root-mean-square of element-wise difference of two 1-D arrays.
   pure function rms_diff_1d(a, b) result(rms)
      real(dp), intent(in) :: a(:), b(:)
      real(dp) :: rms
      rms = norm2(a - b) / sqrt(real(size(a), dp))
   end function rms_diff_1d


   ! -----------------------------------------------------------------------
   !> Root-mean-square of element-wise difference of two 2-D arrays.
   pure function rms_diff_2d(a, b) result(rms)
      real(dp), intent(in) :: a(:,:), b(:,:)
      real(dp) :: rms
      rms = norm2(a - b) / sqrt(real(size(a), dp))
   end function rms_diff_2d


   ! -----------------------------------------------------------------------
   !> Write a 2-D numerical solution and exact solution to a gnuplot-compatible
   !> pm3d data file.
   !>
   !> The output has NX * NY data lines grouped in NY blocks separated by blank
   !> lines (one blank line per j-row, as expected by gnuplot's pm3d).
   !>
   !> Arguments:
   !>   NX, NY    — grid sizes in x and y (node counts, including boundaries)
   !>   u_num     — numerical solution, shape (NX, NY)
   !>   u_ex      — exact (reference) solution, shape (NX, NY)
   !>   filename  — output file name
   !>   header    — comment lines written at the top (prefix "# " is added)
   !>   nhr       — number of header lines
   subroutine write_field_output(NX, NY, u_num, u_ex, filename, header, nhr)
      integer, intent(in) :: NX, NY, nhr
      real(dp), intent(in) :: u_num(:,:), u_ex(:,:)
      character(len=*), intent(in) :: filename
      character(len=70), intent(in) :: header(nhr)
      character(len=*), parameter :: datafmt = '(4(1x,es14.6))'
      real(dp) :: hx, hy
      integer :: funit, i, j, ih

      hx = 1.0_dp / real(NX - 1, dp)
      hy = 1.0_dp / real(NY - 1, dp)

      open(newunit=funit, file=filename, status='unknown', action='write')
      do ih = 1, nhr
         write(funit, '("# ",a)') trim(header(ih))
      end do
      write(funit, '(a)') '# Columns: x   y   u_numerical   u_exact'
      do j = 1, NY
         do i = 1, NX
            write(funit, datafmt) real(i - 1, dp) * hx, real(j - 1, dp) * hy, &
                u_num(i, j), u_ex(i, j)
         end do
         write(funit, *)
      end do
      close(funit)
   end subroutine write_field_output

   ! -----------------------------------------------------------------------
   !> Write a three-column 1-D dataset to a text file.
   !>
   !> Arguments:
   !>   x, y_num, y_ex  — vectors of equal length
   !>   filename        — output file name
   !>   header          — comment lines written at the top (prefix "# " is added)
   !>   nhr             — number of header lines
   !>   columns         — column description line (written with "# " prefix)
   subroutine write_three_column_output(x, y_num, y_ex, filename, header, nhr, columns)
      real(dp), intent(in) :: x(:), y_num(:), y_ex(:)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nhr
      character(len=*), intent(in) :: header(nhr)
      character(len=*), intent(in) :: columns

      integer :: file_unit, idx, header_idx

      if (size(x) /= size(y_num) .or. size(x) /= size(y_ex)) then
         write(error_unit, '(a)') 'Error: write_three_column_output size mismatch'
         stop 1
      end if

      open(newunit=file_unit, file=filename, status='unknown', action='write')
      do header_idx = 1, nhr
         write(file_unit, '("# ",a)') trim(header(header_idx))
      end do
      write(file_unit, '("# ",a)') trim(columns)

      do idx = 1, size(x)
         write(file_unit, '(3(1x,es14.6))') x(idx), y_num(idx), y_ex(idx)
      end do
      close(file_unit)
   end subroutine write_three_column_output
 
   ! -----------------------------------------------------------------------
   !> Sort and merge duplicate (row,col) entries in a COO sparse matrix.
   !>
   !> On entry, nz COO triplets (rnr, snr, a) may contain duplicate (i,j)
   !> positions.  On exit the triplets are sorted row-major and duplicates
   !> are summed; nz is updated to the compressed count.
   !>
   !> The arrays rnr, snr, a must be allocated to at least nz elements on
   !> entry (the in-place operation does not grow them).
   !>
   !> Note: y12ma returns IFAIL=11 when two entries share the same (i,j)
   !> position, so this routine must be called before passing a FEM-assembled
   !> COO matrix to y12ma.
   subroutine compress_matrix(nz, rnr, snr, a)
      integer, intent(inout) :: nz
      integer, intent(inout) :: rnr(:), snr(:)
      real(dp), intent(inout) :: a(:)

      integer :: write_ptr, i

      if (nz <= 1) return

      call quicksort_coo(1, nz, rnr, snr, a)

      write_ptr = 1
      do i = 2, nz
         if (rnr(i) == rnr(write_ptr) .and. snr(i) == snr(write_ptr)) then
            a(write_ptr) = a(write_ptr) + a(i)
         else
            write_ptr = write_ptr + 1
            rnr(write_ptr) = rnr(i)
            snr(write_ptr) = snr(i)
            a(write_ptr)   = a(i)
         end if
      end do

      nz = write_ptr

   contains

      recursive subroutine quicksort_coo(low, high, rnr, snr, a)
         integer, intent(in) :: low, high
         integer, intent(inout) :: rnr(:), snr(:)
         real(dp), intent(inout) :: a(:)

         integer :: i, j, pivot_r, pivot_s, temp_i
         real(dp) :: temp_a

         if (low >= high) return

         pivot_r = rnr(high)
         pivot_s = snr(high)
         i = low - 1

         do j = low, high - 1
            if (rnr(j) < pivot_r .or. &
                (rnr(j) == pivot_r .and. snr(j) < pivot_s)) then
               i = i + 1
               temp_i = rnr(i)
               rnr(i) = rnr(j)
               rnr(j) = temp_i
               temp_i = snr(i)
               snr(i) = snr(j)
               snr(j) = temp_i
               temp_a = a(i)
               a(i)   = a(j)
               a(j)   = temp_a
            end if
         end do

         i = i + 1
         temp_i = rnr(i)
         rnr(i) = rnr(high)
         rnr(high) = temp_i
         temp_i = snr(i)
         snr(i) = snr(high)
         snr(high) = temp_i
         temp_a = a(i)
         a(i)   = a(high)
         a(high) = temp_a

         call quicksort_coo(low, i - 1, rnr, snr, a)
         call quicksort_coo(i + 1, high, rnr, snr, a)
      end subroutine quicksort_coo

   end subroutine compress_matrix

end module y12m_example_helpers
