! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: Gemini
!
! biharmonic_solver.f90 - Kirchhoff Plate (Biharmonic) on a unit square
!
! Solves the 2-D Biharmonic equation
!
!   \nabla^4 u = f(x,y)   on  (0,1) x (0,1)
!
! with Simply Supported boundary conditions:
!   u = 0          on boundary
!   \nabla^2 u = 0 on boundary
!
! Discretisation
! --------------
! Total grid N x N (including boundary nodes), spacing h = 1/(N-1).
! Interior nodes: i,j = 2..N-1 (1-indexed); ndof = (N-2)^2.
!
! Features command-line toggles for MMS difficulty and Stencil width.
!
! Solver: y12ma (double precision), no iterative refinement.

! ==============================================================
! Module biharmonic_solver
! ==============================================================
module biharmonic_solver
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
  implicit none
  private
  public :: run

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = acos(-1.0_dp)

  ! 13-point stencil for \nabla^4
  real(dp), parameter :: stencil_13(-2:2, -2:2) = reshape([ &
      0.0_dp,  0.0_dp,  1.0_dp,  0.0_dp,  0.0_dp, &
      0.0_dp,  2.0_dp, -8.0_dp,  2.0_dp,  0.0_dp, &
      1.0_dp, -8.0_dp, 20.0_dp, -8.0_dp,  1.0_dp, &
      0.0_dp,  2.0_dp, -8.0_dp,  2.0_dp,  0.0_dp, &
      0.0_dp,  0.0_dp,  1.0_dp,  0.0_dp,  0.0_dp  &
  ], [5, 5])

  ! 21-point stencil for \nabla^4
  ! Derived from the convolution of the 5-point and 9-point Laplacians.
  ! Multiplied by 1/3 to match the correct scaling of the Biharmonic operator.
  real(dp), parameter :: stencil_21(-2:2, -2:2) = reshape([ &
      0.0_dp,  1.0_dp,  1.0_dp,  1.0_dp,  0.0_dp, &
      1.0_dp, -2.0_dp,-10.0_dp, -2.0_dp,  1.0_dp, &
      1.0_dp,-10.0_dp, 36.0_dp,-10.0_dp,  1.0_dp, &
      1.0_dp, -2.0_dp,-10.0_dp, -2.0_dp,  1.0_dp, &
      0.0_dp,  1.0_dp,  1.0_dp,  1.0_dp,  0.0_dp  &
  ], [5, 5]) / 3.0_dp

contains

  !> Exact manufactured solution (toggle between easy/hard)
  elemental function u_exact(x, y, mms_level) result(u)
    real(dp), intent(in) :: x, y
    integer, intent(in) :: mms_level
    real(dp) :: u
    if (mms_level == 1) then
      u = sin(pi * x) * sin(pi * y)
    else
      u = sin(3.0_dp * pi * x) * sin(4.0_dp * pi * y)
    end if
  end function u_exact

  !> Manufactured Right-Hand Side ( \nabla^4 u )
  elemental function f_rhs(x, y, mms_level) result(f)
    real(dp), intent(in) :: x, y
    integer, intent(in) :: mms_level
    real(dp) :: f
    if (mms_level == 1) then
      f = 4.0_dp * (pi**4) * sin(pi * x) * sin(pi * y)
    else
      f = 625.0_dp * (pi**4) * sin(3.0_dp * pi * x) * sin(4.0_dp * pi * y)
    end if
  end function f_rhs

  ! ---------------------------------------------------------------
  ! Sparse matrix-vector product: y <- y + alpha * A * x  (COO)
  ! ---------------------------------------------------------------
  subroutine dp_coo_gemv(m, nz, alpha, a, ia, ja, x, y)
    integer, intent(in) :: m, nz
    real(dp), intent(in) :: alpha
    real(dp), intent(in) :: a(nz)
    integer, intent(in) :: ia(nz), ja(nz)
    real(dp), intent(in) :: x(m)
    real(dp), intent(inout) :: y(m)
    integer :: p
    do p = 1, nz
      y(ia(p)) = y(ia(p)) + alpha * a(p) * x(ja(p))
    end do
  end subroutine dp_coo_gemv

  ! ---------------------------------------------------------------
  ! Root-mean-square of element-wise difference of two 2-D arrays
  ! ---------------------------------------------------------------
  function rms_diff(a, b) result(rms)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp) :: rms
    rms = norm2(a - b) / sqrt(real(size(a), dp))
  end function rms_diff

  ! ---------------------------------------------------------------
  ! Assemble the stiffness matrix A in COO format.
  ! ---------------------------------------------------------------
  subroutine assemble_matrix(N, stencil_type, a, rnr, snr, nz)
    integer, intent(in) :: N, stencil_type
    real(dp), intent(out) :: a(:)
    integer, intent(out) :: rnr(:), snr(:)
    integer, intent(out) :: nz

    integer :: i, j, is, js, ni, nj, k_row, k_col, ndof_1d
    real(dp) :: weight, final_weight
    real(dp) :: active_stencil(-2:2, -2:2)

    integer :: cols(25), n_in_row, c
    real(dp) :: vals(25)
    logical :: found

    if (stencil_type == 21) then
      active_stencil = stencil_21
    else
      active_stencil = stencil_13
    end if

    ndof_1d = N - 2
    nz = 0

    do j = 2, N - 1
      do i = 2, N - 1
        k_row = (j - 2) * ndof_1d + (i - 2) + 1
        n_in_row = 0

        do js = -2, 2
          do is = -2, 2
            weight = active_stencil(is, js)
            if (weight == 0.0_dp) cycle

            ni = i + is
            nj = j + js
            final_weight = weight

            ! Ghost node mirroring for Simply Supported bounds
            if (ni == 0) then
              ni = 2; final_weight = -final_weight
            else if (ni == N + 1) then
              ni = N - 1; final_weight = -final_weight
            end if

            if (nj == 0) then
              nj = 2; final_weight = -final_weight
            else if (nj == N + 1) then
              nj = N - 1; final_weight = -final_weight
            end if

            if (ni >= 2 .and. ni <= N - 1 .and. nj >= 2 .and. nj <= N - 1) then
              k_col = (nj - 2) * ndof_1d + (ni - 2) + 1

              found = .false.
              do c = 1, n_in_row
                if (cols(c) == k_col) then
                  vals(c) = vals(c) + final_weight
                  found = .true.
                  exit
                end if
              end do

              if (.not. found) then
                n_in_row = n_in_row + 1
                cols(n_in_row) = k_col
                vals(n_in_row) = final_weight
              end if
            end if
          end do
        end do

        do c = 1, n_in_row
          nz = nz + 1
          rnr(nz) = k_row
          snr(nz) = cols(c)
          a(nz) = vals(c)
        end do
      end do
    end do
  end subroutine assemble_matrix

  ! ---------------------------------------------------------------
  ! Compute the right-hand side from MMS f(x,y).
  ! ---------------------------------------------------------------
  subroutine compute_rhs(N, mms_level, b)
    integer, intent(in) :: N, mms_level
    real(dp), intent(inout) :: b(:)
    integer :: i, j, k_row, ndof_1d
    real(dp) :: h, x, y

    h = 1.0_dp / real(N - 1, dp)
    ndof_1d = N - 2

    do j = 2, N - 1
      y = real(j - 1, dp) * h
      do i = 2, N - 1
        x = real(i - 1, dp) * h
        k_row = (j - 2) * ndof_1d + (i - 2) + 1
        b(k_row) = (h**4) * f_rhs(x, y, mms_level)
      end do
    end do
  end subroutine compute_rhs

  ! ---------------------------------------------------------------
  ! Write the solution to a gnuplot pm3d data file.
  ! ---------------------------------------------------------------
  subroutine write_output(N, u_num, u_ex, filename, header, nhr)
    integer, intent(in) :: N
    real(dp), intent(in) :: u_num(N, N), u_ex(N, N)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nhr
    character(len=70), intent(in) :: header(nhr)
    character(len=*), parameter :: datafmt = '(4(1x,es14.6))'
    real(dp) :: h
    integer :: funit, i, j, ih

    h = 1.0_dp / real(N - 1, dp)
    open(newunit=funit, file=filename, status='unknown', action='write')
    do ih = 1, nhr
      write(funit, '("# ",a)') trim(header(ih))
    end do
    write(funit, '(a)') '# Columns: x   y   u_numerical   u_exact'
    write(funit, '(a)') '#'
    do j = 1, N
      do i = 1, N
        write(funit, datafmt) real(i - 1, dp) * h, real(j - 1, dp) * h, &
            u_num(i, j), u_ex(i, j)
      end do
      write(funit, *)
    end do
    close(funit)
  end subroutine write_output

  ! ---------------------------------------------------------------
  ! Main Driver
  ! ---------------------------------------------------------------
  subroutine run(N, outfile, stencil_type, mms_level)
    use y12m, only: y12ma
    integer, intent(in) :: N, stencil_type, mms_level
    character(len=*), intent(in) :: outfile

    integer :: ndof_1d, ndof, nz_max, nz, nn, nn1, iha
    real(dp) :: h
    real(dp), allocatable :: b(:), a(:), pivot(:), a0(:), b0(:), residual(:)
    integer, allocatable :: snr(:), rnr(:), ha(:,:), snr0(:), rnr0(:)
    real(dp) :: aflag(8), err_max, err_rms, res_inf, res_2
    integer :: iflag(10), ifail, i, j, k
    real(dp), allocatable :: u_full(:,:), u_ex_full(:,:)

    ! Timing variables
    integer :: t0, t1, clock_rate
    integer :: t_assemble, t_solve, t_error, t_output
    real(dp) :: s_assm, s_solv, s_err, s_out, s_tot

    ndof_1d = N - 2
    ndof    = ndof_1d * ndof_1d
    h       = 1.0_dp / real(N - 1, dp)

    write(output_unit, '(a,i0,a,i0,a,i0,a)') &
        'Grid: ', N, ' x ', N, '  (', ndof, ' interior DOFs)'
    write(output_unit, '(a,i0,a)') 'Stencil    : ', stencil_type, '-point'
    write(output_unit, '(a,i0)')   'MMS Level  : ', mms_level

    ! Allocate workspace
    nz_max = stencil_type * ndof
    nn     = max(3 * nz_max, 2 * (ndof_1d + 2) * ndof)
    nn1    = nn
    iha    = ndof
    allocate(a(nn), pivot(ndof), b(ndof), snr(nn), rnr(nn1), ha(ndof, 11))
    allocate(a0(nz_max), b0(ndof), snr0(nz_max), rnr0(nz_max), residual(ndof))
    allocate(u_full(N, N), u_ex_full(N, N))

    call system_clock(count_rate=clock_rate)

    ! 1. Assemble A u = b
    call system_clock(t0)
    b  = 0.0_dp
    call assemble_matrix(N, stencil_type, a, rnr, snr, nz)
    call compute_rhs(N, mms_level, b)

    a0(1:nz) = a(1:nz)
    b0 = b
    snr0(1:nz) = snr(1:nz)
    rnr0(1:nz) = rnr(1:nz)
    call system_clock(t1); t_assemble = t1 - t0

    ! 2. Solve with y12ma
    call system_clock(t0)
    ifail = 0
    call y12ma(ndof, nz, a, snr, nn, rnr, nn1, pivot, ha, iha, &
        aflag, iflag, b, ifail)
    if (ifail /= 0) then
      write(error_unit, '(a,i0)') 'ERROR: y12ma returned IFAIL = ', ifail
      stop 1, quiet = .true.
    end if
    call system_clock(t1); t_solve = t1 - t0

    ! 3. Error vs exact solution and discrete residual
    call system_clock(t0)
    u_full = 0.0_dp
    do j = 2, N - 1
      do i = 2, N - 1
        k = (j - 2) * ndof_1d + (i - 2) + 1
        u_full(i, j) = b(k)
      end do
    end do

    do j = 1, N
      do i = 1, N
        u_ex_full(i, j) = u_exact(real(i - 1, dp) * h, real(j - 1, dp) * h, mms_level)
      end do
    end do

    err_max = maxval(abs(u_full(2:N-1, 2:N-1) - u_ex_full(2:N-1, 2:N-1)))
    err_rms = rms_diff(u_full(2:N-1, 2:N-1), u_ex_full(2:N-1, 2:N-1))

    residual = b0
    call dp_coo_gemv(ndof, nz, alpha=-1.0_dp, a=a0, ia=rnr0, ja=snr0, x=b, y=residual)
    res_inf = maxval(abs(residual))
    res_2   = norm2(residual)
    call system_clock(t1); t_error = t1 - t0

    write(output_unit, '(a,es12.4)') 'Max pointwise error  : ', err_max
    write(output_unit, '(a,es12.4)') 'RMS pointwise error  : ', err_rms
    write(output_unit, '(a,es12.4)') 'Residual ||r||_inf   : ', res_inf
    write(output_unit, '(a,es12.4)') 'Residual ||r||_2     : ', res_2

    ! 4. Write data file
    call system_clock(t0)
    block
      character(len=70) :: hdr(2)
      if (stencil_type == 21) then
        hdr(1) = '21-point biharmonic stencil, Kirchhoff plate on unit square'
      else
        hdr(1) = '13-point biharmonic stencil, Kirchhoff plate on unit square'
      end if
      if (mms_level == 1) then
        hdr(2) = 'MMS Easy: u(x,y) = sin(pi x) sin(pi y)'
      else
        hdr(2) = 'MMS Hard: u(x,y) = sin(3pi x) sin(4pi y)'
      end if
      call write_output(N, u_full, u_ex_full, outfile, hdr, size(hdr))
    end block
    call system_clock(t1); t_output = t1 - t0

    ! 5. Profiling Calculations & Formatting
    s_assm = real(t_assemble, dp) / real(clock_rate, dp)
    s_solv = real(t_solve, dp) / real(clock_rate, dp)
    s_err  = real(t_error, dp) / real(clock_rate, dp)
    s_out  = real(t_output, dp) / real(clock_rate, dp)
    s_tot  = s_assm + s_solv + s_err + s_out

    ! Fallback to prevent divide-by-zero on extremely fast / small grid runs
    if (s_tot == 0.0_dp) s_tot = 1e-12_dp

    write(output_unit, '(/,a)') 'Timing summary:'
    write(output_unit, '(a,g14.6,a,f7.2,a)') '  Assembly     : ', s_assm, ' s  (', (s_assm/s_tot)*100.0_dp, ' %)'
    write(output_unit, '(a,g14.6,a,f7.2,a)') '  Solve        : ', s_solv, ' s  (', (s_solv/s_tot)*100.0_dp, ' %)'
    write(output_unit, '(a,g14.6,a,f7.2,a)') '  Error+resid  : ', s_err,  ' s  (', (s_err /s_tot)*100.0_dp, ' %)'
    write(output_unit, '(a,g14.6,a,f7.2,a)') '  Output       : ', s_out,  ' s  (', (s_out /s_tot)*100.0_dp, ' %)'
    write(output_unit, '(a,g14.6,a,f7.2,a)') '  Total Time   : ', s_tot,  ' s  (', 100.00_dp, ' %)'

  end subroutine run
end module biharmonic_solver

! ==============================================================
! Main program: parse CLI, call run.
! ==============================================================
program biharmonic_driver
  use biharmonic_solver, only: run
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none

  integer :: N, ios, iarg, nargs, iarg_pos
  character(len=256) :: arg, outfile
  integer :: stencil_type, mms_level

  ! Defaults
  N            = 22
  outfile      = 'biharmonic.dat'
  stencil_type = 13
  mms_level    = 1
  nargs        = command_argument_count()

  iarg = 1
  iarg_pos = 1

  do while (iarg <= nargs)
    call get_command_argument(iarg, arg)

    if (trim(arg) == '--help' .or. trim(arg) == '-h') then
      write(*, '(a)') 'Usage: biharmonic_driver [--help] [--stencil 13|21] [--mms easy|hard] [N] [output_file]'
      stop 0
    else if (trim(arg) == '--stencil') then
      iarg = iarg + 1
      call get_command_argument(iarg, arg)
      read(arg, *, iostat=ios) stencil_type
      if (ios /= 0 .or. (stencil_type /= 13 .and. stencil_type /= 21)) then
        write(error_unit, '(a)') 'ERROR: --stencil must be 13 or 21.'
        stop 1
      end if
    else if (trim(arg) == '--mms') then
      iarg = iarg + 1
      call get_command_argument(iarg, arg)
      if (trim(arg) == 'hard') then
        mms_level = 2
      else
        mms_level = 1
      end if
    else
      ! Positional Arguments
      if (iarg_pos == 1) then
        read(arg, *, iostat=ios) N
        if (ios /= 0 .or. N < 3) then
          write(error_unit, '(a)') 'ERROR: N must be an integer >= 3'
          stop 1
        end if
        iarg_pos = 2
      else if (iarg_pos == 2) then
        outfile = trim(arg)
        iarg_pos = 3
      end if
    end if

    iarg = iarg + 1
  end do

  call run(N, trim(outfile), stencil_type, mms_level)
end program biharmonic_driver