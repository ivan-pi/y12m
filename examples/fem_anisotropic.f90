! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: Gemini
!
! fem_anisotropic.f90 - 2D Anisotropic Diffusion via FEM
!
! Solves the equation:
!   -\nabla \cdot ( K \nabla u ) = f(x,y)   on  (0,1) x (0,1)
!
! Features:
!   - Bilinear (Q1) quadrilaterals on a structured grid.
!   - Independent NX and NY for testing rectangular elements.
!   - 2x2 Gauss Quadrature for exact stiffness integration.
!   - 3x3 Gauss Quadrature for highly accurate force integration.
!   - Comprehensive timing summary.
!
! Solver: y12ma (double precision sparse direct).

module fem_solver
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
  implicit none
  private
  public :: run

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = acos(-1.0_dp)

  ! Diffusion Tensor Components
  real(dp), parameter :: kxx = 10.0_dp
  real(dp), parameter :: kyy =  1.0_dp
  real(dp), parameter :: kxy =  3.0_dp

contains

  !> Exact manufactured solution
  elemental function u_exact(x, y) result(u)
    real(dp), intent(in) :: x, y
    real(dp) :: u
    u = exp(x) * sin(pi * y)
  end function u_exact

  !> Manufactured Right-Hand Side: -div(K grad u)
  elemental function f_rhs(x, y) result(f)
    real(dp), intent(in) :: x, y
    real(dp) :: f
    f = -exp(x) * ( (kxx - kyy * pi**2)*sin(pi * y) + 2.0_dp * kxy * pi * cos(pi * y) )
  end function f_rhs

  ! ---------------------------------------------------------------
  ! Root-mean-square error
  ! ---------------------------------------------------------------
  function rms_diff(a, b) result(rms)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: rms
    rms = norm2(a - b) / sqrt(real(size(a), dp))
  end function rms_diff

  ! ---------------------------------------------------------------
  ! Calculates the 4x4 element stiffness matrix Ke.
  ! ---------------------------------------------------------------
  subroutine get_element_stiffness(hx, hy, Ke)
    real(dp), intent(in) :: hx, hy
    real(dp), intent(out) :: Ke(4, 4)

    real(dp) :: xi(2), wt(2)
    real(dp) :: dNdxi(4), dNdeta(4), dNdx(4), dNdy(4)
    real(dp) :: detJ, w, inv_hx, inv_hy
    integer :: g1, g2, a, b

    ! 2-point Gauss rule
    xi = [-1.0_dp/sqrt(3.0_dp), 1.0_dp/sqrt(3.0_dp)]
    wt = [1.0_dp, 1.0_dp]

    detJ = (hx * hy) / 4.0_dp
    inv_hx = 2.0_dp / hx
    inv_hy = 2.0_dp / hy
    Ke = 0.0_dp

    do g1 = 1, 2
      do g2 = 1, 2
        w = wt(g1) * wt(g2) * detJ

        dNdxi(1) = -0.25_dp * (1.0_dp - xi(g2)); dNdeta(1) = -0.25_dp * (1.0_dp - xi(g1))
        dNdxi(2) =  0.25_dp * (1.0_dp - xi(g2)); dNdeta(2) = -0.25_dp * (1.0_dp + xi(g1))
        dNdxi(3) =  0.25_dp * (1.0_dp + xi(g2)); dNdeta(3) =  0.25_dp * (1.0_dp + xi(g1))
        dNdxi(4) = -0.25_dp * (1.0_dp + xi(g2)); dNdeta(4) =  0.25_dp * (1.0_dp - xi(g1))

        dNdx = dNdxi * inv_hx
        dNdy = dNdeta * inv_hy

        do a = 1, 4
          do b = 1, 4
            Ke(a,b) = Ke(a,b) + w * ( &
               dNdx(a) * (kxx * dNdx(b) + kxy * dNdy(b)) + &
               dNdy(a) * (kxy * dNdx(b) + kyy * dNdy(b)) )
          end do
        end do
      end do
    end do
  end subroutine get_element_stiffness

  ! ---------------------------------------------------------------
  ! FEM Assembly Routine
  ! ---------------------------------------------------------------
  subroutine assemble_system(NX, NY, hx, hy, a, rnr, snr, nz, F)
    integer, intent(in) :: NX, NY
    real(dp), intent(in) :: hx, hy
    real(dp), intent(out) :: a(:), F(:)
    integer, intent(out) :: rnr(:), snr(:)
    integer, intent(out) :: nz

    real(dp) :: Ke(4, 4), Fe(4)
    integer :: ndof_1d_x, ndof_1d_y, ndof
    integer :: i, j, a_idx, b_idx, row, col, elem_nodes(4)
    real(dp) :: x_base, y_base

    integer, allocatable :: row_nnz(:), row_cols(:,:)
    real(dp), allocatable :: row_vals(:,:)
    integer :: max_nz_per_row, r, c
    logical :: found

    real(dp) :: xi3(3), wt3(3), N_shape(4), x_g, y_g, detJ
    integer :: g1, g2

    ndof_1d_x = NX - 2
    ndof_1d_y = NY - 2
    ndof = ndof_1d_x * ndof_1d_y

    max_nz_per_row = 9
    allocate(row_nnz(ndof), row_cols(max_nz_per_row, ndof), row_vals(max_nz_per_row, ndof))
    row_nnz = 0
    row_vals = 0.0_dp
    F = 0.0_dp

    call get_element_stiffness(hx, hy, Ke)

    xi3 = [-sqrt(0.6_dp), 0.0_dp, sqrt(0.6_dp)]
    wt3 = [5.0_dp/9.0_dp, 8.0_dp/9.0_dp, 5.0_dp/9.0_dp]
    detJ = (hx * hy) / 4.0_dp

    do j = 1, NY - 1
      y_base = real(j - 1, dp) * hy
      do i = 1, NX - 1
        x_base = real(i - 1, dp) * hx

        elem_nodes(1) = (j - 1) * NX + i
        elem_nodes(2) = (j - 1) * NX + i + 1
        elem_nodes(3) = j * NX + i + 1
        elem_nodes(4) = j * NX + i

        Fe = 0.0_dp
        do g1 = 1, 3
          do g2 = 1, 3
            N_shape(1) = 0.25_dp * (1.0_dp - xi3(g1)) * (1.0_dp - xi3(g2))
            N_shape(2) = 0.25_dp * (1.0_dp + xi3(g1)) * (1.0_dp - xi3(g2))
            N_shape(3) = 0.25_dp * (1.0_dp + xi3(g1)) * (1.0_dp + xi3(g2))
            N_shape(4) = 0.25_dp * (1.0_dp - xi3(g1)) * (1.0_dp + xi3(g2))

            x_g = x_base + (1.0_dp + xi3(g1)) * hx / 2.0_dp
            y_g = y_base + (1.0_dp + xi3(g2)) * hy / 2.0_dp

            Fe(:) = Fe(:) + wt3(g1) * wt3(g2) * detJ * N_shape(:) * f_rhs(x_g, y_g)
          end do
        end do

        do a_idx = 1, 4
          row = get_dof_idx(elem_nodes(a_idx), NX, NY)
          if (row > 0) then
            F(row) = F(row) + Fe(a_idx)

            do b_idx = 1, 4
              col = get_dof_idx(elem_nodes(b_idx), NX, NY)
              if (col > 0) then
                found = .false.
                do c = 1, row_nnz(row)
                  if (row_cols(c, row) == col) then
                    row_vals(c, row) = row_vals(c, row) + Ke(a_idx, b_idx)
                    found = .true.; exit
                  end if
                end do
                if (.not. found) then
                  row_nnz(row) = row_nnz(row) + 1
                  row_cols(row_nnz(row), row) = col
                  row_vals(row_nnz(row), row) = Ke(a_idx, b_idx)
                end if
              else
                F(row) = F(row) - Ke(a_idx, b_idx) * u_exact( &
                   get_x(elem_nodes(b_idx), NX, hx), get_y(elem_nodes(b_idx), NX, hy))
              end if
            end do
          end if
        end do
      end do
    end do

    nz = 0
    do r = 1, ndof
      do c = 1, row_nnz(r)
        nz = nz + 1
        rnr(nz) = r
        snr(nz) = row_cols(c, r)
        a(nz)   = row_vals(c, r)
      end do
    end do
  end subroutine assemble_system

  pure function get_dof_idx(node, NX, NY) result(dof)
    integer, intent(in) :: node, NX, NY
    integer :: dof, ix, iy
    iy = (node - 1) / NX + 1
    ix = mod(node - 1, NX) + 1
    if (ix == 1 .or. ix == NX .or. iy == 1 .or. iy == NY) then
      dof = -1
    else
      dof = (iy - 2) * (NX - 2) + (ix - 2) + 1
    end if
  end function get_dof_idx

  pure function get_x(node, NX, hx) result(x)
    integer, intent(in) :: node, NX
    real(dp), intent(in) :: hx
    real(dp) :: x
    x = real(mod(node - 1, NX), dp) * hx
  end function get_x

  pure function get_y(node, NX, hy) result(y)
    integer, intent(in) :: node, NX
    real(dp), intent(in) :: hy
    real(dp) :: y
    y = real((node - 1) / NX, dp) * hy
  end function get_y

  ! ---------------------------------------------------------------
  ! Write Output
  ! ---------------------------------------------------------------
  subroutine write_output(NX, NY, u_num, u_ex, filename, header, nhr)
    integer, intent(in) :: NX, NY
    real(dp), intent(in) :: u_num(NX, NY), u_ex(NX, NY)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nhr
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
        write(funit, datafmt) real(i - 1, dp) * hx, real(j - 1, dp) * hy, u_num(i, j), u_ex(i, j)
      end do
      write(funit, *)
    end do
    close(funit)
  end subroutine write_output

  ! ---------------------------------------------------------------
  ! Main Driver
  ! ---------------------------------------------------------------
  subroutine run(NX, NY, outfile)
    use y12m, only: y12ma
    integer, intent(in) :: NX, NY
    character(len=*), intent(in) :: outfile

    integer :: ndof_1d_x, ndof_1d_y, ndof, nz_max, nz, nn, nn1, iha
    real(dp) :: hx, hy
    real(dp), allocatable :: F(:), a(:), pivot(:)
    integer, allocatable :: snr(:), rnr(:), ha(:,:)
    real(dp) :: aflag(8), err_max, err_rms
    integer :: iflag(10), ifail, i, j, k
    real(dp), allocatable :: u_full(:,:), u_ex_full(:,:), u_interior(:), u_ex_int(:)

    ! Timing variables
    integer :: t0, t1, clock_rate
    integer :: t_assemble, t_solve, t_error, t_output
    real(dp) :: s_assm, s_solv, s_err, s_out, s_tot

    ndof_1d_x = NX - 2
    ndof_1d_y = NY - 2
    ndof      = ndof_1d_x * ndof_1d_y
    hx        = 1.0_dp / real(NX - 1, dp)
    hy        = 1.0_dp / real(NY - 1, dp)

    write(output_unit, '(a,i0,a,i0,a,i0,a)') 'Grid: ', NX, ' x ', NY, '  (', ndof, ' interior DOFs)'

    nz_max = 9 * ndof
    nn     = max(3 * nz_max, 2 * (max(ndof_1d_x, ndof_1d_y) + 2) * ndof)
    nn1    = nn
    iha    = ndof
    allocate(a(nn), pivot(ndof), F(ndof), snr(nn), rnr(nn1), ha(ndof, 11))
    allocate(u_full(NX, NY), u_ex_full(NX, NY), u_interior(ndof), u_ex_int(ndof))

    call system_clock(count_rate=clock_rate)

    ! 1. Assemble
    call system_clock(t0)
    call assemble_system(NX, NY, hx, hy, a, rnr, snr, nz, F)
    call system_clock(t1); t_assemble = t1 - t0

    ! 2. Solve with y12ma
    call system_clock(t0)
    ifail = 0
    call y12ma(ndof, nz, a, snr, nn, rnr, nn1, pivot, ha, iha, aflag, iflag, F, ifail)
    if (ifail /= 0) then
      write(error_unit, '(a,i0)') 'ERROR: y12ma returned IFAIL = ', ifail
      stop 1
    end if
    call system_clock(t1); t_solve = t1 - t0

    ! 3. Map back and evaluate errors
    call system_clock(t0)
    u_full = 0.0_dp
    do j = 1, NY
      do i = 1, NX
        u_ex_full(i, j) = u_exact(real(i - 1, dp) * hx, real(j - 1, dp) * hy)
        if (i == 1 .or. i == NX .or. j == 1 .or. j == NY) then
           u_full(i,j) = u_ex_full(i,j)
        else
           k = (j - 2) * ndof_1d_x + (i - 2) + 1
           u_full(i, j) = F(k)
        end if
      end do
    end do

    k = 0
    do j = 2, NY - 1
      do i = 2, NX - 1
        k = k + 1
        u_interior(k) = u_full(i, j)
        u_ex_int(k)   = u_ex_full(i, j)
      end do
    end do

    err_max = maxval(abs(u_interior - u_ex_int))
    err_rms = rms_diff(u_interior, u_ex_int)
    call system_clock(t1); t_error = t1 - t0

    write(output_unit, '(a,es12.4)') 'Max pointwise error  : ', err_max
    write(output_unit, '(a,es12.4)') 'RMS pointwise error  : ', err_rms

    ! 4. Write data file
    call system_clock(t0)
    block
      character(len=70) :: hdr(2)
      hdr(1) = 'FEM Q1 Anisotropic Diffusion on unit square'
      hdr(2) = 'MMS: u(x,y) = exp(x) * sin(pi * y)'
      call write_output(NX, NY, u_full, u_ex_full, outfile, hdr, size(hdr))
    end block
    call system_clock(t1); t_output = t1 - t0

    ! 5. Profiling Calculations
    s_assm = real(t_assemble, dp) / real(clock_rate, dp)
    s_solv = real(t_solve, dp) / real(clock_rate, dp)
    s_err  = real(t_error, dp) / real(clock_rate, dp)
    s_out  = real(t_output, dp) / real(clock_rate, dp)
    s_tot  = s_assm + s_solv + s_err + s_out

    if (s_tot == 0.0_dp) s_tot = 1e-12_dp

    write(output_unit, '(/,a)') 'Timing summary:'
    write(output_unit, '(a,g14.6,a,f7.2,a)') '  Assembly     : ', s_assm, ' s  (', (s_assm/s_tot)*100.0_dp, ' %)'
    write(output_unit, '(a,g14.6,a,f7.2,a)') '  Solve        : ', s_solv, ' s  (', (s_solv/s_tot)*100.0_dp, ' %)'
    write(output_unit, '(a,g14.6,a,f7.2,a)') '  Error        : ', s_err,  ' s  (', (s_err /s_tot)*100.0_dp, ' %)'
    write(output_unit, '(a,g14.6,a,f7.2,a)') '  Output       : ', s_out,  ' s  (', (s_out /s_tot)*100.0_dp, ' %)'
    write(output_unit, '(a,g14.6,a,f7.2,a)') '  Total Time   : ', s_tot,  ' s  (', 100.00_dp, ' %)'

  end subroutine run
end module fem_solver

! ==============================================================
! Main program
! ==============================================================
program fem_anisotropic
  use fem_solver, only: run
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none

  integer :: NX, NY, ios, iarg, nargs, iarg_pos
  character(len=256) :: arg, outfile

  NX = 22
  NY = 22
  outfile = 'fem_anisotropic.dat'
  nargs = command_argument_count()

  iarg = 1
  iarg_pos = 1

  do while (iarg <= nargs)
    call get_command_argument(iarg, arg)

    if (trim(arg) == '--help' .or. trim(arg) == '-h') then
      write(*, '(a)') 'Usage: fem_anisotropic [--help] [NX] [NY] [output_file]'
      write(*, '(a)') '  NX, NY       total grid size in X and Y (>= 3)'
      write(*, '(a)') '  output_file  output data file'
      stop 0
    else
      ! Positional Arguments
      if (iarg_pos == 1) then
        read(arg, *, iostat=ios) NX
        if (ios /= 0 .or. NX < 3) then
          write(error_unit, '(a)') 'ERROR: NX must be an integer >= 3'
          stop 1
        end if
        ! If only one N is provided, default NY to match NX (square grid backward compatibility)
        NY = NX
        iarg_pos = 2
      else if (iarg_pos == 2) then
        read(arg, *, iostat=ios) NY
        if (ios /= 0 .or. NY < 3) then
          ! If it's not a number, it must be the output file
          outfile = trim(arg)
          iarg_pos = 4 ! Skip NY
        else
          iarg_pos = 3
        end if
      else if (iarg_pos == 3) then
        outfile = trim(arg)
        iarg_pos = 4
      end if
    end if
    iarg = iarg + 1
  end do

  call run(NX, NY, trim(outfile))
end program fem_anisotropic