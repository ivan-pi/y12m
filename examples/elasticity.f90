! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: Gemini
!
! elasticity_2d.f90 - 2D Linear Elasticity (Plane Stress) on a Cantilever Beam
!
! Solves the Navier-Cauchy equations using Q1 Finite Elements:
!   mu * Lap(u) + (lam + mu) * Grad(Div(u)) = 0
!
! Exact Analytical Solution (Timoshenko Beam):
!   Evaluated to prescribe Dirichlet boundaries and measure L2 errors.
!
! =========================================================================
! DESIGN NOTE: Short-Edge (Y-Axis First) Node Numbering
! =========================================================================
! Direct sparse solvers like Y12M rely on Markowitz pivoting for stability,
! but they lack global topological re-ordering algorithms (like Reverse
! Cuthill-McKee or AMD). Therefore, the initial bandwidth of the input
! matrix strictly dictates how much memory the LU fill-in will consume.
!
! If we number nodes along the long edge first (X-axis, NX=240), two vertically
! adjacent nodes will have IDs separated by 240. Because each node has 2 DOFs,
! the matrix bandwidth becomes at least 480. During Gaussian elimination, this
! wide front causes massive fill-in, quickly exhausting RAM and throwing IFAIL=5.
!
! By numbering nodes along the short edge first (Y-axis, NY=40), the maximum
! ID difference between any connected nodes drops to 40 (bandwidth ~80).
! This topological trick tightly clusters the non-zeros around the diagonal,
! exponentially reducing LU fill-in and allowing massive grids to be solved.
!
! Formula mapping (i,j) to node ID: node = (i - 1) * NY + j
! =========================================================================
!
module elasticity_solver
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
  implicit none
  private
  public :: run

  integer, parameter :: dp = kind(1.0d0)

contains

  ! === Analytical Solutions for Plane Stress Cantilever ===
  elemental function exact_ux(x, y, P, L, D, E, nu, Inertia) result(ux)
    real(dp), intent(in) :: x, y, P, L, D, E, nu, Inertia
    real(dp) :: ux
    ux = -P * y / (6.0_dp * E * Inertia) * &
         (3.0_dp * (x**2 - L**2) - (2.0_dp + nu) * y**2 + 6.0_dp * (1.0_dp + nu) * (D/2.0_dp)**2)
  end function exact_ux

  ! Note: Depth (D) removed from signature as it is not used in the y-displacement formula
  elemental function exact_uy(x, y, P, L, E, nu, Inertia) result(uy)
    real(dp), intent(in) :: x, y, P, L, E, nu, Inertia
    real(dp) :: uy
    ! Standard analytical v-displacement for cantilever fixed at x=L
    uy = P / (6.0_dp * E * Inertia) * &
         (3.0_dp * nu * x * y**2 + x**3 - 3.0_dp * L**2 * x + 2.0_dp * L**3)
  end function exact_uy

  elemental function shear_traction(y, P, D, Inertia) result(tau_xy)
    real(dp), intent(in) :: y, P, D, Inertia
    real(dp) :: tau_xy
    tau_xy = -P / (2.0_dp * Inertia) * ((D/2.0_dp)**2 - y**2)
  end function shear_traction

  ! ---------------------------------------------------------------
  ! Q1 Element Local Stiffness Matrix Integration (Plane Stress)
  ! ---------------------------------------------------------------
  subroutine get_element_stiffness(dx, dy, E, nu, Ke)
    real(dp), intent(in) :: dx, dy, E, nu
    real(dp), intent(out) :: Ke(8,8)

    real(dp) :: Dmat(3,3), B(3,8)
    real(dp) :: xi, eta, weight, detJ
    integer :: q_i, q_j
    real(dp), parameter :: inv_sqrt3 = 1.0_dp / sqrt(3.0_dp)
    real(dp) :: q_pts(2) = [-inv_sqrt3, inv_sqrt3]
    real(dp) :: dN_dxi(4), dN_deta(4)
    integer :: i

    ! Plane stress constitutive matrix
    Dmat = 0.0_dp
    Dmat(1,1) = 1.0_dp;       Dmat(1,2) = nu
    Dmat(2,1) = nu;           Dmat(2,2) = 1.0_dp
    Dmat(3,3) = (1.0_dp - nu) / 2.0_dp
    Dmat = Dmat * (E / (1.0_dp - nu**2))

    Ke = 0.0_dp
    detJ = (dx / 2.0_dp) * (dy / 2.0_dp)

    ! 2x2 Gauss Quadrature
    do q_i = 1, 2
      do q_j = 1, 2
        xi = q_pts(q_i); eta = q_pts(q_j)

        ! Shape function derivatives wrt local coords
        dN_dxi = 0.25_dp * [-(1.0_dp - eta), (1.0_dp - eta), (1.0_dp + eta), -(1.0_dp + eta)]
        dN_deta = 0.25_dp * [-(1.0_dp - xi), -(1.0_dp + xi), (1.0_dp + xi), (1.0_dp - xi)]

        ! Assemble Strain-Displacement Matrix B
        B = 0.0_dp
        do i = 1, 4
          B(1, 2*i-1) = dN_dxi(i) / (dx / 2.0_dp)  ! dN/dx
          B(2, 2*i)   = dN_deta(i) / (dy / 2.0_dp) ! dN/dy
          B(3, 2*i-1) = B(2, 2*i)
          B(3, 2*i)   = B(1, 2*i-1)
        end do

        ! Ke += B^T * D * B * detJ * weight
        weight = 1.0_dp ! 1.0 * 1.0 for Gauss
        Ke = Ke + matmul(transpose(B), matmul(Dmat, B)) * detJ * weight
      end do
    end do
  end subroutine get_element_stiffness

  ! ---------------------------------------------------------------
  ! Global Assembly (Masked Dirichlet & Optimized Boundary Loops)
  ! ---------------------------------------------------------------
  subroutine assemble_system(NX, NY, L, D, P, E, nu, Inertia, a, rnr, snr, b, nz)
    use y12m_example_helpers, only: compress_matrix
    integer, intent(in) :: NX, NY
    real(dp), intent(in) :: L, D, P, E, nu, Inertia
    real(dp), intent(out) :: a(:), b(:)
    integer, intent(out) :: rnr(:), snr(:), nz

    integer :: ex, ey, i, j, row, col, node
    integer :: el_nodes(4), dof_map(8)
    real(dp) :: dx, dy, Ke(8,8), y
    real(dp) :: local_f
    logical, allocatable :: is_dirichlet(:)

    dx = L / real(NX - 1, dp)
    dy = D / real(NY - 1, dp)

    allocate(is_dirichlet(2 * NX * NY))
    is_dirichlet = .false.
    nz = 0
    b = 0.0_dp

    ! 1. Flag Dirichlet DOFs (Right Boundary x = L)
    ! i = NX, so node = (NX - 1) * NY + j
    do j = 1, NY
      node = (NX - 1) * NY + j
      is_dirichlet(2*node - 1) = .true. ! u_x is fixed
      is_dirichlet(2*node)     = .true. ! u_y is fixed
    end do

    ! 2. Assemble Stiffness Matrix (Loop over elements)
    call get_element_stiffness(dx, dy, E, nu, Ke)

    do ey = 1, NY - 1
      do ex = 1, NX - 1
        ! Counter-clockwise element node mapping (Y-axis first numbering)
        ! ex maps to X-index (i), ey maps to Y-index (j)
        el_nodes(1) = (ex - 1) * NY + ey
        el_nodes(2) = ex * NY + ey
        el_nodes(3) = ex * NY + (ey + 1)
        el_nodes(4) = (ex - 1) * NY + (ey + 1)

        ! Map nodes to global DOFs (x=odd, y=even)
        do i = 1, 4
          dof_map(2*i - 1) = 2 * el_nodes(i) - 1
          dof_map(2*i)     = 2 * el_nodes(i)
        end do

        ! Insert 8x8 block into sparse matrix
        do i = 1, 8
          row = dof_map(i)

          ! MASKING: Skip if Dirichlet boundary
          if (is_dirichlet(row)) cycle

          do j = 1, 8
            col = dof_map(j)
            if (abs(Ke(i,j)) > 0) call add(row, col, Ke(i,j))
          end do
        end do
      end do
    end do

    ! 3. RIGHT BOUNDARY (x = L) : Dirichlet Exact Displacement
    do j = 1, NY
       y = -D/2.0_dp + real(j - 1, dp) * dy
       node = (NX - 1) * NY + j

       ! Equation: E * u_x = E * exact_ux
       call add(2*node - 1, 2*node - 1, E)
       b(2*node - 1) = E * exact_ux(L, y, P, L, D, E, nu, Inertia)

       ! Equation: E * u_y = E * exact_uy (D removed)
       call add(2*node, 2*node, E)
       b(2*node) = E * exact_uy(L, y, P, L, E, nu, Inertia)
    end do

    ! 4. LEFT BOUNDARY (x = 0) : Neumann Shear Traction
    do j = 1, NY
       y = -D/2.0_dp + real(j - 1, dp) * dy
       ! i = 1, so node = (1 - 1) * NY + j = j
       node = j

       ! Add nodal force integration
       if (j == 1 .or. j == NY) then
           local_f = shear_traction(y, P, D, Inertia) * (dy / 2.0_dp)
       else
           local_f = shear_traction(y, P, D, Inertia) * dy
       end if

       b(2*node) = b(2*node) + local_f
    end do

    ! Compress duplicates and update the final 'nz' count
    call compress_matrix(nz, rnr, snr, a)

  contains
    subroutine add(r, c, val)
      integer, intent(in) :: r, c
      real(dp), intent(in) :: val
      nz = nz + 1
      rnr(nz) = r
      snr(nz) = c
      a(nz) = val
    end subroutine add
  end subroutine assemble_system

  ! ---------------------------------------------------------------
  ! Solve
  ! ---------------------------------------------------------------
  subroutine run(NX, NY, outfile)
    use y12m, only: y12ma
    integer, intent(in) :: NX, NY
    character(len=*), intent(in) :: outfile

    ! Physics constants (Aluminum 7075)
    real(dp) :: L = 30.0_dp, D = 5.0_dp, P = 1000.0_dp
    real(dp) :: E = 72.1e9_dp, nu = 0.33_dp, Inertia

    integer :: ndof, nz_max, nz, nn, ifail
    real(dp), allocatable :: a(:), b(:), pivot(:)
    integer, allocatable :: snr(:), rnr(:), ha(:,:)
    real(dp) :: aflag(8), err_max_x, err_max_y
    real(dp) :: x_val, y_val, ux_ex, uy_ex
    integer :: iflag(10), i, j, node, funit

    Inertia = (D**3) / 12.0_dp
    ndof = 2 * NX * NY

    ! Extremely generous workspace for Block-Sparse fill-in
    nz_max = 36 * ndof
    nn     = 200 * ndof
    allocate(a(nn), pivot(ndof), b(ndof), snr(nn), rnr(nn), ha(ndof, 11))

    aflag = 0.0_dp
    iflag = 0
    ifail = 0

    write(output_unit, '(a,i0,a,i0,a,i0,a)') 'Grid: ', NX, ' x ', NY, ' (', ndof, ' DOFs)'

    call assemble_system(NX, NY, L, D, P, E, nu, Inertia, a, rnr, snr, b, nz)
    call y12ma(ndof, nz, a, snr, nn, rnr, nn, pivot, ha, ndof, aflag, iflag, b, ifail)

    if (ifail /= 0) then
      write(error_unit, '(a,i0)') 'ERROR: y12ma returned IFAIL = ', ifail
      stop 1
    end if

    ! Compute Errors & Write Output
    err_max_x = 0.0_dp; err_max_y = 0.0_dp

    open(newunit=funit, file=outfile, status='unknown', action='write')
    write(funit, '(a)') '# Columns: x   y   ux_num   uy_num   ux_exact   uy_exact'

    do j = 1, NY
      y_val = -D/2.0_dp + real(j-1,dp)*D/(NY-1)
      do i = 1, NX
        x_val = real(i-1,dp)*L/(NX-1)
        ! Y-axis first numbering for error check
        node = (i - 1) * NY + j

        ux_ex = exact_ux(x_val, y_val, P, L, D, E, nu, Inertia)
        uy_ex = exact_uy(x_val, y_val, P, L, E, nu, Inertia)

        err_max_x = max(err_max_x, abs(b(2*node - 1) - ux_ex))
        err_max_y = max(err_max_y, abs(b(2*node)     - uy_ex))

        write(funit, '(6(1x,es14.6))') x_val, y_val, b(2*node - 1), b(2*node), ux_ex, uy_ex
      end do
      write(funit, *) ! Empty line for gnuplot grid block formatting
    end do

    close(funit)

    write(output_unit, '(a,es12.4)') 'Max Error u_x : ', err_max_x
    write(output_unit, '(a,es12.4)') 'Max Error u_y : ', err_max_y
    write(output_unit, '(a)') 'Solution written to ' // trim(outfile)

  end subroutine run
end module elasticity_solver

! ==============================================================
! Main program: parse CLI
! ==============================================================
program elasticity_2d
  use elasticity_solver, only: run
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none

  integer :: NX, NY, ios, iarg, nargs, pos_count
  character(len=256) :: arg, outfile

  ! Defaults for a 30m x 5m beam
  ! (61x11 nodes yields 60x10 perfectly square elements)
  NX      = 61
  NY      = 11
  outfile = 'elasticity.dat'

  nargs = command_argument_count()
  iarg = 1
  pos_count = 0

  ! Parse arguments
  do while (iarg <= nargs)
    call get_command_argument(iarg, arg)

    if (trim(arg) == '--help' .or. trim(arg) == '-h') then
      write(*, '(a)') 'Usage: elasticity_2d [--help] [NX] [NY] [output_file]'
      write(*, '(a)') '  NX : Number of nodes in the x-direction (Length)'
      write(*, '(a)') '  NY : Number of nodes in the y-direction (Depth)'
      stop 0
    else
      ! Positional arguments
      pos_count = pos_count + 1
      if (pos_count == 1) read(arg, *, iostat=ios) NX
      if (pos_count == 2) read(arg, *, iostat=ios) NY
      if (pos_count == 3) outfile = trim(arg)
    end if

    iarg = iarg + 1
  end do

  if (NX < 2 .or. NY < 2) then
    write(error_unit, '(a)') 'ERROR: Grid dimensions must be at least 2x2'
    stop 1
  end if

  call run(NX, NY, trim(outfile))

end program elasticity_2d
