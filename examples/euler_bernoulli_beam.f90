! SPDX-License-Identifier: GPL-2.0-only
!
! euler_bernoulli_beam.f90 — 1D Hermite cubic FEM for the Euler-Bernoulli beam
!
! Solves the clamped-free (cantilever) Euler-Bernoulli beam with linearly
! tapered bending stiffness on the unit interval [0, 1]:
!
!   d^2/dx^2 [ EI(x) d^2w/dx^2 ] = q(x)   on (0, 1),
!   w(0) = 0,  w'(0) = 0                     (clamped at x = 0),
!   EI(1) w''(1) = 0,  d/dx[EI w'']|_1 = 0   (free tip at x = 1).
!
! The bending stiffness varies linearly from root to tip:
!   EI(x) = 1 + x   (normalised units)
!
! Method of Manufactured Solution (MMS):
!   The exact solution is chosen as the degree-5 polynomial
!     w(x)     = 20 x^2 - 10 x^3 + x^5
!     theta(x) = w'(x) = 40 x - 30 x^2 + 5 x^4
!   which satisfies all four boundary conditions identically.
!   The corresponding distributed load is
!     q(x) = (EI w'')'' = -120 + 120 x + 240 x^2.
!
! Discretisation:
!   Hermite cubic shape functions, 2 DOFs per node (deflection w and
!   rotation theta).  The 4x4 element stiffness matrix is computed by a
!   2-point Gauss-Legendre rule (exact for the linear-EI integrand) and the
!   load vector by a 3-point rule (exact for the quadratic-in-x load).
!   For N elements the global stiffness matrix is (2 N) x (2 N) (after
!   applying the clamped BC at x=0) and has half-bandwidth 3.  Adjacent DOFs
!   carry different physical meanings (deflection/rotation coupling), making
!   the sparsity pattern distinct from purely nodal discretisations.
!
! Node/DOF numbering (N=4 elements example):
!
!   Physical node:  1[C]        2           3           4           5
!   Global DOF  :   1    2      3    4      5    6      7    8      9   10
!                   w    t      w    t      w    t      w    t      w    t
!   Clamped [C] :  [x    x] <- removed; reduced index  r = global_DOF - 2
!   Reduced DOF :                r1   r2     r3   r4     r5   r6     r7   r8
!
! Block-tridiagonal sparsity of K after clamping w=t=0 at node 1 (N=4):
!
!              r1   r2 | r3   r4 | r5   r6 | r7   r8
!   r1 (w2) [  x    x |  x    x |  .    . |  .    . ]
!   r2 (t2) [  x    x |  x    x |  .    . |  .    . ]
!           [---------+---------+---------+---------]
!   r3 (w3) [  x    x |  x    x |  x    x |  .    . ]
!   r4 (t3) [  x    x |  x    x |  x    x |  .    . ]
!           [---------+---------+---------+---------]
!   r5 (w4) [  .    . |  x    x |  x    x |  x    x ]
!   r6 (t4) [  .    . |  x    x |  x    x |  x    x ]
!           [---------+---------+---------+---------]
!   r7 (w5) [  .    . |  .    . |  x    x |  x    x ]
!   r8 (t5) [  .    . |  .    . |  x    x |  x    x ]
!
!   Each 2x2 sub-block couples deflection (w) and rotation (t=theta) at a
!   node or between adjacent nodes.  Block half-bandwidth = 1.
!
! Assembly:
!   Element contributions are appended directly to COO triplet arrays
!   (Phase 1).  Because adjacent elements share boundary-node DOFs the raw
!   list contains duplicate (row, col) pairs; a counting-sort pass compresses
!   and sums these before the solver call (Phase 2).
!
! Solver: y12ma (double-precision sparse direct).
!
! Expected nodal convergence on a uniform mesh (superconvergence for 1D
! Hermite elements with smooth variable coefficients):
!   max |w_h   - w|   = O(h^4)
!   max |theta_h - w'| = O(h^4)
!
! Note: the beam stiffness matrix has condition number O(h^{-4}), so for
! very fine meshes the rounding errors from direct LU factorisation begin
! to dominate before the FEM discretisation error.  Clean O(h^4) convergence
! is observed for the mesh sizes used in the convergence study (N <= 32).
!
! Usage
! -----
!   euler_bernoulli_beam [--help] [--conv] [N] [output_file]
!
!     N            number of Hermite cubic elements (default 16)
!     output_file  output data file (default: euler_bernoulli_beam.dat)
!     --conv       grid-refinement convergence study (N=2,4,8,16,32)

module beam_solver
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
  implicit none
  private
  public :: run

  integer, parameter :: dp = kind(1.0d0)

  ! 2-point Gauss-Legendre quadrature on [0,1]
  integer, parameter :: ng2 = 2
  real(dp), parameter :: sq3inv = 1.0_dp / sqrt(3.0_dp)
  real(dp), parameter :: gp2(2) = [0.5_dp - 0.5_dp*sq3inv, &
                                    0.5_dp + 0.5_dp*sq3inv]
  real(dp), parameter :: gw2(2) = [0.5_dp, 0.5_dp]

  ! 3-point Gauss-Legendre quadrature on [0,1]
  integer, parameter :: ng3 = 3
  real(dp), parameter :: sq35 = sqrt(3.0_dp / 5.0_dp)
  real(dp), parameter :: gp3(3) = [0.5_dp - 0.5_dp*sq35, &
                                    0.5_dp,                &
                                    0.5_dp + 0.5_dp*sq35]
  real(dp), parameter :: gw3(3) = [5.0_dp/18.0_dp, &
                                    8.0_dp/18.0_dp, &
                                    5.0_dp/18.0_dp]

contains

  !> Bending stiffness: EI(x) = 1 + x  (linear taper, normalised units).
  elemental function ei_func(x) result(EI)
    real(dp), intent(in) :: x
    real(dp) :: EI
    EI = 1.0_dp + x
  end function ei_func

  !> Exact deflection (MMS): w(x) = 20 x^2 - 10 x^3 + x^5.
  elemental function exact_w(x) result(w)
    real(dp), intent(in) :: x
    real(dp) :: w
    w = x**2 * (20.0_dp - 10.0_dp*x + x**3)
  end function exact_w

  !> Exact rotation (MMS): theta(x) = w'(x) = 40 x - 30 x^2 + 5 x^4.
  elemental function exact_theta(x) result(theta)
    real(dp), intent(in) :: x
    real(dp) :: theta
    theta = x * (40.0_dp - 30.0_dp*x + 5.0_dp*x**3)
  end function exact_theta

  !> Distributed load derived from MMS: q(x) = (EI w'')'' = -120 + 120 x + 240 x^2.
  elemental function load_q(x) result(q)
    real(dp), intent(in) :: x
    real(dp) :: q
    q = -120.0_dp + 120.0_dp*x + 240.0_dp*x**2
  end function load_q

  !> Evaluate the four Hermite shape functions at local coordinate xi in [0,1].
  !> N(1) = H1 = 1-3xi^2+2xi^3,  N(2) = h*xi*(1-xi)^2,
  !> N(3) = H3 = 3xi^2-2xi^3,    N(4) = h*xi^2*(xi-1).
  pure subroutine shape_N(xi, h, N)
    real(dp), intent(in) :: xi, h
    real(dp), intent(out) :: N(4)
    N(1) = 1.0_dp - 3.0_dp*xi**2 + 2.0_dp*xi**3
    N(2) = h * xi * (1.0_dp - xi)**2
    N(3) = 3.0_dp*xi**2 - 2.0_dp*xi**3
    N(4) = h * xi**2 * (xi - 1.0_dp)
  end subroutine shape_N

  !> Second derivatives d^2N_i/dxi^2 scaled by h^2 so that
  !> d^2N_i/dx^2 = D2N(i) / h^2.
  pure subroutine shape_D2N(xi, h, D2N)
    real(dp), intent(in) :: xi, h
    real(dp), intent(out) :: D2N(4)
    D2N(1) = -6.0_dp + 12.0_dp*xi
    D2N(2) = h * (-4.0_dp + 6.0_dp*xi)
    D2N(3) =  6.0_dp - 12.0_dp*xi
    D2N(4) = h * (-2.0_dp + 6.0_dp*xi)
  end subroutine shape_D2N

  !> 4x4 element stiffness matrix for variable EI, computed by 2-point
  !> Gauss-Legendre quadrature (exact for EI linear in x times quadratic
  !> second-derivative products, i.e., for any polynomial EI of degree 1).
  !>
  !> Ke(i,j) = integral_0^h  EI(x) * N_i''(x) * N_j''(x) dx.
  subroutine elem_stiffness(h, x_e, Ke)
    real(dp), intent(in) :: h, x_e
    real(dp), intent(out) :: Ke(4,4)
    real(dp) :: xi, EI_loc, D2N(4)
    integer :: g, i, j
    Ke = 0.0_dp
    do g = 1, ng2
      xi     = gp2(g)
      EI_loc = ei_func(x_e + xi*h)
      call shape_D2N(xi, h, D2N)
      ! Ke(i,j) += (gw * EI / h^3) * D2N_i * D2N_j
      do j = 1, 4
        do i = 1, 4
          Ke(i,j) = Ke(i,j) + gw2(g) * (EI_loc / h**3) * D2N(i) * D2N(j)
        end do
      end do
    end do
  end subroutine elem_stiffness

  !> Consistent equivalent nodal load vector for q(x) = -120+120x+240x^2,
  !> computed by 3-point Gauss-Legendre quadrature (exact for degree 5).
  !>
  !> fe(i) = integral_0^h  N_i(x) * q(x) dx.
  subroutine elem_force(h, x_e, fe)
    real(dp), intent(in) :: h, x_e
    real(dp), intent(out) :: fe(4)
    real(dp) :: xi, q_loc, N(4)
    integer :: g, i
    fe = 0.0_dp
    do g = 1, ng3
      xi    = gp3(g)
      q_loc = load_q(x_e + xi*h)
      call shape_N(xi, h, N)
      do i = 1, 4
        fe(i) = fe(i) + gw3(g) * h * N(i) * q_loc
      end do
    end do
  end subroutine elem_force

  !> Solve the tapered cantilever beam with N Hermite cubic elements.
  !>
  !> On return, err_w and err_theta are the maximum absolute nodal errors
  !> in deflection and rotation, respectively.
  !> If outfile is non-empty the nodal solution is written there.
  subroutine run(N, outfile, err_w, err_theta)
    use y12m, only: y12ma
    integer, intent(in) :: N
    character(len=*), intent(in) :: outfile
    real(dp), intent(out) :: err_w, err_theta

    integer :: ndof, nz_raw_max, nz_raw, nz, nn, ifail
    integer :: e, li, lj, gi, gj, ri, rj, nd, funit
    real(dp) :: h, x, x_e, w_ex, theta_ex
    real(dp) :: Ke(4,4), fe_vec(4)
    real(dp), allocatable :: val_raw(:), a(:), b(:), pivot(:)
    integer, allocatable :: row_raw(:), col_raw(:), snr(:), rnr(:), ha(:,:)
    real(dp) :: aflag(8)
    integer :: iflag(10)

    h    = 1.0_dp / real(N, dp)
    ndof = 2 * N   ! free DOFs: clamped BC at x=0 removes 2 DOFs

    ! Upper bound on raw COO entries before compression: each element
    ! contributes at most 4*4=16 triplets (fewer for element 1, where 2 DOFs
    ! are clamped).  Adjacent elements share a boundary node, so the
    ! compressed entry count is strictly less than the raw count.
    nz_raw_max = 16 * N

    ! y12ma needs extra space for LU factors; compressed non-zeros for the
    ! block-tridiagonal pattern are at most 12*N.  A 10x buffer is consistent
    ! with other examples (e.g. darcy_flow.f90).
    nn = max(10 * 12 * N, 100 * ndof)
    allocate(a(nn), pivot(ndof), b(ndof), snr(nn), rnr(nn), ha(ndof, 11))
    allocate(val_raw(nz_raw_max), row_raw(nz_raw_max), col_raw(nz_raw_max))

    b     = 0.0_dp
    aflag = 0.0_dp
    iflag = 0
    ifail = 0

    ! Phase 1 — sparse element assembly.
    ! Global DOF numbering: node k owns DOFs 2k-1 (w) and 2k (theta).
    ! Node 1 is clamped: DOFs 1 and 2 are eliminated.
    ! Reduced (free) index: r = global_DOF - 2,  valid for global_DOF > 2.
    ! Adjacent elements share the DOFs of their common boundary node, so the
    ! raw COO list may contain duplicate (ri,rj) pairs at shared nodes.
    nz_raw = 0
    do e = 1, N
      x_e = real(e - 1, dp) * h
      call elem_stiffness(h, x_e, Ke)
      call elem_force(h, x_e, fe_vec)
      do li = 1, 4
        gi = 2*e - 2 + li   ! global DOF index for local DOF li
        if (gi <= 2) cycle   ! skip clamped DOFs at node 1
        ri = gi - 2          ! reduced (free) row index
        b(ri) = b(ri) + fe_vec(li)
        do lj = 1, 4
          gj = 2*e - 2 + lj
          if (gj <= 2) cycle
          rj = gj - 2
          nz_raw = nz_raw + 1
          row_raw(nz_raw) = ri
          col_raw(nz_raw) = rj
          val_raw(nz_raw) = Ke(li, lj)
        end do
      end do
    end do

    ! Phase 2 — compress: sort by (row,col) and sum duplicate entries.
    ! y12ma returns IFAIL=11 if any two triplets share the same (i,j) position.
    call coo_compress(ndof, nz_raw, row_raw, col_raw, val_raw, nz, rnr, snr, a)
    deallocate(val_raw, row_raw, col_raw)

    call y12ma(ndof, nz, a, snr, nn, rnr, nn, pivot, ha, ndof, aflag, iflag, b, ifail)

    if (ifail /= 0) then
      write(error_unit, '(a,i0)') 'ERROR: y12ma returned IFAIL = ', ifail
      stop 1
    end if

    ! Maximum nodal errors over free DOFs (nodes 2 to N+1).
    ! Clamped node 1 removed DOFs 1 and 2; reduced index = global_DOF - 2.
    !   deflection at node nd : reduced index 2*(nd-1) - 1  (= 2*nd - 3)
    !   rotation   at node nd : reduced index 2*(nd-1)      (= 2*nd - 2)
    err_w     = 0.0_dp
    err_theta = 0.0_dp
    do nd = 2, N + 1
      x        = real(nd - 1, dp) * h
      w_ex     = exact_w(x)
      theta_ex = exact_theta(x)
      err_w     = max(err_w,     abs(b(2*(nd-1)-1) - w_ex))   ! deflection
      err_theta = max(err_theta, abs(b(2*(nd-1))   - theta_ex)) ! rotation
    end do

    if (len_trim(outfile) > 0) then
      open(newunit=funit, file=outfile, status='unknown', action='write')
      write(funit, '(a)') '# Euler-Bernoulli beam: tapered, Hermite cubic FEM'
      write(funit, '(a)') '# EI(x)=1+x, q(x)=-120+120x+240x^2; exact: w=20x^2-10x^3+x^5'
      write(funit, '(a)') '# Columns: x   w_num   theta_num   w_exact   theta_exact'
      write(funit, '(5(1x,es14.6))') 0.0_dp, 0.0_dp, 0.0_dp, &
           exact_w(0.0_dp), exact_theta(0.0_dp)
      do nd = 2, N + 1
        x = real(nd - 1, dp) * h
        write(funit, '(5(1x,es14.6))') x, &
             b(2*(nd-1)-1), b(2*(nd-1)), &   ! deflection, rotation
             exact_w(x), exact_theta(x)
      end do
      close(funit)
      write(output_unit, '(a,a)') 'Solution written to ', trim(outfile)
    end if

  end subroutine run

  !> Sort and sum duplicate (row,col) entries in a COO triplet list.
  !>
  !> Uses a counting-sort pass (keyed on (row-1)*ndof + col) followed by a
  !> single linear sweep to accumulate values for repeated positions.
  !> On entry nz_in triplets (row_in, col_in, val_in) may contain duplicates;
  !> on exit (row_out, col_out, val_out) is a duplicate-free, row-major sorted
  !> list of length nz_out.  The output arrays must each have at least nz_in
  !> elements (worst case: no duplicates, nz_out == nz_in).
  !>
  !> Note: y12ma returns IFAIL=11 when two entries share the same (i,j)
  !> position, so compression is mandatory for element-by-element FEM assembly.
  subroutine coo_compress(ndof, nz_in, row_in, col_in, val_in, &
                          nz_out, row_out, col_out, val_out)
    integer, intent(in) :: ndof, nz_in
    integer, intent(in) :: row_in(nz_in), col_in(nz_in)
    real(dp), intent(in) :: val_in(nz_in)
    integer, intent(out) :: nz_out
    integer, intent(out) :: row_out(:), col_out(:)
    real(dp), intent(out) :: val_out(:)

    integer :: k, key, last_key, nkeys
    integer, allocatable :: cnt(:), sorted_idx(:)

    ! Counting sort: key = (row-1)*ndof + col  (1-based, row-major order)
    nkeys = ndof * ndof
    allocate(cnt(nkeys), sorted_idx(nz_in))
    cnt = 0
    do k = 1, nz_in
      key = (row_in(k) - 1) * ndof + col_in(k)
      cnt(key) = cnt(key) + 1
    end do
    do key = 2, nkeys
      cnt(key) = cnt(key) + cnt(key - 1)
    end do
    do k = nz_in, 1, -1
      key = (row_in(k) - 1) * ndof + col_in(k)
      sorted_idx(cnt(key)) = k
      cnt(key) = cnt(key) - 1
    end do

    ! Linear sweep: compress duplicate (row,col) pairs by summing values
    nz_out = 0
    last_key = -1
    do k = 1, nz_in
      key = (row_in(sorted_idx(k)) - 1) * ndof + col_in(sorted_idx(k))
      if (key /= last_key) then
        nz_out = nz_out + 1
        row_out(nz_out) = row_in(sorted_idx(k))
        col_out(nz_out) = col_in(sorted_idx(k))
        val_out(nz_out) = val_in(sorted_idx(k))
        last_key = key
      else
        val_out(nz_out) = val_out(nz_out) + val_in(sorted_idx(k))
      end if
    end do

    deallocate(cnt, sorted_idx)
  end subroutine coo_compress

end module beam_solver

! ==============================================================
! Main program: run a single solve or a convergence study
! ==============================================================
program euler_bernoulli_beam
  use beam_solver, only: run
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  integer :: N, ios, iarg, nargs, pos_count, i
  character(len=256) :: arg, outfile
  real(dp) :: err_w, err_theta, h, h_prev, err_w_prev, err_theta_prev
  real(dp) :: rate_w, rate_theta
  integer, parameter :: n_refine = 5
  integer :: ns(n_refine)
  logical :: do_conv

  N       = 16
  outfile = 'euler_bernoulli_beam.dat'
  do_conv = .false.
  ns      = [2, 4, 8, 16, 32]

  nargs     = command_argument_count()
  iarg      = 1
  pos_count = 0

  do while (iarg <= nargs)
    call get_command_argument(iarg, arg)

    if (trim(arg) == '--help' .or. trim(arg) == '-h') then
      write(*, '(a)') 'Usage: euler_bernoulli_beam [--help] [--conv] [N] [output_file]'
      write(*, '(a)') '  N            number of Hermite cubic elements (default 16)'
      write(*, '(a)') '  output_file  output data file (default: euler_bernoulli_beam.dat)'
      write(*, '(a)') '  --conv       grid-refinement convergence study (N=2,4,8,16,32)'
      stop 0
    else if (trim(arg) == '--conv') then
      do_conv = .true.
    else
      pos_count = pos_count + 1
      if (pos_count == 1) then
        read(arg, *, iostat=ios) N
        if (ios /= 0 .or. N < 1) then
          write(error_unit, '(a)') 'ERROR: N must be a positive integer'
          stop 1
        end if
      else if (pos_count == 2) then
        outfile = trim(arg)
      end if
    end if

    iarg = iarg + 1
  end do

  if (do_conv) then
    write(output_unit, '(a)') &
         'Euler-Bernoulli beam — Hermite cubic FEM convergence study'
    write(output_unit, '(a)') &
         'Tapered cantilever: EI(x)=1+x, q(x)=-120+120x+240x^2'
    write(output_unit, '(a)') &
         'Exact: w=20x^2-10x^3+x^5,  theta=40x-30x^2+5x^4'
    write(output_unit, *)
    write(output_unit, '(a)') &
         '   N        h        max|e_w|      max|e_t|     rate_w  rate_t'
    write(output_unit, '(a)') &
         '  ---  ----------  -----------  -----------  -------  -------'

    h_prev         = 0.0_dp
    err_w_prev     = 0.0_dp
    err_theta_prev = 0.0_dp

    do i = 1, n_refine
      call run(ns(i), '', err_w, err_theta)
      h = 1.0_dp / real(ns(i), dp)
      if (i == 1) then
        write(output_unit, '(i5,f12.6,2(2x,es11.4),2(9x,a1))') &
             ns(i), h, err_w, err_theta, '-', '-'
      else
        rate_w     = log(err_w_prev     / err_w)     / log(h_prev / h)
        rate_theta = log(err_theta_prev / err_theta) / log(h_prev / h)
        write(output_unit, '(i5,f12.6,2(2x,es11.4),2f9.2)') &
             ns(i), h, err_w, err_theta, rate_w, rate_theta
      end if
      h_prev         = h
      err_w_prev     = err_w
      err_theta_prev = err_theta
    end do

    ! Write the finest-mesh nodal solution to file
    call run(ns(n_refine), trim(outfile), err_w, err_theta)
  else
    call run(N, trim(outfile), err_w, err_theta)
    write(output_unit, '(a,i0,a)') 'N = ', N, ' Hermite cubic elements'
    write(output_unit, '(a,es12.4)') 'Max nodal error in deflection : ', err_w
    write(output_unit, '(a,es12.4)') 'Max nodal error in rotation   : ', err_theta
  end if

end program euler_bernoulli_beam
