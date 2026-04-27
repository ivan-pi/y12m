! SPDX-License-Identifier: GPL-2.0-only
!
! y12m_mm.f90 -- Matrix Market test driver for the y12m sparse solver.
!
! Reads a sparse matrix from a Matrix Market file (.mtx), constructs the
! right-hand side b = A * [1, 1, ..., 1]^T (row sums), solves Ax = b with
! the y12ma high-level interface (double precision), then verifies the
! solution against the known exact answer x = [1, 1, ..., 1]^T.
!
! Usage:
!   y12m_mm [options] <input.mtx>
!   y12m_mm [options] -              (read from standard input)
!
! Options:
!   -o, --output FILE    Write results to FILE instead of standard output.
!   -v, --verbose        Also print solver diagnostics (growth factor, min
!                        pivot, garbage-collection counts, etc.).
!   -h, --help           Print this help and exit.
!
! Supported Matrix Market formats:
!   - representation: coordinate
!   - field:          real, integer
!   - symmetry:       general, symmetric, skew-symmetric
!
! For symmetric/skew-symmetric matrices the stored triangle is expanded to
! the full matrix before solving, so the row-sum RHS is correct.
!
! Note: y12ma always uses fixed internal defaults for AFLAG(1-4) and
! IFLAG(2-5).  To tune solver parameters use the lower-level
! y12mb/y12mc/y12md interface instead.
!

! =============================================================================
! Module y12m_mm_driver
!
! Contains all computational logic: reading the Matrix Market file, expanding
! symmetry, building the RHS, calling the solver, and reporting results.
! The main program is a thin wrapper that handles only command-line parsing.
! =============================================================================
module y12m_mm_driver
  use y12m, only: y12ma
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none
  private

  public :: run
  public :: cli_opts

  integer, parameter :: dp = kind(1.0d0)

  ! NIST MMIO interface block (mminfo, mmread, mmwrite).
  ! Placed in the module specification part so the interfaces are visible to
  ! all procedures in the module.  A bare interface block is only valid in a
  ! specification part, not inside a procedure body.
  include "mmio.fi"

  !> Options collected from the command line.
  type :: cli_opts
    character(len=512) :: infile      = ''
    character(len=512) :: outfile     = ''
    logical            :: to_file     = .false.
    logical            :: verbose     = .false.
    logical            :: print_all   = .false.  ! print full x and r vectors
    logical            :: use_stdin   = .false.
    integer            :: fill_factor = 3        ! nn = fill_factor * z
  end type cli_opts

contains

  ! ---------------------------------------------------------------------------
  !> Sparse COO matrix-vector product: y := y + alpha * A * x
  !>
  !> A is stored in coordinate format with nz entries.
  !> No assumptions are made about ordering.
  !>
  !> @param[in]     m      Number of rows (length of y).
  !> @param[in]     nz     Number of non-zero entries.
  !> @param[in]     alpha  Scalar multiplier.
  !> @param[in]     val    Non-zero values,   val(1:nz).
  !> @param[in]     row    Row indices,        row(1:nz), 1-based.
  !> @param[in]     col    Column indices,     col(1:nz), 1-based.
  !> @param[in]     x      Input vector,       x(1:*).
  !> @param[inout]  y      Accumulation vector, y(1:m).
  subroutine sparse_gemv(m, nz, alpha, val, row, col, x, y)
    integer,  intent(in)    :: m, nz
    real(dp), intent(in)    :: alpha
    real(dp), intent(in)    :: val(nz)
    integer,  intent(in)    :: row(nz), col(nz)
    real(dp), intent(in)    :: x(*)
    real(dp), intent(inout) :: y(m)
    integer :: k
    do k = 1, nz
      y(row(k)) = y(row(k)) + alpha * val(k) * x(col(k))
    end do
  end subroutine sparse_gemv

  ! ---------------------------------------------------------------------------
  !> Main driver: read, validate, solve, and report.
  !>
  !> @param[in]  opts  Command-line options (see cli_opts).
  subroutine run(opts)
    type(cli_opts), intent(in) :: opts

    ! ---- Matrix Market metadata ----
    character(len=10) :: rep
    character(len=7)  :: field
    character(len=19) :: symm
    integer           :: mm_rows, mm_cols, mm_entries

    ! ---- Temporary read buffers (mmread requires these to be allocatable) ----
    integer,           allocatable :: indx_tmp(:), jndx_tmp(:), ival_tmp(:)
    real(dp),          allocatable :: rval_tmp(:)
    complex,           allocatable :: cval_tmp(:)  ! plain complex, as declared in mmio.f
    integer :: nnzmax   ! buffer capacity for mmread; kept separate from
                        ! mm_entries to avoid aliasing the IN and OUT arguments

    ! ---- Full COO arrays after symmetry expansion ----
    integer,  allocatable :: row_coo(:), col_coo(:)
    real(dp), allocatable :: val_coo(:)
    integer :: z   ! non-zeros after expansion

    ! ---- Saved originals for residual computation ----
    integer,  allocatable :: row_orig(:), col_orig(:)
    real(dp), allocatable :: val_orig(:), b_orig(:)

    ! ---- Solver arrays ----
    real(dp), allocatable :: a(:), pivot(:), b(:), aflag(:)
    integer,  allocatable :: snr(:), rnr(:), ha(:,:), iflag(:)

    ! ---- Post-solve diagnostics ----
    real(dp), allocatable :: resid(:)
    real(dp) :: res_inf, res_2norm, err_inf, err_2norm, elapsed

    ! ---- Miscellaneous ----
    integer :: n, nn, nn1, iha, ifail
    integer :: inp_unit, out_unit, ios
    integer :: i, t1, t2, clock_rate

    ! =========================================================================
    ! 1. Open input
    ! =========================================================================
    if (opts%use_stdin) then
      inp_unit = 5
    else
      open(newunit=inp_unit, file=trim(opts%infile), status='old', &
           action='read', iostat=ios)
      if (ios /= 0) then
        write(error_unit,'(3a)') 'Error: cannot open "', trim(opts%infile), '"'
        stop 1
      end if
    end if

    ! =========================================================================
    ! 2. Read Matrix Market header and print identification immediately
    ! =========================================================================
    call mminfo(inp_unit, rep, field, symm, mm_rows, mm_cols, mm_entries)

    write(*, '(a,a)')       'File:     ', trim(opts%infile)
    write(*, '(a,i0,a,i0)') 'Size:     ', mm_rows, ' x ', mm_cols
    write(*, '(a,i0)')      'Entries:  ', mm_entries
    write(*, '(a,a)')       'Field:    ', trim(field)
    write(*, '(a,a,a,a)')   'Format:   ', trim(rep), ',  symmetry: ', trim(symm)
    write(*, '(a)') ''

    ! =========================================================================
    ! 3. Validate header
    ! =========================================================================
    if (trim(rep) /= 'coordinate') then
      write(error_unit,'(3a)') &
          'Error: only "coordinate" format is supported, got "', trim(rep), '".'
      stop 1
    end if

    if (trim(field) /= 'real' .and. trim(field) /= 'integer') then
      write(error_unit,'(3a)') &
          'Error: only "real" and "integer" fields are supported, got "', &
          trim(field), '".'
      stop 1
    end if

    if (mm_rows /= mm_cols) then
      write(error_unit,'(a,i0,a,i0,a)') &
          'Error: matrix is not square (', mm_rows, ' x ', mm_cols, &
          '); y12ma requires a square system.'
      stop 1
    end if

    n = mm_rows
    if (n < 2) then
      write(error_unit,'(a,i0,a)') &
          'Error: matrix dimension n=', n, &
          ' is less than 2 (y12ma requires n >= 2).'
      stop 1
    end if

    ! =========================================================================
    ! 4. Read non-zero triplets
    ! =========================================================================
    nnzmax = mm_entries
    allocate(indx_tmp(nnzmax), jndx_tmp(nnzmax), ival_tmp(nnzmax), &
             rval_tmp(nnzmax), cval_tmp(nnzmax))
    indx_tmp = 0
    jndx_tmp = 0
    ival_tmp = 0
    rval_tmp = 0.0_dp

    ! mmread rewinds and re-reads the file from the beginning.
    call mmread(inp_unit, rep, field, symm, mm_rows, mm_cols, mm_entries, &
                nnzmax, indx_tmp, jndx_tmp, ival_tmp, rval_tmp, cval_tmp)

    if (.not. opts%use_stdin) close(inp_unit)

    deallocate(ival_tmp, cval_tmp)

    ! =========================================================================
    ! 5. Expand symmetry to full matrix
    ! =========================================================================
    call expand_symmetry(symm, mm_entries, indx_tmp, jndx_tmp, rval_tmp, &
                         row_coo, col_coo, val_coo, z)

    deallocate(indx_tmp, jndx_tmp, rval_tmp)

    write(*, '(a,i0)') 'NNZ (after symmetry expansion): ', z
    write(*, '(a)') ''

    ! =========================================================================
    ! 6. Build b = A * 1  (row sums) — exact solution is x = 1
    ! =========================================================================
    allocate(b(n), b_orig(n), pivot(n))
    b(1:n) = 0.0_dp
    call sparse_gemv(n, z, 1.0_dp, val_coo, row_coo, col_coo, &
                     [(1.0_dp, i = 1, n)], b)
    b_orig = b

    ! =========================================================================
    ! 7. Save original COO for residual computation after y12ma overwrites a
    ! =========================================================================
    allocate(val_orig(z), row_orig(z), col_orig(z))
    val_orig = val_coo
    row_orig = row_coo
    col_orig = col_coo

    ! =========================================================================
    ! 8. Allocate solver arrays and call y12ma
    !    NN >= 2*Z (recommended 2*Z to 3*Z); NN1 >= Z (recommended 2*Z to 3*Z)
    ! =========================================================================
    nn  = opts%fill_factor * z
    nn1 = opts%fill_factor * z
    iha = n

    allocate(a(nn), snr(nn), rnr(nn1), ha(iha, 11), aflag(8), iflag(10))

    a(1:z)   = val_coo
    snr(1:z) = col_coo
    rnr(1:z) = row_coo

    deallocate(row_coo, col_coo, val_coo)

    ifail = 0
    call system_clock(t1, clock_rate)
    call y12ma(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, &
               aflag, iflag, b, ifail)
    call system_clock(t2)

    elapsed = real(t2 - t1, dp) / real(clock_rate, dp)

    ! =========================================================================
    ! 9. Open output stream
    ! =========================================================================
    if (opts%to_file) then
      open(newunit=out_unit, file=trim(opts%outfile), status='replace', &
           action='write', iostat=ios)
      if (ios /= 0) then
        write(error_unit,'(3a)') &
            'Error: cannot open output file "', trim(opts%outfile), '"'
        stop 1
      end if
    else
      out_unit = 6
    end if

    ! =========================================================================
    ! 10. Report solver failure if any
    ! =========================================================================
    if (ifail /= 0) then
      write(out_unit,'(a,i0)') 'IFAIL = ', ifail
      write(out_unit,'(a)')    'Solver failed; no solution available.'
      if (opts%to_file) close(out_unit)
      stop 1
    end if

    ! =========================================================================
    ! 11. Compute residual r = b_orig - A * x  via sparse_gemv
    !     b(1:n) now holds the computed solution x.
    ! =========================================================================
    allocate(resid(n))
    resid = b_orig
    call sparse_gemv(n, z, -1.0_dp, val_orig, row_orig, col_orig, b, resid)

    res_inf   = maxval(abs(resid))
    res_2norm = norm2(resid)
    err_inf   = maxval(abs(b(1:n) - 1.0_dp))
    err_2norm = norm2(b(1:n) - 1.0_dp)

    ! =========================================================================
    ! 12. Print results
    ! =========================================================================
    call print_results(out_unit, opts%verbose, opts%print_all, n, b, resid, &
                       res_inf, res_2norm, err_inf, err_2norm, &
                       elapsed, aflag, iflag, ifail)

    if (opts%to_file) close(out_unit)

  end subroutine run

  ! ---------------------------------------------------------------------------
  !> Expand a Matrix Market COO triplet set to the full (unsymmetric) matrix.
  !>
  !> For "general" symmetry the arrays are copied unchanged.
  !> For "symmetric" each off-diagonal entry (i,j) gains a mirror at (j,i)
  !> with the same value.  For "skew-symmetric" the mirror is negated.
  !> Diagonal entries are never duplicated.
  subroutine expand_symmetry(symm, n_in, indx, jndx, rval, &
                              row_out, col_out, val_out, nz_out)
    character(len=*), intent(in)  :: symm
    integer,          intent(in)  :: n_in
    integer,          intent(in)  :: indx(n_in), jndx(n_in)
    real(dp),         intent(in)  :: rval(n_in)
    integer,  allocatable, intent(out) :: row_out(:), col_out(:)
    real(dp), allocatable, intent(out) :: val_out(:)
    integer,               intent(out) :: nz_out

    logical :: is_skew
    integer :: p, idx, n_offdiag

    select case (trim(symm))
    case ('general')
      nz_out = n_in
      allocate(row_out(nz_out), col_out(nz_out), val_out(nz_out))
      row_out = indx
      col_out = jndx
      val_out = rval

    case ('symmetric', 'skew-symmetric')
      is_skew = (trim(symm) == 'skew-symmetric')

      n_offdiag = 0
      do p = 1, n_in
        if (indx(p) /= jndx(p)) n_offdiag = n_offdiag + 1
      end do
      nz_out = n_in + n_offdiag

      allocate(row_out(nz_out), col_out(nz_out), val_out(nz_out))

      idx = 0
      do p = 1, n_in
        idx = idx + 1
        row_out(idx) = indx(p)
        col_out(idx) = jndx(p)
        val_out(idx) = rval(p)
        if (indx(p) /= jndx(p)) then
          idx = idx + 1
          row_out(idx) = jndx(p)
          col_out(idx) = indx(p)
          val_out(idx) = merge(-rval(p), rval(p), is_skew)
        end if
      end do

    case default
      write(error_unit,'(3a)') &
          'Error: unsupported symmetry "', trim(symm), &
          '"; supported: general, symmetric, skew-symmetric.'
      stop 1
    end select
  end subroutine expand_symmetry

  ! ---------------------------------------------------------------------------
  !> Write the solution, residual, norms, timing, and optional diagnostics
  !> to out_unit.
  !>
  !> x and r are printed side-by-side, one row per index.  When n > 6 and
  !> print_all is .false., only the first and last 3 elements are shown with
  !> a ":" separator; set print_all = .true. to print every element.
  subroutine print_results(out_unit, verbose, print_all, n, x, resid, &
                           res_inf, res_2norm, err_inf, err_2norm, &
                           elapsed, aflag, iflag, ifail)
    integer,  intent(in) :: out_unit, n, iflag(10), ifail
    logical,  intent(in) :: verbose, print_all
    real(dp), intent(in) :: x(n), resid(n)
    real(dp), intent(in) :: res_inf, res_2norm, err_inf, err_2norm
    real(dp), intent(in) :: elapsed, aflag(8)

    integer :: i, k
    integer, parameter :: NHALF = 3
    logical :: full

    full = print_all .or. n <= 2 * NHALF

    write(out_unit,'(a)') &
        '     i       x[i]              err(x)          r[i]'
    write(out_unit,'(a)') &
        '  ----  ----------------  ------------  ----------------'

    do i = 1, n
      ! In preview mode skip the middle block; print a single ":" instead.
      if (.not. full .and. i == NHALF + 1) then
        write(out_unit,'(a)') '     :'
      end if
      if (full .or. i <= NHALF .or. i > n - NHALF) then
        k = i
        write(out_unit,'(2x,i6,2x,es16.8,2x,es12.4,2x,es16.8)') &
            k, x(k), x(k) - 1.0_dp, resid(k)
      end if
    end do

    write(out_unit,'(a)') ''
    write(out_unit,'(a,es14.6)') '||r||_inf   = ', res_inf
    write(out_unit,'(a,es14.6)') '||r||_2     = ', res_2norm
    write(out_unit,'(a,es14.6)') '||x - 1||_inf = ', err_inf
    write(out_unit,'(a,es14.6)') '||x - 1||_2   = ', err_2norm

    write(out_unit,'(a)') ''
    write(out_unit,'(a,f14.6,a)') 'Wall time: ', elapsed, ' s'

    if (verbose) then
      write(out_unit,'(a)') ''
      write(out_unit,'(a)') 'Solver diagnostics (set by y12ma):'
      write(out_unit,'(2x,a,es14.6)') 'Stability factor   AFLAG(1) = ', aflag(1)
      write(out_unit,'(2x,a,es14.6)') 'Drop tolerance     AFLAG(2) = ', aflag(2)
      write(out_unit,'(2x,a,es14.6)') 'Growth limit       AFLAG(3) = ', aflag(3)
      write(out_unit,'(2x,a,es14.6)') 'Pivot tolerance    AFLAG(4) = ', aflag(4)
      write(out_unit,'(2x,a,es14.6)') 'Growth factor      AFLAG(5) = ', aflag(5)
      write(out_unit,'(2x,a,es14.6)') 'Max element in A   AFLAG(6) = ', aflag(6)
      write(out_unit,'(2x,a,es14.6)') 'Max element in LU  AFLAG(7) = ', aflag(7)
      write(out_unit,'(2x,a,es14.6)') 'Min pivot          AFLAG(8) = ', aflag(8)
      write(out_unit,'(2x,a,i0)')     'Row garbage coll.  IFLAG(6) = ', iflag(6)
      write(out_unit,'(2x,a,i0)')     'Col garbage coll.  IFLAG(7) = ', iflag(7)
      write(out_unit,'(2x,a,i0)')     'Max nnz in A       IFLAG(8) = ', iflag(8)
    end if

    write(out_unit,'(a,i0)') 'IFAIL:     ', ifail
  end subroutine print_results

end module y12m_mm_driver


! =============================================================================
! Main program — command-line parsing only.
! =============================================================================
program y12m_mm
  use y12m_mm_driver, only: run, cli_opts
  implicit none

  type(cli_opts) :: opts
  call parse_args(opts)
  call run(opts)

contains

  !> Parse command-line arguments and populate opts.
  subroutine parse_args(opts)
    type(cli_opts), intent(out) :: opts

    integer :: iarg, nargs
    character(len=512) :: arg

    nargs = command_argument_count()
    iarg  = 1
    do while (iarg <= nargs)
      call get_command_argument(iarg, arg)
      select case (trim(arg))
      case ('-o', '--output')
        iarg = iarg + 1
        if (iarg > nargs) then
          write(*,'(a)') 'Error: -o/--output requires a filename argument.'
          stop 1
        end if
        call get_command_argument(iarg, opts%outfile)
        opts%to_file = .true.
      case ('-v', '--verbose')
        opts%verbose = .true.
      case ('-p', '--print-all')
        opts%print_all = .true.
      case ('-f', '--fill-factor')
        iarg = iarg + 1
        if (iarg > nargs) then
          write(*,'(a)') 'Error: -f/--fill-factor requires an integer argument.'
          stop 1
        end if
        call get_command_argument(iarg, arg)
        read(arg, *, iostat=ios) opts%fill_factor
        if (ios /= 0 .or. opts%fill_factor < 2) then
          write(*,'(a)') 'Error: --fill-factor must be an integer >= 2.'
          stop 1
        end if
      case ('-h', '--help')
        call print_help()
        stop 0
      case ('-')
        opts%use_stdin = .true.
        opts%infile = '<stdin>'
      case default
        if (arg(1:1) == '-') then
          write(*,'(2a)') 'Error: unknown option: ', trim(arg)
          write(*,'(a)')  'Try: y12m_mm --help'
          stop 1
        end if
        if (len_trim(opts%infile) == 0 .or. trim(opts%infile) == '<stdin>') then
          opts%infile = trim(arg)
        else
          write(*,'(2a)') 'Warning: ignoring unexpected argument: ', trim(arg)
        end if
      end select
      iarg = iarg + 1
    end do

    if (len_trim(opts%infile) == 0) then
      write(*,'(a)') 'Error: no input file specified.'
      write(*,'(a)') 'Try: y12m_mm --help'
      stop 1
    end if

    if (trim(opts%infile) == '-') opts%use_stdin = .true.
  end subroutine parse_args

  !> Print the help message.
  subroutine print_help()
    write(*,'(a)') 'Usage: y12m_mm [options] <input.mtx>'
    write(*,'(a)') '       y12m_mm [options] -      (read from stdin)'
    write(*,'(a)') ''
    write(*,'(a)') 'Solve Ax = b where b = A*[1,...,1]^T (row sums).'
    write(*,'(a)') 'The exact solution is therefore x = [1,...,1]^T.'
    write(*,'(a)') 'Uses the y12ma high-level solver (double precision).'
    write(*,'(a)') ''
    write(*,'(a)') 'Options:'
    write(*,'(a)') '  -o, --output FILE   Write results to FILE (default: stdout)'
    write(*,'(a)') '  -v, --verbose       Print solver diagnostics (AFLAG, IFLAG)'
    write(*,'(a)') '  -p, --print-all     Print full x and r vectors (default: head/tail preview)'
    write(*,'(a)') '  -f, --fill-factor N Over-provision solver arrays: nn = N*nnz (default: 3)'
    write(*,'(a)') '  -h, --help          Show this help and exit'
    write(*,'(a)') ''
    write(*,'(a)') 'Supported Matrix Market formats:'
    write(*,'(a)') '  representation: coordinate'
    write(*,'(a)') '  field:          real, integer'
    write(*,'(a)') '  symmetry:       general, symmetric, skew-symmetric'
    write(*,'(a)') ''
    write(*,'(a)') 'Note: y12ma sets AFLAG(1-4) and IFLAG(2-5) to fixed internal'
    write(*,'(a)') 'defaults.  Use the y12mb/y12mc/y12md API to tune these values.'
  end subroutine print_help

end program y12m_mm
