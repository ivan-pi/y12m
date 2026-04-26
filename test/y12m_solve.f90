! SPDX-License-Identifier: GPL-2.0-only
! Assisted-by: GitHub Copilot:claude-sonnet-4.5
!
! y12m_solve.f90 -- File-based test driver for the y12m sparse solver.
!
! This program reads a sparse linear system Ax=b from an ASCII text file (or
! standard input), solves it with the y12ma high-level interface (resolving to
! the double-precision subroutine y12maf), then prints the solution vector, the
! residual r = b - A*x, the 1-norm and 2-norm of r, and the wall-clock time
! spent in the solver.
!
! Usage:
!   y12m_solve [options] <input_file>
!   y12m_solve [options] -              (read from standard input)
!
! Options:
!   -o, --output FILE    Write results to FILE instead of standard output.
!   -v, --verbose        Also print solver diagnostics (growth factor, min
!                        pivot, garbage-collection counts, etc.).
!   -h, --help           Print this help and exit.
!
! Input file format:
!   Line 1:  <name>  --  <description>
!   Line 2:  <n> <kind>             (n = matrix dimension; kind = real|double)
!   Lines:   <row> <col> <value>    (non-zero triplets, any order)
!            0 0 0.0                (end-of-matrix sentinel)
!   Lines:   <b_1>                  (right-hand-side vector, one entry per line)
!            <b_2>
!            ...
!            <b_n>
!
! Note: y12ma always uses fixed internal defaults for AFLAG(1-4) and IFLAG(2-5)
! regardless of user input.  All values reported below are those set by y12ma
! itself.  To tune solver parameters, use the lower-level y12mb/y12mc/y12md
! interface instead.
!
program y12m_solve
  use, intrinsic :: iso_fortran_env, only: input_unit
  use y12m
  implicit none

  ! ---- Upper bound on non-zeros accepted from file ----
  integer, parameter :: MAXTMP = 100000

  ! ---- Temporary read buffers (fixed size, no heap allocation needed) ----
  double precision :: a_tmp(MAXTMP)
  integer          :: row_tmp(MAXTMP), col_tmp(MAXTMP)

  ! ---- Allocatable solver arrays ----
  double precision, allocatable :: a(:), pivot(:), b(:), aflag(:)
  integer,          allocatable :: snr(:), rnr(:), ha(:,:), iflag(:)

  ! ---- Saved copies for residual computation ----
  double precision, allocatable :: a_orig(:), b_orig(:), resid(:)
  integer,          allocatable :: row_orig(:), col_orig(:)

  ! ---- Scalars ----
  integer          :: n, z, nn, nn1, iha, ifail
  integer          :: inp_unit, out_unit, ios
  integer          :: i, k, iarg, nargs
  integer          :: row_i, col_i
  double precision :: val_d, res_1norm, res_2norm, elapsed

  ! system_clock counters
  integer :: t1, t2, clock_rate

  logical :: to_file, verbose, use_stdin

  character(len=512) :: infile, outfile, arg
  character(len=512) :: header_line
  character(len=64)  :: kind_str

  ! ---- Default values ----
  infile    = ''
  outfile   = ''
  to_file   = .false.
  verbose   = .false.
  use_stdin = .false.

  ! ---- Parse command-line arguments ----
  nargs = command_argument_count()
  iarg  = 1
  do while (iarg <= nargs)
    call get_command_argument(iarg, arg)
    select case (trim(arg))
    case ('-o', '--output')
      iarg = iarg + 1
      call get_command_argument(iarg, outfile)
      to_file = .true.
    case ('-v', '--verbose')
      verbose = .true.
    case ('-h', '--help')
      call print_help()
      stop 0
    case default
      if (len_trim(infile) == 0) then
        infile = trim(arg)
      else
        write(*,'(2a)') 'Warning: ignoring unexpected argument: ', trim(arg)
      end if
    end select
    iarg = iarg + 1
  end do

  if (len_trim(infile) == 0) then
    write(*,'(a)') 'Error: no input file specified (use - for stdin).'
    write(*,'(a)') 'Try: y12m_solve --help'
    stop 1
  end if

  ! ---- Open input ----
  if (trim(infile) == '-') then
    use_stdin = .true.
    inp_unit  = input_unit
  else
    open(newunit=inp_unit, file=trim(infile), status='old', &
         action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(3a)') 'Error: cannot open "', trim(infile), '"'
      stop 1
    end if
  end if

  ! ---- Read header line ----
  read(inp_unit, '(a)', iostat=ios) header_line
  if (ios /= 0) then
    write(*,'(a)') 'Error: failed to read header line.'
    stop 1
  end if

  ! ---- Read matrix dimension and precision kind ----
  read(inp_unit, *, iostat=ios) n, kind_str
  if (ios /= 0) then
    write(*,'(a)') 'Error: failed to read matrix dimension/kind line.'
    stop 1
  end if
  if (n < 2) then
    write(*,'(a,i0,a)') 'Error: matrix dimension n=', n, &
        ' is less than 2 (y12ma requires n >= 2).'
    stop 1
  end if

  ! ---- Read non-zero triplets ----
  z = 0
  do
    read(inp_unit, *, iostat=ios) row_i, col_i, val_d
    if (ios /= 0) then
      write(*,'(a)') 'Error: failed while reading triplets.'
      stop 1
    end if
    if (row_i == 0 .and. col_i == 0) exit   ! end-of-matrix sentinel
    z = z + 1
    if (z > MAXTMP) then
      write(*,'(a,i0,a)') 'Error: more than ', MAXTMP, &
          ' non-zeros; recompile with a larger MAXTMP.'
      stop 1
    end if
    row_tmp(z) = row_i
    col_tmp(z) = col_i
    a_tmp(z)   = val_d
  end do

  if (z == 0) then
    write(*,'(a)') 'Error: no non-zero entries found in matrix.'
    stop 1
  end if

  ! ---- Allocate saved copies ----
  allocate(a_orig(z), b_orig(n), row_orig(z), col_orig(z), resid(n))
  a_orig(1:z)   = a_tmp(1:z)
  row_orig(1:z) = row_tmp(1:z)
  col_orig(1:z) = col_tmp(1:z)

  ! ---- Read RHS vector ----
  allocate(b(n), pivot(n))
  do i = 1, n
    read(inp_unit, *, iostat=ios) b(i)
    if (ios /= 0) then
      write(*,'(a,i0)') 'Error: failed reading RHS entry ', i
      stop 1
    end if
  end do
  b_orig(1:n) = b(1:n)

  if (.not. use_stdin) close(inp_unit)

  ! ---- Allocate solver arrays ----
  ! NN >= 2*Z; recommended 2*Z to 3*Z.
  ! NN1 >= Z;  recommended 2*Z to 3*Z.
  nn  = 3 * z
  nn1 = 3 * z
  iha = n

  allocate(a(nn), snr(nn), rnr(nn1), ha(iha,11), aflag(8), iflag(10))

  ! Copy triplets into solver arrays.
  a(1:z)   = a_orig(1:z)
  snr(1:z) = col_orig(1:z)
  rnr(1:z) = row_orig(1:z)

  ! ---- Solve ----
  call system_clock(t1, clock_rate)
  call y12ma(n, z, a, snr, nn, rnr, nn1, pivot, ha, iha, aflag, iflag, b, ifail)
  call system_clock(t2)
  if (clock_rate > 0) then
    elapsed = dble(t2 - t1) / dble(clock_rate)
  else
    elapsed = 0.0d0
  end if

  ! ---- Open output ----
  if (to_file) then
    out_unit = 20
    open(unit=out_unit, file=trim(outfile), status='replace', &
         action='write', iostat=ios)
    if (ios /= 0) then
      write(*,'(3a)') 'Error: cannot open output file "', trim(outfile), '"'
      stop 1
    end if
  else
    out_unit = 6
  end if

  ! ---- Print matrix identification ----
  write(out_unit,'(a,a)') 'Matrix:  ', trim(header_line)
  write(out_unit,'(a,i0,a,i0,a,i0,3a)') &
      'Size:    ', n, ' x ', n, ',  nnz: ', z, &
      ',  kind: ', trim(kind_str), ' (solved in double precision)'

  ! ---- Check for solver failure ----
  if (ifail /= 0) then
    write(out_unit,'(a,i0)') 'IFAIL = ', ifail
    write(out_unit,'(a)')    'Solver failed; no solution available.'
    if (to_file) close(out_unit)
    stop 1
  end if

  ! ---- Compute residual r = b_orig - A_orig * x ----
  ! After y12ma returns, b(1:n) holds the solution x.
  resid(1:n) = b_orig(1:n)
  do k = 1, z
    resid(row_orig(k)) = resid(row_orig(k)) - a_orig(k) * b(col_orig(k))
  end do

  res_1norm = sum(abs(resid(1:n)))
  res_2norm = sqrt(sum(resid(1:n)**2))

  ! ---- Print solution ----
  write(out_unit,'(a)') ''
  write(out_unit,'(a)') 'Solution x:'
  do i = 1, n
    write(out_unit,'(2x,a,i0,a,es22.14)') 'x[', i, '] = ', b(i)
  end do

  ! ---- Print residual ----
  write(out_unit,'(a)') ''
  write(out_unit,'(a)') 'Residual r = b - A*x:'
  do i = 1, n
    write(out_unit,'(2x,a,i0,a,es22.14)') 'r[', i, '] = ', resid(i)
  end do

  ! ---- Print norms ----
  write(out_unit,'(a)') ''
  write(out_unit,'(a,es14.6)') '||r||_1 = ', res_1norm
  write(out_unit,'(a,es14.6)') '||r||_2 = ', res_2norm

  ! ---- Print timing ----
  write(out_unit,'(a)') ''
  write(out_unit,'(a,f14.6,a)') 'Wall time: ', elapsed, ' s'

  ! ---- Print solver diagnostics ----
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

  if (to_file) close(out_unit)

contains

  subroutine print_help()
    write(*,'(a)') 'Usage: y12m_solve [options] <input_file>'
    write(*,'(a)') '       y12m_solve [options] -      (read from stdin)'
    write(*,'(a)') ''
    write(*,'(a)') 'Solve the sparse linear system Ax=b defined in <input_file>'
    write(*,'(a)') 'using the y12ma high-level solver (double precision).'
    write(*,'(a)') ''
    write(*,'(a)') 'Options:'
    write(*,'(a)') '  -o, --output FILE   Write results to FILE (default: stdout)'
    write(*,'(a)') '  -v, --verbose       Print solver diagnostics (AFLAG, IFLAG)'
    write(*,'(a)') '  -h, --help          Show this help and exit'
    write(*,'(a)') ''
    write(*,'(a)') 'Input file format:'
    write(*,'(a)') '  Line 1:  <name>  --  <description>'
    write(*,'(a)') '  Line 2:  <n> <kind>          (n = dimension; kind = real|double)'
    write(*,'(a)') '  Lines:   <row> <col> <value>  (sparse triplets, any order)'
    write(*,'(a)') '           0 0 0.0               (end-of-matrix sentinel)'
    write(*,'(a)') '  Lines:   <b_1> ... <b_n>      (one RHS entry per line)'
    write(*,'(a)') ''
    write(*,'(a)') 'Note: y12ma sets AFLAG(1-4) and IFLAG(2-5) to fixed internal'
    write(*,'(a)') 'defaults.  Use the y12mb/y12mc/y12md API to tune these values.'
  end subroutine print_help

end program y12m_solve
