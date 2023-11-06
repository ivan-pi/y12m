! test_y12m.f -- Solve a 3-by-3 sparse matrix
!
! This test is derived from the Scipy documentation at
! https://docs.scipy.org/doc/scipy/reference/sparse.html#matrix-vector-product
!
program test_y12m

implicit none

integer, parameter :: n = 3
integer, parameter :: nnz = 5

real, allocatable :: a(:), b(:), pivot(:), aa(:,:), x(:)
integer, allocatable :: snr(:), rnr(:), ha(:,:)

integer :: iflag(10), ifail, nn, nn1, iha, z, i
real :: aflag(8)

interface y12ma
    subroutine y12mae(n, z, a, snr, nn, rnr, nn1, &
        pivot, ha, iha, aflag, iflag, b, ifail)
      integer, intent(in) :: n
      integer, intent(inout) :: z
      real, intent(inout) :: a(nn)
      integer, intent(inout) :: snr(nn)
      integer, intent(in) :: nn
      integer, intent(inout) :: rnr(nn1)
      integer, intent(in) :: nn1
      real, intent(inout) :: pivot(n)
      integer, intent(inout) :: ha(iha,11)
      integer, intent(in) :: iha
      real, intent(inout) :: aflag(8)
      integer, intent(inout) :: iflag(10)
      real, intent(inout) :: b(n)
      integer, intent(out) :: ifail
    end subroutine
end interface

iflag = 0
aflag = 0

z = nnz
nn  = 2*nnz
nn1 = 2*nnz
iha = n

! real arrays
allocate(a(nn),b(n),x(n),pivot(n), aa(n,n))

! integer arrays
allocate(snr(nn),rnr(nn1),ha(iha,11))

! store nonzeros
a(1:z) = [1.0, 2.0,      &
                    3.0, &
          4.0,      5.0]

! column indices
snr(1:z) = [1,2,3,1,3]

! row indices (must be colum ordered)
rnr(1:z) = [1,1,2,3,3]

! dense matrix
aa = 0
do i = 1, z
   aa(rnr(i),snr(i)) = a(i)
end do

! Right-hand side
b = [1.0,-3.0,-1.0]
x = b

call y12mae(n,z,a,snr,nn,rnr,nn1,pivot,ha,iha,aflag,iflag,x,ifail)
!call display_flags(aflag,iflag)

if (ifail /= 0) then
   write(*,*) "IFAIL = ", ifail
   write(*,*) "AFLAG(5) - Growth factor: ", aflag(5)
   write(*,*) "AFLAG(8) - Minimal pivot: ", aflag(8)
   error stop ifail
end if

! Expected solution is [1., 0., -1.]

write(*,*) "Solution is", x

write(*,*) "2-norm   = ", norm2(matmul(aa,x) - b)
write(*,*) "Max-norm = ", maxval(abs(matmul(aa,x) - b))


deallocate(a,b,pivot,snr,rnr,ha)

contains

   subroutine display_flags(aflag,iflag)
      real, intent(in) :: aflag(8)
      integer, intent(in) :: iflag(10)
      integer :: i 
      write(*,*) "-------------------------------------"
      do i = 1, 8
         write(*,'("AFLAG(",I0,") = ",G0)') i, aflag(i)
      end do
      write(*,*) "-------------------------------------"
      do i = 1, 10
         write(*,'("IFLAG(",I2,") = ",I0)') i, iflag(i)
      end do
      write(*,*) "-------------------------------------"
   end subroutine

end program