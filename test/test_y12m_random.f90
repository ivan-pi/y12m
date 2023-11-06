! test_y12m_random.f -- Test solve of a sparse matrix with random values
!
! Only the diagonal and a section of the first row are filled
!
! The test is derived from the Scipy documentation at
! https://docs.scipy.org/doc/scipy/reference/sparse.html#example-1
!
program test_y12m_random

implicit none

integer, parameter :: dp = kind(1.0d0)

integer, parameter :: n = 1000
integer, parameter :: nnz = 1200

real(dp), allocatable :: a(:), b(:), x(:), pivot(:)
integer, allocatable :: snr(:), rnr(:),  ha(:,:)

integer :: iflag(10), ifail, z, nn, nn1, iha, i
real(dp) :: aflag(8)

external :: y12mae  ! for single precision
external :: y12maf  ! for double precision

iflag = 0
aflag = 0

z = nnz
nn  = 2*nnz
nn1 = 2*nnz
iha = n

! real arrays
allocate(a(nn),b(n),x(n),pivot(n))

! integer arrays
allocate(snr(nn),rnr(nn1),ha(iha,11))

! generate non-zeros
call random_number(a(1:nnz))

! row 1
do i = 1, 200
   snr(i) = 100 + i
   rnr(i) = 1
end do

! diagonal
do i = 1, 1000
   snr(i + 200) = i 
   rnr(i + 200) = i 
end do

! solution vector
x(1:n) = 1

! right-hand side
call mmult(n,z,a,rnr,snr,x,b)
print *, maxval(b), minval(b)

! solve system
call y12maf(n,z,a,snr,nn,rnr,nn1,pivot,ha,iha,aflag,iflag,b,ifail)
if (ifail /= 0) then
   write(*,*) "Smallest pivot: ", aflag(8)
   error stop ifail
end if

write(*,*) "2-norm of difference   = ", norm2(x - b)
write(*,*) "max-norm of difference = ", maxval(abs(x - b))

deallocate(a,b,x,pivot,snr,rnr,ha)

contains

   ! matrix multiplication in coordinate format
   subroutine mmult(n,nnz,a,ia,ja,x,b)
      implicit none
      integer, intent(in) :: n, nnz
      real(dp), intent(in) :: a(nnz), x(n)
      real(dp), intent(out) :: b(n)
      integer, intent(in) :: ia(nnz), ja(nnz)

      integer :: i, j

      b = 0
      do i = 1, nnz
         b(ia(i)) = b(ia(i)) + a(i)*x(ja(i))
      end do

   end subroutine

end program