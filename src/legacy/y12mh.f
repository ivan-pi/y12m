      subroutine y12mhe(n,nz,a,snr,work,anorm)
c
c   purpose.
c   -------
c
c   this subroutine computes the    one-norm   of a sparse
c   matrix   a.   all parameters  (except    anorm )   have  the
c   same meaning as in the other subroutines in  package   y12m.
c   on exit the   one-norm   of matrix   a   will be stored in
c   anorm.
c
c
c  declaration of the global variables and arrays.
c
      integer n, nz
      integer   snr
      real a, work, anorm
      dimension a(nz), snr(nz), work(n)
c
c  declaration of the internal variables.
c
      integer l
      intrinsic abs
c
c  set all locations of array     work     equal to zero.
c
      do 10 i=1,n
      work(i)=0.0e0
   10 continue
c
c  calculate the sums of the absolute values of the non-zero
c  elements in each row of matrix     a .     store these sums
c  in array     work .
c
      do 20 i=1,nz
      l=snr(i)
      work(l)=work(l)+abs(a(i))
   20 continue
c
c  calculate the one-norm of matrix     a .
c
      anorm=0.0e0
      do 30 i=1,n
      if(work(i).gt.anorm) anorm=work(i)
   30 continue
c
c  stop the computations.
c
      return
      end
c
      subroutine y12mhf(n,nz,a,snr,work,anorm)
c
c   purpose.
c   -------
c
c   this subroutine computes the    one-norm   of a sparse
c   matrix   a.   all parameters  (except    anorm )   have  the
c   same meaning as in the other subroutines in  package   y12m.
c   on exit the   one-norm   of matrix   a   will be stored in
c   anorm.
c
c
c  declaration of the global variables and arrays.
c
      integer n, nz
      integer   snr
      double precision a, work, anorm
      dimension a(nz), snr(nz), work(n)
c
c  declaration of the internal variables.
c
      integer l
      intrinsic dabs
c
c  set all locations of array     work     equal to zero.
c
      do 10 i=1,n
      work(i)=0.0d0
   10 continue
c
c  calculate the sums of the absolute values of the non-zero
c  elements in each row of matrix     a .     store these sums
c  in array     work .
c
      do 20 i=1,nz
      l=snr(i)
      work(l)=work(l)+dabs(a(i))
   20 continue
c
c  calculate the one-norm of matrix     a .
c
      anorm=0.0d0
      do 30 i=1,n
      if(work(i).gt.anorm) anorm=work(i)
   30 continue
c
c  stop the computations.
c
      return
      end
c
