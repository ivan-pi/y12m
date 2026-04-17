c     y12mhf.f - complete version of the double-precision 1-norm routine.
c     The file src/legacy/y12mhf.f is missing the final END statement;
c     this file provides the complete, buildable implementation.
c     The algorithm is identical to the original.
c
      subroutine y12mhf(n,nz,a,snr,work,anorm)
c
c
c   purpose.
c   -------
c
c
c   this subroutine computes the    one-norm   of a sparse
c   matrix   a.   all parameters  (except    anorm )   have  the
c   same meaning as in the other subroutines in  package   y12m.
c   on exit the   one-norm   of matrix   a   will be stored in
c   anorm.
c
c
c
c
c  declaration of the global variables and arrays.
c
c
      integer n, nz
      integer   snr
      double precision a, work, anorm
      dimension a(nz), snr(nz), work(n)
c
c
c  declaration of the internal variables.
c
c
      integer l, i
c
c
c  set all locations of array     work     equal to zero.
c
c
      do 10 i=1,n
   10 work(i)=0.0d0
c
c
c  calculate the sums of the absolute values of the non-zero
c  elements in each column of matrix     a .     store these sums
c  in array     work .
c
c
      do 20 i=1,nz
      l=snr(i)
   20 work(l)=work(l)+abs(a(i))
c
c
c  calculate the one-norm of matrix     a .
c
c
      anorm=0.0d0
      do 30 i=1,n
   30 if(work(i).gt.anorm) anorm=work(i)
c
c
c  stop the computations.
c
c
      return
      end
