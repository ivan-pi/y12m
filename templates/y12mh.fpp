#:set postfix = ['e','f']
#:set type = ['real', 'double precision']
#:set literal_zero = ['0.0e0', '0.0d0']
#:set abs_func = ['abs','dabs']
#:for pf, real_t, zero, abs in zip(postfix,type,literal_zero,abs_func)
      subroutine y12mh${pf}$(n,nz,a,snr,work,anorm)
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
      ${real_t}$ a, work, anorm
      dimension a(nz), snr(nz), work(n)
c
c  declaration of the internal variables.
c
      integer l
      intrinsic ${abs}$
c
c  set all locations of array     work     equal to zero.
c
      do 10 i=1,n
      work(i)=${zero}$
   10 continue
c
c  calculate the sums of the absolute values of the non-zero
c  elements in each row of matrix     a .     store these sums
c  in array     work .
c
      do 20 i=1,nz
      l=snr(i)
      work(l)=work(l)+${abs}$(a(i))
   20 continue
c
c  calculate the one-norm of matrix     a .
c
      anorm=${zero}$
      do 30 i=1,n
      if(work(i).gt.anorm) anorm=work(i)
   30 continue
c
c  stop the computations.
c
      return
      end
c
#:endfor
