      parameter (iha = 10,nn = 200,nn1 = 100)
      implicit real (a-b,g,p,t-v),integer(c,f,h-n,r-s,z)
      real a(nn),pivot(iha),b(iha),aflag(8)
      integer snr(nn),rnr(nn1),ha(iha,11),iflag(10)
      data nin/5/, nout/6/
c
c     initialization of the parameters.
c
      read (nin,101)n,z
101   format(2i4)
c
c     initialize the non-zero elements of matrix
c     a  in arbitrary order.
c
      do 120 k = 1,z
      read(nin,110) rnr(k), snr(k), a(k)
110   format(2i4,f12.6)
120   continue
c
c
c     initialize the components of the right-hand side vector  b.
c
      read(nin,130) (b(k),k = 1,n)
130   format(6f12.6)
c
c     call the subroutine  y12mae.
c
      call y12mae(n,z,a,snr,nn,rnr,nn1,pivot,ha,iha,aflag,
     1 iflag,b,ifail)
c
c     print the results.
c
      write(nout,1)ifail
1     format('1the error diagnostic parameter  ifail',
     1 '  is equal to',i4)
      if(ifail.gt.0)go to 5
      write(nout,2)
    2 format('0the solution vector is given below.')
      do 4 i = 1,6
      write(nout,3)i,b(i)
    3 format(' ',i10,f20.5)
    4 continue
c
c     print the auxiliary information about the solution.
c     this is optionally.
c
      write(nout,11)aflag(6)
11    format('0the largest element in the original',
     1 ' matrix is:'f12.5)
      write(nout,12)aflag(7)
   12 format('0the largest element found in the elimination',
     1f12.5)
      write(nout,13)aflag(8)
   13 format('0the minimal(in absolute value)pivotal element',
     1 f12.5)
      write(nout,14)aflag(5)
   14 format('0the growth factor is:',1pd12.2)
      write(nout,15)iflag(6)
   15 format('0the number of collections in the row list',i5)
      write(nout,16)iflag(7)
   16 format('0the number of collections in the column list',i5)
      write(nout,17)iflag(8)
   17 format('0the largest number of elements found in array a',
     1 i9)
    5 continue
      stop
      end
