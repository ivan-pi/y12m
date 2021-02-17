      parameter (iha=10,nn=200,nn1=100,iha1=60)
      implicit real(a-b,g,p,t-v),integer(c,f,h-n,r-s,z)
      real a(nn),y(iha),b(iha),aflag(11),b1(iha),x(iha),a1(iha1)
      integer snr(nn),rnr(nn1),ha(iha,13),iflag(12),sn(iha1)
      data nin/5/, nout/6/
c
c     initialization of the parameters.
c
      aflag(1)=128.0
      aflag(2)=1.e-3
      aflag(3)=1.e+16
      aflag(4)=1.e-12
      iflag(2)=2
      iflag(3)=1
      iflag(4)=1
      iflag(5)=2
      iflag(11)=25
      read (nin,101)n,z
101   format(2i4)
c
c     initialize the non-zero elements of matrix  a  in arbitrary order.
c
      do 120 k=1,z
      read(nin,110) rnr(k), snr(k), a(k)
110   format(2i4,f12.6)
120   continue
c
c
c     initialize the components of the right-hand side vector  b.
c
      read(nin,130) (b(k),k=1,n)
130   format(6f12.6)
c
c     call the subroutine  y12mfe.
c
      call y12mfe(n,a,snr,nn,rnr,nn1,a1,sn,z,ha,iha,b,b1,x,y,
     1 aflag,iflag,ifail)
c
c     print the results.
c
      write(nout,1)ifail
    1 format('1the error diagnostic parameter  ifail  is equal to',i4)
      if(ifail.gt.0)go to 5
      write(nout,2)
    2 format('0the solution vector is given below.')
      do 4 i=1,6
      write(nout,3)i,x(i)
    3 format(' ',i10,f20.5)
    4 continue
c
c     print the auxiliary information about the solution.
c     this is optionally.
c
      write(nout,11)aflag(6)
   11 format('0the largest element in the original matrix is:'f12.5)
      write(nout,12)aflag(7)
   12 format('0the largest element found in the elimination',f12.5)
      write(nout,13)aflag(8)
   13 format('0the minimal(in absolute value)pivotal element',f12.5)
      write(nout,14)aflag(5)
   14 format('0the growth factor is:',1pe12.2)
      write(nout,15)iflag(6)
   15 format('0the number of collections in the row list',i5)
      write(nout,16)iflag(7)
   16 format('0the number of collections in the column list',i5)
      write(nout,17)iflag(8)
   17 format('0the largest number of elements found in array a',i9)
      write(nout,18) iflag(12)
   18 format('0the number of iterations is:,'i5)
      write(nout,19) aflag(9)
   19 format('0estimation of the absolute error is:',1pe10.2)
    5 continue
      stop
      end
