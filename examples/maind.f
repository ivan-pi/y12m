      parameter (ibdim=1000, iadim=25*ibdim, icdim=20*ibdim)
c     parameter (ibdim=600, iadim=25*ibdim, icdim=20*ibdim)
c     implicit real(a-b,g,p,t-v),integer(c,f,h-n,r-s,z)
c     real a(15000),pivot(600),b(600)
c     integer*2 snr(15000),rnr(12000),ha(600,11)
c     integer snr(15000),rnr(12000),ha(600,11)
c     ibdim=600
c     iadim=25*ibdim
c     icdim=20*ibdim
c     call main0(ibdim,iadim,icdim,a,pivot,b,snr,rnr,ha)
c     stop
c     end
c     subroutine main0(ibdim,iadim,icdim,a,pivot,b,snr,rnr,ha)
c
c     this program tests the subroutine of package y12m by
c     means of matrices of class d(n,c) and/or class e(n,c).
c     the program uses data consisting of 7 integers.
c     the first of these integers (in the program it is de-
c     noted by indexp) should be equal to:
c         1     if the matrices of class d(n,c) has to be used.
c         2     if the matrices of class e(n,c) has to be used.
c         3     if the matrices of both classes has to be used.
c     the next three integers (called nstart,nincr and nend
c     in the program) are used to determine the range of the
c     parameter n. nstart should be larger than 22. nincr
c     should be larger than zero. nend should be larger than
c     or equal to nstart. if these conditions are satisfied
c     then the range of n for the systems solved is:
c            n = nstart(nincr)nend
c     the last three integers (called cstart,cincr and cend
c     in the program) are used in the program to determine
c     the range of parameter c (the sparsity of the matrix
c     is changed by c). cstart should be chosen between 2 and
c     n-13. cincr should be larger than 1. cend should be
c     smaller than n-13. if these conditions are satisfied,
c     then the range of c for the systems solved is given as:
c            c = cstart(cincr)cend
c     in the version which is distributed with the package
c     the data read in line 12 are:
c         3   250  50  600  4  40  204
c     if these data are changed,then the parameter ibdim
c     skould be set equal to or larger than the value for
c     nend.
c     parameter ibdim=600, iadim=25*ibdim, icdim=20*ibdim
      implicit real(a-b,g,p,t-v),integer(c,f,h-n,r-s,z)
      real a(iadim),pivot(ibdim),b(ibdim),aflag(8)
      real time, y12cck
c     integer*2 snr(iadim),rnr(icdim),ha(ibdim,11),iflag(10)
      integer snr(iadim),rnr(icdim),ha(ibdim,11),iflag(10)
      data nin/5/, nout /6/
c
c     call y12cqe('y12mae      ')
      iha=ibdim
      write(nout,50)
50    format(' indexp,nstart,nincr,nend,cstart,cincr,cend:')
      read(nin,*)indexp,nstart,nincr,nend,cstart,cincr,cend
      write(nout,*)indexp,nstart,nincr,nend,cstart,cincr,cend
52    format()
c52    format(10(1x,i5))
      if(indexp.gt.0.and.indexp.lt.4)go to 301
      write(nout,302)
302   format(' indexp is out of range')
      go to 54
301   continue
      if(nstart.gt.21)go to 303
      write(nout,304)
304   format(' nstart is out of range')
      go to 54
303   continue
      if(nincr.ge.1)go to 305
      write(nout,306)
306   format(' nincr is out of range')
      go to 54
305   continue
      if(nend.ge.nstart)go to 307
      write(nout,308)
308   format(' nend is out of range')
      go to 54
307   continue
      if(cstart.gt.1.and.cstart.lt.nend-13)go to 309
      write(nout,310)
310   format(' cstart is out of range')
      go to 54
309   continue
      if(cend.ge.cstart.and.cend.lt.nend-13)go to 311
      write(nout,312)
312   format(' cend is out of range')
      go to 54
311   continue
      if(cincr.gt.0)go to 313
      write(nout,314)
314   format(' cincr is out of range')
      go to 54
313   continue
      n=nstart
      c=cstart
      if(indexp.eq.2)index=2
      if(indexp.eq.1)index=1
      if(indexp.eq.3)index=1
      nn=iadim
      nn1=icdim
   51 continue
      if(n.gt.ibdim) go to 410
      if(index.eq.1)print 41,n,c
   41 format('1','matrix of class d,with parameters n=',i4,'  c=',i4)
      if(index.eq.2)print 33,n,c
   33 format('1matrix of class e: n=',i4,' c='i4)
      if(index.eq.1) call matrd1(n,z,c,nn,nn1,a,snr,rnr)
       if(index.eq.2) call matre1(n,z,c,nn,nn1,a,snr,rnr)
      write(nout,53)z
53    format('0the number of non-zero elements in the original',
     1 ' matrix is equal to:',10x,i8)
      do 1 i=1,n
      b(i)=0.0
1     continue
      do 2 i=1,z
      lrow=rnr(i)
      lcol=snr(i)
2     b(lrow)=b(lrow)+a(i)
      ifail1=0
      time=y12cck(1,ifail1)
      call y12mae(n,z,a,snr,nn,rnr,nn1,pivot,ha,iha,aflag,iflag,b,ifail)
      ifail1=0
      time=-time+y12cck(1,ifail1)
      if(ifail.ne.0)go to 3
      write(nout,110)time
110   format(' subs after y12mae: ',f12.2)
      t=0.0
      do 4 i=1,n
      tt=abs(b(i)-1.)
      if(tt.gt.t)t=tt
4     continue
      write(nout,102) t
102   format('0the largest error is: ',1pe10.2)
3     continue
      write(nout,104)iflag(8)
  104 format('0the largest number of elements in a:',i8)
      write(nout,105)iflag(6)
  105 format('0the number of collections in row list',i4)
      write(nout,106)iflag(7)
  106 format('0the number of collections in column list',i4)
      write(nout,107)aflag(6)
  107 format('0the largest element in original matrix',1pe10.2)
      write(nout,108)aflag(7)
  108 format('0the largest element in lu-matrix',1pe10.2)
      write(nout,109)aflag(5)
  109 format('0the growth factor is:',1pe10.2)
      write(nout,100)aflag(8)
  100 format('0the minimal pivotal element',1pe10.2)
      write(nout,204) aflag(2),aflag(1)
204   format('0the drop tolerance: ',1pe10.2,
     1 ',     the stability factor: ',1pe10.2)
      write(nout,103) ifail
103   format('0the error diagnostic parameter = ',i4)
      c=c+cincr
      if(c.le.cend)go to 51
      c=cstart
      n=n+nincr
      if(n.le.nend)go to 51
      if(indexp.ne.3)go to 54
      c=cstart
      n=nstart
      index=index+1
      if(index.eq.2)go to 51
54    continue
c     call y12cqe('term.       ')
      stop
410   continue
      write(nout,420) ibdim,n
420   format(' maxn = ',i4,' is greater than n = ',i4)
      stop
      end
