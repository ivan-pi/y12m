      parameter (ibdim=2000,iadim=50*ibdim,icdim=50*ibdim,
     * iddim=10*ibdim,maxsto=1)
      implicit real(a-b,g,p,t-y),integer(c,f,h-n,r-s,z)
      real a(iadim),y(ibdim),b(ibdim),aflag(11),b1(ibdim),w(ibdim)
      real a1(iddim),x(ibdim),tm,anorm,rcond,cond
c     integer*2 snr(iadim),rnr(icdim),ha(ibdim,13),iflag(12)
      integer snr(iadim),rnr(icdim),ha(ibdim,13),iflag(12)
c     integer*2 sn(iddim)
      integer sn(iddim)
      data nin/5/,nout/6/
      lstop=0
      iter=0
      m=200
      n=200
      c=20
      index=16
      nn=iadim
      nn1=icdim
      alpha=1.0
      aflag(1)=4.0
      aflag(2)=1.0
      aflag(3)=1.e+16
      aflag(4)=1.e-12
      iflag(2)=3
      iflag(3)=1
      iflag(4)=1
      write(nout,50)
50    format(' m, n, c, index, alpha:')
      read(nin,*) m, n, c, index, alpha
      write(nout,*) m, n, c, index, alpha
52    format()
c52    format(10(1x,i5))
      write(nout,60)
60    format(' (aflag(j),j=1,4),(iflag(j),j=2,4):')
      read(nin,*)(aflag(j),j=1,4),(iflag(j),j=2,4)
520   format(4(1x,f8.2),5i5)
521   format(4(1x,e8.2),5i5)
      write(nout,521)(aflag(j),j=1,4),(iflag(j),j=2,4)
   51 continue
      iflag(5)=2
      call matrf2(m,n,c,index,alpha,nn,nn1,nz,a,snr,rnr,ifejlm)
      if(ifejlm.gt.0)go to 531
      do 1 i=1,m
      b(i)=0.0
1     continue
      do 2 i=1,nz
      lrow=rnr(i)
      lcol=snr(i)
2     b(lrow)=b(lrow)+a(i)
      iflag(11)=64
      if(nz.le.32000)go to 33
      ifejlm=8
      go to 53
33    continue
      nstart=1
      if(m.eq.n) go to 300
      do 301 i=1,nz
      a(nz+i)=a(i)
      snr(nz+i)=rnr(i)
      snr(i)=snr(i)+m
  301 rnr(nz+i)=snr(i)
      nz=2*nz
      do 302 i=1,m
      a(nz+i)=1.0
      snr(nz+i)=i
  302 rnr(nz+i)=i
      nz=nz+m
      n=n+m
      nstart=m+1
      do 303 i=nstart,n
  303 b(i)=0.0
  300 call y12mhe(n,nz,a,snr,w,anorm)
      iha=ibdim
      print 37,m,n,c,index,alpha
   37 format('1','matrix of class     f2     with parameters   m=',i6,'
     1      n=',i6,'       c=',i6,'    index=',i6,'    alpha=',f16.1)
      if(iter.eq.0) print 38
   38 format('0','the system is solved directly')
      if(iter.eq.1) print 39
   39 format('0','the system is solved by iterative refinement')
      print 55,nz
   55 format('0','the number of non-zeros before the elimination:',i12)
      print 56,anorm
   56 format('0','the 1-norm of matrix    a   is equal to:',1pe12.2)
      call time(itime1)
      if(iter.eq.1) go to 31
      call y12mbe(n,nz,a,snr,nn,rnr,nn1,ha,iha,aflag,iflag,ifail)
      call y12mce(n,nz,a,snr,nn,rnr,nn1,y,b,ha,iha,
     *            aflag,iflag,ifail)
      call y12mde(n,a,nn,b,y,snr,ha,iha,iflag,ifail)
      go to 32
   31 call y12mfe(n,a,snr,nn,rnr,nn1,a1,sn,nz,ha,iha,b,b1,x,y,
     1 aflag,iflag,ifail)
   32 call time(itime2)
      if(ifail.ne.0)go to 3
      call y12mge(n,nn,a,snr,w,y,anorm,rcond,iha,ha,iflag,ifail)
      tm=float(itime2-itime1)/1.e+6
      write(6,210)tm
  210 format('0','the computing time is equal to:',5x,f9.2)
      cond=1.0/rcond
      write(6,57) cond
   57 format('0','the eestimate of the condition number of matrix    a
     1  is       cond =',1pe12.2)
      t=0.0
      do 4 i=nstart,n
      if(iter.eq.1) tt=abs(x(i)-1.)
      if(iter.eq.0) tt=abs(b(i)-1.)
      if(tt.gt.t)t=tt
4     continue
      write(6,102) t
102   format('0the largest error is: ',1pe10.2)
3     continue
      write(6,104)iflag(8)
  104 format('0the largest number of elements in a:',i8)
      write(6,105)iflag(6)
  105 format('0the number of collections in row list',i4)
      write(6,106)iflag(7)
  106 format('0the number of collections in column list',i4)
      write(6,107)aflag(6)
  107 format('0the largest element in original matrix',1pe10.2)
      write(6,108)aflag(7)
  108 format('0the largest element in lu-matrix',1pe10.2)
      write(6,109)aflag(5)
  109 format('0the growth factor is:',1pe10.2)
      write(6,100)aflag(8)
  100 format('0the minimal pivotal element',1pe10.2)
      if(iter.eq.0) go to 401
      write(6,202)iflag(12)
  202 format('0the number of iterations is:',i4)
      write(6,902)iflag(11)
  902 format('0the number of allowed iterations is:',i4)
      write(6,203)aflag(9),aflag(10),aflag(11)
  203 format('0the norms of the last correction , the residual  and
     1the solution vector:',3(4x,1pe14.6))
  401 write(6,204) aflag(2),aflag(1)
204   format('0the drop tolerance: ',1pe10.2,
     1 ',     the stability factor: ',1pe10.2)
      write(6,103) ifail
103   format('0the error diagnostic parameter = ',i4)
531   continue
      write(6,903) ifejlm
  903 format('0the error diagnostic parameter for,the matrix generator
     1is          ifejlm =',i3)
53    continue
      if(n.gt.m) n=n-m
      iter=iter+1
      if(iter.lt.2) go to 51
c      iter=0
c      aflag(2)=aflag(2)/10.0
c      lstop=lstop+1
c      if(lstop.eq.maxsto) aflag(2)=0.0
c      if(lstop.le.maxsto) go to 51
      stop
      end
