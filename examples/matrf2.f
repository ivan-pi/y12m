      subroutine matrf2(m,n,c,index,alpha,nn,nn1,nz,a,snr,rnr,ifejlm)           
c                                                                               
c                                                                               
c   purpose.                                                                    
c   - - - - -                                                                   
c                                                                               
c                                                                               
c   the subroutine generates sparse  (rectangular or square)  matrices.         
c   the dimensions of the matrix and the average number of non-zero             
c   elements per row can be specified by the user.     moreover,                
c   the user can also change the sparsity pattern and the condition num-        
c   ber of the matrix.   the non-zero elements of the desired matrix            
c   will be accumulated  (in an arbitrary order)  in the first   nz             
c   positions of array   a.   the column and the row numbers of the non-        
c   zero element stored in   a(i),   i=1(1)nz,   will be found in               
c   snr(i)   and   rnr(i),   respectively.                                      
c                                                                               
c                                                                               
c                                                                               
c                                                                               
c   input parameters.                                                           
c                                                                               
c                                                                               
c   m       - integer.  the number of the rows in the desired matrix.           
c             n < m+1 < maxint+1    must be specified.                          
c                                                                               
c   n       - integer.  the number of the columns in the desired matrix.        
c             21 < n < maxint+1   must be specified.                            
c                                                                               
c   c       - integer.  the sparsity pattern can be changed by means of         
c             this parameter.   10 < c < n-10   must be specified.              
c                                                                               
c   index   - integer.   the average number of non-zero elements per row        
c             in the matrix will be equal to   index.                           
c             1 < index < n-c-8    must be specified.                           
c                                                                               
c   alpha   - real.   the condition number of the matrix can be changed         
c             by this parameter.   alpha > 0.0   must be specified.             
c             if   alpha   is approximately equal to   1.0   then the           
c             generated matrix is well-conditioned.   large values of           
c             alpha   will usually produce ill-conditioned matrices.            
c             note that no round-off errors during the computations in          
c             this subroutine are made if   alpha = 2**i   (where   i           
c             is an arbitrary integer which produces numbers in the             
c             machine range).                                                   
c                                                                               
c   nn      - integer.   the length of arrays   a   and   snr   (see            
c             below). index*m+109 < nn < maxint+1    must be specified.         
c                                                                               
c   nn1     - integer.   the length of array   rnr   (see below).               
c             index*m+109 < nn1 < maxint+1   must be specified.                 
c                                                                               
c                                                                               
c                                                                               
c                                                                               
c                                                                               
c   output parameters.                                                          
c                                                                               
c                                                                               
c   nz      - integer.   the number of non-zero elements in the matrix.         
c                                                                               
c   a(nn)   - real array. the non-zero elements of the matrix generated         
c             are accumulated in the first   nz   locations of array  a.        
c                                                                               
c   snr(nn) - integer*2 array.   the column number of the non-zero ele-         
c             ment kept in   a(i),   i=1(1)nz,   is stored in   snr(i).         
c                                                                               
c   rnr(nn1)- integer*2 array.   the row number of the non-zero element         
c             kept in   a(i),   i=1(1)nz,   is stored in   rnr(i).              
c                                                                               
c   ifejlm   - integer.   ifejlm=0  indicates that the call is successful.      
c             error diagnostics are given by means of positive values of        
c             this parameter as follows:                                        
c                  ifejlm=1   -   n       is out of range.                      
c                  ifejlm=2   -   m       is out of range.                      
c                  ifejlm=3   -   c       is out of range.                      
c                  ifejlm=4   -   index   is out of range.                      
c                  ifejlm=5   -   nn      is out of range.                      
c                  ifejlm=6   -   nn1     is out of range.                      
c                  ifejlm=7   -   alpha   is out of range.                      
c                                                                               
c                                                                               
c                                                                               
c                                                                               
      real a,alpha,alpha1                                                       
      integer m,n,nz,c,nn,nn1,ifejlm,m1,nz1,rr1,rr2,rr3,k,m2,n1,n2              
c     integer*2 snr,rnr                                                         
      integer snr,rnr                                                           
      dimension a(nn),snr(nn),rnr(nn1)                                          
c     maxint -- the larget integer such that          maxint                    
c     and      -maxint     are representable on the computer                    
c     maxint=32768                                                              
      maxint=int(2.0**31-1.0)                                                   
      m1=m                                                                      
      ifejlm=0                                                                  
      nz1=index*m+110                                                           
      k=1                                                                       
      alpha1=alpha                                                              
      index1=index-1                                                            
c                                                                               
c  check the parameters.                                                        
c                                                                               
      if(n.ge.22)go to 1                                                        
    2 ifejlm=1                                                                  
      return                                                                    
    1 if(n.gt.maxint)go to 2                                                    
      if(m.ge.n)go to 3                                                         
    4 ifejlm=2                                                                  
      return                                                                    
    3 if(m.gt.maxint)go to 4                                                    
      if(c.lt.11)go to 6                                                        
      if(n-c.ge.11)go to 5                                                      
    6 ifejlm=3                                                                  
      return                                                                    
    5 if(index.lt.2)go to 12                                                    
      if(n-c-index.ge.9)go to 13                                                
   12 ifejlm=4                                                                  
   13 if(nn.ge.nz1)go to 7                                                      
    8 ifejlm=5                                                                  
      return                                                                    
    7 if(nn.gt.maxint)go to 8                                                   
      if(nn1.ge.nz1)go to 9                                                     
   10 ifejlm=6                                                                  
      return                                                                    
    9 if(nn1.gt.maxint)go to 10                                                 
      if(alpha.gt.0.)go to 11                                                   
      ifejlm=7                                                                  
      return                                                                    
   11 continue                                                                  
c                                                                               
c  end of the error check.   begin to generate the non-zero elements of         
c  the required matrix.                                                         
c                                                                               
      do 20 i=1,n                                                               
      a(i)=1.                                                                   
      snr(i)=i                                                                  
   20 rnr(i)=i                                                                  
      nz=n                                                                      
      j1=1                                                                      
      do 21 j=1,index1                                                          
      j1=-j1                                                                    
      do 22 i=1,n                                                               
      a(nz+i)=j1*j*i                                                            
      if(i+c+j-1.le.n)snr(nz+i)=i+c+j-1                                         
      if(i+c+j-1.gt.n)snr(nz+i)=c+i+j-1-n                                       
   22 rnr(nz+i)=i                                                               
   21 nz=nz+n                                                                   
      rr1=10                                                                    
      rr2=nz                                                                    
      rr3=1                                                                     
   25 continue                                                                  
      do 26 i=1,rr1                                                             
      a(rr2+i)=alpha*i                                                          
      snr(rr2+i)=n-rr1+i                                                        
      rnr(rr2+i)=rr3                                                            
   26 continue                                                                  
      if(rr1.eq.1)go to 27                                                      
      rr2=rr2+rr1                                                               
      rr1=rr1-1                                                                 
      rr3=rr3+1                                                                 
      go to 25                                                                  
   27 nz=nz+55                                                                  
   29 m1=m1-n                                                                   
      alpha=1./alpha                                                            
      if(m1.le.0)go to 28                                                       
      n2=k*n                                                                    
      if(m1.ge.n)m2=n                                                           
      if(m1.lt.n)m2=m1                                                          
      l=0                                                                       
      if(m2.le.10)l=11-m2                                                       
      do 30 i=1,m2                                                              
      a(nz+i)=alpha*(k+1)                                                       
      snr(nz+i)=i+l                                                             
   30 rnr(nz+i)=n2+i                                                            
      nz=nz+m2                                                                  
      j1=1                                                                      
      do 41 j=1,index1                                                          
      j1=-j1                                                                    
      do 42 i=1,m2                                                              
      a(nz+i)=alpha*j*j1*((k+1)*i+1.)                                           
      if(i+c+j-1.le.n)snr(nz+i)=i+c+j-1                                         
      if(i+c+j-1.gt.n)snr(nz+i)=c+i+j-1-n                                       
   42 rnr(nz+i)=n2+i                                                            
   41 nz=nz+m2                                                                  
      k=k+1                                                                     
      go to 29                                                                  
   28 continue                                                                  
      alpha=1./alpha1                                                           
      rr1=1                                                                     
      rr2=nz                                                                    
   35 continue                                                                  
      do 36 i=1,rr1                                                             
      a(rr2+i)=alpha*(rr1+1-i)                                                  
      snr(rr2+i)=i                                                              
      rnr(rr2+i)=m-10+rr1                                                       
   36 continue                                                                  
      if(rr1.eq.10)go to 34                                                     
      rr2=rr2+rr1                                                               
      rr1=rr1+1                                                                 
      go to 35                                                                  
   34 nz=nz+55                                                                  
      alpha=alpha1                                                              
      return                                                                    
      end                                                                       
