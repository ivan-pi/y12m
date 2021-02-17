      subroutine matre1(n,z,c,nn,nn1,a,snr,rnr)                                 
      implicit real(a-b,g,t-v),integer(c,f,h-n,r-s,z)                           
      real a(nn)                                                                
c     integer*2 snr(nn),rnr(nn1)                                                
      integer snr(nn),rnr(nn1)                                                  
      do 1 i=1,n                                                                
      a(i)=4.                                                                   
      snr(i)=i                                                                  
    1 rnr(i)=i                                                                  
      r=n-1                                                                     
      do 2 i=1,r                                                                
      r1=n+i                                                                    
      a(r1)=-1.                                                                 
      snr(r1)=i+1                                                               
    2 rnr(r1)=i                                                                 
      r2=2*n-1                                                                  
      do 3 i=1,r                                                                
      r1=r2+i                                                                   
      a(r1)=-1.                                                                 
      snr(r1)=i                                                                 
    3 rnr(r1)=i+1                                                               
      r=n-c                                                                     
      r2=3*n-2                                                                  
      do 4 i=1,r                                                                
      r1=r2+i                                                                   
      a(r1)=-1.                                                                 
      snr(r1)=i+c                                                               
    4 rnr(r1)=i                                                                 
      r2=4*n-2-c                                                                
      do 5 i=1,r                                                                
      r1=r2+i                                                                   
      a(r1)=-1.                                                                 
      snr(r1)=i                                                                 
    5 rnr(r1)=i+c                                                               
      z=5*n-2*c-2                                                               
      return                                                                    
      end                                                                       
