#! SPDX-License-Identifier: GPL-2.0-only
#! Assisted-by: Claude Code
#:set postfix = ['e', 'f']
#:set type = ['real', 'double precision']
#:set stab = ['16.0e0', '16.0d0']
#:set drop = ['1.0e-12', '1.0d-12']
#:set grow = ['1.0e+16', '1.0d+16']
#:for pf, real_t, stab, drop, grow in zip(postfix,type,stab,drop,grow)
      subroutine y12ma${pf}$(n, z, a, snr, nn, rnr, nn1, pivot, ha,
     1iha,aflag,iflag,b,ifail)
      implicit ${real_t}$ (a-b,g,p,t-y), integer (c,f,h-n,r-s,z)
      ${real_t}$ a(nn), pivot(n), aflag(8), b(n)
      integer snr(nn), rnr(nn1), ha(iha,11), iflag(10)
      external y12mb${pf}$, y12mc${pf}$, y12md${pf}$
      aflag(1)=${stab}$
      aflag(2)=${drop}$
      aflag(3)=${grow}$
      aflag(4)=${drop}$
      iflag(2)=2
      iflag(3)=1
      iflag(4)=0
      iflag(5)=1
      call y12mb${pf}$(n,z,a,snr,nn,rnr,nn1,ha,iha,aflag,iflag,ifail)
      if(ifail.ne.0)go to 1
      call y12mc${pf}$(n,z,a,snr,nn,rnr,nn1,pivot,b,ha,iha,aflag,iflag,
     1 ifail)
      if(ifail.ne.0)go to 1
      call y12md${pf}$(n,a,nn,b,pivot,snr,ha,iha,iflag,ifail)
1     return
      end
c
#:endfor
