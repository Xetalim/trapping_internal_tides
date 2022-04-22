      subroutine tsadvc(m,n,mm,nn,k1m,k1n)
c
c --- version 2.7
      implicit none
c
      include 'dimensions.h'
      include 'common_blocks.h'
      include 'stmt_funcs.h'
c
      real pold,pmid,pnew,wlist1(idm,jdm),wlist2(idm,jdm)
c
c
      do 1 j=1,jj
      do 1 i=1,ii
      wlist1(i,j)=wts1
      wlist2(i,j)=wts2
      if (abs(dp(i,j,k1n)-dp(i,j,k1m))/
     .    max(dp(i,j,k1n)+dp(i,j,k1m),epsil).gt..25) then
        wlist1(i,j)=1.
        wlist2(i,j)=0.
      end if
 1    continue
c
      do 461 k=1,kk
      km=k+mm
      kn=k+nn
c
c --- ------------------------------------
c --- advection of thermodynamic variables in mixed layer 
c --- ------------------------------------
c
      do 461 j=1,jj
      do 461 i=1,ii
      pold=max(0.,dpold(i,j,k))
      pmid=max(0.,dp(i,j,km))
      pnew=max(0.,dp(i,j,kn))
      dp(i,j,km)=pmid*wlist1(i,j)+(pold+pnew)*wlist2(i,j)

 461  continue
c
      return
      end
