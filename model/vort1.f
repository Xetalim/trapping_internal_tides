      subroutine vort1(m,mm,nn)
c
c --- version 2.8.1 (biharmonic)
      implicit none
c
      include 'dimensions.h'
      include 'common_blocks.h'
c
      real dpmx(idm,jdm),cutoff,vort(idm,jdm)
      real hh2(idm,jdm)
      integer ib,jb
c
      cutoff=onem
c
      do 9 k=1,kk
      km=k+mm
      kn=k+nn
c
c --- store total (barotropic plus baroclinic) flow at old and mid time in
c --- -utotn,vtotn- and -utotm,vtotm- respectively. store minimum thickness
c --- values for use in pot.vort. calculation in -dpmx-.
c
      do 812 j=1,jj
      do 812 i=1,ii
      vort(i,j)=huge			!  diagnostic use
 812  dpmx(i,j)=2.*cutoff
c
       do j=1,jj
       do i=1,ii
      utotm(i,j)=u(i,j,km)+ubavg(i,j,m)
      vtotm(i,j)=v(i,j,km)+vbavg(i,j,m)
      enddo
      enddo
      
      do 807 j=2,jj
      do 807 i=2,ii
      dpmx(i,j  )=max(dpmx(i,j  ),dp(i,j,km)+dp(i-1,j,km))
 807  dpmx(i,j)=max(dpmx(i,j),dp(i,j-1,km)+dp(i-1,j-1,km))
c
      do 808 j=2,jj
      do 808 i=2,ii
      dpmx(i  ,j)=max(dpmx(i  ,j),dp(i,j,km)+dp(i,j-1,km))
 808  dpmx(i,j)=max(dpmx(i,j),dp(i-1,j,km)+dp(i-1,j-1,km))
c
      do 707 i=2,ii
 707  dpmx(i,1  )=max(dpmx(i,1  ),dp(i,1,km)+dp(i-1,1,km))
c
      do 708 i=2,ii
      dpmx(i  ,1)=max(dpmx(i  ,1),dp(i,1,km)+dp(i,1,km))
 708  dpmx(i,1)=max(dpmx(i,1),dp(i-1,1,km)+dp(i-1,1,km))
c
      do 607 j=2,jj
      dpmx(1,j  )=max(dpmx(1,j  ),dp(1,j,km))
 607  dpmx(1,j)=max(dpmx(1,j),dp(1,j-1,km))
c
      do 608 j=2,jj
 608  dpmx(1  ,j)=max(dpmx(1  ,j),dp(1,j,km)+dp(1,j-1,km))
c

c --- vorticity, pot.vort., defor. at interior points (incl. promontories)
      do 63 j=2,jj
      do 63 i=2,ii
  63  vort(i,j)=(vtotm(i,j)*scuv(i,j)-vtotm(i-1,j)*scuv(i-1,j)
     .          -utotm(i,j)*scuv(i,j)+utotm(i,j-1)*scuv(i,j-1))
     .          *sc2i(i,j)
      do 64 j=2,jj
      jb=j+1
      if(j.eq.jj) jb=jj
      do 64 i=2,ii
      ib=i+1
      if(i.eq.ii) ib=ii
      hh2(i,j)=0.125*max(8.*cutoff,2.*(dp(i,j,km)+dp(i-1,j,km)+
     .                      dp(i,j-1,km)+dp(i-1,j-1,km))
     .   ,dpmx(i,j),dpmx(i-1,j),dpmx(ib,j),dpmx(i,j-1),dpmx(i,jb))
 64   potvor(i,j,k)=(vort(i,j)+corio(i,j)) / hh2(i,j)

      do 65 i=2,ii
      ib=i+1
      if(i.eq.ii) ib=ii
      hh2(i,1)=0.125*max(8.*cutoff,4.*(dp(i,1,km)+dp(i-1,1,km))
     .   ,dpmx(i,1),dpmx(i-1,1),dpmx(ib,1),dpmx(i,2))
 65   potvor(i,1,k)=corio(i,1) / hh2(i,1)

      do 66 j=2,jj
      jb=j+1
      if(j.eq.jj) jb=jj
      vort(1,j)=(vtotm(1,j)*scuv(1,j))
     .          *sc2i(1,j)
      hh2(1,j)=0.125*max(8.*cutoff,2.*(dp(1,j,km)+
     .                      dp(1,j-1,km))
     .   ,dpmx(1,j),dpmx(2,j),dpmx(1,j-1),dpmx(1,jb))
 66   potvor(1,j,k)=(vort(1,j)+corio(1,j)) / hh2(1,j)

  9   continue

      return
      end
