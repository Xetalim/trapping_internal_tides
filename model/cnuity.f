      subroutine cnuity(m,mm,nn)
c
c --- version 2.8
      implicit none
c
      include 'dimensions.h'
      include 'common_blocks.h'
c
      integer ia,ib,ja,jb,ka,kb
      real q,clip,flxhi,flxlo
      real p0(idm,jdm,kdm+1),p1(idm,jdm,kdm+1)
      real pmax,pmin
      real uflux4(idm,jdm),vflux4(idm,jdm)
      real uflux3(idm,jdm),vflux3(idm,jdm)
      real ph(idm,jdm,kdm+1)
c
c --- ------------------------------------------------------
c --- continuity equation (flux-corrected transport version)
c --- ------------------------------------------------------
c
      do 41 j=1,jj
c
      do 74 i=1,ii
 74   util3(i,j)=0.
c
      do 40 i=1,ii
      uflux(i,j)=0.
      uflux2(i,j)=0.
      uflux4(i,j)=0.
      uflux3(i,j)=0.
      utotm(i,j)=0.
 40   utotn(i,j)=0.
c
      do 41 i=1,ii
      p(i,j,1)=0.
      p0(i,j,1)=0.
      p1(i,j,1)=0.
      ph(i,j,1)=0.
      vflux(i,j)=0.
      vflux2(i,j)=0.
      vflux4(i,j)=0.
      vflux3(i,j)=0.
      vtotm(i,j)=0.
 41   vtotn(i,j)=0.
c
      do 76 k=1,kk
      km=k+mm
      kn=k+nn
c
c --- uflux/vflux = low-order (diffusive) mass fluxes at old time level.
c --- uflux2/vflux2 = 'antidiffusive' fluxes, defined as high-order minus low-
c --- order fluxes. high-order fluxes are second-order in space, time-centered.
c
      do 11 j=1,jj1
      do 11 i=2,ii1
      utotm(i,j)=(u(i,j,km)+ubavg(i,j,m))*scuv(i,j)
      if (utotm(i,j).ge.0.) then
        q=min(dp(i-1,j,kn),max(0.,depthu(i,j)-util3(i-1,j)))
      else
        q=min(dp(i  ,j,kn),max(0.,depthu(i,j)-util3(i  ,j)))
      end if
      uflux(i,j)=utotm(i,j)*q
      uflux4(i,j)=ubavg(i,j,m)*dpu(i,j,km)*scuv(i,j) 
      uflux3(i,j)=u(i,j,km)*dpu(i,j,km) 
 11   uflux2(i,j)=utotm(i,j)*dpu(i,j,km)-uflux(i,j)
c
      do 12 j=2,jj1
      do 12 i=1,ii1
      vtotm(i,j)=(v(i,j,km)+vbavg(i,j,m))*scuv(i,j)
      if (vtotm(i,j).ge.0.) then
        q=min(dp(i,j-1,kn),max(0.,depthv(i,j)-util3(i,j-1)))
      else
        q=min(dp(i,j  ,kn),max(0.,depthv(i,j)-util3(i,j  )))
      end if
      vflux(i,j)=vtotm(i,j)*q
      vflux4(i,j)=vbavg(i,j,m)*dpv(i,j,km)*scuv(i,j) 
      vflux3(i,j)=v(i,j,km)*dpv(i,j,km) 
 12   vflux2(i,j)=vtotm(i,j)*dpv(i,j,km)-vflux(i,j)
c
      do 32 i=1,ii1
      vtotm(i,1)=(v(i,1,km)+vbavg(i,1,m))*scuv(i,1)
        q=min(dp(i,1  ,kn),max(0.,depthv(i,1)-util3(i,1  )))
      if (vtotm(i,j).ge.0.) then
        q=min(2.*dp(i,1,kn)-dp(i,2,kn),max(0.,depthv(i,1)-
     .(2.*util3(i,1)-util3(i,2))))
      else
        q=min(dp(i,1  ,kn),max(0.,depthv(i,1)-util3(i,1  )))
      endif
      vflux(i,1)=vtotm(i,1)*q
      vflux3(i,1)=v(i,1,km)*dpv(i,1,km) 
      vflux4(i,1)=vbavg(i,1,m)*dpv(i,1,km)*scuv(i,1) 
 32   vflux2(i,1)=vtotm(i,1)*dpv(i,1,km)-vflux(i,1)
c
c --- advance -dp- field using low-order (diffusive) flux values
c
      do 19 j=1,jj1
      do 19 i=1,ii1
      dpold(i,j,k)=dp(i,j,kn)
      util3(i,j)=util3(i,j)+dp(i,j,kn)
 19   dp(i,j,kn)=dp(i,j,kn)-(uflux(i+1,j)-uflux(i,j)
     .                      +vflux(i,j+1)-vflux(i,j))*delt1*sc2i(i,j)
c
c --- at each grid point, determine the ratio of the largest permissible
c --- pos. (neg.) change in -dp- to the sum of all incoming (outgoing) fluxes
c
      do 26 j=1,jj1
      ja=max(1,j-1)
      jb=min(jj1,j+1)
      do 26 i=1,ii1
      ia=max(1,i-1)
      ib=min(ii1,i+1)
      util1(i,j)=max(dp(i,j,kn),dp(ia,j,kn),dp(ib,j,kn),
     .                          dp(i,ja,kn),dp(i,jb,kn))
      util2(i,j)=max(0.,
     .           min(dp(i,j,kn),dp(ia,j,kn),dp(ib,j,kn),
     .                          dp(i,ja,kn),dp(i,jb,kn)))
c
      util1(i,j)=(util1(i,j)-dp(i,j,kn))
     ./((max(0.,uflux2(i,j))-min(0.,uflux2(i+1,j))
     .  +max(0.,vflux2(i,j))-min(0.,vflux2(i,j+1))+epsil)
     .*delt1*sc2i(i,j))
c
 26   util2(i,j)=(util2(i,j)-dp(i,j,kn))
     ./((min(0.,uflux2(i,j))-max(0.,uflux2(i+1,j))
     .  +min(0.,vflux2(i,j))-max(0.,vflux2(i,j+1))-epsil)
     .*delt1*sc2i(i,j))
c
c --- limit antidiffusive fluxes
c --- (keep track in -utotn,vtotn- of discrepancy between high-order
c --- fluxes and the sum of low-order and clipped antidiffusive fluxes.
c --- this will be used later to restore nondivergence of barotropic flow)
c
      do 28 j=1,jj1
      do 28 i=2,ii1
      if (uflux2(i,j).ge.0.) then
        clip=min(1.,util1(i,j),util2(i-1,j))
      else
        clip=min(1.,util2(i,j),util1(i-1,j))
      end if
      utotn(i,j)=utotn(i,j)+uflux2(i,j)*(1.-clip)
 28   uflux(i,j)=uflux2(i,j)*clip
c
      do 49 i=1,ii-1
      if (vflux2(i,1).ge.0.) then
        clip=min(1.,util1(i,1))
      else
        clip=min(1.,util2(i,1))
      end if
      vtotn(i,1)=vtotn(i,1)+vflux2(i,1)*(1.-clip)
 49   vflux(i,1)=vflux2(i,1)*clip
c
      do 29 j=2,jj-1
      do 29 i=1,ii-1
      if (vflux2(i,j).ge.0.) then
        clip=min(1.,util1(i,j),util2(i,j-1))
      else
        clip=min(1.,util2(i,j),util1(i,j-1))
      end if
      vtotn(i,j)=vtotn(i,j)+vflux2(i,j)*(1.-clip)
 29   vflux(i,j)=vflux2(i,j)*clip
c
c --- evaluate effect of antidiffusive fluxes on -dp- field
c
      do 15 j=1,jj1
      do 15  i=1,ii1
      dp(i,j,kn)=dp(i,j,kn)-(uflux(i+1,j)-uflux(i,j)
     .                      +vflux(i,j+1)-vflux(i,j))*delt1*sc2i(i,j)
      delp5(i,j,kn)=delp5(i,j,kn)-(uflux4(i+1,j)-uflux4(i,j)
     .                      +vflux4(i,j+1)-vflux4(i,j))*delt1*sc2i(i,j)
 15   p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
c
 76   continue
c
      call boundp(p,2)
      call boundp(dp,nn+1)
      call boundp(delp5,nn+1)
      call boundp(dpold,1)

c --- restore nondivergence of vertically integrated mass flow by
c --- recovering fluxes lost in the flux limiting process.
c --- treat these fluxes as an 'upstream' barotropic correction to
c --- the sum of diffusive and antidiffusive fluxes obtained so far.
c
      do 77 k=1,kk
      km=k+mm
      kn=k+nn
c
      do 44 j=1,jj1
      do 44 i=2,ii1
      if (utotn(i,j).ge.0.) then
        q=dp(i-1,j,kn)/p(i-1,j,kk+1)
      else
        q=dp(i  ,j,kn)/p(i  ,j,kk+1)
      end if
 44   uflux(i,j)=utotn(i,j)*q
c
      do 45 j=2,jj1
      do 45 i=1,ii1
      if (vtotn(i,j).ge.0.) then
        q=dp(i,j-1,kn)/p(i,j-1,kk+1)
      else
        q=dp(i,j  ,kn)/p(i,j  ,kk+1)
      end if
 45   vflux(i,j)=vtotn(i,j)*q
c
      do 55 i=1,ii1
      if (vtotn(i,1).ge.0.) then
      q=(2.*dp(i,1,kn)-dp(i,2,kn))/(2.*p(i,1,kk+1)-p(i,2,kk+1))
      else
        q=dp(i,1  ,kn)/p(i,1  ,kk+1)
      endif
 55   vflux(i,1)=vtotn(i,1)*q
c
      do 14 j=1,jj1
      do 14 i=1,ii1
      dp(i,j,kn)=dp(i,j,kn)-(uflux(i+1,j)-uflux(i,j)
     .                      +vflux(i,j+1)-vflux(i,j))*delt1*sc2i(i,j)
 14   p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)

c
 77   continue
      call boundp(p,2)
      call boundp(dp,nn+1)

c --- add bottom-pressure restoring term arising from split-explicit treatment
c --- of continuity equation (step 4 in appendix b to 1992 brhs paper)
c
      do 39 j=1,jj1
      do 39 k=1,kk
      kn=k+nn
      do 39 i=1,ii1
      delp5(i,j,kn)=delp5(i,j,kn)+dp(i,j,kn)*(pbot(i,j)/p(i,j,kk+1)-1.)
      dp(i,j,kn)=dp(i,j,kn)*pbot(i,j)/p(i,j,kk+1)
 39   p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)

      call boundp(p,2)
      call boundp(dp,nn+1)
      call boundp(delp5,nn+1)
c
c
c --- ---------------------------------------------------------------------
c --- biharmonic thickness diffusion (literally, interface depth diffusion)
c --- ---------------------------------------------------------------------
c
c     if (thkdff.eq.0.) return
c
      do 23 j=1,jj
c
      do 22 i=1,ii
 22   uflux(i,j)=0.
c
      do 23 i=1,ii
 23   vflux(i,j)=0.
c
c
      do 13 k=2,kk
c
      do 141 j=1,jj1
c --- limit fluxes to avoid intertwining interfaces
      do 141 i=2,ii1
      flxhi= .25*(p(i  ,j,kk+1)-p(i  ,j,k))*scu2(i  ,j)
      flxlo=-.25*(p(i-1,j,kk+1)-p(i-1,j,k))*scu2(i-1,j)
c
 141  uflux(i,j)=min(flxhi,max(flxlo,
     .   delt1*thkdff*(p(i-1,j,k)-p(i,j,k))*scuv(i,j)))
c
      do 151 i=1,ii1
c --- limit fluxes to avoid intertwining interfaces
      do 151 j=2,jj1
      flxhi= .25*(p(i,j  ,kk+1)-p(i,j  ,k))*scu2(i,j  )
      flxlo=-.25*(p(i,j-1,kk+1)-p(i,j-1,k))*scu2(i,j-1)
c
 151  vflux(i,j)=min(flxhi,max(flxlo,
     .   delt1*thkdff*(p(i,j-1,k)-p(i,j,k))*scuv(i,j)))
c
      do 18 j=1,jj1
      do 18 i=1,ii1
 18   p(i,j,k)=p(i,j,k)-(uflux(i+1,j)-uflux(i,j)
     .                  +vflux(i,j+1)-vflux(i,j))*sc2i(i,j)
c
 13   continue
      call boundp(p,2)

      do k=1,kk
      kn=k+nn
      km=k+mm
      do i=1,ii
      do j=1,jj
      p1(i,j,k+1)=p1(i,j,k)+dp(i,j,km)
      p0(i,j,k+1)=p0(i,j,k)+dpold(i,j,k)
      enddo
      enddo
      enddo
      call boundpne1(p,p1,p0)

      do k=1,kk
      kn=k+nn
      do i=1,ii
      do j=1,jj 
      dp(i,j,kn)=p(i,j,k+1)-p(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
