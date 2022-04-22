      subroutine momeq2(itime,m,n,mm,nn,k1n,k1m)
c
c --- version 2.8.1 (biharmonic)
      implicit none
c
      include 'dimensions.h'
      include 'common_blocks.h'
      include 'stmt_funcs.h'
c
      real botstr(idm,jdm),visc(idm,jdm),
     .     wgtia(idm,jdm),wgtib(idm,jdm),wgtja(idm,jdm),wgtjb(idm,jdm),
     .     dpxy,dpia,dpib,dpja,dpjb,visca,viscb,ptopl,pbotl,cutoff,q,
     .     dt1inv,vibb,verm,uold(idm,jdm,kdm),vold(idm,jdm,kdm),ujbb
      real dl2u(idm,jdm),dl2uja(idm,jdm),dl2ujb(idm,jdm),
     .     dl2v(idm,jdm),dl2via(idm,jdm),dl2vib(idm,jdm)
      real pu(idm,jdm,kdm+1),pv(idm,jdm,kdm+1)
      real defor1(idm,jdm),defor2(idm,jdm)
      real uflux1(idm,jdm),vflux1(idm,jdm)
      real uflux3(idm,jdm),vflux3(idm,jdm)
      real pgfx(idm,jdm),pgfy(idm,jdm)
      real uja(idm,jdm),ujb(idm,jdm)
      real via(idm,jdm),vib(idm,jdm)

      integer kan,jc,ia,ja,ih,itime
c
      cutoff=onem
      ih=ii/2
c     goto 832
c
      do 804 j=1,jj
      do 804 i=1,ii
c --- spatial weighting function for pressure gradient calculation:
      pu(i,j,1)=0.
      pv(i,j,1)=0.
      defor1(i,j)=0.
      defor2(i,j)=0.
        via(i,j)=0.
        vib(i,j)=0.
        uja(i,j)=0.
        ujb(i,j)=0.
       pgfx(i,j)=0.
       pgfy(i,j)=0.
       dl2v(i,j)=0.
       dl2u(i,j)=0.
       dl2via(i,j)=0.
       dl2vib(i,j)=0.
       dl2uja(i,j)=0.
       dl2ujb(i,j)=0.
      util1(i,j)=0.
 804  util2(i,j)=0.
c
      do j=1,jj
      vpfl1(j)=0.
      vpfl5(j)=0.
      enddo

      do 9 k=1,kk
      km=k+mm
      kn=k+nn
c
c --- store total (barotropic plus baroclinic) flow at old and mid time in
c --- -utotn,vtotn- and -utotm,vtotm- respectively. 
c
      do 807 j=1,jj
      do 807 i=1,ii
      utotm(i,j)=u(i,j,km)+ubavg(i,j,m)
      utotn(i,j)=u(i,j,kn)+ubavg(i,j,n)
      uflux(i,j)=utotm(i,j)*max(dpu(i,j,km),cutoff)
 807  pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,km)
c
      do 808 j=1,jj
      do 808 i=1,ii
      vtotm(i,j)=v(i,j,km)+vbavg(i,j,m)
      vtotn(i,j)=v(i,j,kn)+vbavg(i,j,n)
      vflux(i,j)=vtotm(i,j)*max(dpv(i,j,km),cutoff)
 808  pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,km)
c
      do i=1,ii1
      do j=1,jj1
      vpfl1(j)=vpfl1(j)+0.5*(v(i,j,km)*dpv(i,j,km)+
     .v(i,j+1,km)*dpv(i,j+1,km))*pre(i,j,k)
      vpfl5(j)=vpfl5(j)+0.5*(v(i,j,km)*dpv(i,j,km)+
     .v(i,j+1,km)*dpv(i,j+1,km))*pr5(i,j,k)
      enddo
      enddo

c --- define auxiliary velocity fields (via,vib,uja,ujb) to implement
c --- sidewall friction along near-vertical bottom slopes. wgtja,wgtjb,wgtia,
c --- wgtib indicate the extent to which a sidewall is present.
c
      do 885 j=1,jj1
      do 885 i=1,ii
      wgtjb(i,j)=max(0.,min(1.,(pu(i,j,k+1)-depthu(i,j+1))
     .          /max(pu(i,j,k+1)-pu(i,j,k),epsil)))
 885  ujb(i,j)=(1.-wgtjb(i,j))*utotn(i,j+1)+wgtjb(i,j)*slip*utotn(i,j)
      continue
c
      do 805 j=2,jj
      do 805 i=1,ii
      wgtja(i,j)=max(0.,min(1.,(pu(i,j,k+1)-depthu(i,j-1))
     .          /max(pu(i,j,k+1)-pu(i,j,k),epsil)))
 805  uja(i,j)=(1.-wgtja(i,j))*utotn(i,j-1)+wgtja(i,j)*slip*utotn(i,j)
c
      do i=1,ii
      wgtja(i,1)=max(0.,min(1.,(pu(i,1,k+1)-depthu(i,1))
     .          /max(pu(i,1,k+1)-pu(i,1,k),epsil)))
      uja(i,1)=(1.-wgtja(i,1))*utotn(i,1)+wgtja(i,1)*slip*utotn(i,1)
      enddo
c
      do j=1,jj
      do i=1,ii
      dl2u(i,j)=utotn(i,j)
      enddo
      enddo
c
      do 886 i=1,ii1
      do 886 j=1,jj
      wgtib(i,j)=max(0.,min(1.,(pv(i,j,k+1)-depthv(i+1,j))
     .          /max(pv(i,j,k+1)-pv(i,j,k),epsil)))
 886  vib(i,j)=(1.-wgtib(i,j))*vtotn(i+1,j)+wgtib(i,j)*slip*vtotn(i,j)

      do 806 i=2,ii
      do 806 j=1,jj
      wgtia(i,j)=max(0.,min(1.,(pv(i,j,k+1)-depthv(i-1,j))
     .          /max(pv(i,j,k+1)-pv(i,j,k),epsil)))
 806  via(i,j)=(1.-wgtia(i,j))*vtotn(i-1,j)+wgtia(i,j)*slip*vtotn(i,j)

      do j=1,jj
      wgtia(1,j)=max(0.,min(1.,(pv(1,j,k+1))
     .          /max(pv(1,j,k+1)-pv(1,j,k),epsil)))
      via(1,j)=wgtia(1,j)*slip*vtotn(1,j)
      enddo

      do i=1,ii
      do j=1,jj
      dl2v(i,j)=vtotn(i,j)
       enddo
       enddo

      do 63 j=1,jj1
      do 63 i=1,ii1
 63   defor1(i,j)=((utotn(i+1,j)*scuv(i+1,j)-utotn(i,j)*scuv(i,j))
     .            -(vtotn(i,j+1)*scuv(i,j+1)-vtotn(i,j)*scuv(i,j)))**2
     .            *sc2i(i,j)
c
      do 64 j=2,jj
      do 64 i=2,ii
 64   defor2(i,j)=(vib(i-1,j)*scuv(i,j)-via(i,j)*scuv(i-1,j)
     .            +ujb(i,j-1)*scuv(i,j)-uja(i,j)*scuv(i,j-1))**2
     .            *sc2i(i,j)
c
      do i=2,ii
      ujbb=utotn(i,1)
      defor2(i,1)=(vib(i-1,1)*scuv(i,1)-via(i,1)*scuv(i-1,1)
     .            +ujbb*scuv(i,1)-uja(i,1)*scuv(i,1))**2
     .            *sc2i(i,1)
      enddo
c
      do j=2,jj
      vibb=vtotn(1,j)
      defor2(1,j)=(vibb*scuv(1,j)-via(1,j)*scuv(1,j)
     .            +ujb(1,j-1)*scuv(1,j)-uja(1,j)*scuv(1,j-1))**2
     .            *sc2i(1,j)
      enddo
c
      vibb=vtotn(1,1)
      ujbb=utotn(1,1)
      defor2(1,1)=(vibb*scuv(1,1)-via(1,1)*scuv(1,1)
     .            +ujbb*scuv(1,1)-uja(1,1)*scuv(1,1))**2
     .            *sc2i(1,1)
c
c --- define auxiliary del2 fields (dl2via,dl2vib,dl2uja,dl2ujb) to imple-
c --- ment biharmonic sidewall friction along near-vertical bottom slopes.
c
      do 905 j=2,jj1
      do 905 i=2,ii1
 905  dl2uja(i,j)=(1.-wgtja(i,j))*dl2u(i,j-1)+wgtja(i,j)*slip*dl2u(i,j)
c
      do i=2,ii1
      dl2uja(i,1)=(1.-wgtja(i,1))*dl2u(i,1)+wgtja(i,1)*slip*dl2u(i,1)
      enddo
c
      do 705 j=1,jj1
      do 705 i=2,ii1
 705  dl2ujb(i,j)=(1.-wgtjb(i,j))*dl2u(i,j+1)+wgtjb(i,j)*slip*dl2u(i,j)
c
      do 906 i=2,ii1
      do 906 j=2,jj1
 906  dl2via(i,j)=(1.-wgtia(i,j))*dl2v(i-1,j)+wgtia(i,j)*slip*dl2v(i,j)
c
      do j=2,jj1
      dl2via(1,j)=wgtia(1,j)*slip*dl2v(1,j)
      enddo
c
      do 706 i=1,ii1
      do 706 j=2,jj1
 706  dl2vib(i,j)=(1.-wgtib(i,j))*dl2v(i+1,j)+wgtib(i,j)*slip*dl2v(i,j)
c
c --- ----------
c --- u equation
c --- ----------
c
c --- deformation-dependent eddy viscosity coefficient
c
      do 37 j=1,jj1
      do 37 i=2,ii1
 37   visc(i,j)=max(veldff,viscos*
     .sqrt(.5*(defor1(i,j)+defor1(i-1,j)+defor2(i,j)+defor2(i,j+1))))
c
      do 830 j=1,jj1
      do 824 i=2,ii-2
      verm=sqrt(800./min(800.,max(depthu(i,j),depthu(i+1,j))
     ./onem))
 824  uflux1(i,j)=(visc(i,j)+visc(i+1,j))*(dl2u(i,j)-dl2u(i+1,j))
     .             *hfharm(max(dpu(i  ,j,km),onemm),
     .                     max(dpu(i+1,j,km),onemm))
     .       *verm *scu2(i,j)*2./(scuv(i,j)+scuv(i+1,j))
c
      verm=sqrt(800./min(800.,max(depthu(1,j),depthu(2,j))
     ./onem))
      uflux1(1,j)=(visc(2,j)+visc(2,j))*(dl2u(1,j)-dl2u(2,j))
     .             *hfharm(max(dpu(1,j,km),onemm),
     .                     max(dpu(2,j,km),onemm))
     .     *verm   *scu2(1,j)*2./(scuv(1,j)+scuv(2,j))
c
      verm=sqrt(800./min(800.,max(depthu(ii1,j),depthu(ii,j))
     ./onem))
      uflux1(ii1,j)=(visc(ii1,j)+visc(ii1,j))*(dl2u(ii1,j)-dl2u(ii,j))
     .             *hfharm(max(dpu(ii1,j,km),onemm),
     .                     max(dpu(ii,j,km),onemm))
     .     *verm   *scu2(ii1,j)*2./(scuv(ii1,j)+scuv(ii,j))
c
c --- lateral turb. momentum flux (at vorticity points)
c --- (left and right fluxes are evaluated separately because of sidewalls)
 830  continue
c
      do 840 i=2,ii1
      do 820 j=2,jj1
      dpxy=max(dpu(i,j ,km),onemm)
      dpja=max(dpu(i,j-1,km),onemm)
      verm=sqrt(800./min(800.,max(depthu(i,j),depthu(i,j-1))
     ./onem))
c
        visca=visc(i,j-1)
 820  uflux2(i,j)=(visc(i,j)+visca)*(dl2uja(i,j)-dl2u(i,j))
     .            *hfharm(dpja+wgtja(i,j)*(dpxy-dpja),dpxy)
     .*verm       *scu2(i,j )*2./(scuv(i,j)+scuv(i,j-1))
      do 822 j=1,jj-2
      dpxy=max(dpu(i,j ,km),onemm)
      dpjb=max(dpu(i,j+1,km),onemm)
      verm=sqrt(800./min(800.,max(depthu(i,j),depthu(i,j+1))
     ./onem))
c
        viscb=visc(i,j+1)
 822  uflux3(i,j)=(visc(i,j)+viscb)*(dl2u(i,j)-dl2ujb(i,j))
     .            *hfharm(dpjb+wgtjb(i,j)*(dpxy-dpjb),dpxy)
     . *verm      *scu2(i,j+1)*2./(scuv(i,j)+scuv(i,j+1))
      dpxy=max(dpu(i,jj1 ,km),onemm)
      dpjb=max(dpu(i,jj,km),onemm)
      verm=sqrt(800./min(800.,max(depthu(i,jj),depthu(i,jj1))
     ./onem))
      uflux3(i,jj1)=(visc(i,jj1)+visc(i,jj1))*(dl2u(i,jj1)-
     .dl2ujb(i,jj1))
     .            *hfharm(dpjb+wgtjb(i,jj1)*(dpxy-dpjb),dpxy)
     . *verm*scu2(i,jj)*2./(scuv(i,jj1)+scuv(i,jj))
      verm=sqrt(800./min(800.,depthu(i,1)
     ./onem))
      dpxy=max(dpu(i,1 ,km),onemm)
      dpja=max(dpu(i,1,km),onemm)
      uflux2(i,1)=(visc(i,1)+visc(i,1))*(dl2uja(i,1)-dl2u(i,1))
     .            *hfharm(dpja+wgtja(i,1)*(dpxy-dpja),dpxy)
     .*verm      *scu2(i,1 )*2./(scuv(i,1)+scuv(i,1))
 840  continue
c
      if (k.gt.1) then
c
c --- pressure force in x direction
c --- ('scheme 2' from appendix -a- in bleck-smith paper)
c
      do 96 j=1,jj1
      do 96 i=2,ii1
      util1(i,j)=max(0.,min(depthu(i,j)-pu(i,j,k),
     .           pu(i,j,k+1)-pu(i,j,2),h1))
 96   pgfx(i,j)=(montg(i,j,k)-montg(i-1,j,k))
     .  *util1(i,j)
c
      do 98 j=2,jj1
      do 98 i=2,ii1
c
 98   gradx(i,j)=(pgfx(i,j)
     .+(h1-util1(i,j))*
     .  (pgfx (i-1,j)+pgfx (i+1,j)+pgfx (i,j-1)+pgfx (i,j+1))/
     .  (util1(i-1,j)+util1(i+1,j)+util1(i,j-1)+util1(i,j+1)+epsil))/h1
 
      do 88 i=2,ii1
 
 88   gradx(i,1)=(pgfx(i,1)
     .+(h1-util1(i,1))*
     .  (pgfx (i-1,1)+pgfx (i+1,1)+pgfx (i,2))/
     .  (util1(i-1,1)+util1(i+1,1)+util1(i,2)+epsil))/h1
c
      end if				!  k > 1
c
      do 16 j=1,jj1
      do 16 i=2,ii1
c
      ptopl=min(depthu(i,j),.5*(p(i,j,k  )+p(i-1,j,k  )))
      pbotl=min(depthu(i,j),.5*(p(i,j,k+1)+p(i-1,j,k+1)))
c
      botstr(i,j)=0.
      if(bdrag) then
c --- bottom boundary layer stress. stress profile is assumed linear
      botstr(i,j)=-utotn(i,j)*.5*(drag(i,j)+drag(i-1,j))*
     .     (max(depthu(i,j)-thkbot,          pbotl       )
     .     -max(depthu(i,j)-thkbot,min(ptopl,pbotl-onemm)))
     .     /max(dpu(i,j,km),onemm)
c
      endif
c --- time smoothing of -u- field  (part 1)
      u(i,j,km)=u(i,j,km)*(wuv1*dpu(i,j,km)+onemm)
     .         +u(i,j,kn)* wuv2*dpu(i,j,kn)
c
 16   continue
      do 6 j=1,jj1
      do 6 i=2,ii1

      uold(i,j,k)=u(i,j,kn)
      u(i,j,kn)=u(i,j,kn)+delt1*(-scui(i,j)*(gradx(i,j)
     .+.25*((utotm(i+1,j)**2-utotm(i-1,j)**2)+(vtotm(i,j)**2-
     .vtotm(i-1,j)**2)+(vtotm(i,j+1)**2-vtotm(i-1,j+1)**2)))
     .+.125*(vflux(i  ,j)+vflux(i  ,j+1)+vflux(i-1,j)+vflux(i-1,j+1))
     .     *(potvor(i,j,k)+potvor(i,j+1,k)) - ubrhs(i,j) 
     .+ botstr(i,j)
     .-(uflux1(i,j)-uflux1(i-1,j)
     . +uflux3(i,j)-uflux2(i  ,j))/(scu2(i,j)*max(dpu(i,j,km),onemm)))
c
   6  continue
c     if(k.eq.2.and.itime.ge.72) then
c     print *, u(2,14,kn),u(ii1,14,kn), 'unew'
c     print *, uold(2,14,k),uold(ii1,14,k), 'uold'
c     print *, gradx(2,14),gradx(ii1,14), 'gradx'
c     print *,(utotm(3,14)**2-utotm(1,14)**2)+(vtotm(2,14)**2-
c    .vtotm(1,14)**2)+(vtotm(2,15)**2-vtotm(1,15)**2),
c    .(utotm(ii,14)**2-utotm(ii-2,14)**2)+(vtotm(ii1,14)**2-
c    .vtotm(ii-2,14)**2)+(vtotm(ii1,15)**2-vtotm(ii-2,15)**2),
c    .'nonl'
c     print *,vflux(2,14)+vflux(2,15)+vflux(1,14)+vflux(1,15),
c    .vflux(ii1,14)+vflux(ii1,15)+vflux(ii-2,14)+vflux(ii-2,15),
c    .'vflux'
c     print *,potvor(2,14,k)+potvor(2,15,k),
c    .potvor(ii1,14,k)+potvor(ii1,15,k),'potvor'
c     print *, ubrhs(2,14),ubrhs(ii1,14), 'ubrhs'
c     print *,uflux1(2,14),uflux1(ii1,14), 'uflux1'
c     print *,uflux1(1,14),uflux1(ii-2,14), 'uflux1b'
c     print *,uflux3(2,14),uflux3(ii1,14), 'uflux3'
c     print *,uflux2(2,14),uflux2(ii1,14), 'uflux2'
c     endif
c --- ----------
c --- v equation
c --- ----------
c
c --- deformation-dependent eddy viscosity coefficient
c
      do 38 j=2,jj1
      do 38 i=1,ii1
 38   visc(i,j)=max(veldff,viscos*
     .sqrt(.5*(defor1(i,j)+defor1(i,j-1)+defor2(i,j)+defor2(i+1,j))))
c
      do i=1,ii1
      visc(i,1)=max(veldff,viscos*
     .sqrt(.5*(defor1(i,1)+defor2(i,1)+defor2(i+1,1))))
      enddo
c
      do 826 i=1,ii1
      do 826 j=1,jj-2
      verm=sqrt(800./min(800.,max(depthv(i,j),depthv(i,j+1))
     ./onem))
            vflux1(i,j)=(visc(i,j)+visc(i,j+1))*(dl2v(i,j)-dl2v(i,j+1))
     .             *hfharm(max(dpv(i,j  ,km),onemm),
     .                     max(dpv(i,j+1,km),onemm))
     .  *verm      *scu2(i,j)*2./(scuv(i,j)+scuv(i,J+1))
 826  continue
c
      do i=1,ii1
      verm=sqrt(800./min(800.,max(depthv(i,jj),depthv(i,jj1))
     ./onem))
      vflux1(i,jj1)=(visc(i,jj1)+visc(i,jj1))*(dl2v(i,jj1)-dl2v(i,jj))
     .             *hfharm(max(dpv(i,jj1  ,km),onemm),
     .                     max(dpv(i,jj,km),onemm))
     .   *verm     *scu2(i,jj1)*2./(scuv(i,jj1)+scuv(i,jj))
      enddo   
c --- lateral turb. momentum flux (at vorticity points)
c --- (left and right fluxes are evaluated separately because of sidewalls)
c
      do 821 i=2,ii1
      do 821 j=2,jj1
      dpxy=max(dpv(i ,j,km),onemm)
      dpia=max(dpv(i-1,j,km),onemm)
      dpib=max(dpv(i+1,j,km),onemm)
c
        visca=visc(i-1,j)
        viscb=visc(i+1,j)

      verm=sqrt(800./min(800.,max(depthv(i,j),depthv(i-1,j))
     ./onem))
      vflux2(i,j)=(visc(i,j)+visca)*(dl2via(i,j)-dl2v(i,j))
     .            *hfharm(dpia+wgtia(i,j)*(dpxy-dpia),dpxy)
     .  *verm     *scu2(i ,j)*2./(scuv(i,j)+scuv(i-1,j))
      verm=sqrt(800./min(800.,max(depthv(i,j),depthv(i+1,j))
     ./onem))
 821  vflux3(i,j)=(visc(i,j)+viscb)*(dl2v(i,j)-dl2vib(i,j))
     .            *hfharm(dpib+wgtib(i,j)*(dpxy-dpib),dpxy)
     . *verm      *scu2(i+1,j)*2./(scuv(i,j)+scuv(i+1,j))
c
      do 921 j=2,jj1
      dpxy=max(dpv(1 ,j,km),onemm)
      dpib=max(dpv(2,j,km),onemm)
c
        viscb=visc(2,j)
      dpia=onemm
      verm=sqrt(800./min(800.,max(depthv(1,j),depthv(1,j))
     ./onem))

      vflux2(1,j)=(visc(1,j)+visc(1,j))*(dl2via(1,j)-dl2v(1,j))
     .            *hfharm(dpia+wgtia(1,j)*(dpxy-dpia),dpxy)
     .*verm       *scu2(1 ,j)*2./(scuv(1,j)+scuv(1,j))
      verm=sqrt(800./min(800.,max(depthv(1,j),depthv(2,j))
     ./onem))
 921  vflux3(1,j)=(visc(1,j)+viscb)*(dl2v(1,j)-dl2vib(1,j))
     .            *hfharm(dpib+wgtib(1,j)*(dpxy-dpib),dpxy)
     . *verm      *scu2(2,j)*2./(scuv(1,j)+scuv(2,j))
c
      do j=2,jj1
      dpxy=max(dpv(ii1,j,km),onemm)
      dpib=onemm
        viscb=visc(ii1,j)

      verm=sqrt(800./min(800.,max(depthv(ii1,j),depthv(ii,j))
     ./onem))
      vflux3(ii1,j)=(visc(ii1,j)+viscb)*(dl2v(ii1,j)-dl2vib(ii1,j))
     .            *hfharm(dpib+wgtib(ii1,j)*(dpxy-dpib),dpxy)
     .  *verm     *scu2(ii,j)*2./(scuv(ii1,j)+scuv(ii,j))
      enddo
c
      if (k.gt.1) then
c
c --- pressure force in y direction
c --- ('scheme 2' from appendix -a- in bleck-smith paper)
c
      do 97 j=2,jj1
      do 97 i=1,ii1
      util2(i,j)=max(0.,min(depthv(i,j)-pv(i,j,k),
     .           pv(i,j,k+1)-pv(i,j,2),h1))
 97   pgfy(i,j)=(montg(i,j,k)-montg(i,j-1,k))
     .  *util2(i,j)
c

      do 89 j=2,jj1
c
 89   grady(1,j)=(pgfy(1,j)
     .+(h1-util2(1,j))*
     .  (pgfy (2,j)+pgfy (1,j-1)+pgfy (1,j+1))/
     .  (util2(2,j)+util2(1,j-1)+util2(1,j+1)+epsil))/h1
c
      do 99 i=2,ii1
      do 99 j=2,jj1
c
 99   grady(i,j)=(pgfy(i,j)
     .+(h1-util2(i,j))*
     .  (pgfy (i-1,j)+pgfy (i+1,j)+pgfy (i,j-1)+pgfy (i,j+1))/
     .  (util2(i-1,j)+util2(i+1,j)+util2(i,j-1)+util2(i,j+1)+epsil))/h1
c
      end if				!  k > 1
c
      do 17 j=2,jj1
      do 17 i=1,ii1
c
      ptopl=min(depthv(i,j),.5*(p(i,j,k  )+p(i,j-1,k  )))
      pbotl=min(depthv(i,j),.5*(p(i,j,k+1)+p(i,j-1,k+1)))
c
c --- bottom boundary layer stress. stress profile is assumed linear
      botstr(i,j)=0.
      if(bdrag) then
      botstr(i,j)=-vtotn(i,j)*.5*(drag(i,j)+drag(i,j-1))*
     .     (max(depthv(i,j)-thkbot,          pbotl       )
     .     -max(depthv(i,j)-thkbot,min(ptopl,pbotl-onemm)))
     .     /max(dpv(i,j,km),onemm)
c
      endif
c --- time smoothing of -v- field  (part 1)
      v(i,j,km)=v(i,j,km)*(wuv1*dpv(i,j,km)+onemm)
     .         +v(i,j,kn)* wuv2*dpv(i,j,kn)
c
 17   continue
      do 7 j=2,jj1
      do 7 i=1,ii1
      vold(i,j,k)=v(i,j,kn)
 7    v(i,j,kn)=v(i,j,kn)+delt1*(-scui(i,j)*(grady(i,j)
     .+.25*((vtotm(i,j+1)**2-vtotm(i,j-1)**2)+(utotm(i,j)**2-
     .utotm(i,j-1)**2)+(utotm(i+1,j)**2-utotm(i+1,j-1)**2)))
     .-.125*(uflux(i,j  )+uflux(i+1,j  )+uflux(i,j-1)+uflux(i+1,j-1))
     .     *(potvor(i,j,k)+potvor(i+1,j,k)) - vbrhs(i,j) 
     .+ botstr(i,j)
     .-(vflux1(i,j)-vflux1(i,j-1)
     . +vflux3(i,j)-vflux2(i,j  ))/(scu2(i,j)*max(dpv(i,j,km),onemm)))
c
 9    continue
c
      call boundu(u,uold,nn+1)
      call boundv(v,vold,nn+1)

      dt1inv = 1./delt1
c
      do 14 k=1,kk
      kn=k+nn
      do 14 j=1,jj
c
      do 12 i=1,ii
 12   p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
c
c --- compute new -dpu,dpv- field. save old -dpu,dpv- values in -pu,pv-.
c
      do 13 i=1,ii
 13   pu(i,j,k+1)=dpu(i,j,kn)
c
      do 14 i=1,ii
 14   pv(i,j,k+1)=dpv(i,j,kn)
c
      call dpudpv(p,depthu,depthv,dpu(1,1,k1n),dpv(1,1,k1n))
c
c --- substitute depth-weighted averages for (u,v) at massless grid points.
c --- (scan layers in top-down direction to save time.)
c --- extract barotropic velocities generated during most recent baroclinic
c --- time step and use them to force barotropic flow field.
c
      do j=1,jj
      do i=1,ii
      utotn(i,j)=0.
      vtotn(i,j)=0.
      enddo
      enddo

      do 31 j=1,jj1
      do 32 i=1,ii
 32   utotn(i,j)=0.
      do 33 k=1,kk
      km=k+mm
      kn=k+nn
      kan=max(1,k-1)+nn
      do 33 i=2,ii1
      q=min(dpu(i,j,km),dpu(i,j,kn),onem)
      u(i,j,kn)=(u(i,j,kn)*q+u(i,j,kan)*(onem-q))/onem
      utotn(i,j)=utotn(i,j)+u(i,j,kn)*dpu(i,j,kn)
 33   continue
      do 31 i=2,ii1
 31   utotn(i,j)=utotn(i,j)/depthu(i,j)
c
      do 30 j=2,jj1
      do 34 i=1,ii
 34   vtotn(i,j)=0.
      do 35 k=1,kk
      km=k+mm
      kn=k+nn
      kan=max(1,k-1)+nn
      do 35 i=1,ii1
      q=min(dpv(i,j,km),dpv(i,j,kn),onem)
      v(i,j,kn)=(v(i,j,kn)*q+v(i,j,kan)*(onem-q))/onem
 35   vtotn(i,j)=vtotn(i,j)+v(i,j,kn)*dpv(i,j,kn)
      do 30 i=1,ii1
 30   vtotn(i,j)=vtotn(i,j)/depthv(i,j)

c --- time smoothing of -u,v- fields  (part 2)
c
      do 22 k=1,kk
      km=k+mm
      kn=k+nn
c
      do 24 j=1,jj1
c
      do 24 i=2,ii1
      u(i,j,kn)=u(i,j,kn)-utotn(i,j)
 24   u(i,j,km)=(u(i,j,km)+u(i,j,kn)*wuv2*dpu(i,j,kn))/
     .   (wuv1*dpu(i,j,km)+onemm+wuv2*(pu(i,j,k+1)+dpu(i,j,kn)))
c
      do 22 j=2,jj1
      do 22 i=1,ii1
      v(i,j,kn)=v(i,j,kn)-vtotn(i,j)
 22   v(i,j,km)=(v(i,j,km)+v(i,j,kn)*wuv2*dpv(i,j,kn))/
     .   (wuv1*dpv(i,j,km)+onemm+wuv2*(pv(i,j,k+1)+dpv(i,j,kn)))
c
 832  continue
      kn=kk+nn
      do 867 j=1,jj
c
      do 865 i=1,ii
      utotn(i,j)=utotn(i,j)*dt1inv
 865  ubavg(i,j,n)=ubavg(i,j,m)
c
      do 866 i=1,ii
      vtotn(i,j)=vtotn(i,j)*dt1inv
 866  vbavg(i,j,n)=vbavg(i,j,m)
c
      do 867 i=1,ii
 867  pbavg(i,j,n)=pbavg(i,j,m)
c
      call bounduu(u,nn+1)
      call boundvv(v,nn+1)
      call bounduu(u,mm+1)
      call boundvv(v,mm+1)

      return
      end
