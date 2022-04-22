      subroutine barotp(itime,n)
c
c --- version 2.8
      implicit none
c
      include 'dimensions.h'
      include 'common_blocks.h'
c
      real utndcy(idm,jdm),vtndcy(idm,jdm),umax,uabs,hulp(idm)
      real botstrub(idm,jdm),botstrvb(idm,jdm),pbo(idm,jdm)
      real pbten(idm,jdm),vbo(idm,jdm),ubo(idm,jdm)
      integer lll,ml,nl,mn,ll,itime,ia,ja
      logical vthenu
c
c --- ------------------------------------------------------------------------
c --- advance barotropic equations from baroclinic time level -m- to level -n-
c --- ------------------------------------------------------------------------
c
      ml=n
      nl=3
      
c --- explicit time integration of barotropic flow (forward-backward scheme)
c --- in order to combine forward-backward scheme with leapfrog treatment of
c --- coriolis term, v-eqn must be solved before u-eqn every other time step
      vthenu=.false.
c
      do 840 lll=1,lstep
c
c --- continuity equation
c
      do 843 j=1,jj1
      do 843 i=1,ii1
       pbo(i,j)=pbavg(i,j,nl)
      pbten(i,j)=((ubavg(i+1,j,ml)*depthu(i+1,j)*scuv(i+1,j)
     .                 -ubavg(i  ,j,ml)*depthu(i  ,j)*scuv(i  ,j))
     .                 +(vbavg(i,j+1,ml)*depthv(i,j+1)*scuv(i,j+1)
     .                 -vbavg(i,j  ,ml)*depthv(i,j  )*scuv(i,j )))
     .  *sc2i(i,j)
 843  pbavg(i,j,nl)=(1.-wbaro)*pbavg(i,j,ml)+wbaro*pbavg(i,j,nl)
     . -(1.+wbaro)*dlt*pbten(i,j)
c
      do i=1,ii1
c     hulp(i)=pbavg(i,1,nl)
      enddo

c     call boundpb(itime,pbavg,nl,lll)

      mn=ml
      if (vthenu) go to 901
c
c --- u momentum equation
c
 900  continue
c
      do 841 j=1,jj1
      do 841 i=2,ii1
      ubo(i,j)=ubavg(i,j,nl)
      utndcy(i,j)=-thref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scui(i,j)
     .+(vbavg(i  ,j,mn)*depthv(i  ,j)+vbavg(i  ,j+1,mn)*depthv(i  ,j+1)
     . +vbavg(i-1,j,mn)*depthv(i-1,j)+vbavg(i-1,j+1,mn)*depthv(i-1,j+1))
     . *(pvtrop(i,j)+pvtrop(i,j+1))*.125
      
c

c --- bottom boundary layer stress. stress profile is assumed linear
      botstrub(i,j)=-ubavg(i,j,nl)*.5*(dragb(i,j)+dragb(i-1,j))*
     .     thkbot/depthu(i,j)

      ubavg(i,j,nl)=(1.-wbaro)*ubavg(i,j,ml)+wbaro*ubavg(i,j,nl)
     . +(1.+wbaro)*dlt*(utndcy(i,j)+utotn(i,j)+botstrub(i,j))
c
 841  continue
      call boundub(ubavg,utndcy,botstrub,utotn,nl)
c
      mn=nl
      if (vthenu) go to 902
c
c --- v momentum equation
c
 901  continue
c
      do 842 i=1,ii1
      do 842 j=2,jj1
      vbo(i,j)=vbavg(i,j,nl)
      vtndcy(i,j)=-thref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scui(i,j)
     .-(ubavg(i,j  ,mn)*depthu(i,j  )+ubavg(i+1,j  ,mn)*depthu(i+1,j  )
     . +ubavg(i,j-1,mn)*depthu(i,j-1)+ubavg(i+1,j-1,mn)*depthu(i+1,j-1))
     . *(pvtrop(i,j)+pvtrop(i+1,j))*.125
       
      botstrvb(i,j)=-vbavg(i,j,nl)*.5*(dragb(i,j)+dragb(i,j-1))*
     .     thkbot/depthv(i,j)

      vbavg(i,j,nl)=(1.-wbaro)*vbavg(i,j,ml)+wbaro*vbavg(i,j,nl)
     . +(1.+wbaro)*dlt*(vtndcy(i,j)+vtotn(i,j)+botstrvb(i,j))
c
 842  continue
c
c     call boundvb(vbavg,vtndcy,botstrvb,vtotn,nl)
      call boundvb(itime,vbavg,nl,lll)

      mn=nl
      if (vthenu) go to 900
c
c --- switch order in which -u,v- equations are solved
 902  vthenu=.not.vthenu
c
      ll=ml
      ml=nl
      nl=ll
c
 840  continue

      return
      end
