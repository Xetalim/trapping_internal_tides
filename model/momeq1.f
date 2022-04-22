      subroutine momeq1(itime,m,n,mm,nn,k1m)
c
c --- version 2.8
      implicit none
c
      include 'dimensions.h'
      include 'common_blocks.h'
      include 'stmt_funcs.h'
c
      real phi,plo,ubot,vbot,stresx,stresy,pmax,pmin
      integer itime,ih,ia,ja,ib,jb
      real ph(idm,jdm,kdm+1)
c
c
c --- --------------------
c --- hydrostatic equation
c --- --------------------
c
      ih=ii/2
      do j=1,jj
      do i=1,ii
      p(i,j,1)=0.
      ph(i,j,1)=0.
      enddo
      enddo
c
      do 8 j=1,jj
      do 8 k=1,kk
      do 8 i=1,ii
      montg(i,j,k)=0.
 8    p(i,j,k+1)=p(i,j,k)+dp(i,j,k+mm)
c
       do 80 j=1,jj1
      do 80 i=1,ii1
c --- store (1+eta) (= p_total/p_prime) in -util1-
      util1(i,j)=1.+pbavg(i,j,m)/p(i,j,kk+1)
c
c --- save p_total at mxlyr.bottom in -util3-
      util3(i,j)=p(i,j,2)*util1(i,j)
c
c --- m_prime in lowest layer:
      montg(i,j,kk)=psikk(i,j)
     .             -pbavg(i,j,m)*(theta(kk)+thbase)*thref
      pr5(i,j,kk)=0.
 80   pre(i,j,kk)=0.
c
c --- m_prime in remaining layers:
      if(kk.ge.3) then
      do 83 k=kk-1,2,-1
      do 83 i=1,ii1
      do 83 j=1,jj1
      montg(i,j,k)=montg(i,j,k+1)+p(i,j,k+1)*util1(i,j)
     .             *(theta(k+1)-theta(k))*thref
 83   pre(i,j,k)=pre(i,j,k+1)+(p(i,j,k+1)-pin(i,j,k+1))
     .             *(theta(k+1)-theta(k))*thref
      endif

      do 13 j=1,jj1
      do 13 i=1,ii1
      montg(i,j,1)=montg(i,j,2)+p(i,j,2)*util1(i,j)
     .             *(theta(2)-theta(1))*thref
 13   pre(i,j,1)=pre(i,j,2)+(p(i,j,2)-pin(i,j,2))
     .             *(theta(2)-theta(1))*thref
c
      do j=1,jj
      do k=1,kk
      do i=1,ii
      ph(i,j,k+1)=ph(i,j,k)+delp5(i,j,k+mm)
      enddo
      enddo
      enddo
c
      do k=kk-1,1,-1
      do i=1,ii1
      do j=1,jj1
      pr5(i,j,k)=pr5(i,j,k+1)+ph(i,j,k+1)
     .             *(theta(k+1)-theta(k))*thref
      enddo
      enddo
      enddo
c
      call dpudpv(p,depthu,depthv,dpu(1,1,k1m),dpv(1,1,k1m))
c
c +++ ++++++++++++++++++
c +++ momentum equations
c +++ ++++++++++++++++++
c
      do i=1,ii
      do j=1,jj
      drag(i,j)=0.
      enddo
      enddo

      bdrag=.true.
      if(bdrag) then
c --- bottom drag (standard bulk formula)
c
      do 804 j=1,jj1
      do 800 i=1,ii1
      util1(i,j)=0.
 800  util2(i,j)=0.
 
      do 801 k=1,kk
      kn=k+nn
      do 801 i=1,ii1
      phi=max(p(i,j,k+1),p(i,j,kk+1)-thkbot)
      plo=max(p(i,j,k  ),p(i,j,kk+1)-thkbot)
      util1(i,j)=util1(i,j)+(u(i,j,kn)+u(i+1,j,kn))*(phi-plo)
 801  util2(i,j)=util2(i,j)+(v(i,j,kn)+v(i,j+1,kn))*(phi-plo)
c
      do 804 i=1,ii1
      ubot=ubavg(i,j,n)+ubavg(i+1,j,n)+util1(i,j)/thkbot
      vbot=vbavg(i,j,n)+vbavg(i,j+1,n)+util2(i,j)/thkbot
      drag(i,j)=0.003*
     .(.25*sqrt(ubot*ubot+vbot*vbot)+cbar)
     .          *g/(thref*thkbot)
 804  dragb(i,j)=pbot(i,j)/pbot(ih,9)*0.3*exp(-float(itime)/240.)*
     .(.25*sqrt(ubot*ubot+vbot*vbot)+cbar)
     .          *g/(thref*thkbot)
c
      endif
c --- store r.h.s. of barotropic u/v eqn. in -ubrhs,vbrhs-
c --- store layer 1 pressure gradient plus wind forcing in -gradx,grady-
c
      do 69 j=1,jj1
      do 69 i=2,ii1
      ubrhs(i,j)=
     . (vbavg(i  ,j,m)*depthv(i  ,j)+vbavg(i  ,j+1,m)*depthv(i  ,j+1)
     . +vbavg(i-1,j,m)*depthv(i-1,j)+vbavg(i-1,j+1,m)*depthv(i-1,j+1))
     . *(pvtrop(i,j)+pvtrop(i,j+1))*.125
c
 69   gradx(i,j)=montg(i,j,1)-montg(i-1,j,1)
c
      do 70 j=2,jj1
      do 70 i=1,ii1
      vbrhs(i,j)=
     .-(ubavg(i,j  ,m)*depthu(i,j  )+ubavg(i+1,j  ,m)*depthu(i+1,j  )
     . +ubavg(i,j-1,m)*depthu(i,j-1)+ubavg(i+1,j-1,m)*depthu(i+1,j-1))
     . *(pvtrop(i,j)+pvtrop(i+1,j))*.125
c
 70   grady(i,j)=montg(i,j,1)-montg(i,j-1,1)
c
      return
      end
