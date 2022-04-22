      subroutine geopar
c
c --- set up model parameters related to geography
c
c --- version 2.8
      implicit none
c
      include 'dimensions.h'
      include 'common_blocks.h'
      include 'stmt_funcs.h'
c
      real realat
      integer iha
c
c --- express vorticity points in rotated grid in terms of old (lat/lon) grid
c
      iha=(ii+1)/2
c --- define coriolis parameter and grid size
      do 56 i=1,ii
        realat=alat(xpivn-float(i)+.5,gridn)
      do 56 j=1,jj
c     corio(i,j)=0.0*sin(realat)*4.*pi/86400.
      corio(i,j)=1.e-4
c
c --- scux,scuy,scvx,scvy: grid scale at u,v points in x,y dir. respectively
      scuv(i,j)=111.2e5*(alat(xpivn-(iha-.5),gridn)
     .                  -alat(xpivn-(iha+.5),gridn))*radian
c --- size of grid cells (length x width) at u,v,p,q points resp.
      scu2(i,j)=scuv(i,j)*scuv(i,j)
c
      scui(i,j)=1./scuv(i,j)
      sc2i(i,j)=1./scu2(i,j)
   56   continue
c
c --- initialize some arrays
c
      do 209 k=1,kk
      do 209 j=1,jj
      do 209 i=1,ii
 209   p(i,j,k+1)=0.
c
      do 210 j=1,jj
      do 210 i=1,ii
      pbot(i  ,j  )=0.
      p(i  ,j  ,1)=0.
      do 210 k=1,kk
      dp(i  ,j  ,k   )=0.
 210  dp(i  ,j  ,k+kk)=0.

c
c --- initialize  u,ubavg,utotm,uflx,uflux,uflux2/3,uja,ujb  at points
c --- located upstream and downstream (in i direction) of p points.
c --- initialize  depthu,dpu,utotn,pgfx  upstream and downstream of p points
c --- as well as at lateral neighbors of interior u points.
c
      do 156 j=1,jj
      do 156 i=1,ii
      depthu(i,j)=0.
      utotn (i,j)=0.
c
      do 156 k=1,kk
      dpu(i,j,k   )=0.
      dpu(i,j,k+kk)=0.
c
 156  continue
c
      do 158 j=1,jj
      do 158 i=1,ii
      ubavg(i,j,1)=0.
      ubavg(i,j,2)=0.
      ubavg(i,j,3)=0.
      utotm (i,j)=0.
      uflux (i,j)=0.
      uflux2(i,j)=0.
c
      do 158 k=1,kk
      dpu(i,j,k   )=0.
      dpu(i,j,k+kk)=0.
      u(i,j,k   )=0.
 158  u(i,j,k+kk)=0.
c
c --- initialize  v,vbavg,vtotm,vflx,vflux,vflux2/3,via,vib  at points
c --- located upstream and downstream (in j direction) of p points.
c --- initialize  depthv,dpv,vtotn,pgfy  upstream and downstream of p points
c --- as well as at lateral neighbors of interior v points.
c
      do 166 i=1,ii
      do 166 j=1,jj
      depthv(i,j)=0.
      vtotn (i,j)=0.
c
      do 166 k=1,kk
      dpv(i,j,k   )=0.
      dpv(i,j,k+kk)=0.
c
 166  continue
c
      do 168 i=1,ii
      do 168 j=1,jj
      vbavg(i,j,1)=0.
      vbavg(i,j,2)=0.
      vbavg(i,j,3)=0.
      vtotm (i,j)=0.
      vflux (i,j)=0.
      vflux2(i,j)=0.
c
      do 168 k=1,kk
      dpv(i,j,k   )=0.
      dpv(i,j,k+kk)=0.
      v(i,j,k   )=0.
 168  v(i,j,k+kk)=0.

      return
      end
