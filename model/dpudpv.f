      subroutine dpudpv(p,depthu,depthv,dpu,dpv)
c
c --- version 2.7
      implicit none
c
c --- ----------------------------------
c --- define layer depth at  u,v  points
c --- ----------------------------------
c
      include 'dimensions.h'
c
      real p(idm,jdm,kdm+1),dpu(idm,jdm,kdm),dpv(idm,jdm,kdm),
     .     depthu(idm,jdm),depthv(idm,jdm)
c
      do 180 k=1,kk
c
      do 154 j=1,jj1
      do 154 i=2,ii1
 154  dpu(i,j,k)=max(0.,
     .           min(depthu(i,j),.5*(p(i,j,k+1)+p(i-1,j,k+1)))-
     .           min(depthu(i,j),.5*(p(i,j,k  )+p(i-1,j,k  ))))
c
      do 155 j=2,jj1
      do 155 i=1,ii1
 155  dpv(i,j,k)=max(0.,
     .           min(depthv(i,j),.5*(p(i,j,k+1)+p(i,j-1,k+1)))-
     .           min(depthv(i,j),.5*(p(i,j,k  )+p(i,j-1,k  ))))
c
      do i=1,ii1
      dpv(i,1,k)=max(0.,
     .           min(depthv(i,1),1.5*p(i,1,k+1)-.5*p(i,2,k+1))-
     .           min(depthv(i,1),1.5*p(i,1,k)-.5*p(i,2,k)   ))
      enddo
c
 180  continue
      return
      end
