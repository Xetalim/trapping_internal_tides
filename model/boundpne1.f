      subroutine boundpne1(p2,p1,p0)
c
      implicit none
c
      include 'dimensions.h'
      include 'common_blocks.h'
c
      real p0(idm,jdm,kdm+1),cx(kdm)
      real epps
      real p2(idm,jdm,kdm+1),p1(idm,jdm,kdm+1)
c
       epps=1.e-30

      do j=14,1,-1
      do i=1,ii-1
      do k=2,kk
      cx(k)=-(p2(i,j,k)-p0(i,j,k))
     &/(p2(i,j,k)+p0(i,j,k)-2.*p1(i,j+1,k)+epps)
      if(cx(k).lt.0.) then   
      cx(k)=-(p2(i,j+1,k)-p0(i,j+1,k))*scuv(i,j+1)
     &/(p1(i,j+2,k)-0.5*(p2(i,j+1,k)+p0(i,j+1,k))+epps)/delt1
      cx(k)=max(cx(k),-2.*scuv(i,j+1)/delt1)
      cx(k)=min(0.,cx(k))
      p2(i,j,k)=(p0(i,j,k)-cx(k)*delt1/scuv(i,j+1)*(p1(i,j+1,k)-
     &0.5*p0(i,j,k)))
     &/(1.-0.5*cx(k)*delt1/scuv(i,j+1))
      endif
      enddo
      do k=2,kk
      if(p2(i,j,k).lt.p2(i,j,k-1)) p2(i,j,k)=p2(i,j,k-1)
      enddo
      enddo

      enddo
 
      return
      end
