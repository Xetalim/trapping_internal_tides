      subroutine boundub(q,qtn,qbo,qto,l)
c
      implicit none
c
      include 'dimensions.h'
      include'common_blocks.h'
c
      real q(idm,jdm,2*kdm),z(idm,jdm),qtn(idm,jdm)
     .,qbo(idm,jdm),qto(idm,jdm),alfa
c
      do j=1,jj
      do i=1,ii
       z(i,j)=q(i,j,l)
      enddo
      enddo

      do j=1,jj
       z(ii,j)=0.     
       z(1,j)=0.
      enddo
      do i=2,ii1
       z(i,jj)=0.    
      do j=1,13
      alfa=((14.-j)/13.)**2
       z(i,j)=z(i,j)-alfa*(1.+wbaro)*dlt*(qtn(i,j)+qbo(i,j)+qto(i,j))
     
      enddo
      enddo

      do j=1,jj
      do i=1,ii
       q(i,j,l)=z(i,j)
      enddo
      enddo

      return
      end
