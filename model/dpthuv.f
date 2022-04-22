      subroutine dpthuv
c
c --- define water depth (bottom pressure) at  u,v  points and barotp.pot.vort.
c
c --- version 2.7
      implicit none
c
      include 'dimensions.h'
      include 'common_blocks.h'
c
      real uvdep,a,b
c
c --- function for determining depth at u,v points
      uvdep(a,b)=min(a,b)
c
      do j=1,jj
      do i=1,ii
      depthu(i,j)=0.
      depthv(i,j)=0.
      pvtrop(i,j)=0.
      enddo
      enddo
    
      do 151 j=1,jj1
      do 151 i=2,ii1
 151  depthu(i,j)=uvdep(pbot(i,j),pbot(i-1,j))
c
      do 152 i=1,ii1
      do 152 j=2,jj1
 152  depthv(i,j)=uvdep(pbot(i,j),pbot(i,j-1))
c
      do i=1,ii1
      depthv(i,1)=pbot(i,1)
      enddo
c
      do 153 j=2,jj1
      do 153 i=2,ii1
 153  pvtrop(i,j)=corio(i,j)*4./(pbot(i,j  )+pbot(i-1,j  )
     .                          +pbot(i,j-1)+pbot(i-1,j-1))
c
      do i=2,ii1
      pvtrop(i,1)=corio(i,1)*2./(pbot(i,1)+pbot(i-1,1))
      pvtrop(i,jj)=corio(i,jj)*2./(pbot(i,jj1)+pbot(i-1,jj1))
      enddo

      do j=2,jj1
      pvtrop(1,j)=corio(1,j)*2./(pbot(1,j)+pbot(1,j-1))
      pvtrop(ii,j)=corio(ii,j)*2./(pbot(ii1,j)+pbot(ii1,j-1))
      enddo
      pvtrop(1,1)=corio(1,1)/pbot(1,1)
      pvtrop(1,jj)=corio(1,jj)/pbot(1,jj1)
      pvtrop(ii,1)=corio(ii,1)/pbot(ii1,1)
      pvtrop(ii,jj)=corio(ii,jj)/pbot(ii1,jj1)
      return
      end
