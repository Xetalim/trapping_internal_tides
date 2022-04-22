      subroutine boundp(q,l)
c
      implicit none
c
      include 'dimensions.h'
c
      real q(idm,jdm,2*kdm),p(idm,jdm,kdm)
c
      do k=1,kk
      m=l+k-1
      do j=1,jj
      do i=1,ii
       p(i,j,k)=q(i,j,m)
      enddo
      enddo
      enddo

      do k=1,kk
      do j=1,jj
       p(ii,j,k)=0.       
      enddo
      do i=1,ii
       p(i,jj,k)=0.
      enddo
      enddo

      do k=1,kk
      m=l+k-1
      do j=1,jj
      do i=1,ii
       q(i,j,m)=p(i,j,k)
      enddo
      enddo
      enddo

      return
      end
