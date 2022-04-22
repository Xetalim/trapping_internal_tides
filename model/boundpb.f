      subroutine boundpb(itime,q,nl,lll)
c
      implicit none
c
      include 'dimensions.h'
      include 'common_blocks.h'
c
      integer itime,lll,nl,iha
      real q(idm,jdm,3),z(idm,jdm),tsec,alfa
c
      do j=1,jj
      do i=1,ii
       z(i,j)=q(i,j,nl)
      enddo
      enddo

      do j=1,jj
       z(ii,j)=0.     
      enddo
      do i=1,ii
       z(i,jj)=0.    
      enddo
      tsec =max((itime-2)*baclin+lll*batrop*2.,0.)
      iha=(ii-2)/2
      do j=1,13
      do i=1,ii-1
       alfa=((14.-j)/13.)**2
c      z(i,j)=z(i,j)*(1.-(1.-0.1*(j-1))**2)+(1.-0.1*(j-1))**2*
       z(i,j)=z(i,j)*(1.-alfa)+alfa*
     .onecm*1.*sin(tsec*2.*pi/(12.*3600.))
     .*(1.-1.*exp(-float(itime)/240.))
      enddo
      enddo

      do j=1,jj
      do i=1,ii
       q(i,j,nl)=z(i,j)
      enddo
      enddo

      return
      end
