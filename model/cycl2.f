      subroutine cycl2
c
c --- version 2.6
      implicit none
c
      include 'dimensions.h'
      include 'common_blocks.h'
      include 'stmt_funcs.h'

      real pb
      real poos(idm,jdm),pmer(idm,jdm)
      integer ih

      ih=ii/2

      pb=4300.

      do k=1,kk-1
      do i=1,ii-1
      do j=1,jj-1
      p(i,j,k+1)=100.*float(k)
      dpini(i,j,k)=p(i,j,k+1)-p(i,j,k)
      enddo
      enddo
      enddo

      do i=1,ii-1
      do j=1,jj-1
      p(i,j,kk+1)=pb
      poos(i,j)=pb
      pmer(i,j)=pb
      enddo
      enddo
      
      do i=1,ii-1
      do j=jj-34,jj
      poos(i,j)=-.5*pb*sin(pi*0.5*(j+17.-jj)/17.)+0.5*pb
      if(poos(i,j).lt.50.) poos(i,j)=50.
      enddo
      enddo
      
      do j=1,jj-1
      do i=1,ii-1
      pmer(i,j)=pb-pb*((ih+.5-i)/(ih+.5))**2
      if(pmer(i,j).lt.50.) pmer(i,j)=50.
      enddo
      enddo
      
      do i=1,ii1
      do j=1,jj1
      p(i,j,kk+1)=poos(i,j)
      if(pmer(i,j).lt.poos(i,j)) p(i,j,kk+1)=pmer(i,j)
       enddo
       enddo

      do k=1,kk-1
      do i=1,ii-1
      do j=1,jj-1
      if (p(i,j,k+1).gt.p(i,j,kk+1)) then
      p(i,j,k+1)=p(i,j,kk+1)
      dpini(i,j,k)=p(i,j,k+1)-p(i,j,k)
      endif
      enddo
      enddo
      enddo

      do j=1,jj-1
      do i=1,ii-1
       dpini(i,j,kk)=p(i,j,kk+1)-p(i,j,kk)
      enddo
      enddo

c
      do j=1,jj1
      do i=1,ii1
        do k=1,kk
          dpini(i,j,k)=dpini(i,j,k)*onem/(thref*(1-theta(k)))
        enddo
      enddo
      enddo

      do  j=1,jj
      do  i=1,ii
      pin(i,j,1)=0.
      enddo
      enddo
      	
      do 18 k=1,kk
      do 18 j=1,jj1
      do 18 i=1,ii1
  18  pin(i,j,k+1)=pin(i,j,k)+dpini(i,j,k)
 
      do k=1,kk
      do j=1,jj
      do i=1,ii
      delp5(i,j,k)=0.
      delp5(i,j,k+kk)=0.
      enddo
      enddo
      enddo

      end
