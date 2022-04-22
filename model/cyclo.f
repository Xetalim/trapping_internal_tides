      subroutine cyclo
c
c --- version 2.6
      implicit none
c
      include 'dimensions.h'
      include 'common_blocks.h'
      include 'stmt_funcs.h'

      real pb,pmax,pmin
      real poos(idm,jdm),pmer(idm,jdm)
      integer ih

      ih=ii/2

      pb=4300.

      do k=1,kk-1
      do i=1,ii-1
      do j=1,jj-1
      p(i,j,k+1)=100.*float(k)
      dp(i,j,k)=p(i,j,k+1)-p(i,j,k)
      enddo
      enddo
      enddo

      do k=2,kk
      pmax=0.
      pmin=1.e12
      do i=1,ii1
      do j=1,jj1
      if(p(i,j,k).gt.pmax) pmax=p(i,j,k)
      if(p(i,j,k).lt.pmin) pmin=p(i,j,k)
       enddo
       enddo
       print *, pmax, pmin, k , 'depth'
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
c     if(i.eq.1) print*, poos(i,j), j, 'bottom-oos'
      enddo
      enddo
      
      do j=1,jj-1
      do i=1,ii-1
      pmer(i,j)=pb-pb*((ih+.5-i)/(ih+.5))**2
      if(pmer(i,j).lt.50.) pmer(i,j)=50.
      if(j.eq.101) print*, pmer(i,j), i, 'bottom-mer'
      enddo
      enddo
      
      do i=1,ii1
      do j=1,jj1
      p(i,j,kk+1)=poos(i,j)
      if(pmer(i,j).lt.poos(i,j)) p(i,j,kk+1)=pmer(i,j)
       enddo
       enddo

      pmax=0.
      pmin=1.e12
      do i=1,ii1
      do j=1,jj1
      if(p(i,j,kk+1).gt.pmax) pmax=p(i,j,kk+1)
      if(p(i,j,kk+1).lt.pmin) pmin=p(i,j,kk+1)
       enddo
       enddo
       print *, pmax, pmin, 'depth-bottom'

      do k=1,kk-1
      do i=1,ii-1
      do j=1,jj-1
      if (p(i,j,k+1).gt.p(i,j,kk+1)) then
      p(i,j,k+1)=p(i,j,kk+1)
      dp(i,j,k)=p(i,j,k+1)-p(i,j,k)
      endif
      enddo
      enddo
      enddo

      do j=1,jj-1
      do i=1,ii-1
       dp(i,j,kk)=p(i,j,kk+1)-p(i,j,kk)
      enddo
      enddo

c
      do j=1,jj1
      do i=1,ii1
        do k=1,kk
          dp(i,j,k)=dp(i,j,k)*onem/(thref*(1-theta(k)))
        enddo
      enddo
      enddo
      
      do  j=1,jj
      do  i=1,ii
      p(i,j,1)=0.
      pin(i,j,1)=0.
      enddo
      enddo
      	
      do 18 k=1,kk
      do 18 j=1,jj1
      do 18 i=1,ii1
      p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
  18  pin(i,j,k+1)=p(i,j,k+1)
 
      do k=2,kk+1
      pmax=0.
      pmin=1.e12
      do i=1,ii1
      do j=1,jj1
      if(p(i,j,k).gt.pmax) pmax=p(i,j,k)
      if(p(i,j,k).lt.pmin) pmin=p(i,j,k)
       enddo
       enddo
       print *, pmax, pmin, k , 'depth'
       enddo


c --- Calculate montgomery potentials
      do j=1,jj1
          do i=1,ii1
            montg(i,j,1)=0.0     
	    do k=1,kk-1
              montg(i,j,k+1)=montg(i,j,k)
     .              -p(i,j,k+1)*(theta(k+1)-theta(k))*thref
	    enddo
            psikk(i,j)=montg(i,j,kk)
          enddo
      enddo

c    do k=1,kk
c     do j=1,jj
c         do i=1,ii
c             montg(i,j,k)=montg(i,j,k)-1.*psikk(i,j)
c    enddo
c         enddo
c     enddo
c
        mntin=montg(ih,1,1)

       do i=1,idm
       do j=1,jdm
       pbavg(i,j,1)=0.
       enddo
       enddo

      do j=1,jj1
          do i=1,ii1
c           psikk(i,j)=montg(i,j,kk)
            pbot(i,j)=p(i,j,kk+1)
      enddo
      enddo

      pmax=0.
      pmin=1.e12
      do i=1,ii1
      do j=1,jj1
      if(pbot(i,j).gt.pmax) pmax=pbot(i,j)
      if(pbot(i,j).lt.pmin) pmin=pbot(i,j)
       enddo
       enddo
       print *, pmax, pmin, 'pbot'

	    do k=1,kk
      do j=1,jj
          do i=1,ii
      u(i,j,k)=0.
      v(i,j,k)=0.
      potvor(i,j,k)=0.
      delp5(i,j,k)=0.
      enddo
      enddo
      enddo

      call waterdepth
c

      do k=1,3
      do i=1,idm
      do j=1,jdm
      ubavg(i,j,k)=0.
      vbavg(i,j,k)=0.
      enddo
      enddo
      enddo
c
      call timeslot

      end

      subroutine waterdepth
      include 'dimensions.h'
      include 'common_blocks.h'

c
      do i=1,ii
      do j=1,jj
      depthu(i,j)=0.
      depthv(i,j)=0.
      enddo
      enddo

      do k=1,2*kk
      do i=1,ii
      do j=1,jj
      dpu(i,j,k)=0.
      dpv(i,j,k)=0.
      enddo
      enddo
      enddo

      do 152 j=1,jj1
      do 152 i=2,ii1
152   depthu(i,j)=min(p(i,j,kk+1),p(i-1,j,kk+1))
c
      do 153 j=2,jj1
      do 153 i=1,ii1
153   depthv(i,j)=min(p(i,j,kk+1),p(i,j-1,kk+1))
c
      do i=1,ii1
      depthv(i,1)=p(i,1,kk+1)
      enddo
c
c --- define layer depth at  u,v  points
c
      do 180 k=1,kk
c
      do 157 j=1,jj1
      do 157 i=2,ii1
      dpu(i,j,k)=max(0.,
     .           min(depthu(i,j),.5*(p(i,j,k+1)+p(i-1,j,k+1)))-
     .           min(depthu(i,j),.5*(p(i,j,k  )+p(i-1,j,k  ))))
 157   dpu(i,j,k+kk)=dpu(i,j,k)
c
      do 155 j=2,jj1
      do 155 i=1,ii1
      dpv(i,j,k)=max(0.,
     .           min(depthv(i,j),.5*(p(i,j,k+1)+p(i,j-1,k+1)))-
     .           min(depthv(i,j),.5*(p(i,j,k  )+p(i,j-1,k  ))))
 155  dpv(i,j,k+kk)=dpv(i,j,k)
c
      do i=1,ii1
      dpv(i,1,k)=max(0.,
     .           min(depthv(i,1),p(i,1,k+1))-
     .           min(depthv(i,1),p(i,1,k  )))
      dpv(i,1,k+kk)=dpv(i,1,k)
      enddo
c
 180  continue
      return
      end
c
      subroutine timeslot
      include 'dimensions.h'
      include 'common_blocks.h'

c --- copy variables into second time slot
c
      do k=1,kdm
      do j=1,jdm
      do i=1,idm
        u(i,j,k+kdm)=u(i,j,k)
        v(i,j,k+kdm)=v(i,j,k)
        dp(i,j,k+kdm)=dp(i,j,k)
        delp5(i,j,k+kdm)=delp5(i,j,k)
        dpini(i,j,k)=dp(i,j,k)
      end do
      end do
      end do
      do k=2,3
        do j=1,jdm
          do i=1,idm
            pbavg(i,j,k)=pbavg(i,j,1)
            ubavg(i,j,k)=ubavg(i,j,1)
            vbavg(i,j,k)=vbavg(i,j,1)
          enddo
        enddo
      enddo
c      
      end
