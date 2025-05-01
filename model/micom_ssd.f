      program micom
c
c --- --------------------------------------
c --- miami isopycnic coordinate ocean model
c
c ---         v e r s i o n    3.0 ; imau version
c --- --------------------------------------
c
      implicit none
c
      include 'dimensions.h'
      include 'common_blocks.h'
      include 'stmt_funcs.h'
c
      real day1,day2,x,x1,time0,timav,cold
      integer nstep0,ni,no,ia,ib,ja,jb,ka,kb,itijd,
     .mo0,mo1,mo2,mo3,itime,iaa,jaa,ih
      real umax,umin,pmax,pmin,uabs,tstp,tsec
      real start_time, current_time, elapsed, estimated_total, eta
      integer total_steps
c
      data nstep0,time0/0,0./

      namelist /micom_nml/ day1, day2, trcrin, trcout, lpout
c
      call init_block_data 

c
      ih=ii/2

      open(10, file='micom.nml')
      read(10, nml=micom_nml)
      close(10)
      write (lp,101) thkdff,temdff,veldff,viscos,diapyc,vertmx
 101  format (' turb. flux parameters:',1p/
     .  ' thkdff,temdff,veldff =',3e9.2/
     .  ' viscos,diapyc,vertmx =',3e9.2)
c
      realism=.false.
      forc=.false.
      thc=.false.
      tstp=86400./baclin
      print *, tstp, 'tstp'
c
c --- 'lstep'= number of barotropic steps per baroclinic time step.
c --- lstep   m u s t   be even.
c
      lstep=baclin/batrop
      lstep=2*((lstep+1)/2)
      dlt=baclin/lstep
      write (lp,'(i4,'' barotropic steps per baroclinic time step'')')
     .  lstep
c
c --- set up parameters defining the geographic environment
c
	
      call geopar
      	
c
c --- subtract constant 'thbase' from theta to reduce roundoff errors
c
      do 14 k=1,kk
 14   theta(k)=sigma(k)/1000.-thbase
c
c --- model is to be integrated from time step 'nstep1' to 'nstep2'
      nstep1=day1*86400./baclin + .0001
      nstep2=day2*86400./baclin + .0001
c
      IF (nstep1.le.0) THEN
c
c --- set up initial conditions
c
            call cyclo
c        goto 6666
         delt1=baclin
c
      ELSE  !  nstep1 > 0
c
c --- start from restart file
c
          call restin(nstep0,time0)
c
         nstep0=time0*86400./baclin+.0001
         write (lp,111) nstep0,nstep1
 111  format (9x,'time step in restart file -',i9,5x,'wanted -',i9)
         IF (nstep0.lt.nstep1) STOP 'wrong time in restart file.'
c
c if restart file allows restart at a time greater than nstep1, use it
         nstep2=nstep0+nstep2-nstep1
         nstep1=nstep0
c
         call cycl2
         delt1=baclin+baclin
      ENDIF !  nstep1 > 0  or  = 0
c
      time=0.
c --- set barotp.pot.vort. and layer thickness (incl.bottom pressure) 
c --- at u,v points
c
      call dpthuv
c
      do 16 m=1,1
      mm=(m-1)*kk
c
      do 18 j=1,jj
      do 18 k=1,kk
      do 18 i=1,ii
 18   p(i,j,k+1)=p(i,j,k)+dp(i,j,k+mm)
c
      call dpudpv(p,depthu,depthv,dpu(1,1,1+mm),dpv(1,1,1+mm))
 16   continue
c
      nstep=nstep1
      write (lp,'(/2(a,f8.1),2(a,i9),a/)') 'model starts at day',
     .   time0,', goes to day',time0+day2-day1,'   (steps',nstep1,
     .   ' --',nstep2,')'
      timav=time0
c
c --- ---------------------
c --- main loop starts here
c --- ---------------------
c
      call vort1(1,0,0)
      do 705 j=1,jj
      do 705 i=1,ii
      ssh(i,j)=montg(i,j,1)+thref*pbavg(i,j,1)
 705  continue
c     call sshini
c     goto 5555
c --- letter 'm' refers to mid-time level (example: dp(i,j,km) )
c --- letter 'n' refers to old and new time level
c
      iaa=2
      jaa=2

      do k=1,kdm
      do j=1,jdm
      do i=1,idm
      dpm(i,j,k)=0.
      enddo
      enddo
      enddo

c     do j=1,jdm
c     do i=1,idm
c     ssm(i,j)=0.
c     enddo
c     enddo

      do j=1,jdm
      vp1(j)=0.
      vp5(j)=0.
      enddo
      total_steps = max(1, nstep2-nstep1)
      DO itijd = 1, total_steps
c     do itijd=1,80
c     do itijd=1,480
c     do itime = 1,140
      itime=itijd+nstep0
      if (itijd.eq.1) print *, itijd, itime, ' tijd-start'
      m=mod(nstep  ,2)+1
      n=mod(nstep+1,2)+1
      mm=(m-1)*kk
      nn=(n-1)*kk
      k1m=1+mm
      k1n=1+nn
c
      nstep=nstep+1
      time=time0+(nstep-nstep0)*baclin/86400.
      diagno=.false.
      if (mod(time+.0001,8*diagfq).lt..0002) diagno=.true.
      if (nstep.ge.nstep2) diagno=.true.
c     if (itime.le.4) diagno=.true.
c
      if (diagno) then 
       write(*,*)'itime=',itime
      endif

      IF (MOD((itime * baclin ), 86400.) < baclin) THEN
      do k=1,kdm
      do j=1,jdm
      do i=1,idm
      delp5(i,j,k)=0.
      delp5(i,j,k+kdm)=0.
      enddo
      enddo
      enddo
      ENDIF

      if(mod(itime,120).eq.0) then
      call cpu_time(current_time)
      elapsed = current_time - start_time
      estimated_total = elapsed / (itime) * total_steps
      eta = estimated_total - elapsed
            print '(A,I4,A,I4,A,F6.2,A)',
     &     "itime ", itime, "/", total_steps,
     &     " | ETA: ", eta, " s"

!     print *, itime, 'itime'
      endif
      call cnuity(m,mm,nn)
c     print *, dp(1,240,30+nn),dp(52,240,30+nn), 'dp-cn'

      if(mod(itime,120).eq.0) then
      pmax=-onem*4000.
      pmin=onem*4000.
      do k=1,kk
      do j=1,jj1
       do i=1,ii1
      if(dp(i,j,k+nn)-dpini(i,j,k).gt.pmax) then 
       pmax=dp(i,j,k+nn)-dpini(i,j,k)    
       ia=i
       ja=j
       ka=k
       endif
      if(dp(i,j,k+nn)-dpini(i,j,k).lt.pmin) then
       pmin=dp(i,j,k+nn)-dpini(i,j,k)   
       ib=i
       jb=j
       kb=k
       endif
      enddo
      enddo
      enddo
      if(pmax/onem*pmax/onem.gt.pmin/onem*pmin/onem) then
      print *, pmax/onecm, ia, ja, ka, 'dpmax'
      print *, dp(ia,ja,ka+nn)/onem, 
     .dpini(ia,ja,ka)/onem, 'dp'
      else
      print *, pmin/onecm, ib, jb, kb, 'dpmin'
      print *, dp(ib,jb,kb+nn)/onem, 
     .dpini(ib,jb,kb)/onem, 'dp'
      endif
      endif

      if(mod(itime,120).eq.0) then
      umax=-1.e12
      umin=1.e12
      do j=241,jj
       do i=1,ii
      if(vbavg(i,j,n).gt.umax) then 
       umax=vbavg(i,j,n)
       ia=i
       ja=j
       endif
      if(vbavg(i,j,n).lt.umin) then 
       umin=vbavg(i,j,n)
       ib=i
       jb=j
       endif
      enddo
      enddo
      print *, 1.e-2*umax, ia, ja, 'vtr'
      print *, 1.e-2*umin, ib, jb, 'vtr'
      endif

c     print *, 'tsadvc'
      call tsadvc(m,n,mm,nn,k1m,k1n)
c     print*, 'vort'
      call vort1(m,mm,nn)
c     print*, 'momeq1'
      call momeq1(itime,m,n,mm,nn,k1m)
c     print*, 'momeq2'
      call momeq2(itime,m,n,mm,nn,k1n,k1m)
c     print *, 'barotp'

      if(mod(itime,120).eq.0) then
      umax=0.
      do k=1,kk
      do j=1,jj
       do i=1,ii
       uabs=sqrt(u(i,j,k+nn)**2+v(i,j,k+nn)**2)
      if(uabs.gt.umax) then 
       umax=uabs
       ia=i
       ja=j
       ka=k
       endif
      enddo
      enddo
      enddo
      print *, umax, ia, ja, ka, 'ucl'
      endif

      call barotp(itime,n)

      if(mod(itime,120).eq.0) then
      pmax=-onem*4000.
      pmin=onem*4000.
      do k=1,3
      do j=1,jj1
       do i=1,ii1
      if(pbavg(i,j,k).gt.pmax) then 
       pmax=pbavg(i,j,k)
       ia=i
       ja=j
       endif
      if(pbavg(i,j,k).lt.pmin) then
       pmin=pbavg(i,j,k)
       ib=i
       jb=j
       endif
      enddo
      enddo
      enddo
c     if(pmax/onecm*pmax/onecm.gt.pmin/onecm*pmin/onecm) then
      print *, pmax/onecm, ia, ja, 'zmax'
c     else
      print *, pmin/onecm, ib, jb, 'zmin'
c     endif
      tsec =max((itime-2)*baclin+lstep*batrop*2.,0.)
      print *, (pbavg(ih,1,n)+montg(ih,1,1)-mntin)/onecm, 
     .sin(tsec*2.*pi/(12.*3600.)), n, 'zfor'

      umax=0.
      do k=1,3
      do j=1,jj
       do i=1,ii
       uabs=sqrt(ubavg(i,j,k)**2)
      if(uabs.gt.umax) then 
       umax=uabs
       ia=i
       ja=j
       endif
       enddo
       enddo
      enddo
      print *, umax, ia, ja, 'ubavgmax'
      iaa=ia
      jaa=ja

      umax=0.
      do k=1,3
      do j=1,jj
       do i=1,ii
       uabs=sqrt(vbavg(i,j,k)**2)
      if(uabs.gt.umax) then 
       umax=uabs
       ia=i
       ja=j
       endif
      enddo
      enddo
      enddo
      print *, umax, ia, ja, 'vbavgmax'
      endif
c
c     if(mod(itime,120).eq.0) then
c     print *, pbavg(1,14,n),pbavg(ii1,14,n), 'pbavg'
c     print *, vbavg(1,14,n),vbavg(ii1,14,n), 'vbavg'
c     print *, ubavg(2,14,n),ubavg(ii1,14,n), 'ubavg'
c     print *, dp(1,14,2+nn)-dpini(1,14,2),dp(ii1,14,2+nn)
c    .-dpini(ii1,14,2), 'dp'
c     print *, v(1,14,2+nn),v(ii1,14,2+nn), 'v'
c     print *, u(2,14,2+nn),u(ii1,14,2+nn), 'u'
c     endif
      if(mod(itime,120).eq.0) then
c -----------------------------------------------------------------
c
c --- output and diagnostic calculations
      do 604 j=1,jj1
      do 604 i=1,ii1
      ssh(i,j)=(montg(i,j,1)+thref*pbavg(i,j,n))/onecm
 604  continue
      print *, itime, ssh(ih,1), 'ssh'
      print *, ssh(ih,4),ssh(ih,7),ssh(ih,10)
      print *, ssh(ih,14),ssh(ih,20),ssh(ih,30)
      endif
c
c -----------------------------------------------------------------
c
      if (diagno) then
c
      write (lp,100) nstep,int((time+.001)/360.),mod(time+.001,360.)
 100  format (' time step',i9,9x,'y e a r',i6,9x,'d a y',f9.2)
c
      timav=time
c
c --- output to restart file
c
      call restout(itime)
c
      write (lp,105) nstep
 105  format (' step',i9,' -- archiving completed --')
c
      ENDIF
      IF (itime * baclin.gt.3. * 86400.) THEN
   
      do k=1,kdm
      do j=1,jdm
      do i=1,idm
      dpm(i,j,k)=dpm(i,j,k)+(dp(i,j,k+nn)-delp5(i,j,k+nn))/tstp
      enddo
      enddo
      enddo

      do j=1,jdm
      vp1(j)=vp1(j)+vpfl1(j)/tstp
      vp5(j)=vp5(j)+vpfl5(j)/tstp
      enddo

c     do j=1,jdm
c     do i=1,idm
c     ssm(i,j)=ssm(i,j)+ssh(i,j)/tstp
c     enddo
c     enddo

c     endif
      IF (MOD((itime * baclin ), 86400.) < baclin) THEN
      call avgout
      do k=1,kdm
      do j=1,jdm
      do i=1,idm
      dpm(i,j,k)=0.
      enddo
      enddo
      enddo

c     do j=1,jdm
c     do i=1,idm
c     ssm(i,j)=0.
c     enddo
c     enddo

      do j=1,jdm
      vp1(j)=0.
      vp5(j)=0.
      enddo

      ENDIF
      endif

c     IF (itime * baclin.gt.9. * 86400..and.
c    &itime*baclin.le.10. * 86400.) THEN
c     IF (MOD (itime * baclin, diagfq * 86400./2.) < baclin) THEN
c        CALL BinOut (k1n,n)
c     ENDIF
c     endif

      IF (itime * baclin.gt.4. * 86400..and.
     &itime*baclin.le.5. * 86400.) THEN
      IF (MOD (itime * baclin, diagfq * 86400./2.) < baclin) THEN
      do k=1,kdm
      do j=1,jdm
      do i=1,idm
      dp5(i,j,k)=dp(i,j,k+nn)-delp5(i,j,k+nn)
      enddo
      enddo
      enddo
         CALL BinOut (k1n,n)
      ENDIF
      endif

c     IF (itime * baclin.gt.9. * 86400..and.
c    &itime*baclin.le.10. * 86400.) THEN
c     IF (MOD (itime * baclin, diagfq * 86400./2.) < baclin) THEN
c     do k=1,kdm
c     do j=1,jdm
c     do i=1,idm
c     dp5(i,j,k)=dp(i,j,k+nn)-delp5(i,j,k+nn)
c     enddo
c     enddo
c     enddo
c        CALL BinOut (k1n,n)
c     ENDIF
c     endif

c     IF (itime * baclin.gt.11. * 86400..and.
c    &itime*baclin.le.12. * 86400.) THEN
c     IF (MOD (itime * baclin, diagfq * 86400./2.) < baclin) THEN
c     do k=1,kdm
c     do j=1,jdm
c     do i=1,idm
c     dp5(i,j,k)=dp(i,j,k+nn)-delp5(i,j,k+nn)
c     enddo
c     enddo
c     enddo
c        CALL BinOut (k1n,n)
c     ENDIF
c     endif

      IF (itime * baclin.gt.15. * 86400..and.
     &itime*baclin.le.16. * 86400.) THEN
      IF (MOD (itime * baclin, diagfq * 86400./2.) < baclin) THEN
      do k=1,kdm
      do j=1,jdm
      do i=1,idm
      dp5(i,j,k)=dp(i,j,k+nn)-delp5(i,j,k+nn)
      enddo
      enddo
      enddo
         CALL BinOut (k1n,n)
      ENDIF
      endif

c     IF (itime * baclin.gt.39. * 86400..and.
c    &itime*baclin.le.40. * 86400.) THEN
c     IF (MOD (itime * baclin, diagfq * 86400./2.) < baclin) THEN
c        CALL BinOut (k1n,n)
c     ENDIF
c     endif

c     IF (itime * baclin.gt.9. * 86400.) THEN
c     do 704 j=1,jj
c     do 704 i=1,ii
c     ssh(i,j)=montg(i,j,1)+thref*pbavg(i,j,n)
c704  continue
c     endif
c     IF (itime * baclin.gt.9.5 * 86400.) THEN
c     IF (MOD (itime * baclin, diagfq * 86400./5.) < baclin) THEN
c        CALL sshout(2) 
c     ENDIF
c     endif
c     IF (itime * baclin.gt.15.5 * 86400..and.
c    .    itime * baclin.le.16.0 * 86400) THEN
c     IF (MOD (itime * baclin, diagfq * 86400./5.) < baclin) THEN
c     do 504 j=1,jj
c     do 504 i=1,ii
c     ssh(i,j)=montg(i,j,1)+thref*pbavg(i,j,n)
c504  continue
c        CALL sshout(2) 
c     ENDIF
c     endif
      delt1=baclin+baclin
      ENDDO ! itime-loop
 5555 continue
c     call avgout
 6666 continue
      WRITE (6,*) '(normal end of MICOM)'

      end
