      subroutine restout(itime)
c
      include 'dimensions.h'
      include 'common_blocks.h'
c

c --- output to restart file
      character flnm*13, time_str*5
c
      write (time_str, '(i5.5)') NINT (time)

         no=12
      flnm(:8) = 'restart_:'
      flnm(9:13)=time_str

c
       if(time.le.5.) then
       OPEN (unit=no,file=flnm(:13),status='unknown',form='unformatted')
c           rewind no
               write (no) nstep,time,u
               write (no) nstep,time,v
               write (no) nstep,time,dp
               write (no) nstep,time,ubavg,vbavg,pbavg
               write (no) nstep,time,pbot,psikk
            CLOSE (unit=no)
       endif


       if(time.gt.5.) then
       OPEN (unit=no,file=flnmrso,status='unknown',form='unformatted')
            rewind no
               write (no) nstep,time,u
               write (no) nstep,time,v
               write (no) nstep,time,dp
               write (no) nstep,time,ubavg,vbavg,pbavg
               write (no) nstep,time,pbot,psikk
            CLOSE (unit=no)
      endif

      return
      end
