      subroutine restin(nstep0,time0)
c
      include 'dimensions.h'
      include 'common_blocks.h'
c

c --- start from restart file
c
         ni=11
c
            OPEN (unit=ni,file=flnmrsi,status='old',form='unformatted')
               READ (ni) nstep0,time0,u
               READ (ni) nstep0,time0,v
               READ (ni) nstep0,time0,dp
               READ (ni) nstep0,time0,ubavg,vbavg,pbavg
               read (ni) nstep0,time0,pbot,psikk
            CLOSE (unit=ni)

      return
      end
