      subroutine BinOut (k1n,n)
c
      include 'dimensions.h'
      include 'common_blocks.h'
c
      character flnm*14, time_str*5

! schrijf weg:
      if(time.gt.4..and.time.le.5.) then
      write (time_str, '(i5.5)') NINT (16*time-16*4)
      endif

      if(time.gt.9..and.time.le.10.) then
      write (time_str, '(i5.5)') NINT (16*time-16*8)
      endif

      if(time.gt.11..and.time.le.12.) then
      write (time_str, '(i5.5)') NINT (16*time-16*9)
      endif

      if(time.gt.15..and.time.le.16.) then
      write (time_str, '(i5.5)') NINT (16*time-16*12)
      endif

      flnm(:9) = 'analysis_:'
      flnm(10:14)=time_str
      open (unit = 99, file = flnm(:14), form = 'unformatted')
!! 3D velden:
         write (99) dp5
      close (99)

      end subroutine BinOut
