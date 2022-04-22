      subroutine avgout 
c
      include 'dimensions.h'
      include 'common_blocks.h'
c
      character flnm*13, time_str*5

! schrijf weg:
      write (time_str, '(i5.5)') NINT (time)

      flnm(:8) = 'average_:'
      flnm(9:13)=time_str
      open (unit = 99, file = flnm(:13), form = 'unformatted')
!! 3D velden:
         write (99) dpm
         write (99) vp1
         write (99) vp5
      close (99)

      end
