      subroutine sshini 
c
      include 'dimensions.h'
      include 'common_blocks.h'
c
      character flnm*13, time_str*5

! schrijf weg:
      write (time_str, '(i5.5)') NINT (0.)

      flnm(:8) = 'sshoutp_:'
      flnm(9:13)=time_str
      open (unit = 99, file = flnm(:13), form = 'unformatted')
         write (99) ssh
      close (99)

      end subroutine sshini
