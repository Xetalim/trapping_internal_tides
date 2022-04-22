       program read016

       double precision t(53,322,43)
       double precision y(53,322,43)
       dimension x(53,322,43,16),v(53,322,42)
       dimension a(53,322,42),fi(53,322,42)
       integer ia,ja,ka,ib,jb,kb
       real vmx,tmx,tmn

       character path*44

       path = '/Users/4302001/surfdrive/WGF/2022/Model/env/'

       open (99,file=
     &path // "average_00004"
     &,form="unformatted")
       read (99) t
       close (99)

       onem=98060.

       open (99,file=
     &path // "analyse_00001"
     &,form="unformatted")
       read (99) y
       close (99)

       do i=1,53
       do j=1,322
       do k=1,43
       x(i,j,k,1)=(y(i,j,k)-t(i,j,k))/980.6
       enddo
       enddo
       enddo

       open (99,file=
     &path // "analyse_00002"
     &,form="unformatted")
       read (99) y
       close (99)

       do i=1,53
       do j=1,322
       do k=1,43
       x(i,j,k,2)=(y(i,j,k)-t(i,j,k))/980.6
       enddo
       enddo
       enddo


       open (99,file=
     &path // "analyse_00003"
     &,form="unformatted")
       read (99) y
       close (99)

       do i=1,53
       do j=1,322
       do k=1,43
       x(i,j,k,3)=(y(i,j,k)-t(i,j,k))/980.6
       enddo
       enddo
       enddo


       open (99,file=
     &path // "analyse_00004"
     &,form="unformatted")
       read (99) y
       close (99)

       do i=1,53
       do j=1,322
       do k=1,43
       x(i,j,k,4)=(y(i,j,k)-t(i,j,k))/980.6
       enddo
       enddo
       enddo


       open (99,file=
     &path // "analyse_00005"
     &,form="unformatted")
       read (99) y
       close (99)

       do i=1,53
       do j=1,322
       do k=1,43
       x(i,j,k,5)=(y(i,j,k)-t(i,j,k))/980.6
       enddo
       enddo
       enddo


       open (99,file=
     &path // "analyse_00006"
     &,form="unformatted")
       read (99) y
       close (99)

       do i=1,53
       do j=1,322
       do k=1,43
       x(i,j,k,6)=(y(i,j,k)-t(i,j,k))/980.6
       enddo
       enddo
       enddo


       open (99,file=
     &path // "analyse_00007"
     &,form="unformatted")
       read (99) y
       close (99)

       do i=1,53
       do j=1,322
       do k=1,43
       x(i,j,k,7)=(y(i,j,k)-t(i,j,k))/980.6
       enddo
       enddo
       enddo


       open (99,file=
     &path // "analyse_00008"
     &,form="unformatted")
       read (99) y
       close (99)

       do i=1,53
       do j=1,322
       do k=1,43
       x(i,j,k,8)=(y(i,j,k)-t(i,j,k))/980.6
       enddo
       enddo
       enddo


       open (99,file=
     &path // "analyse_00009"
     &,form="unformatted")
       read (99) y
       close (99)

       do i=1,53
       do j=1,322
       do k=1,43
       x(i,j,k,9)=(y(i,j,k)-t(i,j,k))/980.6
       enddo
       enddo
       enddo


       open (99,file=
     &path // "analyse_00010"
     &,form="unformatted")
       read (99) y
       close (99)

       do i=1,53
       do j=1,322
       do k=1,43
       x(i,j,k,10)=(y(i,j,k)-t(i,j,k))/980.6
       enddo
       enddo
       enddo


       open (99,file=
     &path // "analyse_00011"
     &,form="unformatted")
       read (99) y
       close (99)

       do i=1,53
       do j=1,322
       do k=1,43
       x(i,j,k,11)=(y(i,j,k)-t(i,j,k))/980.6
       enddo
       enddo
       enddo


       open (99,file=
     &path // "analyse_00012"
     &,form="unformatted")
       read (99) y
       close (99)

       do i=1,53
       do j=1,322
       do k=1,43
       x(i,j,k,12)=(y(i,j,k)-t(i,j,k))/980.6
       enddo
       enddo
       enddo


       open (99,file=
     &path // "analyse_00013"
     &,form="unformatted")
       read (99) y
       close (99)

       do i=1,53
       do j=1,322
       do k=1,43
       x(i,j,k,13)=(y(i,j,k)-t(i,j,k))/980.6
       enddo
       enddo
       enddo


       open (99,file=
     &path // "analyse_00014"
     &,form="unformatted")
       read (99) y
       close (99)

       do i=1,53
       do j=1,322
       do k=1,43
       x(i,j,k,14)=(y(i,j,k)-t(i,j,k))/980.6
       enddo
       enddo
       enddo


       open (99,file=
     &path // "analyse_00015"
     &,form="unformatted")
       read (99) y
       close (99)

       do i=1,53
       do j=1,322
       do k=1,43
       x(i,j,k,15)=(y(i,j,k)-t(i,j,k))/980.6
       enddo
       enddo
       enddo


       open (99,file=
     &path // "analyse_00016"
     &,form="unformatted")
       read (99) y
       close (99)

       do i=1,53
       do j=1,322
       do k=1,43
       x(i,j,k,16)=(y(i,j,k)-t(i,j,k))/980.6
       enddo
       enddo
       enddo

      do l=1,16
      do k=2,42
      do i=1,53
      do j=1,322
      x(i,j,k,l)=x(i,j,k,l)+x(i,j,k-1,l)
      enddo
      enddo
      enddo
      enddo

      pi=4.*atan(1.)

      do k=1,42
      do i=1,53
      do j=1,322
      xm=0.
      xsq=0.
      do l=1,16
      xm=xm+x(i,j,k,l)
      xsq=xsq+x(i,j,k,l)*x(i,j,k,l)
      enddo
      xm=xm/16.
      xc=0.
      xs=0.
      v(i,j,k)=0.
      do l=1,16
      xc=xc+(x(i,j,k,l)-xm)*cos(l*pi/4)
      xs=xs+(x(i,j,k,l)-xm)*sin(l*pi/4)
      enddo
      xc=xc/8.
      xs=xs/8.
      do l=1,16
      v(i,j,k)=v(i,j,k)+(xm+xc*cos(l*pi/4)+xs*sin(l*pi/4)-x(i,j,k,l))*
     &(xm+xc*cos(l*pi/4)+xs*sin(l*pi/4)-x(i,j,k,l))
      enddo
      v(i,j,k)=v(i,j,k)/xsq
      a(i,j,k)=sqrt(xc**2+xs**2)
      fi(i,j,k)=atan2(xs,xc)
      enddo
      enddo
      enddo

      do k=38,28,-1
      print *, k, fi(36,265,k),fi(35,265,k),fi(37,265,k)
      enddo
      do k=38,28,-1
      print *, k, fi(38,265,k),fi(34,265,k)
      enddo

      open (unit = 99, file='iwa038-04-harmonie', form='unformatted')
      write (99) a
      write (99) fi
      write (99) v
      close (99)
 666  continue

       end
