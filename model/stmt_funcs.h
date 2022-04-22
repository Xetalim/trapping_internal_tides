c-----------------------------------------------------------------------------
      real dist1,alat1,grid,a,b,athird,c1,c2,c3,c4,c5,c6,c7,r,s,t
      real alat,dist,sig,dsigdt,dsigds,tofsig,sofsig
      real a0,a1,a2,cubr,cubq,cuban,cubrl,cubim,harmon,hfharm
c
      parameter (athird=1./3.)
      data c1,c2,c3,c4,c5,c6,c7/-7.2169e-2,4.9762e-2,8.0560e-1,
     .     -7.5911e-3,-3.0063e-3,3.5187e-5,3.7297e-5/
c
c --- formulae relating latitude to grid distance from equator
      alat(dist1,grid)=(2.*atan(exp(dist1*grid/radian))-pi/2.)
      dist(alat1,grid)=alog(tan((2.*alat1+pi)/4.))*radian/grid
c
c --- harmonic mean
      harmon(a,b)=2.*a*b/(a+b)
c
c --- harmonic mean divided by 2
      hfharm(a,b)=a*b/(a+b)
c
c --- -----------------
c --- equation of state
c --- -----------------
c
c --- sigma-theta as a function of temp (deg c) and salinity (mil)
c --- (friedrich-levitus 3rd degree polynomial fit)
c
      sig(t,s)=(c1+c3*s+t*(c2+c5*s+t*(c4+c7*s+c6*t)))*1.e-3
c
c --- d(sig)/dt
      dsigdt(t,s)=(c2+c5*s+2.*t*(c4+c7*s+1.5*c6*t))*1.e-3
c
c --- d(sig)/ds
      dsigds(t,s)=(c3+t*(c5+t*c7))*1.e-3
c
c --- auxiliary statements for finding root of 3rd degree polynomial
      a0(s)=(c1+c3*s)/c6
      a1(s)=(c2+c5*s)/c6
      a2(s)=(c4+c7*s)/c6
      cubq(s)=athird*a1(s)-(athird*a2(s))**2
      cubr(r,s)=athird*(.5*a1(s)*a2(s)-1.5*(a0(s)-1.e3*r/c6))
     .   -(athird*a2(s))**3
c --- if q**3+r**2>0, water is too dense to yield real root at given
c --- salinitiy. setting q**3+r**2=0 in that case is equivalent to
c --- lowering sigma until a double real root is obtained.
      cuban(r,s)=athird*atan2(sqrt(max(0.,
     .   -(cubq(s)**3+cubr(r,s)**2))),cubr(r,s))
      cubrl(r,s)=sqrt(-cubq(s))*cos(cuban(r,s))
      cubim(r,s)=sqrt(-cubq(s))*sin(cuban(r,s))
c
c --- temp (deg c) as a function of sigma and salinity (mil)
      tofsig(r,s)=-cubrl(r,s)+sqrt(3.)*cubim(r,s)-athird*a2(s)
c
c --- salinity (mil) as a function of sigma and temperature (deg c)
      sofsig(r,t)=(r*1.e3-c1-t*(c2+t*(c4+c6*t)))/(c3+t*(c5+c7*t))
c
c --- Note: changing the above-defined equation of state also requires 
c --- modification of statement functions embedded in -mxlayr.f-
c-----------------------------------------------------------------------------
