c-----------------------------------------------------------------------------
      c o m m o n  /micom1/
     . u(idm,jdm,2*kdm),v(idm,jdm,2*kdm)      ! velocity components
     .,dp(idm,jdm,2*kdm),dpold(idm,jdm,kdm)   ! layer thickness
     .,dpu(idm,jdm,2*kdm),dpv(idm,jdm,2*kdm)  ! layer thickness at u,v points
     .,delp5(idm,jdm,2*kdm)
     .,p(idm,jdm,kdm+1),dpini(idm,jdm,kdm)     ! interface pressure
     .,corio(idm,jdm),dp5(idm,jdm,kdm)
     .,psikk(idm,jdm)
     .,potvor(idm,jdm,kdm)                        ! potential vorticity
     .,dpm(idm,jdm,kdm)
     .,ssh(idm,jdm),ssm(idm,jdm)
     .,pre(idm,jdm,kdm)
     .,pr5(idm,jdm,kdm)
     .,vp1(jdm),vp5(jdm)
     .,pin(idm,jdm,kdm+1)
c
      real u,v,dp,dpold,dpu,dpv,p,corio,psikk,delp5
     .,potvor,dpm,ssh,pre,dpini,pin,vp1,ssm
     .,dp5,pr5,vp5
c
      c o m m o n  /micom2/
     . montg(idm,jdm,kdm)  ! montgomery potential
     .,ubavg(idm,jdm,3),vbavg(idm,jdm,3)      ! barotropic velocity
     .,pbavg(idm,jdm,3)                       ! barotropic pressure
     .,ubrhs(idm,jdm),vbrhs(idm,jdm)          ! rhs of barotropic u,v eqns.
     .,utotm(idm,jdm),vtotm(idm,jdm)          ! total (barotrop.+baroclin.)..
     .,utotn(idm,jdm),vtotn(idm,jdm)          ! ..velocities at 2 time levels
     .,uflux(idm,jdm),vflux(idm,jdm)          ! horizontal mass fluxes
     .,uflux2(idm,jdm),vflux2(idm,jdm)        ! more mass fluxes
     .,vpfl1(jdm),vpfl5(jdm)
c
      real montg,ubavg,vbavg,pbavg,ubrhs,vbrhs,utotm,
     .     vtotm,utotn,vtotn,uflux,vflux,uflux2,vflux2,
     .     vpfl1,vpfl5
c
      c o m m o n  /micom3/
     . util1(idm,jdm),util2(idm,jdm)          ! arrays for temporary storage
     .,util3(idm,jdm)          ! arrays for temporary storage
     .,scuv(idm,jdm)            ! mesh size at u pts in x,y dir.
     .,scu2(idm,jdm)            ! grid box size at u,v pts
     .,scui(idm,jdm)          ! inverses of scux,scvy
     .,sc2i(idm,jdm)          ! inverses of scp2,scq2
     .,gradx(idm,jdm),grady(idm,jdm)          ! horiz. presssure gradient
     .,depthu(idm,jdm),depthv(idm,jdm)        ! bottom pres. at u,v points
     .,pvtrop(idm,jdm)                        ! pot.vort. of barotropic flow
     .,drag(idm,jdm),dragb(idm,jdm)           ! bottom drag
     .,pbot(idm,jdm)                          ! bottom pressure at t=0
c
      real util1,util2,util3,scuv,scui,
     .     scu2,sc2i,gradx,pbot,
     .     grady,depthu,depthv,pvtrop,drag,dragb
c  
c ---  s w i t c h e s    (if set to .true., then...)
c --- diagno      output model fields and diagnostic messages
c --- rotat       rotate poles in mercator projection by 90 deg.
c --- thermo      use thermodynamic forcing functions
c --- windf       include wind stress in forcing functions
c --- relax       activate lateral boundary nudging
c --- trcrin      initialize tracer from restart file
c --- trcout      advect tracer and save results in history/restart file
c --- realism     include realistic time dependent forcing
c --- forc        includes buoyancy forcing
c --- bdrag       includes bottom drag
c --- thc
c
      logical diagno,rotat,thermo,windf,relax,trcrin,trcout,realism
     .,forc,bdrag,thc
      common/swtchs/diagno,rotat,thermo,windf,relax,trcrin,trcout
     .,realism,forc,bdrag,thc
c
      common/varbls/time,delt1,dlt,w0,w1,w2,w3,ws0,ws1,ws2,ws3,
     .      area,watcum,empcum,mntin,
     .      nstep,nstep1,nstep2,lstep,l0,l1,l2,l3,ls0,ls1,ls2,ls3
c
      real time,delt1,dlt,w0,w1,w2,w3,ws0,ws1,ws2,ws3,
     .     area,watcum,empcum,mntin
      integer nstep,nstep1,nstep2,lstep,l0,l1,l2,l3,ls0,ls1,ls2,ls3
c
c --- 'baclin' = baroclinic time step
c --- 'batrop' = barotropic time step
c --- 'thkdff' = diffusion velocity (cm/s) for thickness diffusion
c --- 'veldff' = diffusion velocity (cm/s) for momentum dissipation
c --- 'temdff' = diffusion velocity (cm/s) for temp/salin. mixing
c --- 'viscos' is nondimensional, used in deformation-dependent viscosity
c --- 'diapyc' = diapycnal diffusivity times buoyancy freq. (cm**2/s**2)
c --- 'vertmx' = diffusion velocity (cm/s) for mom.mixing across mix.layr.base
c --- 'mixfrq' = number of time steps between diapycnal mixing calculations
c --- 'h1'     = depth interval used in lateral weighting of hor.pres.grad.
c --- slip = +1  for free-slip boundary cond., slip = -1  for non-slip cond.
c --- 'cbar'   = rms flow speed (cm/s) for linear bottom friction law
c --- 'diagfq' = number of days between model diagnostics (incl.output)
c --- 'ntracr' = number of time steps between tracer transport
c --- 'wuv1/2' = weights for time smoothing of u,v field
c --- 'wts1/2' = weights for time smoothing of t,s field
c --- 'wbaro'  = weight for time smoothing of barotropic u,v,p field
c --- 'thkmin' = minimum mixed-layer thickness (m)
c --- 'thkbot' = thickness of bottom boundary layer (pressure units)
c --- 'sigjmp' = minimum density jump at mixed-layer bottom
c ---' salmin' = minimum salinity allowed in an isopycnic layer
c
      common/parms1/sigma(kdm),theta(kdm),thbase,baclin,batrop,thkdff,
     .              veldff,temdff,viscos,diapyc,vertmx,h1,slip,cbar,
     .              diagfq,wuv1,wuv2,wts1,wts2,wbaro,thkmin,
     .              thkbot,sigjmp,salmin(kdm),mixfrq,ntracr
c
      real sigma,theta,thbase,baclin,batrop,thkdff,veldff,temdff,viscos,
     .     diapyc,vertmx,h1,slip,cbar,diagfq,wuv1,wuv2,wts1,wts2,
     .     wbaro,thkmin,thkbot,sigjmp,salmin
      integer mixfrq,ntracr
c
c --- 'tenm,onem,...' = pressure thickness values corresponding to 10m,1m,...
c --- 'g'      = gravity acceleration
c --- 'csubp'  = specific heat of air at constant pressure (j/g/deg)
c --- 'spcifh' = specific heat of sea water (j/g/deg)
c --- 'cd'     = drag coefficient
c --- 'ct'     = thermal transfer coefficient
c --- 'airdns' = air density at sea level (g/cm**3)
c --- 'evaplh' = latent heat of evaporation (j/g)
c --- 'thref'  = reference value of specific volume (cm**3/g)
c --- 'epsil'  = small nonzero number used to prevent division by zero
c
      common/consts/tenm,onem,tencm,onecm,onemm,g,csubp,spcifh,cd,ct,
     .              airdns,evaplh,thref,epsil,huge,radian,pi
c
      real tenm,onem,tencm,onecm,onemm,g,csubp,spcifh,cd,ct,airdns,
     .     evaplh,thref,epsil,huge,radian,pi
c
c --- grid point where detailed diagnostics are desired:
      common/testpt/itest,jtest
c
      integer itest,jtest
c
c --- two   m a p   p r o j e c t i o n s   are supported:
c --- rotat = .false.: conventional mercator projection
c ---                  (i and x pointing south, j and y pointing east)
c --- rotat = .true.:  mercator projection with poles rotated +90 deg.
c ---                  (i and x pointing east, j and y pointing north)
c
c --- if rotat=.true., (xpivo,ypivo) and (xpivn,ypivn), respectively, mark
c --- the location of the intersection of the true and model equator ('pivot
c --- point') in the lat/lon and rotated model grid.
c --- if rotat=.false., xpivn specifies the i-index of the equator (xpivo,
c --- ypivo,ypivn not used).
c --- grido = mesh size of lat/lon grid in degrees
c --- gridn = mesh size of actual model grid (rotated or not) in deg. longitude
c
      common/pivot/xpivo,ypivo,grido,xpivn,ypivn,gridn
c
      real xpivo,ypivo,grido,xpivn,ypivn,gridn
c
c --- locations of lat/lon grid points and of associated x-coordinate
c --- direction as viewed in rotated model grid:
c
      common/newgrd/ xnew(nlato,nlongo), ynew(nlato,nlongo),
     .              dxnew(nlato,nlongo),dynew(nlato,nlongo)
c
      real xnew, ynew,dxnew,dynew
c
c --- locations of rotated model grid points and of associated x-coordinate 
c --- direction as viewed in lat/lon grid:
c
      common/oldgrd/ xold(nlatn,nlongn), yold(nlatn,nlongn),
     .              dxold(nlatn,nlongn),dyold(nlatn,nlongn)
c
      real xold,yold,dxold,dyold
c
      character*48 flnmdep,flnmrsi,flnmrso,flnmarc
     .,flnmrss,flnmrse
c
      common/iovars/flnmdep,flnmrsi,flnmrso,flnmarc
     .,flnmrss,flnmrse
c
c
c> Revision history:
c>
c> July 1997 - eliminated 3-D arrays -uold,vold- (used in time smoothing)
c-----------------------------------------------------------------------------
