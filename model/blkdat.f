      SUBROUTINE init_block_data 
c
      include 'dimensions.h'
      include 'common_blocks.h'

c --- 'baclin' = baroclinic time step
c --- 'batrop' = barotropic time step
c --- 'diagfq' = number of days between model diagnostics (incl.output)
      data baclin,batrop/180.,3./,diagfq/.125/
c
c --- 'thkdff' = diffusion velocity (cm/s) for thickness diffusion
c --- 'veldff' = diffusion velocity (cm/s) for momentum dissipation
c --- 'temdff' = diffusion velocity (cm/s) for temp/salin. mixing
c --- 'viscos' is nondimensional, used in deformation-dependent viscosity
      data thkdff/0.5/,veldff/1.0/,temdff/0.0/,viscos/0.2/
c
c --- 'diapyc' = diapycnal diffusivity times buoyancy freq. (cm**2/s**2)
c --- 'vertmx' = diffusion velocity (cm/s) for mom.mixing across mix.layr.base
c --- 'mixfrq' = number of time steps between diapycnal mixing calculations
c --- 'h1'     = depth interval used in lateral weighting of hor.pres.grad.
c --- 'thkmin' = minimum mixed-layer thickness (m)
      data diapyc/0./,vertmx/0./,mixfrq/5/    
      data h1/980600./,thkmin/20./
c
c --- slip=+1  for free-slip boundary cond., slip=-1  for non-slip cond.
c --- 'cbar'   = rms flow speed (cm/s) for linear bottom friction law
c --- 'thkbot' = thickness of bottom boundary layer (pressure units)
c --- 'sigjmp' = minimum density jump at mixed-layer bottom (theta units)
      data slip/-1./,cbar/10./,thkbot/980600./,sigjmp/2.e-5/
c
c --- weights for time smoothing
      data wuv1,wuv2/.75,.125/
      data wts1,wts2/.875,.0625/
      data wbaro/.125/
c
c --- layer thicknesses in units of pressure:
      data tenm,onem,tencm,onecm,onemm/980600.,98060.,9806.,980.6,98.06/
      data radian/57.2957795/,pi/3.1415926536/
c
c --- 'g'      = gravitational acceleration
c --- 'csubp'  = specific heat of air at constant pressure (j/g/deg)
c --- 'spcifh' = specific heat of sea water (j/g/deg)
c --- 'cd'     = drag coefficient
c --- 'ct'     = thermal transfer coefficient
c --- 'airdns' = air density at sea level (g/cm**3)
c --- 'evaplh' = latent heat of evaporation (j/g)
c --- 'thref'  = reference value of specific volume (cm**3/g)
c --- 'epsil'  = small nonzero number used to prevent division by zero
      data g/980.6/,csubp/1.0057/,spcifh/3.99/,cd/.0013/,ct/.0012/
      data airdns/1.2e-3/,evaplh/2.47e3/,thref/1./,epsil/1.e-11/
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
c
c --- the original lat/lon grid points are specified in terms of
c ---      nlato,nlongo  --  grid dimension in lat/lon (i,j) direction,
c ---      grido         --  grid size (degrees),
c ---      xpivo,ypivo   --  the i,j index of the equatorial pivot point.
c
c --- the grid points in the rotated mercator grid are specified in terms of
c ---      nlatn,nlongn  --  grid dimension in (rotated) i,j direction,
c ---      gridn         --  grid size in degrees longitude,
c ---      xpivn,ypivn   --  the i,j index of the equatorial pivot point.
c
      data xpivn,ypivn,gridn/1200. ,0., 0.0375/
c
c ---  s w i t c h e s    (if set to .true., then...)
c --- rotat       rotate poles in mercator projection by 90 deg.
c --- thermo      use thermodynamic forcing functions
c --- windf       use wind stress forcing function
c --- relax       activate lateral boundary nudging
c
      data thermo/.false./,rotat/.false./,windf/.false./,relax/.false./
c
c --- 'lp' = logical unit number for printer output
      data lp/6/
c
c --- use 'huge' to initialize array portions that the code should never access
      data huge/1.e33/
c
c --- i/o file names
c
      data flnmrsi/ 'restart_in'/
      data flnmrso/ 'restart_out'/
      data flnmrss/ 'restar6_out'/
      data flnmrse/ 'restar8_out'/
      data thbase/.0/

c --- layer densities (sigma units):
      bin_restart = .true.
      OPEN (99, file = "thetas")
         READ (99,*) theta
      CLOSE (99)
     
       sigma = 1000.0 *theta
 
      end
