c===================== include file "scalar.h" =======================
#ifdef hcomments
c
c  various scalar quantities:
c
c     dtts   = time step for density & tracers (in seconds)
c     dtuv   = time step for baroclinic velocity (in seconds)
c     dtbt   = time step for barotropic velocity (in seconds)
c              used by free surface model
c     c2dtts = 2*dtts
c     c2dtuv = 2*dtuv
c     c2dtbt = 2*dtbt
c     area   = surface area of ocean (cm**2)
c     volume = volume of ocean (cm**3)
c     omega  = earth's rotation rate (radians/sec)
c     radius = earth's radius (cm)
c     grav   = earth's gravitational acceleration (cm/sec**2)
c     cdbot  = bottom drag coefficient
c
c     ncon   = number of  passes through convective code in tracer
c     ntbt   = number of free surface time steps per baroclinic 
c              velocity time steps.
c
c  non dimensional quantities:
c
c     radian = degrees per radian
c     pi     = something good to eat (or something else)
c---------------------------------------------------------------------
#endif

#ifdef dttsacc
      real dtacc(km)
#endif
      real    dtts, dtuv, dtbt, c2dtts, c2dtuv, c2dtbt
      real    area, volume, omega, radius, grav, cdbot,
     &        seglength
      integer ncon, ntbt, ntbt2, ntbtp0

      common  /scalar/
     &        dtts, dtuv, dtbt, c2dtts, c2dtuv, c2dtbt,
     &        area, volume, omega, radius, grav, cdbot,
     &        ncon, ntbt, ntbt2, ntbtp0,seglength
#ifdef dttsacc
     &  ,dtacc
#endif
      real    radian, pi,rho0,rho0r,rho_zero,oceanheat,oceanheat1
      common  /ndcon/
     &        radian, pi,rho0,rho0r,rho_zero,oceanheat,oceanheat1
