c============================== domain.h =============================
#ifdef hcomments
c---------------------------------------------------------------------
c     stlon = longitudinal start reference of domain
c     stlat = latitudinal start reference of domain
c     dxdeg = horizontal resolution in longitude direction
c     dydeg = horizontal resolution in latitude direction
c     iall  = number of longitudinal grid points in whole domain
c             (including one outer halo for cyclic/closed boundaries)
c     jall  = number of latitudinal grid points in whole domain
c             (including one outer halo for closed boundaries)
c     km    = number of vertical depth levels in model domain
c     nt    = number of tracers (temperature, salinity, etc.)
c     maxp  = maximum number of processors likely to be used
c---------------------------------------------------------------------
c     For a user defined model set all parameters as required.
c---------------------------------------------------------------------
#endif

      real      stlon, stlat, dxdeg, dydeg
      integer   iall, jall, km, nt, maxp
      integer   ntx
c set nt here to avoid having to set for each model resolution
#if defined ptcore || defined ptmonthly
      parameter(ntx=1)
#endif
#if defined ptisotope || defined ptisotopesalinity
      parameter(ntx=3)
#else
      parameter(ntx=2)
#endif
      parameter ( km = 19, nt = ntx, maxp = 1 )
#ifdef vvcoarse
      parameter ( iall  = 22, jall  = 23 )
#endif
#ifdef vcoarse
      parameter ( iall  = 47, jall  = 52 )
#endif
#ifdef coarse
      parameter ( iall  = 92, jall  = 105 )
#endif
#ifdef coarse_latlon
      parameter ( iall  = 122, jall  = 80 )
#endif
#ifdef fine
      parameter ( iall  = 182, jall  = 211 )
#endif
#ifdef vcoarse_latlon
      parameter ( iall  = 62, jall  = 40 )
#endif
#ifdef fine_latlon
      parameter ( iall  = 242, jall  = 160 )
#endif
