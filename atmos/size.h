c====================== include file "size.h" ==========================
c                                                                       
c-----------------------------------------------------------------------
c     USER INPUT: 
c-----------------------------------------------------------------------
c
c     imta    = number of grid points in the longitudinal direction
c              (calculated points are from 2 through imt-1. End points
c               are boundaries)          
c
c     jmta    = number of grid points (latitude rows) in the latitudinal
c              direction (calculated points are from 2 through jmt-1.           
c              End points are boundaries)   
c                                
c     km     = number of grid points in the vertical direction          
c              (calculated points are from 1 through km)
c
c     nt     = number of tracers (temperature, salinity, ...)
c
c     mnisle = maximum number of islands (unconnected land masses)
c
c     maxipp = maximum number of all island perimeter points
c-----------------------------------------------------------------------
c                                                                       
#ifdef vvcoarse
      parameter ( imta  = 22, jmta  = 23 )
#endif
#ifdef vcoarse
      parameter ( imta  = 47, jmta  = 52 )
#endif
#ifdef vcoarse_latlon
      parameter ( imta  = 62, jmta  = 40 )
#endif
#ifdef coarse
      parameter ( imta  = 92, jmta  = 105 )
#endif
#ifdef coarse_latlon
      parameter ( imta  = 122, jmta  = 80 )
#endif
#ifdef fine
      parameter ( imta  = 182, jmta  = 211 )
#endif
#ifdef fine_latlon
      parameter ( imta  = 242, jmta  = 160 )
#endif

      parameter (nat=2)
