c===================== include file "frees.h" ========================
#ifdef hcomments
c
c     variables for free surface model calculation of external mode
c
c     h0     =  free surface model sea surface height
c     u0     =  free surface model velocity east
c     v0     =  free surface model velocity north
c     h      = depth over "u,v" points
c     hr     = reciprocal depth over "u,v" points
c     zu     = vertically averaged zonal forcing 
c     zv     = vertically averaged meridional forcing 
#ifndef free_eb
c     freeav = time average of free surface fields
#endif
c---------------------------------------------------------------------
#endif

      real    h0(fimin:fimax, fjmin:fjmax, 3)
      real    u0(fimin:fimax, fjmin:fjmax, 3)
      real    v0(fimin:fimax, fjmin:fjmax, 3)
      real    h (fimin:fimax, fjmin:fjmax)
      real    hr(fimin:fimax, fjmin:fjmax)
      real    zu(fimin:fimax, fjmin:fjmax)
      real    zv(fimin:fimax, fjmin:fjmax)
#ifndef free_eb
      real    freeav(3, fimin:fimax, fjmin:fjmax)
#endif
#ifdef de_checkbd
#ifdef varts
      real    dchkbd(imt,jmt),dchkbdx
#else
      real    dchkbd
#endif
#endif

      common  /freesurf/
     &        h0, u0, v0, h, hr, zu, zv
#ifndef free_eb
     &,       freeav
#endif
#ifdef de_checkbd
#ifdef varts
     &,       dchkbd,dchkbdx
#else
     &,       dchkbd
#endif
#endif

