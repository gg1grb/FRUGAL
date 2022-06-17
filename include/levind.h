c===================== include file "levind.h" =======================
#ifdef hcomments
c
c     vertical level indicators which define model geometry & 
c     bottom topography:
c
c     kmt = number of vertical boxes over "t" points
c     kmu = number of vertical boxes over "u,v" points
c---------------------------------------------------------------------
#endif

      integer kmt(0:imax+1,0:jmax+1), kmu(0:imax+1,0:jmax+1)
#ifdef obc
     &   ,kmtopen(imt,jmt)
#endif
      common  /levind/
     &        kmt, kmu
#ifdef obc
     &   ,kmtopen
#endif

