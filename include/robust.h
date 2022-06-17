c------------------------------------------------------------------
c contains the common blocks for data for robust runs
c------------------------------------------------------------------

#ifdef robust
      real rtracer(imt,jmt,km,12),gammar(km)
#ifdef seasonal
      real rtracersea(imt,jmt,kmsea,2,12) 
#endif
#endif
#if defined ptisotoperobust || defined ptisotopepointrobust
      real d18orobust(imt,jmt,km,nt-2)
#endif
      common /robustdata/ 
#ifdef robust
     & rtracer,gammar
#ifdef seasonal
     & ,rtracersea
#endif
#endif
#if defined ptisotoperobust || defined ptisotopepointrobust
     & d18orobust
#endif

