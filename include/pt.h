c------------------------"pt.h"-----------------------------
c  contains the arrays to store data for passive tracer runs
c-----------------------------------------------------------

#ifdef ptcore
      real upt(imt,jmt,km,12),vpt(imt,jmt,km,12),
     $     tptx(imt,jmt,km,2,12)
      integer convpt(imt,jmt,km,12)
#endif
#ifdef ptmonthly
      real upt(imt,jmt,km,2),vpt(imt,jmt,km,2),
     $ tptx(imt,jmt,km,2)
#endif
      real spt(imt,jmt,nt),tpt(km,imt,jmt,2)
      integer conv(imt,jmt,km)
      common /ptx/ upt,vpt,tptx,tpt,conv,spt
#ifdef ptcore
     $, convpt
#endif
