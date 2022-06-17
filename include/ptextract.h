c----------------------------'ptextract.h"--------------------
c  contains variables associated with extracting the 
c  circulation for passive tracer runs
c-------------------------------------------------------------

      real smonth
#ifdef daily
     & ,sday
#endif
      integer msaheadx,msahead
      integer conv(imt,jmt,km)

      common /ptex/smonth,msaheadx,msahead,conv
#ifdef daily
     & ,sday
#endif
      
