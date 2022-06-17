c====================== include file "param.h" =========================
c
c     main parameter file which sets ocean characteristics:
c
c     "size.h" sets model size "imta", "jmta"
c
#include "size.h"
c
c     nvar   = number of prognostic variables
c     lseg   = maximum number of longitudinal stream function segments
c     nlatpr = maximum number of latitudes for matrix printouts
c              on diagnostic time steps
c     nhreg  = number of regions in the horizontal used for averaging
c              tracers.
c     nvreg  = number of regions in the vertical used for term balance
c              calculations. note "nvreg" isn`t used for tracer averages
c     numreg = total number of regions ( = product of nhreg & nvreg)
c              used for term balance calculations
c
      parameter (nhreg=5, nvreg=2, numreg=nhreg*nvreg)
      parameter (nvar=nt+2)
c
c
      parameter (imtam1=imta-1)                                
      parameter (imtap1=imta+1, imtam2=imta-2
     &,          jmtap1=jmta+1, jmtam1=jmta-1, jmtam2=jmta-2
     &,          jscana=jmtam2)
c
#ifdef mrw_rivers
      real dxyscale,dtriver
      integer imtr,jmtr,nsega,nlandp
      parameter(dxyscale=12,imtr=(imta-2)*12+2,jmtr=jmta*12,
     &          nsega=5,nlandp=266760)
#endif
