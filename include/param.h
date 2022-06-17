c===================== include file "param.h" ========================
#ifdef hcomments
c
c  positional indicies used by the majority of subroutines
c     i, j, k, l, m, n
c
c  processor identification variables
c     npe    = number of processors (always 1 if sequential mode)
c     me     = processor id (always 0 if sequential mode)
c     master = processor 0 (always true if sequential mode)
c
c  ocean domain characteristics:
c     nt     = number of tracers
c     imt    = number of grid points in longitudinal direction
c              this varies per sub domain ... imax-1 <= imt <= imax
c     jmt    = number of grid points in latitudinal direction
c              this varies per sub domain ... jmax-1 <= jmt <= jmax
c     ioff   = longitudinal offset of a sub domain
c     joff   = latitudinal offset of a sub domain
c---------------------------------------------------------------------
#endif

#include "domain.h"
      integer   imax, jmax, bmax
      parameter (imax=  iall, jmax=  jall)
      integer   ewsg, nssg
      parameter (ewsg=   1, nssg=   1)
      integer    i, j, k, l, m, n

      integer    npe, me, pcyc
      logical    master
      common    /idtag/ 
     &           npe, me, pcyc, master

      integer    imt, jmt, ioff, joff
      integer    imtp1, imtm1, imtm2, imtm3
      integer    jmtp1, jmtm1, jmtm2, jmtm3
      integer    imu, imum1, imum2

      integer   fimin, fimax, fjmin, fjmax
#ifdef seasonal
      integer kmsea
#endif
#if defined coupled || defined iceonly
      integer imta,jmta,nat
      integer nhreg,nvreg,numreg,nvar
      integer imtam1,imtap1,jmtam1,jmtap1,jscana
      integer imtam2,jmtam2
c      integer maxsbc,numsbc,numosbc,numasbc
#endif
#ifdef varts
c  subtsmax is the maximum number of sub-time steps
c  pps - advective and diffusive terms require different time levels
c        where they cross a time step boundary. the arrays uea,uwa,...,
c        tea,twa..... contain the time levels to be used when calculating
c        these terms. the size of these arrays is (-pps:pps,-pps:pps)
c        pps=1 for normal advection scheme, but more for quick scheme
      integer   subtsmax,pps
#endif
      parameter (fimin=1, fimax=imax, fjmin=1, fjmax=jmax)

      parameter (imt=imax, jmt=jmax, ioff=0, joff=0)
      parameter (imtp1=imt+1, imtm1=imt-1, imtm2=imt-2, 
     &           jmtp1=jmt+1, jmtm1=jmt-1, jmtm2=jmt-2, 
     &           imu=imt, imum1=imu-1, imum2=imu-2)
#ifdef varts
#ifdef vvcoarse
      parameter(subtsmax=4,pps=1)
#endif
#ifdef vcoarse
      parameter(subtsmax=8,pps=1)
#endif
#ifdef coarse
      parameter(subtsmax=8,pps=1)
#endif
#ifdef fine 
      parameter(subtsmax=128,pps=1)
#endif
#ifdef vcoarse_latlon
      parameter(subtsmax=4,pps=1)
#endif
#ifdef coarse_latlon
      parameter(subtsmax=8,pps=1)
#endif
#ifdef fine_latlon
      parameter(subtsmax=128,pps=1)
#endif
#endif
#ifdef seasonal
#ifdef level19
      parameter(kmsea=10)
#endif
#ifdef level36
      parameter(kmsea=17)
#endif
#endif
#include "pconst.h"

#if defined coupled || defined iceonly
#include "parama.h"
#endif
c ds is the no. of degrees west of 0E of the centre of the grid.
c must be between 0 and 360
      integer ds
      parameter(ds=40)
#ifdef adjustrelax
      integer yearadjustrelaxstart
      parameter(yearadjustrelaxstart=10.)
#endif
