c====================== include file "config.h" ====================== 
#ifdef hcomments 
c 
c     declarations of variables used by the pre processing routines
c     to generate appropriate sub-files for the main SEA program
c--------------------------------------------------------------------- 
#endif 

#include "domain.h"
#include "pconst.h"

      integer    imt, jmt, imtm1, jmtm1, imtp1, jmtp1, imtm2
      integer    imtnb, jmtnb
      parameter (imt=iall, jmt=jall, imtm1=imt-1, jmtm1=jmt-1)
      parameter (imtp1=imt+1, jmtp1=jmt+1, imtm2=imt-2)
      parameter (imtnb=imt-2, jmtnb=jmt-2)

      integer    max_ewsg, max_nssg, ewsg, nssg
      parameter (max_ewsg=iall/4, max_nssg=jall/4)

      integer    kmt(0:imtp1,0:jmtp1), kmu(0:imtp1,0:jmtp1)
      integer    sd_kmt(0:imtp1,0:jmtp1), sd_kmu(0:imtp1,0:jmtp1)
      common    /levind/ 
     &           kmt, kmu, sd_kmt, sd_kmu

      integer    map(0:max_ewsg, 0:max_nssg)
      integer    imtsd(0:max_ewsg, 0:max_nssg)
      integer    jmtsd(0:max_ewsg, 0:max_nssg)
      integer    ioffsd(0:max_ewsg, 0:max_nssg)
      integer    joffsd(0:max_ewsg, 0:max_nssg)
      common    /subdomains/
     &           map, imtsd, jmtsd, ioffsd, joffsd, ewsg, nssg

      real       dx, dxr, dx2r, dx4r, dy, dyr, dy2r, dy4r
      real       csu(jmt), csur(jmt), cst(jmt), cstr(jmt)
      real       phi(jmt), phit(jmt), sine(jmt), tng(jmt)
      real       fcor(jmt), tanra(jmt)
      real       dz(km), dzr(km), dz2r(km), c2dz(km)
      real       dzw(0:km), dzwr(0:km), dzw2r(0:km)
      common    /grdvar/ 
     &           dx, dxr, dx2r, dx4r, dy, dyr, dy2r, dy4r,
     &           csu, csur, cst, cstr, phi, phit, sine, tng,
     &           fcor, tanra, 
     &           dz, dzr, dz2r, c2dz, dzw, dzwr, dzw2r
