
      subroutine sumflux
c
c=======================================================================
c     sum fluxes over atmospheric time steps
c     based on code by a.fanning
c
c     author:   m.eby   e-mail: eby@uvic.ca
c=======================================================================
c
      implicit none

#include "param.h"
#include "cembm.h"
#include "atm.h"
#include "levinda.h"
#ifdef mrw_rivers
#include "riv.h"
#include "grdvara.h"
#endif

      real dta

      dta = dtatm


      do j = 2,jmtam1
        do i = 2,imtam1
          if (kmta(i,j) .eq. 0) then

            flux(i,j,1) = flux(i,j,1) + dts*(dnswr(i,j) - uplwr(i,j) 
     &                  - upltnt(i,j) - upsens(i,j))
 
           endif

	    flux(i,j,2) = flux(i,j,2) + dts*(precip(i,j) - evap(i,j))



#ifdef mrw_rivers
       flux(i,j,2)=flux(i,j,2)+runoff(i,j)/
     &    (dxta(i,j)*dyta(i,j))
#endif
# ifdef uvic_embm_seasonal_average
          atbar(i,j) = atbar(i,j) + dta*at(i,j,2,1)
          aqbar(i,j) = aqbar(i,j) + dta*at(i,j,2,2)
#endif

	enddo
      enddo
c
# ifdef uvic_embm_seasonal
      totaltime = totaltime + dta
# endif
      avetime = avetime + dts
c
      return
      end




      subroutine aveflux
      implicit none
c
c=======================================================================
c     calculate the average fluxes for next ocean time step
c     based on code by a.fanning
c
c     author:   m.eby   e-mail: eby@uvic.ca
c=======================================================================
c
#include "param.h"
#include "cembm.h"
#include "atm.h"
#include "levinda.h"
#include "csbc.h"
#include "mapsbc.h" 
#include "iounit.h"
#include "ctmngr.h"
C
      real fa,fb
      integer num

c
c-----------------------------------------------------------------------
c     calculate average net fluxes. convert heat flux to cal/cm**2/s 
c     from kW/m**2 and fresh water flux (cm/s) to an apparent salt 
c     flux (g/cm**2/s) using global ocean average salinity, socn
c-----------------------------------------------------------------------
c
c put heat flux into Wm-2
c      fa = 2.389e-8/avetime
       fa=1.e-3/avetime
c put fresh water flux into mm/day
      fb = 864000./avetime
c      fb = -socn*rhoocn/avetime
      num = 0
      do j=2,jmtam1
        do i=2,imtam1
	  if (kmta(i,j) .eq. 0) then
# ifdef uvic_embm_ice
            sbcocn(i,j,ihflx) = fa*flux(i,j,1)
# else
c           check we are not taking more heat than in first layer
c           (for use without an ice model, to prevent freezing sst)
            sbcocn(i,j,ihflx) = max(fa*flux(i,j,1), sbcocn(i,j,ihflx))
# endif
c was ipme not ifflx
            sbcocn(i,j,ifflx) = fb*flux(i,j,2)
          endif
	enddo
      enddo	
c
      return
      end

