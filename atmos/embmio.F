      subroutine initatmos
      implicit none

c=======================================================================
c     initialize variables for the energy-moisture balance model
c     based on code by a.fanning
c
c     author:   m.eby   e-mail: eby@uvic.ca
c=======================================================================

# include "param.h"
# include "slap.h"
# include "atm.h"
# ifdef uvic_embm_ice
#  include "ice.h"
# endif
# include "cembm.h"
# include "csbc.h"
# include "levinda.h"
# include "mapsbc.h"
# include "switch.h"
# include "ctmngr.h"

      real val
      integer index

      character*10 field

c-----------------------------------------------------------------------
c     set initial or default restart conditions
c-----------------------------------------------------------------------
# if !defined atmos_start
      if(itt.eq.1)then
      nats = namix
      irstdy = 0
      msrsdy = 0
      do j=1,jmta
        do i=1,imta
          do n=1,nat
            val = 0.0
            if (n .eq. 1) val = 10.0
            if (n .eq. 2) val = 0.01
            at(i,j,1,n) = val
            at(i,j,2,n) = val
	  enddo
# ifdef uvic_embm_ice
	  hice(i,j,1) = 0.0
	  hice(i,j,2) = 0.0
	  aice(i,j,1) = 0.0
	  aice(i,j,2) = 0.0
#endif
   	enddo
      enddo
      endif
# else
      nats = namix
      irstdy = 0
      msrsdy = 0
      do j=1,jmta
        do i=1,imta
          do n=1,nat
            val = 0.0
            if (n .eq. 1) val = 10.0
            if (n .eq. 2) val = 0.01
            at(i,j,1,n) = val
            at(i,j,2,n) = val
	  enddo
# ifdef uvic_embm_ice
	  hice(i,j,1) = 0.0
	  hice(i,j,2) = 0.0
	  aice(i,j,1) = 0.0
	  aice(i,j,2) = 0.0
#endif
   	enddo
      enddo
#endif

c-----------------------------------------------------------------------
c     zero initial flux accumulators
c-----------------------------------------------------------------------

      do n=1,nat
        do j=1,jmta
          do i=1,imta
            flux(i,j,n) = 0.0
          enddo
        enddo
      enddo          

c-----------------------------------------------------------------------
c     set initial approximation for vector of unknowns
c-----------------------------------------------------------------------

      do n=1,nat
	do index=1,orda
	  slaplhs(index,n) = at(gridcol(index),gridrow(index),2,n)
	enddo
 	call embmbc (at(1,1,2,n))
      enddo

      return
      end 



