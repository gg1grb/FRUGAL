      subroutine fluxes
      implicit none

c=======================================================================
c     calculate energy and moisture fluxes
c     based on code by a.fanning
c
c     Note: evaporation and precipitation are in cm/s
c           and humidities are in g/g
c
c     for Thompson and Warren outgoing radiation (see: Thommpson S.J.,
c     and S.G. Warren 'parameterization of outgoing ...'J. Atmos. Sci.,
c     39, 2667-2680, 1982.
c
c     for cloud parameterization (1) see:
c     Berliand, M.E. and T.G. Berliand, 'determining the net ...', 
c     Isv. Akad. Nauk. SSSR Ser. Geofiz No. 1, 1952.
c     or
c     Esbensen S.K., and Y. Kushnir, 'the heat budget of the global ...'
c     Clim. Res. Inst., Rep 29., Dep. Atmos. Sci., Oregon State Univ.,
c     Oregon, 1981.
c
c     author:   m.eby   e-mail: eby@uvic.ca
c=======================================================================

# include "param.h"
# include "cembm.h"
# include "atm.h"
# include "levinda.h"
# include "csbc.h"
#ifdef landevap
# include "riv.h"
# include "grdvara.h"
#endif
#ifdef cmeans3
#include "cmipmean3.h"
#include "fields.h"
#include "ctmngr.h"
#endif
#ifdef enhancedso
#include "coorda.h"
#endif
      real epsln,fa,fb,fc,fd,fe,b00,b10,b20,b01,b11,b21,b02,b12,b22,
     & b03,b13,b23,dswr,rh,rhrh,dalt,rp,dt,ws,cf,temp,
     & storage,spdevap,aquifer
      integer ii,jj

      epsln=1.e-20

c-----------------------------------------------------------------------
c     set appropriate constants
c-----------------------------------------------------------------------

      fa = 3.80116e03/ssp
      fb = rhoatm/rhoocn
      fc = vlocn*rhoocn
      fd = 0.94*rhoatm*cpatm

# ifndef uvic_embm_planetary
c----------------------------------------------------------------------
c     Set Thomson and Warren constants
c----------------------------------------------------------------------

      b00 = 2.43414e2
      b10 = -3.47968e1
      b20 = 1.02790e1
      b01 = 2.60065
      b11 = -1.62064
      b21 = 6.34856e-1
      b02 = 4.40272e-3
      b12 = -2.26092e-2
      b22 = 1.12265e-2
      b03 = -2.05237e-5
      b13 = -9.67e-5
      b23 = 5.62925e-5

# endif

      do j=2,jmtam1

c----------------------------------------------------------------------
c  calculate surface solar flux
c----------------------------------------------------------------------

        do i=2,imtam1
        dswr = (1.0 - scatter(i,j))*solins(i,j)*coalbedo(i,j)
c-----------------------------------------------------------------------
c         calculate saturation specific humidity at surface air temp
c-----------------------------------------------------------------------
 
# ifdef uvic_embm_mountains
          temp = at(i,j,2,1)-elev(i,j)/14285.
          ssh(i,j) = fa*exp(17.67*temp/(temp + 243.5))
# else
          ssh(i,j) = fa*exp(17.67*at(i,j,2,1)/(at(i,j,2,1) + 243.5))
# endif
# if defined uvic_embm_clouds || !defined uvic_embm_planetary
          rh = max(0.0, min(1.0, at(i,j,2,2)/(ssh(i,j) + epsln)))
# endif
#ifdef cmeans3
      if(itt.ge.cmipoutput(3,1).and.itt.le.cmipoutput(3,2))then
      a2d3(i,j,9)=a2d3(i,j,9)+rh
      endif
#endif
 
c-----------------------------------------------------------------------
c         calculate outgoing longwave radiation
c-----------------------------------------------------------------------
 
# ifdef uvic_embm_planetary
          outlwr(i,j) = esplnt(i,j)*(at(i,j,2,1) + 273.15)**4
# else
          rhrh = rh*rh
          outlwr(i,j) = 1.0e3*(b00 + b10*rh + b20*rhrh
     &                + (b01 + b11*rh + b21*rhrh)*at(i,j,2,1)
     &                + (b02 + b12*rh + b22*rhrh)*at(i,j,2,1)**2
     &                + (b03 + b13*rh + b23*rhrh)*at(i,j,2,1)**3)
#endif
# if defined uvic_embm_CO2 || defined lgm
     &                - anthro
# endif

c-----------------------------------------------------------------------
c         set the incoming short wave
c-----------------------------------------------------------------------

          dnswr(i,j) = dswr
c-----------------------------------------------------------------------
c         calculations only for ocean points
c-----------------------------------------------------------------------

          if (kmta(i,j) .eq. 0) then

            ws = spd(i,j)

# ifdef uvic_embm_astress
     &         + aspd(i,j)
# endif
            dt = surf(i,j,1) - at(i,j,2,1)
      dalt = 1.0022e-3 + 8.22e-5*(dt) + 2.66e-6*ws
      dalt =  max(dalt, 6.0e-05)
      dalt =  min(dalt, 2.190e-03)
#ifdef xlth
      if(dalt .lt. 2.19e-03)then
       fe=fb*1.5
      else
       fe=fb
      endif
#endif

c-----------------------------------------------------------------------
c           calculate evaporation or sublimation (ensure it is positive)
c-----------------------------------------------------------------------

#ifdef xlth
          evap(i,j)=max(0.0, fe*dalt*ws*(surf(i,j,2) - at(i,j,2,2)))
#else
           evap(i,j) = max(0.0, fb*dalt*ws*(surf(i,j,2) - at(i,j,2,2)))
#endif
           upltnt(i,j) = fc*evap(i,j)
#ifdef enhancedso
           if(yta(i,j).le.-60.)then
           evap(i,j)=evap(i,j)*10.0
           upltnt(i,j) = upltnt(i,j)*10.0
           endif
#endif

c-----------------------------------------------------------------------
c           calculate upward sensible heat flux
c-----------------------------------------------------------------------

            upsens(i,j) = fd*dalt*ws*(dt)
#ifdef enhancedso
            if(yta(i,j).le.-60.)then
            upsens(i,j) =             upsens(i,j)*10.
            endif
#endif
c-----------------------------------------------------------------------
c           calculate upward longwave re-radiation
c-----------------------------------------------------------------------

# ifdef uvic_embm_clouds
            cf = 11.5 - 22.8*rh + 11.7*rh**2
c prevent negative implied cloud cover
c           cf=amax1(cf,0.4)
#ifdef cmeans3
           a2d3(i,j,11)=a2d3(i,j,11)+cf
#endif
            uplwr(i,j) = esocn*(at(i,j,2,1) + 273.15)**4 
     &                   *(0.39 - 0.05*(at(i,j,2,2)*1.61e3)**.5)*cf
     &                 + 4.0*esocn*(at(i,j,2,1) + 273.15)**3*dt               
# else
           uplwr(i,j) = esocn*(surf(i,j,1) + 273.15)**4
     &                 - esatm(i,j)*(at(i,j,2,1) + 273.15)**4
# endif
	  endif

#ifdef landevap
c---------------------------------------------------------------------
c  calculate evaporation over land
c---------------------------------------------------------------------

      if(kmta(i,j).eq.1)then

#ifdef mrw_rivers
      storage=0.0
      do jj=1,12
      do ii=2,13
      storage=storage+hriver((i-2)*12+ii,(j-1)*12+jj,npriv)
      end do
      end do
      storage=storage/144.
#endif
      if(j.ge.25.and.j.le.36)then
      aquifer=2.5+500*cos(5.*(j-30.5)*180./60.*2.*3.141592654/360.)**2
      elseif(j.ge.13.and.j.le.48)then
      aquifer=10.0
      else
      aquifer=2.5
      endif

      storage=amin1(storage/aquifer,c1)

      ws=amax1(200.,spd(i,j))*10.0

      dalt = 1.0022e-3 + 2.66e-6*ws
      dalt =  max(dalt, 6.0e-05)
      dalt =  min(dalt, 2.190e-03)

#ifdef xlth
      evap(i,j)=
     &  max(0.0,fe*storage*dalt*ws*(surf(i,j,2)-at(i,j,2,2)))
#else
      evap(i,j)=
     &  max(0.0,fb*storage*dalt*ws*(surf(i,j,2)-at(i,j,2,2)))
#endif
      upltnt(i,j) = fc*evap(i,j)
      endif
#endif

        enddo
      enddo
	
      return
      end


      subroutine precipitate
      implicit none

c=======================================================================
c     calculate precipitation explicitly and update humidity
c     based on code by a.fanning
c
c     Note: evaporation and precipitation are in cm/s
c           and humidities are in g/g
c
c     author:   m.eby   e-mail: eby@uvic.ca
c=======================================================================

# include "param.h"
# include "cembm.h"
# include "atm.h"
# include "levinda.h"
# include "iounit.h"
#ifdef mrw_rivers
# include "riv.h"
#endif
#include "flux.h"
      real fa,qmax,fx,f1,f2,f3,f4

c-----------------------------------------------------------------------
c     set appropriate constants
c-----------------------------------------------------------------------

      fa = rhoatm*hq/(rhoocn*dts)

      do j=2,jmtam1
        do i=2,imtam1

c-----------------------------------------------------------------------
c         check if specific humidity is greater than saturation 
c-----------------------------------------------------------------------

#ifndef mrw_rivers
        qmax = 0.85*ssh(i,j)
#else
#ifdef uvic_embm_mountains
        if((at(i,j,2,1)-elev(i,j)/14285.).lt.0.0)then
#else
        if(at(i,j,2,1).lt.0.0)then
#endif 
        qmax=0.95*ssh(i,j)
        else
        qmax=0.85*ssh(i,j)
        endif
#endif
          if (at(i,j,2,2) .gt. qmax) then
            precip(i,j) = fa*(at(i,j,2,2) - qmax)
            at(i,j,2,2) = qmax
          else 
            precip(i,j) = 0.0
            if (at(i,j,2,2) .lt. 0.0) then
              at(i,j,2,2) = 0.0
c martin
c              write(stdout,*) '=>Warning: negative q(i,j) at ', i, j
            endif
          endif


#if defined peforrivers || defined pforrivers
         precip(i,j)=peflux(i,j,3)/864000.
#endif
        enddo
      enddo

#ifdef cyclic
       do j=1,jmta
       precip(1,j)=precip(imta-1,j)
       precip(imta,j)=precip(2,j)
       end do
#endif
c-----------------------------------------------------------------------
c     calculate river runoff from continents
c-----------------------------------------------------------------------

#ifdef mrw_rivers
      call rivmodel
# endif
	
      return
      end


