c====================== include file "atm.h" ===========================
c
c     arrays for the energy-moisture balance model
c
c     note: units for heat flux are in g/cm**2/s**3/cm**2 or (mW/m**2)
c           units for fresh water flux are in cm/s
c           downward is into the surface (ocean, ice or land)
c           upward is into the bottom of the atmosphere
c           outward is out of the top of the atmosphere
c           inward is into the top of the atmosphere
c
c     outlwr(imta,jmta)	 = outgoing atmosphere longwave flux
c     at(imta,jmta,2,nat)	 = tracers (previous and most recent)
c     surf(imta,jmta,nat)	 = surface values (1=temperature, 2=humidity)
c     flux(imta,jmta,nat)	 = downward flux (1=heat, 2=fresh water)
c     precip(imta,jmta)	 = precipitation
c     evap(imta,jmta)	 = evaporation
c     dnswr(imta,jmta)  	 = downward surface shortwave flux
c     uplwr(imta,jmta)  	 = upward surface longwave flux
c     upsens(imta,jmta)	 = upward surface sensible heat flux
c     upltnt(imta,jmta)    = upward surface latent heat flux
c     ssh(imta,jmta) 	 = saturation specific humidity
c     spd(imta,jmta) 	 = wind speed (cm/sec)
c     da(jmta,nat)	 = tracer diffusivity parameterization
c     coalbedo(jmta)	 = planetary coalbedo
c     solins(jmta)	 = solar insolation (units of solarconst)
c     esatm(jmta)	 = atmosphere emissivity times Stefan's constant
c     scatter(jmta)	 = solar scattering into atmosphere
#ifdef uvic_embm_planetary
c     esplnt(jmta)	 = planetary emissivity times Stefan's constant
#endif
#ifdef uvic_embm_seasonal_average
c     rtbar(imta,jmta)	 = running average atmospheric temperature
c     rqbar(imta,jmta)	 = running average specific humidity
c     atbar(imta,jmta)	 = average temperature accumulator
c     aqbar(imta,jmta)	 = average humidity accumulator
#endif
#ifdef uvic_embm_mountains
c     elev(imta,jmta)	 = elevations (in cm) or factor for temperature
#endif
#ifdef uvic_embm_astress
c     atx(imta,jmta)	 = anomalous taux
c     aty(imta,jmta)	 = anomalous tauy
c     aspd(imta,jmta)	 = anomalous wind speed
c     tbar(imta,jmta)	 = average atmospheric temperature
c     qbar(imta,jmta)	 = average specific humidity
#endif
c     cori(jmta,2)	 = Coriolis parameter
#if defined uvic_embm_astress || defined uvic_embm_gstress
c     mask(imta,jmta)	= mask for stress feedback calculations
c
      integer maskw(imta,jmta)
c     logical lastatmos
      common /atm/ maskw
c     common /atm/ lastatmos
#endif
c
      real outlwr(imta,jmta),at(imta,jmta,2,nat),surf(imta,jmta,nat)
      real flux(imta,jmta,nat), precip(imta,jmta), evap(imta,jmta)
      real dnswr(imta,jmta), uplwr(imta,jmta), upsens(imta,jmta)
      real upltnt(imta,jmta),ssh(imta,jmta),spd(imta,jmta),
     &     da(imta,jmta,nat)
      real toasolar(imta,jmta),toatot(imta,jmta)
      real coalbedo(imta,jmta), solins(imta,jmta), 
     &      esatm(imta,jmta), scatter(imta,jmta)

#ifdef uvic_embm_planetary
      real esplnt(imta,jmta)
#endif
#ifdef uvic_embm_seasonal_average
      real rtbar(imta,jmta), rqbar(imta,jmta)
      real atbar(imta,jmta), aqbar(imta,jmta)
#endif
#if defined uvic_embm_mountains 
      real elev(imta,jmta)
#endif
#if defined uvic_embm_astress || defined uvic_embm_gstress
      real atx(imta,jmta), aty(imta,jmta), aspd(imta,jmta)
      real tbar(imta,jmta), qbar(imta,jmta),zat(jmta),zah(jmta)
     &     ,zwx(jmta),zwy(jmta)
#endif
      real cori(imta,jmta,2),totsurheatflux
c#ifdef advectq
      real spdu(imta,jmta),spdv(imta,jmta)
c#endif
#ifdef adjustheatdiffusivity
      real heatdiffusivityfbnorth,heatdiffusivityfbsouth
#endif
#ifdef adjustwvdiffusivity
      real wvdiffusivityfbnorth,wvdiffusivityfbsouth
#endif
c
      common /atm/ outlwr, at, surf, flux, precip, evap, dnswr, uplwr
     &,            upsens, upltnt, spd, ssh, da, coalbedo, solins
     &,            esatm, scatter,toasolar,toatot
#ifdef uvic_embm_planetary
     &,	           esplnt
#endif
#ifdef uvic_embm_seasonal_average
     &,	           rtbar, rqbar, atbar, aqbar
#endif
#if defined uvic_embm_mountains 
     &,            elev
#endif
#ifdef uvic_embm_astress
     &,            atx, aty, aspd, tbar, qbar,zat,zah,zwx,zwy
#endif
     &,            cori,totsurheatflux
#ifdef advectq
     &, spdu,spdv
#endif
#ifdef adjustheatdiffusivity
     & ,heatdiffusivityfbnorth,heatdiffusivityfbsouth
#endif
#ifdef adjustwvdiffusivity
     & ,wvdiffusivityfbnorth,wvdiffusivityfbsouth
#endif
c

