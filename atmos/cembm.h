c====================== include file "cembm.h" =========================
c
c     parameters for use in the energy balance model (also see atm.h)
c
c     addflux		= logical flag for adding only even mode fluxes
c     nats		= number of atmospheric time steps since mixing
c     namix		= time steps between mixing (set in atmos.in)
c     lf		= time step flag (1=>leapfrog, 2=>forward)
c     pyear		= calendar year (set in atmos.in, -/+ = BC/AD)
c     dtatm 		= time step for atmosphere (set in atmos.in)
c     dts		= time step (2*dtatm=>leapfrog, dtatm=>forward)
c     co2ccn		= CO2 concentration
c     anthro		= radiative forcing by CO2
c     co2ccni		= initial CO2 concentration
c     co2yri		= last year of initial CO2 concentration
c     co2ccnf		= final CO2 concentration
c     co2yrf		= first year of final CO2 concentration
c     solarconst	= solar constant (g/s**3)
c     ssp		= sea surface pressure (g/cm/s**2)
c     cdatm 		= drag coefficient (dimensionless)
c     cpatm		= atmospheric heat capacity (cm**2/s**2/K)
c     ht		= scale height for temperature
c     hq		= scale height for specific humidity
c     rhoatm		= density of air at sea surface (g/cm**3)
c     rhoocn		= representative sea surface density
c     esocn		= ocean emissivity times Stefan's constant 
c     vlocn		= latent heat of vaporization of water
c     socn		= constant ocean salinity
c     cdice 		= drag coefficient (dimensionless)
c     rhoice		= ice density (g/cm**3)
c     esice		= ice emissivity times Stefan's constant 
c     slice		= latent heat of sublimation of ice
c     flice		= latent heat of fusion of ice (cm**2/s**2)
c     condice		= ice conductivity (g*cm/s**3/K)
c     kadv		= time step divider for ice advection
c     ntats		= time step counter for time averaging
c     avetime		= time for averaging fluxes
c     totaltime		= total time for long term averages
c
      logical addflux
      integer nats, namix, lf, kadv, ntats
      real pyear, dtatm, dts, co2ccn, anthro, co2ccni, co2yri
     &,              co2ccnf, co2yrf, solarconst, ssp, cdatm, cpatm, ht
     &,              hq, rhoatm, rhoocn, esocn, vlocn, socn, cdice
     &,              rhoice, esice, slice, flice, condice, avetime
     &,              totaltime,hadvect,co2startyear,co2endyear
     &,              co2ccnv(29)
c
      real function Dalton
      common /cembm/ addflux
      common /cembm/ nats, namix, lf, kadv, ntats
      common /cembm/ pyear, dtatm, dts, co2ccn, anthro, co2ccni, co2yri
     &,              co2ccnf, co2yrf, solarconst, ssp, cdatm, cpatm, ht
     &,              hq, rhoatm, rhoocn, esocn, vlocn, socn, cdice
     &,              rhoice, esice, slice, flice, condice, avetime 
     &,              totaltime,hadvect,co2startyear,co2endyear,co2ccnv
      common /dal/ function Dalton
c
