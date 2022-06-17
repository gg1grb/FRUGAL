c-----------------------------------------------------------------
c data for bulk formula fluxes
c-----------------------------------------------------------------
       real windx(imt,jmt,12),windy(imt,jmt,12),
     &      wvx(imt,jmt),wvy(imt,jmt)
#if defined nordic_tracer || defined nordic_heateof
     &     ,eofws(2,imt,jmt)
#endif
#ifdef coupled
     &     ,heatflux(imt,jmt)
#ifndef coupledfwclim
     &     ,freshflux(imt,jmt)
#else
     &     ,freshflux(imt,jmt,3)
#endif
#endif
#ifdef surheatflux
     &     ,heatflux(imt,jmt,3)
#endif
#ifdef surfreshflux
     &     ,freshflux(imt,jmt,3)
#endif
#if defined peforrivers || defined pforrivers
     &  ,peflux(imt,jmt,3)
#endif
#ifdef tstar
     &  ,couplest
#endif
#ifdef natdiag
     & ,naqcon(4)
#endif
       common /fluxdata/ windx,windy,wvx,wvy
#ifdef nordic_tracer
     &     ,eofws
#endif
#ifdef surheatflux
     &     ,heatflux
#endif
#ifdef surfreshflux
     &     ,freshflux
#endif
#ifdef coupled
     &     ,heatflux
#ifndef coupledfwclim
     &     ,freshflux
#endif
#endif
#if defined peforrivers || defined pforrivers
     &  ,peflux
#endif
#ifdef tstar
     &  ,couplest
#endif
#ifdef natdiag
     & ,naqcon
#endif

