c===================== include file "ctmngr.h" =======================
#ifdef hcomments
c
c     time manager clock parameters 
c
c     itt     = current time step number
c     itbt    = current barotropic timestep
c     itbtp   = current pass of free surface model timestepping scheme
c     totsec  = accumulated time in seconds from jan 1, year 0
c     totday  = accumulated time in days
c     years   = accumulated time in years
c     stamp   = date & time stamp corresponding to itt
c     totsec_start  = time in days run started
c                     (to fix rounding errors in tmngr.F)
c---------------------------------------------------------------------
#endif

      character stamp*34,filenameyears*8
      real    totsec, totday, years, totsec_start
     & ,rmonth,fahead,fbehind
#ifdef daily
     & ,rday,dayone
#endif
#if defined nordic_tracer || defined nordic_heateof
     & ,dayone
#endif
#ifdef insolation_vary
     & ,modelstartyear
#endif
#ifdef era40_6hourly
     & ,hourbegin,tothour,hourone
#endif
      integer itt, itbt, itbtp, countloop
     & ,mahead,mbehind,maheadx
      common  /ctmngr/ 
     &        totsec, totday, years, totsec_start, itt, itbt, itbtp
     &,rmonth,fahead,fbehind,mahead,maheadx,mbehind
#ifdef daily
     & ,rday,dayone
#endif
#if defined nordic_tracer || defined nordic_heateof
     & ,dayone
#endif
#ifdef insolation_vary
     & ,modelstartyear
#endif
#ifdef era40_6hourly
     & ,hourbegin,tothour,hourone
#endif

      common  /ctstmp/
     &        stamp,filenameyears
#if defined coupled || defined iceonly
      integer ntimes,nfulltimes
      parameter (ntimes = 100, nfulltimes = 20)
      integer iday(ntimes), msday(ntimes)
      integer nextfulltime, nexttime
      integer initial, iref, irunstart, itime, iuser
      integer iruntime, imodeltime, ireftime, iusertime
      integer ihalfstep, isunday
      integer itemptime,itemptime2,itmptime,itmptime2,itmptime3
      integer idt, idtd2
      integer itt0, itta, irstdy, msrsdy
      real dayoyr
      common /tmngrii/ nextfulltime, nexttime
      common /tmngrii/ initial, iref, irunstart, itime, iuser
      common /tmngrii/ iruntime, imodeltime, ireftime, iusertime
      common /tmngrii/ ihalfstep, isunday
      common /tmngrii/ itemptime,itemptime2,itmptime,itmptime2,itmptime3
      common /tmngrii/ idt, idtd2
      common /tmngrii/ iday, msday
      common /tmngrii/ itt0, itta, irstdy, msrsdy,countloop
c
      common /tmngrr/ dayoyr
#endif
