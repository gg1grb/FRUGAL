      real cmean3
#ifdef extractsurheatflux
      real surheatfluxmean3(imta,jmta)
#endif
      real h0mean3(imta,jmta),u0mean3(imta,jmta),v0mean3(imta,jmta)
     & ,wsxmean3(imta,jmta),wsymean3(imta,jmta)
#ifdef ocean_ice
      real aicemean3(imta,jmta),hicemean3(imta,jmta),
     & ticemean3(imta,jmta)
#ifdef uvic_ice_cf
      real uicemean3(imt,jmt),vicemean3(imt,jmt)
#endif
#endif
      real umean3(km,imt,jmt),vmean3(km,imt,jmt),
     &     tmean3(km,imt,jmt,nt)
#ifdef savew
      real storewmean3(km,imt,jmt)
#endif
#if defined coupled
      real atmosmean3(imta,jmta,40)
#endif
      common /fields3mean/ymean3,
#ifdef extractsurheatflux
     & surheatfluxmean3,
#endif
     & h0mean3,u0mean3,v0mean3,wsxmean3,wsymean3,
#ifdef ocean_ice
     & aicemean3,hicemean3,ticemean3,
#ifdef uvic_ice_cf
     & uicemean3,vicemean3,
#endif
#endif
     & umean3,vmean3,
     &     tmean3,
#ifdef savew
     & storewmean3,
#endif
#if defined coupled
     & atmosmean3
#endif
      
