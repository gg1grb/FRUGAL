      real cmean2
#ifdef extractsurheatflux
      real surheatfluxmean2(imta,jmta)
#endif
      real h0mean2(imta,jmta),u0mean2(imta,jmta),v0mean2(imta,jmta)
     & ,wsxmean2(imta,jmta),wsymean2(imta,jmta)
#ifdef ocean_ice
      real aicemean2(imta,jmta),hicemean2(imta,jmta),
     & ticemean2(imta,jmta)
#ifdef uvic_ice_cf
      real uicemean2(imt,jmt),vicemean2(imt,jmt)
#endif
#endif
      real umean2(km,imt,jmt),vmean2(km,imt,jmt),
     &     tmean2(km,imt,jmt,nt)
#ifdef savew
      real storewmean2(km,imt,jmt)
#endif
#if defined coupled
      real atmosmean2(imta,jmta,40)
#endif
      common /fields2mean/ymean2,
#ifdef extractsurheatflux
     & surheatfluxmean2,
#endif
     & h0mean2,u0mean2,v0mean2,wsxmean2,wsymean2,
#ifdef ocean_ice
     & aicemean2,hicemean2,ticemean2,
#ifdef uvic_ice_cf
     & uicemean2,vicemean2,
#endif
#endif
     & umean2,vmean2,
     &     tmean2,
#ifdef savew
     & storewmean2,
#endif
#if defined coupled
     & atmosmean2
#endif
      
