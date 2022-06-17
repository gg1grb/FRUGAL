      real cmean1
#ifdef extractsurheatflux
      real surheatfluxmean1(imta,jmta)
#endif
      real h0mean1(imta,jmta),u0mean1(imta,jmta),v0mean1(imta,jmta)
     & ,wsxmean1(imta,jmta),wsymean1(imta,jmta)
#ifdef ocean_ice
      real aicemean1(imta,jmta),hicemean1(imta,jmta),
     & ticemean1(imta,jmta)
#ifdef uvic_ice_cf
      real uicemean1(imt,jmt),vicemean1(imt,jmt)
#endif
#endif
      real umean1(km,imt,jmt),vmean1(km,imt,jmt),
     &     tmean1(km,imt,jmt,nt)
#ifdef savew
      real storewmean1(km,imt,jmt)
#endif
#if defined coupled
      real atmosmean1(imta,jmta,40)
#endif
      common /fields1mean/ymean1,
#ifdef extractsurheatflux
     & surheatfluxmean1,
#endif
     & h0mean1,u0mean1,v0mean1,wsxmean1,wsymean1,
#ifdef ocean_ice
     & aicemean1,hicemean1,ticemean1,
#ifdef uvic_ice_cf
     & uicemean1,vicemean1,
#endif
#endif
     & umean1,vmean1,
     &     tmean1,
#ifdef savew
     & storewmean1,
#endif
#if defined coupled
     & atmosmean1
#endif
      
