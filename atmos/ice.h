c====================== include file "ice.h" ===========================
c
c     arrays for the ice model
#if defined uvic_embm_ice || defined ocean_ice
c
c     hice(imt,jmt,3)	 = ice thickness (tau-1, tau, tau+1)
c     aice(imt,jmt,3)	 = ice area (tau-1, tau, tau+1)
c     tice(imt,jmt)	 = ice surface temperature
c     frzpt(imt,jmt)	 = freezing temperature of sea water
c     coalice(jmt)	 = planetary coalbedo with ice
# ifdef time_averages
c     ta_hi(imt,jmt)	 = time average ice thickness
c     ta_ai(imt,jmt)	 = time average ice area
c     ta_tice(imt,jmt)	 = time average ice surface temperature
c     ta_gice(imt,jmt)   = time averaged ice growth rate
# endif
# ifdef prescribelandice
c      icemask		 = mask for constant ice albedo over land
# endif
# ifdef uvic_ice_cf
c     uice(imt,jmt,2)	 = u ice velocity (previous and most recent)
c     vice(imt,jmt,2)	 = v ice velocity (previous and most recent)
c     presur(imt,jmt)	 = pressure due to internal stress
c     uvmx(imt,jmt)	 = mask for u C-grid ice velocities
c     uvmy(imt,jmt)	 = mask for v C-grid ice velocities
c     kfd		 = iteration counter for free drift solver
c     kcf		 = iteration counter for cavitating fluid solver
c
      integer kfd,kcf
      common /ice/ kfd, kcf
c
      real uice(imt,jmt,2), vice(imt,jmt,2), presur(imt,jmt)
      real uvmx(imt,jmt), uvmy(imt,jmt)
      common /ice/ uice, vice, presur, uvmx, uvmy
# endif
# ifdef uvic_ice_adv
c     kfd		 = iteration counter for free drift solver
      integer kfd
      common /ice/ kfd
c
#endif
c
      real hice(imt,jmt,3), aice(imt,jmt,3), tice(imt,jmt)
      real frzpt(imt,jmt),
     &    coalice(imt,jmt)
#ifdef mrw_ice_adv
      real uice(imt,jmt),vice(imt,jmt)
#endif
c
      common /ice/ hice, aice, tice, frzpt, coalice
#ifdef mrw_ice_adv
     & ,uice,vice
#endif
# ifdef time_averages
     &,            ta_hi, ta_ai, ta_tice, ta_gice
# endif
#endif
# ifdef prescribelandice
      integer icemask(imt,jmt)
      common /icemask/ icemask
# endif
c
