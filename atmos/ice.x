c====================== include file "ice.h" ===========================
c
c     arrays for the ice model
#ifdef uvic_embm_ice
c
c     hice(imta,jmta,3)	 = ice thickness (tau-1, tau, tau+1)
c     aice(imta,jmta,3)	 = ice area (tau-1, tau, tau+1)
c     tice(imta,jmta)	 = ice surface temperature
c     frzpt(imta,jmta)	 = freezing temperature of sea water
c     coalice(jmta)	 = planetary coalbedo with ice
# ifdef uvic_embm_landice
c      icemask		 = mask for constant ice albedo over land
# endif
# ifdef uvic_ice_cf
c     uice(imta,jmta,2)	 = u ice velocity (previous and most recent)
c     vice(imta,jmta,2)	 = v ice velocity (previous and most recent)
c     presur(imta,jmta)	 = pressure due to internal stress
c     uvmx(imta,jmta)	 = mask for u C-grid ice velocities
c     uvmy(imta,jmta)	 = mask for v C-grid ice velocities
c     kfd		 = iteration counter for free drift solver
c     kcf		 = iteration counter for cavitating fluid solver
c
      common /ice/ kfd, kcf
c
      real uice(imta,jmta,2), vice(imta,jmta,2), presur(imta,jmta)
      real uvmx(imta,jmta), uvmy(imta,jmta)
      common /ice/ uice, vice, presur, uvmx, uvmy
# endif
c
      real hice(imta,jmta,3), aice(imta,jmta,3), tice(imta,jmta)
      real frzpt(imta,jmta), coalice(jmta)
# ifdef uvic_embm_landice
      integer icemask(imta,jmta)
# endif
c
      common /ice/ hice, aice, tice, frzpt, coalice
# ifdef uvic_embm_landice
      common /icemask/ icemask
# endif
#endif
c
