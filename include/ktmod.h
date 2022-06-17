c==================== include file "ktmod.h" =========================
#ifdef hcomments
c
c   variables used in the mixed layer model
c
c     pdz(km)    = equivalent layer thicknesses (distorted spinup)
c     eqvdep(km) = equivalent layer depth       (  ditto  )
c
c     decay(km)  = exponential defining the depth decay of wind mixing
c     wndh       = wind mixing scale depth (10,000cm)
c
c     xlamda     = efficiency of the wind mixing (0.7)
c     xepsilon   = efficiency of buoyancy driven mixing (0.15)
c
c     rho_zero   = reference density of seawater
c     rho_air    = density of air
c
c     delpsf     = potential energy change per unit
c                  'non-solar density change' (cm3/s2)
c---------------------------------------------------------------------
#endif

      real   pdz(km), eqvdep(km), decay(km)
      real   wndh, xlamda, xepsilon, rho_zero, rho_air, delpsf
      common /ktmod/ 
     &       pdz, eqvdep, decay, wndh, xlamda, xepsilon,
     &       rho_zero, rho_air, delpsf
