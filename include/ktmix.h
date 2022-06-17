c==================== include file "ktmix.h" =========================
#ifdef hcomments
c
c   variables used in the kraus-turner mixed layer module
c
c     kt_lamda   = efficiency of the wind mixing (0.7)
c     kt_epsln   = efficiency of buoyancy driven mixing (0.15)
c     kt_wndh    = wind mixing scale depth (10,000cm)
c     decay(km)  = exponential defining the depth decay of wind mixing
c     rho_zero   = reference density of seawater
c     delpsf     = potential energy change per unit
c                  'non-solar density change' (cm3/s2)
c
c   coefficients for mixing the layers
c     xcz, xc1, xc2, xc3
c---------------------------------------------------------------------
#endif

      real     decay(km), xcz(km), xc1(km), xc2(km), xc3(km)
      real     kt_wndh, kt_lamda, kt_epsln, delpsf
      common  /ktmod/ 
     &         decay, xcz, xc1, xc2, xc3,
     &         kt_wndh, kt_lamda, kt_epsln, delpsf
      real   pdz(km), eqvdep(km)
      real   wndh, xlamda, xepsilon, rho_air
      common /ktmod/
     &       wndh, xlamda, xepsilon,
     &       rho_air


