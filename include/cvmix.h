c===================== include file "cvmix.h" ========================
#ifdef hcomments
c
c  variables used for vertical diffusion
c     fkph   = vertical diffusion coefficient (cm**2/sec)
c     fkpm   = vertical viscosity coefficient (cm**2/sec)
c     bvdc   = background vertical diffusion constant
c     bvvc   = background vertical viscosity constant
c     aidif  = coefficient for implicit time differencing for
c              vertical diffusion. aidif=1 gives the fully implicit
c              case. aidif=0 gives the fully explicit case
#if defined ppvmix || defined implicitvmix
c     vvclim = vertical viscosity coefficient limit
c     vdclim = vertical diffusion coefficient limit
#endif
#if defined mixing_alpha || defined full_alpha
c     alpham gives enhanced mixing coefficient
#endif
c
c  derived quantities:
c     vvc  = vertical viscosity coeff
c     vdc  = vertical diffusion coeff
c---------------------------------------------------------------------
#endif

      real    vvc(km), vdc(km),vdcbg(km)
      real    fkph, fkpm, bvdc, bvvc, aidif
#if defined ppvmix || defined implicitvmix
      real    vvclim, vdclim
#endif
#if defined mixing_alpha || defined full_alpha
      real    alpham
#endif
      common  /cvmix/
     &        vvc, vdc, fkph, fkpm, bvdc, bvvc, aidif,vdcbg
#if defined ppvmix || defined implicitvmix
     &,       vvclim, vdclim
#endif
#if defined mixing_alpha || defined full_alpha
     &,       alpham
#endif

