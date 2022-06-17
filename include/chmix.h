c===================== include file "chmix.h" ========================
#ifdef hcomments
c
c     horizontal mixing coefficients
c
c     am  = lateral eddy viscosity (cm**2/sec)
c     ah  = lateral eddy diffusivity (cm**2/sec)
c     aq  = lateral turbulent kinetic energy diffusivity 
c           (cm**2/sec)
c---------------------------------------------------------------------
#endif

      real    am, ah(km)
      real    aau(imax,jmax),bbu(imax,jmax), 
     &        ccu(imax,jmax), ddu(imax,jmax), 
     &        ggu(imax,jmax),hhu(imax,jmax),
     &        aat(km,imax,jmax), bbt(km,imax,jmax), 
     &        cct(km,imax,jmax), ddt(km,imax,jmax)
      common  /chmix/
     &        am, ah, bbu, ccu, ddu, ggu, hhu, bbt, cct, ddt
     & ,aat,aau
