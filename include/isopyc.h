c====================== include file "isopyc.h" ======================
c
c  variables for isopycnal diffusion scheme
c
c    ahisop = along isopycnal diffusion coefficient for mixing of
c               tracers  [cm**2/sec]
c               ("ah" = in "chmix.h is used as background horizontal
c               diffusion to prevent checkerboarding)
c    slmxr  = maximum slope of isopycnal allowed in mixing scheme
c               limited to prevent excessively large vertical mixing
c               that could create numerical instabilities
c               a value of 200 translates to a slope of 1:200
c
c    fk#  = isopycnal diffusion tensor components (compute five)
c          fk1(3) is resolved at center of eastern face of t,s box
c          fk2(3) is resolved at center of northern face of t,s box
c          fk3(1-3) are resolved at center of top face of t,s box
c     components (1,1) & (2,2) of tensor are assumed = 1.0
c     components (1,2) & (2,1) of tensor are assumed = 0.0
c
c    (in isop0) 
c      e#  = density gradients (1=zonal, 2=meridional, 3=vertical)
c    (in isop1) 
c      e1  = diffusion at (i,j) east face due to sloping isopycnal
c      e1m = diffusion at (i+1,j) east face due to sloping isopycnal
c      e2  = diffusion at (i,j) north face due to sloping isopycnal
c      e2m = diffusion at (i,j+1) north face due to sloping isopycnal
c      e3  = diffusion at (i,j) top face due to sloping isopycnal
c=====================================================================

      real   fk1(km,imax,jmax,3:3), fk2(km,imax,jmax,3:3),
     &       fk3(km,imax,jmax,3)
      real   e1(km), e1m(km), e2(km), e2m(km), e3(km+1)
      real   ahisop, slmxr

      common /isopyc/ 
     &       fk1, fk2, fk3, e1, e1m, e2, e2m, e3, ahisop, slmxr

