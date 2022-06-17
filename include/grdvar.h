c===================== include file "grdvar.h" =======================
#ifdef hcomments
c
c     @(#) SCCS module: grdvar.h, version 1.12
c
c  This version assumes all dx and dy boxes are of equal size when
c  measured in degrees.  Variables dxt, dxu etc replaced by dx, 
c  dyt,dyu etc by dy and dzt by dz.
c
c     dx     = longitudinal width of "t" grid box (in cm)
c     dxr    = reciprocal of "dx"
c     dx2r   = reciprocal of "2*dx"
c     dx4r   = reciprocal of "4*dx"
c
c     dy     = latitudinal height of "t" grid box (in cm)
c     dyr    = reciprocal of "dy"
c     dy2r   = reciprocal of "2*dy"
c     dy4r   = reciprocal of "4*dy"
c
c     csu     = cosine of "u,v" grid point latitude
c     csur    = reciprocal of "csu"
c     cst     = cosine of "t" grid point latitude
c     cstr    = reciprocal of "cst"
c     phi     = latitude of "u,v" grid point in radians
c     phit    = latitude of "t" grid point in radians
c     sine    = sine of "u,v" grid point latitude
c     tng     = tan of "u,v" grid point latitude
c     fcor    = 2*omega*sine(j)
c
c     dz(k)  = level thickness of "t" and "u,v" grid boxes (in cm)
c               (vertical separation between "w" velocity points)
c     c2dz(k)= "2*dzt"
c     dzr(k) = reciprocal of dz
c     dz2r(k)= reciprocal of "2*dz"
c     dzw(k)  = vertical separation between centers of levels k & k+1
c     dzwr(k) = reciprocal of dzw
c     dzw2r(k)= reciprocal of "2*dzw"
c     dzur(k)= upper diffusion grid factor = 1.0/(dzw(k-1)*dz(k))
c     dzlr(k)= lower diffusion grid factor = 1.0/(dzw(k)*dz(k))
c     dzwur(k)= upper diffusion grid factor = 1.0/(dz(k)*dzw(k))
c     dzwlr(k)= lower diffusion grid factor = 1.0/(dz(k+1)*dzw(k))
c     tanra(j)= tan(j)/radius on the "u,v" grid
c---------------------------------------------------------------------
#endif

      real dxt(imt,jmt),dxtr(imt,jmt),dxt2r(imt,jmt),dxt4r(imt,jmt),
     &     dyt(imt,jmt),dytr(imt,jmt),dyt2r(imt,jmt),dyt4r(imt,jmt),
     &     dxu(imt,jmt),dxur(imt,jmt),dxu2r(imt,jmt),dxu4r(imt,jmt),
     &     dyu(imt,jmt),dyur(imt,jmt),dyu2r(imt,jmt),dyu4r(imt,jmt),
#ifdef uvic_ice_cf
     &     csu(imt,jmt),cstr(imt,jmt),
#endif
     &     phi(imt,jmt),phit(imt,jmt),
     &     sine(imt,jmt),tng(imt,jmt),
     &     fcor(imt,jmt),tanra(imt,jmt)
      real    dz(km), dzr(km), dz2r(km), c2dz(km), dzur(km), dzlr(km)
      real    dzw(0:km), dzwr(0:km), dzw2r(0:km)
      common  /grdvar/
     &     dxt,dxtr,dxt2r,dxt4r,
     &     dyt,dytr,dyt2r,dyt4r,
     &     dxu,dxur,dxu2r,dxu4r,
     &     dyu,dyur,dyu2r,dyu4r,
#ifdef uvic_ice_cf
     &     csu,cstr,
#endif
     &     phi,phit,
     &     sine,tng,
     &     fcor,tanra,
     &     dz, dzr, dz2r, c2dz, dzw, dzwr, dzw2r, dzur, dzlr
