c====================== include file "grdvar.h" ========================
c
c     variables which are functions of the grid defined by "coord.h"
c
c     dxt     = longitudinal width of "t" grid box at the 
c               equator (in cm)
c     dxtr    = reciprocal of "dxt"
c     dxt2r   = reciprocal of "2*dxt"
c     dxt4r   = reciprocal of "4*dxt"
c     dxu     = longitudinal width of "u,v" grid box at the 
c               equator (in cm)
c     dxur    = reciprocal of "dxu"
c     dxu2r   = reciprocal of "2*dxu"
c     dxu4r   = reciprocal of "4*dxu"
c     dxmetr  = reciprocal of "(dxt(i)+dxt(i+1))"
c     duw     = xu(i) - xt(i)
c     due     = xt(i+1) - xu(i)
c     dus     = yu(jrow) - yt(jrow)
c     dun     = yt(jrow+1) - yu(jrow)
c
c     dyt     = latitudinal height of "t" grid box (in cm)
c     dytr    = reciprocal of "dyt"
c     dyt2r   = reciprocal of "2*dyt"
c     dyt4r   = reciprocal of "4*dyt"
c     dyu     = latitudinal height of "u,v" grid box (in cm)
c     dyur    = reciprocal of "dyu"
c     dyu2r   = reciprocal of "2*dyu"
c     dyu4r   = reciprocal of "4*dyu"
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
      real dxta,  dxtra,  dxt2ra, dxua
      real dxura, dxu2ra, dxu4ra, dxt4ra
      real dyta,  dytra,  dyt2ra, dyua
      real dyura, dyu2ra, dyu4ra, dyt4ra
      real csua,  csura,  csta,   cstra
      real cstdytra, cstdyt2ra
      real csudyura, csudyu2ra
      real cst_dytra, csu_dyura
      real phia,  phita,  sinea,  tnga
      real fcora
      real dxmetra, duwa, duea
      real duna, dusa

      common /grdvara/ dxta(imta,jmta),dxtra(imta,jmta),
     &                 dxt2ra(imta,jmta),dxua(imta,jmta)
      common /grdvara/ dxura(imta,jmta),dxu2ra(imta,jmta),
     &                 dxu4ra(imta,jmta),dxt4ra(imta,jmta)
      common /grdvara/ dyta(imta,jmta),dytra(imta,jmta),
     &                 dyt2ra(imta,jmta),dyua(imta,jmta)
      common /grdvara/ dyura(imta,jmta),dyu2ra(imta,jmta),
     &                 dyu4ra(imta,jmta),dyt4ra(imta,jmta)
      common /grdvara/ csua(imta,jmta),csura(imta,jmta),
     &                 csta(imta,jmta),cstra(imta,jmta)
      common /grdvara/ cstdytra(imta,jmta),cstdyt2ra(imta,jmta)
      common /grdvara/ csudyura(imta,jmta),csudyu2ra(imta,jmta)
      common /grdvara/ cst_dytra(imta,jmta),csu_dyura(imta,jmta)
      common /grdvara/ phia(imta,jmta),phita(imta,jmta),
     &                 sinea(imta,jmta),tnga(imta,jmta)
      common /grdvara/ fcora(imta,jmta) 
      common /grdvara/ dxmetra(imta,jmta),duwa(imta,jmta),
     &                 duea(imta,jmta)
      common /grdvara/ duna(imta,jmta),dusa(imta,jmta)
