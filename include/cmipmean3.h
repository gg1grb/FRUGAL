      integer no3d3,no2d3,na2d3,noza3,nozb3
      real cmean3
      parameter(no3d3=16,no2d3=8,na2d3=18,noza3=4,nozb3=8)
      real o3d3(imt,jmt,0:km,no3d3),o2d3(imt,jmt,no2d3),
     & a2d3(imta,jmta,na2d3),oza3(-90:90,km,noza3),
     & ozb3(-90:90,nozb3)

      common /cmip3/ o3d3,o2d3,a2d3,oza3,ozb3,cmean3

c o3d3indices
c  1 pot temp
c  2 salinity
c  3 u vel
c  4 v vel
c  5 w vel
c  6 hor kinetic energy
c  7 speed of sound
c  8 density
c  9 northward heat transport advective
c 10 northward heat transport diffusive
c 11 eastward heat transport advective
c 12 eastward heat transport diffusive
c 13 northward salt transport advective
c 14 northward salt transport diffusive
c 15 eastward salt transport advective
c 16 eastward salt transport diffusive

c o2d3indices
c  1 ssh
c  2 sea ice thickness
c  3 sea ice fraction
c  4 sea ice u
c  5 sea ice v
c  6 fresh water into ocean
c  7 heat into ocean

c a2d3 indices
c  1 air pressure at sea level
c  2 total precip
c  3 temp
c  4 wsx
c  5 wsy
c  6 surface latent heat flux
c  7 surface sensible heat flux
c  8 net downward surface longwave
c  9 net down surface shortwave
c  10 surface temp
c  11 atmospheric water vapour content - in winds.F
c  12 runoff flux
c  13 u wind
c  14 v wind
c  15 specific humidity
c  16 incoming toa shortwave
c  17 outgoing toa shortwave
c  18 toa outgoing longwave

c oza3 indices
c  1 meridional volume transport atlantic
c  2 meridional volume transport pacific
c  3 meridional volume transport pacific/indian
c  4 meridional volume transport global

c ozb3 indices
c  1 meridional heat transport atlantic
c  2 meridional heat transport pacific
c  3 meridional heat transport pacific/indian
c  4 meridional heai transport global
c  5 meridional salt transport atlantic
c  6 meridional salt transport pacific
c  7 meridional salt transport pacific/indian
c  8 meridional salt transport global

