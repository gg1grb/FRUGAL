c----------------------include "topog.h"-------------------------
c
c  contains topography array at 1/12th degree for ocean topography,
c  land topography,and runoff topography
c
      integer*2 topog(-2256:2256,-1152:1152)
      common /topogr/ topog
