      integer no3d2,no2d2,na2d2
      real cmean2
      parameter(no3d2=5,no2d2=1,na2d2=1)
      real o3d2(imt,jmt,0:km,no3d2),o2d2(imt,jmt,no2d2),
     & a2d2(imta,jmta,na2d2)

      common /cmip2/ o3d2,o2d2,a2d2,cmean2

c o3d2 indices
c  1 pot temp
c  2 salinity
c  3 u vel
c  4 v vel
c  5 w vel

c o2d2 indices
c  1 ssh

