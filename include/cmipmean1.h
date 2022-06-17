      integer no3d1,no2d1,na2d1
      real cmean1
      parameter(no3d1=5,no2d1=1,na2d1=1)
      real o3d1(imt,jmt,0:km,no3d1),o2d1(imt,jmt,no2d1),
     & a2d1(imta,jmta,na2d1)

      common /cmip1/ o3d1,o2d1,a2d1,cmean1

c o3d1 indices
c  1 pot temp
c  2 salinity
c  3 u vel
c  4 v vel
c  5 w vel

c o2d1 indices
c  1 ssh

