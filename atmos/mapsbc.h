c====================== include file "mapsbc.h" ========================
c
c     indicies for surface boundary condition arrays.  names must be
c     unique and are set at the begining of "driver"
c
c     itaux is the x component of wind stress
c     itauy is the y component of wind stress
c     ihflx is heat flux
c     ipme is for precipitation minus evaporation
c     ipsw is for penitrating shortwave (into ocean first layer)
c     isst is the ocean model SST
c     isss is the ocean model SSS
c     isu is the x component of the ocean model first layer velocity
c     isv is the y component of the ocean model first layer velocity
c     it10 is the 10 metre atmopheric temperature
c     iq10 is the 10 metre atmopheric humidity
c     ilwr is incoming longwave radiation
c     iswr is incoming shortwave radiation
c     igu is the x component of the ocean model second layer velocity
c     igv is the y component of the ocean model second layer velocity
c
      integer itaux, itauy, ihflx, ifflx,ipme, ipsw, isst,
     &            isss, isu, isv, it10, iq10, ilwr, iswr, igu, igv

      common /mapsbc/ itaux, itauy, ihflx, ifflx,ipme, ipsw, isst,
     &            isss, isu, isv, it10, iq10, ilwr, iswr, igu, igv
c
