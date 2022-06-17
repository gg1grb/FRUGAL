#ifdef mrw_rivers
c====================== include file "riv.h" ===========================
c
c
      real dxriv(imtr,jmtr),dyriv(imtr,jmtr),cstriv(imtr,jmtr)
      real runoff(imta,jmta)
      real hriver(imtr,jmtr,2)
      integer ltop(imtr,jmtr)
      integer nriv,npriv 

      common /riv/ hriver,runoff
      common /riv/ dxriv,dyriv,cstriv
      common /riv/ ltop,nriv,npriv
#endif
