c====================== include file "levind.h" ========================
c
c     vertical level indicators which define model geometry & bottom
c     topography:
c
c     kmt = number of vertical boxes over "t" points
c     kmu = number of vertical boxes over "u,v" points
c
      integer kmta,kmua
      common /levinda/ kmta(imta,jmta), kmua(imta,jmta)
c
