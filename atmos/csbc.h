c====================== include file "csbc.h" ==========================
c
c                    surface boundary conditions (S.B.C.)
c
c     numosbc = number of ocean S.B.C..
c     numasbc = number of atmospheric S.B.C..
c     maxsbc  = maximum number of allowable S.B.C. fields.
c     numsbc  = total number of S.B.C. this is the sum of the 
c               atmosphere S.B.C. and ocean S.B.C. fields. it must
c
c     sbcocn  = S.B.C. fields on the ocean grid.
c               n=1..numosbc are for the incoming ocean S.B.C.
c               n=numosbc+1..numosbc+numasbc are for the outgoing S.B.C.
c     sbcname  = S.B.C. character string names (eg: "SST", "taux", etc)
c               note that these character strings should be changed with
c               care since they are used to identify fields for special
c               treatment
c     dunits  = dimensional units for S.B.C. (eg: "deg C")
c     mapsbc  = indices for relating how S.B.C. fields are stored
c               within arrays "sbcocn" and "sbcatm"  
c               eg: m = mapsbc(n) where sbcocn(,,m) is S.B.C. "n"
c
c
      integer numosbc,numasbc,maxsbc,numsbc
c
      parameter (numosbc=5
     &,          numasbc=4
c     &,          maxsbc=numosbc+numasbc, numsbc=maxsbc)
     &,          maxsbc=6, numsbc=6)
      logical asegs, asege, afirst, alast
      character*10 sbcname(maxsbc)
      character*15 dunits(maxsbc)
c
      real sbcocn(imta,jmta,-1:numsbc)
      common /csbc/ sbcocn 
      real coabc(maxsbc), crits(maxsbc)
      common /csbc/ coabc, crits
      integer mapsbc(maxsbc), numpas, iombc, ntspos
      common /csbc1/ mapsbc, numpas, iombc, ntspos
      common /csbc2/ sbcname, dunits
      common /csbc3/ asegs, asege, afirst, alast 
c
