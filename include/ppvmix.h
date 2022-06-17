c====================== include file "ppmix.h" =======================
#ifdef hcomments
c
c  variables for vertical mixing scheme of pacanowski & philander 
c     fricmx = max vertical mixing coefficient
c     wndmix = min vertical mixing in level 1 to simulate wind mixing
c
c  derived quantities:
c     drhom  = vertical density difference on "t" grid box bottoms
c     riu    = richardson number for "u" grid box bottoms
c---------------------------------------------------------------------
#endif

      real    riu(km,imax,jmax)
      real    wndmix, fricmx, epsln
      parameter (epsln=1.0e-20)
      common  /cppmix/ 
     &        riu, wndmix, fricmx

