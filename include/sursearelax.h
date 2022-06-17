c----------------------------------------------------------------
c  contains arrays to store seasonally varying 
c  surface fields for surface relaxation
c---------------------------------------------------------------
      real rstracer(imt,jmt,nt,12)
      common /surrelax/ rstracer
#ifdef adjustrelax
      real rstraceradjust(imt,jmt,2,12)
      common /surrelax/ rstraceradjust
#endif
