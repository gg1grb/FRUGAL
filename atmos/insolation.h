c=================== include file "insolation.h" =======================
c
c     note: this solution of berger 1978 is valid only for 1,000,000 
c           years centered on present-day. for longer period the 
c           solution 1990 must be used (contact Berger for this)                
c                                                              
c                                                              
c     refer to :                                      
c     Berger A. 1978. a simple algorithm to compute long term
c                     variations of daily or monthly insolation
c                     contr. 18  Inst of Astronomy and Geophysics
c                     Universite Catholique de Louvain
c                     Louvain-la-Neuve  Belgium              
c                                                                 
c     Berger A. 1978. long term variations of daily insolation 
c                     and quaternary climatic changes
c                     J. of Atmospheric Sciences  35  2362-2367
c                              
c     nef  nob  nop have been reduced to  19  18  9  
c 
c     these are the expansion coefficients for minimal efficiency, they 
c     can be increased by referring to Berger, 1978:
c     eccentricity           ecc   table 1
c     obliquity              xob   table 2
c     longitude perihelion   perh  table 3
c                    
      integer nef,nob,nop
      parameter	(nef=19, nob=18, nop=9)
c
      real ae(nef), be(nef), ce(nef), ye(nef), ze(nef)
     &,	   aob(nob), bob(nob), cob(nob), yob(nob), zob(nob)
     &,	   aop(nop), bop(nop), cop(nop), yop(nop), zop(nop)
      real pir, pirr, step
     &,              test, xod, xop, prm, nd, ecc, perh, pre, xob
c
      common /insol/ ae, be, ce, ye, ze, aob, bob, cob, yob, zob
     &,              aop, bop, cop, yop, zop, pir, pirr, step
     &,              test, xod, xop, prm, nd, ecc, perh, pre, xob




