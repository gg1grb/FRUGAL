c====================== include file "slap.h" ==========================
c
c     variables needed for the energy balance model solution technique
c
c     5-point stencil
c     max number non-zeros: total, upper and lower triangles (with diag)
c
      integer maxnelt, maxneltu, maxneltl
      parameter  (maxnelt = (jmta+2)*(3*imta) + 2*(imta+2)*(jmta+1))
      parameter  (maxneltu = 3*(jmta+2)*(imta+2) - imta - jmta)
      parameter  (maxneltl = maxneltu)
c
c     number direction vectors to orthogonalize against
c
      integer nsave
      parameter  (nsave = 10)
c
c     work space
c
      integer lenw,leniw
      parameter  (lenw = 1 + (imta)*(jmta)*(nsave+7) +
     &	         nsave*(nsave+3) + maxneltu + maxneltl)
      parameter  (leniw =  (imta)*(jmta)*4+maxneltu + maxneltl + 32)
      real       rwork(lenw)
      integer    iwork(leniw)
c
c     the tracer matrix
c
      integer    slapia(maxnelt,2,nat), slapja(maxnelt,2,nat)
      integer    nelt, orda
c
      real       slapam(maxnelt,2,nat), slaprhs(imta*jmta,nat)
      real       slaplhs(imta*jmta,nat), sltol(nat)
c
      common /slapmats/ slapia, slapja, nelt, orda
     &,	                slapam, slaprhs, slaplhs, sltol
c
c     index maps
c
      integer    indexgrid(imta,jmta), gridrow(imta*jmta)
      integer    gridcol(imta*jmta)
c
      common /slapmaps/ indexgrid, gridrow, gridcol
c
