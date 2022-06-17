#ifdef mrw_rivers
      subroutine rivinit
      implicit none
#include "param.h"
#include "grdvara.h"
#include "riv.h"
#include "levinda.h"
      real nseg1,nseg2,fx
      integer ii,jj,ll,it,testpoint,testfield,nhole

c------------------------------------------------------------
c set up runoff grid with a resolution 1/dxyscale
c that of the atmospheric model
c------------------------------------------------------------

      do i=2,imtr-1
      ii=int((i-2)/dxyscale)+2
      dxriv(i)=dxta(ii)/dxyscale
      end do
#ifdef cyclic
      dxriv(1)=dxriv(imtr-1)
      dxriv(imtr)=dxriv(2)
#endif
 
      do j=1,jmtr
      jj=int((j-1)/dxyscale)+1
      dyriv(j)=dyta(jj)/dxyscale
      end do

      do j=1,jmtr
      csuriv(j)=cosd((j-jmtr/2.)/jmtr*180.)
      cstriv(j)=cosd((j-0.5-jmtr/2)/jmtr*180.)
      end do

c----------------------------------------------------------------
c  read in land topography with a resolution of 1/3x1/4 degree
c  note values are in m but only relative height needed
c----------------------------------------------------------------

      open(46,file='data/landtopog',status='old')
      do j=jmtr,1,-1
      read(46,78)(ltop(i,j),i=3*(imtr-2)/4+2,imtr-1),
     &           (ltop(i,j),i=2,3*(imtr-2)/4+1)
      end do

   78 format(13i6)

c---------------------------------------------------------------
c  set points at or below sea level to 1m
c  (points which are actally ocean in the OGCM are then set to 0m)
c---------------------------------------------------------------

      do j=1,jmtr
      do i=1,imtr
      if(ltop(i,j).le.0)ltop(i,j)=1
      end do
      end do

c--------------------------------------------------------------
c set points defined as ocean by kmta to 0m
c--------------------------------------------------------------

      do j=1,jmtr
      do i=2,imtr-1
      ii=int((i-2)/dxyscale)+2
      jj=int((j-1)/dxyscale)+1
      if(kmta(ii,jj).eq.0)ltop(i,j)=0
      end do
      end do

c----------------------------------------------------------------
c fill in holes
c----------------------------------------------------------------

      it=0
  500 continue
      nhole=0
      it=it+1
      testfield=0

      do j=2,jmtr-1
      do i=2,imtr-1

      if(ltop(i,j).gt.0)then

      testpoint=0
      do l=-1,1
      do ll=-1,1
      if(l.eq.0.and.ll.eq.0)then
      else
      if(ltop(i+l,j+ll).lt.ltop(i,j))then
      testpoint=1
      endif
      endif
      end do
      end do

      if(testpoint.eq.0)then
      ltop(i,j)=ltop(i,j)+10
      testfield=1
      nhole=nhole+1
      endif

      endif
      end do
      end do

      if(testfield.eq.1)then
      if(mod(nhole,10).eq.0)
     & print *,'iteration ',it,'no. of holes is ',nhole
      goto 500
      endif

#ifdef cyclic
c--------------------------------------------------------------
c  set cyclic boundary conditions
c--------------------------------------------------------------

      do j=1,jmtr
      ltop(1,j)=ltop(imtr-1,j)
      ltop(imtr,j)=ltop(2,j)
      end do
#endif

c     do j=2,jmta-1
c     do i=2,imta-1
c     fx=0
c     do l=-1,1
c     do ll=-1,1
c     if(kmta(i+l,j+ll).eq.1)fx=1
c     end do
c     end do
c     if(kmta(i,j).eq.0.and.fx.eq.1)dd(i,j)=1
c     end do
c     end do
c     do j=1,jmta
c     write(95,45)(dd(i,j),i=1,imta)
c     end do
c  45 format(80i1)      
c     do j=1,9
c     do i=1,imta
c     if(dd(i,j).eq.1)write(95,46)i,j
c     end do
c     end do
c  46 format('      rivername(',i2,',',i2,')='ANTARCTICA   6'')
c     do j=10,jmta
c     do i=1,imta
c     if(dd(i,j).eq.1)write(95,45)i,j
c     end do
c     end do
c  45 format('      rivername(',i2,',',i2,')=')

c--------------------------------------------------------------
c  set time levels
c--------------------------------------------------------------

      nriv=1
      npriv=2


      return

      end



      subroutine rivmodel
      implicit none
#include "param.h"
#include "riv.h"
#include "grdvara.h"
#include "atm.h"
#include "cembm.h"
#include "levinda.h"

      real triver,qriver,areaubox,arealbox,
     &     drivn,drivne,drive,drivse,drivs
     & ,drivsw,drivw,drivnw,gradient,fx,vriver
      real tothriver,dtsriv,dh
      integer ii,jj,il,jl,ntemp


c------------------------------------------------------------
c  move dowm time levels and set time step
c------------------------------------------------------------

      ntemp=nriv
      nriv=npriv
      npriv=ntemp

      dtsriv = dtatm


c------------------------------------------------------------
c  initialise runoff array to zero
c------------------------------------------------------------

      do j=1,jmta
      do i=1,imta
      runoff(i,j)=0.0
      end do
      end do

c------------------------------------------------------------
c  set hriver(i,j,npriv) equal to hriver(i,j,nriv)
c  so that flow direction is based on nriv time level
c  and npriv time level is used to save new values of hriver
c------------------------------------------------------------

      do j=1,jmtr
      do i=1,imtr
      hriver(i,j,npriv)=hriver(i,j,nriv)
      end do
      end do

c------------------------------------------------------------
c  set cyclic points at next time level 
c  to zero water height so as to only
c  accumulate next time step's water
c------------------------------------------------------------

      do j=1,jmtr
      hriver(1,j,npriv)=0.0
      hriver(imtr,j,npriv)=0.0
      end do      

c------------------------------------------------------------
c  for each land grid box
c------------------------------------------------------------

      do j=2,jmtr-1
      do i=2,imtr-1
     
      if(ltop(i,j).gt.0)then

c------------------------------------------------------------
c  find lowest surrounding box including in water depth
c------------------------------------------------------------

      drive=(dxriv(i)+dxriv(i+1))*0.5*cstriv(j)
      drivw=(dxriv(i)+dxriv(i-1))*0.5*cstriv(j)
      drivn=(dyriv(j)+dyriv(j+1))*0.5
      drivs=(dyriv(j)+dyriv(j-1))*0.5
      drivne=sqrt(drive*drive+drivn*drivn)
      drivnw=sqrt(drivw*drivw+drivn*drivn)
      drivse=sqrt(drive*drive+drivs*drivs)
      drivsw=sqrt(drivw*drivw+drivs*drivs)

      gradient=(ltop(i,j)-ltop(i+1,j))/drive
      il=1
      jl=0

      fx=(ltop(i,j)-ltop(i-1,j))/drivw
      if(fx.gt.gradient)then
      il=-1
      jl=0
      gradient=fx
      endif

      fx=(ltop(i,j)-ltop(i,j+1))/drivn
      if(fx.gt.gradient)then
      il=0
      jl=1
      gradient=fx
      endif

      fx=(ltop(i,j)-ltop(i,j-1))/drivs
      if(fx.gt.gradient)then
      il=0
      jl=-1
      gradient=fx
      endif

      fx=(ltop(i,j)-ltop(i+1,j+1))/drivne
      if(fx.gt.gradient)then
      il=1
      jl=1
      gradient=fx
      endif
 
      fx=(ltop(i,j)-ltop(i+1,j-1))/drivse
      if(fx.gt.gradient)then
      il=1
      jl=-1
      gradient=fx
      endif

      fx=(ltop(i,j)-ltop(i-1,j-1))/drivsw
      if(fx.gt.gradient)then
      il=-1
      jl=-1
      gradient=fx
      endif

      fx=(ltop(i,j)-ltop(i-1,j+1))/drivnw
      if(fx.gt.gradient)then
      il=-1
      jl=1
      gradient=fx 
      endif

c-------------------------------------------------------------
c  convert gradient from m/cm to cm/cm
c-------------------------------------------------------------

      gradient=gradient*100.

c-------------------------------------------------------------
c  calculate velocity of river and set upper and lower limits
c-------------------------------------------------------------

      vriver=35.*sqrt(gradient/5.e-5)
      if(vriver.lt.15.)vriver=15.
      if(vriver.gt.50.)vriver=50.

c-------------------------------------------------------------
c  calculate time scale to drain grid box and volume loss
c-------------------------------------------------------------

      areaubox=dxriv(i)*dyriv(j)*cstriv(j)
      arealbox=dxriv(i+il)*dyriv(j+jl)*cstriv(j+jl)

      triver=sqrt(areaubox)/vriver
      if(triver.lt.abs(dtsriv))triver=abs(dtsriv)
      qriver=areaubox*hriver(i,j,nriv)/triver

c-------------------------------------------------------------
c  if lowest box is land drain grid box into its lowest neighbour
c-------------------------------------------------------------

      if(ltop(i+il,j+jl).gt.0.0)then

      hriver(i+il,j+jl,npriv)=hriver(i+il,j+jl,npriv) + dtsriv*
     &    qriver/arealbox    
      hriver(i,j,npriv)=hriver(i,j,npriv) - dtsriv*
     &    qriver/areaubox

c------------------------------------------------------------
c  else accumulate in runoff array with atmos model resolution 
c------------------------------------------------------------

      else

      hriver(i,j,npriv)=hriver(i,j,npriv) - dtsriv*
     &    qriver/areaubox

      ii=int((i+il-2)/dxyscale)+2
      jj=int((j+jl-1)/dxyscale)+1

      runoff(ii,jj)=runoff(ii,jj)+qriver*dtsriv 
      endif

c------------------------------------------------------------
c  add p-e
c------------------------------------------------------------

      ii=int((i-2)/dxyscale)+2
      jj=int((j-1)/dxyscale)+1

      hriver(i,j,npriv)=hriver(i,j,npriv)+
#ifdef landevap
     &  dtsriv*(precip(ii,jj)-evap(ii,jj))
#else
     &  dtsriv*precip(ii,jj)
#endif

c-----------------------------------------------------------
c  ensure hriver is not negative
c-----------------------------------------------------------

      if(hriver(i,j,npriv).lt.0.0)hriver(i,j,npriv)=0.0

c-----------------------------------------------------------
c  next point.........
c-----------------------------------------------------------

      endif

      end do
      end do

#ifdef cyclic
c-----------------------------------------------------------
c  set cyclic boundary conditions
c-----------------------------------------------------------

      do j=1,jmta
c     runoff(2,j)=runoff(2,j)+runoff(imta,j)
c     runoff(imta-1,j)=runoff(imta-1,j)+runoff(1,j)
      runoff(1,j)=runoff(imta-1,j)
      runoff(imta,j)=runoff(2,j)
      end do

      do j=1,jmtr
      hriver(imtr-1,j,npriv)=hriver(imtr-1,j,npriv)+hriver(1,j,npriv)
      hriver(2,j,npriv)=hriver(2,j,npriv)+hriver(imtr,j,npriv)
      hriver(1,j,npriv)=hriver(imtr-1,j,npriv)
      hriver(imtr,j,npriv)=hriver(2,j,npriv)
      end do
#endif
c martin
      tothriver=0.0
      do i=2,imtr-1 
      do j=2,jmtr-1
      ii=int((i-2)/dxyscale)+2
      jj=int((j-1)/dxyscale)+1
      if(kmta(ii,jj).eq.1)tothriver=tothriver+hriver(i,j,npriv)*
     &  dxriv(i)*dyriv(j)*cstriv(j)
      end do
      end do
      write(93,*)tothriver

      return
      end 

#endif

