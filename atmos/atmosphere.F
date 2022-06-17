      subroutine atmosphere
      implicit none
#include "param.h"
#include "atm.h"
#include "coorda.h"
#include "csbc.h"
#include "iounit.h"
#include "levinda.h"
#include "switch.h"
#include "ctmngr.h"
#include "mapsbc.h"
#include "cembm.h"
#include "scalar.h"
#include "grdvara.h"
#ifdef mrw_rivers
#include "riv.h"
#endif
#include "topog.h"
#include "levind.h"
#include "transform.h"
#ifdef uvic_embm_ice
#include "ice.h"
#endif
#include "embmdata.h"
      real xntspas,fx,runofftot
      integer mm,ntspas,loop,ix(imta)

c--------------------------------------------------------------------
c set up grid at start of run
c--------------------------------------------------------------------

      if(first)then
c#ifdef uvic_ice_cf
c      stop 'uvic_ice_cf not compatible with vargrid yet'
c#endif
c      call grida

c--------------------------------------------------------------------
c read in zonal mean wind stress, temp and humidity
c--------------------------------------------------------------------
#ifdef uvic_embm_astress
      stop 'read in of zonally averaged quantities
     &   may need modifiying if using variable grid'
      open(56,file='data/zaatemp0kam',status='old')
      do j=1,jmta
      read(56,*)zat(i,j)
      end do
      close(56)
      open(56,file='data/zahumid0kam',status='old')
      do j=1,jmta
      read(56,*)zah(i,j)
      end do
      close(56)
      open(56,file='data/zawindstressx0kam',status='old')
      do j=1,jmta
      read(56,*)zwx(i,j)
      end do
      close(56)
      open(56,file='data/zawindstressy0kam',status='old')
      do j=1,jmta
      read(56,*)zwy(i,j)
      end do
      close(56)
#endif
  
c-----------------------------------------------------------------------
c     read land mask kmt
c-----------------------------------------------------------------------

      do j=1,jmta
      do i=1,imta
      kmta(i,j)=kmt(i,j)
      end do
      end do
      do i=1,imta
      kmta(i,jmta)=0
      end do
      do j=1,jmta
      do i=1,imta
      if(kmta(i,j).gt.0)then
      fx=0
      else
      fx=1
      endif
      kmta(i,j)=fx
      end do
      end do
#ifdef cyclic
      do j=1,jmta
      kmta(1,j)=kmta(imta-1,j)
      kmta(imta,j)=kmta(2,j)
      end do
#endif

#ifdef uvic_embm_mountains
c-----------------------------------------------------------------------
c  read in elevation 
c-----------------------------------------------------------------------

      call transformtopog(3)
      do j=1,jmta
      do i=1,imta
      elev(i,j)=datamodel(i,j)*c100
      end do
      end do
#ifdef cyclic
      do j=1,jmta
      elev(1,j)=elev(imta-1,j)
      elev(imta,j)=elev(2,j)
      end do
#endif
      do j=1,jmta
      do i=1,imta
      if(kmta(i,j).ne.1)elev(i,j)=0.0
      if(kmta(i,j).eq.1.and.elev(i,j).lt.0.0)elev(i,j)=0.0
      end do
      end do
#endif
#ifdef mrw_rivers
c-----------------------------------------------------------------------
c     initialize the river model
c-----------------------------------------------------------------------

      call rivinit

#endif


#ifdef uvic_embm_ice
#ifdef coldstart
      if(init)then
      do i=1,imta
      do j=1,jmta
      print *,i,j,sbcocn(i,j,isst)
      if(sbcocn(i,j,isst).le.0.0)then
      do mm=1,3
      hice(i,j,mm)=200.0
      aice(i,j,mm)=1.0
      end do
      tice(i,j)=-30.
      endif
      end do
      end do
      endif
#endif
#endif
      endif

c
c
c-----------------------------------------------------------------------
c     Initialize S.B.C. names and related items.
c     
c     Dimensional units for the S.B.C. are for the model where they
c     will be used. eg: Heat flux will be used as a S.B.C. for the
c     ocean in cal/cm**2/sec.
c-----------------------------------------------------------------------

      do n=1,maxsbc
        mapsbc(n)  = 0
	coabc(n)   = c0
	crits(n)   = c0
	sbcname(n) = 'unknown'
	dunits(n)  = 'unknown'
      enddo
c
      if(first)then
      write (stdout,*) ' '
      write (stdout,*)
     & ' Surface Boundary Conditions are defined as follows:'
c
      write (stdout,*) ' '
      write (stdout,*) ' for the ocean (applied to MOM):'
      endif
c-------------------------------------------------------------------
c    define order of sbcs
c-------------------------------------------------------------------

      isu=-1
      isv=0
      itaux=1
      itauy=2
      ihflx=3
      ifflx=4
      isst =5
      isss =6
c

      do n=1,numsbc
c

c       order S.B.C.'s such that S.B.C. "n" has index "m"=mapsbc(n)
c       within arrays "sbcocn" and "sbcatm"
c
        mapsbc(n) = n
        m = mapsbc(n)
	if (m .eq. 1) then
	  sbcname(m) = ' taux     '
	  dunits(m)  = ' dynes/cm**2   '
	  coabc(m)   = 1.0
	  crits(m)   = 1.e-3
	else if (m .eq. 2) then
	  sbcname(m) = ' tauy     '
	  dunits(m)  = ' dynes/cm**2   '
	  coabc(m)   = 1.0
	  crits(m)   = 1.e-3
	else if (m .eq. 3) then
	  sbcname(m) = ' heat flux'
	  dunits(m)  = ' cal/cm**2/sec '
	  coabc(m)   = 1.0
	  crits(m)   = 1.e-5
	else if (m .eq. 4) then
	  sbcname(m) = ' salt flux'
	  dunits(m)  = ' gm/cm**2/sec  '
	  coabc(m)   = 1.0
	  crits(m)   = 1.e-6
	else if (m .eq. 5) then
	  sbcname(m) = ' SST      '
	  dunits(m)  = ' deg C         '
	  coabc(m)   = 0.0
	  crits(m)   = 0.05
	else if (m .eq. 6) then
	  sbcname(m) = ' SSS      '
	  dunits(m)  = ' (S-35.0)/1000 '
	  coabc(m)   = 1.0
	  crits(m)   = 1.e-5
	endif
c
      if(first)then
        if (n .eq. numosbc+1) then
          write (stdout,*) ' '
          write (stdout,*) ' for the atmosphere (from MOM):'
	endif
	if (n .le. numosbc) then
	  write (stdout,*) ' S.B.C. # ',m,' is ',sbcname(m)
     &,', atmos native units are ', dunits(m)

        else
	  write (stdout,*) ' S.B.C. # ',m,' is ',sbcname(m)
     &,', MOM native units are ', dunits(m)

	endif
        endif
      enddo


c-----------------------------------------------------------------------
c     do the introductory atmosphere setup once per run
c-----------------------------------------------------------------------
c     "setatm" must do the following:
c       1) set up the atmospheric S.B.C. grid definition
c       2) define the atmosphere land/sea mask
c       3) set the atmosphere time step "dtatm" {seconds}
c----------------------------------------------------------------------
      call setembm (dtatm)

c-----------------------------------------------------------------------
c     compute the number of atmospheric time steps "numots" for this run and 
c     in the segment 
c-----------------------------------------------------------------------

      xntspas =seglength/dtatm
      ntspas=seglength/dtatm
      if(xntspas-ntspas.ne.0.0)
     &stop 'atmos segment length not multiple of atmospheric timestep' 


c-----------------------------------------------------------------------
c     S T A R T    S E G M E N T    L O O P
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
c       get the atmospheric S.B.C. from MOM
c-----------------------------------------------------------------------
        call gasbc

c
c-----------------------------------------------------------------------
c       call the atmospheric model once for each time step until one
c       segment of "segtim" days is complete. hold atmos S.B.C. fixed
c       during each segment and predict average S.B.C. for MOM
c-----------------------------------------------------------------------

        do loop=1,ntspas
          call embm
        enddo

c----------------------------------------------------------------------
c  calculate toa energy balance
c----------------------------------------------------------------------

        do j=1,jmta
        do i=1,imta
        toatot(i,j)=toasolar(i,j)-outlwr(i,j)
        end do
        end do
c martin
#ifdef mrw_rivers
        runofftot=0.0
        do j=1,jmta
        do i=1,imta
        runofftot=runofftot+runoff(i,j)
        end do
        end do
        write(92,*)itt,dts,runofftot
        
c-----------------------------------------------------------------------
c       update wind stress anomalies
c-----------------------------------------------------------------------

	call windstress
# endif

c-----------------------------------------------------------------------
c     E N D    S E G M E N T    L O O P
c-----------------------------------------------------------------------


8800  format (/,/,1x,'The model has been configured to run for'
     &,g14.7,' days in ',i10,' segments of ',g14.7,' days each.'
     &,/1x,'There will be ',i6,' ocean time steps per segment'
     &,' and ',i6,' atmosphere time steps per segment.'/
     &,/1x,' The ocean "dts" =', g14.7, ' seconds'
     &,', and the atmosphere "dtatm" =', g14.7,' seconds'
     &,/,/)
9000  format(/,/,1x,'==== E N T E R I N G    ',a9,'    S E G M E N T  #'
     &,       i4, ',  time = ',a32,' ===='/,/)
9100  format (10x,'==> Finished reading MOM S.B.C. data from file'
     &,    ' restrt_sbc.dta at time = ',a32)
9200  format (10x,'==> Finished writing MOM S.B.C. data to file '
     &,    'restrt_sbc.dta at time = ',a32)
c      stop
      return
      end


#include "setembm.F"
#include "gasbc.F"
#include "embm.F"
#include "embmio.F"
c#include "blkdtaa.F"
#include "grida.F"
#include "insolation.F"
#include "solve.F"
#include "slap.F"
#include "couple.F"
#include "rivmodel.F"
#include "fluxes.F"
#include "embmbc.F"
#include "icemodel.F"
#include "winds.F"
