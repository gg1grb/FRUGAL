      subroutine insolation (phia, solins)
      implicit none
c=======================================================================
c     calculate daily insolation
c     based on code by a.fanning and a.berger
c
c     input:
c       dayoyr = days + fractional days since start of calendar year
c       phi    = latitude in degrees
c     output:
c       solins = insolation (units as for solarconst)
c
c     author:   m.eby   e-mail: eby@uvic.ca
c=======================================================================
c

# include "insolation.h"
# include "scalar.h"
# include "cembm.h"
# include "ctmngr.h"
# include "iounit.h"
      real xl,xllp,xee,xse,xlam,dlamm,anm,ranm,xec,ranv,
     & anv,tls,rphi,s,sind,cosd,rdelta,aphi,adelta,dayl,
     & ww,phid,sinp,cosp,tp,stp,rdayl,phia,delta,solins,
     & t,xes,arg,rp,prg

      nd = nint(dayoyr)

c--------------------------------------------------------------------
c     set longitude of perihelion relative to vernal equinox
c--------------------------------------------------------------------

      xl = perh + 180.0

c--------------------------------------------------------------------
c     calendar date  month (ma) and day (ja)
c     nd  number of this day in a year of 365 days
c     xlam = mean long. sun for true long. = 0
c     dlamm = mean long. sun for ma-ja
c--------------------------------------------------------------------

      xllp = xl*pir
      xee = ecc**2.0
      xse = sqrt(1.0-xee)
      xlam = (ecc/2.0 + ecc*xee/8.0)*(1.0+xse)*sin(xllp)
     &     - xee/4.0*(1.5 + xse)*sin(2.0*xllp)
     &     + ecc*xee/8.0*(1.0/3.0 + xse)*sin(3.20*xllp)
      xlam = 2.0*xlam/pir
      dlamm = xlam + (nd - 80)*step
      anm = dlamm - xl
      ranm = anm*pir
      xec = xee*ecc
      ranv = ranm + (2.0*ecc - xec/4.0)*sin(ranm)
     &     + 5.0/4.0*ecc**2*sin(2.0*ranm)
     &     + 13.0/12.0*xec*sin(3.0*ranm)
      anv = ranv/pir

c-------------------------------------------------------------------
c     true longitude of the earth
c-------------------------------------------------------------------

      tls = anv + xl

      rphi = phia*pir
      s = 1.0/(pi*((1.0-ecc*ecc)/(1.0+ecc*cos((tls-xl)*pir)))**2)
      sind = sin(xob*pir)*sin(tls*pir)
      cosd = sqrt(1.0-sind**2)
      rdelta = atan(sind/cosd)
      delta = rdelta/pir
      aphi = abs(phia)
      adelta = abs(delta)
c      if (years.ge.5064.9) then
c#ifdef insolation_vary
c        open(71,file='insolv.1',status='replace')
c#else
c        open(71,file='insolc.1',status='replace')
c#endif
c      endif
c------------------------------------------------------------------
c     singularity for aphi = 90 and delta = 0 => polar night
c------------------------------------------------------------------

      if ((abs(aphi - 90.0) .gt. test) .or. (adelta .gt. test)) then
        if(adelta .le. test) then
          dayl = 12.0
          ww = s*cos(rphi)
c          if (years.ge.5064.9) then
c         write(71,*) '2',ww,s,rphi
c          endif
        else
          if(aphi .le. test) then
            dayl = 12.0
            ww = s*cos(rdelta)
c             if (years.ge.5064.9) then
c            write(71,*) '3',ww,s,rdelta
c            endif
          else
            phid = phia*delta
            if ((aphi .gt. 90.0 - adelta) .and. (phid .ne. 0.0)) then
              if (phid .gt. 0.0) then
c               polar continual day
                dayl = 24.00
	        ww = s*sind*sin(rphi)*pi
c                 if (years.ge.5064.9) then
c                 write(71,*) '4',ww,s,rphi,sind,pi
c                 endif
              else
c               polar continual night
                dayl = 0.0
                ww = 0.0
c                 if (years.ge.5064.9) then
c                 write(71,*) '5',ww
c                 endif
              endif
            else
c             daily sunrise and sunset
              sinp = sind*sin(rphi)
              cosp = cosd*cos(rphi)
              tp = -sinp/cosp
              if (tp.gt.1.) then
                tp = 1.0
              elseif (tp.lt.-1.) then
                tp = -1.0
              endif
              stp = sqrt(1.0-tp*tp)
              rdayl = acos(tp)
              dayl = 24.0*rdayl/pi
              ww = s*(rdayl*sinp+cosp*stp)
            endif
          endif
        endif
      else
c       polar continual night
        dayl = 0.0
        ww = 0.0
      endif

      solins = solarconst*ww


      return
      end


      subroutine berger
      implicit none
c=======================================================================
c     calculate eccentricity, precession and obliquity for any year
c     based on code by a.fanning and a.berger
c
c     author:   m.eby   e-mail: eby@uvic.ca
c=======================================================================

#include "param.h"
# include "insolation.h"
# include "scalar.h"
# include "cembm.h"
#include "iounit.h"
#include "switch.h"


      real t,xes,xec,arg,prg,rp
      logical done

      pir = pi/180.0
      pirr = pir/3600.0
      step = 360.0/365.25
      test = 1.0e-4
      xod = 23.320556
      xop = 3.392506
      prm = 50.439273

c-------------------------------------------------------------------      
c     eccentricity                                       
c     number of terms to keep for expansion solution of eccentricity
c-------------------------------------------------------------------
       
      do i=1,nef
        be(i) = ye(i)*pirr
        ce(i) = ze(i)*pir
      enddo

c-------------------------------------------------------------------
c     obliquity relative to mean ecliptic of date
c     number of terms to keep for expansion solution of obliquity
c-------------------------------------------------------------------

      do i=1,nob
        bob(i) = yob(i)*pirr
        cob(i) = zob(i)*pir
      enddo

c-------------------------------------------------------------------
c     general precession in longitude
c     number of terms to keep in series expansion of solution of 
c     general precession
c-------------------------------------------------------------------

      do i=1,nop
        bop(i) = yop(i)*pirr
        cop(i) = zop(i)*pir
      enddo

c------------------------------------------------------------------
c     Berger's solution is referenced to 1950.
c     convert t such that zero is referenced about calendar year 0
c------------------------------------------------------------------

      t = pyear - 1950
c      write(6,*) 'in berger:t=',t,pyear
c-----------------------------------------------------------------------
c      calculate eccentricity (ecc)
c-----------------------------------------------------------------------

      xes = 0.0
      xec = 0.0
      do i=1, nef
        arg = be(i)*t + ce(i)
        xes = xes + ae(i)*sin(arg)
        xec = xec + ae(i)*cos(arg)
      enddo

      ecc = sqrt(xes*xes + xec*xec)

c-----------------------------------------------------------------------
c     calculate precession (pre)
c-----------------------------------------------------------------------

      if (abs(xec) .gt. 1.0e-08) then
        rp = atan(xes/xec)
        if (xec .lt. 0.0) rp = rp + pi
        if (xec .gt. 0.0 .and. xes .lt. 0.0) rp = rp + 2.0*pi
      else
        if (xes .lt. 0.0) rp = 1.5*pi
        if (xes .eq. 0.0) rp = 0.0
        if (xes .gt. 0.0) rp = pi/2.0
      endif
      perh = rp/pir

      prg = prm*t
      do i=1, nop
        arg = bop(i)*t + cop(i)
        prg = prg + aop(i)*sin(arg)
      enddo
      prg = (prg/3600.0) + xop
      perh = perh + prg
      if (abs(perh) .gt. 360.) perh = perh-float(int(perh/360.0))*360.0
      if (perh .lt. 0.0) perh = perh + 360.0

      pre = ecc*sin(perh*pir)

c-----------------------------------------------------------------------
c     calculate obliquity (xob)
c-----------------------------------------------------------------------

      xob = xod
      do i=1,nob
        arg = bob(i)*t + cob(i)
        xob = xob+aob(i)/3600.0*cos(arg)
      enddo

c-----------------------------------------------------------------------
c     write summary
c-----------------------------------------------------------------------

      if(first)then
      write(stdout,*) 
      write(stdout,*) 'Insolation parameters'
      write(stdout,*) 'year = ', t + 1950, '(not always model year)'
      write(stdout,*) 'eccentricity = ', ecc
      write(stdout,*) 'obliquity = ', xob
      write(stdout,*) 'precession factor = ', perh+180.0
      write(stdout,*) 'precession = ', pre
      endif

	return
	end
