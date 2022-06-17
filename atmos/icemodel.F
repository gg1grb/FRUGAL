#ifdef uvic_embm_ice
      subroutine therm
c=======================================================================
c     thermodynamic ice model
c
c     Note: if run with embm this routine must be called after "fluxes" 
c           and before "solve"
c
c     calculates ice and open water growth rates based on surface 
c     energy budget. see Parkinson and Washington, JGR, Vol.84, C1, 
c     311-337, 1979 and Hibler, JPO, Vol.9, 815-846, 1979
c
c     heat and fresh water fluxes between the ocean and atmosphere 
c     are adjusted depending on ice growth or melt. ice thickness is 
c     changed by the amount of growth or melt and sublimation
c
c     author:   m.eby   e-mail: eby@uvic.ca
c=======================================================================
c
# include "param.h"
# include "csbc.h"
# include "cembm.h"
# include "atm.h"
# include "ice.h"
# include "levind.h"
# include "grdvar.h"
# include "mapsbc.h"
# include "iounit.h"
#include "fields.h"
#include "timelv.h"
c
      real epsln,fa,fb,fc,fd,fe,ff,fg,fh,fi,ho,amin,amax,scale,tol,
     &     errmax,dswr,ai,ftopo,fbot,ao,dho,ws,tcdh,ti,tiold,
     &     dalt,fpt,fm,fn,fl,dusens,dt,qi,ultnt,dultnt,usens,ulwr,f,
     &     df,delta,sub,ftopi,dulwr,dha,dhi,dh,hmax
      integer mm,maxit,index,imaxi,jmaxi,iter

      epsln=1.e-20
      fa = dts/(rhoice*flice)
      fb = 0.94*rhoatm*cpatm
      fc = 3.80116e3/ssp
      fd = rhoatm/rhoice
      fe = rhoatm*slice
      ff = rhoice*flice/41859.
      fg = rhoice/rhoocn
      fh = 21.8746*265.5
      ho = 1.0 
      amin = 0.15
      amax = 0.99
      scale = 10.0
      maxit = 10
      tol = 0.01
      errmax = 0.0
      imaxi = 0
      jmaxi = 0
# ifdef uvic_embm_ice_leads
      amax = 0.99
      hmax = 1000.0
# endif
# ifdef uvic_ice_cf
      index = 3
# else
      index = lf
# endif
c
      do j=2,jmtm1
        do i=2,imtm1
c      if(kmt(i,j).gt.0)then

        dswr = (1.0 - scatter(i,j))*solins(i,j)*coalice(i,j)
c
          ai = aice(i,j,2)

	  if (kmt(i,j) .gt. 0) then
          ao = 1.0 - ai
c
c-----------------------------------------------------------------------
c           calculate fluxes to and from the ocean (without ice)
c-----------------------------------------------------------------------
c
	    ftopo = dnswr(i,j) - uplwr(i,j) - upsens(i,j) - upltnt(i,j)
c           fbot = sbcocn(i,j,ihflx)/2.389e-8
            fbot = sbcocn(i,j,ihflx)/1.e-3
c           calculate growth of ice in open water areas
            dho = fa*(fbot - ftopo)
c
c martin
c           if (ai .ne. 0.0)then
            if (ai .ne. 0.0.or.(t(1,i,j,1,nct(i,j)).le.frzpt(i,j).and.
     &      dho.gt.0.0)) then
            if(ai.eq.0.0)then
            ai=0.1
            do mm=1,3
            hice(i,j,mm)=1.
            end do
            endif
            if(ai.gt.c0.and.hice(i,j,2).eq.c0)then
            do mm=1,3
            hice(i,j,mm)=1.
            end do
            endif
c
c-----------------------------------------------------------------------
c             find ice temperature by balancing the surface heat budget:
c               tcdh*(ti - fpt) = dswr - ultnt - usens - ulwr
c             using Newton's method:
c               t(i+1) = t(i) - f(t(i))/df(t(i)) 
c             where:
c               f(t(i)) = dswr - ultnt - usens - ulwr - tcdh*(ti - fpt)
c               -df(t(i)) = dultnt + dusens + dulwr + tcdh
c-----------------------------------------------------------------------
c
              ws = spd(i,j)
# ifdef uvic_embm_astress
     &           + aspd(i,j)
# endif
              tcdh = condice/hice(i,j,2)
              ti = tice(i,j) 
              tiold = tice(i,j)

              dalt = 1.0022e-3+8.22e-5*(ti-at(i,j,2,1))+2.66e-6*ws
              dalt =  max(dalt, 6.0e-05)
              dalt =  min(dalt, 2.190e-03)
#ifdef xlth
              if(dalt .lt. 2.19e-03)then
              fi=fe*1.5
              else
              fi=fe
              endif
#endif

              fpt = frzpt(i,j)
# ifdef uvic_embm_clouds
              rh =  at(i,j,2,2)/ssh(i,j)
              dulwr = 4.0*esice*(at(i,j,2,1) + 273.15)**3
              cf = 11.5 - 22.8*rh + 11.7*rh**2
c prevent negative implied could cover
c           cf=amax1(cf,0.4)
              fm = esice*(at(i,j,2,1) + 273.15)**4 
     &             *(0.39 - 0.05*(at(i,j,2,2)*1.61e3)**.5)*cf
# else
              fm = esatm(i,j)*(at(i,j,2,1) + 273.15)**4
              fn = 4.0*esice
# endif
c
#ifdef xlth
                  fl = fi*dalt*ws
#else
                  fl = fe*dalt*ws
#endif    
                 dusens = fb*dalt*ws
                  do iter=1,maxit
c      write(97,*)iter,hice(i,j,3),ti,at(i,j,2,1),fc
               if(ti.gt.0.0.or.ti.lt.-80.)goto 889
                    dt = ti - at(i,j,2,1)
	            qi = fc*exp(21.8746*ti/(ti + 265.5))
                    ultnt = fl*(qi - at(i,j,2,2))
                    dultnt = fl*qi*fh/(ti + 265.5)**2
                    usens = dusens*dt
# ifdef uvic_embm_clouds
                    ulwr = fm + dulwr*dt
# else
                    ulwr = esice*(ti + 273.15)**4 - fm
                    dulwr = fn*(ti + 273.15)**3
# endif
                    f = dswr - ultnt - usens - ulwr - tcdh*(ti - fpt)
                    df = dultnt + dusens + dulwr + tcdh
                    delta = f/df
                    ti = ti + delta
                    if (abs(delta) .le. tol) goto 100
                  enddo
c                 if temperature below freezing and did not converge 
c                 set to last converged temperature
                  if (abs(delta) .gt. errmax .and. ti .lt. fpt) then
                    errmax = abs(delta)
                    imaxi = i
                    jmaxi = j
                    ti = tiold
                  endif
c
c-----------------------------------------------------------------------
c                 set maximum tice to freezing and calculate fluxes
c-----------------------------------------------------------------------
c
  889        ti=tiold
100               ti = min(ti, fpt)
                  dt = ti - at(i,j,2,1)
	          qi = fc*exp(21.8746*ti/(ti + 265.5))
                  sub = max(0.0, dalt*ws*(qi - at(i,j,2,2)))
#ifdef xlth
                  ultnt = fi*sub
#else
                  ultnt = fe*sub
#endif
                  sub = fd*sub
                  usens = dusens*dt
# ifdef uvic_embm_clouds
                  ulwr = fm + dulwr*dt
# else
                  ulwr = esice*(ti + 273.15)**4 - fm
# endif
c
c-----------------------------------------------------------------------
c                 add ice covered area fluxes to ocean area fluxes
c-----------------------------------------------------------------------
c
                  tice(i,j) = ti
c martin
                  if(ai.eq.0.0)tice(i,j)=c0

                  dnswr(i,j) = dnswr(i,j)*ao + dswr*ai
                  upltnt(i,j) = upltnt(i,j)*ao + ultnt*ai
                  upsens(i,j) = upsens(i,j)*ao + usens*ai
                  uplwr(i,j) = uplwr(i,j)*ao + ulwr*ai
                  ftopi = dswr - ulwr - usens - ultnt
                
c
c-----------------------------------------------------------------------
c                 calculate change in ice volume due to sublimation of  
c                 ice (dha). adjust evaporation to the atmosphere to  
c                 account for sublimation from ice. subtract this 
c                 adjustment from the ocean freshwater flux
c-----------------------------------------------------------------------
c
                  dha = -dts*sub
                  sub = sub*ai*fg
c martin
c                  flux(i,j,2) = flux(i,j,2) + dts*sub
c
                  evap(i,j) = evap(i,j)*ao + sub
                else
                  ftopi = 0.0
                  dha = 0.0
                endif
c
c-----------------------------------------------------------------------
c           calculate total change in ice volume (dh)
c-----------------------------------------------------------------------
c
c           total growth of ice from the ocean
            dhi = ai*fa*(fbot - ftopi) + ao*dho
c           total growth (loss limited to total amount)
            dha = dha * ai
            dh = -min(hice(i,j,index), -(dhi + dha))
c
c-----------------------------------------------------------------------
c           calculate new area and thickness from thermodynamics
c-----------------------------------------------------------------------
c
c           use minimum area (amin) of open water
            ai = max(amin, aice(i,j,index))
            aice(i,j,3) = aice(i,j,index) + ((1.0 - ai)*max(0.0, dho)/ho
     &         + 0.5*min(0.0, dhi)*ai/(hice(i,j,index) + epsln))
            hice(i,j,3) = hice(i,j,index) + dh
c
# ifdef uvic_embm_ice_leads
c           leads parameterization
            ai = amax + hice(i,j,3)*(1.0 - amax)/hmax
            aice(i,j,3) = min(aice(i,j,3), ai)
# endif
            aice(i,j,3) = max(0.0, min(1.0, aice(i,j,3)))
c           zero ice area where thickness is zero (use h scaled to mm)
            aice(i,j,3) = min(aice(i,j,3), hice(i,j,3)*scale)
c
c-----------------------------------------------------------------------
c           adjust fluxes to the ocean due to ice melt or growth
c-----------------------------------------------------------------------
c
            dh = dh - dha
            flux(i,j,1) = flux(i,j,1) + ff*dh
c martin
c            flux(i,j,2) = flux(i,j,2) - fg*dh
c
          else
c
c-----------------------------------------------------------------------
c           adjust shortwave flux for land ice and set new ice area
c-----------------------------------------------------------------------
c
            ai = aice(i,j,2)
           dnswr(i,j) = dnswr(i,j)*(1.0 - ai) + dswr*ai
#ifdef prescribelandice
          aice(i,j,3)=icemask(i,j)
#else
# ifdef uvic_embm_mountains
           if (at(i,j,2,1)-elev(i,j)/14286. .le. -10.0) then
# else
           if (at(i,j,2,1) .le. -10.0) then
# endif
             aice(i,j,3) = 1.0
           else
             aice(i,j,3) = 0.0
           endif
#endif
c
	  endif
c
c endif ocean only
c      endif

        enddo
      enddo
c
      if (errmax .gt. 0.0) write (stdout,*) 
     &  '==> Warning: ice temperature not converging: errmax, i, j:'
     &, errmax, imaxi, jmaxi
#if mrw_ice_adv
c
c-----------------------------------------------------------------------
c     advect ice thickness and area with last ocean velocities
c-----------------------------------------------------------------------
c
      call advice (hice)
      call advice (aice)
#endif
c limit ice thickness to 4m
      do j=1,jmt
      do i=1,imt
      hice(i,j,3)=amin1(hice(i,j,3),400.)
      end do
      end do
c
c-----------------------------------------------------------------------
c     set boundary conditions
c-----------------------------------------------------------------------
c
      call setbcx (hice(1,1,3), imt, jmt)
      call setbcx (aice(1,1,3), imt, jmt)
c
c-----------------------------------------------------------------------
c     shuffle time levels
c-----------------------------------------------------------------------
c
      do j=1,jmt
        do i=1,imt
          hice(i,j,1) = hice(i,j,2)
          hice(i,j,2) = hice(i,j,3)
          aice(i,j,1) = aice(i,j,2)
          aice(i,j,2) = aice(i,j,3)
        enddo
      enddo
c
      return
      end

#ifdef uvic_ice_cf
      subroutine cavfld
c=======================================================================
c     calculate ice velocities using a "cavitating fluid" rheology
c     based on code by g.flato
c
c     see Flato and Hibler, JPO, Vol.22, 626-651, 1992
c
c     author:   m.eby   e-mail: eby@uvic.ca
c=======================================================================
c
#  include "param.h"
#  include "cembm.h"
#  include "csbc.h"
#  include "ice.h"
#  include "grdvara.h"
#  include "levind.h"
#  include "mapsbc.h"
#  include "atm.h"     
c
      real cw(imt,jmt), Ax(imt,jmt), Ay(imt,jmt), Bx(imt,jmt)
      real By(imt,jmt), XX(imt,jmt), YY(imt,jmt), pmax(imt,jmt)
      real epsln
c
      integer jsi(jmt), jei(jmt)
c
c     cw    = drag coefficient between ice and water
c     Ax, Ay, Bx, By, XX, YY = momentum equation coefficients
c     pmax  = maximum ice strength
c     jsi   = starting indices for ice segments
c     jei   = ending indices for ice segments
c
c     pstar = constant used to calculate maximum strength (g/s**2/cm)
c     sinth = sin of turning angle (25 degrees)
c     costh = cos of turning angle (25 degrees)
c
      epsln=1.e-20
      pstar = 2.75e5
      sinth = 0.4226
      costh = 0.9063
      cwmin = 0.055
c
c     set convergence criteria and maximum iterations for solvers
      tolfd = 0.01
      mxitfd = 150
      tolcf = 0.5
      mxitcf = 1500
c
      cdrho = cdice*rhoocn
c
c-----------------------------------------------------------------------
c     calculate starting and ending latitudes for velocity calculations
c-----------------------------------------------------------------------
c
      nseg = 0
      inc = 1
      do j=2,jmtm1
        ice = 0
        jm2 = max(j-2,1)
        jp2 = min(j+2,jmt)
        do i=2,imtm1
          hi = hice(i,jm2,2) + hice(i,j-1,2) + hice(i,j,2) 
     &       + hice(i,j+1,2) + hice(i,jp2,2)
          if (hi .ne. 0.0) ice = 1
        enddo
        if (ice .eq. 1) then
          nseg = nseg + inc
          if (inc .eq. 1) jsi(nseg) = j
          if (j .eq. jmtm1) jei(nseg) = jmtm1
          inc = 0 
        else
          do i=1,imt
            uice(i,j,2) = 0.0
            vice(i,j,2) = 0.0
          enddo
          if (inc .eq. 0) jei(nseg) = j
          inc = 1
        endif
      enddo
      if (nseg .eq. 0) return
c
c-----------------------------------------------------------------------
c     calculate compressive strength and zero initial arrays
c-----------------------------------------------------------------------
c
      do j=1,jmt
        do i=1,imt
c         calculate strength and zero pressure
          pmax(i,j) = pstar*exp(-20.0*(1.0 - aice(i,j,2)))*hice(i,j,2)
          presur(i,j) = 0.0
          uice(i,j,1) = 0.0
          vice(i,j,1) = 0.0
          cw(i,j) = 0.0
          Ax(i,j) = 0.0
          Ay(i,j) = 0.0
          Bx(i,j) = 0.0
          By(i,j) = 0.0
          XX(i,j) = 0.0
          YY(i,j) = 0.0
        enddo
      enddo
c
      do kme=1,2
c
c-----------------------------------------------------------------------
c       average velocities after the first Euler step
c-----------------------------------------------------------------------
c
        if (kme .ge. 2) then
c         set different convergence criteria after the first pass
c          tolfd = 0.001
c          tolcf = 0.1
          do n=1,nseg
            do j=jsi(n),jei(n)
              do i=2,imtm1
                uice(i,j,2) = 0.5*(uice(i,j,1) + uice(i,j,2))
                vice(i,j,2) = 0.5*(vice(i,j,1) + vice(i,j,2))
              enddo
            enddo
          enddo
        endif
c
c-----------------------------------------------------------------------
c       calculate drag between the ice and ocean (averaged to B-grid)
c-----------------------------------------------------------------------
c
        do n=1,nseg
          do j=jsi(n),jei(n)
            do i=2,imtm1
c             save velocity for averaging on the correction step
              uice(i,j,1) = uice(i,j,2)
              vice(i,j,1) = vice(i,j,2)
              dub = 0.5*(uice(i,j,2) + uice(i,j+1,2)) - sbcocn(i,j,isu)
              dvb = 0.5*(vice(i,j,2) + vice(i+1,j,2)) - sbcocn(i,j,isv)
              cw(i,j) = cdrho*sqrt(dub**2 + dvb**2)
            enddo
          enddo
        enddo
        call setbcx (uice(1,1,1), imt, jmt)
        call setbcx (vice(1,1,1), imt, jmt)
        call setbcx (cw, imt, jmt)
c       set minimum drag coefficient
        do j=1,jmt
          do i=1,imt
            cw(i,j) = max(cwmin, cw(i,j))
          enddo
        enddo
c
c-----------------------------------------------------------------------
c       calculate forces and free drift velocities on B-grid
c-----------------------------------------------------------------------
c
        do n=1,nseg
          do j=jsi(n),jei(n)
            do i=2,imtm1
              gwu = sbcocn(i,j,igu)
              gwv = sbcocn(i,j,igv)
              Aij = cw(i,j)*costh
              Bij = 0.25*(hice(i+1,j+1,2) + hice(i,j+1,2) + hice(i,j,2)
     &            + hice(i+1,j,2))*rhoice*cori(i,j,1) + cw(i,j)*sinth
              XX(i,j) = Aij*gwu - Bij*gwv + sbcocn(i,j,itaux)
              YY(i,j) = Aij*gwv + Bij*gwu + sbcocn(i,j,itauy)
              detAB = (Aij)**2 + (Bij)**2
              if (abs(detAB) .gt. epsln) then
                uice(i,j,2) = (Aij*XX(i,j) + Bij*YY(i,j))/detAB
                vice(i,j,2) = (Aij*YY(i,j) - Bij*XX(i,j))/detAB
              else
                uice(i,j,2) = 0.0
                vice(i,j,2) = 0.0
              endif
            enddo
          enddo
        enddo
        call setbcx (XX, imt, jmt)
        call setbcx (YY, imt, jmt)
        call setbcx (uice(1,1,2), imt, jmt)
        call setbcx (vice(1,1,2), imt, jmt)
c
c-----------------------------------------------------------------------
c       average forces and velocities to C-grid
c-----------------------------------------------------------------------
c
        do n=1,nseg
          do j=jsi(n),jei(n)
            do i=2,imtm1
              cwx = 0.5*(cw(i,j-1) + cw(i,j))
              Ax(i,j) = cwx*costh
              Bx(i,j) = 0.25*(hice(i,j,2) + hice(i+1,j,2))
     &            *rhoice*(cori(i,j-1,1) + cori(i,j,1)) + cwx*sinth
              cwy = 0.5*(cw(i-1,j) + cw(i,j))
              Ay(i,j) = cwy*costh
              By(i,j) = 0.5*(hice(i,j,2) + hice(i,j+1,2))
     &                  *rhoice*cori(i,j,1) + cwy*sinth
            enddo
          enddo
          do j=jei(n),jsi(n),-1
            do i=imtm1,2,-1
              XX(i,j) = 0.5*(XX(i,j-1) + XX(i,j))
              YY(i,j) = 0.5*(YY(i-1,j) + YY(i,j))
              uice(i,j,2) = 0.5*(uice(i,j-1,2) + uice(i,j,2))*uvmx(i,j)
              vice(i,j,2) = 0.5*(vice(i-1,j,2) + vice(i,j,2))*uvmy(i,j)
             enddo
          enddo
        enddo
        call setbcx (Ax, imt, jmt)
        call setbcx (Ay, imt, jmt)
        do j=1,jmt
          do i=1,imt
            Ax(i,j) = max(cwmin*costh, Ax(i,j))
            Ay(i,j) = max(cwmin*costh, Ay(i,j))
          enddo
        enddo
        call setbcx (XX, imt, jmt)
        call setbcx (YY, imt, jmt)
        call setbcx (uice(1,1,2), imt, jmt)
        call setbcx (vice(1,1,2), imt, jmt)
c
c-----------------------------------------------------------------------
c       calculate free drift velocity on C-grid by under-relaxation
c-----------------------------------------------------------------------
c
        w = 0.25
        do kfd=1,mxitfd
          emax = 0.0
          do n=1,nseg
            do j=jsi(n),jei(n)
              do i=2,imtm1
                dpdx = (presur(i+1,j) - presur(i,j))*dxura(i,j)
                Bvbar = 0.25*Bx(i,j)*(vice(i,j,2) + vice(i+1,j,2)
     &                + vice(i,j-1,2) + vice(i+1,j-1,2))
                du = (uice(i,j,2) + (dpdx - XX(i,j) - Bvbar)/Ax(i,j))
     &               *w*uvmx(i,j)
                uice(i,j,2) = uice(i,j,2) - du
                emax = max(abs(du), emax)
              enddo
            enddo
          enddo
          call setbcx (uice(1,1,2), imt, jmt)
          do n=1,nseg
            do j=jsi(n),jei(n)
              do i=2,imtm1
                dpdy = (presur(i,j+1) - presur(i,j))*dyura(i,j)
                Bubar = 0.25*By(i,j)*(uice(i-1,j,2) + uice(i,j,2)
     &                + uice(i-1,j+1,2) + uice(i,j+1,2))
                dv = (vice(i,j,2) + (dpdy - YY(i,j) + Bubar)/Ay(i,j))
     &               *w*uvmy(i,j)
                vice(i,j,2) = vice(i,j,2) - dv
                emax = max(abs(dv), emax)
              enddo
            enddo
          enddo
          call setbcx (vice(1,1,2), imt, jmt)
c         compare error with tolerance
          if (emax .le. tolfd) go to 10
        enddo
        if (kme .ne. 1)
     &    write(*,*) '=> Warning: iterations exceeded in calculating',
     &    ' free drift. iterations, maximum error:', kfd, emax
c
c-----------------------------------------------------------------------
c       calculate cavitating fluid correction to free drift velocity
c-----------------------------------------------------------------------
c
10      do n=1,nseg
          do j=jsi(n),jei(n)
            do i=2,imtm1
c             reuse arrays XX and YY to hold updated velocities
              XX(i,j) = uice(i,j,2)
              YY(i,j) = vice(i,j,2)
            enddo
          enddo
        enddo
        do kcf=1,mxitcf
          emax = 0.0
          do n=1,nseg
            do j=jsi(n),jei(n)
              do i=2,imt
                if (kmt(i,j) .gt. 0) then
c                 calculate divergence times cst(j) (leave cst in 
c                 divergence since it is needed to calculate dp)
               D = (vice(i,j,2)*csua(i,j) - vice(i,j-1,2)*csua(i,j-1))
     &         *dytra(i,j) + (uice(i,j,2) - uice(i-1,j,2))*dxtra(i,j)
                  if (D .lt. 0.0 .or. presur(i,j) .gt. 0.0) then
c                   calculate pressure increment (dp)
                    dp = -D/((uvmx(i-1,j)/Ax(i-1,j)
     &                 + uvmx(i,j)/Ax(i,j))*cstra(i,j)*dxtra(i,j)**2 
     &                 + (csua(i,j-1)*uvmy(i,j-1)/Ay(i,j-1)
     &                 + csua(i,j)*uvmy(i,j)/Ay(i,j))*dytra(i,j)**2)
c                   limit dp by strength for convergence or pressure 
c                   for divergence
                    dp = min(dp, pmax(i,j) - presur(i,j))
                    dp = max(dp, -presur(i,j))
                    presur(i,j) = presur(i,j) + dp
c                   calculate and apply velocity corrections
                   du = uvmx(i-1,j)*dp*dxtra(i-1,j)*cstra(i,j)/Ax(i-1,j)
                    XX(i-1,j) = XX(i-1,j) - du
                    du = uvmx(i,j)*dp*dxtra(i,j)*cstra(i,j)/Ax(i,j)
                    XX(i,j) = XX(i,j) + du
                    dv = uvmy(i,j-1)*dp*dytra(i,j-1)/Ay(i,j-1)
                    YY(i,j-1) = YY(i,j-1) - dv
                    dv = uvmy(i,j)*dp*dytra(i,j)/Ay(i,j)
                    YY(i,j) = YY(i,j) + dv
                    emax = max(abs(du), emax)
                    emax = max(abs(dv), emax)
                  endif
                endif
              enddo
            enddo
          enddo
          do n=1,nseg
            do j=jsi(n),jei(n)
              do i=2,imtm1
                uice(i,j,2) = XX(i,j)
                vice(i,j,2) = YY(i,j)
              enddo
            enddo
          enddo
          call setbcx (presur, imt, jmt)
          call setbcx (uice(1,1,2), imt, jmt)
          call setbcx (vice(1,1,2), imt, jmt)
c         compare error with tolerance
          if (emax .le. tolcf) goto 20
        enddo
        write(*,*) '=> Warning: iterations exceeded in calculating',
     &    ' cf correction. iterations, maximum error:', kcf, emax
c
20    enddo
c
c-----------------------------------------------------------------------
c     advect ice thickness and area on C-grid with velocities at tau
c-----------------------------------------------------------------------
c
      call advice (uice(1,1,2), vice(1,1,2), hice, jsi, jei, nseg)
      call advice (uice(1,1,2), vice(1,1,2), aice, jsi, jei, nseg)
c
      return
      end


      subroutine advice (u, v, t, jsi, jei, nseg)
c=======================================================================
c     upstream advection of tracers with C-grid velocities
c     forward step kadv times using dt = dts/kadv
c
c     input:
c       u    = C-grid u component of ice velocity
c       v    = C-grid v component of ice velocity
c       t    = tracer to be advected (level lf)
c       jsi  = starting indices for ice segments
c       jei  = ending indices for ice segments
c       nseg = number of latitude segments for ice advection
c     output:
c       t = advected tracer (level 3)
c
c     author:   m.eby   e-mail: eby@uvic.ca
c=======================================================================
c
#  include "param.h"
#  include "cembm.h"
#  include "grdvara.h"
c
      real u(imt,jmt), v(imt,jmt), t(imt,jmt,3), afs(imt), tmp(imt,jmt)
c
      integer jsi(jmt), jei(jmt)
c
      do n=1,nseg
        do j=jsi(n),jei(n)
          do i=1,imt
            tmp(i,j) = t(i,j,lf)
          enddo
        enddo
      enddo
c
      dt = dts/float(kadv)
      do k=1, kadv
c
        do i=1,imt
          afs(i) = 0.0
        enddo
        do n=1,nseg
          do j=jsi(n),jei(n)
            afw = u(1,j)*(tmp(1,j) + tmp(2,j))
     &          + abs(u(1,j))*(tmp(1,j) - tmp(2,j))
            do i=2,imtm1
              afe = u(i,j)*(tmp(i,j) + tmp(i+1,j))
     &            + abs(u(i,j))*(tmp(i,j) - tmp(i+1,j))
              afn = v(i,j)*(tmp(i,j) + tmp(i,j+1))
     &            + abs(v(i,j))*(tmp(i,j) - tmp(i,j+1))
          t(i,j,3) = tmp(i,j) - dt*cstra(i,j)*((afe - afw)*dXt2ra(i,j)
     &             + (afn*csua(i,j) - afs(i)*csua(i,j-1))*dyt2ra(i,j))
              afw = afe
              afs(i) = afn
            enddo
          enddo
        enddo
c
        call setbcx (t(1,1,3), imt, jmt)
c
        do n=1,nseg
          do j=jsi(n),jei(n)
            do i=1,imt
              tmp(i,j) = t(i,j,3)
            enddo
          enddo
        enddo
c
      enddo
c
      return
      end
#endif
#if mrw_ice_adv
      subroutine advice (t)
c=======================================================================
c     upstream advection of tracers with C-grid velocities
c     forward step kadv times using dt = dts/kadv
c
c     input:
c       uice    = C-grid u component of ice velocity
c       vice    = C-grid v component of ice velocity
c       t    = tracer to be advected (level lf)
c       jsi  = starting indices for ice segments
c       jei  = ending indices for ice segments
c       nseg = number of latitude segments for ice advection
c     output:
c       t = advected tracer (level 3)
c
c     author:   m.eby   e-mail: eby@uvic.ca - adapted by g.bigg
c=======================================================================
c
#  include "param.h"
#  include "timelv.h"
#  include "cembm.h"
#  include "grdvara.h"
#  include "ice.h"
#  include "csbc.h"
#  include "levind.h"
c
      real t(imt,jmt,3), tmp(imt,jmt)
      real frac,fn,fs,fe,fw
c
c----------------------------------------------------------------------
c      horizontal advection
c----------------------------------------------------------------------
c
      isu=-1
      isv=0

      do j=2,jmtm1
      do i=2,imtm1
      im1=i-1
      jp1=j+1
      ip1=i+1
      jm1=j-1

      fn=(dxua(i,j)+dxua(ip1,j))*
     & (sbcocn(i,j,isv)+sbcocn(im1,j,isv))*p25
      fs=(dxua(i,jm1)+dxua(ip1,jm1))*
     & (sbcocn(i,jm1,isv)+sbcocn(im1,jm1,isv))*p25
      fe=(dyua(i,j)+dyua(i,jm1))*
     & (sbcocn(i,j,isu)+sbcocn(i,jm1,isu))*p25
      fw=(dyua(im1,j)+dyua(im1,jm1))*
     & (sbcocn(im1,j,isu)+sbcocn(im1,jm1,isu))*p25

      frac=0.01
      uice(i,j)=sbcocn(i,j,isu)*frac
      vice(i,j)=sbcocn(i,j,isv)*frac

            tmp(i,j)= frac/(dxua(i,j)*dyua(i,j))*(
c       u advection (zonal)
     &      fw * (t(i-1,j,3) + t(i,j,3))*p5
     &    - fe * (t(i+1,j,3) + t(i,j,3))*p5
c  v-advection (meridional
     &    + fs * (t(i,j-1,3) + t(i,j,3))*p5
     &    - fn * (t(i,j+1,3) + t(i,j,3))*p5  )
c
      end do
      end do

      do j=1,jmt
      do i=1,imt
      if(kmt(i,j).eq.0) then
      t(i,j,3)=t(i,j,3)
      else
      t(i,j,3)=t(i,j,3)+dts*tmp(i,j)
      endif
      end do
      end do

      return
      end
#endif
#endif
