      subroutine solve (n)
      implicit none

c
c=======================================================================
c     solve for tracer distribution after diffusion using SLAP
c     based on code by a.fanning
c
c     input:
c       n    = tracer number
c
c     author:   m.eby   e-mail: eby@uvic.ca
c=======================================================================
c
# include "param.h"
# include "slap.h"
# include "atm.h"
# include "cembm.h"
#include "levinda.h"
#ifdef advectq
#include "coorda.h"
#include "grdvara.h"
#endif

c
      real fa,tol,err
      integer isym,itol,itmax,iunit,iter,ierr
#ifdef advectq
      real fvn,fvs,fue,fuw
      integer im,ip,jm,jp
#endif
c
c-----------------------------------------------------------------------
c     move previous solution in time
c-----------------------------------------------------------------------
c
      do j=1,jmta
        do i=1,imta
          at(i,j,1,n) = at(i,j,2,n)
          at(i,j,2,n) = 0.0
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     if solving for temperature add in latent heat from precipitation
c-----------------------------------------------------------------------
c
      if (n .eq. 1) then 
        fa = dts*vlocn*rhoocn/(cpatm*rhoatm*ht)
        do j=2,jmtam1
          do i=2,imtam1
             k = indexgrid(i,j)
               slaprhs(k,n) = slaprhs(k,n) + fa*precip(i,j)
              if(kmta(i,j).eq.1)
     &   slaprhs(k,n)=slaprhs(k,n)-dts/(cpatm*rhoatm*ht)*upltnt(i,j)
          enddo
        enddo
      endif


#ifdef advectq
c----------------------------------------------------------------------
c  calculate change in humidity due to advection by wind
c----------------------------------------------------------------------




      if (n .eq. 2) then

c ensure cyclic conditions are set

      do j=1,jmta
      spdu(1,j)=spdu(imta-1,j)
      spdu(imta,j)=spdu(2,j)
      spdv(1,j)=spdv(imta-1,j)
      spdv(imta,j)=spdv(2,j)
      at(1,j,2,2)=at(imta-1,j,2,2)
      at(imta,j,2,2)=at(2,j,2,2)
      dyua(1,j)=dyua(imta-1,j)
      dyua(imta,j)=dyua(2,j)
      dxua(1,j)=dxua(imta-1,j)
      dxua(imta,j)=dxua(2,j)
      end do


        do j=2,jmt-1
        do i=2,imt-1

        im=i-1
        ip=i+1
        jm=j-1
        jp=j+1

        fvn =
     & (spdv(i,j)+spdv(im,j)) *
     & (dxua(i,j)+dxua(im,j))*p25
        fvs =
     & (spdv(i,jm)+spdv(im,jm)) *
     & (dxua(i,jm)+dxua(im,jm))*p25
        fue =
     & (spdu(i,j)+spdu(i,jm)) *
     & (dyua(i,j)+dyua(i,jm))*p25
        fuw =
     & (spdu(im,j)+spdu(im,jm)) *
     & (dyua(im,j)+dyua(im,jm))*p25

c use time level 2 as temporary storage

         k = indexgrid(i,j)

         slaprhs(k,n)=slaprhs(k,n) +dts*(
c       at(i,j,2,2) = at(i,j,1,2) + dts* (
c       at(i,j,2,2) = at(i,j,1,2) + 0.5*dts* (
c  u-advection (zonal)
     &      fuw * (at(im,j,1,2) + at(i,j,1,2))*p5
     &    - fue * (at(ip,j,1,2) + at(i,j,1,2))*p5
c  v-advection (meridional)
     &    + fvs * (at(i,jm,1,2) + at(i,j,1,2))*p5
     &    - fvn * (at(i,jp,1,2) + at(i,j,1,2))*p5  )

     &   /(dxta(i,j)*dyta(i,j))

        enddo
        enddo

c put back into apropriate time level

c       do j=2,jmt-1
c         do i=2,imt-1
c           at(i,j,1,2)=at(i,j,2,2)
c        end do
c      end do
cifdef advectqpositive
c       do j=2,jmt-1
c         do i=2,imt-1
c           if(at(i,j,1,2).lt.0.0)at(i,j,1,2)=0.0
c        end do
c      end do
cendif

c  set cylic conditions

c     do j=1,jmta
c     at(1,j,1,2)=at(imta-1,j,1,2)
c     at(imta,j,1,2)=at(2,j,1,2)
c     end do

      endif

#endif


c
      tol = sltol(n)
c
      isym = 0      !zero for non-symmetric, one for symmetric
      itol = 0      !stopping criterion selector
      itmax = 900   !max number of iterations
      iunit = 0     !output unit number: zero for no monitoring
c
c-----------------------------------------------------------------------
c     solve for tracer using "slap"
c-----------------------------------------------------------------------
c
       call sslugm(orda, slaprhs(1,n), slaplhs(1,n), nelt
     &,           slapia(1,lf,n), slapja(1,lf,n), slapam(1,lf,n)
     &,           isym, nsave, itol, tol, itmax, iter, err
     &,           ierr, iunit, rwork, lenw, iwork, leniw)
 
      if (ierr .ne. 0) write(*,*) 'iter, ierr, n =', iter, ierr, n
c
c-----------------------------------------------------------------------
c     copy new solution from left hand side
c-----------------------------------------------------------------------
c
      do k=1,orda
        i = gridcol(k)
        j = gridrow(k)
        at(i,j,2,n) = slaplhs(k,n)
      enddo
c
c-----------------------------------------------------------------------
c     set boundary conditions
c-----------------------------------------------------------------------
c
      call embmbc (at(1,1,2,n))
c
      return 
      end


      subroutine indexmaps
      implicit none
c
c=======================================================================
c     produces maps from interior tracer points on grid to vector of
c     unknowns and back again. cyclic and closed N/S conditions handled
c     implicitly so only use 2:imtm1, 2:jmtm1     
c     based on code by a.fanning
c
c     author:   m.eby   e-mail: eby@uvic.ca
c=======================================================================
c
# include "param.h"
# include "slap.h"
c
c
      orda = 0
      do j=2,jmtam1
        do i=2,imtam1
          orda = orda + 1
          indexgrid(i,j) = orda
          gridrow(orda) = j
          gridcol(orda) = i
        enddo
      enddo
c
      return
      end


      subroutine newrhs
      implicit none
c
c=======================================================================
c     calculate new right hand side
c     based on code by a.fanning
c
c     author:   m.eby   e-mail: eby@uvic.ca
c=======================================================================
c
# include "param.h"
# include "slap.h"
# include "cembm.h"
# include "atm.h"
# include "levinda.h"

      real fa,fc
      integer kk
c
c-----------------------------------------------------------------------
c     air temperature
c-----------------------------------------------------------------------
c
      fa = dts/(cpatm*rhoatm*ht)
      do j=2,jmtam1
c       incoming swr absorbed by the atmosphere is fc*dnswr
        do i=2,imtam1
        fc = scatter(i,j)/(1.0 - scatter(i,j))
          k = indexgrid(i,j)
          if (kmta(i,j) .eq. 0) then
            slaprhs(k,1) = at(i,j,lf,1) + fa*(fc*dnswr(i,j)
     &                   - outlwr(i,j) + uplwr(i,j) + upsens(i,j))
          else
            slaprhs(k,1) = at(i,j,lf,1) + fa*(fc*dnswr(i,j) 
     &                   - outlwr(i,j) + dnswr(i,j))
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     specific humidity
c-----------------------------------------------------------------------
c
      fa = dts*rhoocn/(rhoatm*hq)
      do j=2,jmtam1
        do i=2,imtam1
          k = indexgrid(i,j)
          if (kmta(i,j) .eq. 0) then
            slaprhs(k,2) = at(i,j,lf,2) + fa*evap(i,j)
          else
            slaprhs(k,2) = at(i,j,lf,2)
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     add other tracer source and sink terms here
c-----------------------------------------------------------------------
c
      do n=3,nat
        do j=2,jmtam1
          do i=2,imtam1
            k = indexgrid(i,j)
            slaprhs(k,n) = at(i,j,lf,n)
          enddo
        enddo
      enddo


c
      return
      end


      subroutine coef (n)
      implicit none
c
c=======================================================================
c     compute matrix coefficients
c       coefficient storage
c       non-zero element A(i,j) is stored in AM(pos)
c       row index i             is stored in ia(pos)
c       column index j          is stored in ja(pos)
c       coefficients for a single grid point are stored as a group
c	these form one row of the coefficient matrix
c     based on code by a.fanning
c
c     input:
c       n    = tracer number
c
c     author:   m.eby   e-mail: eby@uvic.ca
c=======================================================================
c
# include "param.h"
# include "slap.h"
# include "grdvara.h"
# include "cembm.h"
# include "atm.h"
c
      real fa2,fc2,fb3,fb1,fb2
      integer index,itmp
c
      index = 0
c
      do j=2,jmtam1
c
        do i=2,imtam1
c
        fa2 = dts*(da(i,j,n) + da(i,j-1,n))*dyura(i,j-1)*
     &        dyt2ra(i,j)
        fc2 = dts*(da(i,j,n) + da(i,j+1,n))*dyura(i,j)*
     &        dyt2ra(i,j)
          fb3 = dts*da(i,j,n)*dxura(i,j)*dxtra(i,j)
          fb1 = dts*da(i,j,n)*dxura(i-1,j)*dxtra(i,j)
          fb2 = 1.0 - fa2 - fc2 - fb1 - fb3
c
c         closed N/S bdys:
c
          if (j .eq. 2) fb2 = fb2 + fa2
          if (j .eq. jmtam1) fb2 = fb2 + fc2
c
c         Store the relevant coefficients
c         and their row/col indices in am, ia, ja respect.
c
c         main block diagonal
          itmp = indexgrid(i,j)        ! row and col index of b2
          index = index + 1
          slapam(index,lf,n) = fb2     ! diagonal of block and matrix
          slapia(index,lf,n) = itmp
          slapja(index,lf,n) = itmp
c
          if (i .gt. 2) then           ! first sub diagonal
            index = index + 1
            slapam(index,lf,n) = fb1
            slapia(index,lf,n) = itmp
            slapja(index,lf,n) = itmp-1
          else 
            if (i .eq. 2) then         ! cyclic domain i-1 maps to nx 
              index = index + 1
              slapam(index,lf,n) = fb1
              slapia(index,lf,n) = itmp
              slapja(index,lf,n) = indexgrid(imtam1,j)
            endif
          endif
c
          if (i .lt. imtam1) then       ! first super diagonal
            index = index + 1
            slapam(index,lf,n) = fb3
            slapia(index,lf,n) = itmp
            slapja(index,lf,n) = itmp+1
          else 
            if (i .eq. imtam1) then     ! cyclic domain nx+1 maps to 1 
              index = index + 1
              slapam(index,lf,n) = fb3
              slapia(index,lf,n) = itmp
              slapja(index,lf,n) = indexgrid(2,j)
            endif
          endif
c
c         block sub-diagonal
          if (j .gt. 2) then           ! element in sub-diagonal block
            index = index + 1
            slapam(index,lf,n) = fa2   ! diagonal element (of block)
            slapia(index,lf,n) = itmp
            slapja(index,lf,n) = indexgrid(i,j-1)
          endif
c
c         block super-diagonal
          if (j .lt. jmtam1) then       ! element in sub-diagonal block
            index = index + 1
            slapam(index,lf,n) = fc2   ! diagonal element (of block)
            slapia(index,lf,n) = itmp
            slapja(index,lf,n) = indexgrid(i,j+1)
          endif
c
        enddo
      enddo
c
      nelt = index

c
      return
      end
