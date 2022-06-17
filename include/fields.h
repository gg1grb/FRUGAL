c===================== include file "fields.h" =======================
#ifdef hcomments
c
c  dynamic field arrays:
c     t = tracer quantites. the 1st is temperature & 
c                           the 2nd is salinity
c         if nt > 2 then other tracers are allowed.
c     u = zonal component of velocity.
c     v = meridional component of velocity.
c---------------------------------------------------------------------
#endif

#ifdef quick
      real    t(km,0:imax+1,0:jmax+1,nt,3)
#else
      real    t(km,imax,jmax,nt,3)
#endif
#if defined ptmonthly || defined ptcore
      real    u(km,imax,jmax), v(km,imax,jmax)
#else
      real    u(km,imax,jmax,3), v(km,imax,jmax,3)
#endif
#ifdef ptaccel
      real ptstore(km,imt,jmt),ptstorem(km,imt,jmt)
#endif
      common  /slabs/
     &        t, u, v 
#ifdef ptaccel
     &,       ptstore,ptstorem
#endif
#ifdef hcomments
c---------------------------------------------------------------------
c  the working arrays are
c 
c     uclin  = internal mode u (used for advective fluxes &
c               diagnostics)
c     vclin  = internal mode v 
c     rhoo   = density at (i  ,j  ) on the "t" grid
c     rhpo   = density at (i+1,j  ) on the "t" grid
c     rhop   = density at (i  ,j+1) on the "t" grid
c     rhpp   = density at (i+1,j+1) on the "t" grid
c     dpdx   = zonal gradient of pressure on "u,v" grid
c     dpdy   = meridional gradient of pressure on "u,v" grid
c     fue    = advective coeff for eastern face of "u,v" grid box
c              in "clinic" & "t" grid box in "tracer"
c     fuw    = advective coeff for western face of grid box
c     fvn    = advective coeff for northern face of grid box
c     fvs    = advective coeff for southern face of grid box
c     fm     = (0,1) over "t" grid (land,ocean) points
c     gm     = (0,1) over "u,v" grid (land,ocean) points
c     vmf    = array (1) used for vertical differences of u
c     vmf    = array (2) used for vertical differences of v
c     vtf    = array used for vertical tracer flux
c     fw     = vertical velocity defined at the bottom of "u,v" 
c              boxes in "clinic" and "t" boxes in "tracer"
c     fwb1   = w * (quantity 1) defined at the bottoms of the boxes
c     fwb2   = w * (quantity 2) defined at the bottoms of the boxes
c     restr  = term for restoring surface tracers to prescribed 
c              values via newtonain damping
c     rests  = time scale for restoring surface tracers (days)
c
c     maskoo = land mask at (i,j) on the 't' grid
c     maskpo = land mask at (i+1,j) on the 't' grid
c     maskmo = land mask at (i-1,j) on the 't' grid
c     maskop = land mask at (i,j+1) on the 't' grid
c     maskom = land mask at (i,j-1) on the 't' grid
c---------------------------------------------------------------------
#endif

      real    dpdx(km), dpdy(km), vmf(0:km,2), vtf(0:km)
      real    fue(km), fuw(km), fvn(km), fvs(km)
      real    fw(0:km), fwb1(0:km), fwb2(0:km)
      real    rhoo(km), rhpo(km), rhpp(km), rhop(km)
      real    totalsalt,fwadjustatlantic,fwadjustsouth
      real    tempz(imt,jmt)
      integer maskoo(km)
      integer maskpo(km), maskmo(km), maskop(km), maskom(km)
      integer mask(imt,jmt)
      real    atover,southover,btdpkeep
      common  /work/
     &        dpdx, dpdy, fue, fuw, fvn, fvs, vmf, vtf,
     &        fw, fwb1, fwb2, rhoo, rhpo, rhpp, rhop,
     &        maskoo, maskpo, maskmo, maskop, maskom,totalsalt,
     &        mask,atover,southover,fwadjustatlantic,fwadjustsouth
     & ,btdpkeep
      real volumeflux(10)
      common /tsdatax/ volumeflux
      common /tsdatax/ tempz
#if defined relaxsurface && defined fullfs
      real relaxfwflux(imax,jmax)
      common /relaxfw/ relaxfwflux
#endif
#ifdef brine
      real ice(imt,jmt)
      common /icebrine/ ice
#endif
#if defined control || defined transient || defined crad_th2o || defined trad_ch2o
      real fwstore(imt,jmt,0:360*141)
      common /cmipstore/ fwstore
#endif
#if defined hosep1 || defined hosec1
      real maskhose(imt,jmt),areahose,freshfluxhose
      common /cmiphose/ maskhose,areahose,freshfluxhose
#endif
#if defined control || defined transient || defined crad_th2o || defined trad_ch2o || defined hosep1 || defined hosec1
c cmipoutput(?,1) start time for averaging
c cmipoutput(?,2) end time for averaging
c cmipoutput(1,?) 1 day means
c cmipoutput(2,?) 5 day means
c cmipoutput(3,?) 30 day means
      integer cmipoutput(3,2)
      common /cmipout/ cmipoutput
#endif
#ifndef coupled
c heat and fw fluxes for ocean-only model
      real oco_hf(imt,jmt),oco_fw(imt,jmt)
      common /oconly/ oco_hf, oco_fw
#endif
