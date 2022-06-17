c===================== include file "cvbc.h" =========================
#ifdef hcomments
c
c     global sea surface forcing variables:
c       sst = sea surface temperatures
c       sss = sea surface salinities
c       wsx = wind stress (in longitudinal (w->e) direction) 
c       wsy = wind stress (in latitudinal (s->n) direction)
c
c     vertical boundary condition variables:
c       smf = surface momentum flux   (dynes/cm**2)
c         1 => zonal wind stress
c         2 => meridional wind stress
c       bmf = bottom momentum flux
c         1 => zonal bottom drag
c         2 => meridional bottom drag 
c       stf = surface tracer flux    
c         1 => surface heat flux (cal/cm**2/sec = cm*degC/sec
c               = ly/sec) (assuming rho*cp = 1 cal/degC/cm**3)
c         2 => surface water flux (ppt/cm**2/sec)
c       btf = bottom tracer flux (for consistency but normally zero!)
c         1 => bottom heat flux
c         2 => bottom water flux
c
c     gamma  = relaxation timescale
c---------------------------------------------------------------------
#endif

      real    sst(imax,jmax), sss(imax,jmax)
#if defined ptisotope || defined ptisotopesalinity
      real    ssd18o(imax,jmax)
#endif
      real    wsx(imax,jmax), wsy(imax,jmax)
      real    smf(2), bmf(2), stf(nt), btf(nt)
      real    gamma
#ifdef obc
      real    ssh(imt,jmt),sshx(imt,jmt,12)
#endif
      common  /cvbc/
     &        sst, sss, wsx, wsy,
     &        smf, bmf, stf, btf
     &,       gamma
#if defined ptisotope || defined ptisotopesalinity
     &,       ssd18o
#endif
#ifdef obc
     &,       ssh,sshx
#endif
