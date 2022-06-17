c====================== include file "iounit.h" ========================
c
c     i/o units and related variables
c
c     taum1disk = disk pointer for tau-1 latitude rows
c     taudisk   = disk pointer for tau   latitude rows
c     taup1disk = disk pointer for tau+1 latitude rows
c     kflds     = disk unit used for two dimensional fields
c     latdisk   = disk units for latitude rows (alternately pointed to 
c                by taum1disk, taudisk, and taup1disk)
c     wide_open_mw = logical to indicate that the MW is fully opened.
c              if .true. then jmw = jmta and there are no latitude rows
c              on disk. Instead, they are all in the MW.
c              if .false. then jmw < jmta and all latitude rows are on
c              disk so they must be transferred between the MW and disk.
c              
c     iodoc  = unit for documentation
c     iostab = unit for stability testing
c     iotim  = unit for time means
c     iotim1 = scratch disk (SSD) unit for accumulating time means
c     ionew1 = unit for reading sponge layer data
c     ionew2 = mirror unit of sponge layer data on SSD
c     for the following, a control # < 0 implies that unformatted data
c     will be written to a unit selected by the i/o manager "iomngr.F"
c     and given a hardwired name (grep getunit *.F to see names) 
c     and formatted data (to stdout) will be written. if a # > 0 and
c      # <> stdout, only unformatted data will be written.
c
c     iotavg = control # for tracer averages
c     iotmb  = control # for writing tracer meridional budget.
c     iotrmb = control # for term balances for tracer and momentum
c     ioglen = control # for writing global energetics integrals
c     iovmsf = control # for writing meridional stream function
c     iogyre = control # for writing gyre transport.
c     ioprxz = control # for writing x-z sections from latitudes
c     ioext  = control # for writing external mode (stream function)
c     iodsp  = control # for writing diagnostic surface pressure
c     iotsi  = control # for writing time step integrals
c     ioxbt  = control # for writing time averaged xbt data
c     iozmbc = control # for writing zonal mean surf boundary conditons
c
      integer iotavg
      integer ioglen, iovmsf, iogyre, ioprxz, ioext, iodsp
      integer iotsi, iozmbc, ionew1, ionew2, ioxbt

      common /iounit/ iotavg
      common /iounit/ ioglen, iovmsf, iogyre, ioprxz, ioext, iodsp
      common /iounit/ iotsi, iozmbc, ionew1, ionew2, ioxbt
c
c     iotext  = 80 character text string for describing the details
c                of the next unformatted data record.
c     expnam  = 60 character text string for the experiment name
c
      character iotext*80, expnam*60
      common /iounitc/ iotext, expnam
