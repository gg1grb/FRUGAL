c===================== include file "iounit.h" =======================
#ifdef hcomments
c
c     i/o units & related variables
c
c     fnrest  = main name of archive file used to start run
c     dir     = directory indicating restarts or snapshots
c     afile   = archive file name with appropriate extension
c
c     iorest  = unit number for general input/output 
c     unit_s  = unit number for slave sub domain restarts
c     ufmu    = unit number for unformatted archive
c     dacu    = unit number for direct access archive
c
c     stdin   = standard input unit number
c     stdout  = standard output unit number
c     stderr  = standard error unit number
c
c     dlev    = depth level at which a horizontal cross-section of
c               the 3D arrays is taken when writing out snapshots.
c---------------------------------------------------------------------
#endif

      integer    dlev
      character  fnrest*10, dir*10, afile*30, oplabel*1
      common /iofile/ 
     &           dlev, fnrest, dir, afile, oplabel

      integer    iorest, unit_s, ufmu, dacu
      parameter  (iorest=10, unit_s=11, ufmu=12, dacu=12)

      integer    stdin, stdout, stderr
      parameter  (stdin=5, stdout=6, stderr=6)

