c====================== include file "ident.h" =======================
#ifdef hcomments
c
c     definitions used to initialise and archive the dynamic arrays
c
c---------------------------------------------------------------------
#endif

      integer    n2dv, n3dv
      parameter (n2dv=3, n3dv=2+nt)

      real     fill2(n2dv), fill3(n3dv)
      character*40
     &         vid2(n2dv), vunit2(n2dv), vid3(n3dv), vunit3(n3dv)
      logical  vtgrd2(n2dv), vtgrd3(n3dv)

      common   /ident_r/  fill2, fill3
      common   /ident_c/  vid2, vid3, vunit2, vunit3
      common   /ident_l/  vtgrd2, vtgrd3

