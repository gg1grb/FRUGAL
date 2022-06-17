c===================== include file "coord.h" ========================
#ifdef hcomments
c
c          model grid point coordinates
c
c     grid definition:
c
c     the model uses a staggard arakawa "b" grid. horizontally, 
c     tracer quantities are defined at the centers of "t" grid 
c     boxes and velocities are defined at the centers of "u,v" grid 
c     boxes.  the centers of "u,v" grid boxes are located at the 
c     northeast corners of "t" grid boxes. the 1st "t" grid box is 
c     located in the southwest corner of the "t" grid. this grid 
c     system is replicated and stacked vertically one on top of 
c     another from the surface of the ocean downward. vertically,
c     tracers & velocities are defined at the centers of their 
c     corresponding boxes and are at the same depths.
c
c     primary grid information:
c
c     xt(i) = longitude of the ith "t" point in degrees.   i=1..imt
c     xu(i) = longitude of the ith "u,v" point in degrees. i=1..imt
c     yt(j) = latitude of the jth "t" point in degrees.    j=1..jmt
c     yu(j) = latitude of the jth "u,v" point in degrees.  j=1..jmt
c     zt(k) = distance (cm) from surface down to center of level k 
c             (depth of "t" and "u,v" grid points)
c     zw(k) = distance (cm) from surface down to bottom of level k 
c---------------------------------------------------------------------
#endif

      real xt(imt,0:jmt+1),yt(imt,0:jmt+1),xu(imt,0:jmt),yu(imt,0:jmt)
      real    zw(0:km), zt(0:km+1)
      common  /coord/
     &        xt, yt, xu, yu, zw, zt

