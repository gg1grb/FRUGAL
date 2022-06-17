c===================== include file "state.h" ========================
#ifdef hcomments
c
c---------------------------------------------------------------------
#endif

      real    tq, sq
      integer kk

      real    to(km), so(km), c(km,9),
     &        to1(km), so1(km), co1(km,9),
     &        to2(km), so2(km), co2(km,9)
      common  /rstate/
     &        to, so, c, to1, so1, co1, to2, so2, co2
