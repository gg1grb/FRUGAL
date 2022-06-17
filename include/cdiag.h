c===================== include file "cdiag.h" ========================
#ifdef hcomments
c
c     variables used for computing diagnostics:
c
c     ektot = total kinetic energy per unit volume
c     tddt  = rate of change of tracer per unit volume
c     dtabs = avg. per unit vol. of modulus of rate of change of tracer 
c     tvar  = rate of change of tracer variance per unit volume
c---------------------------------------------------------------------
#endif

      real    tddt(nt), dtabs(nt), tvar(nt), ektot
      common  /cdiag/  tddt, dtabs, tvar, ektot

