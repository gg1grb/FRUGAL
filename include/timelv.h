c===================== include file "timelv.h" =======================
#ifdef hcomments
c
c  time level pointers (indices) 
c
c  tracer pointers
c     tm = time level pointer (index) corresponding to "itt-2"
c     tc = time level pointer (index) corresponding to "itt-1"
c     tp = time level pointer (index) corresponding to "itt"
c       where "itt" is the baroclinic model timestep.
c
c  baroclinic pointers
c     um = time level pointer (index) corresponding to "itt-2"
c     uc = time level pointer (index) corresponding to "itt-1"
c     up = time level pointer (index) corresponding to "itt"
c       where "itt" is the baroclinic model timestep.
c
c  free surface model 
c     nm0 = time level pointer to "itt0-2"
c     nc0 = time level pointer to "itt0-1"
c     np0 = time level pointer to "itt0"
c       where "itt0" is the free surface model timestep.
c
c     ?nm  = next nm value (similarly for ?nc, ?np, nnm0, nnc0, nnp0)
c---------------------------------------------------------------------
#endif
#ifdef varts
      integer
     & nmu(imt,jmt),ncu(imt,jmt),npu(imt,jmt),
     & nmt(imt,jmt),nct(imt,jmt),npt(imt,jmt),
     & nmh(imt,jmt),nch(imt,jmt),nph(imt,jmt)
      integer nnmu(imt,jmt),nncu(imt,jmt),nnpu(imt,jmt),
     & nnmt(imt,jmt),nnct(imt,jmt),nnpt(imt,jmt),
     & nnmh(imt,jmt),nnch(imt,jmt),nnph(imt,jmt)
      integer ntst(imt,jmt),ntsu(imt,jmt),ntsh(imt,jmt),
     & difftstew(imt,jmt),difftstns(imt,jmt),
     & difftsuew(imt,jmt),difftsuns(imt,jmt),
     & difadind(subtsmax,subtsmax)
      integer 
     & uwa(-pps:pps,-pps:pps),uea(-pps:pps,-pps:pps),
     & usa(-pps:pps,-pps:pps),una(-pps:pps,-pps:pps),
     & uwd(-pps:pps,-pps:pps),ued(-pps:pps,-pps:pps),
     & usd(-pps:pps,-pps:pps),und(-pps:pps,-pps:pps)
      integer twa(-pps:pps,-pps:pps),tea(-pps:pps,-pps:pps),
     & tsa(-pps:pps,-pps:pps),tna(-pps:pps,-pps:pps),
     & twd(-pps:pps,-pps:pps),ted(-pps:pps,-pps:pps),
     & tsd(-pps:pps,-pps:pps),tnd(-pps:pps,-pps:pps)
      integer maxband,maxsubts,subts,
     & um,uc,up,unm,unc,unp,tm,tc,tp,tnm,tnc,tnp,
     & bm,bc,bp,bnm,bnc,bnp
      integer  nm0, nc0, np0
      integer  nnm0, nnc0, nnp0
      common  /ctlevp/ nmu,ncu,npu,nmt,nct,npt,nmh,nch,nph,
     & nnmu,nncu,nnpu,nnmt,nnct,nnpt,nnmh,nnch,nnph,
     & difftstew,difftstns,difftsuew,difftsuns,
     & difadind,ntst,ntsu,ntsh
      common /ctlevpa/
     & uwa,uea,usa,una,uwd,ued,usd,und,
     & twa,tea,tsa,tna,twd,ted,tsd,tnd
      common /ctlevpb/
     & maxband,maxsubts,subts
      common  /ctlevpc/
     & um,uc,up,unm,unc,unp,tm,tc,tp,tnm,tnc,tnp,
     & bm,bc,bp,bnm,bnc,bnp,
     &          nm0, nc0, np0,
     &         nnm0, nnc0, nnp0
#endif
#ifndef varts
      integer  tm, tc ,tp, um, uc, up, nm0, nc0, np0
      integer  tnm, tnc, tnp, unm, unc, unp, nnm0, nnc0, nnp0
      common  /ctlevp/
     &         tm, tc, tp, um, uc, up, nm0, nc0, np0, 
     &         tnm, tnc, tnp, unm, unc, unp, nnm0, nnc0, nnp0
#endif

