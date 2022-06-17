c======================== include file "isopyc.h" ======================
#ifdef isopycmixgriffies
c
c     isopycnal diffusion variables:
c
c     ahisop = isopycnal tracer mixing coefficient (cm**2/sec) 
c     drodx  = d(rho)/dx local to east face of T cell 
c     drody  = d(rho)/dy local to north face of T cell
c     drodz  = d(rho)/dz local to bottom face of T cell
c     Ai_e   = diffusion coefficient on eastern face of T cell
c     Ai_n   = diffusion coefficient on northern face of T cell
c     Ai_bx  = diffusion coefficient on bottom face of T cell
c     Ai_by  = diffusion coefficient on bottom face of T cell
c
c     fzisop = vertical structure function for the isopycnal diffusion
c              coefficient.
c     slmxr  = reciprocal of maximum allowable slope of isopycnals for
c              small angle approximation
c
      real K11, K22, K33
      real alphai,betai,ddxt,ddyt,ddzt,Ai_ez,Ai_nz,Ai_bx,Ai_by,
     &     ahisop,fzisop,slmxr

      common /cisop/ alphai(km,imt,jmt), betai(km,imt,jmt)

      common /cisop/ ddxt(km,imt,jmt,2)
      common /cisop/ ddyt(km,imt,jmt,2)
      common /cisop/ ddzt(0:km,imt,jmt,2)

      common /cisop/ Ai_ez(km,imt,jmt,0:1,0:1)
      common /cisop/ Ai_nz(km,imt,jmt,0:1,0:1)
      common /cisop/ Ai_bx(km,imt,jmt,0:1,0:1)
      common /cisop/ Ai_by(km,imt,jmt,0:1,0:1)
      common /cisop/ K11(km,imt,jmt)
      common /cisop/ K22(km,imt,jmt)
      common /cisop/ K33(km,imt,jmt)
      common /cisop/ fzisop(km), slmxr 
#ifdef england
      common /cisop/ ahisop(km)
#else
      common /cisop/ ahisop
#endif
      real tmask,tmaske,tmasks,tmaskw,tmaskn
      common /cisop/ tmask(km,0:imt,0:jmt),tmaskn(km,0:imt,0:jmt),
     & tmasks(km,0:imt,0:jmt),tmaskw(km,0:imt,0:jmt),
     & tmaske(km,0:imt,0:jmt)
# ifdef gent_mcwilliams
c     adv_vetiso = zonal isopycnal mixing velocity computed at the 
c                  center of the eastern face of the "t" cells
c     adv_vntiso = meridional isopycnal mixing velocity computed at
c                  the center of the northern face of the "t" cells
c                  (Note: this includes the cosine as in "adv_vnt")
c     adv_vbtiso = vertical isopycnal mixing velocity computed at the
c                  center of the top face of the "t" cells
c     adv_fbiso  = "adv_vbtiso" * (tracer) evaluated at the center of
c                  the bottom face of the "t" cells
c     athkdf = isopycnal thickness diffusivity (cm**2/sec)
      real adv_vetiso,adv_vntiso,adv_vbtiso,adv_fbiso
      common /cisop/ adv_vetiso(km,imt,jmt)
      common /cisop/ adv_vntiso(km,imt,jmt)
      common /cisop/ adv_vbtiso(0:km,imt,jmt)   
      common /cisop/ adv_fbiso(0:km)   
# endif
      integer ip,jq,kr,i_delta,j_delta,k_delta
      real slmx,athkdf,del_dm,s_dm,s_dmr,delta_iso,ft,
     & delta1,delta2,s_minus,s_plus,dtxsqr,tprime,sprime,sc,dzt4r,
     & Ai0,sumz,sxe,syn,fzisopb,sumx,sxb,sumy,facty,syb,flux_x,flux_y
      real diff_fe,diff_fw,diff_fn,diff_fs,diff_fbiso
      common /cisop1/ slmx,athkdf,del_dm,s_dm,s_dmr,delta_iso,ft,
     & delta1,delta2,s_minus,s_plus,dtxsqr(km),tprime,sprime,sc,dzt4r,
     & Ai0,sumz,sxe,syn,fzisopb,sumx,sxb,sumy,facty,syb,flux_x,flux_y
      common /cisop1/ip,jq,kr,i_delta,j_delta,k_delta
      common /cisop/ diff_fe(km),diff_fw(km),diff_fn(km),diff_fs(km)
     &              ,diff_fbiso(0:km)
      

      real drodxe
      real drodze
      real drodyn
      real drodzn
      real drodxb
      real drodyb
      real drodzb
#ifdef full_tensor
      real drodye
      real drodxn
#endif
      common /cisop2/ drodxe(km,imt,jmt,0:1)
     &               ,drodze(km,imt,jmt,0:1,0:1)
     &               ,drodyn(km,imt,jmt,0:1)
     &               ,drodzn(km,imt,jmt,0:1,0:1)
     &               ,drodxb(km,imt,jmt,0:1,0:1)
     &               ,drodyb(km,imt,jmt,0:1,0:1)
     &               ,drodzb(km,imt,jmt,0:1)
#ifdef full_tensor
     &               ,drodye(km,imt,jmt,0:1,0:1)
     &               ,drodxn(km,imt,jmt,0:1,0:1)
#endif
      integer iy,jy
      common /cisop2/ iy,jy
#endif
