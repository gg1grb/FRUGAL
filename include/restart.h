c===================== include file "restart.h" ======================
#ifdef hcomments
c
c     global definitions and variables for archive reading/writing
c
c---------------------------------------------------------------------
#endif

      integer    i, j, k, l, ns, ew
      integer    filetype, srl, cri
      common    /restart_scalars/
     &           filetype, srl, cri

      real       full2d(iall,jall), full3d(km,iall,jall)
      real       sub2d(imax,jmax+km), sub3d(km,imax,jmax)
      integer    map(0:ewsg-1, 0:nssg-1)
      integer    imtsd(0:ewsg-1, 0:nssg-1)
      integer    jmtsd(0:ewsg-1, 0:nssg-1)
      integer    ioffsd(0:ewsg-1, 0:nssg-1)
      integer    joffsd(0:ewsg-1, 0:nssg-1)
      common    /restart_arrays/
     &           full2d, full3d, sub2d, sub3d,
     &           map, imtsd, jmtsd, ioffsd, joffsd

      character*80 ifile(0:ewsg-1, 0:nssg-1)
      character*80 tfile(0:ewsg-1, 0:nssg-1)
      common    /restart_files/
     &           ifile, tfile

#define LEVITUS         1
#define RESTART         2
#define SNAPSHOT        3

#ifdef hdf

#define OBSERVATIONAL   1
#define MODELLED        2

#define RESOLUTION_4    1
#define RESOLUTION_6    2

#define NON_ROTATED     1
#define ROTATED         2

#define TGRID           1

#define FMODE_ATTR      'FMODE'
#define FTYPE_ATTR      'FTYPE'
#define TGRID_ATTR      'T_GRID'
#define LEVEL_ATTR      'LEVEL'
#define TSTEP_ATTR      'TIMESTEP'
#define MYEAR_ATTR      'YEAR'
#define MDAY_ATTR       'DAY'
#define MSEC_ATTR       'SECONDS'
#define RESOL_ATTR      'RESOLUTION'
#define DATE_ATTR       'DATE'
#define VERSION_ATTR    'VERSION'
#define ROTATION_ATTR   'ROTATION'

#define   INT_HDF   DFNT_INT32
# define  REAL_HDF  DFNT_FLOAT32

      integer*4  sfstart, sfcreate, sfendacc, sfend,
     &           sfsdtstr, sfwdata, sfsfill, sfdimid,
     &           sfsdmname, sfsdmstr, sfsdscale, sfsattr,
     &           sfn2index, sfselect, sfrdata
      external   sfstart, sfcreate, sfendacc, sfend,
     &           sfsdtstr, sfwdata, sfsfill, sfdimid,
     &           sfsdmname, sfsdmstr, sfsdscale, sfsattr,
     &           sfn2index, sfselect, sfrdata

      integer*4  dimslon, dimslat, dimslev
      parameter  (dimslon=iall, dimslat=jall, dimslev=km)

      character*40
     &           lonlabel(2), latlabel(2),
     &           label2(2), unit2(2), fmat2(2),
     &           label3(3), unit3(3), fmat3(3)
      integer*4  sdsid, ret, hdfu, index, dimid1, dimid2, dimid3
      real       xt(iall), xu(iall), yt(jall), yu(jall), zt(km+1)

      common /hdf_c/
     &        lonlabel, latlabel,
     &        label2, unit2, fmat2, label3, unit3, fmat3
      common /hdf_r/
     &        xt, xu, yt, yu, zt
      common /hdf_i/
     &        sdsid, ret, hdfu, index, dimid1, dimid2, dimid3

#endif
