#include <define.h>
! ======================================================
! Aggreate/screen high-resolution urban dataset
! to a lower resolutioin/subset data, suitable for running
! regional or point cases.
! ======================================================
SUBROUTINE makeurbandata( casename,dir_rawdata,dir_srfdata, &
                          lc_year,edgen,edgee,edges,edgew )

   USE precision
   USE netcdf
   USE ncio
   USE omp_lib
 
   IMPLICIT NONE

   CHARACTER(LEN=256), intent(in) :: casename
   CHARACTER(len=256), intent(in) :: dir_rawdata
   CHARACTER(len=256), intent(in) :: dir_srfdata

   INTEGER , intent(in) :: lc_year
   REAL(r8), intent(in) :: edgen
   REAL(r8), intent(in) :: edgee
   REAL(r8), intent(in) :: edges
   REAL(r8), intent(in) :: edgew

   CHARACTER(len=*), parameter :: DATASRC = "MOD"
   CHARACTER(len=*), parameter :: Title   = "Land surface model input urban data"
   CHARACTER(len=*), parameter :: Authors = "Dai YJ group at Sun Yat-sen University"
   CHARACTER(len=*), parameter :: Address = "School of Atmospheric Sciences, Sun Yat-sen University, Zhuhai, China"
   CHARACTER(len=*), parameter :: Email   = "yuanh25@mail.sysu.edu.cn"

#ifdef USE_LCZ
   INTEGER, parameter :: hxy  = 6000
   INTEGER, parameter :: nlcz = 10
#else
   INTEGER, parameter :: den_clss = 3
#endif
   INTEGER, parameter :: nxy  = 1200
   INTEGER, parameter :: bxy  = 600
   INTEGER, parameter :: mon  = 12
   INTEGER, parameter :: rid  = 33
   INTEGER, parameter :: ns   = 2
   INTEGER, parameter :: nr   = 2
   INTEGER, parameter :: ulev = 10

! input variables
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: hlat
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: hlats
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: hlatn
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: hlon
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: hlonw
   REAL(r8), ALLOCATABLE, DIMENSION(:)   :: hlone
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: gfcc_tc
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: gedi_th
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: modur
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: gl30_wt
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: harea
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: wtrf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: htrf

#ifdef USE_LCZ
   INTEGER , ALLOCATABLE, DIMENSION(:,:) :: lcz
#else
   INTEGER , ALLOCATABLE, DIMENSION(:,:) :: urrgid
   INTEGER , ALLOCATABLE, DIMENSION(:,:) :: urden
#endif

   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: hlai, hsai

#ifndef USE_LCZ
   REAL(r8), DIMENSION(den_clss,rid)  :: hwrcan, wtrd, emrf, emwl, ncar_wt
   REAL(r8), DIMENSION(den_clss,rid)  :: emimrd, emperd, ncar_ht
   REAL(r8), DIMENSION(den_clss,rid)  :: thrf, thwl, tbmin, tbmax
   
   REAL(r8), DIMENSION(den_clss,rid,ulev ):: cvrf, cvwl, cvimrd, &
                                             tkrf, tkwl, tkimrd
   REAL(r8), DIMENSION(den_clss,rid,nr,ns):: albrf, albwl, albimrd, albperd
#endif
   
   ! output variables
   REAL(r8), ALLOCATABLE, DIMENSION(:)     :: latso
   REAL(r8), ALLOCATABLE, DIMENSION(:)     :: lonso
   REAL(r8), ALLOCATABLE, DIMENSION(:,:)   :: area
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: wgt_top
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: tc
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: urwt
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: htop
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: pct_tc
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: pct_urwt
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: htop_ur
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: ht_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: wt_rf
#ifndef USE_LCZ
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: hwr_can
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: wt_rd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: em_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: em_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: em_imrd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: em_perd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: th_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: th_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: tb_min
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: tb_max

   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: cv_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: cv_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: cv_imrd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: tk_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: tk_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: tk_imrd

   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: alb_rf
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: alb_wl
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: alb_imrd
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: alb_perd
#endif

   REAL(r8), ALLOCATABLE, DIMENSION(:,:)     :: hgt, avg
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:)   :: ur_dc
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:)   :: pct_ur
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: ur_lai
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: ur_sai
   
   REAL(r8), ALLOCATABLE, DIMENSION(:,:,:,:) :: wgt_lai, wgt_sai   
!-----------------------------------------------------------------

   CHARACTER(len=256) lndname
   CHARACTER(len=256) cyear
   CHARACTER(len=256) :: reg1, reg2, reg3, reg4

   INTEGER :: reg(4)
   ! variable ids
   INTEGER :: lat_dimid, lon_dimid, mon_dimid, den_dimid
   INTEGER :: ncid, lat_vid, lon_vid, urlat_vid, urlon_vid, mon_vid, den_vid, &
              pct_tcvid, pct_urvid, pct_urwtvid, ur_laivid, ur_saivid, htop_urvid
   INTEGER :: upftvid, gedi_thvid, gfcc_tcvid, gl30_wtvid, ur_clssvid, &
              laivid, saivid, wt_rfvid, ht_rfvid 
#ifndef USE_LCZ
   INTEGER :: ns_dimid, nr_dimid, ulev_dimid
   INTEGER :: ns_vid, nr_vid, lev_vid
   INTEGER :: ur_rgvid, hwr_canvid, wt_rdvid, em_rfvid, em_wlvid, em_imrdvid, em_perdvid
   INTEGER :: cv_rfvid, cv_wlvid, cv_imrdvid
   INTEGER :: th_rfvid, th_wlvid, tbminvid, tbmaxvid
   INTEGER :: tk_rfvid, tk_wlvid, tk_imrdvid
   INTEGER :: alb_rfvid, alb_imrdvid, alb_perdvid, alb_wlvid
   INTEGER :: uxid
#endif

   REAL(r8) :: ldelta
   REAL(r8) :: pi, deg2rad, re, dx, dy, sumur
   REAL(r8) :: x_delta, y_delta
   
   INTEGER  :: i, j, k, io, jo, m, n, ii, jj, mm, cont, sumth, p, inx
   INTEGER  :: n_ns(2), n_nr(2), n_den(3), n_ulev(10), n_mon(12)
   INTEGER  :: XY2D(2), XY3D(3), XY4D(4), UR3D(3), UL3D(4), XY5D(5)

   INTEGER  :: ei, ej, si, sj
   INTEGER  :: reglat,reglon,reglon_,sreglat,sreglon,ereglat,ereglon
   LOGICAL  :: fileExists

   pi = 4._r8*atan(1.)
   deg2rad = pi/180._r8
   re = 6.37122e6 * 0.001

   write(cyear,'(i4.4)') lc_year

   print*, ">>> allocating memory..."
   ! IF (USE_LCZ) THEN
#ifdef USE_LCZ
   allocate( hlat  (hxy) )
   allocate( hlats (hxy) )
   allocate( hlatn (hxy) )
   allocate( hlon  (hxy) )
   allocate( hlonw (hxy) )
   allocate( hlone (hxy) )

   allocate( harea(hxy, hxy) )
   allocate( lcz  (hxy, hxy) )

   allocate( hlai (hxy, hxy, mon) )
   allocate( hsai (hxy, hxy, mon) )

   allocate( ur_dc     (lon_points, lat_points, nlcz ) )
   allocate( pct_ur    (lon_points, lat_points, nlcz ) )
   allocate( wgt_top   (lon_points, lat_points, nlcz ) )
   allocate( tc        (lon_points, lat_points, nlcz ) )
   allocate( urwt      (lon_points, lat_points, nlcz ) )
   allocate( htop      (lon_points, lat_points, nlcz ) )
   allocate( pct_tc    (lon_points, lat_points, nlcz ) )
   allocate( pct_urwt  (lon_points, lat_points, nlcz ) )
   allocate( htop_ur   (lon_points, lat_points, nlcz ) )

   allocate( ur_sai    (lon_points, lat_points, nlcz, mon ) )
   allocate( ur_lai    (lon_points, lat_points, nlcz, mon ) )
   allocate( wgt_sai   (lon_points, lat_points, nlcz, mon ) )
   allocate( wgt_lai   (lon_points, lat_points, nlcz, mon ) )
   ! ELSE
#else
   allocate( hlat  (nxy) )
   allocate( hlats (nxy) )
   allocate( hlatn (nxy) )
   allocate( hlon  (nxy) )
   allocate( hlonw (nxy) )
   allocate( hlone (nxy) )

   allocate( harea   (nxy, nxy) )
   allocate( urrgid  (nxy, nxy) )
   allocate( urden   (nxy, nxy) )

   allocate( hlai    (nxy, nxy, mon) )
   allocate( hsai    (nxy, nxy, mon) )

   
   allocate( ur_dc     (lon_points, lat_points, den_clss ) )
   allocate( pct_ur    (lon_points, lat_points, den_clss ) )
   allocate( wgt_top   (lon_points, lat_points, den_clss ) )
   allocate( tc        (lon_points, lat_points, den_clss ) )
   allocate( urwt      (lon_points, lat_points, den_clss ) )
   allocate( htop      (lon_points, lat_points, den_clss ) )
   allocate( pct_tc    (lon_points, lat_points, den_clss ) )
   allocate( pct_urwt  (lon_points, lat_points, den_clss ) )
   allocate( htop_ur   (lon_points, lat_points, den_clss ) )
   allocate( hwr_can   (lon_points, lat_points, den_clss ) )
   allocate( wt_rf     (lon_points, lat_points, den_clss ) )
   allocate( wt_rd     (lon_points, lat_points, den_clss ) )
   allocate( em_rf     (lon_points, lat_points, den_clss ) )
   allocate( em_wl     (lon_points, lat_points, den_clss ) )
   allocate( em_imrd   (lon_points, lat_points, den_clss ) )
   allocate( em_perd   (lon_points, lat_points, den_clss ) )
   allocate( ht_rf     (lon_points, lat_points, den_clss ) )
   allocate( th_rf     (lon_points, lat_points, den_clss ) )
   allocate( th_wl     (lon_points, lat_points, den_clss ) )
   allocate( tb_min    (lon_points, lat_points, den_clss ) )
   allocate( tb_max    (lon_points, lat_points, den_clss ) )
   allocate( ur_sai    (lon_points, lat_points, den_clss, mon ) )
   allocate( ur_lai    (lon_points, lat_points, den_clss, mon ) )
   allocate( wgt_sai   (lon_points, lat_points, den_clss, mon ) )
   allocate( wgt_lai   (lon_points, lat_points, den_clss, mon ) )
   allocate( cv_rf     (lon_points, lat_points, den_clss, ulev) )
   allocate( cv_wl     (lon_points, lat_points, den_clss, ulev) )
   allocate( cv_imrd   (lon_points, lat_points, den_clss, ulev) )
   allocate( tk_rf     (lon_points, lat_points, den_clss, ulev) )
   allocate( tk_wl     (lon_points, lat_points, den_clss, ulev) )
   allocate( tk_imrd   (lon_points, lat_points, den_clss, ulev) )
   allocate( alb_rf    (lon_points, lat_points, den_clss, nr, ns) )
   allocate( alb_wl    (lon_points, lat_points, den_clss, nr, ns) )
   allocate( alb_imrd  (lon_points, lat_points, den_clss, nr, ns) )
   allocate( alb_perd  (lon_points, lat_points, den_clss, nr, ns) )

   hwrcan  (:,:) = 0.
   !wtrf    (:,:) = 0.
   wtrd    (:,:) = 0.
   emrf    (:,:) = 0.
   emwl    (:,:) = 0.
   emimrd  (:,:) = 0.
   emperd  (:,:) = 0.
   !htrf    (:,:) = 0.
   thrf    (:,:) = 0.
   thwl    (:,:) = 0.
   tbmin   (:,:) = 0.
   tbmax   (:,:) = 0.

   hwr_can  (:,:,:) = 0.
   wt_rf    (:,:,:) = 0.
   wt_rd    (:,:,:) = 0.
   em_rf    (:,:,:) = 0.
   em_wl    (:,:,:) = 0.
   em_imrd  (:,:,:) = 0.
   em_perd  (:,:,:) = 0.
   ht_rf    (:,:,:) = 0.
   th_rf    (:,:,:) = 0.
   th_wl    (:,:,:) = 0.
   tb_min   (:,:,:) = 0.
   tb_max   (:,:,:) = 0.
   tkrf     (:,:,:) = 0.
   tkwl     (:,:,:) = 0.
   tkimrd   (:,:,:) = 0.
   cvrf     (:,:,:) = 0.
   cvwl     (:,:,:) = 0.
   cvimrd   (:,:,:) = 0.

   tk_rf    (:,:,:,:) = 0.
   tk_wl    (:,:,:,:) = 0.
   tk_imrd  (:,:,:,:) = 0.
   cv_rf    (:,:,:,:) = 0.
   cv_wl    (:,:,:,:) = 0.
   cv_imrd  (:,:,:,:) = 0.
   albrf    (:,:,:,:) = 0.
   albwl    (:,:,:,:) = 0.
   albimrd  (:,:,:,:) = 0.

   alb_rf   (:,:,:,:,:) = 0.
   alb_wl   (:,:,:,:,:) = 0.
   alb_imrd (:,:,:,:,:) = 0.
   alb_perd (:,:,:,:,:) = 0.     
   ! ENDIF
#endif

   allocate( gfcc_tc (nxy, nxy) )
   allocate( gedi_th (nxy, nxy) )
   allocate( gl30_wt (nxy, nxy) )
   allocate( modur   (nxy, nxy) )
   allocate( wtrf    (nxy, nxy) )
   allocate( htrf    (nxy, nxy) )

   allocate( latso     (lat_points) )
   allocate( lonso     (lon_points) )
   allocate( area      (lon_points, lat_points) )
   allocate( hgt       (lon_points, lat_points) )
   allocate( avg       (lon_points, lat_points) )
   ! initialization
   hgt      (:,:)   = 0.
   wtrf     (:,:)   = 0.
   htrf     (:,:)   = 0.
   avg      (:,:)   = 0.
   tc       (:,:,:) = 0.
   urwt     (:,:,:) = 0.
   htop     (:,:,:) = 0.
   pct_tc   (:,:,:) = 0.
   pct_ur   (:,:,:) = 0.
   pct_urwt (:,:,:) = 0.
   htop_ur  (:,:,:) = 0.
   ur_dc    (:,:,:) = 0.

   ur_lai   (:,:,:,:) = 0.
   ur_sai   (:,:,:,:) = 0.
   wgt_lai  (:,:,:,:) = 0.
   wgt_sai  (:,:,:,:) = 0.

   cont = 0.
   sumth= 0.

   x_delta = (edgee-edgew)*1._r8/lon_points*1._r8
   y_delta = (edgen-edges)*1._r8/lat_points*1._r8

   DO i = 1, lon_points
      lonso(i) = edgew + (i-0.5)*x_delta
      IF (lonso(i) > 180) THEN
         lonso(i) = lonso(i) - 360
      ENDIF
   ENDDO

   DO i = 1, lat_points
      latso(i) = edgen - (i-0.5)*y_delta
   ENDDO

#ifdef USE_LCZ
   ldelta = 5._r8/(hxy*1._r8)
#else
   ldelta = 5._r8/(nxy*1._r8)
#endif

#ifndef USE_LCZ
   CALL nccheck( nf90_open(TRIM(dir_rawdata)//'urban_5x5/urban_properties.nc', nf90_nowrite, ncid) )

   CALL nccheck( nf90_inq_varid(ncid, "CANYON_HWR"         , hwr_canvid  ) )
   CALL nccheck( nf90_inq_varid(ncid, "WTLUNIT_ROOF"       , wt_rfvid    ) )
   CALL nccheck( nf90_inq_varid(ncid, "WTROAD_PERV"        , wt_rdvid    ) )
   CALL nccheck( nf90_inq_varid(ncid, "EM_ROOF"            , em_rfvid    ) )
   CALL nccheck( nf90_inq_varid(ncid, "EM_WALL"            , em_wlvid    ) )
   CALL nccheck( nf90_inq_varid(ncid, "EM_IMPROAD"         , em_imrdvid  ) )
   CALL nccheck( nf90_inq_varid(ncid, "EM_PERROAD"         , em_perdvid  ) )
   CALL nccheck( nf90_inq_varid(ncid, "ALB_ROOF"           , alb_rfvid   ) )
   CALL nccheck( nf90_inq_varid(ncid, "ALB_WALL"           , alb_wlvid   ) )
   CALL nccheck( nf90_inq_varid(ncid, "ALB_IMPROAD"        , alb_imrdvid ) )
   CALL nccheck( nf90_inq_varid(ncid, "ALB_PERROAD"        , alb_perdvid ) )
   CALL nccheck( nf90_inq_varid(ncid, "HT_ROOF"            , ht_rfvid    ) )
   CALL nccheck( nf90_inq_varid(ncid, "TK_ROOF"            , tk_rfvid    ) )
   CALL nccheck( nf90_inq_varid(ncid, "TK_WALL"            , tk_wlvid    ) )
   CALL nccheck( nf90_inq_varid(ncid, "TK_IMPROAD"         , tk_imrdvid  ) )
   CALL nccheck( nf90_inq_varid(ncid, "CV_ROOF"            , cv_rfvid    ) )
   CALL nccheck( nf90_inq_varid(ncid, "CV_WALL"            , cv_wlvid    ) )
   CALL nccheck( nf90_inq_varid(ncid, "CV_IMPROAD"         , cv_imrdvid  ) )
   CALL nccheck( nf90_inq_varid(ncid, "THICK_ROOF"         , th_rfvid    ) )
   CALL nccheck( nf90_inq_varid(ncid, "THICK_WALL"         , th_wlvid    ) )
   CALL nccheck( nf90_inq_varid(ncid, "T_BUILDING_MIN"     , tbminvid    ) )
   CALL nccheck( nf90_inq_varid(ncid, "T_BUILDING_MAX"     , tbmaxvid    ) )

   CALL nccheck( nf90_get_var(ncid, hwr_canvid  , hwrcan  ) )
   CALL nccheck( nf90_get_var(ncid, wt_rfvid    , ncar_wt ) )
   CALL nccheck( nf90_get_var(ncid, wt_rdvid    , wtrd    ) )
   CALL nccheck( nf90_get_var(ncid, em_rfvid    , emrf    ) )
   CALL nccheck( nf90_get_var(ncid, em_wlvid    , emwl    ) )
   CALL nccheck( nf90_get_var(ncid, em_imrdvid  , emimrd  ) )
   CALL nccheck( nf90_get_var(ncid, em_perdvid  , emperd  ) )
   CALL nccheck( nf90_get_var(ncid, alb_rfvid   , albrf   ) )
   CALL nccheck( nf90_get_var(ncid, alb_wlvid   , albwl   ) )
   CALL nccheck( nf90_get_var(ncid, alb_imrdvid , albimrd ) )
   CALL nccheck( nf90_get_var(ncid, alb_perdvid , albperd ) )
   CALL nccheck( nf90_get_var(ncid, ht_rfvid    , ncar_ht ) )
   CALL nccheck( nf90_get_var(ncid, tk_rfvid    , tkrf    ) )
   CALL nccheck( nf90_get_var(ncid, tk_wlvid    , tkwl    ) )
   CALL nccheck( nf90_get_var(ncid, tk_imrdvid  , tkimrd  ) )
   CALL nccheck( nf90_get_var(ncid, cv_rfvid    , cvrf    ) )
   CALL nccheck( nf90_get_var(ncid, cv_wlvid    , cvwl    ) )
   CALL nccheck( nf90_get_var(ncid, cv_imrdvid  , cvimrd  ) )
   CALL nccheck( nf90_get_var(ncid, th_rfvid    , thrf    ) )
   CALL nccheck( nf90_get_var(ncid, th_wlvid    , thwl    ) )
   CALL nccheck( nf90_get_var(ncid, tbminvid    , tbmin   ) )
   CALL nccheck( nf90_get_var(ncid, tbmaxvid    , tbmax   ) )

   CALL nccheck( nf90_close(ncid) )
#endif

   sreglat = 90 - int((90.-edgen)/5.)*5
   ereglat = 90 - int((90.-edges-0.5*ldelta)/5.)*5

   sreglon = -180 + int((edgew+180)/5.)*5
   ereglon = -180 + int((edgee+180-0.5*ldelta)/5.)*5
   
   DO reglat = sreglat, ereglat, -5
      DO reglon_ = sreglon, ereglon, 5

         reglon = reglon_
         IF (reglon_ > 180) reglon = reglon_ - 360
         reg(1) = reglat
         reg(2) = reglon
         reg(3) = reglat - 5
         reg(4) = reglon + 5

         write(reg1, "(i4)") reg(1)
         write(reg2, "(i4)") reg(2)
         write(reg3, "(i4)") reg(3)
         write(reg4, "(i4)") reg(4)

         ! PRINT*, "Processing ",TRIM(creg(1))//'_'//TRIM(creg(2))//'_'//TRIM(creg(3))//'_'//TRIM(creg(4))
         lndname = TRIM(dir_rawdata)//'srf_5x5/RG_'//TRIM(adjustL(reg1))//'_'//&
                   TRIM(adjustL(reg2))//'_'//TRIM(adjustL(reg3))//'_'//TRIM(adjustL(reg4))//'.MOD'//TRIM(cyear)//'.nc'

         inquire (file=lndname, exist=fileExists)
         IF (fileExists) THEN
            ! read 500m modis urban data
            
            CALL nccheck( nf90_open(trim(lndname), nf90_nowrite    , ncid       ) )

            CALL nccheck( nf90_inq_varid(ncid, "PCT_URBAN"     , upftvid    ) )
            CALL nccheck( nf90_inq_varid(ncid, "HTOP"          , gedi_thvid ) )

            CALL nccheck( nf90_get_var(ncid, upftvid    , modur   ) )
            CALL nccheck( nf90_get_var(ncid, gedi_thvid , gedi_th ) )

            CALL nccheck( nf90_close(ncid) )
         ELSE
            print*, 'Please check ', lndname
            CYCLE
         ENDIF

         lndname = TRIM(dir_rawdata)//'urban_5x5/RG_'//TRIM(adjustL(reg1))//'_'//&
                    TRIM(adjustL(reg2))//'_'//TRIM(adjustL(reg3))//'_'//TRIM(adjustL(reg4))//'.SRF'//TRIM(cyear)//'.nc' !TRIM(year)//'.nc'

         CALL nccheck( nf90_open(TRIM(lndname), nf90_nowrite, ncid) )

         CALL nccheck( nf90_inq_varid(ncid, "PCT_Tree"   , gfcc_tcvid) )
         !CALL nccheck( nf90_inq_varid(ncid, "Hgt_Tree"   , gedi_thvid) )
         CALL nccheck( nf90_inq_varid(ncid, "PCT_Water"  , gl30_wtvid) )

         CALL nccheck( nf90_get_var(ncid, gfcc_tcvid, gfcc_tc) )
         !CALL nccheck( nf90_get_var(ncid, gedi_thvid, gedi_th) )
         CALL nccheck( nf90_get_var(ncid, gl30_wtvid, gl30_wt) )

         CALL nccheck( nf90_close(ncid) )

         ! read urban LAI
         lndname = TRIM(dir_rawdata)//'building_5x5/RG_'//TRIM(adjustL(reg1))//'_'//&
                    TRIM(adjustL(reg2))//'_'//TRIM(adjustL(reg3))//'_'//TRIM(adjustL(reg4))//'.BLD.nc'
         CALL nccheck( nf90_open(lndname, nf90_nowrite, ncid) )

         CALL nccheck( nf90_inq_varid(ncid, "HT_ROOF"     , ht_rfvid) )
         CALL nccheck( nf90_get_var  (ncid, ht_rfvid      , htrf    ) )
         CALL nccheck( nf90_inq_varid(ncid, "WTLUNIT_ROOF", wt_rfvid) )
         CALL nccheck( nf90_get_var  (ncid, wt_rfvid      , wtrf    ) )

         CALL nccheck( nf90_close(ncid) )

#ifdef USE_LCZ
         lndname = TRIM(dir_rawdata)//'lcz_5x5/RG_'//TRIM(adjustL(reg1))//'_'//&
                    TRIM(adjustL(reg2))//'_'//TRIM(adjustL(reg3))//'_'//TRIM(adjustL(reg4))//'.LCZ.nc'

         PRINT*, ">>> Processing file ", trim(lndname), "..."
         CALL nccheck( nf90_open(lndname, nf90_nowrite, ncid) )

         CALL nccheck( nf90_inq_varid(ncid, "LCZ"     , ur_clssvid) )
         CALL nccheck( nf90_get_var  (ncid, ur_clssvid, lcz       ) )

         CALL nccheck( nf90_close(ncid) )

         lndname = TRIM(dir_rawdata)//'lai_5x5/RG_'//TRIM(adjustL(reg1))//'_'//&
                    TRIM(adjustL(reg2))//'_'//TRIM(adjustL(reg3))//'_'//TRIM(adjustL(reg4))//'.LCZLAI_v1'//TRIM(cyear)//'.nc'
         CALL nccheck( nf90_open(lndname , nf90_nowrite    , ncid  ) )

         CALL nccheck( nf90_inq_varid(ncid, "URBAN_TREE_LAI"  , laivid) )
         CALL nccheck( nf90_inq_varid(ncid, "URBAN_TREE_SAI"  , saivid) )
         CALL nccheck( nf90_get_var  (ncid, laivid            , hlai  ) )
         CALL nccheck( nf90_get_var  (ncid, saivid            , hsai  ) )

         CALL nccheck( nf90_close(ncid) )

         DO i = 1, hxy
            hlats(i) = reg(1) - i*ldelta
            hlatn(i) = reg(1) - (i-1)*ldelta
            hlonw(i) = reg(2) + (i-1)*ldelta
            hlone(i) = reg(2) + i*ldelta
         ENDDO

         DO i = 1, hxy
            dx = (hlone(1)-hlonw(1))*deg2rad
            dy = sin(hlatn(i)*deg2rad) - sin(hlats(i)*deg2rad)
            harea(:,i) = dx*dy*re*re
         ENDDO

         si = 1; ei = hxy
         sj = 1; ej = hxy

         IF (reglat == sreglat) si = int((reg(1)-edgen)*hxy/5)+1
         IF (reglon == sreglon) sj = int((edgew-reg(2))*hxy/5)+1
         IF (reglat == ereglat) ei = int((reg(1)-edges)*hxy/5)
         IF (reglon == ereglon) ej = int((edgee-reg(2))*hxy/5)
         IF (ei < si) ei = si
         IF (ej < sj) ej = sj

         DO i = si, ei
            DO j = sj, ej
               ! calculate io, jo
               !io = NINT((hlats(i)+ldelta/2+ 90.)/y_delta+0.5)
               !io = nyo+1-io
               ii = CEILING(i*1./5)
               jj = CEILING(j*1./5)

               IF (x_delta>0 .and. y_delta>0) THEN
                  io = NINT((90.-(hlats(i)+ldelta/2))/y_delta+0.5)
                  jo = NINT((hlonw(j)+ldelta/2+180.)/x_delta+0.5)
               ELSE
                  io = 1
                  jo = 1
               ENDIF

               IF (gedi_th(jj,ii) > 0) THEN
                  hgt(jo,io) = hgt(jo,io) + gedi_th(jj,ii)*harea(j,i)
                  avg(jo,io) = avg(jo,io) + harea(j,i)
               ENDIF

               IF (lcz(j,i)>0 .and. lcz(j,i)<=10) THEN
                  IF (modur(jj,ii)<=0 .and. lcz(j,i) >0) THEN
                     print*, 'lcz= ',lcz(j,i)
                  ENDIF
                  IF (modur(jj,ii) > 0) THEN                           
                     inx = int(lcz(j,i))
                     ur_dc(jo,io,inx) = ur_dc(jo,io,inx) + harea(j,i)
                     
                     IF (htrf(j,i) > 0) THEN
                        ht_rf(jo,io,3) = ht_rf(jo,io,3) + htrf(j,i)*harea(j,i)
                     ENDIF

                     IF (wtrf(j,i) > 0) THEN
                        wt_rf(jo,io,3) = wt_rf(jo,io,3) + wtrf(j,i)*harea(j,i)
                     ENDIF
                     ! 加权：
                     ! 粗网格城市水体(植被)覆盖度=粗网格城市水体(植被)覆盖度+500m城市格点水体(植被)覆盖度*500m城市格点面积
                     ! 加权系数；粗网格城市格点面积
                     IF (gl30_wt(jj,ii) > 0.) THEN
                        urwt (jo,io,inx) = urwt (jo,io,inx) + gl30_wt(jj,ii)*harea(j,i)
                     ENDIF
                     IF (gfcc_tc(jj,ii) > 0) THEN
                        tc   (jo,io,inx) = tc   (jo,io,inx) + gfcc_tc(jj,ii)*harea(j,i)
                        ! 树高加权
                        ! 粗网格城市树高=粗网格城市树高+500m城市格点植被覆盖度*城市格点树高*城市格点面积
                        ! 加权系数:城市格点植被覆盖度*城市面积
                        DO mm = 1, 12
                           IF (hlai(j,i,mm) > 0) THEN
                              ur_lai (jo,io,inx,mm) = ur_lai (jo,io,inx,mm) + hlai(j,i,mm)*harea(j,i)*gfcc_tc(jj,ii)
                              wgt_lai(jo,io,inx,mm) = wgt_lai(jo,io,inx,mm) + harea(j,i)*gfcc_tc(jj,ii)
                           ENDIF

                           IF (hsai(j,i,mm) > 0) THEN
                              ur_sai (jo,io,inx,mm) = ur_sai (jo,io,inx,mm) + hsai(j,i,mm)*harea(j,i)*gfcc_tc(jj,ii)
                              wgt_sai(jo,io,inx,mm) = wgt_sai(jo,io,inx,mm) + harea(j,i)*gfcc_tc(jj,ii)
                           ENDIF
                        ENDDO
                        IF (gedi_th(jj,ii) > 0) THEN
                           htop(jo,io,inx) = htop(jo,io,inx) + gedi_th(jj,ii)*gfcc_tc(jj,ii)*harea(j,i)
                           wgt_top(jo,io,inx) = wgt_top(jo,io,inx) + gfcc_tc(jj,ii)*harea(j,i)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
#else
         ! read NCAR 1km urban class and region id
         lndname = TRIM(dir_rawdata)//'urban_5x5/RG_'//TRIM(adjustL(reg1))//'_'//&
                   TRIM(adjustL(reg2))//'_'//TRIM(adjustL(reg3))//'_'//TRIM(adjustL(reg4))//'.NCAR.nc'
         
         PRINT*, ">>> Processing file ", trim(lndname), "..." 
         
         CALL nccheck( nf90_open(lndname, nf90_nowrite, ncid) )

         CALL nccheck( nf90_inq_varid(ncid, "URBAN_DENSITY_CLASS", ur_clssvid) )
         CALL nccheck( nf90_inq_varid(ncid, "REGION_ID"          , ur_rgvid  ) )

         CALL nccheck( nf90_get_var(ncid, ur_clssvid  , urden ) )
         CALL nccheck( nf90_get_var(ncid, ur_rgvid    , urrgid) )

         CALL nccheck( nf90_close(ncid) )

         lndname = TRIM(dir_rawdata)//'lai_5x5/RG_'//TRIM(adjustL(reg1))//'_'//&
                    TRIM(adjustL(reg2))//'_'//TRIM(adjustL(reg3))//'_'//TRIM(adjustL(reg4))//'.UrLAI_v3'//TRIM(cyear)//'.nc'

         CALL nccheck( nf90_open(lndname  , nf90_nowrite    , ncid  ) )

         CALL nccheck( nf90_inq_varid(ncid, "URBAN_TREE_LAI", laivid) )
         CALL nccheck( nf90_inq_varid(ncid, "URBAN_TREE_SAI", saivid) )
         CALL nccheck( nf90_get_var  (ncid, laivid          , hlai  ) )
         CALL nccheck( nf90_get_var  (ncid, saivid          , hsai  ) )

         CALL nccheck( nf90_close(ncid) )

         lndname = TRIM(dir_rawdata)//'building_5x5/RG_'//TRIM(adjustL(reg1))//'_'//&
                    TRIM(adjustL(reg2))//'_'//TRIM(adjustL(reg3))//'_'//TRIM(adjustL(reg4))//'.BLD.nc'
         CALL nccheck( nf90_open(lndname, nf90_nowrite, ncid) )

         CALL nccheck( nf90_inq_varid(ncid, "HT_ROOF"     , ht_rfvid) )
         CALL nccheck( nf90_get_var  (ncid, ht_rfvid      , htrf    ) )
         CALL nccheck( nf90_inq_varid(ncid, "WTLUNIT_ROOF", wt_rfvid) )
         CALL nccheck( nf90_get_var  (ncid, wt_rfvid      , wtrf    ) )

         CALL nccheck( nf90_close(ncid) )

         ! calculate the edge of small grids(500m)
         DO i = 1, nxy
            hlats(i) = reg(1) - i*ldelta
            hlatn(i) = reg(1) - (i-1)*ldelta
            hlonw(i) = reg(2) + (i-1)*ldelta
            hlone(i) = reg(2) + i*ldelta
         ENDDO

         DO i = 1, nxy
            dx = (hlone(1)-hlonw(1))*deg2rad
            dy = sin(hlatn(i)*deg2rad) - sin(hlats(i)*deg2rad)
            harea(:,i) = dx*dy*re*re
         ENDDO

         si = 1; ei = nxy
         sj = 1; ej = nxy

         IF (reglat == sreglat) si = int((reg(1)-edgen)*nxy/5)+1
         IF (reglon == sreglon) sj = int((edgew-reg(2))*nxy/5)+1
         IF (reglat == ereglat) ei = int((reg(1)-edges)*nxy/5)
         IF (reglon == ereglon) ej = int((edgee-reg(2))*nxy/5)
         IF (ei < si) ei = si
         IF (ej < sj) ej = sj

#ifdef USE_POINT_DATA
         DO i = 1, nxy
            DO j= 1, nxy
               IF (gedi_th(j,i)>0 .and. modur(j,i)>0) THEN
                  hgt(1,1) = hgt(1,1) + gedi_th(j,i)*harea(j,i)*modur(j,i)/100
                  avg(1,1) = avg(1,1) + harea(j,i)*modur(j,i)/100
               ENDIF
            ENDDO
         ENDDO
#endif
         DO i = si, ei
            DO j = sj, ej
               ! calculate io, jo
               !io = NINT((hlats(i)+ldelta/2+ 90.)/y_delta+0.5)
               !io = nyo+1-io

               IF (x_delta>0 .and. y_delta>0) THEN
                  io = NINT((90.-(hlats(i)+ldelta/2))/y_delta+0.5)
                  jo = NINT((hlonw(j)+ldelta/2+180.)/x_delta+0.5)
               ELSE
                  print*, 'Site point is ', (hlats(i)+ldelta/2), (hlonw(j)+ldelta/2)  ! for debug
                  io = 1
                  jo = 1
               ENDIF
               ! 聚合Simard树高
               ! 加权：
               ! 粗网格树高=粗网格树高+simard树高*1km格点面积
               ! 加权系数：粗网格格点面积
#ifndef USE_POINT_DATA
               IF (gedi_th(j,i) > 0) THEN
                  hgt(jo,io) = hgt(jo,io) + gedi_th(j,i)*harea(j,i)
                  avg(jo,io) = avg(jo,io) + harea(j,i)
               ENDIF
#endif
#ifndef USE_POINT_DATA
               IF (modur(j,i) > 0) THEN
#endif
                  IF (urden(j,i) > 0) THEN
                     ! IF (urden(j,i) == 1) THEN
                  
                     ! 加权：
                     ! 粗网格城市水体(植被)覆盖度=粗网格城市水体(植被)覆盖度+500m城市格点水体(植被)覆盖度*500m城市格点面积
                     ! 加权系数；粗网格城市格点面积
                     ! check for FI-Torni
                     IF (modur(j,i)<=0 .and. urden(j,i)>0) THEn
                        print*, 'MODIS and NCAR inconsistance'
                     ENDIF
                     inx = int(urden(j,i))
                     IF (gl30_wt(j,i) > 0.) THEN
                        urwt (jo,io,inx) = urwt (jo,io,inx) + gl30_wt(j,i)*harea(j,i)
                     ENDIF
                     IF (gfcc_tc(j,i) > 0.) THEN

                        tc(jo,io,inx) = tc(jo,io,inx) + gfcc_tc(j,i)*harea(j,i)

                        DO ii = 1, 12
                           IF (hlai(j,i,ii) > 0) THEN
                              ur_lai (jo,io,inx,ii) = ur_lai (jo,io,inx,ii) + hlai(j,i,ii)*harea(j,i)*gfcc_tc(j,i)
                              wgt_lai(jo,io,inx,ii) = wgt_lai(jo,io,inx,ii) + harea(j,i)*gfcc_tc(j,i)
                           ENDIF

                           IF (hsai(j,i,ii) > 0) THEN
                              ur_sai (jo,io,inx,ii) = ur_sai (jo,io,inx,ii) + hsai(j,i,ii)*harea(j,i)*gfcc_tc(j,i)
                              wgt_sai(jo,io,inx,ii) = wgt_sai(jo,io,inx,ii) + harea(j,i)*gfcc_tc(j,i)
                           ENDIF
                        ENDDO
                        ! 树高加权
                        ! 粗网格城市树高=粗网格城市树高+500m城市格点植被覆盖度*城市格点树高*城市格点面积
                        ! 加权系数:城市格点植被覆盖度*城市面积
                        IF (gedi_th(j,i) > 0) THEN
                           htop   (jo,io,inx) = htop   (jo,io,inx) + gedi_th(j,i)*gfcc_tc(j,i)*harea(j,i)
                           wgt_top(jo,io,inx) = wgt_top(jo,io,inx) + gfcc_tc(j,i)*harea(j,i)
                        ENDIF
                        ! ENDIF
                     ENDIF

                     uxid             = urrgid(j,i)
                     ur_dc(jo,io,inx) = ur_dc(jo,io,inx) + harea(j,i)
                     
                     IF (htrf(j,i) > 0) THEN
                        ht_rf(jo,io,inx) = ht_rf(jo,io,inx) + htrf(j,i)*harea(j,i)
                     ELSE
                        ht_rf(jo,io,inx) = ht_rf(jo,io,inx) + ncar_ht(inx,uxid)*harea(j,i)
                     ENDIF

                     IF (wtrf(j,i) > 0) THEN
                        wt_rf(jo,io,inx) = wt_rf(jo,io,inx) + wtrf(j,i)*harea(j,i)
                     ELSE
                        wt_rf(jo,io,inx) = wt_rf(jo,io,inx) + ncar_wt(inx,uxid)*harea(j,i)
                     ENDIF
                     
                     hwr_can  (jo,io,inx) = hwr_can  (jo,io,inx) + hwrcan  (inx,uxid)*harea(j,i)
                     !wt_rf    (jo,io,inx) = wt_rf    (jo,io,inx) + wtrf    (inx,uxid)*harea(j,i)
                     wt_rd    (jo,io,inx) = wt_rd    (jo,io,inx) + wtrd    (inx,uxid)*harea(j,i)
                     em_rf    (jo,io,inx) = em_rf    (jo,io,inx) + emrf    (inx,uxid)*harea(j,i)
                     em_wl    (jo,io,inx) = em_wl    (jo,io,inx) + emwl    (inx,uxid)*harea(j,i)
                     em_imrd  (jo,io,inx) = em_imrd  (jo,io,inx) + emimrd  (inx,uxid)*harea(j,i)
                     em_perd  (jo,io,inx) = em_perd  (jo,io,inx) + emperd  (inx,uxid)*harea(j,i)
                     th_rf    (jo,io,inx) = th_rf    (jo,io,inx) + thrf    (inx,uxid)*harea(j,i)
                     th_wl    (jo,io,inx) = th_wl    (jo,io,inx) + thwl    (inx,uxid)*harea(j,i)
                     tb_min   (jo,io,inx) = tb_min   (jo,io,inx) + tbmin   (inx,uxid)*harea(j,i)
                     tb_max   (jo,io,inx) = tb_max   (jo,io,inx) + tbmax   (inx,uxid)*harea(j,i)
                     !ht_rf    (jo,io,inx) = ht_rf    (jo,io,inx) + htrf    (inx,uxid)*harea(j,i)

                     alb_rf  (jo,io,inx,:,:) = alb_rf  (jo,io,inx,:,:) + albrf  (inx,uxid,:,:)*harea(j,i)
                     alb_wl  (jo,io,inx,:,:) = alb_wl  (jo,io,inx,:,:) + albwl  (inx,uxid,:,:)*harea(j,i)
                     alb_imrd(jo,io,inx,:,:) = alb_imrd(jo,io,inx,:,:) + albimrd(inx,uxid,:,:)*harea(j,i)
                     alb_perd(jo,io,inx,:,:) = alb_perd(jo,io,inx,:,:) + albperd(inx,uxid,:,:)*harea(j,i)

                     tk_rf  (jo,io,inx,:) = tk_rf  (jo,io,inx,:) + tkrf  (inx,uxid,:)*harea(j,i)
                     tk_wl  (jo,io,inx,:) = tk_wl  (jo,io,inx,:) + tkwl  (inx,uxid,:)*harea(j,i)
                     DO m = 1, 10
                        ! tkimrd与cvimrd有缺省值，计算需要跳过
                        IF (tkimrd(inx,uxid,m) .ne. -999) THEN
                           tk_imrd(jo,io,inx,m) = tk_imrd(jo,io,inx,m) + tkimrd(inx,uxid,m)*harea(j,i)
                        ENDIF
                        IF (cvimrd(inx,uxid,m) .ne. -999.) THEN
                           cv_imrd(jo,io,inx,m) = cv_imrd(jo,io,inx,m) + cvimrd(inx,uxid,m)*harea(j,i)
                        ENDIF
                     ENDDO
                     cv_rf  (jo,io,inx,:) = cv_rf  (jo,io,inx,:) + cvrf  (inx,uxid,:)*harea(j,i)
                     cv_wl  (jo,io,inx,:) = cv_wl  (jo,io,inx,:) + cvwl  (inx,uxid,:)*harea(j,i)
                  ENDIF

                  ! 根据MODIS城市覆盖对城市格点补充，并将其归类为MD urban
                  !IF (modur(j,i)>0 .and. urden(j,i)<=0) THEN
                  IF (urden(j,i) <= 0) THEN
                     ! print*, 'mod Processing'
                     IF (gl30_wt(j,i) > 0.) THEN
                        urwt (jo,io,3) = urwt (jo,io,3) + gl30_wt(j,i)*harea(j,i)*modur(j,i)/100
                     ENDIF
                     IF (gfcc_tc(j,i) > 0.) THEN
                        tc(jo,io,3) = tc(jo,io,3) + gfcc_tc(j,i)*harea(j,i)*modur(j,i)/100

                        DO ii = 1, 12
                           IF (hlai(j,i,ii) > 0) THEN
                              ur_lai (jo,io,3,ii) = ur_lai (jo,io,3,ii) + hlai(j,i,ii)*harea(j,i)*gfcc_tc(j,i)
                              wgt_lai(jo,io,3,ii) = wgt_lai(jo,io,3,ii) + harea(j,i)*gfcc_tc(j,i)
                           ENDIF

                           IF (hsai(j,i,ii) > 0) THEN
                              ur_sai (jo,io,3,ii) = ur_sai (jo,io,3,ii) + hsai(j,i,ii)*harea(j,i)*gfcc_tc(j,i)
                              wgt_sai(jo,io,3,ii) = wgt_sai(jo,io,3,ii) + harea(j,i)*gfcc_tc(j,i)
                           ENDIF
                        ENDDO

                        IF (gedi_th(j,i) > 0) THEN
                           htop(jo,io,3) = htop(jo,io,3) + gedi_th(j,i)*gfcc_tc(j,i)*harea(j,i)*modur(j,i)/100
                           wgt_top(jo,io,3) = wgt_top(jo,io,3) + harea(j,i)*gfcc_tc(j,i)*modur(j,i)/100
                        ENDIF
                     ENDIF

                     !ur_dc(jo,io,3) = ur_dc(jo,io,3) + harea(j,i)*modur(j,i)/100
                     uxid           = urrgid(j,i)
                     ! 部分格点MODIS与NCAR不一致(NCAR没有城市ID)，因此通过距离MODIS格点最近的NCAR城市ID赋值
                     IF (reg(1)==-45 .and. reg(3)==-50 .and. reg(2)==65 .and. reg(4)==70) THEN
                        uxid = 30
                     ENDIF
                     ! 城市建筑属性聚合
                     ! 加权：
                     ! 粗网格城市属性=粗网格城市属性+细网格城市属性*细网格面积*MODIS_PCT_URBAN
                     ! 加权系数：细网格面积(ur_dc)
#ifdef USE_POINT_DATA
                     print*, 'Urban class = ', urden(j,i)
                     ur_dc(jo,io,3) = ur_dc(jo,io,3) + harea(j,i)
                     
                     IF (htrf(j,i) > 0) THEN
                        ht_rf(jo,io,3) = ht_rf(jo,io,3) + htrf(j,i)*harea(j,i)
                     ELSE
                        ht_rf(jo,io,3) = ht_rf(jo,io,3) + ncar_ht(inx,3)*harea(j,i)
                     ENDIF

                     IF (wtrf(j,i) > 0) THEN
                        wt_rf(jo,io,3) = wt_rf(jo,io,3) + wtrf(j,i)*harea(j,i)
                     ELSE
                        wt_rf(jo,io,3) = wt_rf(jo,io,3) + ncar_wt(inx,3)*harea(j,i)
                     ENDIF
                     
                     hwr_can  (jo,io,3) = hwr_can  (jo,io,3) + hwrcan  (3,uxid)*harea(j,i)!*modur(j,i)/100
                     !wt_rf    (jo,io,3) = wt_rf    (jo,io,3) + wtrf    (3,uxid)*harea(j,i)*modur(j,i)/100
                     wt_rd    (jo,io,3) = wt_rd    (jo,io,3) + wtrd    (3,uxid)*harea(j,i)!*modur(j,i)/100
                     em_rf    (jo,io,3) = em_rf    (jo,io,3) + emrf    (3,uxid)*harea(j,i)!*modur(j,i)/100
                     em_wl    (jo,io,3) = em_wl    (jo,io,3) + emwl    (3,uxid)*harea(j,i)!*modur(j,i)/100
                     em_imrd  (jo,io,3) = em_imrd  (jo,io,3) + emimrd  (3,uxid)*harea(j,i)!*modur(j,i)/100
                     em_perd  (jo,io,3) = em_perd  (jo,io,3) + emperd  (3,uxid)*harea(j,i)!*modur(j,i)/100
                     th_rf    (jo,io,3) = th_rf    (jo,io,3) + thrf    (3,uxid)*harea(j,i)!*modur(j,i)/100
                     th_wl    (jo,io,3) = th_wl    (jo,io,3) + thwl    (3,uxid)*harea(j,i)!*modur(j,i)/100
                     tb_min   (jo,io,3) = tb_min   (jo,io,3) + tbmin   (3,uxid)*harea(j,i)!*modur(j,i)/100
                     tb_max   (jo,io,3) = tb_max   (jo,io,3) + tbmax   (3,uxid)*harea(j,i)!*modur(j,i)/100
                     !ht_rf    (jo,io,3) = ht_rf    (jo,io,3) + htrf    (3,uxid)*harea(j,i)*modur(j,i)/100

                     alb_rf  (jo,io,3,:,:) = alb_rf  (jo,io,3,:,:) + albrf  (3,uxid,:,:)*harea(j,i)!*modur(j,i)/100
                     alb_wl  (jo,io,3,:,:) = alb_wl  (jo,io,3,:,:) + albwl  (3,uxid,:,:)*harea(j,i)!*modur(j,i)/100
                     alb_imrd(jo,io,3,:,:) = alb_imrd(jo,io,3,:,:) + albimrd(3,uxid,:,:)*harea(j,i)!*modur(j,i)/100
                     alb_perd(jo,io,3,:,:) = alb_perd(jo,io,3,:,:) + albperd(3,uxid,:,:)*harea(j,i)!*modur(j,i)/100

                     tk_rf  (jo,io,3,:) = tk_rf  (jo,io,3,:) + tkrf  (3,uxid,:)*harea(j,i)!*modur(j,i)/100
                     tk_wl  (jo,io,3,:) = tk_wl  (jo,io,3,:) + tkwl  (3,uxid,:)*harea(j,i)!*modur(j,i)/100
                     DO m = 1, 10
                        IF (tkimrd(3,uxid,m) .ne. -999.) THEN
                           tk_imrd(jo,io,3,m) = tk_imrd(jo,io,3,m) + tkimrd(3,uxid,m)*harea(j,i)!*modur(j,i)/100
                        ENDIF
                        IF (cvimrd(3,uxid,m) .ne. -999.) THEN
                           cv_imrd(jo,io,3,m) = cv_imrd(jo,io,3,m) + cvimrd(3,uxid,m)*harea(j,i)!*modur(j,i)/100
                        ENDIF
                     ENDDO
                     cv_rf  (jo,io,3,:) = cv_rf  (jo,io,3,:) + cvrf  (3,uxid,:)*harea(j,i)!*modur(j,i)/100
                     cv_wl  (jo,io,3,:) = cv_wl  (jo,io,3,:) + cvwl  (3,uxid,:)*harea(j,i)!*modur(j,i)/100
#else
                     ur_dc(jo,io,3) = ur_dc(jo,io,3) + harea(j,i)*modur(j,i)/100

                     IF (htrf(j,i) > 0) THEN
                        ht_rf(jo,io,3) = ht_rf(jo,io,3) + htrf(j,i)*harea(j,i)*modur(j,i)/100
                     ELSE
                        ht_rf(jo,io,3) = ht_rf(jo,io,3) + ncar_ht(inx,3)*harea(j,i)*modur(j,i)/100
                     ENDIF

                     IF (wtrf(j,i) > 0) THEN
                        wt_rf(jo,io,3) = wt_rf(jo,io,3) + wtrf(j,i)*harea(j,i)
                     ELSE
                        wt_rf(jo,io,3) = wt_rf(jo,io,3) + ncar_wt(inx,3)*harea(j,i)
                     ENDIF

                     hwr_can  (jo,io,3) = hwr_can  (jo,io,3) + hwrcan  (3,uxid)*harea(j,i)*modur(j,i)/100
                     !wt_rf    (jo,io,3) = wt_rf    (jo,io,3) + wtrf    (3,uxid)*harea(j,i)*modur(j,i)/100
                     wt_rd    (jo,io,3) = wt_rd    (jo,io,3) + wtrd    (3,uxid)*harea(j,i)*modur(j,i)/100
                     em_rf    (jo,io,3) = em_rf    (jo,io,3) + emrf    (3,uxid)*harea(j,i)*modur(j,i)/100
                     em_wl    (jo,io,3) = em_wl    (jo,io,3) + emwl    (3,uxid)*harea(j,i)*modur(j,i)/100
                     em_imrd  (jo,io,3) = em_imrd  (jo,io,3) + emimrd  (3,uxid)*harea(j,i)*modur(j,i)/100
                     em_perd  (jo,io,3) = em_perd  (jo,io,3) + emperd  (3,uxid)*harea(j,i)*modur(j,i)/100
                     th_rf    (jo,io,3) = th_rf    (jo,io,3) + thrf    (3,uxid)*harea(j,i)*modur(j,i)/100
                     th_wl    (jo,io,3) = th_wl    (jo,io,3) + thwl    (3,uxid)*harea(j,i)*modur(j,i)/100
                     tb_min   (jo,io,3) = tb_min   (jo,io,3) + tbmin   (3,uxid)*harea(j,i)*modur(j,i)/100
                     tb_max   (jo,io,3) = tb_max   (jo,io,3) + tbmax   (3,uxid)*harea(j,i)*modur(j,i)/100
                     !ht_rf    (jo,io,3) = ht_rf    (jo,io,3) + htrf    (3,uxid)*harea(j,i)*modur(j,i)/100

                     alb_rf  (jo,io,3,:,:) = alb_rf  (jo,io,3,:,:) + albrf  (3,uxid,:,:)*harea(j,i)*modur(j,i)/100
                     alb_wl  (jo,io,3,:,:) = alb_wl  (jo,io,3,:,:) + albwl  (3,uxid,:,:)*harea(j,i)*modur(j,i)/100
                     alb_imrd(jo,io,3,:,:) = alb_imrd(jo,io,3,:,:) + albimrd(3,uxid,:,:)*harea(j,i)*modur(j,i)/100
                     alb_perd(jo,io,3,:,:) = alb_perd(jo,io,3,:,:) + albperd(3,uxid,:,:)*harea(j,i)*modur(j,i)/100

                     tk_rf  (jo,io,3,:) = tk_rf  (jo,io,3,:) + tkrf  (3,uxid,:)*harea(j,i)*modur(j,i)/100
                     tk_wl  (jo,io,3,:) = tk_wl  (jo,io,3,:) + tkwl  (3,uxid,:)*harea(j,i)*modur(j,i)/100
                     DO m = 1, 10
                        IF (tkimrd(3,uxid,m) .ne. -999.) THEN
                           tk_imrd(jo,io,3,m) = tk_imrd(jo,io,3,m) + tkimrd(3,uxid,m)*harea(j,i)*modur(j,i)/100
                        ENDIF
                        IF (cvimrd(3,uxid,m) .ne. -999.) THEN
                           cv_imrd(jo,io,3,m) = cv_imrd(jo,io,3,m) + cvimrd(3,uxid,m)*harea(j,i)*modur(j,i)/100
                        ENDIF
                     ENDDO
                     cv_rf  (jo,io,3,:) = cv_rf  (jo,io,3,:) + cvrf  (3,uxid,:)*harea(j,i)*modur(j,i)/100
                     cv_wl  (jo,io,3,:) = cv_wl  (jo,io,3,:) + cvwl  (3,uxid,:)*harea(j,i)*modur(j,i)/100
#endif
                  ENDIF
#ifndef USE_POINT_DATA
               ENDIF
#endif
            ENDDO
         ENDDO
#endif
      ENDDO
   ENDDO

   ! IF (USE_LCZ) THEN
#ifdef USE_LCZ
   DO i = 1, lat_points 
      DO j = 1, lon_points
         DO k = 1, 10
            IF (ur_dc(j,i,k) > 0) THEN
               ! calculate urban tree cover
               pct_tc  (j,i,k) = tc  (j,i,k) / ur_dc(j,i,k) !* 100
               ! calculate urban water cover
               pct_urwt(j,i,k) = urwt(j,i,k) / ur_dc(j,i,k) !* 100
               ht_rf    (j,i,k) = ht_rf    (j,i,k) / ur_dc(j,i,k)
               wt_rf    (j,i,k) = wt_rf    (j,i,k) / ur_dc(j,i,k)
               IF (wgt_top(j,i,k) > 0.) THEN
               ! calculate urban tree height
                  htop_ur (j,i,k) = htop(j,i,k) / wgt_top(j,i,k)!tc   (j,i,k-1)
               ENDIF
               DO ii = 1, 12
                  IF (wgt_lai(j,i,k,ii) > 0.) THEN
                  ! calculate urban tree height
                     ur_lai (j,i,k,ii) = ur_lai(j,i,k,ii) / wgt_lai(j,i,k,ii)!tc   (j,i,k-1)
                  ENDIF
                  IF (wgt_sai(j,i,k,ii) > 0.) THEN
                  ! calculate urban tree height
                     ur_sai (j,i,k,ii) = ur_sai(j,i,k,ii) / wgt_sai(j,i,k,ii)!tc   (j,i,k-1)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

         IF (avg(j,i) > 0) THEN
            hgt(j,i) = hgt(j,i) / avg(j,i)
         ENDIF

         sumur = sum(ur_dc(j,i,1:10)) !ur_dc(j,i,1) + ur_dc(j,i,2) + ur_dc(j,i,3)
         IF (sumur > 0.) THEN
            ! calculate 3 types urban cover, sum(pct_ur(j,i,:))=100
            DO k = 1, 10
               pct_ur(j,i,k) = ur_dc(j,i,k) / sumur * 100.
            ENDDO
         ENDIF
         ! 检查LCZ是否超过100%
         IF (sum(pct_ur(j,i,1:10)) > 1e-6 .and. abs(sum(pct_ur(j,i,1:10))-100) > 1e-3) THEN
            PRINT *, 'urban_pct > 100'
            PRINT *, pct_ur(j,i,1:10)
         ENDIF
      ENDDO
   ENDDO
   PRINT *, "********************************"
   DO i = 1, lat_points
      DO j = 1, lon_points
         DO k =1, 10
            ! nccheck for htop_ur of urban grid
            ! 如果该城市格点有植被覆盖却没有树高，则将聚合过程中生成的该格点的htop数据
            ! 赋值为城市树高
            IF (pct_tc(j,i,k) > 0 .and. htop_ur(j,i,k) == 0) THEN
               htop_ur(j,i,k) = hgt(j,i)
            ENDIF

            ! 如果经过上一步仍没有树高数据
            ! 则将该点同纬度的所有树高求平均赋值城市树高
            IF (pct_tc(j,i,k) > 0 .and. htop_ur(j,i,k) == 0) THEN
               DO p = 1, lon_points
                  IF (hgt(p,i) > 0) THEN
                     sumth = sumth + hgt(p,i)
                     cont  = cont  + 1
                  ENDIF
               ENDDO

               IF (cont > 0) THEN
                  htop_ur(j,i,k) = sumth/cont
               ENDIF
               cont  = 0
               sumth = 0
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   DO i=1, 12
      n_mon(i) = i
   ENDDO

   lndname = trim(dir_srfdata)//trim(cyear)//'/'//'urban_0.5x0.5.MOD.nc'
   print *, ">>> writing out surface data file ", trim(lndname)," ..."

   CALL nccheck( nf90_create (trim(lndname)   , NF90_NETCDF4, ncid      ) )
   CALL nccheck( nf90_def_dim(ncid, "lat"     , lat_points  , lat_dimid ) )
   CALL nccheck( nf90_def_dim(ncid, "lon"     , lon_points  , lon_dimid ) )
   CALL nccheck( nf90_def_dim(ncid, "LCZ_type", nlcz        , den_dimid ) )
   CALL nccheck( nf90_def_dim(ncid, "mon"     , mon         , mon_dimid ) )

   CALL nccheck( nf90_def_var(ncid, "lat"     , NF90_DOUBLE, lat_dimid , lat_vid ) )
   CALL nccheck( nf90_def_var(ncid, "lon"     , NF90_DOUBLE, lon_dimid , lon_vid ) )
   CALL nccheck( nf90_def_var(ncid, "month"   , NF90_INT   , mon_dimid , mon_vid ) )

   CALL nccheck( nf90_put_att(ncid, lat_vid , "long_name", "Latitude"        ) )
   CALL nccheck( nf90_put_att(ncid, lat_vid , "units"    , "degrees_north"   ) )
   CALL nccheck( nf90_put_att(ncid, lon_vid , "long_name", "Longitude"       ) )
   CALL nccheck( nf90_put_att(ncid, lon_vid , "units"    , "degrees_east"    ) )
   CALL nccheck( nf90_put_att(ncid, mon_vid , "long_name", "month"           ) )
   CALL nccheck( nf90_put_att(ncid, mon_vid , "units"    , "month"           ) )

   XY3D = (/lon_dimid, lat_dimid, den_dimid/)
   XY4D = (/lon_dimid, lat_dimid, den_dimid, mon_dimid/)

   CALL nccheck( nf90_def_var(ncid, "LCZ_TREE_PCT"  , NF90_DOUBLE, XY3D, pct_tcvid  ) )
   CALL nccheck( nf90_def_var(ncid, "LCZ_WATER_PCT" , NF90_DOUBLE, XY3D, pct_urwtvid) )
   CALL nccheck( nf90_def_var(ncid, "LCZ_TREE_TOP"  , NF90_DOUBLE, XY3D, htop_urvid ) )
   CALL nccheck( nf90_def_var(ncid, "LCZ_PCT"       , NF90_DOUBLE, XY3D, pct_urvid  ) )
   CALL nccheck( nf90_def_var(ncid, "LCZ_WT_ROOF"   , NF90_DOUBLE, XY3D, wt_rfvid   ) )
   CALL nccheck( nf90_def_var(ncid, "LCZ_HT_ROOF"   , NF90_DOUBLE, XY3D, ht_rfvid   ) )
   CALL nccheck( nf90_def_var(ncid, "LCZ_TREE_LAI"  , NF90_DOUBLE, XY4D, ur_laivid  ) )
   CALL nccheck( nf90_def_var(ncid, "LCZ_TREE_SAI"  , NF90_DOUBLE, XY4D, ur_saivid  ) )

   CALL nccheck( nf90_put_att(ncid, pct_urvid , "units"     , "%"                           ) )
   CALL nccheck( nf90_put_att(ncid, pct_urvid , "long_name" , "Percentage of each LCZ type" ) )
   CALL nccheck( nf90_put_att(ncid, pct_urvid , "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, pct_tcvid, "units"     , "%"                            ) )
   CALL nccheck( nf90_put_att(ncid, pct_tcvid, "long_name" , "Urban percent tree cover"     ) )
   CALL nccheck( nf90_put_att(ncid, pct_tcvid, "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, htop_urvid, "units"     , "m"                           ) )
   CALL nccheck( nf90_put_att(ncid, htop_urvid, "long_name" , "Urban tree top height"       ) )
   CALL nccheck( nf90_put_att(ncid, htop_urvid, "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, pct_urwtvid, "units"     , "%"                          ) )
   CALL nccheck( nf90_put_att(ncid, pct_urwtvid, "long_name" , "Urban percent water cover"  ) )
   CALL nccheck( nf90_put_att(ncid, pct_urwtvid, "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, wt_rfvid  , "units"     , "unitless"                    ) )
   CALL nccheck( nf90_put_att(ncid, wt_rfvid  , "long_name" , "fraction of roof"            ) )
   CALL nccheck( nf90_put_att(ncid, wt_rfvid  , "_FillValue", -999.) )

    CALL nccheck( nf90_put_att(ncid, ht_rfvid  , "units"     , "meters"                     ) )
   CALL nccheck( nf90_put_att(ncid, ht_rfvid  , "long_name" , "height of roof"              ) )
   CALL nccheck( nf90_put_att(ncid, ht_rfvid  , "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, ur_laivid, "units"     , "m^2/m^2"                      ) )
   CALL nccheck( nf90_put_att(ncid, ur_laivid, "long_name" , "Urban tree monthly lai"       ) )
   CALL nccheck( nf90_put_att(ncid, ur_laivid, "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, ur_saivid, "units"     , "m^2/m^2"                      ) )
   CALL nccheck( nf90_put_att(ncid, ur_saivid, "long_name" , "Urban tree monthly sai"       ) )
   CALL nccheck( nf90_put_att(ncid, ur_saivid, "_FillValue", -999.) )

   CALL nccheck( nf90_enddef(ncid) )

   CALL nccheck( nf90_inq_varid(ncid, "lat"          , urlat_vid  ) )
   CALL nccheck( nf90_put_var  (ncid, urlat_vid      , latso      ) )

   CALL nccheck( nf90_inq_varid(ncid, "lon"          , urlon_vid  ) )
   CALL nccheck( nf90_put_var  (ncid, urlon_vid      , lonso      ) )

   CALL nccheck( nf90_inq_varid(ncid, "month"        , mon_vid    ) )
   CALL nccheck( nf90_put_var  (ncid, mon_vid        , n_mon      ) )

   CALL nccheck( nf90_inq_varid(ncid, "LCZ_TREE_PCT" , pct_tcvid  ) )
   CALL nccheck( nf90_put_var  (ncid, pct_tcvid      , pct_tc     ) )

   CALL nccheck( nf90_inq_varid(ncid, "LCZ_WATER_PCT", pct_urwtvid) )
   CALL nccheck( nf90_put_var  (ncid, pct_urwtvid    , pct_urwt   ) )

   CALL nccheck( nf90_inq_varid(ncid, "LCZ_TREE_TOP" , htop_urvid ) )
   CALL nccheck( nf90_put_var  (ncid, htop_urvid     , htop_ur    ) )

   CALL nccheck( nf90_inq_varid(ncid, "LCZ_WT_ROOF"  , wt_rfvid   ) )
   CALL nccheck( nf90_put_var  (ncid, wt_rfvid       , wt_rf      ) )

   CALL nccheck( nf90_inq_varid(ncid, "LCZ_HT_ROOF"  , ht_rfvid   ) )
   CALL nccheck( nf90_put_var  (ncid, ht_rfvid       , ht_rf      ) )

   CALL nccheck( nf90_inq_varid(ncid, "LCZ_TREE_LAI" , ur_laivid  ) )
   CALL nccheck( nf90_put_var  (ncid, ur_laivid      , ur_lai     ) )

   CALL nccheck( nf90_inq_varid(ncid, "LCZ_TREE_SAI" , ur_saivid  ) )
   CALL nccheck( nf90_put_var  (ncid, ur_saivid      , ur_sai     ) )

   CALL nccheck( nf90_inq_varid(ncid, "LCZ_PCT"      , pct_urvid  ) )
   CALL nccheck( nf90_put_var  (ncid, pct_urvid      , pct_ur     ) )

   CALL nccheck( nf90_close(ncid) )
#else
   DO i = 1, lat_points
      DO j = 1, lon_points
         DO k = 1, 3
               IF (ur_dc(j,i,k) > 0) THEN
                  ! calculate urban tree cover
                  pct_tc  (j,i,k) = tc  (j,i,k) / ur_dc(j,i,k) !* 100
                  ! calculate urban water cover
                  pct_urwt(j,i,k) = urwt(j,i,k) / ur_dc(j,i,k) !* 100
                  IF (wgt_top(j,i,k) > 0.) THEN
                  ! calculate urban tree height
                     htop_ur (j,i,k) = htop(j,i,k) / wgt_top(j,i,k)!tc   (j,i,k-1)
                  ENDIF
                  DO ii = 1, 12
                     IF (wgt_lai(j,i,k,ii) > 0.) THEN
                     ! calculate urban tree height
                        ur_lai (j,i,k,ii) = ur_lai(j,i,k,ii) / wgt_lai(j,i,k,ii)!tc   (j,i,k-1)
                     ENDIF
                     IF (wgt_sai(j,i,k,ii) > 0.) THEN
                     ! calculate urban tree height
                        ur_sai (j,i,k,ii) = ur_sai(j,i,k,ii) / wgt_sai(j,i,k,ii)!tc   (j,i,k-1)
                     ENDIF
                  ENDDO
               
                  hwr_can  (j,i,k) = hwr_can  (j,i,k) / ur_dc(j,i,k)
                  wt_rf    (j,i,k) = wt_rf    (j,i,k) / ur_dc(j,i,k)
                  wt_rd    (j,i,k) = wt_rd    (j,i,k) / ur_dc(j,i,k)
                  em_rf    (j,i,k) = em_rf    (j,i,k) / ur_dc(j,i,k)
                  em_wl    (j,i,k) = em_wl    (j,i,k) / ur_dc(j,i,k)
                  em_imrd  (j,i,k) = em_imrd  (j,i,k) / ur_dc(j,i,k)
                  em_perd  (j,i,k) = em_perd  (j,i,k) / ur_dc(j,i,k)
                  th_rf    (j,i,k) = th_rf    (j,i,k) / ur_dc(j,i,k)
                  th_wl    (j,i,k) = th_wl    (j,i,k) / ur_dc(j,i,k)
                  tb_min   (j,i,k) = tb_min   (j,i,k) / ur_dc(j,i,k)
                  tb_max   (j,i,k) = tb_max   (j,i,k) / ur_dc(j,i,k)
                  ht_rf    (j,i,k) = ht_rf    (j,i,k) / ur_dc(j,i,k)
   
                  alb_rf  (j,i,k,:,:) = alb_rf  (j,i,k,:,:) / ur_dc(j,i,k)
                  alb_wl  (j,i,k,:,:) = alb_wl  (j,i,k,:,:) / ur_dc(j,i,k)
                  alb_imrd(j,i,k,:,:) = alb_imrd(j,i,k,:,:) / ur_dc(j,i,k)
                  alb_perd(j,i,k,:,:) = alb_perd(j,i,k,:,:) / ur_dc(j,i,k)
   
                  tk_rf  (j,i,k,:) = tk_rf  (j,i,k,:) / ur_dc(j,i,k)
                  tk_wl  (j,i,k,:) = tk_wl  (j,i,k,:) / ur_dc(j,i,k)
                  DO m = 1, 10
                     IF (tk_imrd(j,i,k,m) >= 0.) THEN
                        tk_imrd(j,i,k,m) = tk_imrd(j,i,k,m) / ur_dc(j,i,k)
                     ENDIF
                     IF (cv_imrd(j,i,k,m) >= 0.) THEN
                        cv_imrd(j,i,k,m) = cv_imrd(j,i,k,m) / ur_dc(j,i,k)
                     ENDIF
                  ENDDO
                  cv_rf  (j,i,k,:) = cv_rf  (j,i,k,:) / ur_dc(j,i,k)
                  cv_wl  (j,i,k,:) = cv_wl  (j,i,k,:) / ur_dc(j,i,k)
               ENDIF
         ENDDO

         IF (avg(j,i) > 0) THEN
            hgt(j,i) = hgt(j,i) / avg(j,i)
         ENDIF 

         sumur = ur_dc(j,i,1) + ur_dc(j,i,2) + ur_dc(j,i,3)
         IF (sumur > 0.) THEN
            ! calculate 3 types urban cover, sum(pct_ur(j,i,:))=100
            pct_ur(j,i,1) = ur_dc(j,i,1) / sumur * 100.
            pct_ur(j,i,2) = ur_dc(j,i,2) / sumur * 100.
            pct_ur(j,i,3) = ur_dc(j,i,3) / sumur * 100.
         ENDIF

         ! 检查3类城市是否超过100%
         IF (sum(pct_ur(j,i,1:3)) > 1e-6 .and. abs(sum(pct_ur(j,i,1:3))-100) > 1e-3) THEN
            PRINT *, 'urban_pct > 100'
            PRINT *, pct_ur(j,i,1:3)
         ENDIF
      ENDDO
   ENDDO   

   PRINT *, "********************************"   
   DO i = 1, lat_points
      DO j = 1, lon_points
         DO k =1, 3
            !nccheck for htop_ur of urban grid
            ! 如果该城市格点有植被覆盖却没有树高，则将聚合过程中生成的该格点的htop数据
            ! 赋值为城市树高
            IF (pct_tc(j,i,k) > 0 .and. htop_ur(j,i,k) == 0) THEN
               htop_ur(j,i,k) = hgt(j,i)
            ENDIF

            ! 如果经过上一步仍没有树高数据
            ! 则将该点同纬度的所有树高求平均赋值城市树高
            IF (pct_tc(j,i,k) > 0 .and. htop_ur(j,i,k) == 0) THEN
               DO p = 1, lon_points
                  IF (hgt(p,i) > 0) THEN
                     sumth = sumth + hgt(p,i)
                     cont  = cont  + 1
                  ENDIF
               ENDDO

               IF (cont > 0) THEN
                  htop_ur(j,i,k) = sumth/cont
               ENDIF

               cont  = 0
               sumth = 0
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   DO i = 1, 3
      n_den(i) = i
   ENDDO

   DO i = 1, 2
      n_nr(i) = i
   ENDDO

   DO i = 1, 2
      n_ns(i) = i
   ENDDO

   DO i = 1, 10
      n_ulev(i) = i
   ENDDO

   DO i = 1, 12
      n_mon(i) = i
   ENDDO

   lndname = trim(dir_srfdata)//trim(cyear)//'/urban_0.5x0.5.MOD.nc'
   print *, ">>> writing out surface data file ", trim(lndname), " ..."

   CALL nccheck( nf90_create(trim(lndname), NF90_NETCDF4, ncid) )

   CALL nccheck( nf90_def_dim(ncid, "lat"     , lat_points     , lat_dimid ) )
   CALL nccheck( nf90_def_dim(ncid, "lon"     , lon_points     , lon_dimid ) )
   CALL nccheck( nf90_def_dim(ncid, "numsolar", ns      , ns_dimid  ) )
   CALL nccheck( nf90_def_dim(ncid, "numrad"  , nr      , nr_dimid  ) )
   CALL nccheck( nf90_def_dim(ncid, "ulev"    , ulev    , ulev_dimid) )
   CALL nccheck( nf90_def_dim(ncid, "density" , den_clss, den_dimid ) )
   CALL nccheck( nf90_def_dim(ncid, "month"   , mon     , mon_dimid ) )

   CALL nccheck( nf90_def_var(ncid, "lat"     , NF90_DOUBLE, lat_dimid , lat_vid ) )
   CALL nccheck( nf90_def_var(ncid, "lon"     , NF90_DOUBLE, lon_dimid , lon_vid ) )
   CALL nccheck( nf90_def_var(ncid, "numsolar", NF90_INT  , ns_dimid  , ns_vid  ) )
   CALL nccheck( nf90_def_var(ncid, "numrad"  , NF90_INT  , nr_dimid  , nr_vid  ) )
   CALL nccheck( nf90_def_var(ncid, "ulev"    , NF90_INT  , ulev_dimid, lev_vid ) )
   CALL nccheck( nf90_def_var(ncid, "density" , NF90_INT  , den_dimid , den_vid ) )
   CALL nccheck( nf90_def_var(ncid, "month"   , NF90_INT  , mon_dimid , mon_vid ) )

   CALL nccheck( nf90_put_att(ncid, lat_vid , "long_name", "Latitude"        ) )
   CALL nccheck( nf90_put_att(ncid, lat_vid , "units"    , "degrees_north"   ) )
   CALL nccheck( nf90_put_att(ncid, lon_vid , "long_name", "Longitude"       ) )
   CALL nccheck( nf90_put_att(ncid, lon_vid , "units"    , "degrees_east"    ) )
   CALL nccheck( nf90_put_att(ncid, ns_vid  , "long_name", "solar band"      ) )
   CALL nccheck( nf90_put_att(ncid, ns_vid  , "units"    , "1-dir,2-diff"    ) )
   CALL nccheck( nf90_put_att(ncid, nr_vid  , "long_name", "radiation band"  ) )
   CALL nccheck( nf90_put_att(ncid, nr_vid  , "units"    , "1-vis,2-nir"     ) )
   CALL nccheck( nf90_put_att(ncid, lev_vid , "long_name", "urban layer"     ) )
   CALL nccheck( nf90_put_att(ncid, lev_vid , "units"    , "urban layer"     ) )
   CALL nccheck( nf90_put_att(ncid, den_vid , "long_name", "urban density"   ) )
   CALL nccheck( nf90_put_att(ncid, den_vid , "units"    , "1-tbd,2-hd,3-md" ) )
   CALL nccheck( nf90_put_att(ncid, mon_vid , "long_name", "month"           ) )
   CALL nccheck( nf90_put_att(ncid, mon_vid , "units"    , "month"           ) )

   XY2D = (/ lon_dimid, lat_dimid /)
   XY3D = (/ lon_dimid, lat_dimid, den_dimid /)
   CALL nccheck( nf90_def_var(ncid, "URBAN_TREE_PCT" , NF90_DOUBLE, XY3D, pct_tcvid  ) )
   CALL nccheck( nf90_def_var(ncid, "URBAN_WATER_PCT", NF90_DOUBLE, XY3D, pct_urwtvid) )
   CALL nccheck( nf90_def_var(ncid, "URBAN_TREE_TOP" , NF90_DOUBLE, XY3D, htop_urvid ) )
   CALL nccheck( nf90_def_var(ncid, "CANYON_HWR"     , NF90_DOUBLE, XY3D, hwr_canvid ) )
   CALL nccheck( nf90_def_var(ncid, "WTLUNIT_ROOF"   , NF90_DOUBLE, XY3D, wt_rfvid   ) )
   CALL nccheck( nf90_def_var(ncid, "WTROAD_PERV"    , NF90_DOUBLE, XY3D, wt_rdvid   ) )
   CALL nccheck( nf90_def_var(ncid, "EM_ROOF"        , NF90_DOUBLE, XY3D, em_rfvid   ) )
   CALL nccheck( nf90_def_var(ncid, "EM_WALL"        , NF90_DOUBLE, XY3D, em_wlvid   ) )
   CALL nccheck( nf90_def_var(ncid, "EM_IMPROAD"     , NF90_DOUBLE, XY3D, em_imrdvid ) )
   CALL nccheck( nf90_def_var(ncid, "EM_PERROAD"     , NF90_DOUBLE, XY3D, em_perdvid ) )
   CALL nccheck( nf90_def_var(ncid, "HT_ROOF"        , NF90_DOUBLE, XY3D, ht_rfvid   ) )
   CALL nccheck( nf90_def_var(ncid, "THICK_ROOF"     , NF90_DOUBLE, XY3D, th_rfvid   ) )
   CALL nccheck( nf90_def_var(ncid, "THICK_WALL"     , NF90_DOUBLE, XY3D, th_wlvid   ) )
   CALL nccheck( nf90_def_var(ncid, "T_BUILDING_MIN" , NF90_DOUBLE, XY3D, tbminvid   ) )
   CALL nccheck( nf90_def_var(ncid, "T_BUILDING_MAX" , NF90_DOUBLE, XY3D, tbmaxvid   ) )

   XY4D = (/ lon_dimid, lat_dimid, den_dimid, ulev_dimid /)
   UR3D = (/ lon_dimid, lat_dimid, den_dimid  /)
   CALL nccheck( nf90_def_var(ncid, "URBAN_PCT" , NF90_DOUBLE, UR3D, pct_urvid ) )
   CALL nccheck( nf90_def_var(ncid, "TK_ROOF"   , NF90_DOUBLE, XY4D, tk_rfvid  ) )
   CALL nccheck( nf90_def_var(ncid, "TK_WALL"   , NF90_DOUBLE, XY4D, tk_wlvid  ) )
   CALL nccheck( nf90_def_var(ncid, "TK_IMPROAD", NF90_DOUBLE, XY4D, tk_imrdvid) )
   CALL nccheck( nf90_def_var(ncid, "CV_ROOF"   , NF90_DOUBLE, XY4D, cv_rfvid  ) )
   CALL nccheck( nf90_def_var(ncid, "CV_WALL"   , NF90_DOUBLE, XY4D, cv_wlvid  ) )
   CALL nccheck( nf90_def_var(ncid, "CV_IMPROAD", NF90_DOUBLE, XY4D, cv_imrdvid) )

   UL3D = (/ lon_dimid, lat_dimid, den_dimid, mon_dimid /)
   CALL nccheck( nf90_def_var(ncid, "URBAN_TREE_LAI", NF90_DOUBLE, UL3D, ur_laivid ) )
   CALL nccheck( nf90_def_var(ncid, "URBAN_TREE_SAI", NF90_DOUBLE, UL3D, ur_saivid ) )

   XY5D = (/ lon_dimid, lat_dimid, den_dimid, nr_dimid, ns_dimid /)
   CALL nccheck( nf90_def_var(ncid, "ALB_ROOF"   , NF90_DOUBLE, XY5D, alb_rfvid   ) )
   CALL nccheck( nf90_def_var(ncid, "ALB_WALL"   , NF90_DOUBLE, XY5D, alb_wlvid   ) )
   CALL nccheck( nf90_def_var(ncid, "ALB_IMPROAD", NF90_DOUBLE, XY5D, alb_imrdvid ) )
   CALL nccheck( nf90_def_var(ncid, "ALB_PERROAD", NF90_DOUBLE, XY5D, alb_perdvid ) )

   CALL nccheck( nf90_put_att(ncid, pct_urvid , "units"     , "%"                                          ) )
   CALL nccheck( nf90_put_att(ncid, pct_urvid , "long_name" , "Percentage of each urban type (density)"    ) )
   CALL nccheck( nf90_put_att(ncid, pct_urvid , "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, pct_tcvid, "units"     , "%"                                           ) )
   CALL nccheck( nf90_put_att(ncid, pct_tcvid, "long_name" , "Urban percent tree cover"                    ) )
   CALL nccheck( nf90_put_att(ncid, pct_tcvid, "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, htop_urvid, "units"     , "m"                                          ) )
   CALL nccheck( nf90_put_att(ncid, htop_urvid, "long_name" , "Urban tree top height"                      ) )
   CALL nccheck( nf90_put_att(ncid, htop_urvid, "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, pct_urwtvid, "units"     , "%"                                         ) )
   CALL nccheck( nf90_put_att(ncid, pct_urwtvid, "long_name" , "Urban percent water cover"                 ) )
   CALL nccheck( nf90_put_att(ncid, pct_urwtvid, "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, ur_laivid, "units"     , "m^2/m^2"                                     ) )
   CALL nccheck( nf90_put_att(ncid, ur_laivid, "long_name" , "Urban tree monthly lai"                      ) )
   CALL nccheck( nf90_put_att(ncid, ur_laivid, "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, ur_saivid, "units"     , "m^2/m^2"                                     ) )
   CALL nccheck( nf90_put_att(ncid, ur_saivid, "long_name" , "Urban tree monthly sai"                      ) )
   CALL nccheck( nf90_put_att(ncid, ur_saivid, "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, hwr_canvid, "units"     , "unitless"                                   ) )
   CALL nccheck( nf90_put_att(ncid, hwr_canvid, "long_name" , "canyon height to width ratio"               ) )
   CALL nccheck( nf90_put_att(ncid, hwr_canvid, "_FillValue", -999.   ) )

   CALL nccheck( nf90_put_att(ncid, wt_rfvid  , "units"     , "unitless"                                   ) )
   CALL nccheck( nf90_put_att(ncid, wt_rfvid  , "long_name" , "fraction of roof"                           ) )
   CALL nccheck( nf90_put_att(ncid, wt_rfvid  , "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, wt_rdvid  , "units"     , "unitless"                                   ) )
   CALL nccheck( nf90_put_att(ncid, wt_rdvid  , "long_name" , "fraction of pervious road"                  ) )
   CALL nccheck( nf90_put_att(ncid, wt_rdvid  , "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, em_rfvid  , "units"     , "unitless"                                   ) )
   CALL nccheck( nf90_put_att(ncid, em_rfvid  , "long_name" , "emissivity of roof"                         ) )
   CALL nccheck( nf90_put_att(ncid, em_rfvid  , "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, em_wlvid  , "units"     , "unitless"                                   ) )
   CALL nccheck( nf90_put_att(ncid, em_wlvid  , "long_name" , "emissivity of wall"                         ) )
   CALL nccheck( nf90_put_att(ncid, em_wlvid  , "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, em_imrdvid, "units"     , "unitless"                                   ) )
   CALL nccheck( nf90_put_att(ncid, em_imrdvid, "long_name" , "emissivity of impervious road"              ) )
   CALL nccheck( nf90_put_att(ncid, em_imrdvid, "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, em_perdvid, "units"     , "unitless"                                   ) )
   CALL nccheck( nf90_put_att(ncid, em_perdvid, "long_name" , "emissivity of pervious road"                ) )
   CALL nccheck( nf90_put_att(ncid, em_perdvid, "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, ht_rfvid  , "units"     , "meters"                                     ) )
   CALL nccheck( nf90_put_att(ncid, ht_rfvid  , "long_name" , "height of roof"                             ) )
   CALL nccheck( nf90_put_att(ncid, ht_rfvid  , "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, th_rfvid  , "units"     , "m"                                          ) )
   CALL nccheck( nf90_put_att(ncid, th_rfvid  , "long_name" , "thickness of roof"                          ) )
   CALL nccheck( nf90_put_att(ncid, th_rfvid  , "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, th_wlvid  , "units"     , "m"                                          ) )
   CALL nccheck( nf90_put_att(ncid, th_wlvid  , "long_name" , "thickness of wall"                          ) )
   CALL nccheck( nf90_put_att(ncid, th_wlvid  , "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, tbminvid  , "units"     , "K"                                          ) )
   CALL nccheck( nf90_put_att(ncid, tbminvid  , "long_name" , "minimum interior building temperature"      ) )
   CALL nccheck( nf90_put_att(ncid, tbminvid  , "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, tbmaxvid  , "units"     , "K"                                          ) )
   CALL nccheck( nf90_put_att(ncid, tbmaxvid  , "long_name" , "maximum interior building temperature"      ) )
   CALL nccheck( nf90_put_att(ncid, tbmaxvid  , "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, tk_rfvid  , "units"     , "W/m*K"                                      ) )
   CALL nccheck( nf90_put_att(ncid, tk_rfvid  , "long_name" , "thermal conductivity of roof"               ) )
   CALL nccheck( nf90_put_att(ncid, tk_rfvid  , "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, tk_wlvid  , "units", "W/m*K"                                           ) )
   CALL nccheck( nf90_put_att(ncid, tk_wlvid  , "long_name", "thermal conductivity of wall"                ) )
   CALL nccheck( nf90_put_att(ncid, tk_wlvid  , "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, tk_imrdvid, "units"     , "W/m*K"                                      ) )
   CALL nccheck( nf90_put_att(ncid, tk_imrdvid, "long_name" , "thermal conductivity of impervious road"    ) )
   CALL nccheck( nf90_put_att(ncid, tk_imrdvid, "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, cv_rfvid  , "units"     , "J/m^3*K"                                    ) )
   CALL nccheck( nf90_put_att(ncid, cv_rfvid  , "long_name" , "volumetric heat capacity of roof"           ) )
   CALL nccheck( nf90_put_att(ncid, cv_rfvid  , "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, cv_wlvid  , "units"     , "J/m^3*K"                                    ) )
   CALL nccheck( nf90_put_att(ncid, cv_wlvid  , "long_name" , "volumetric heat capacity of wall"           ) )
   CALL nccheck( nf90_put_att(ncid, cv_wlvid  , "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, cv_imrdvid, "units"     , "J/m^3*K"                                    ) )
   CALL nccheck( nf90_put_att(ncid, cv_imrdvid, "long_name" , "volumetric heat capacity of impervious road") )
   CALL nccheck( nf90_put_att(ncid, cv_imrdvid, "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, alb_rfvid , "units"     , "unitless"                                   ) )
   CALL nccheck( nf90_put_att(ncid, alb_rfvid , "long_name" , "albedo of roof"                             ) )
   CALL nccheck( nf90_put_att(ncid, alb_rfvid , "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, alb_wlvid , "units"     , "unitless"                                   ) )
   CALL nccheck( nf90_put_att(ncid, alb_wlvid , "long_name" , "albedo of wall"                             ) )
   CALL nccheck( nf90_put_att(ncid, alb_wlvid , "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, alb_imrdvid, "units"     , "unitless"                                  ) )
   CALL nccheck( nf90_put_att(ncid, alb_imrdvid, "long_name" , "albedo of impervious road"                 ) )
   CALL nccheck( nf90_put_att(ncid, alb_imrdvid, "_FillValue", -999.) )

   CALL nccheck( nf90_put_att(ncid, alb_perdvid, "units"     , "unitless"                                  ) )
   CALL nccheck( nf90_put_att(ncid, alb_perdvid, "long_name" , "albedo of pervious road"                   ) )
   CALL nccheck( nf90_put_att(ncid, alb_perdvid, "_FillValue", -999.) )

   CALL nccheck( nf90_enddef(ncid) )
   
   CALL nccheck( nf90_inq_varid(ncid, "lat"            , urlat_vid  ) )
   CALL nccheck( nf90_put_var  (ncid, urlat_vid        , latso     ) )

   CALL nccheck( nf90_inq_varid(ncid, "lon"            , urlon_vid  ) )
   CALL nccheck( nf90_put_var  (ncid, urlon_vid        , lonso      ) )

   call nccheck( nf90_inq_varid(ncid, "numsolar"       , ns_vid     ) )
   CALL nccheck( nf90_put_var  (ncid, ns_vid           , n_ns       ) )

   CALL nccheck( nf90_inq_varid(ncid, "numrad"         , nr_vid     ) )
   CALL nccheck( nf90_put_var  (ncid, nr_vid           , n_nr       ) )

   CALL nccheck( nf90_inq_varid(ncid, "month"          , mon_vid    ) )
   CALL nccheck( nf90_put_var  (ncid, mon_vid          , n_mon      ) )

   CALL nccheck( nf90_inq_varid(ncid, "density"        , den_vid    ) )
   CALL nccheck( nf90_put_var  (ncid, den_vid          , n_den      ) )

   CALL nccheck( nf90_inq_varid(ncid, "ulev"           , lev_vid    ) )
   CALL nccheck( nf90_put_var  (ncid, lev_vid          , n_ulev     ) )

   CALL nccheck( nf90_inq_varid(ncid, "URBAN_TREE_PCT" , pct_tcvid  ) )
   CALL nccheck( nf90_put_var  (ncid, pct_tcvid        , pct_tc    ) )

   CALL nccheck( nf90_inq_varid(ncid, "URBAN_WATER_PCT", pct_urwtvid) )
   CALL nccheck( nf90_put_var  (ncid, pct_urwtvid      , pct_urwt  ) )

   CALL nccheck( nf90_inq_varid(ncid, "URBAN_TREE_TOP" , htop_urvid ) )
   CALL nccheck( nf90_put_var  (ncid, htop_urvid       , htop_ur   ) )

   CALL nccheck( nf90_inq_varid(ncid, "URBAN_TREE_LAI" , ur_laivid  ) )
   CALL nccheck( nf90_put_var  (ncid, ur_laivid        , ur_lai    ) )

   CALL nccheck( nf90_inq_varid(ncid, "URBAN_TREE_SAI" , ur_saivid  ) )
   CALL nccheck( nf90_put_var  (ncid, ur_saivid        , ur_sai    ) )

   CALL nccheck( nf90_inq_varid(ncid, "CANYON_HWR"     , hwr_canvid ) )
   CALL nccheck( nf90_put_var  (ncid,hwr_canvid        , hwr_can   ) )

   CALL nccheck( nf90_inq_varid(ncid, "WTLUNIT_ROOF"   , wt_rfvid   ) )
   CALL nccheck( nf90_put_var  (ncid, wt_rfvid         , wt_rf     ) )

   CALL nccheck( nf90_inq_varid(ncid, "WTROAD_PERV"    , wt_rdvid   ) )
   CALL nccheck( nf90_put_var  (ncid, wt_rdvid         , wt_rd     ) )

   CALL nccheck( nf90_inq_varid(ncid, "EM_ROOF"        , em_rfvid   ) )
   CALL nccheck( nf90_put_var  (ncid, em_rfvid         , em_rf     ) )

   CALL nccheck( nf90_inq_varid(ncid, "EM_WALL"        , em_wlvid   ) )
   CALL nccheck( nf90_put_var  (ncid, em_wlvid         , em_wl     ) )

   CALL nccheck( nf90_inq_varid(ncid, "EM_IMPROAD"     , em_imrdvid ) )
   CALL nccheck( nf90_put_var  (ncid, em_imrdvid       , em_imrd   ) )

   CALL nccheck( nf90_inq_varid(ncid, "EM_PERROAD"     , em_perdvid ) )
   CALL nccheck( nf90_put_var  (ncid, em_perdvid       , em_perd   ) )

   CALL nccheck( nf90_inq_varid(ncid, "HT_ROOF"        , ht_rfvid   ) )
   CALL nccheck( nf90_put_var  (ncid, ht_rfvid         , ht_rf     ) )

   CALL nccheck( nf90_inq_varid(ncid, "THICK_ROOF"     , th_rfvid   ) )
   CALL nccheck( nf90_put_var  (ncid, th_rfvid         , th_rf     ) )

   CALL nccheck( nf90_inq_varid(ncid, "THICK_WALL"     , th_wlvid   ) )
   CALL nccheck( nf90_put_var  (ncid, th_wlvid         , th_wl     ) )

   CALL nccheck( nf90_inq_varid(ncid, "T_BUILDING_MIN" , tbminvid   ) )
   CALL nccheck( nf90_put_var  (ncid, tbminvid         , tb_min    ) )

   CALL nccheck( nf90_inq_varid(ncid, "T_BUILDING_MAX" , tbmaxvid   ) )
   CALL nccheck( nf90_put_var  (ncid, tbmaxvid         , tb_max    ) )

   CALL nccheck( nf90_inq_varid(ncid, "URBAN_PCT"      , pct_urvid  ) )
   CALL nccheck( nf90_put_var  (ncid, pct_urvid        , pct_ur    ) )

   CALL nccheck( nf90_inq_varid(ncid, "TK_ROOF"        , tk_rfvid   ) )
   CALL nccheck( nf90_put_var  (ncid, tk_rfvid         , tk_rf     ) )

   CALL nccheck( nf90_inq_varid(ncid, "TK_WALL"        , tk_wlvid   ) )
   CALL nccheck( nf90_put_var  (ncid, tk_wlvid         , tk_wl     ) )

   CALL nccheck( nf90_inq_varid(ncid, "TK_IMPROAD"     , tk_imrdvid ) )
   CALL nccheck( nf90_put_var  (ncid, tk_imrdvid       , tk_imrd   ) )

   CALL nccheck( nf90_inq_varid(ncid, "CV_ROOF"        , cv_rfvid   ) )
   CALL nccheck( nf90_put_var  (ncid, cv_rfvid         , cv_rf     ) )

   CALL nccheck( nf90_inq_varid(ncid, "CV_WALL"        , cv_wlvid   ) )
   CALL nccheck( nf90_put_var  (ncid, cv_wlvid         , cv_wl     ) )

   CALL nccheck( nf90_inq_varid(ncid, "CV_IMPROAD"     , cv_imrdvid ) )
   CALL nccheck( nf90_put_var  (ncid, cv_imrdvid       , cv_imrd   ) )

   CALL nccheck( nf90_inq_varid(ncid, "ALB_ROOF"       , alb_rfvid  ) )
   CALL nccheck( nf90_put_var  (ncid, alb_rfvid        , alb_rf    ) )

   CALL nccheck( nf90_inq_varid(ncid, "ALB_WALL"       , alb_wlvid  ) )
   CALL nccheck( nf90_put_var  (ncid, alb_wlvid        , alb_wl    ) )

   CALL nccheck( nf90_inq_varid(ncid, "ALB_IMPROAD"    , alb_imrdvid) )
   CALL nccheck( nf90_put_var  (ncid, alb_imrdvid      , alb_imrd  ) )

   CALL nccheck( nf90_inq_varid(ncid, "ALB_PERROAD"    , alb_perdvid) )
   CALL nccheck( nf90_put_var  (ncid, alb_perdvid      , alb_perd  ) )

   CALL nccheck( nf90_close(ncid) )
#endif

   PRINT*, "*** SUCCESS write surface file ***"
   
   deallocate( hlat    )
   deallocate( hlats   )
   deallocate( hlatn   )
   deallocate( hlon    )
   deallocate( hlonw   )
   deallocate( hlone   )
   deallocate( harea   )
   deallocate( gfcc_tc )
   deallocate( gedi_th )
   deallocate( gl30_wt )
   deallocate( hlai    )
   deallocate( hsai    )

#ifndef USE_LCZ
   deallocate( urrgid  )
   deallocate( hwr_can   )
   deallocate( wt_rf     )
   deallocate( wt_rd     )
   deallocate( em_rf     )
   deallocate( em_wl     )
   deallocate( em_imrd   )
   deallocate( em_perd   )
   deallocate( ht_rf     )
   deallocate( th_rf     )
   deallocate( th_wl     )
   deallocate( tb_min    )
   deallocate( tb_max    )
   deallocate( cv_rf     )
   deallocate( cv_wl     )
   deallocate( cv_imrd   )
   deallocate( tk_rf     )
   deallocate( tk_wl     )
   deallocate( tk_imrd   )
   deallocate( alb_rf    )
   deallocate( alb_wl    )
   deallocate( alb_imrd  )
   deallocate( alb_perd  )
   ! ENDIF
#endif

   deallocate( latso     )
   deallocate( lonso     )
   deallocate( area      )
   deallocate( ur_dc     )
   deallocate( pct_ur    )

   deallocate( pct_tc    )
   deallocate( pct_urwt  )
   deallocate( htop_ur   )
   deallocate( ur_lai    )
   deallocate( ur_sai    )
END SUBROUTINE
