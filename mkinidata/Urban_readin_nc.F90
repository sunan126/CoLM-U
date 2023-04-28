#include <define.h>

SUBROUTINE Urban_readin_nc (dir_srfdata,dir_atmdata,nam_urbdata,nam_atmdata,lc_year)

! ===========================================================
! Read in the Urban dataset
! ===========================================================

      USE precision
      USE GlobalVars
      USE LC_Const
      USE MOD_TimeVariables
      USE MOD_TimeInvariants
      USE MOD_UrbanTimeInvars
      USE ncio
      USE omp_lib
      USE UrbanLCZ_Const

      IMPLICIT NONE

      INTEGER, intent(in) :: lc_year    ! which year of land cover data used
      CHARACTER(LEN=256), intent(in) :: dir_srfdata
      CHARACTER(LEN=256), intent(in) :: nam_urbdata
      CHARACTER(LEN=256), intent(in) :: dir_atmdata
      CHARACTER(LEN=256), intent(in) :: nam_atmdata

      CHARACTER(LEN=256) :: lndname
      CHARACTER(len=256) :: cyear

      INTEGER :: ncid
      INTEGER :: wtlunitroof_vid, htroof_vid, canyonhwr_vid, wtroadperv_vid, wtroadimperv_vid
      INTEGER :: urbanwaterpct_vid, urbantreepct_vid, urbantreetop_vid
      INTEGER :: albroof_vid, albwall_vid, albimproad_vid, albperroad_vid
      INTEGER :: emroof_vid, emwall_vid, emimproad_vid, emperroad_vid
      INTEGER :: cvroof_vid, cvwall_vid, cvimproad_vid
      INTEGER :: tkroof_vid, tkwall_vid, tkimproad_vid
      INTEGER :: tbuildingmax_vid, tbuildingmin_vid
      INTEGER :: thickroof_vid, thickwall_vid
      INTEGER :: reg_vid, pop_vid, veh_vid, wed_vid, weh_vid, wdh_vid, met_vid, hol_vid
      INTEGER :: urbanpop_vid, urbanlucy_vid
      INTEGER :: i, j, u, t, l, m, npatch, k, lucy_id

      REAL(r8) :: thick_roof, thick_wall

      ! define variables for reading in
      ! -------------------------------------------
      ! 城市形态结构参数
#ifdef USE_POINT_DATA
#ifdef USE_OBS_PARA
      REAL(r8) :: rfwt, rfht, tpct, wpct, hw_point, htop_point, prwt
#endif
#endif
      ! parameters for LUCY
      INTEGER , allocatable :: urbanlucy(:,:)! LUCY region id

      REAL(r8), allocatable :: urbanpop (:,:,:)! population density

      REAL(r8):: lweek_holiday(231,7)   , &  ! week holidays
                 lvehc_prof   (231,24,2), &  ! diurnal traffic profile
                 lhum_prof    (231,24 ) , &  ! diurnal metabolize profile
                 lfix_holiday (231,365) , &  ! Fixed public holidays, holiday(0) or workday(1)
                 lvehicle     (231,3)        ! vehicle numbers per thousand people

      ! morphological parameter
      REAL(r8), allocatable :: wtlunitroof   (:,:,:)
      REAL(r8), allocatable :: htroof        (:,:,:)
      REAL(r8), allocatable :: canyonhwr     (:,:,:)
      REAL(r8), allocatable :: wtroadperv    (:,:,:)
      REAL(r8), allocatable :: urbanwaterpct (:,:,:)
      REAL(r8), allocatable :: urbantreepct  (:,:,:)
      REAL(r8), allocatable :: urbantreetop  (:,:,:)
      

      ! albedo
      REAL(r8), allocatable :: albroof   (:,:,:,:,:)
      REAL(r8), allocatable :: albwall   (:,:,:,:,:)
      REAL(r8), allocatable :: albimproad(:,:,:,:,:)
      REAL(r8), allocatable :: albperroad(:,:,:,:,:)

      ! emissivity
      REAL(r8), allocatable :: emroof        (:,:,:)
      REAL(r8), allocatable :: emwall        (:,:,:)
      REAL(r8), allocatable :: emimproad     (:,:,:)
      REAL(r8), allocatable :: emperroad     (:,:,:)

      ! thermal pars of roof, wall, impervious
      REAL(r8), allocatable :: cvroof      (:,:,:,:)
      REAL(r8), allocatable :: cvwall      (:,:,:,:)
      REAL(r8), allocatable :: cvimproad   (:,:,:,:)

      ! thermal pars of roof, wall, impervious
      REAL(r8), allocatable :: tkroof      (:,:,:,:)
      REAL(r8), allocatable :: tkwall      (:,:,:,:)
      REAL(r8), allocatable :: tkimproad   (:,:,:,:)

      ! room maximum and minimum temperature
      REAL(r8), allocatable :: tbuildingmax  (:,:,:)
      REAL(r8), allocatable :: tbuildingmin  (:,:,:)

      ! thickness of roof and wall
      REAL(r8), allocatable :: thickroof     (:,:,:)
      REAL(r8), allocatable :: thickwall     (:,:,:)

      allocate ( wtlunitroof   (1:lon_points,1:lat_points,1:N_URB) )
      allocate ( htroof        (1:lon_points,1:lat_points,1:N_URB) )
      allocate ( canyonhwr     (1:lon_points,1:lat_points,1:N_URB) )
      allocate ( wtroadperv    (1:lon_points,1:lat_points,1:N_URB) )
      allocate ( urbanwaterpct (1:lon_points,1:lat_points,1:N_URB) )
      allocate ( urbantreepct  (1:lon_points,1:lat_points,1:N_URB) )
      allocate ( urbantreetop  (1:lon_points,1:lat_points,1:N_URB) )
      allocate ( albroof       (1:lon_points,1:lat_points,1:N_URB,2,2) )
      allocate ( albwall       (1:lon_points,1:lat_points,1:N_URB,2,2) )
      allocate ( albimproad    (1:lon_points,1:lat_points,1:N_URB,2,2) )
      allocate ( albperroad    (1:lon_points,1:lat_points,1:N_URB,2,2) )
      allocate ( emroof        (1:lon_points,1:lat_points,1:N_URB) )
      allocate ( emwall        (1:lon_points,1:lat_points,1:N_URB) )
      allocate ( emimproad     (1:lon_points,1:lat_points,1:N_URB) )
      allocate ( emperroad     (1:lon_points,1:lat_points,1:N_URB) )
      allocate ( cvroof        (1:lon_points,1:lat_points,1:N_URB,nl_roof) )
      allocate ( cvwall        (1:lon_points,1:lat_points,1:N_URB,nl_wall) )
      allocate ( cvimproad     (1:lon_points,1:lat_points,1:N_URB,nl_soil) )
      allocate ( tkroof        (1:lon_points,1:lat_points,1:N_URB,nl_roof) )
      allocate ( tkwall        (1:lon_points,1:lat_points,1:N_URB,nl_wall) )
      allocate ( tkimproad     (1:lon_points,1:lat_points,1:N_URB,nl_soil) )
      allocate ( tbuildingmax  (1:lon_points,1:lat_points,1:N_URB) )
      allocate ( tbuildingmin  (1:lon_points,1:lat_points,1:N_URB) )
      allocate ( thickroof     (1:lon_points,1:lat_points,1:N_URB) )
      allocate ( thickwall     (1:lon_points,1:lat_points,1:N_URB) )
      allocate ( urbanpop      (1:lon_points,1:lat_points,1:N_URB) )
      allocate ( urbanlucy     (1:lon_points,1:lat_points) )

#ifdef USE_LCZ

      write(cyear,'(i4.4)') lc_year
      lndname = trim(dir_srfdata)//trim(cyear)//'/'//trim(nam_urbdata)
      print*,trim(lndname)
      CALL nccheck( nf90_open(trim(lndname), nf90_nowrite, ncid) )

      CALL nccheck( nf90_inq_varid(ncid, "LCZ_LUCY_id"  , urbanlucy_vid ) )
      CALL nccheck( nf90_get_var  (ncid, urbanlucy_vid  , urbanlucy     ) )

      CALL nccheck( nf90_inq_varid(ncid, "LCZ_POP_DEN"  , urbanpop_vid  ) )
      CALL nccheck( nf90_get_var  (ncid, urbanpop_vid   , urbanpop      ) )

      CALL nccheck( nf90_inq_varid(ncid, "LCZ_WT_ROOF"  , wtlunitroof_vid) )
      CALL nccheck( nf90_get_var  (ncid, wtlunitroof_vid, wtlunitroof    ) )

      CALL nccheck( nf90_inq_varid(ncid, "LCZ_HT_ROOF"  , htroof_vid     ) )
      CALL nccheck( nf90_get_var  (ncid, htroof_vid     , htroof_vid     ) )

      CALL nccheck( nf90_close(ncid) )
      !TODO: change 360, 720 to parameters
      do i=1, 360
         do j=1,720
            canyonhwr   (j,i,:) = h2w     (:)
            wtroadperv  (j,i,:) = perfrac (:)
            emroof      (j,i,:) = roofem  (:)
            emwall      (j,i,:) = wallem  (:)
            emimproad   (j,i,:) = roadem  (:)
            emperroad   (j,i,:) = perem   (:)
            tbuildingmax(j,i,:) = 297.65  !TODO: check from WRF TBL, =TARGTEMP(294.15)+GAPTEMP(3.5)
            tbuildingmin(j,i,:) = 290.65  !TODO: check from WRF TBL, =TARGTEMP(294.15)-GAPTEMP(3.5)
            thickroof   (j,i,:) = rooftk  (:)
            thickwall   (j,i,:) = walltk  (:)

            albroof     (j,i,:,1,1) = roofalb(:)
            albroof     (j,i,:,1,2) = roofalb(:)
            albroof     (j,i,:,2,1) = roofalb(:)
            albroof     (j,i,:,2,2) = roofalb(:)
            albwall     (j,i,:,1,1) = wallalb(:)
            albwall     (j,i,:,1,2) = wallalb(:)
            albwall     (j,i,:,2,1) = wallalb(:)
            albwall     (j,i,:,2,2) = wallalb(:)
            albimproad  (j,i,:,1,1) = roadalb(:)
            albimproad  (j,i,:,1,2) = roadalb(:)
            albimproad  (j,i,:,2,1) = roadalb(:)
            albimproad  (j,i,:,2,2) = roadalb(:)
            albperroad  (j,i,:,1,1) = peralb (:)
            albperroad  (j,i,:,1,2) = peralb (:)
            albperroad  (j,i,:,2,1) = peralb (:)
            albperroad  (j,i,:,2,2) = peralb (:)

            do k=1,10
               cvroof   (j,i,k,:) = roofcv(k)
               cvwall   (j,i,k,:) = wallcv(k)
               cvimproad(j,i,k,:) = roadcv(k)
               tkroof   (j,i,k,:) = throof(k)
               tkwall   (j,i,k,:) = thwall(k)
               tkimproad(j,i,k,:) = throad(k)
            enddo
         enddo
      enddo
#else
! READ in urban data
      write(cyear,'(i4.4)') lc_year
      lndname = trim(dir_srfdata)//trim(cyear)//'/'//trim(nam_urbdata)
      print*,trim(lndname)
      CALL nccheck( nf90_open(trim(lndname), nf90_nowrite, ncid) )

      CALL nccheck( nf90_inq_varid(ncid, "URBAN_LUCY_id",   urbanlucy_vid    ) )
      CALL nccheck( nf90_inq_varid(ncid, "URBAN_POP_DEN",   urbanpop_vid     ) )
      CALL nccheck( nf90_inq_varid(ncid, "WTLUNIT_ROOF",    wtlunitroof_vid  ) )
      CALL nccheck( nf90_inq_varid(ncid, "HT_ROOF",         htroof_vid       ) )
      CALL nccheck( nf90_inq_varid(ncid, "CANYON_HWR",      canyonhwr_vid    ) )
      CALL nccheck( nf90_inq_varid(ncid, "WTROAD_PERV",     wtroadperv_vid   ) )
      CALL nccheck( nf90_inq_varid(ncid, "URBAN_WATER_PCT", urbanwaterpct_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "URBAN_TREE_PCT",  urbantreepct_vid ) )
      CALL nccheck( nf90_inq_varid(ncid, "URBAN_TREE_TOP",  urbantreetop_vid ) )
      CALL nccheck( nf90_inq_varid(ncid, "ALB_ROOF",        albroof_vid      ) )
      CALL nccheck( nf90_inq_varid(ncid, "ALB_WALL",        albwall_vid      ) )
      CALL nccheck( nf90_inq_varid(ncid, "ALB_IMPROAD",     albimproad_vid   ) )
      CALL nccheck( nf90_inq_varid(ncid, "ALB_PERROAD",     albperroad_vid   ) )
      CALL nccheck( nf90_inq_varid(ncid, "EM_ROOF",         emroof_vid       ) )
      CALL nccheck( nf90_inq_varid(ncid, "EM_WALL",         emwall_vid       ) )
      CALL nccheck( nf90_inq_varid(ncid, "EM_IMPROAD",      emimproad_vid    ) )
      CALL nccheck( nf90_inq_varid(ncid, "EM_PERROAD",      emperroad_vid    ) )
      CALL nccheck( nf90_inq_varid(ncid, "CV_ROOF",         cvroof_vid       ) )
      CALL nccheck( nf90_inq_varid(ncid, "CV_WALL",         cvwall_vid       ) )
      CALL nccheck( nf90_inq_varid(ncid, "CV_IMPROAD",      cvimproad_vid    ) )
      CALL nccheck( nf90_inq_varid(ncid, "TK_ROOF",         tkroof_vid       ) )
      CALL nccheck( nf90_inq_varid(ncid, "TK_WALL",         tkwall_vid       ) )
      CALL nccheck( nf90_inq_varid(ncid, "TK_IMPROAD",      tkimproad_vid    ) )
      CALL nccheck( nf90_inq_varid(ncid, "T_BUILDING_MAX",  tbuildingmax_vid ) )
      CALL nccheck( nf90_inq_varid(ncid, "T_BUILDING_MIN",  tbuildingmin_vid ) )
      CALL nccheck( nf90_inq_varid(ncid, "THICK_ROOF",      thickroof_vid    ) )
      CALL nccheck( nf90_inq_varid(ncid, "THICK_WALL",      thickwall_vid    ) )

      CALL nccheck( nf90_get_var(ncid, urbanlucy_vid  ,   urbanlucy    ) )
      CALL nccheck( nf90_get_var(ncid, urbanpop_vid   ,   urbanpop     ) )
      CALL nccheck( nf90_get_var(ncid, wtlunitroof_vid,   wtlunitroof  ) )
      CALL nccheck( nf90_get_var(ncid, htroof_vid,        htroof       ) )
      CALL nccheck( nf90_get_var(ncid, canyonhwr_vid,     canyonhwr    ) )
      CALL nccheck( nf90_get_var(ncid, wtroadperv_vid,    wtroadperv   ) )
      CALL nccheck( nf90_get_var(ncid, urbanpop_vid,      urbanpop     ) )
      CALL nccheck( nf90_get_var(ncid, urbanwaterpct_vid, urbanwaterpct) )
      CALL nccheck( nf90_get_var(ncid, urbantreepct_vid,  urbantreepct ) )
      CALL nccheck( nf90_get_var(ncid, urbantreetop_vid,  urbantreetop ) )
      CALL nccheck( nf90_get_var(ncid, albroof_vid,       albroof      ) )
      CALL nccheck( nf90_get_var(ncid, albwall_vid,       albwall      ) )
      CALL nccheck( nf90_get_var(ncid, albimproad_vid,    albimproad   ) )
      CALL nccheck( nf90_get_var(ncid, albperroad_vid,    albperroad   ) )
      CALL nccheck( nf90_get_var(ncid, emroof_vid,        emroof       ) )
      CALL nccheck( nf90_get_var(ncid, emwall_vid,        emwall       ) )
      CALL nccheck( nf90_get_var(ncid, emimproad_vid,     emimproad    ) )
      CALL nccheck( nf90_get_var(ncid, emperroad_vid,     emperroad    ) )
      CALL nccheck( nf90_get_var(ncid, cvroof_vid,        cvroof       ) )
      CALL nccheck( nf90_get_var(ncid, cvwall_vid,        cvwall       ) )
      CALL nccheck( nf90_get_var(ncid, cvimproad_vid,     cvimproad    ) )
      CALL nccheck( nf90_get_var(ncid, tkroof_vid,        tkroof       ) )
      CALL nccheck( nf90_get_var(ncid, tkwall_vid,        tkwall       ) )
      CALL nccheck( nf90_get_var(ncid, tkimproad_vid,     tkimproad    ) )
      CALL nccheck( nf90_get_var(ncid, tbuildingmax_vid,  tbuildingmax ) )
      CALL nccheck( nf90_get_var(ncid, tbuildingmin_vid,  tbuildingmin ) )
      CALL nccheck( nf90_get_var(ncid, thickroof_vid,     thickroof    ) )
      CALL nccheck( nf90_get_var(ncid, thickwall_vid,     thickwall    ) )

      CALL nccheck( nf90_close(ncid) )

      lndname = trim("/stu01/dongwz/data/CLMrawdata/urban_5x5/LUCY_rawdata.nc")
      CALL nccheck( nf90_open(trim(lndname), nf90_nowrite, ncid) )

      CALL nccheck( nf90_inq_varid(ncid, "vehicle"    , veh_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "weekendday" , wed_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "weekendhour", weh_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "weekdayhour", wdh_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "metabolism" , met_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "holiday"    , hol_vid) )

      CALL nccheck( nf90_get_var(ncid, veh_vid, lvehicle) )
      CALL nccheck( nf90_get_var(ncid, wed_vid, lweek_holiday) )
      CALL nccheck( nf90_get_var(ncid, weh_vid, lvehc_prof(:,:,2)) )
      CALL nccheck( nf90_get_var(ncid, wdh_vid, lvehc_prof(:,:,1)) )
      CALL nccheck( nf90_get_var(ncid, met_vid, lhum_prof) )
      CALL nccheck( nf90_get_var(ncid, hol_vid, lfix_holiday) )

      CALL nccheck( nf90_close(ncid) )
#ifdef USE_POINT_DATA
#ifdef USE_OBS_PARA

      lndname = trim(dir_atmdata)//'/'//trim(nam_atmdata)
      print*, lndname
      CALL nccheck( nf90_open(lndname, nf90_nowrite, ncid) )

      CALL nccheck( nf90_inq_varid(ncid, "impervious_area_fraction" , wtroadimperv_vid ) ) ! imperivous area fraciton
      CALL nccheck( nf90_inq_varid(ncid, "tree_area_fraction"       , urbantreepct_vid ) ) ! tree area fraction
      CALL nccheck( nf90_inq_varid(ncid, "water_area_fraction"      , urbanwaterpct_vid) ) ! water area fraction
      CALL nccheck( nf90_inq_varid(ncid, "roof_area_fraction"       , wtlunitroof_vid  ) ) ! roof area fraction
      CALL nccheck( nf90_inq_varid(ncid, "building_mean_height"     , htroof_vid       ) ) ! building mean height
      CALL nccheck( nf90_inq_varid(ncid, "tree_mean_height"         , urbantreetop_vid ) ) ! tree mean height
      CALL nccheck( nf90_inq_varid(ncid, "canyon_height_width_ratio", canyonhwr_vid    ) ) ! H2W

      CALL nccheck( nf90_get_var(ncid, wtroadimperv_vid,  prwt      ) )
      CALL nccheck( nf90_get_var(ncid, urbantreepct_vid,  tpct      ) )
      CALL nccheck( nf90_get_var(ncid, urbanwaterpct_vid, wpct      ) )
      CALL nccheck( nf90_get_var(ncid, wtlunitroof_vid,   rfwt      ) )
      CALL nccheck( nf90_get_var(ncid, htroof_vid,        rfht      ) )
      CALL nccheck( nf90_get_var(ncid, urbantreetop_vid,  htop_point) )
      CALL nccheck( nf90_get_var(ncid, canyonhwr_vid,     hw_point  ) )

      CALL nccheck( nf90_close(ncid) )

      wtroadperv   (1,1,:) = 1 - (prwt-rfwt)/(1-rfwt-wpct) !1. - prwt
      urbantreepct (1,1,:) = tpct*100
      urbanwaterpct(1,1,:) = wpct*100
      wtlunitroof  (1,1,:) = rfwt
      htroof       (1,1,:) = rfht
      urbantreetop (1,1,:) = htop_point
      canyonhwr    (1,1,:) = hw_point
      !albroof(:,:,:,:,:) = 0.4
#endif
#endif
#endif

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i, j, u, t, l, m, npatch, lucy_id) &
!$OMP PRIVATE(thick_roof, thick_wall)
#endif
      DO u = 1, numurban

         npatch = urb2patch(u)
         i = patch2lon(npatch)
         j = patch2lat(npatch)
         m = patchclass(npatch)
         t = urbclass(u)

         lucy_id         = urbanlucy     (i,j)
         IF (lucy_id > 0) THEN
            vehicle(u,:)      = lvehicle(lucy_id,:)
            week_holiday(u,:) = lweek_holiday(lucy_id,:)
            weh_prof(u,:)     = lvehc_prof(lucy_id,:,2)
            wdh_prof(u,:)     = lvehc_prof(lucy_id,:,1)
            hum_prof(u,:)     = lhum_prof(lucy_id,:)
            fix_holiday(u,:)  = lfix_holiday(lucy_id,:)
         ENDIF

         popcell(u)      = urbanpop      (i,j,t)
         froof(u)        = wtlunitroof   (i,j,t) !roof fractional cover
         hroof(u)        = htroof        (i,j,t) !average building height
         hwr(u)          = canyonhwr     (i,j,t) !average building height to their distance
         fgper(u)        = wtroadperv    (i,j,t) !pervious fraction to ground area

         alb_roof(:,:,u) = albroof   (i,j,t,:,:) !albedo of roof
         alb_wall(:,:,u) = albwall   (i,j,t,:,:) !albedo of walls
         alb_gimp(:,:,u) = albimproad(i,j,t,:,:) !albedo of impervious
         alb_gper(:,:,u) = albperroad(i,j,t,:,:) !albedo of pervious road

         em_roof(u)      = emroof        (i,j,t) !emissivity of roof
         em_wall(u)      = emwall        (i,j,t) !emissiviry of wall
         em_gimp(u)      = emimproad     (i,j,t) !emissivity of impervious
         em_gper(u)      = emperroad     (i,j,t) !emissivity of pervious

         cv_roof(:,u)    = cvroof      (i,j,t,:) !heat capacity of roof [J/(m2 K)]
         cv_wall(:,u)    = cvwall      (i,j,t,:) !heat capacity of wall [J/(m2 K)]
         cv_gimp(:,u)    = cvimproad   (i,j,t,:) !heat capacity of impervious [J/(m2 K)]

         tk_roof(:,u)    = tkroof      (i,j,t,:) !thermal conductivity of roof [W/m-K]
         tk_wall(:,u)    = tkwall      (i,j,t,:) !thermal conductivity of wall [W/m-K]
         tk_gimp(:,u)    = tkimproad   (i,j,t,:) !thermal conductivity of impervious [W/m-K]

         thick_roof      = thickroof     (i,j,t) !thickness of roof [m]
         thick_wall      = thickwall     (i,j,t) !thickness of wall [m]

#ifdef URBAN_BEM
         t_roommax(u)    = tbuildingmax  (i,j,t) !maximum temperature of inner room [K]
         t_roommin(u)    = tbuildingmin  (i,j,t) !minimum temperature of inner room [K]
#else
         t_roommax(u)    = 373.16                !maximum temperature of inner room [K]
         t_roommin(u)    = 180.00                !minimum temperature of inner room [K]
#endif

#ifdef URBAN_WATER
         flake(u) = urbanwaterpct(i,j,t)/100. !urban water fractional cover
#else
         flake(u) = 0.
#endif

#ifdef URBAN_TREE
         ! set tree fractional cover (<= 1.-froof)
         ! 植被覆盖占非水体面积部分的比例
         fveg(npatch) = urbantreepct(i,j,t)/100. !urban tree percent
         IF (flake(u) < 1.) THEN
            fveg(npatch) = fveg(npatch)/(1.-flake(u))
         ELSE
            fveg(npatch) = 0.
         ENDIF
         ! 假设树的覆盖比例小于等于地面比例(屋顶没有树)
         fveg(npatch) = min(fveg(npatch), 1.-froof(u))
#else
         fveg(npatch) = 0.
#endif

         ! set urban tree crown top and bottom [m]
         htop(npatch) = min(hroof(u), urbantreetop(i,j,t))
         htop(npatch) = max(2., htop(npatch))
         hbot(npatch) = htop(npatch)*hbot0(m)/htop0(m)
         hbot(npatch) = max(1., hbot(npatch))

         ! roof and wall layer depth
         DO l=1, nl_roof
            z_roof(l,u) = (l-0.5)*(thick_roof/nl_roof)
         ENDDO

         DO l=1, nl_wall
            z_wall(l,u) = (l-0.5)*(thick_wall/nl_wall)
         ENDDO

         dz_roof(1,u) = 0.5*(z_roof(1,u)+z_roof(2,u))
         DO l = 2, nl_roof-1
            dz_roof(l,u) = 0.5*(z_roof(l+1,u)-z_roof(l-1,u))
         ENDDO
         dz_roof(nl_roof,u) = z_roof(nl_roof,u)-z_roof(nl_roof-1,u)

         dz_wall(1,u) = 0.5*(z_wall(1,u)+z_wall(2,u))
         DO l = 2, nl_wall-1
            dz_wall(l,u) = 0.5*(z_wall(l+1,u)-z_wall(l-1,u))
         ENDDO
         dz_wall(nl_wall,u) = z_wall(nl_wall,u)-z_wall(nl_wall-1,u)

         ! lake depth and layer depth
         !lakedepth(npatch) = 1.
         !dz_lake(:,npatch) = lakedepth(npatch) / nl_lake

      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      deallocate ( wtlunitroof   )
      deallocate ( htroof        )
      deallocate ( canyonhwr     )
      deallocate ( wtroadperv    )
      deallocate ( urbanwaterpct )
      deallocate ( urbantreetop  )
      deallocate ( albroof       )
      deallocate ( albwall       )
      deallocate ( albimproad    )
      deallocate ( albperroad    )
      deallocate ( emroof        )
      deallocate ( emwall        )
      deallocate ( emimproad     )
      deallocate ( emperroad     )
      deallocate ( cvroof        )
      deallocate ( cvwall        )
      deallocate ( cvimproad     )
      deallocate ( tkroof        )
      deallocate ( tkwall        )
      deallocate ( tkimproad     )
      deallocate ( tbuildingmax  )
      deallocate ( tbuildingmin  )
      deallocate ( thickroof     )
      deallocate ( thickwall     )

!      CALL nccheck( nf90_close(ncid) )
END SUBROUTINE Urban_readin_nc
