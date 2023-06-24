#include <define.h>

SUBROUTINE LuLccEnergyMassConserve
! -------------------------------
! Created by Wanyi Lin and Hua Yuan, 06/2023
! -------------------------------

   USE precision
   USE GlobalVars
   USE MOD_TimeInvariants
   USE MOD_PFTimeInvars
   USE MOD_PCTimeInvars
   USE MOD_UrbanTimeInvars
   USE MOD_LuLccTimeInvars
   USE MOD_TimeVariables
   USE MOD_PFTimeVars
   USE MOD_PCTimeVars
   USE MOD_UrbanTimeVars
   USE MOD_LuLccTimeVars
   USE MOD_LuLccTransferMatrix
   IMPLICIT NONE


   INTEGER, parameter   :: nigbp = 17
   INTEGER, allocatable :: frnp_(:) !index of source patches

   INTEGER  :: k,ilc,num,itmpnp_
   INTEGER  :: i, j, np, np_, selfnp_
   REAL(r8) :: lccpcttmp(nigbp)
   REAL(r8) :: diff,wgt(1:nl_soil)

! #ifdef OPENMP
! print *, 'OPENMP enabled, threads num = ', OPENMP
! !$OMP PARALLEL DO NUM_THREADS(OPENMP) &
! !$OMP PRIVATE(i,j,np,np_)
! #endif

   DO j = 1, lat_points
      DO i = 1, lon_points
         np = grid_patch_s (i,j) !初始化后的编号
         np_= grid_patch_s_(i,j) !保存了上一年的编号

         IF (np.le.0) CYCLE  ! 不考虑海洋

         DO WHILE (np.le.grid_patch_e(i,j))
         ! 找出每个新的patch来源的去年的patch

#ifdef IGBP_CLASSIFICATION

            IF ( patchtype(np) .ne. 0 ) THEN
               ! 非土壤patch,使用restart值
               itmpnp_ = np_
               DO WHILE (itmpnp_ .le. grid_patch_e_(i,j))
                  IF (patchclass_(itmpnp_) .eq. patchclass(np)) THEN
                     wliq_soisno (:,np)   = wliq_soisno_(:,itmpnp_)
                     wice_soisno (:,np)   = wice_soisno_(:,itmpnp_)
                     t_soisno    (:,np)   = t_soisno_   (:,itmpnp_)
                     z_sno       (:,np)   = z_sno_      (:,itmpnp_)
                     dz_sno      (:,np)   = dz_sno_     (:,itmpnp_)
                     t_grnd      (np)     = t_grnd_     (itmpnp_)
                     tleaf       (np)     = tleaf_      (itmpnp_)
                     ldew        (np)     = ldew_       (itmpnp_)
                     sag         (np)     = sag_        (itmpnp_)
                     scv         (np)     = scv_        (itmpnp_)
                     snowdp      (np)     = snowdp_     (itmpnp_)
                     fveg        (np)     = fveg_       (itmpnp_)
                     fsno        (np)     = fsno_       (itmpnp_)
                     sigf        (np)     = sigf_       (itmpnp_)
                     green       (np)     = green_      (itmpnp_)
                     alb         (:,:,np) = alb_        (:,:,itmpnp_)
                     ssun        (:,:,np) = ssun_       (:,:,itmpnp_)
                     ssha        (:,:,np) = ssha_       (:,:,itmpnp_)
                     thermk      (np)     = thermk_     (itmpnp_)
                     extkb       (np)     = extkb_      (itmpnp_)
                     extkd       (np)     = extkd_      (itmpnp_)
                     zwt         (np)     = zwt_        (itmpnp_)
                     wa          (np)     = wa_         (itmpnp_)

                     t_lake       (:,np)= t_lake_      (:,itmpnp_)
                     lake_icefrac (:,np)= lake_icefrac_(:,itmpnp_)

                     EXIT
                  ENDIF
                  itmpnp_ = itmpnp_ + 1
               ENDDO
               np = np + 1
               CYCLE
            ENDIF


            ! ===============土壤patch==================
            ! 计算来源patch的数目,即转换百分比>0
            lccpcttmp(:) = lccpct(np,1:nigbp,patchclass(np))
            num = count(lccpcttmp .gt. 0)

            IF (num .gt. 0) THEN
               allocate ( frnp_ (num) )

               ! case1:存在与自己类型不同的来源patch
               IF ( (sum(lccpcttmp) - lccpcttmp(patchclass(np))) .gt. 0 ) THEN

                  ! 获取来源patch的索引 np_, 存储为frnp_
                  k = 0
                  DO ilc = 1, nigbp
                     IF (lccpcttmp(ilc) .gt. 0) THEN
                        k = k + 1
                        itmpnp_ = np_
                        DO WHILE (itmpnp_ .le. grid_patch_e_(i,j))

                           ! 如果存在与今年的patch类型相同的来源patch，存储为selfnp_
                           IF (patchclass_(itmpnp_) .eq. patchclass(np)) THEN
                              selfnp_ = itmpnp_
                           ENDIF

                           IF (patchclass_(itmpnp_) .eq. ilc) THEN
                              frnp_(k) = itmpnp_
                              EXIT
                           ENDIF
                           itmpnp_ = itmpnp_ + 1
                        ENDDO

                     ELSE
                       CYCLE
                     ENDIF
                  ENDDO

                  ! 初始化变量
                  wliq_soisno  (:,np)                = 0  !liquid water in layers [kg/m2]
                  wice_soisno  (:,np)                = 0  !ice lens in layers [kg/m2]
                  t_soisno     (maxsnl+1:nl_soil,np) = 0  !soil + snow layer te ！mperature [K]
                  z_sno        (:,np)                = 0  !node depth [m]
                  dz_sno       (:,np)                = 0  !interface depth [m]
                  ldew         (np)                  = 0  !depth of water on foliage [mm]
                  sag          (np)                  = 0  !non dimensional snow age [-]
                  scv          (np)                  = 0  !snow cover, water equivalent [mm]
                  snowdp       (np)                  = 0  !snow depth [meter]
                  fveg         (np)                  = 0  !fraction of vegetation cover
                  fsno         (np)                  = 0  !fraction of snow cover on ground
                  sigf         (np)                  = 0  !fraction of veg cover, excluding snow-covered veg [-]
                  green        (np)                  = 0  !leaf greenness
                  zwt          (np)                  = 0  !the depth to water table [m]
                  wa           (np)                  = 0  !water storage in aquifer [mm]
                  tleaf        (np)                  = 0  !leaf temperature [K]
                  alb          (:,:,np)              = 0  !averaged albedo [-]
                  ssun         (:,:,np)              = 0  !sunlit canopy absorption for solar radiation (0-1)
                  ssha         (:,:,np)              = 0  !shaded canopy absorption for solar radiation (0-1)
                  thermk       (np)                  = 0  !canopy gap fraction for tir radiation
                  extkb        (np)                  = 0  !(k, g(mu)/mu) direct solar extinction coefficient
                  extkd        (np)                  = 0  !diffuse and scattered diffuse PAR extinction      coefficient


                  ! 获取来源patch中最多的雪层
                  itmpnp_ = frnp_(1)
                  IF ( num .ge. 2 ) THEN
                     DO k = 2, num
                        IF ( sum(z_sno_(:,frnp_(k))) .gt. sum(z_sno_(:,frnp_(k-1))) ) THEN
                           itmpnp_ = frnp_(k)
                        ENDIF
                     ENDDO
                  ENDIF

                  ! 雪层变量调整
                  t_soisno (:0,np) = t_soisno_(:0,itmpnp_)
                  z_sno    (: ,np) = z_sno_   (: ,itmpnp_)
                  dz_sno   (: ,np) = dz_sno_  (: ,itmpnp_)


                  ! 计算温度调整权重
                  wgt(1:nl_soil) = 0
                  DO k = 1, num
                     wgt(1:nl_soil) = wgt(1:nl_soil) + cvsoil_(:,frnp_(k))*lccpcttmp(k)
                  ENDDO

                  ! 其他变量调整
                  DO k = 1, num
                     wliq_soisno (:,np) = wliq_soisno (:,np) + wliq_soisno_(:,frnp_(k))*lccpcttmp(k)/sum(lccpcttmp(1:num))
                     wice_soisno (:,np) = wice_soisno (:,np) + wice_soisno_(:,frnp_(k))*lccpcttmp(k)/sum(lccpcttmp(1:num))
                     t_soisno (1:nl_soil,np) = t_soisno (1:nl_soil,np) + t_soisno_(1:nl_soil,frnp_(k))*cvsoil_(:,frnp_(k))*lccpcttmp(k)/wgt
                     ldew  (np) = ldew  (np) + ldew_   (frnp_(k))*lccpcttmp(k)/sum(lccpcttmp(1:num))
                     sag   (np) = sag   (np) + sag_    (frnp_(k))*lccpcttmp(k)/sum(lccpcttmp(1:num))
                     scv   (np) = scv   (np) + scv_    (frnp_(k))*lccpcttmp(k)/sum(lccpcttmp(1:num))
                     snowdp(np) = snowdp(np) + snowdp_ (frnp_(k))*lccpcttmp(k)/sum(lccpcttmp(1:num))
                     fsno  (np) = fsno  (np) + fsno_   (frnp_(k))*lccpcttmp(k)/sum(lccpcttmp(1:num))
                     sigf  (np) = sigf  (np) + sigf_   (frnp_(k))*lccpcttmp(k)/sum(lccpcttmp(1:num))
                     zwt   (np) = zwt   (np) + zwt_    (frnp_(k))*lccpcttmp(k)/sum(lccpcttmp(1:num))
                     wa    (np) = wa    (np) + wa_     (frnp_(k))*lccpcttmp(k)/sum(lccpcttmp(1:num))

                     IF ( green_(frnp_(k)) .eq. 1 ) THEN
                     ! 任意一个来源patch的绿度是1，则设置为1
                        green (np) = green_(frnp_(k))
                     ELSE
                        print*, 'test: green_=',green_(frnp_(k)),'patchclass_=',patchclass_(frnp_(k))
                     ENDIF

                  ENDDO


                  ! 对以下变量，如果该patch在去年也存在，用restart赋值
                  IF (lccpcttmp(patchclass(np)) .gt. 0) THEN
                     tleaf    (np) = tleaf_    (selfnp_)
                     fveg     (np) = fveg_     (selfnp_)
                     alb  (:,:,np) = alb_  (:,:,selfnp_)
                     ssun (:,:,np) = ssun_ (:,:,selfnp_)
                     ssha (:,:,np) = ssha_ (:,:,selfnp_)
                     thermk   (np) = thermk_   (selfnp_)
                     extkb    (np) = extkb_    (selfnp_)
                     extkd    (np) = extkd_    (selfnp_)

                  ! 否则用面积加权平均
                  ELSE
                     print*, 'new patch, set area weighted average', 'i=',i, ',j=',j,';patchclass=',patchclass(np), &
                             ' ,patchfrac=',patchfrac(np)
                     DO k = 1, num
                        tleaf    (np) = tleaf    (np) + tleaf_    (frnp_(k)) *lccpcttmp(k)/sum(lccpcttmp(1:num))
                        fveg     (np) = fveg     (np) + fveg_     (frnp_(k)) *lccpcttmp(k)/sum(lccpcttmp(1:num))
                        alb  (:,:,np) = alb  (:,:,np) + alb_  (:,:,frnp_(k)) *lccpcttmp(k)/sum(lccpcttmp(1:num))
                        ssun (:,:,np) = ssun (:,:,np) + ssun_ (:,:,frnp_(k)) *lccpcttmp(k)/sum(lccpcttmp(1:num))
                        ssha (:,:,np) = ssha (:,:,np) + ssha_ (:,:,frnp_(k)) *lccpcttmp(k)/sum(lccpcttmp(1:num))
                        thermk   (np) = thermk   (np) + thermk_   (frnp_(k)) *lccpcttmp(k)/sum(lccpcttmp(1:num))
                        extkb    (np) = extkb    (np) + extkb_    (frnp_(k)) *lccpcttmp(k)/sum(lccpcttmp(1:num))
                        extkd    (np) = extkd    (np) + extkd_    (frnp_(k)) *lccpcttmp(k)/sum(lccpcttmp(1:num))
                     ENDDO

                  ENDIF

                  ! 地面温度赋值
                  IF ( sum( z_sno(:,np) ) .eq. 0 )  THEN
                  ! 雪层节点求和为0，则取土壤第一层的温度
                     t_grnd(np) = t_soisno(1,np)
                  ELSE
                  ! 存在雪层的话，取最上层的雪的温度
                     DO k = maxsnl+1, 0
                        IF ( z_sno(k,np) .le. 0 ) THEN
                           t_grnd(np) = t_soisno(k,np)
                           EXIT
                        ENDIF
                     ENDDO
                  ENDIF

                  deallocate ( frnp_ )

               ELSE
                  ! case2: 只有自己变成了自己(新一年的这种类型减少了或者占比没变)
                  itmpnp_ = np_
                  DO WHILE (itmpnp_ .le. grid_patch_e_(i,j))
                     IF (patchclass_(itmpnp_) .eq. patchclass(np)) THEN
                        wliq_soisno (:,np) = wliq_soisno_ (:,itmpnp_)
                        wice_soisno (:,np) = wice_soisno_ (:,itmpnp_)
                        t_soisno    (:,np) = t_soisno_    (:,itmpnp_)
                        z_sno       (:,np) = z_sno_       (:,itmpnp_)
                        dz_sno      (:,np) = dz_sno_      (:,itmpnp_)
                        t_grnd        (np) = t_grnd_        (itmpnp_)
                        tleaf         (np) = tleaf_         (itmpnp_)
                        ldew          (np) = ldew_          (itmpnp_)
                        sag           (np) = sag_           (itmpnp_)
                        scv           (np) = scv_           (itmpnp_)
                        snowdp        (np) = snowdp_        (itmpnp_)
                        fveg          (np) = fveg_          (itmpnp_)
                        fsno          (np) = fsno_          (itmpnp_)
                        sigf          (np) = sigf_          (itmpnp_)
                        green         (np) = green_         (itmpnp_)
                        alb       (:,:,np) = alb_       (:,:,itmpnp_)
                        ssun      (:,:,np) = ssun_      (:,:,itmpnp_)
                        ssha      (:,:,np) = ssha_      (:,:,itmpnp_)
                        thermk        (np) = thermk_        (itmpnp_)
                        extkb         (np) = extkb_         (itmpnp_)
                        extkd         (np) = extkd_         (itmpnp_)
                        zwt           (np) = zwt_           (itmpnp_)
                        wa            (np) = wa_            (itmpnp_)
                        EXIT
                     ENDIF

                     itmpnp_ = itmpnp_ + 1
                  ENDDO

                  deallocate ( frnp_ )
               ENDIF

            ELSE
               ! ! case3：如果这个patch的类型并不存在于新一年的地表数据中，用初始化的值。目前测试年份没有出现这种情况。
               ! print*, 'error: patchclass(np) not exist in new srfdata; use mkinitial value! stop'
               ! print*, 'i=',i, ',j=',j, ';patchclass=',patchclass(np), ',patchfrac=',patchfrac(np)
            ENDIF
#endif

            np = np + 1

         ENDDO

      ENDDO
   ENDDO

END SUBROUTINE LuLccEnergyMassConserve
! ---------- EOP ------------
