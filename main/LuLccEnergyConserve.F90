#include <define.h>

  SUBROUTINE LuLccEnergyConserve
! -------------------------------
! Created by Hua Yuan, 04/2022
! May be renamed as "LuLccEnergyMassConserve"
! -------------------------------

    USE precision
    USE GlobalVars
    USE MOD_TimeInvariants
    USE MOD_PFTimeInvars
    USE MOD_PCTimeInvars
    USE MOD_LuLccTimeInvars
    USE MOD_TimeVariables
    USE MOD_PFTimeVars
    USE MOD_PCTimeVars
    USE MOD_LuLccTimeVars
    USE MOD_LuLccTMatrix


    IMPLICIT NONE
!TODO: need coding below...

    INTEGER, allocatable :: tmp(:) !某个patch的lc类型在to序列中的位置
    INTEGER, allocatable :: tmpfr(:) !某个patch的lc来源
    REAL(r8), allocatable :: lccpct(:) !某个patch的lc来源百分比，包括自己
    INTEGER, allocatable :: tmpnp_(:) !某个patch的来源patch的索引
    INTEGER :: fr(nnum) !lc旧序列
    INTEGER :: to(nnum) !lc新序列
    REAL(r8) :: pctt(nnum) !lcc百分比
    REAL(r8) :: pctself(nigbp) !保持自我的百分比
    REAL(r8) :: diff,wgt(1:nl_soil)
    INTEGER :: k,num,itmpnp_

    INTEGER i, j, np, np_

! #ifdef OPENMP
! print *, 'OPENMP enabled, threads num = ', OPENMP
! !$OMP PARALLEL DO NUM_THREADS(OPENMP) &
! !$OMP PRIVATE(i,j,np,np_)
! #endif

    DO j = 1, lat_points 
      DO i = 1, lon_points
        np = grid_patch_s (i,j) !初始化后的编号
        np_= grid_patch_s_(i,j) !保存了上一年的编号

        IF (np.le.0 .or. np_.le.0) CYCLE

        DO WHILE (np.le.grid_patch_e(i,j)) 
          ! 找出每个新的patch来源的去年的patch

#ifdef IGBP_CLASSIFICATION
          fr = lcfr(:,np)
          to = lcto(:,np)
          pctt = dchg(:,np)
          pctself = selfchg(:,np)
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

          ! 土壤patch；计数新一年ptach的来源patch
          num = count(to .eq. patchclass(np))
          IF (pctself(patchclass(np)) .gt. 0) THEN
            num = num + 1
          ENDIF

          IF (num .gt. 0) THEN           
            allocate ( tmp    (num) )
            allocate ( tmpfr  (num) )
            allocate ( lccpct (num) )
            allocate ( tmpnp_ (num) )

            ! 统计: to序列中lc==np的个数,并获取索引号 ->继而获取from序列的lc、LCC百分比
            IF ( count(to .eq. patchclass(np)) .gt. 0) THEN
              ! case1:来源patch中存在patchclass类型与自己不同的情况
              IF (pctself(patchclass(np)) .gt. 0) THEN           
                itmpnp_ = np_
                DO WHILE (itmpnp_ .le. grid_patch_e_(i,j))
                  IF (patchclass_(itmpnp_) .eq. patchclass(np)) THEN                 
                    tmpnp_(num) = itmpnp_                  
                    EXIT
                  ENDIF
                  itmpnp_ = itmpnp_ + 1
                ENDDO
                lccpct(num) = pctself(patchclass(np))
                num = num - 1
              ENDIF

              tmp(1:1) = findloc(to, patchclass(np))
              tmpfr(1) = fr(tmp(1))
              lccpct(1) = pctt(tmp(1))
                            
              DO k = 1, num
                
                IF (k .ge. 2) THEN
                  tmp(k:k) = findloc(to(tmp(k-1)+1:nnum),patchclass(np)) + tmp(k-1)
                  tmpfr(k) = fr(tmp(k))
                  lccpct(k) = pctt(tmp(k))
                ENDIF
                
                ! TO DO:检查有没有 在np_中找不到fr的情况
                itmpnp_ = np_
                DO WHILE (itmpnp_ .le. grid_patch_e_(i,j))
                  ! 根据 tmpfr获取该类型的索引，然后存到一个动态数组
                  IF (patchclass_(itmpnp_) .eq. tmpfr(k)) THEN
                    tmpnp_(k) = itmpnp_
                    EXIT
                  ENDIF
                  itmpnp_ = itmpnp_ + 1
                ENDDO
              ENDDO
              
              IF (pctself(patchclass(np)) .gt. 0) THEN 
                num = num + 1
              ENDIF

              ! ! CHECK===================
              ! print*, wliq_soisno(:,np)
              ! print*, wliq_soisno(1:10,np)
              ! print*, 'cv test:', cvsoil_(1,np)
              ! IF (sum(z_sno_(:,np_)) < 0) THEN
              !   print*,'np=',np,'lat=',gridlatd(j),'lon=',gridlond(i)
              !   print*, 'old layer z_sno_',z_sno_(:,np_)             
              !   print*, 'new layer z_sno',z_sno(:,np)             
              !   print*, 'old layer dz_sno_',dz_sno_(:,np_)             
              !   print*, 'new layer dz_sno',dz_sno(:,np)                 
              ! ENDIF
              ! ! CHECK===================

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
              extkd        (np)                  = 0  !diffuse and scattered diffuse PAR extinction coefficient

              wgt(1:nl_soil) = 0   ! 温度调整权重
              itmpnp_ = tmpnp_(1)  ! 获取来源patch中最多的雪层
              DO k = 1, num
                wgt(1:nl_soil) = wgt(1:nl_soil) + cvsoil_(:,tmpnp_(k))*lccpct(k)
                IF ( (k .ge. 2) .AND. (num .ge. 2) ) THEN
                  IF ( sum(z_sno_(:,tmpnp_(k))) .lt. sum(z_sno_(:,tmpnp_(k-1))) ) THEN
                    itmpnp_ = tmpnp_(k)
                  ENDIF
                ENDIF
              ENDDO


              DO k = 1, num
                wliq_soisno (:,np) = wliq_soisno (:,np) + wliq_soisno_(:,tmpnp_(k))*lccpct(k)/sum(lccpct(1:num))
                wice_soisno (:,np) = wice_soisno (:,np) + wice_soisno_(:,tmpnp_(k))*lccpct(k)/sum(lccpct(1:num))
                t_soisno (1:nl_soil,np) = t_soisno (1:nl_soil,np) + t_soisno_(1:nl_soil,tmpnp_(k))*cvsoil_(:,tmpnp_(k))*lccpct(k)/wgt
                ldew  (np) = ldew  (np) + ldew_   (tmpnp_(k))*lccpct(k)/sum(lccpct(1:num))
                sag   (np) = sag   (np) + sag_    (tmpnp_(k))*lccpct(k)/sum(lccpct(1:num))
                scv   (np) = scv   (np) + scv_    (tmpnp_(k))*lccpct(k)/sum(lccpct(1:num))
                snowdp(np) = snowdp(np) + snowdp_ (tmpnp_(k))*lccpct(k)/sum(lccpct(1:num))
                fsno  (np) = fsno  (np) + fsno_   (tmpnp_(k))*lccpct(k)/sum(lccpct(1:num))
                sigf  (np) = sigf  (np) + sigf_   (tmpnp_(k))*lccpct(k)/sum(lccpct(1:num))
                zwt   (np) = zwt   (np) + zwt_    (tmpnp_(k))*lccpct(k)/sum(lccpct(1:num))
                wa    (np) = wa    (np) + wa_     (tmpnp_(k))*lccpct(k)/sum(lccpct(1:num))

                IF ( green_(tmpnp_(k)).NE. 1 ) THEN
                  print*, 'test: green=',green(np),'patchlass=',patchclass(np)
                ELSE
                  green (np) = green_(tmpnp_(k))
                ENDIF

              ENDDO

              ! 对以下变量，如果该patch在去年也存在，用去年的节点和温度赋值；
              IF (pctself(patchclass(np)) .gt. 0) THEN 
                t_soisno(:0,np)  = t_soisno_(:0,tmpnp_(num))
                z_sno   (:,np)   = z_sno_   (:, tmpnp_(num))
                dz_sno  (:,np)   = dz_sno_  (:, tmpnp_(num))
                tleaf   (np)     = tleaf_   (tmpnp_(num))
                fveg    (np)     = fveg_    (tmpnp_(num))
                alb     (:,:,np) = alb_     (:,:,tmpnp_(num))
                ssun    (:,:,np) = ssun_    (:,:,tmpnp_(num))
                ssha    (:,:,np) = ssha_    (:,:,tmpnp_(num))
                thermk  (np)     = thermk_  (tmpnp_(num))
                extkb   (np)     = extkb_   (tmpnp_(num))
                extkd   (np)     = extkd_   (tmpnp_(num))

                IF (.NOT. (coszen(np) == coszen_(tmpnp_(num))) ) THEN
                  print*, 'error: restart coszen should == mkini coszen'
                ENDIF

              ! 否则雪的层数和温度用来源patch中最多的层赋值，其他变量用加权平均
              ELSE                
                ! print*, 'new patch, set largest available snow layers, set area weighted tleaf', 'i=',i, ',j=',j,';patchclass=',patchclass(np), ',patchfrac=',patchfrac(np)
                t_soisno (:0,np) = t_soisno_(:0,itmpnp_)
                z_sno (:,np) = z_sno_ (:,itmpnp_)
                dz_sno(:,np) = dz_sno_(:,itmpnp_)
                DO k = 1, num
                  tleaf (np)     = tleaf (np)     + tleaf_  (tmpnp_(k))     *lccpct(k)/sum(lccpct(1:num))
                  fveg  (np)     = fveg  (np)     + fveg_   (tmpnp_(k))     *lccpct(k)/sum(lccpct(1:num))
                  alb   (:,:,np) = alb   (:,:,np) + alb_    (:,:,tmpnp_(k)) *lccpct(k)/sum(lccpct(1:num))
                  ssun  (:,:,np) = ssun  (:,:,np) + ssun_   (:,:,tmpnp_(k)) *lccpct(k)/sum(lccpct(1:num))
                  ssha  (:,:,np) = ssha  (:,:,np) + ssha_   (:,:,tmpnp_(k)) *lccpct(k)/sum(lccpct(1:num))
                  thermk(np)     = thermk(np)     + thermk_ (tmpnp_(k))     *lccpct(k)/sum(lccpct(1:num))
                  extkb (np)     = extkb (np)     + extkb_  (tmpnp_(k))     *lccpct(k)/sum(lccpct(1:num))
                  extkd (np)     = extkd (np)     + extkd_  (tmpnp_(k))     *lccpct(k)/sum(lccpct(1:num))
                ENDDO

              ENDIF

              ! 地面温度赋值
              IF ( sum( z_sno(:,np) ) .eq. 0 )  THEN
                t_grnd(np) = t_soisno(1,np)
              ELSE
                DO k = maxsnl+1, 0
                  IF ( z_sno(k,np) .le. 0 ) THEN
                    t_grnd(np) = t_soisno(k,np)
                    EXIT
                  ENDIF
                ENDDO
              ENDIF


              ! CHECK==========================
              ! IF (pctself(patchclass(np)) .gt. 0) THEN 
              !   diff = abs(wliq_soisno (1,np) - wliq_soisno_(1,tmpnp_(num)))
              !   IF ( (diff .gt. wliq_soisno_(1,tmpnp_(num))*0.5) .and. (wliq_soisno_(1,tmpnp_(num)) .gt. 0) ) THEN
              !     print*, 'np=',np, ', j=',j, ', i=',i
              !     print*, 'diff = ',diff, 'new=',wliq_soisno (1,np), 'old=', wliq_soisno_(1,tmpnp_(num))
              !   ENDIF
              ! ENDIF

              IF (np .eq. 290858) THEN 
                print*, 'np=',np, ', j=',j, ', i=',i,'lat=',gridlatd(j),'lon=',gridlond(i)
                print*, 'patchclass=',patchclass(np), ',patchfrac=',patchfrac(np),'tmpfr=',tmpfr
                print*, 'new wliq_soisno=', wliq_soisno (:,np)  
                print*, 'old =           ', wliq_soisno_(:,tmpnp_(num))
                print*, 'new wice_soisno=', wice_soisno (:,np)  
                print*, 'old =           ', wice_soisno_(:,tmpnp_(num))
                print*, 'new t_soisno   =', t_soisno    (:,np)  
                print*, 'old =           ', t_soisno_   (:,tmpnp_(num))
                print*, 'new z_sno      =', z_sno       (:,np)  
                print*, 'old =           ', z_sno_      (:,tmpnp_(num))
                print*, 'new dz_sno     =', dz_sno      (:,np)  
                print*, 'old =           ', dz_sno_     (:,tmpnp_(num))  
                print*, 'new t_grnd     =', t_grnd      (np)    
                print*, 'old =           ', t_grnd_     (tmpnp_(num))
                print*, 'new tleaf      =', tleaf       (np)    
                print*, 'old =           ', tleaf_      (tmpnp_(num))
                print*, 'new ldew       =', ldew        (np)    
                print*, 'old =           ', ldew_       (tmpnp_(num))
                print*, 'new sag        =', sag         (np)    
                print*, 'old =           ', sag_        (tmpnp_(num))
                print*, 'new scv        =', scv         (np)    
                print*, 'old =           ', scv_        (tmpnp_(num))  
                print*, 'new snowdp     =', snowdp      (np)    
                print*, 'old =           ', snowdp_     (tmpnp_(num))
                print*, 'new fveg       =', fveg        (np)    
                print*, 'old =           ', fveg_       (tmpnp_(num))
                print*, 'new fsno       =', fsno        (np)    
                print*, 'old =           ', fsno_       (tmpnp_(num))
                print*, 'new sigf       =', sigf        (np)    
                print*, 'old =           ', sigf_       (tmpnp_(num))
                print*, 'new green      =', green       (np)    
                print*, 'old =           ', green_      (tmpnp_(num))  
                print*, 'new alb        =', alb         (:,:,np)
                print*, 'old =           ', alb_        (:,:,tmpnp_(num))
                print*, 'new ssun       =', ssun        (:,:,np)
                print*, 'old =           ', ssun_       (:,:,tmpnp_(num))
                print*, 'new ssha       =', ssha        (:,:,np)
                print*, 'old =           ', ssha_       (:,:,tmpnp_(num))
                print*, 'new thermk     =', thermk      (np)    
                print*, 'old =           ', thermk_     (tmpnp_(num))
                print*, 'new extkb      =', extkb       (np)    
                print*, 'old =           ', extkb_      (tmpnp_(num))  
                print*, 'new extkd      =', extkd       (np)    
                print*, 'old =           ', extkd_      (tmpnp_(num))
                print*, 'new zwt        =', zwt         (np)    
                print*, 'old =           ', zwt_        (tmpnp_(num))
                print*, 'new wa         =', wa          (np)    
                print*, 'old =           ', wa_         (tmpnp_(num))              
              ENDIF

              ! CHECK===========================

              deallocate ( tmp    )
              deallocate ( tmpfr  )
              deallocate ( lccpct )
              deallocate ( tmpnp_ )

            ELSE

              ! case2: 只有自己变成了自己(新一年的这种类型减少了或者占比没变)
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

                  ! print*, 'unchanged wliq_soisno'
                  EXIT
                ENDIF

                ! CHECK===========================
                ! IF (num .ge. 1) THEN
                ! IF (tmpnp_(num) == 0) THEN
                !   print*, tmpnp_
                !   print*, num
                ! ENDIF
                IF (.NOT. (coszen(np) == coszen_(itmpnp_)) ) THEN
                  print*, 'error: restart coszen should == mkini coszen'
                ENDIF
                ! ENDIF
                ! CHECK===========================

                itmpnp_ = itmpnp_ + 1
              ENDDO

              deallocate ( tmp    )
              deallocate ( tmpfr  )
              deallocate ( lccpct )
              deallocate ( tmpnp_ )
            ENDIF

          ELSE
            ! case3：原始地表数据中并没有这种patch，可能是mksrf过程中凭空产生的；占比都很少。
            ! print*, 'error: patchclass(np) not exist in list <to>; use mkinitial value' !stop
            ! print*, 'i=',i, ',j=',j, ';patchclass=',patchclass(np), ',patchfrac=',patchfrac(np)
          ENDIF
#endif

          np = np + 1
        
        ENDDO

        ! ! CHECK===================
        ! print*,'np=',grid_patch_s(i,j),'to',grid_patch_e(i,j)
        ! print*,'old patchclass_=',patchclass_(grid_patch_s_(i,j):grid_patch_e_(i,j)),&
        ! ',wliq_=',wliq_soisno_(1,grid_patch_s_(i,j):grid_patch_e_(i,j)),&
        ! ',wice_=',wice_soisno_(1,grid_patch_s_(i,j):grid_patch_e_(i,j))
        ! print*,'new patchclass=',patchclass(grid_patch_s(i,j):grid_patch_e(i,j)),&
        ! ',wliq=',wliq_soisno(1,grid_patch_s(i,j):grid_patch_e(i,j)),&
        ! ',wice=',wice_soisno(1,grid_patch_s(i,j):grid_patch_e(i,j))
        ! ! CHECK===================

      ENDDO
    ENDDO        

  END SUBROUTINE LuLccEnergyConserve
! ---------- EOP ------------
