#include <define.h>

SUBROUTINE CLMDRIVER (idate,deltim,dolai,doalb,dosst,oro)


!=======================================================================
!
! CoLM MODEL DRIVER (with URBAN)
!
!=======================================================================

 USE precision
 USE PhysicalConstants, only: tfrz, rgas, vonkar
 USE GlobalVars
 USE LC_Const
 USE MOD_TimeInvariants
 USE MOD_TimeVariables
 USE MOD_UrbanTimeInvars
 USE MOD_UrbanTimeVars
 USE MOD_1D_Forcing
 USE MOD_1D_Fluxes
 USE omp_lib

 IMPLICIT NONE

  INTEGER,  intent(in) :: idate(3) !model calendar for next time step (year, julian day, seconds)
  REAL(r8), intent(in) :: deltim   !seconds in a time-step

  LOGICAL,  intent(in) :: dolai    !true if time for time-varying vegetation paramter
  LOGICAL,  intent(in) :: doalb    !true if time for surface albedo calculation
  LOGICAL,  intent(in) :: dosst    !true if time for update sst/ice/snow

  REAL(r8), intent(inout) :: oro(numpatch)  !ocean(0)/seaice(2)/ flag

  INTEGER :: i, m, u
  LOGICAL :: run_urban_model = .false.      !don't run urban model in default

! ======================================================================

#ifdef URBAN_MODEL
  run_urban_model = .true.
#endif

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i, m, u) &
!$OMP SCHEDULE(STATIC, 1)
#endif
  DO i = 1, numpatch

     m = patchclass(i)
     !print *, "patch:", i, "patchclass:", m

#ifdef URBAN_ONLY
     IF (m .ne. URBAN) CYCLE
#endif

     ! For slab urban or not urban patches
     IF (.not.run_urban_model .or. m.ne.URBAN) THEN

        CALL CLMMAIN (i, idate,           coszen(i),       deltim,          &
        patchlonr(i),    patchlatr(i),    patchclass(i),   patchtype(i),    &
        doalb,           dolai,           dosst,           oro(i),          &

      ! SOIL INFORMATION AND LAKE DEPTH
        soil_s_v_alb(i), soil_d_v_alb(i), soil_s_n_alb(i), soil_d_n_alb(i), &
        porsl(1:,i),     psi0(1:,i),      bsw(1:,i),       hksati(1:,i),    &
        csol(1:,i),      dksatu(1:,i),    dkdry(1:,i),     rootfr(1:,m),    &
        lakedepth(i),    dz_lake(1:,i),                                     &

      ! VEGETATION INFORMATION
        htop(i),         hbot(i),         sqrtdi(m),                        &
        effcon(m),       vmax25(m),       slti(m),         hlti(m),         &
        shti(m),         hhti(m),         trda(m),         trdm(m),         &
        trop(m),         gradm(m),        binter(m),       extkn(m),        &
        chil(m),         rho(1:,1:,m),    tau(1:,1:,m),                     &

      ! ATMOSPHERIC FORCING
        forc_pco2m(i),   forc_po2m(i),    forc_us(i),      forc_vs(i),      &
        forc_t(i),       forc_q(i),       forc_prc(i),     forc_prl(i),     &
        forc_rain(i),    forc_snow(i),    forc_psrf(i),    forc_pbot(i),    &
        forc_sols(i),    forc_soll(i),    forc_solsd(i),   forc_solld(i),   &
        forc_frl(i),     forc_hgt_u(i),   forc_hgt_t(i),   forc_hgt_q(i),   &
        forc_rhoair(i),                                                     &

      ! LAND SURFACE VARIABLES REQUIRED FOR RESTART
        z_sno(maxsnl+1:,i),               dz_sno(maxsnl+1:,i),              &
        t_soisno(maxsnl+1:,i),            wliq_soisno(maxsnl+1:,i),         &
        wice_soisno(maxsnl+1:,i),                                           &

        t_grnd(i),       tleaf(i),        ldew(i),         &
        sag(i),          scv(i),          snowdp(i),       fveg(i),         &
        fsno(i),         sigf(i),         green(i),        lai(i),          &
        sai(i),          alb(1:,1:,i),    ssun(1:,1:,i),   ssha(1:,1:,i),   &
        thermk(i),       extkb(i),        extkd(i),                         &

        zwt(i),          wa(i),                                             &
        t_lake(1:,i),    lake_icefrac(1:,i),                                &

      ! SNICAR snow model related
        snw_rds(:,i),    ssno(:,:,:,i),                                     &
        mss_bcpho(:,i),  mss_bcphi(:,i),  mss_ocpho(:,i),  mss_ocphi(:,i),  &
        mss_dst1(:,i),   mss_dst2(:,i),   mss_dst3(:,i),   mss_dst4(:,i),   &

      ! additional diagnostic variables for output
        laisun(i),       laisha(i),                                         &
        rstfac(i),       h2osoi(1:,i),    wat(i),                           &

      ! FLUXES
        taux(i),         tauy(i),         fsena(i),        fevpa(i),        &
        lfevpa(i),       fsenl(i),        fevpl(i),        etr(i),          &
        fseng(i),        fevpg(i),        olrg(i),         fgrnd(i),        &
        trad(i),         tref(i),         tmax(i),         tmin(i),         &
        qref(i),         rsur(i),         rnof(i),         qintr(i),        &
        qinfl(i),        qdrip(i),        rst(i),          assim(i),        &
        respc(i),        sabvsun(i),      sabvsha(i),      sabg(i),         &
        sr(i),           solvd(i),        solvi(i),        solnd(i),        &
        solni(i),        srvd(i),         srvi(i),         srnd(i),         &
        srni(i),         solvdln(i),      solviln(i),      solndln(i),      &
        solniln(i),      srvdln(i),       srviln(i),       srndln(i),       &
        srniln(i),       qcharge(i),      xerr(i),         zerr(i),         &

      ! TUNABLE modle constants
        zlnd,            zsno,            csoilc,          dewmx,           &
        wtfact,          capr,            cnfac,           ssi,             &
        wimp,            pondmx,          smpmax,          smpmin,          &
        trsmx0,          tcrit,                                             &

      ! additional variables required by coupling with WRF model
        emis(i),         z0m(i),          zol(i),          rib(i),          &
        ustar(i),        qstar(i),        tstar(i),                         &
        fm(i),           fh(i),           fq(i) )
     ENDIF

     ! For urban model and urban patches
     IF (run_urban_model .and. m.eq.URBAN) THEN

        u = patch2urb(i)
        !print *, "patch:", i, "urban:", u

        CALL UrbanCLMMAIN ( &
      ! MODEL RUNNING PARAMETERS
        i               ,idate           ,coszen(i)       ,deltim          ,&
        patchlonr(i)    ,patchlatr(i)    ,patchclass(i)   ,patchtype(i)    ,&

      ! URBAN PARAMETERS
        froof(u)        ,flake(u)        ,hroof(u)        ,hwr(u)          ,&
        fgper(u)        ,em_roof(u)      ,em_wall(u)      ,em_gimp(u)      ,&
        em_gper(u)      ,cv_roof(:,u)    ,cv_wall(:,u)    ,cv_gimp(:,u)    ,&
        tk_roof(:,u)    ,tk_wall(:,u)    ,tk_gimp(:,u)    ,z_roof(:,u)     ,&
        z_wall(:,u)     ,dz_roof(:,u)    ,dz_wall(:,u)                     ,&
        lakedepth(i)    ,dz_lake(1:,i)                                     ,&
      ! LUCY输入变量
        fix_holiday(u,:),week_holiday(u,:),hum_prof(u,:)  ,popcell(u)      ,&
        vehicle(u,:)    ,weh_prof(u,:)   ,wdh_prof(u,:)   ,Fahe(u)         ,&
      ! SOIL INFORMATION AND LAKE DEPTH
        porsl(1:,i)     ,psi0(1:,i)      ,bsw(1:,i)       ,hksati(1:,i)    ,&
        csol(1:,i)      ,dksatu(1:,i)    ,dkdry(1:,i)     ,rootfr(1:,m)    ,&
        alb_roof(:,:,u) ,alb_wall(:,:,u) ,alb_gimp(:,:,u) ,alb_gper(:,:,u) ,&

      ! VEGETATION INFORMATION
        htop(i)         ,hbot(i)         ,sqrtdi(m)       ,chil(m)         ,&
        effcon(m)       ,vmax25(m)       ,slti(m)         ,hlti(m)         ,&
        shti(m)         ,hhti(m)         ,trda(m)         ,trdm(m)         ,&
        trop(m)         ,gradm(m)        ,binter(m)       ,extkn(m)        ,&
        rho(1:,1:,m)    ,tau(1:,1:,m)                                      ,&

      ! ATMOSPHERIC FORCING
        forc_pco2m(i)   ,forc_po2m(i)    ,forc_us(i)      ,forc_vs(i)      ,&
        forc_t(i)       ,forc_q(i)       ,forc_prc(i)     ,forc_prl(i)     ,&
        forc_rain(i)    ,forc_snow(i)    ,forc_psrf(i)    ,forc_pbot(i)    ,&
        forc_sols(i)    ,forc_soll(i)    ,forc_solsd(i)   ,forc_solld(i)   ,&
        forc_frl(i)     ,forc_hgt_u(i)   ,forc_hgt_t(i)   ,forc_hgt_q(i)   ,&
        forc_rhoair(i)  ,Fhac(u)         ,Fwst(u)         ,Fach(u)         ,&

      ! LAND SURFACE VARIABLES REQUIRED FOR RESTART
        z_sno_roof  (maxsnl+1:,u)        ,z_sno_gimp  (maxsnl+1:,u)        ,&
        z_sno_gper  (maxsnl+1:,u)        ,z_sno_lake  (maxsnl+1:,u)        ,&
        dz_sno_roof (maxsnl+1:,u)        ,dz_sno_gimp (maxsnl+1:,u)        ,&
        dz_sno_gper (maxsnl+1:,u)        ,dz_sno_lake (maxsnl+1:,u)        ,&
        t_roofsno   (maxsnl+1:,u)        ,t_gimpsno   (maxsnl+1:,u)        ,&
        t_gpersno   (maxsnl+1:,u)        ,t_lakesno   (maxsnl+1:,u)        ,&
        wliq_roofsno(maxsnl+1:,u)        ,wliq_gimpsno(maxsnl+1:,u)        ,&
        wliq_gpersno(maxsnl+1:,u)        ,wliq_lakesno(maxsnl+1:,u)        ,&
        wice_roofsno(maxsnl+1:,u)        ,wice_gimpsno(maxsnl+1:,u)        ,&
        wice_gpersno(maxsnl+1:,u)        ,wice_lakesno(maxsnl+1:,u)        ,&
        z_sno       (maxsnl+1:,i)        ,dz_sno      (maxsnl+1:,i)        ,&
        wliq_soisno (maxsnl+1:,i)        ,wice_soisno (maxsnl+1:,i)        ,&
        t_soisno    (maxsnl+1:,i)        ,t_wallsun   (1:,u)               ,&
        t_wallsha   (1:,u)                                                 ,&

        lai(i)          ,sai(i)          ,fveg(i)         ,sigf(i)         ,&
        green(i)        ,tleaf(i)        ,ldew(i)         ,t_grnd(i)       ,&

        sag_roof(u)     ,sag_gimp(u)     ,sag_gper(u)     ,sag_lake(u)     ,&
        scv_roof(u)     ,scv_gimp(u)     ,scv_gper(u)     ,scv_lake(u)     ,&
        snowdp_roof(u)  ,snowdp_gimp(u)  ,snowdp_gper(u)  ,snowdp_lake(u)  ,&
        fsno_roof(u)    ,fsno_gimp(u)    ,fsno_gper(u)    ,fsno_lake(u)    ,&
        sag(i)          ,scv(i)          ,snowdp(i)       ,fsno(i)         ,&
        alb(1:,1:,i)    ,ssun(1:,1:,i)   ,ssha(1:,1:,i)   ,sroof(1:,1:,u)  ,&
        swsun(1:,1:,u)  ,swsha(1:,1:,u)  ,sgimp(1:,1:,u)  ,sgper(1:,1:,u)  ,&
        slake(1:,1:,u)  ,lwsun(u)        ,lwsha(u)        ,lgimp(u)        ,&
        lgper(u)        ,lveg(u)         ,fwsun(u)        ,dfwsun(u)       ,&
        t_room(u)       ,troof_inner(u)  ,twsun_inner(u)  ,twsha_inner(u)  ,&
        t_roommax(u)    ,t_roommin(u)    ,tafu(u)                          ,&

        zwt(i)          ,wa(i)                                             ,&
        t_lake(1:,i)    ,lake_icefrac(1:,i)                                ,&

      ! additional diagnostic variables for output
        laisun(i)       ,laisha(i)                                         ,&
        rstfac(i)       ,h2osoi(1:,i)    ,wat(i)                           ,&

      ! FLUXES
        taux(i)         ,tauy(i)         ,fsena(i)        ,fevpa(i)        ,&
        lfevpa(i)       ,fsenl(i)        ,fevpl(i)        ,etr(i)          ,&
        fseng(i)        ,fevpg(i)        ,olrg(i)         ,fgrnd(i)        ,&
        trad(i)         ,tref(i)         ,tmax(i)         ,tmin(i)         ,&
        qref(i)         ,rsur(i)         ,rnof(i)         ,qintr(i)        ,&
        qinfl(i)        ,qdrip(i)        ,rst(i)          ,assim(i)        ,&
        respc(i)        ,sabvsun(i)      ,sabvsha(i)      ,sabg(i)         ,&
        sr(i)           ,solvd(i)        ,solvi(i)        ,solnd(i)        ,&
        solni(i)        ,srvd(i)         ,srvi(i)         ,srnd(i)         ,&
        srni(i)         ,solvdln(i)      ,solviln(i)      ,solndln(i)      ,&
        solniln(i)      ,srvdln(i)       ,srviln(i)       ,srndln(i)       ,&
        srniln(i)       ,qcharge(i)      ,xerr(i)         ,zerr(i)         ,&

      ! TUNABLE modle constants
        zlnd            ,zsno            ,csoilc          ,dewmx           ,&
        wtfact          ,capr            ,cnfac           ,ssi             ,&
        wimp            ,pondmx          ,smpmax          ,smpmin          ,&
        trsmx0          ,tcrit                                             ,&

      ! additional variables required by coupling with WRF model
        emis(i)         ,z0m(i)          ,zol(i)          ,rib(i)          ,&
        ustar(i)        ,qstar(i)        ,tstar(i)        ,fm(i)           ,&
        fh(i)           ,fq(i)                                              )
     ENDIF

  ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

END SUBROUTINE CLMDRIVER
! ----------------------------------------------------------------------
! EOP
