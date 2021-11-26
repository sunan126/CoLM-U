#include <define.h>

MODULE UrbanFlux

!-----------------------------------------------------------------------
  USE precision
  USE GlobalVars
  IMPLICIT NONE
  SAVE

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: UrbanOnlyFlux
  PUBLIC :: UrbanVegFlux
  PUBLIC :: dewfraction

! PRIVATE MEMBER FUNCTIONS:
  PRIVATE :: cal_z0_displa

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE UrbanOnlyFlux ( &
        ! 模型运行信息
        ipatch      ,deltim                                ,&
        ! 外强迫
        hu          ,ht          ,hq          ,us          ,&
        vs          ,thm         ,th          ,thv         ,&
        qm          ,psrf        ,rhoair      ,Fhac        ,&
        Fwst        ,Fach                                  ,&
        ! 城市参数
        hroof       ,hlr         ,nurb        ,pondmx      ,&
        fcover                                             ,&
        ! 地面状态
        z0h_g       ,obug        ,ustarg      ,zlnd        ,&
        zsno        ,fsno_roof   ,fsno_gimp   ,fsno_gper   ,&
        htvp_roof   ,htvp_gimp   ,htvp_gper   ,troof       ,&
        twsun       ,twsha       ,tgimp       ,tgper       ,&
        qroof       ,qgimp       ,qgper       ,dqroofdT    ,&
        dqgimpdT    ,dqgperdT                              ,&
        ! 输出
        taux        ,tauy        ,fsenroof    ,fsenwsun    ,&
        fsenwsha    ,fsengimp    ,fsengper    ,fevproof    ,&
        fevpgimp    ,fevpgper    ,croofs      ,cwalls      ,&
        cgrnds      ,croofl      ,cgimpl      ,cgperl      ,&
        croof       ,cgimp       ,cgper       ,tref        ,&
        qref        ,z0m         ,zol         ,rib         ,&
        ustar       ,qstar       ,tstar       ,fm          ,&
        fh          ,fq          ,tafu                      )

!=======================================================================
     USE precision
     USE PhysicalConstants, only: cpair,vonkar,grav
     USE FRICTION_VELOCITY
     IMPLICIT NONE

!----------------------- Dummy argument --------------------------------
     INTEGER, intent(in) :: &
        ipatch     ! patch index [-]

     REAL(r8), intent(in) :: &
        deltim     ! seconds in a time step [second]

     ! atmospherical variables and observational height
     REAL(r8), intent(inout) :: &
        hu,       &! observational height of wind [m]
        ht,       &! observational height of temperature [m]
        hq         ! observational height of humidity [m]

     REAL(r8), intent(in) :: &
        us,       &! wind component in eastward direction [m/s]
        vs,       &! wind component in northward direction [m/s]
        thm,      &! intermediate variable (tm+0.0098*ht) [K]
        th,       &! potential temperature (kelvin)
        thv,      &! virtual potential temperature (kelvin)
        qm,       &! specific humidity at agcm reference height [kg/kg]
        psrf,     &! atmosphere pressure at the surface [pa] [not used]
        rhoair     ! density air [kg/m3]

     REAL(r8), intent(in) :: &
        Fhac,     &! flux from heat or cool AC
        Fwst,     &! waste heat from cool or heat
        Fach       ! flux from air exchange

     ! 城市参数
     INTEGER, intent(in) :: &
        nurb       ! number of aboveground urban components [-]

     REAL(r8), intent(in) :: &
        hroof,    &! average building height [m]
        hlr,      &! average building height to length of side [-]
        pondmx,   &! maximum ponding of roof/impervious [mm]
        fcover(0:4)! coverage of aboveground urban components [-]

     ! 地面状态
     REAL(r8), intent(in) :: &
        z0h_g,    &! roughness length for bare ground, sensible heat [m]
        obug,     &! monin-obukhov length for bare ground (m)
        ustarg,   &! friction velocity for bare ground [m/s]
        zlnd,     &! roughness length for soil [m]
        zsno,     &! roughness length for snow [m]
        fsno_roof,&! fraction of ground covered by snow [-]
        fsno_gimp,&! fraction of ground covered by snow [-]
        fsno_gper,&! fraction of ground covered by snow [-]
        htvp_roof,&! latent heat of vapor of water (or sublimation) [j/kg]
        htvp_gimp,&! latent heat of vapor of water (or sublimation) [j/kg]
        htvp_gper,&! latent heat of vapor of water (or sublimation) [j/kg]

        troof,    &! temperature of roof [K]
        twsun,    &! temperature of sunlit wall [K]
        twsha,    &! temperature of shaded wall [K]
        tgimp,    &! temperature of impervious road [K]
        tgper,    &! pervious ground temperature [K]

        qroof,    &! roof specific humidity [kg/kg]
        qgimp,    &! imperivous road specific humidity [kg/kg]
        qgper,    &! pervious ground specific humidity [kg/kg]
        dqroofdT, &! d(qroof)/dT
        dqgimpdT, &! d(qgimp)/dT
        dqgperdT   ! d(qgper)/dT

     ! 输出
     REAL(r8), intent(out) :: &
        taux,     &! wind stress: E-W [kg/m/s**2]
        tauy,     &! wind stress: N-S [kg/m/s**2]
        fsenroof, &! sensible heat flux from roof [W/m2]
        fsenwsun, &! sensible heat flux from snulit wall [W/m2]
        fsenwsha, &! sensible heat flux from shaded wall [W/m2]
        fsengimp, &! sensible heat flux from impervious road [W/m2]
        fsengper, &! sensible heat flux from pervious ground [W/m2]
        fevproof, &! evaperation heat flux from roof [W/m2]
        fevpgimp, &! evaperation heat flux from impervious road [W/m2]
        fevpgper, &! evaporation heat flux from pervious ground [mm/s]

        croofs,   &! deriv of roof sensible heat flux wrt soil temp [w/m**2/k]
        cwalls,   &! deriv of wall sensible heat flux wrt soil temp [w/m**2/k]
        cgrnds,   &! deriv of soil sensible heat flux wrt soil temp [w/m**2/k]
        croofl,   &! deriv of roof latent heat flux wrt soil temp [w/m**2/k]
        cgimpl,   &! deriv of gimp latent heat flux wrt soil temp [w/m**2/k]
        cgperl,   &! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
        croof,    &! deriv of roof total heat flux wrt soil temp [w/m**2/k]
        cgimp,    &! deriv of gimp total heat flux wrt soil temp [w/m**2/k]
        cgper,    &! deriv of soil total heat flux wrt soil temp [w/m**2/k]

        tref,     &! 2 m height air temperature [kelvin]
        qref,     &! 2 m height air humidity [kg/kg]

        z0m,      &! effective roughness [m]
        zol,      &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,      &! bulk Richardson number in surface layer
        ustar,    &! friction velocity [m/s]
        tstar,    &! temperature scaling parameter
        qstar,    &! moisture scaling parameter
        fm,       &! integral of profile function for momentum
        fh,       &! integral of profile function for heat
        fq,       &! integral of profile function for moisture
        tafu       ! effective urban air temperature (2nd layer, walls)

!------------------------ LOCAL VARIABLES ------------------------------
     INTEGER ::   &
        niters,   &! maximum number of iterations for surface temperature
        iter,     &! iteration index
        nmozsgn    ! number of times moz changes sign

     REAL(r8) ::  &
        beta,     &! coefficient of conective velocity [-]
        dth,      &! diff of virtual temp. between ref. height and surface
        dqh,      &! diff of humidity between ref. height and surface
        dthv,     &! diff of vir. poten. temp. between ref. height and surface
        obu,      &! monin-obukhov length (m)
        obuold,   &! monin-obukhov length from previous iteration
        ram,      &! aerodynamical resistance [s/m]
        rah,      &! thermal resistance [s/m]
        raw,      &! moisture resistance [s/m]
        raih,     &! temporary variable [kg/m2/s]
        raiw,     &! temporary variable [kg/m2/s]
        fh2m,     &! relation for temperature at 2m
        fq2m,     &! relation for specific humidity at 2m
        fm10m,    &! integral of profile function for momentum at 10m
        thvstar,  &! virtual potential temperature scaling parameter
        um,       &! wind speed including the stablity effect [m/s]
        ur,       &! wind speed at reference height [m/s]
        wc,       &! convective velocity [m/s]
        wc2,      &! wc**2
        zeta,     &! dimensionless height used in Monin-Obukhov theory
        zii,      &! convective boundary height [m]
        zldis,    &! reference height "minus" zero displacement heght [m]
        z0mg,     &! roughness length over ground, momentum [m]
        z0hg,     &! roughness length over ground, sensible heat [m]
        z0qg       ! roughness length over ground, latent heat [m]

     REAL(r8) evplwet, evplwet_dtl, elwmax, elwdif

!----------------------- defination for 3d run ------------------------ !

     INTEGER, parameter :: nlay = 3  ! potential layer number

     INTEGER ::   &
        clev,     &! current layer index
        numlay     ! available layer number

     REAL(r8) ::  &
        ktop,     &! K value at a specific height
        utop,     &! u value at a specific height
        fht,      &! integral of profile function for heat at the top layer
        fqt,      &! integral of profile function for moisture at the top layer
        fmtop,    &! fm value at a specific height
        phih,     &! phi(h), similarity function for sensible heat
        displa,   &! displacement height for urban
        displau,  &! displacement height for urban building
        z0mu,     &! roughless length
        z0hu,     &! roughless length for sensible heat
        z0qu,     &! roughless length for latent heat
        tg,       &! ground temperature
        qg         ! ground specific humidity

     REAL(r8) ::  &
        fg,       &! ground fractional cover
        sqrtdragc,&! sqrt(drag coefficient)
        lm,       &! mix length within canopy
        fai,      &! frontal area index
        fwet,     &! fractional wet area
        delta,    &! 0 or 1
        alpha      ! exponential extinction factor for u/k decline within canopy

     REAL(r8), dimension(0:nurb) :: &
        tu,       &! termperature array
        fc,       &! fractional cover array
        canlev,   &! urban canopy layer lookup table
        rb,       &! leaf boundary layer resistance [s/m]
        cfh,      &! heat conductance for leaf [m/s]
        cfw,      &! latent heat conductance for leaf [m/s]
        wtl0,     &! normalized heat conductance for air and leaf [-]
        wtlq0,    &! normalized latent heat cond. for air and leaf [-]

        ei,       &! vapor pressure on leaf surface [pa]
        deidT,    &! derivative of "ei" on "tl" [pa/K]
        qsatl,    &! leaf specific humidity [kg/kg]
        qsatldT    ! derivative of "qsatl" on "tlef"

     REAL(r8), dimension(nlay) :: &
        fah,      &! weight for thermal resistance to upper layer
        faw,      &! weight for moisture resistance to upper layer
        fgh,      &! weight for thermal resistance to lower layer
        fgw,      &! weight for moisture resistance to lower layer
        ueff_lay, &! effective wind speed within canopy layer [m/s]
        ueff_lay_,&! effective wind speed within canopy layer [m/s]
        taf,      &! air temperature within canopy space [K]
        qaf,      &! humidity of canopy air [kg/kg]
        rd,       &! aerodynamic resistance between layers [s/m]
        rd_,      &! aerodynamic resistance between layers [s/m]
        cah,      &! heat conductance for air [m/s]
        cgh,      &! heat conductance for ground [m/s]
        caw,      &! latent heat conductance for air [m/s]
        cgw,      &! latent heat conductance for ground [m/s]
        wtshi,    &! sensible heat resistance for air, grd and leaf [-]
        wtsqi,    &! latent heat resistance for air, grd and leaf [-]
        wta0,     &! normalized heat conductance for air [-]
        wtg0,     &! normalized heat conductance for ground [-]
        wtaq0,    &! normalized latent heat conductance for air [-]
        wtgq0,    &! normalized heat conductance for ground [-]
        wtll,     &! sum of normalized heat conductance for air and leaf
        wtlql      ! sum of normalized heat conductance for air and leaf

     ! temporal
     INTEGER i
     REAL(r8) bee, tmpw1, tmpw2, fact, facq

!-----------------------End Variable List-------------------------------

! initialization
     tu(0) = troof; tu(1) = twsun; tu(2) = twsha

     fc(:)  = fcover(0:nurb)
     fg     = 1 - fcover(0)
     canlev = (/3, 2, 2/)
     numlay = 2

!-----------------------------------------------------------------------
! initial roughness length for z0mg, z0hg, z0qg
! 计算城市仅地面(不含建筑物、植被)的粗糙度
!TODO: 不透水面的粗糙度怎么定义？

     !TODO: change to original
     !z0mg = (1.-fsno)*zlnd + fsno*zsno
     IF (fsno_gper > 0) THEN
        z0mg = zsno
     ELSE
        z0mg = zlnd
     ENDIF
     z0hg = z0mg
     z0qg = z0mg

!-----------------------------------------------------------------------
! initial saturated vapor pressure and humidity and their derivation
!    0: roof, 1: sunlit wall, 2: shaded wall
!-----------------------------------------------------------------------

     qsatl(0) = qroof
     qsatldT(0) = dqroofdT
     DO i = 1, nurb
        CALL qsadv(tu(i),psrf,ei(i),deiDT(i),qsatl(i),qsatldT(i))
     ENDDO

!-----------------------------------------------------------------------
! 计算加权平均的qg, tg
!-----------------------------------------------------------------------

     ! 设定权重
     fah(1) =  1; fah(2) = fg; fah(3) = fg
     faw(1) =  1; faw(2) = fg; faw(3) = fg
     fgh(1) = fg; fgh(2) = fg; fgh(3) = fg
     fgw(1) = fg; fgw(2) = fg; fgw(3) = fg

     ! 加权后的qg, tg
     tg = ( tgimp*fcover(3) + tgper*fcover(4) ) / fgh(3)
     qg = ( qgimp*fcover(3) + qgper*fcover(4) ) / fgw(3)

!-----------------------------------------------------------------------
! initial for fluxes profile
!-----------------------------------------------------------------------

     nmozsgn = 0     !number of times moz changes sign
     obuold  = 0.    !monin-obukhov length from previous iteration
     zii     = 1000. !m (pbl height)
     beta    = 1.    !- (in computing W_*)

!-----------------------------------------------------------------------
! scaling factor bee
!-----------------------------------------------------------------------
!NOTE: bee value, the default is 1
     bee = 1.

!-----------------------------------------------------------------------
! calculate z0m and displa
!-----------------------------------------------------------------------

     ! Macdonald et al., 1998, Eq. (23), A=4.43
     displau = hroof * (1 + 4.43**(-fcover(0))*(fcover(0) - 1))
     fai  = 4/PI*hlr*fcover(0)
     z0mu = (hroof - displau) * &
          exp( -(0.5*1.2/vonkar/vonkar*(1-displau/hroof)*fai)**(-0.5) )

     ! 比较地面和城市的z0m和displa大小，取大者
     ! maximum assumption
     IF (z0mu < z0mg) z0mu = z0mg

     ! roughness length and displacement height for sensible
     ! and latent heat transfer
     z0m = z0mu

     displa  = displau
     displau = max(hroof/2., displa)

!-----------------------------------------------------------------------
! calculate layer decay coefficient
!-----------------------------------------------------------------------

     ! Raupach, 1992
     sqrtdragc = min( (0.003+0.3*fai)**0.5, 0.3 )

     ! Kondo, 1971
     alpha = hroof/(hroof-displa)/(vonkar/sqrtdragc)

!-----------------------------------------------------------------------
! first guess for taf and qaf for each layer
! a large differece from previous schemes
!-----------------------------------------------------------------------

     IF (numlay .eq. 2) THEN
        taf(3) = (tg + 2.*thm)/3.
        qaf(3) = (qg + 2.*qm )/3.
        taf(2) = (2.*tg + thm)/3.
        qaf(2) = (2.*qg + qm )/3.
     ENDIF

! initialization and input values for Monin-Obukhov
     ! have been set before
     z0hu = z0m; z0qu = z0m
     ur   = max(0.1, sqrt(us*us+vs*vs))    !limit set to 0.1
     dth  = thm - taf(3)
     dqh  =  qm - qaf(3)
     dthv = dth*(1.+0.61*qm) + 0.61*th*dqh

     ! 确保观测高度 >= hroof+10.
     hu = max(hroof+10., hu)
     ht = max(hroof+10., ht)
     hq = max(hroof+10., hq)

     zldis = hu - displa

     IF (zldis <= 0.0) THEN
        write(6,*) 'the obs height of u less than the zero displacement heght'
        CALL abort
     ENDIF

     CALL moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mu,um,obu)

     niters=6

! ======================================================================
!    BEGIN stability iteration
! ======================================================================

     ITERATION : DO iter = 1, niters         !begin stability iteration

!-----------------------------------------------------------------------
! Aerodynamical resistances
!-----------------------------------------------------------------------
! Evaluate stability-dependent variables using moz from prior iteration

        !NOTE: displat=hroof, z0mt=0, are set for roof
        ! fmtop is calculated at the same height of fht, fqt
        CALL moninobukm(hu,ht,hq,displa,z0mu,z0hu,z0qu,obu,um, &
           hroof,0.,ustar,fh2m,fq2m,hroof,fmtop,fm,fh,fq,fht,fqt,phih)

! Aerodynamic resistance
        ! 09/16/2017:
        ! note that for ram, it is the resistance from Href to z0mv+displa
        ! however, for rah and raw is only from Href to canopy effective
        ! exchange height.
        ! for Urban: from Href to roof height
        ! so rah/raw is not comparable with that of 1D case
        ram = 1./(ustar*ustar/um)

        ! 05/02/2016: calculate resistance from the top layer (effective exchange
        ! height) to reference height
        ! for Urban: from roof height to reference height
        rah = 1./(vonkar/(fh-fht)*ustar)
        raw = 1./(vonkar/(fq-fqt)*ustar)

        ! update roughness length for sensible/latent heat
        z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.5e-5)**0.45)
        z0qg = z0hg

        z0hu = max(z0hg, z0hu)
        z0qu = max(z0qg, z0qu)

!-----------------------------------------------------------------------
! new method to calculate rd and ueffect
! the kernel part of 3d model
!-----------------------------------------------------------------------

        ! initialization
        rd(:)  = 0.
        rd_(:) = 0.
        ueff_lay(:)  = 0.
        ueff_lay_(:) = 0.

        ! calculate canopy top wind speed (utop) and exchange coefficient (ktop)
        ! need to update each time as obu changed after each iteration
        utop = ustar/vonkar * fmtop
        ktop = vonkar * (hroof-displa) * ustar / phih

        ueff_lay(3) = utop

        !REAL(r8) FUNCTION kintegral(ktop, fc, bee, alpha, z0mg, &
        !      displah, htop, hbot, obu, ustar, ztop, zbot)
        !rd(3)  = kintegral(ktop, 1., bee, alpha, z0mg, displa/hroof, &
        !   hroof, 0., obug, ustarg, hroof, displa+z0m)

        !REAL(r8) FUNCTION frd(ktop, htop, hbot, &
        !      ztop, zbot, displah, z0h, obu, ustar, &
        !      z0mg, alpha, bee, fc)
        rd(3) = frd(ktop, hroof, 0., hroof, displa+z0m, 0., z0h_g, &
           obug, ustarg, z0mg, alpha, bee, 1.)

        !REAL(r8) FUNCTION uintegral(utop, fc, bee, alpha, z0mg, htop, hbot, ztop, zbot)
        !ueff_lay(2)  = uintegral(utop, 1., bee, alpha, z0mg, hroof, 0., hroof, z0mg)

        !REAL(r8) FUNCTION ueffect(utop, htop, hbot, ztop, zbot, z0mg, alpha, bee, fc)
        ueff_lay(2) = ueffect(utop, hroof, 0., hroof, z0mg, z0mg, alpha, bee, 1.)

        !rd(2)  = kintegral(ktop, 1., bee, alpha, z0mg, displa/hroof, &
        !   hroof, 0., obug, ustarg, displa+z0m, z0qg)
        rd(2) = frd(ktop, hroof, 0., displa+z0m, z0qg, 0., z0h_g, &
           obug, ustarg, z0mg, alpha, bee, 1.)

        !print *, "------------------------"
        !print *, "rd :", rd
        !print *, "rd_:", rd_

!-----------------------------------------------------------------------
! Bulk boundary layer resistance of leaves
!-----------------------------------------------------------------------

        rb(:) = 0.

        DO i = 0, nurb
           clev  = canlev(i)
           rb(i) = rhoair * cpair / ( 11.8 + 4.2*ueff_lay(clev) )
        ENDDO

!-----------------------------------------------------------------------
! dimensional and non-dimensional sensible and latent heat conductances
! for canopy and soil flux calculations.
!-----------------------------------------------------------------------

        !NOTE: 0: roof, 1: sunlit wall, 2: shaded wall,
        !      3: impervious road, 4: pervious road, 5: vegetation
        cfh(:) = 0.
        cfw(:) = 0.

        DO i = 0, nurb
           cfh(i) = 1 / rb(i)
           cfw(i) = 1 / rb(i)
        ENDDO

        ! 为了简单处理，墙面没有水交换
        cfw(1:2) = 0.

        ! initialization
        cah(:) = 0.
        caw(:) = 0.
        cgh(:) = 0.
        cgw(:) = 0.

        ! 计算每层的阻抗
        DO i = 3, 2, -1
           IF (i == 3) THEN
              cah(i) = 1. / rah
              caw(i) = 1. / raw
           ELSE
              cah(i) = 1. / rd(i+1)
              caw(i) = 1. / rd(i+1)
           ENDIF

           cgh(i) = 1. / rd(i)
           cgw(i) = 1. / rd(i)
        ENDDO

        ! claculate wtshi, wtsqi
        wtshi(:) = cah(:)*fah(:) + cgh(:)*fgh(:)
        wtsqi(:) = caw(:)*faw(:) + cgw(:)*fgw(:)

        DO i = 0, nurb
           clev = canlev(i)
           wtshi(clev) = wtshi(clev) + fc(i)*cfh(i)
           wtsqi(clev) = wtsqi(clev) + fc(i)*cfw(i)
        ENDDO

        DO i = 3, 3-numlay+1, -1
           wtshi(i) = 1./wtshi(i)
           wtsqi(i) = 1./wtsqi(i)
        ENDDO

        wta0(:)  = cah(:) * wtshi(:) * fah(:)
        wtg0(:)  = cgh(:) * wtshi(:) * fgh(:)

        wtaq0(:) = caw(:) * wtsqi(:) * faw(:)
        wtgq0(:) = cgw(:) * wtsqi(:) * fgw(:)

        ! calculate wtl0, wtll, wtlq0, wtlql
        wtll(:)  = 0.
        wtlql(:) = 0.

        DO i = 0, nurb
           clev = canlev(i)

           wtl0(i)    = cfh(i) * wtshi(clev) * fc(i)
           wtll(clev) = wtll(clev) + wtl0(i)*tu(i)

           wtlq0(i)    = cfw(i) * wtsqi(clev) * fc(i)
           wtlql(clev) = wtlql(clev) + wtlq0(i)*qsatl(i)
        ENDDO

        IF (numlay .eq. 2) THEN

           ! - Equations:
           ! taf(3) = wta0(3)*thm    + wtg0(3)*taf(2) + wtll(3)
           ! taf(2) = wta0(2)*taf(3) + wtg0(2)*tg     + wtll(2)
           !
           ! qaf(3) = wtaq0(3)*qm     + wtgq0(3)*qaf(2) + wtlql(3)
           ! qaf(2) = wtaq0(2)*qaf(3) + wtgq0(2)*qg     + wtlql(2)

           ! 06/20/2021, yuan: 考虑人为热
           tmpw1  = wta0(3)*thm + wtll(3) &
                  + wtshi(3)*1/3.*(Fhac+Fwst)/rhoair/cpair
           fact   = 1. - wta0(2)*wtg0(3)
           ! 06/20/2021, yuan: 考虑人为热
           taf(2) = (wta0(2)*tmpw1 + wtg0(2)*tg + wtll(2) + &
                     wtshi(2)*(2/3.*(Fhac+Fwst)+Fach)/rhoair/cpair) / fact

           tmpw1  = wtaq0(3)*qm + wtlql(3)
           facq   = 1. - wtaq0(2)*wtgq0(3)
           qaf(2) = (wtaq0(2)*tmpw1 + wtgq0(2)*qg + wtlql(2)) / facq

           qaf(3) = wtaq0(3)*qm + wtgq0(3)*qaf(2) + wtlql(3)
           taf(3) = wta0(3)*thm +  wtg0(3)*taf(2) +  wtll(3)
           taf(3) = taf(3) + wtshi(3)*1/3.*(Fhac+Fwst)/rhoair/cpair

        ENDIF

!-----------------------------------------------------------------------
! Update monin-obukhov length and wind speed including the stability effect
!-----------------------------------------------------------------------

        ! 这里使用的是最高层的taf和qaf
        !TODO: 是否合理，运行单点模型测试
        dth = thm - taf(3)
        dqh =  qm - qaf(3)

        tstar = vonkar/(fh-fht)*dth
        qstar = vonkar/(fq-fqt)*dqh

        thvstar = tstar + 0.61*th*qstar
        zeta = zldis*vonkar*grav*thvstar / (ustar**2*thv)
        IF (zeta .ge. 0.) THEN                             !stable
           zeta = min(2.,max(zeta,1.e-6))
        ELSE                                             !unstable
           zeta = max(-100.,min(zeta,-1.e-6))
        ENDIF
        obu = zldis/zeta

        IF (zeta .ge. 0.) THEN
           um  = max(ur,.1)
        ELSE
           wc  = (-grav*ustar*thvstar*zii/thv)**(1./3.)
           wc2 = beta*beta*(wc*wc)
           um  = sqrt(ur*ur+wc2)
        ENDIF

        IF (obuold*obu .lt. 0.) nmozsgn = nmozsgn+1
        IF (nmozsgn >= 4) EXIT

        obuold = obu

     ENDDO ITERATION                         !end stability iteration

! ======================================================================
!    END stability iteration
! ======================================================================

     z0m = z0mu
     zol = zeta
     rib = min(5.,zol*ustar**2/(vonkar**2/fh*um**2))

     ! sensible heat fluxes
     fsenroof = rhoair*cpair*cfh(0)*(troof-taf(3))
     fsenwsun = rhoair*cpair*cfh(1)*(twsun-taf(2))
     fsenwsha = rhoair*cpair*cfh(2)*(twsha-taf(2))

     ! latent heat fluxes
     fevproof = rhoair*cfw(0)*(qsatl(0)-qaf(3))

     ! fact   = 1. - wta0(2)*wtg0(3)
     ! facq   = 1. - wtaq0(2)*wtgq0(3)
     ! deduce: croofs = rhoair*cpair*cfh(0)*(1.-wtg0(3)*wta0(2)*wtl0(0)/fact-wtl0(0))
     croofs = rhoair*cpair*cfh(0)*(1.-wtl0(0)/fact)
     cwalls = rhoair*cpair*cfh(1)*(1.-wtl0(1)/fact)
     ! deduce: croofl = rhoair*cfw(0)*(1.-wtgq0(3)*wtaq0(2)*wtlq0(0)/facq-wtlq0(0))*qsatldT(0)
     croofl = rhoair*cfw(0)*(1.-wtlq0(0)/facq)*qsatldT(0)

     croof = croofs + croofl*htvp_roof

#if(defined CLMDEBUG)
#endif

     tafu = taf(2)

!-----------------------------------------------------------------------
! wind stresses
!-----------------------------------------------------------------------

     taux = - rhoair*us/ram
     tauy = - rhoair*vs/ram

!-----------------------------------------------------------------------
! fluxes from ground to canopy space
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! 计算城市地面各组分的感热、潜热
!-----------------------------------------------------------------------

     fsengimp = cpair*rhoair*cgh(2)*(tgimp-taf(2))
     fsengper = cpair*rhoair*cgh(2)*(tgper-taf(2))

     fevpgimp = rhoair*cgw(2)*(qgimp-qaf(2))
     fevpgper = rhoair*cgw(2)*(qgper-qaf(2))

!-----------------------------------------------------------------------
! Derivative of soil energy flux with respect to soil temperature (cgrnd)
!-----------------------------------------------------------------------

     cgrnds = cpair*rhoair*cgh(2)*(1.-wtg0(2)/fact)
     cgimpl = rhoair*cgw(2)*(1.-wtgq0(2)/facq)*dqgimpdT
     cgperl = rhoair*cgw(2)*(1.-wtgq0(2)/facq)*dqgperdT

     cgimp  = cgrnds + cgimpl*htvp_gimp
     cgper  = cgrnds + cgperl*htvp_gper

!-----------------------------------------------------------------------
! 2 m height air temperature
!-----------------------------------------------------------------------

     tref = thm + vonkar/(fh-fht)*dth * (fh2m/vonkar - fh/vonkar)
     qref =  qm + vonkar/(fq-fqt)*dqh * (fq2m/vonkar - fq/vonkar)

  END SUBROUTINE UrbanOnlyFlux


  SUBROUTINE  UrbanVegFlux ( &
        ! 模型运行信息
        ipatch      ,deltim                                ,&
        ! 外强迫
        hu          ,ht          ,hq          ,us          ,&
        vs          ,thm         ,th          ,thv         ,&
        qm          ,psrf        ,rhoair      ,frl         ,&
        po2m        ,pco2m       ,par         ,sabv        ,&
        rstfac      ,Fhac        ,Fwst        ,Fach        ,&
        ! 城市和植被参数
        hroof       ,hlr         ,nurb        ,pondmx      ,&
        fcover      ,ewall       ,egimp       ,egper       ,&
        ev          ,htop        ,hbot        ,lai         ,&
        sai         ,sqrtdi      ,effcon      ,vmax25      ,&
        slti        ,hlti        ,shti        ,hhti        ,&
        trda        ,trdm        ,trop        ,gradm       ,&
        binter      ,extkd       ,dewmx       ,etrc        ,&
        ! 地面状态
        z0h_g       ,obug        ,ustarg      ,zlnd        ,&
        zsno        ,fsno_roof   ,fsno_gimp   ,fsno_gper   ,&
        htvp_roof   ,htvp_gimp   ,htvp_gper   ,troof       ,&
        twsun       ,twsha       ,tgimp       ,tgper       ,&
        qroof       ,qgimp       ,qgper       ,dqroofdT    ,&
        dqgimpdT    ,dqgperdT    ,sigf        ,tl          ,&
        ldew                                               ,&
        ! 长波辐射
        Ainv        ,B           ,B1          ,dBdT        ,&
        SkyVF       ,VegVF                                 ,&
        ! 输出
        taux        ,tauy        ,fsenroof    ,fsenwsun    ,&
        fsenwsha    ,fsengimp    ,fsengper    ,fevproof    ,&
        fevpgimp    ,fevpgper    ,croofs      ,cwalls      ,&
        cgrnds      ,croofl      ,cgimpl      ,cgperl      ,&
        croof       ,cgimp       ,cgper       ,fsenl       ,&
        fevpl       ,etr         ,rst         ,assim       ,&
        respc       ,lwsun       ,lwsha       ,lgimp       ,&
        lgper       ,lveg        ,lout        ,tref        ,&
        qref        ,z0m         ,zol         ,rib         ,&
        ustar       ,qstar       ,tstar       ,fm          ,&
        fh          ,fq          ,tafu                      )

!=======================================================================

     USE precision
     USE PhysicalConstants, only: vonkar,grav,hvap,cpair,stefnc
     USE FRICTION_VELOCITY
     USE ASSIM_STOMATA_conductance
     IMPLICIT NONE

!-----------------------Arguments---------------------------------------
     INTEGER,  intent(in) :: &
        ipatch     ! patch index

     REAL(r8), intent(in) :: &
        deltim     ! seconds in a time step [second]

     ! 外强迫
     REAL(r8), intent(inout) :: &
        hu,       &! observational height of wind [m]
        ht,       &! observational height of temperature [m]
        hq         ! observational height of humidity [m]

     REAL(r8), intent(in) :: &
        us,       &! wind component in eastward direction [m/s]
        vs,       &! wind component in northward direction [m/s]
        thm,      &! intermediate variable (tm+0.0098*ht)
        th,       &! potential temperature (kelvin)
        thv,      &! virtual potential temperature (kelvin)
        qm,       &! specific humidity at reference height [kg/kg]
        psrf,     &! pressure at reference height [pa]
        rhoair,   &! density air [kg/m**3]

        frl,      &! atmospheric infrared (longwave) radiation [W/m2]
        par,      &! par absorbed per unit sunlit lai [w/m**2]
        sabv,     &! solar radiation absorbed by vegetation [W/m2]
        rstfac,   &! factor of soil water stress to plant physiologocal processes

        po2m,     &! atmospheric partial pressure  o2 (pa)
        pco2m,    &! atmospheric partial pressure co2 (pa)

        Fhac,     &! flux from heat or cool AC
        Fwst,     &! waste heat from cool or heat
        Fach       ! flux from air exchange

     ! 城市和植被参数
     INTEGER,  intent(in) :: &
        nurb       ! number of aboveground urban components [-]

     REAL(r8), intent(in) :: &
        hroof,    &! average building height [m]
        hlr,      &! average building height to length of side [-]
        pondmx,   &! maximum ponding of roof/impervious [mm]
        fcover(0:5)! coverage of aboveground urban components [-]

     REAL(r8), intent(in) :: &
        ewall,    &! emissivity of walls
        egimp,    &! emissivity of impervious road
        egper,    &! emissivity of pervious road
        ev         ! emissivity of vegetation

     REAL(r8), intent(in) :: &
        htop,     &! PFT crown top height [m]
        hbot,     &! PFT crown bottom height [m]
        lai,      &! adjusted leaf area index for seasonal variation [-]
        sai,      &! stem area index  [-]
        sqrtdi,   &! inverse sqrt of leaf dimension [m**-0.5]

        effcon,   &! quantum efficiency of RuBP regeneration (mol CO2 / mol quanta)
        vmax25,   &! maximum carboxylation rate at 25 C at canopy top
                   ! the range : 30.e-6 <-> 100.e-6 (mol co2 m-2 s-1)
        shti,     &! slope of high temperature inhibition function     (s1)
        hhti,     &! 1/2 point of high temperature inhibition function (s2)
        slti,     &! slope of low temperature inhibition function      (s3)
        hlti,     &! 1/2 point of low temperature inhibition function  (s4)
        trda,     &! temperature coefficient in gs-a model             (s5)
        trdm,     &! temperature coefficient in gs-a model             (s6)
        trop,     &! temperature coefficient in gs-a model         (273+25)
        gradm,    &! conductance-photosynthesis slope parameter
        binter,   &! conductance-photosynthesis intercept

        extkd,    &! diffuse and scattered diffuse PAR extinction coefficient
        dewmx,    &! maximum dew
        etrc       ! maximum possible transpiration rate (mm/s)

     ! 地面状态
     REAL(r8), intent(in) :: &
        z0h_g,    &! roughness length for bare ground, sensible heat [m]
        obug,     &! monin-obukhov length for bare ground (m)
        ustarg,   &! friction velocity for bare ground [m/s]
        zlnd,     &! roughness length for soil [m]
        zsno,     &! roughness length for snow [m]
        fsno_roof,&! fraction of ground covered by snow
        fsno_gimp,&! fraction of ground covered by snow
        fsno_gper,&! fraction of ground covered by snow
        htvp_roof,&! latent heat of vapor of water (or sublimation) [j/kg]
        htvp_gimp,&! latent heat of vapor of water (or sublimation) [j/kg]
        htvp_gper,&! latent heat of vapor of water (or sublimation) [j/kg]

        troof,    &! temperature of roof [K]
        twsun,    &! temperature of sunlit wall [K]
        twsha,    &! temperature of shaded wall [K]
        tgimp,    &! temperature of impervious road [K]
        tgper,    &! pervious ground temperature [K]

        qroof,    &! roof specific humidity [kg/kg]
        qgimp,    &! imperivous road specific humidity [kg/kg]
        qgper,    &! pervious ground specific humidity [kg/kg]
        dqroofdT, &! d(qroof)/dT
        dqgimpdT, &! d(qgimp)/dT
        dqgperdT, &! d(qgper)/dT
        sigf       !

     REAL(r8), intent(inout) :: &
        tl,       &! leaf temperature [K]
        ldew       ! depth of water on foliage [mm]

     REAL(r8), intent(in) :: Ainv(5,5)     !Inverse of Radiation transfer matrix
     REAL(r8), intent(in) :: SkyVF(5)      !View factor to sky
     REAL(r8), intent(in) :: VegVF(5)      !View factor to veg
     REAL(r8), intent(inout) :: B(5)       !Vectors of incident radition on each surface
     REAL(r8), intent(inout) :: B1(5)      !Vectors of incident radition on each surface
     REAL(r8), intent(inout) :: dBdT(5)    !Vectors of incident radition on each surface

     REAL(r8), intent(out) :: &
        taux,     &! wind stress: E-W [kg/m/s**2]
        tauy,     &! wind stress: N-S [kg/m/s**2]
        fsenroof, &! sensible heat flux from roof [W/m2]
        fsenwsun, &! sensible heat flux from sunlit wall [W/m2]
        fsenwsha, &! sensible heat flux from shaded wall [W/m2]
        fsengimp, &! sensible heat flux from impervious road [W/m2]
        fsengper, &! sensible heat flux from pervious ground [W/m2]
        fevproof, &! evaporation heat flux from roof [mm/s]
        fevpgimp, &! evaporation heat flux from impervious road [mm/s]
        fevpgper, &! evaporation heat flux from pervious ground [mm/s]

        croofs,   &! deriv of roof sensible heat flux wrt soil temp [w/m**2/k]
        cwalls,   &! deriv of wall sensible heat flux wrt soil temp [w/m**2/k]
        cgrnds,   &! deriv of ground latent heat flux wrt soil temp [w/m**2/k]
        croofl,   &! deriv of roof latent heat flux wrt soil temp [w/m**2/k]
        cgimpl,   &! deriv of impervious latent heat flux wrt soil temp [w/m**2/k]
        cgperl,   &! deriv of soil atent heat flux wrt soil temp [w/m**2/k]
        croof,    &! deriv of roof total flux wrt soil temp [w/m**2/k]
        cgimp,    &! deriv of impervious total heat flux wrt soil temp [w/m**2/k]
        cgper,    &! deriv of soil total heat flux wrt soil temp [w/m**2/k]

        tref,     &! 2 m height air temperature [kelvin]
        qref       ! 2 m height air humidity

     REAL(r8), intent(out) :: &
        fsenl,    &! sensible heat from leaves [W/m2]
        fevpl,    &! evaporation+transpiration from leaves [mm/s]
        etr,      &! transpiration rate [mm/s]
        rst,      &! stomatal resistance
        assim,    &! rate of assimilation
        respc      ! rate of respiration

     REAL(r8), intent(inout) :: &
        lwsun,    &! net longwave radiation of sunlit wall
        lwsha,    &! net longwave radiation of shaded wall
        lgimp,    &! net longwave radiation of impervious road
        lgper,    &! net longwave radiation of pervious road
        lveg,     &! net longwave radiation of vegetation
        lout       ! out-going longwave radiation

     REAL(r8), intent(inout) :: &
        z0m,      &! effective roughness [m]
        zol,      &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,      &! bulk Richardson number in surface layer
        ustar,    &! friction velocity [m/s]
        tstar,    &! temperature scaling parameter
        qstar,    &! moisture scaling parameter
        fm,       &! integral of profile function for momentum
        fh,       &! integral of profile function for heat
        fq,       &! integral of profile function for moisture
        tafu       ! effective urban air temperature (2nd layer, walls)

!-----------------------Local Variables---------------------------------
! assign iteration parameters
     INTEGER, parameter :: itmax  = 40   !maximum number of iteration
     INTEGER, parameter :: itmin  = 6    !minimum number of iteration
     REAL(r8),parameter :: delmax = 3.0  !maximum change in leaf temperature [K]
     REAL(r8),parameter :: dtmin  = 0.01 !max limit for temperature convergence [K]
     REAL(r8),parameter :: dlemin = 0.1  !max limit for energy flux convergence [w/m2]

     REAL(r8) dtl(0:itmax+1)             !difference of tl between two iterative step

     REAL(r8) ::  &
        zldis,    &! reference height "minus" zero displacement heght [m]
        zii,      &! convective boundary layer height [m]
        z0mv,     &! roughness length, momentum [m]
        z0mu,     &! roughness length, momentum [m]
        z0hu,     &! roughness length, sensible heat [m]
        z0qu,     &! roughness length, latent heat [m]
        zeta,     &! dimensionless height used in Monin-Obukhov theory
        beta,     &! coefficient of conective velocity [-]
        wc,       &! convective velocity [m/s]
        wc2,      &! wc**2
        dth,      &! diff of virtual temp. between ref. height and surface
        dthv,     &! diff of vir. poten. temp. between ref. height and surface
        dqh,      &! diff of humidity between ref. height and surface
        obu,      &! monin-obukhov length (m)
        um,       &! wind speed including the stablity effect [m/s]
        ur,       &! wind speed at reference height [m/s]
        uaf,      &! velocity of air within foliage [m/s]
        fh2m,     &! relation for temperature at 2m
        fq2m,     &! relation for specific humidity at 2m
        fm10m,    &! integral of profile function for momentum at 10m
        thvstar,  &! virtual potential temperature scaling parameter
        eah,      &! canopy air vapor pressure (pa)
        pco2g,    &! co2 pressure (pa) at ground surface (pa)
        pco2a,    &! canopy air co2 pressure (pa)

        ram,      &! aerodynamical resistance [s/m]
        rah,      &! thermal resistance [s/m]
        raw,      &! moisture resistance [s/m]
        clai,     &! canopy heat capacity [Jm-2K-1]
        del,      &! absolute change in leaf temp in current iteration [K]
        del2,     &! change in leaf temperature in previous iteration [K]
        dele,     &! change in heat fluxes from leaf [K]
        dele2,    &! change in heat fluxes from leaf [K]
        det,      &! maximum leaf temp. change in two consecutive iter [K]
        dee,      &! maximum leaf temp. change in two consecutive iter [K]

        obuold,   &! monin-obukhov length from previous iteration
        tlbef,    &! leaf temperature from previous iteration [K]
        err,      &! balance error

        rs,       &! sunlit leaf stomatal resistance [s/m]
        rsoil,    &! soil respiration
        gah2o,    &! conductance between canopy and atmosphere
        gdh2o,    &! conductance between canopy and ground
        tprcor     ! tf*psur*100./1.013e5

     INTEGER it, nmozsgn

     REAL(r8) evplwet, evplwet_dtl, etr_dtl, elwmax, elwdif
     REAL(r8) irab, dirab_dtl, fsenl_dtl, fevpl_dtl
     REAL(r8) z0mg, z0hg, z0qg, cint(3)
     REAL(r8) fevpl_bef, fevpl_noadj, dtl_noadj, erre

!----------------------- defination for 3d run ------------------------ !
     INTEGER, parameter :: nlay = 3
     INTEGER, parameter :: avec(5) = (/0,0,0,0,1/) !unit vector

     INTEGER ::   &
        clev,     &! current layer index
        botlay,   &! botom layer index
        numlay     ! available layer number

     REAL(r8) ::  &
        ktop,     &! K value at a specific height
        utop,     &! u value at a specific height
        fht,      &! integral of profile function for heat at the top layer
        fqt,      &! integral of profile function for moisture at the top layer
        fmtop,    &! fm value at a specific height
        phih,     &! phi(h), similarity function for sensible heat
        displa,   &! displacement height for urban
        displau,  &! displacement height for urban building
        displav,  &! displacement height for urban vegetation
        displav_lay,&!displacement height for urban vegetation layer
        z0mv_lay, &! roughless length for vegetation
        ueff_veg, &! effective wind speed within canopy layer [m/s]
        tg,       &! ground temperature
        qg         ! ground specific humidity

     REAL(r8) ::  &
        fg,       &! ground fractional cover
        sqrtdragc,&! sqrt(drag coefficient)
        lm,       &! mix length within canopy
        fai,      &! frontal area index
        lsai,     &! lai+sai
        fwet,     &! fractional wet area
        delta,    &! 0 or 1
        alpha      ! exponential extinction factor for u/k decline within canopy

     REAL(r8) ::  &
        lwsun_bef,&! change of lw for the last time
        lwsha_bef,&! change of lw for the last time
        lgimp_bef,&! change of lw for the last time
        lgper_bef,&! change of lw for the last time
        lveg_bef   ! change of lw for the last time

     REAL(r8), dimension(0:nurb) :: &
        tu,       &! termperature array
        fc,       &! fractional cover array
        canlev,   &! urban canopy layer lookup table
        rb,       &! leaf boundary layer resistance [s/m]
        cfh,      &! heat conductance for leaf [m/s]
        cfw,      &! latent heat conductance for leaf [m/s]
        wtl0,     &! normalized heat conductance for air and leaf [-]
        wtlq0,    &! normalized latent heat cond. for air and leaf [-]

        ei,       &! vapor pressure on leaf surface [pa]
        deidT,    &! derivative of "ei" on "tl" [pa/K]
        qsatl,    &! leaf specific humidity [kg/kg]
        qsatldT    ! derivative of "qsatl" on "tlef"

     REAL(r8), dimension(nlay) :: &
        fah,      &! weight for thermal resistance to upper layer
        faw,      &! weight for moisture resistance to upper layer
        fgh,      &! weight for thermal resistance to lower layer
        fgw,      &! weight for moisture resistance to lower layer
        ueff_lay, &! effective wind speed within canopy layer [m/s]
        ueff_lay_,&! effective wind speed within canopy layer [m/s]
        taf,      &! air temperature within canopy space [K]
        qaf,      &! humidity of canopy air [kg/kg]
        rd,       &! aerodynamic resistance between layers [s/m]
        rd_,      &! aerodynamic resistance between layers [s/m]
        cah,      &! heat conductance for air [m/s]
        cgh,      &! heat conductance for ground [m/s]
        caw,      &! latent heat conductance for air [m/s]
        cgw,      &! latent heat conductance for ground [m/s]
        wtshi,    &! sensible heat resistance for air, grd and leaf [-]
        wtsqi,    &! latent heat resistance for air, grd and leaf [-]
        wta0,     &! normalized heat conductance for air [-]
        wtg0,     &! normalized heat conductance for ground [-]
        wtaq0,    &! normalized latent heat conductance for air [-]
        wtgq0,    &! normalized heat conductance for ground [-]
        wtll,     &! sum of normalized heat conductance for air and leaf
        wtlql      ! sum of normalized heat conductance for air and leaf

     ! temporal
     INTEGER i
     REAL(r8) bee, cf, tmpw1, tmpw2, fact, facq
     REAL(r8) B_5, B1_5, dBdT_5, X(5), dX(5)

!-----------------------End Variable List-------------------------------

! initialization of errors and  iteration parameters
     it    = 1    !counter for leaf temperature iteration
     del   = 0.0  !change in leaf temperature from previous iteration
     dele  = 0.0  !latent head flux from leaf for previous iteration

     dtl   = 0.
     fevpl_bef = 0.

! initial values for z0hg, z0qg

     !TODO: change to original
     !z0mg = (1.-fsno)*zlnd + fsno*zsno
     IF (fsno_gper > 0) THEN
        z0mg = zsno
     ELSE
        z0mg = zlnd
     ENDIF
     z0hg = z0mg
     z0qg = z0mg

!-----------------------------------------------------------------------
! scaling-up coefficients from leaf to canopy
!-----------------------------------------------------------------------

     cint(1) = (1.-exp(-0.110*lai))/0.110
     cint(2) = (1.-exp(-extkd*lai))/extkd
     cint(3) = lai

!-----------------------------------------------------------------------
! initial saturated vapor pressure and humidity and their derivation
!-----------------------------------------------------------------------

     !clai = 4.2 * 1000. * 0.2
     clai = 0.0
     lsai = lai + sai

     ! index 0:roof, 1:sunlit wall, 2:shaded wall, 3: vegetation
     tu(0) = troof; tu(1) = twsun; tu(2) = twsha; tu(3) = tl

     fc(:)  = fcover(0:nurb)
     fc(3)  = fcover(5)
     fg     = 1 - fcover(0)
     canlev = (/3, 2, 2, 1/)

     B_5    = B(5)
     B1_5   = B1(5)
     dBdT_5 = dBdT(5)

     CALL dewfraction (sigf,lai,sai,dewmx,ldew,fwet)

     qsatl(0) = qroof
     qsatldT(0) = dqroofDT
     DO i = 1, nurb
        CALL qsadv(tu(i),psrf,ei(i),deiDT(i),qsatl(i),qsatldT(i))
     ENDDO

     ! 保留上次长波辐射
     lwsun_bef = lwsun
     lwsha_bef = lwsha
     lgimp_bef = lgimp
     lgper_bef = lgper
     lveg_bef  = lveg
!-----------------------------------------------------------------------
! 计算加权平均的qg, tg
!-----------------------------------------------------------------------

     !TODO: no wet
     fah(1) =  1; fah(2) = fg; fah(3) = fg
     faw(1) =  1; faw(2) = fg; faw(3) = fg
     fgh(1) = fg; fgh(2) = fg; fgh(3) = fg
     fgw(1) = fg; fgw(2) = fg; fgw(3) = fg

     ! 加权后的qg, tg
     tg = ( tgimp*fcover(3) + tgper*fcover(4) ) / fgh(3)
     qg = ( qgimp*fcover(3) + qgper*fcover(4) ) / fgw(3)

!-----------------------------------------------------------------------
! initial for fluxes profile
!-----------------------------------------------------------------------

     nmozsgn = 0     !number of times moz changes sign
     obuold  = 0.    !monin-obukhov length from previous iteration
     zii     = 1000. !m (pbl height)
     beta    = 1.    !- (in computing W_*)

!-----------------------------------------------------------------------
! scaling factor bee
!-----------------------------------------------------------------------
!NOTE: bee value, the default is 1
     bee = 1.

!-----------------------------------------------------------------------
! calculate z0m and displa for layers
!-----------------------------------------------------------------------

     ! 计算自身和整个面积的z0和displa (不考虑建筑物的存在)
     CALL cal_z0_displa(lsai, htop, 1., z0mv, displav)
     CALL cal_z0_displa(lsai, htop, fc(3), z0mv_lay, displav_lay)

     ! Macdonald et al., 1998, Eq. (23), A=4.43
     displau = hroof * (1 + 4.43**(-fcover(0))*(fcover(0) - 1))
     fai  = 4/PI*hlr*fcover(0)
     z0mu = (hroof - displau) * &
        exp( -(0.5*1.2/vonkar/vonkar*(1-displau/hroof)*fai)**(-0.5) )

     ! 比较植被、裸地和建筑物的z0m和displa大小，取大者
     ! maximum assumption
     ! 11/26/2021, yuan: remove the below
     !IF (z0mu < z0mv_lay) z0mu = z0mv_lay
     IF (z0mu < z0mg) z0mu = z0mg
     IF (displau < displav_lay) displau = displav_lay

     displa = displau
     z0m    = z0mu

     ! 层次设定
     IF (z0mv+displav > z0mu+displau) THEN
        numlay = 2; botlay = 2; canlev(3) = 2
     ELSE
        numlay = 3; botlay = 1
     ENDIF

!-----------------------------------------------------------------------
! calculate layer decay coefficient
!-----------------------------------------------------------------------

     ! initialization
     sqrtdragc = min( (0.003+0.3*fai)**0.5, 0.3 )

     ! Kondo, 1971
     alpha = hroof/(hroof-displa)/(vonkar/sqrtdragc)

     displau = max(hroof/2., displau)

!-----------------------------------------------------------------------
! first guess for taf and qaf for each layer
! a large differece from previous schemes
!-----------------------------------------------------------------------
     taf(:) = 0.
     qaf(:) = 0.

     IF (numlay .eq. 2) THEN
        taf(3) = (tg + 2.*thm)/3.
        qaf(3) = (qg + 2.*qm )/3.
        taf(2) = (2.*tg + thm)/3.
        qaf(2) = (2.*qg + qm )/3.
     ENDIF

     IF (numlay .eq. 3) THEN
        taf(3) = (tg + 3.*thm)/4.
        qaf(3) = (qg + 3.*qm )/4.
        taf(2) = (tg + thm)/2.
        qaf(2) = (qg + qm )/2.
        taf(1) = (3.*tg + thm)/4.
        qaf(1) = (3.*qg + qm )/4.
     ENDIF

!-----------------------------------------------------------------------
! some environment variables
! how to calculate rsoil and what is its usage?
!-----------------------------------------------------------------------
     pco2a = pco2m
     tprcor = 44.6*273.16*psrf/1.013e5
     rsoil = 0.   !respiration (mol m-2 s-1)
     !rsoil = 1.22e-6*exp(308.56*(1./56.02-1./(tg-227.13)))
     !rsoil = rstfac * 0.23 * 15. * 2.**((tg-273.16-10.)/10.) * 1.e-6
     !rsoil = 5.22 * 1.e-6
     rsoil = 0.22 * 1.e-6

! initialization and input values for Monin-Obukhov
     ! have been set before
     z0hu = z0m; z0qu = z0m
     ur = max(0.1, sqrt(us*us+vs*vs))    !limit set to 0.1
     dth = thm - taf(3)
     dqh =  qm - qaf(3)
     dthv = dth*(1.+0.61*qm) + 0.61*th*dqh

     ! 确保观测高度 >= hroof+10.
     hu = max(hroof+10., hu)
     ht = max(hroof+10., ht)
     hq = max(hroof+10., hq)

     zldis = hu - displa

     IF (zldis <= 0.0) THEN
        write(6,*) 'the obs height of u less than the zero displacement heght'
        CALL abort
     ENDIF

     CALL moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mu,um,obu)

! ======================================================================
!    BEGIN stability iteration
! ======================================================================

     DO WHILE (it .le. itmax)

        tlbef = tl

        del2  = del
        dele2 = dele

!-----------------------------------------------------------------------
! Aerodynamical resistances
!-----------------------------------------------------------------------
! Evaluate stability-dependent variables using moz from prior iteration

        CALL moninobukm(hu,ht,hq,displa,z0mu,z0hu,z0qu,obu,um, &
           hroof,0.,ustar,fh2m,fq2m,hroof,fmtop,fm,fh,fq,fht,fqt,phih)

! Aerodynamic resistance
        ! 09/16/2017:
        ! note that for ram, it is the resistance from Href to z0mu+displa
        ! however, for rah and raw is only from Href to canopy effective
        ! exchange height.
        ! so rah/raw is not comparable with that of 1D case
        ram = 1./(ustar*ustar/um)

        ! 05/02/2016: calculate resistance from the top layer (effective exchange
        ! height) to reference height
        ! for urban, from roof height to reference height
        rah = 1./(vonkar/(fh-fht)*ustar)
        raw = 1./(vonkar/(fq-fqt)*ustar)

! update roughness length for sensible/latent heat
        z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.5e-5)**0.45)
        z0qg = z0hg

        z0hu = max(z0hg, z0hu)
        z0qu = max(z0qg, z0qu)

!-----------------------------------------------------------------------
! new method to calculate rd and ueffect
! the kernel part of 3d model
!-----------------------------------------------------------------------

        ! initialization
        rd(:)  = 0.
        rd_(:) = 0.
        ueff_lay(:)  = 0.
        ueff_lay_(:) = 0.

        ! calculate canopy top wind speed (utop) and exchange coefficient (ktop)
        ! need to update each time as obu changed after each iteration
        utop = ustar/vonkar * fmtop
        ktop = vonkar * (hroof-displa) * ustar / phih

        ueff_lay(3)  = utop
        ueff_lay_(3) = utop

        ! REAL(r8) FUNCTION kintegral(ktop, fc, bee, alpha, z0mg, &
        !      displah, htop, hbot, obu, ustar, ztop, zbot)
        !rd(3)  = kintegral(ktop, 1., bee, alpha, z0mg, displau/hroof, &
        !   hroof, 0., obug, ustarg, hroof, displau+z0mu)

        ! REAL(r8) FUNCTION frd(ktop, htop, hbot, &
        !      ztop, zbot, displah, z0h, obu, ustar, &
        !      z0mg, alpha, bee, fc)
        rd(3) = frd(ktop, hroof, 0., hroof, displau+z0mu, 0., z0h_g, &
           obug, ustarg, z0mg, alpha, bee, 1.)

        ! REAL(r8) FUNCTION uintegral(utop, fc, bee, alpha, z0mg, htop, hbot, ztop, zbot)
        !ueff_lay(2)  = uintegral(utop, 1., bee, alpha, z0mg, hroof, 0., hroof, z0mg)

        ! REAL(r8) FUNCTION ueffect(utop, htop, hbot, &
        !      ztop, zbot, z0mg, alpha, bee, fc)
        ueff_lay(2) = ueffect(utop, hroof, 0., hroof, z0mg, z0mg, alpha, bee, 1.)

        IF (numlay == 3) THEN
           ! REAL(r8) FUNCTION kintegral(ktop, fc, bee, alpha, z0mg, &
           !      displah, htop, hbot, obu, ustar, ztop, zbot)
           !rd(2)  = kintegral(ktop, 1., bee, alpha, z0mg, displau/hroof, &
           !   hroof, 0., obug, ustarg, displau+z0mu, displav+z0mv)
           rd(2) = frd(ktop, hroof, 0., displau+z0mu, displav+z0mv,0., z0h_g, &
              obug, ustarg, z0mg, alpha, bee, 1.)

           !rd(1)  = kintegral(ktop, 1., bee, alpha, z0mg, displau/hroof, &
           !   hroof, 0., obug, ustarg, displav+z0mv, z0qg)
           rd(1) = frd(ktop, hroof, 0., displav+z0mv, z0qg, 0., z0h_g, &
              obug, ustarg, z0mg, alpha, bee, 1.)
        ELSE
           !rd(2)  = kintegral(ktop, 1., bee, alpha, z0mg, displau/hroof, &
           !   hroof, 0., obug, ustarg, displau+z0mu, z0qg)
           rd(2) = frd(ktop, hroof, 0., displau+z0mu, z0qg, 0., z0h_g, &
              obug, ustarg, z0mg, alpha, bee, 1.)
        ENDIF

        !ueff_lay(2)  = uintegral(utop, 1., bee, alpha, z0mg, hroof, 0., hroof, z0mg)
        !print *, "htop/hbot:", htop, hbot  !fordebug
        !ueff_veg  = uintegral(utop, 1., bee, alpha, z0mg, hroof, 0., htop, hbot)

        !ueff_lay_(2) = ueffect(utop, hroof, 0., hroof, z0mg, z0mg, alpha, bee, 1.)
        ueff_veg = ueffect(utop, hroof, 0., htop, hbot, z0mg, alpha, bee, 1.)

        !print *, "ueff_lay :", ueff_lay
        !print *, "ueff_lay_:", ueff_lay_
        !print *, "------------------------"
        !print *, "rd :", rd
        !print *, "rd_:", rd_

!-----------------------------------------------------------------------
! Bulk boundary layer resistance of leaves
!-----------------------------------------------------------------------
        rb(:) = 0.

        DO i = 0, nurb
           IF (i == 3) THEN
              cf = 0.01*sqrtdi*sqrt(ueff_veg)
              rb(i) = 1./cf
              cycle
           ENDIF
           clev = canlev(i)
           rb(i) = rhoair * cpair / ( 11.8 + 4.2*ueff_lay(clev) )
        ENDDO

!-----------------------------------------------------------------------
! stomatal resistances
!-----------------------------------------------------------------------

        IF (lai > 0.) THEN

           rb = rb / lai

           clev = canlev(3)
           eah = qaf(clev) * psrf / ( 0.622 + 0.378 * qaf(clev) )    !pa

!-----------------------------------------------------------------------
! note: calculate resistance for leaves
!-----------------------------------------------------------------------
           CALL stomata (vmax25,effcon ,slti   ,hlti   ,&
              shti    ,hhti    ,trda   ,trdm   ,trop   ,&
              gradm   ,binter  ,thm    ,psrf   ,po2m   ,&
              pco2m   ,pco2a   ,eah    ,ei(3)  ,tu(3)  ,&
              par     ,rb(3)   ,raw    ,rstfac ,cint(:),&
              assim   ,respc   ,rs     )
        ELSE
           rs = 2.e4; assim = 0.; respc = 0.
        ENDIF

! above stomatal resistances are for the canopy, the stomatal rsistances
! and the "rb" in the following calculations are the average for single leaf. thus,
        rs = rs * lai

!-----------------------------------------------------------------------
! dimensional and non-dimensional sensible and latent heat conductances
! for canopy and soil flux calculations.
!-----------------------------------------------------------------------

        cfh(:) = 0.
        cfw(:) = 0.

        DO i = 0, nurb

           IF (i == 3) THEN

              clev = canlev(i)
              delta = 0.0
              IF (qsatl(i)-qaf(clev) .gt. 0.) delta = 1.0

              ! 计算感热阻抗
              cfh(i) = lsai / rb(i)

              ! for building walls, cfw=0., no water transfer
              ! for canopy, keep the same but for one leaf
              ! 计算潜热阻抗
              cfw(i) = (1.-delta*(1.-fwet))*lsai/rb(i) + &
                 (1.-fwet)*delta* ( lai/(rb(i)+rs) )
           ELSE
              cfh(i) = 1 / rb(i)
              cfw(i) = 1 / rb(i)
           ENDIF
        ENDDO

        ! 为了简单处理，墙面没有水交换
        cfw(1:2) = 0.

        ! initialization
        cah(:) = 0.
        caw(:) = 0.
        cgh(:) = 0.
        cgw(:) = 0.

        ! 计算每层的阻抗
        DO i = 3, botlay, -1
           IF (i == 3) THEN
              cah(i) = 1. / rah
              caw(i) = 1. / raw
           ELSE
              cah(i) = 1. / rd(i+1)
              caw(i) = 1. / rd(i+1)
           ENDIF

           cgh(i) = 1. / rd(i)
           cgw(i) = 1. / rd(i)
        ENDDO

        ! claculate wtshi, wtsqi
        wtshi(:) = cah(:)*fah(:) + cgh(:)*fgh(:)
        wtsqi(:) = caw(:)*faw(:) + cgw(:)*fgw(:)

        DO i = 0, nurb
           clev = canlev(i)
           wtshi(clev) = wtshi(clev) + fc(i)*cfh(i)
           wtsqi(clev) = wtsqi(clev) + fc(i)*cfw(i)
        ENDDO

        DO i = 3, 3-numlay+1, -1
           wtshi(i) = 1./wtshi(i)
           wtsqi(i) = 1./wtsqi(i)
        ENDDO

        wta0(:) = cah(:) * wtshi(:) * fah(:)
        wtg0(:) = cgh(:) * wtshi(:) * fgh(:)

        wtaq0(:) = caw(:) * wtsqi(:) * faw(:)
        wtgq0(:) = cgw(:) * wtsqi(:) * fgw(:)

        ! calculate wtl0, wtll, wtlq0, wtlql
        wtll(:)  = 0.
        wtlql(:) = 0.

        DO i = 0, nurb
           clev = canlev(i)

           wtl0(i)  = cfh(i) * wtshi(clev) * fc(i)
           wtll(clev) = wtll(clev) + wtl0(i)*tu(i)

           wtlq0(i) = cfw(i) * wtsqi(clev) * fc(i)
           wtlql(clev) = wtlql(clev) + wtlq0(i)*qsatl(i)
        ENDDO

        ! 根据层数来计算空气温度、湿度
        ! to solve taf(:) and qaf(:)

        IF (numlay .eq. 2) THEN

           ! - Equations:
           ! taf(3) = wta0(3)*thm    + wtg0(3)*taf(2) + wtll(3)
           ! taf(2) = wta0(2)*taf(3) + wtg0(2)*tg     + wtll(2)
           !
           ! qaf(3) = wtaq0(3)*qm     + wtgq0(3)*qaf(2) + wtlql(3)
           ! qaf(2) = wtaq0(2)*qaf(3) + wtgq0(2)*qg     + wtlql(2)

           ! 06/20/2021, yuan: 考虑人为热
           tmpw1  = wta0(3)*thm + wtll(3) &
                  + wtshi(3)*1/3.*(Fhac+Fwst)/rhoair/cpair
           fact   = 1. - wta0(2)*wtg0(3)
           ! 06/20/2021, yuan: 考虑人为热
           taf(2) = (wta0(2)*tmpw1 + wtg0(2)*tg + wtll(2) + &
                     wtshi(2)*(2/3.*(Fhac+Fwst)+Fach)/rhoair/cpair) / fact

           tmpw1  = wtaq0(3)*qm + wtlql(3)
           facq   = 1. - wtaq0(2)*wtgq0(3)
           qaf(2) = (wtaq0(2)*tmpw1 + wtgq0(2)*qg + wtlql(2)) / facq

           qaf(3) = wtaq0(3)*qm + wtgq0(3)*qaf(2) + wtlql(3)
           taf(3) = wta0(3)*thm +  wtg0(3)*taf(2) +  wtll(3)
           taf(3) = taf(3) + wtshi(3)*1/3.*(Fhac+Fwst)/rhoair/cpair

        ENDIF

        IF (numlay .eq. 3) THEN

           ! - Equations:
           ! taf(3) = wta0(3)*thm    + wtg0(3)*taf(2) + wtll(3)
           ! taf(2) = wta0(2)*taf(3) + wtg0(2)*taf(1) + wtll(2)
           ! taf(1) = wta0(1)*taf(2) + wtg0(1)*tg     + wtll(1)
           !
           ! qaf(3) = wtaq0(3)*qm     + wtgq0(3)*qaf(2) + wtlql(3)
           ! qaf(2) = wtaq0(2)*qaf(3) + wtaq0(2)*qaf(1) + wtlql(2)
           ! qaf(1) = wtaq0(1)*qaf(2) + wtaq0(1)*qg     + wtlql(1)

           tmpw1  = wta0(3)*thm + wtll(3) &
                  + wtshi(3)*1/3.*(Fhac+Fwst)/rhoair/cpair
           tmpw2  = wtg0(1)*tg  + wtll(1) &
                  + wtshi(1)*1/3.*(Fhac+Fwst)/rhoair/cpair
           fact   = 1. - wta0(2)*wtg0(3) - wtg0(2)*wta0(1)
           taf(2) = (wta0(2)*tmpw1 + wtg0(2)*tmpw2 + wtll(2) + &
                     wtshi(2)*(1/3.*(Fhac+Fwst)+Fach)/rhoair/cpair) / fact

           tmpw1  = wtaq0(3)*qm + wtlql(3)
           tmpw2  = wtgq0(1)*qg + wtlql(1)
           facq   = 1. - wtaq0(2)*wtgq0(3) - wtgq0(2)*wtaq0(1)
           qaf(2) = (wtaq0(2)*tmpw1 + wtgq0(2)*tmpw2 + wtlql(2)) / facq

           qaf(1) = wtaq0(1)*qaf(2) + wtgq0(1)*qg + wtlql(1)
           taf(1) =  wta0(1)*taf(2) +  wtg0(1)*tg +  wtll(1) &
                  + wtshi(1)*1/3.*(Fhac+Fwst)/rhoair/cpair

           qaf(3) = wtaq0(3)*qm + wtgq0(3)*qaf(2) + wtlql(3)
           taf(3) = wta0(3)*thm +  wtg0(3)*taf(2) +  wtll(3) &
                  + wtshi(3)*1/3.*(Fhac+Fwst)/rhoair/cpair

        ENDIF

!-----------------------------------------------------------------------
! IR radiation, sensible and latent heat fluxes and their derivatives
!-----------------------------------------------------------------------
! the partial derivatives of areodynamical resistance are ignored
! which cannot be determined analtically

        !NOTE: ONLY for vegetation
        i = 3

! sensible heat fluxes and their derivatives
        fsenl = rhoair * cpair * cfh(i) * (tl - taf(botlay))

        ! 09/24/2017: why fact/facq here? bugs? YES
        ! 09/25/2017: re-written, check it clearfully
        ! 11/25/2021: re-written, double check
        IF (botlay == 2) THEN
           fsenl_dtl = rhoair * cpair * cfh(i) * (1.-wtl0(i)/fact)
        ELSE
           fsenl_dtl = rhoair * cpair * cfh(i) * (1.-wta0(1)*wtg0(2)*wtl0(i)/fact-wtl0(i))
        ENDIF


! latent heat fluxes and their derivatives
        etr = rhoair * (1.-fwet) * delta * lai/(rb(i)+rs) &
            * (qsatl(i) - qaf(botlay))

        IF (botlay == 2) THEN
           etr_dtl = rhoair * (1.-fwet) * delta * lai/(rb(i)+rs) &
                   * (1.-wtlq0(i)/facq)*qsatldT(i)
        ELSE
           etr_dtl = rhoair * (1.-fwet) * delta * lai/(rb(i)+rs) &
                   * (1.-wtaq0(1)*wtgq0(2)*wtlq0(i)/facq-wtlq0(i))*qsatldT(i)
        ENDIF

        IF (etr.ge.etrc) THEN
           etr = etrc
           etr_dtl = 0.
        ENDIF

        evplwet = rhoair * (1.-delta*(1.-fwet)) * lsai/rb(i) &
                * (qsatl(i) - qaf(botlay))

        IF (botlay == 2) THEN
           evplwet_dtl = rhoair * (1.-delta*(1.-fwet)) * lsai/rb(i) &
                       * (1.-wtlq0(i)/facq)*qsatldT(i)
        ELSE
           evplwet_dtl = rhoair * (1.-delta*(1.-fwet)) * lsai/rb(i) &
                       * (1.-wtaq0(1)*wtgq0(2)*wtlq0(i)/facq-wtlq0(i))*qsatldT(i)
        ENDIF

        IF (evplwet.ge.ldew/deltim) THEN
           evplwet = ldew/deltim
           evplwet_dtl = 0.
        ENDIF

        fevpl = etr + evplwet
        fevpl_dtl = etr_dtl + evplwet_dtl

        erre = 0.
        fevpl_noadj = fevpl
        IF ( fevpl*fevpl_bef < 0. ) THEN
           erre  = -0.9*fevpl
           fevpl =  0.1*fevpl
        ENDIF

!-----------------------------------------------------------------------
! difference of temperatures by quasi-newton-raphson method for the non-linear system equations
!-----------------------------------------------------------------------

        ! 计算irab, dirab_dtl
        B(5)    = B_5*tl**4
        B1(5)   = B1_5*tl**4
        dBdT(5) = dBdT_5*tl**3
        X  = matmul(Ainv, B)
        ! dBdT前5项为0, dBdT*(0,0,0,0,0,1)
        dX = matmul(Ainv, dBdT*avec)

        ! 每步温度迭代进行计算, 最后一次应为tlbef
        irab = ( (sum(X(1:4)*VegVF(1:4)) + frl*VegVF(5))*ev - B1(5))/fcover(5)*fg
        dirab_dtl = ( sum(dX(1:4)*VegVF(1:4))*ev - dBdT(5) )/fcover(5)*fg

        ! 迭代叶片温度变化
        dtl(it) = (sabv + irab - fsenl - hvap*fevpl) &
           / (lsai*clai/deltim - dirab_dtl + fsenl_dtl + hvap*fevpl_dtl)
        dtl_noadj = dtl(it)

        ! check magnitude of change in leaf temperature limit to maximum allowed value

        IF (it .le. itmax) THEN

           ! put brakes on large temperature excursions
           IF (abs(dtl(it)).gt.delmax) THEN
              dtl(it) = delmax*dtl(it)/abs(dtl(it))
           ENDIF

           IF ((it.ge.2) .and. (dtl(it-1)*dtl(it).le.0.)) THEN
              dtl(it) = 0.5*(dtl(it-1) + dtl(it))
           ENDIF

        ENDIF

        tl = tlbef + dtl(it)
        tu(3) = tl

!-----------------------------------------------------------------------
! square roots differences of temperatures and fluxes for use as the condition of convergences
!-----------------------------------------------------------------------

        del  = sqrt( dtl(it)*dtl(it) )
        dele = dtl(it) * dtl(it) * &
           ( dirab_dtl**2 + fsenl_dtl**2 + hvap*fevpl_dtl**2 )
        dele = sqrt(dele)

!-----------------------------------------------------------------------
!  saturated vapor pressures and canopy air temperature, canopy air humidity
!-----------------------------------------------------------------------
! Recalculate leaf saturated vapor pressure (ei_)for updated leaf temperature
! and adjust specific humidity (qsatl_) proportionately
        CALL qsadv(tu(i),psrf,ei(i),deiDT(i),qsatl(i),qsatldT(i))

! update vegetation/ground surface temperature, canopy air temperature,
! canopy air humidity

        ! calculate wtll, wtlql
        wtll(:)  = 0.
        wtlql(:) = 0.

        DO i = 0, nurb
           clev = canlev(i)
           wtll(clev)  =  wtll(clev) + wtl0(i)*tu(i)
           wtlql(clev) = wtlql(clev) + wtlq0(i)*qsatl(i)
        ENDDO

        IF (numlay .eq. 2) THEN

           ! - Equations:
           ! taf(3) = wta0(3)*thm    + wtg0(3)*taf(2) + wtll(3)
           ! taf(2) = wta0(2)*taf(3) + wtg0(2)*tg     + wtll(2)
           !
           ! qaf(3) = wtaq0(3)*qm     + wtgq0(3)*qaf(2) + wtlql(3)
           ! qaf(2) = wtaq0(2)*qaf(3) + wtgq0(2)*qg     + wtlql(2)

           tmpw1  = wta0(3)*thm + wtll(3) &
                  + wtshi(3)*1/3.*(Fhac+Fwst)/rhoair/cpair
           fact   = 1. - wta0(2)*wtg0(3)
           taf(2) = (wta0(2)*tmpw1 + wtg0(2)*tg + wtll(2) + &
                     wtshi(2)*(2/3.*(Fhac+Fwst)+Fach)/rhoair/cpair) / fact

           tmpw1  = wtaq0(3)*qm + wtlql(3)
           facq   = 1. - wtaq0(2)*wtgq0(3)
           qaf(2) = (wtaq0(2)*tmpw1 + wtgq0(2)*qg + wtlql(2)) / facq

           qaf(3) = wtaq0(3)*qm + wtgq0(3)*qaf(2) + wtlql(3)
           taf(3) = wta0(3)*thm + wtg0 (3)*taf(2) + wtll (3)
           taf(3) = taf(3) + wtshi(3)*1/3.*(Fhac+Fwst)/rhoair/cpair

        ENDIF

        IF (numlay .eq. 3) THEN

           ! - Equations:
           ! taf(3) = wta0(3)*thm    + wtg0(3)*taf(2) + wtll(3)
           ! taf(2) = wta0(2)*taf(3) + wtg0(2)*taf(1) + wtll(2)
           ! taf(1) = wta0(1)*taf(2) + wtg0(1)*tg     + wtll(1)
           !
           ! qaf(3) = wtaq0(3)*qm     + wtgq0(3)*qaf(2) + wtlql(3)
           ! qaf(2) = wtaq0(2)*qaf(3) + wtaq0(2)*qaf(1) + wtlql(2)
           ! qaf(1) = wtaq0(1)*qaf(2) + wtaq0(1)*qg     + wtlql(1)

           tmpw1  = wta0(3)*thm + wtll(3) &
                  + wtshi(3)*1/3.*(Fhac+Fwst)/rhoair/cpair
           tmpw2  = wtg0(1)*tg  + wtll(1) &
                  + wtshi(1)*1/3.*(Fhac+Fwst)/rhoair/cpair
           fact   = 1. - wta0(2)*wtg0(3) - wtg0(2)*wta0(1)
           taf(2) = (wta0(2)*tmpw1 + wtg0(2)*tmpw2 + wtll(2) + &
                     wtshi(2)*(1/3.*(Fhac+Fwst)+Fach)/rhoair/cpair) / fact

           tmpw1  = wtaq0(3)*qm + wtlql(3)
           tmpw2  = wtgq0(1)*qg + wtlql(1)
           facq   = 1. - wtaq0(2)*wtgq0(3) - wtgq0(2)*wtaq0(1)
           qaf(2) = (wtaq0(2)*tmpw1 + wtgq0(2)*tmpw2 + wtlql(2)) / facq

           qaf(1) = wtaq0(1)*qaf(2) + wtgq0(1)*qg + wtlql(1)
           taf(1) = wta0 (1)*taf(2) + wtg0 (1)*tg + wtll (1) &
                  + wtshi(1)*1/3.*(Fhac+Fwst)/rhoair/cpair

           qaf(3) = wtaq0(3)*qm + wtgq0(3)*qaf(2) + wtlql(3)
           taf(3) = wta0 (3)*thm + wtg0(3)*taf(2) + wtll (3) &
                  + wtshi(3)*1/3.*(Fhac+Fwst)/rhoair/cpair

        ENDIF

! update co2 partial pressure within canopy air
        ! 05/02/2016: may have some problem with gdh2o, however,
        ! this variable seems never used here. Different height
        ! level vegetation should have different gdh2o, i.e.,
        ! different rd(layer) values.
        gah2o = 1.0/raw * tprcor/thm                     !mol m-2 s-1
        gdh2o = 1.0/rd(botlay) * tprcor/thm              !mol m-2 s-1

        pco2a = pco2m - 1.37*psrf/max(0.446,gah2o) * &
           (assim - respc - rsoil)

!-----------------------------------------------------------------------
! Update monin-obukhov length and wind speed including the stability effect
!-----------------------------------------------------------------------

        ! 这里使用的是最高层的taf和qaf
        ! 如何进行限制?是不是梯度太大的问题?运行单点模型测试
        dth = thm - taf(3)
        dqh =  qm - qaf(3)

        tstar = vonkar/(fh-fht)*dth
        qstar = vonkar/(fq-fqt)*dqh

        thvstar = tstar + 0.61*th*qstar
        zeta = zldis*vonkar*grav*thvstar / (ustar**2*thv)
        IF (zeta .ge. 0.) THEN                             !stable
           zeta = min(2.,max(zeta,1.e-6))
        ELSE                                             !unstable
           zeta = max(-100.,min(zeta,-1.e-6))
        ENDIF
        obu = zldis/zeta

        IF (zeta .ge. 0.) THEN
           um = max(ur,.1)
        ELSE
           wc = (-grav*ustar*thvstar*zii/thv)**(1./3.)
           wc2 = beta*beta*(wc*wc)
           um = sqrt(ur*ur+wc2)
        ENDIF

        IF (obuold*obu .lt. 0.) nmozsgn = nmozsgn+1
        IF (nmozsgn .ge. 4) obu = zldis/(-0.01)
        obuold = obu

!-----------------------------------------------------------------------
! Test for convergence
!-----------------------------------------------------------------------

        it = it+1

        IF (it .gt. itmin) THEN
           fevpl_bef = fevpl
           det = max(del,del2)
           dee = max(dele,dele2)
           IF (det .lt. dtmin .and. dee .lt. dlemin) EXIT
        ENDIF

     ENDDO

! ======================================================================
!     END stability iteration
! ======================================================================

     z0m = z0mu
     zol = zeta
     rib = min(5.,zol*ustar**2/(vonkar**2/fh*um**2))

! canopy fluxes and total assimilation amd respiration

     IF (lai .gt. 0.001) THEN
        rst = rs/lai
     ELSE
        rs    = 2.0e4
        assim = 0.
        respc = 0.
        rst   = 2.0e4
     ENDIF
     respc = respc + rsoil

! canopy fluxes and total assimilation amd respiration

     fsenl = fsenl + fsenl_dtl*dtl(it-1) &
        ! add the imbalanced energy below due to T adjustment to sensibel heat
        + (dtl_noadj-dtl(it-1)) * (lsai*clai/deltim - dirab_dtl &
        + fsenl_dtl + hvap*fevpl_dtl) &
        ! add the imbalanced energy below due to q adjustment to sensibel heat
        + hvap*erre

     etr     = etr     +     etr_dtl*dtl(it-1)
     evplwet = evplwet + evplwet_dtl*dtl(it-1)
     fevpl   = fevpl_noadj
     fevpl   = fevpl   +   fevpl_dtl*dtl(it-1)

     elwmax  = ldew/deltim

     elwdif  = max(0., evplwet-elwmax)
     evplwet = min(evplwet, elwmax)

     fevpl = fevpl - elwdif
     fsenl = fsenl + hvap*elwdif

!-----------------------------------------------------------------------
! Update dew accumulation (kg/m2)
!-----------------------------------------------------------------------

     ldew = max(0., ldew-evplwet*deltim)

!-----------------------------------------------------------------------
! balance check
!-----------------------------------------------------------------------

     err = sabv + irab + dirab_dtl*dtl(it-1) &
         - fsenl - hvap*fevpl

#if(defined CLMDEBUG)
     IF (abs(err) .gt. .2) &
        write(6,*) 'energy imbalance in UrbanVegFlux.F90', &
        i,it-1,err,sabv,irab,fsenl,hvap*fevpl
#endif

!-----------------------------------------------------------------------
! 植被温度变化计算长波辐射改变
! 包括墙壁、地面吸收和向上长波辐射
!-----------------------------------------------------------------------

     ! 各组分长波辐射吸收值
     lwsun = ( ewall*X(1) - B1(1) ) / (1-ewall)
     lwsha = ( ewall*X(2) - B1(2) ) / (1-ewall)
     lgimp = ( egimp*X(3) - B1(3) ) / (1-egimp)
     lgper = ( egper*X(4) - B1(4) ) / (1-egper)
     lveg  = ( (sum(X(1:4)*VegVF(1:4)) + frl*VegVF(5))*ev - B1(5) )
     lout  = sum( X * SkyVF )

     ! +因叶片温度变化，各组分长波辐射吸收值
     lwsun = lwsun + ( ewall*dX(1) ) / (1-ewall) * dtl(it-1)
     lwsha = lwsha + ( ewall*dX(2) ) / (1-ewall) * dtl(it-1)
     lgimp = lgimp + ( egimp*dX(3) ) / (1-egimp) * dtl(it-1)
     lgper = lgper + ( egper*dX(4) ) / (1-egper) * dtl(it-1)
     lveg  = lveg  + ( sum(dX(1:4)*VegVF(1:4))*ev - dBdT(5) ) * dtl(it-1)
     lout  = lout  + sum( dX * SkyVF * dtl(it-1) )

     ! Energy balance check
     err = lwsun + lwsha + lgimp + lgper + lveg + lout

     IF (abs(err-frl) > 1e-6) THEN
        print *, "Longwave - Energy Balance Check error!", err-frl
     ENDIF

     ! 计算单位面积
     IF (fcover(1) > 0.) lwsun = lwsun / fcover(1) * fg !/ (4*fwsun*HL*fb/fg)
     IF (fcover(2) > 0.) lwsha = lwsha / fcover(2) * fg !/ (4*fwsha*HL*fb/fg)
     IF (fcover(3) > 0.) lgimp = lgimp / fcover(3) * fg !/ fgimp
     IF (fcover(4) > 0.) lgper = lgper / fcover(4) * fg !/ fgper
     IF (fcover(5) > 0.) lveg  = lveg  / fcover(5) * fg !/ fv/fg

     ! 加上上次余量
     lwsun = lwsun + lwsun_bef
     lwsha = lwsha + lwsha_bef
     lgimp = lgimp + lgimp_bef
     lgper = lgper + lgper_bef
     lveg  = lveg  + lveg_bef

     tafu = taf(2)

!-----------------------------------------------------------------------
! wind stresses
!-----------------------------------------------------------------------

     taux = - rhoair*us/ram
     tauy = - rhoair*vs/ram

!-----------------------------------------------------------------------
! fluxes from roof, walls to canopy space
!-----------------------------------------------------------------------

     ! sensible heat fluxes
     fsenroof = rhoair*cpair*cfh(0)*(troof-taf(3))
     fsenwsun = rhoair*cpair*cfh(1)*(twsun-taf(2))
     fsenwsha = rhoair*cpair*cfh(2)*(twsha-taf(2))

     ! latent heat fluxes
     fevproof = rhoair*cfw(0)*(qsatl(0)-qaf(3))

     croofs = rhoair*cpair*cfh(0)*(1.-wtg0(3)*wta0(2)*wtl0(0)/fact-wtl0(0))
     cwalls = rhoair*cpair*cfh(1)*(1.-wtl0(1)/fact)
     croofl = rhoair*cfw(0)*(1.-wtgq0(3)*wtaq0(2)*wtlq0(0)/facq-wtlq0(0))*qsatldT(0)

     croof = croofs + croofl*htvp_roof

!-----------------------------------------------------------------------
! 计算城市地面各组分的感热、潜热
!-----------------------------------------------------------------------

     fsengimp = cpair*rhoair*cgh(botlay)*(tgimp-taf(botlay))
     fsengper = cpair*rhoair*cgh(botlay)*(tgper-taf(botlay))

     fevpgimp = rhoair*cgw(botlay)*(qgimp-qaf(botlay))
     fevpgper = rhoair*cgw(botlay)*(qgper-qaf(botlay))

!-----------------------------------------------------------------------
! Derivative of soil energy flux with respect to soil temperature
!-----------------------------------------------------------------------

     IF (botlay == 2) THEN
        cgrnds = cpair*rhoair*cgh(2)*(1.-wtg0(2)/fact)
        cgimpl = rhoair*cgw(2)*(1.-wtgq0(2)/facq)*dqgimpdT
        cgperl = rhoair*cgw(2)*(1.-wtgq0(2)/facq)*dqgperdT
     ELSE !botlay == 1
        cgrnds = cpair*rhoair*cgh(1)*(1.-wta0(1)*wtg0(2)*wtg0(1)/fact-wtg0(1))
        cgimpl = rhoair*cgw(1)*(1.-wtaq0(1)*wtgq0(2)*wtgq0(1)/facq-wtgq0(1))*dqgimpdT
        cgperl = rhoair*cgw(1)*(1.-wtaq0(1)*wtgq0(2)*wtgq0(1)/facq-wtgq0(1))*dqgperdT
     ENDIF

     cgimp = cgrnds + cgimpl*htvp_gimp
     cgper = cgrnds + cgperl*htvp_gper

!-----------------------------------------------------------------------
! 2 m height air temperature
!-----------------------------------------------------------------------

     tref = thm + vonkar/(fh-fht)*dth * (fh2m/vonkar - fh/vonkar)
     qref =  qm + vonkar/(fq-fqt)*dqh * (fq2m/vonkar - fq/vonkar)

  END SUBROUTINE UrbanVegFlux
!----------------------------------------------------------------------


  SUBROUTINE dewfraction (sigf,lai,sai,dewmx,ldew,fwet)

!=======================================================================
! Original author: Yongjiu Dai, September 15, 1999
!
! determine fraction of foliage covered by water and
! fraction of foliage that is dry and transpiring
!
!=======================================================================

     USE precision
     IMPLICIT NONE

     REAL(r8), intent(in) :: sigf   ! fraction of veg cover, excluding snow-covered veg [-]
     REAL(r8), intent(in) :: lai    ! leaf area index  [-]
     REAL(r8), intent(in) :: sai    ! stem area index  [-]
     REAL(r8), intent(in) :: dewmx  ! maximum allowed dew [0.1 mm]
     REAL(r8), intent(in) :: ldew   ! depth of water on foliage [kg/m2/s]

     REAL(r8), intent(out) :: fwet  ! fraction of foliage covered by water [-]

     REAL(r8) lsai                  ! lai + sai
     REAL(r8) dewmxi                ! inverse of maximum allowed dew [1/mm]
     REAL(r8) vegt                  ! sigf*lsai
!
!-----------------------------------------------------------------------
! Fwet is the fraction of all vegetation surfaces which are wet
! including stem area which contribute to evaporation
     lsai = lai + sai
     dewmxi = 1.0/dewmx
      ! why * sigf? may have bugs
      ! 06/17/2018:
      ! for ONLY one PFT, there may be no problem
      ! but for multiple PFTs, bugs exist!!!
      ! convert the whole area ldew to sigf ldew
     vegt   =  lsai

     fwet = 0
     IF (ldew > 0.) THEN
        fwet = ((dewmxi/vegt)*ldew)**.666666666666

! Check for maximum limit of fwet
        fwet = min(fwet,1.0)

     ENDIF

  END SUBROUTINE dewfraction
!----------------------------------------------------------------------

  REAL(r8) FUNCTION uprofile(utop, fc, bee, alpha, z0mg, htop, hbot, z)

     USE precision
     USE FRICTION_VELOCITY
     IMPLICIT NONE

     REAL(r8), intent(in) :: utop
     REAL(r8), intent(in) :: fc
     REAL(r8), intent(in) :: bee
     REAL(r8), intent(in) :: alpha
     REAL(r8), intent(in) :: z0mg
     REAL(r8), intent(in) :: htop
     REAL(r8), intent(in) :: hbot
     REAL(r8), intent(in) :: z

     REAL(r8) :: ulog,uexp

     ! when canopy LAI->0, z0->zs, fac->1, u->umoninobuk
     ! canopy LAI->large, fac->0 or=0, u->log profile
     ulog = utop*log(z/z0mg)/log(htop/z0mg)
     uexp = utop*exp(-alpha*(1-(z-hbot)/(htop-hbot)))

     uprofile = bee*fc*min(uexp,ulog) + (1-bee*fc)*ulog

     RETURN
  END FUNCTION uprofile

  REAL(r8) FUNCTION kprofile(ktop, fc, bee, alpha, &
                    displah, htop, hbot, obu, ustar, z)

     USE precision
     USE FRICTION_VELOCITY
     IMPLICIT NONE

     REAL(r8), parameter :: com1 = 0.4
     REAL(r8), parameter :: com2 = 0.08

     REAL(r8), intent(in) :: ktop
     REAL(r8), intent(in) :: fc
     REAL(r8), intent(in) :: bee
     REAL(r8), intent(in) :: alpha
     REAL(r8), intent(in) :: displah
     REAL(r8), intent(in) :: htop
     REAL(r8), intent(in) :: hbot
     REAL(r8), intent(in) :: obu
     REAL(r8), intent(in) :: ustar
     REAL(r8), intent(in) :: z

     REAL(r8) :: kexp
     REAL(r8) :: klin, klins
     REAL(r8) :: kcob
     REAL(r8) :: fac

     klin = ktop*z/htop

     ! 02/07/2018: changed combination
     fac  = 1. / (1.+exp(-(displah-com1)/com2))
! 05/29/2021, yuan: bug. not initialized
     !TODO: 检查fac的设定，为什么设置为0
     fac  = 0.
     kcob = 1. / (fac/klin + (1.-fac)/kmoninobuk(0.,obu,ustar,z))

     kexp     = ktop*exp(-alpha*(htop-z)/(htop-hbot))
     kprofile = 1./( bee*fc/min(kexp,kcob) + (1-bee*fc)/kcob )

     RETURN

  END FUNCTION kprofile

  REAL(r8) FUNCTION uintegral(utop, fc, bee, alpha, z0mg, &
                    htop, hbot, ztop, zbot)

     USE precision
     IMPLICIT NONE

     REAL(r8), intent(in) :: utop
     REAL(r8), intent(in) :: fc
     REAL(r8), intent(in) :: bee
     REAL(r8), intent(in) :: alpha
     REAL(r8), intent(in) :: z0mg
     REAL(r8), intent(in) :: htop
     REAL(r8), intent(in) :: hbot
     REAL(r8), intent(in) :: ztop
     REAL(r8), intent(in) :: zbot

     INTEGER  :: i, n
     REAL(r8) :: dz, z, u

     ! 09/26/2017: change fixed n -> fixed dz
     dz = 0.001 ! fordebug
     n  = int( (ztop-zbot) / dz ) + 1

     uintegral = 0.

     DO i = 1, n
        IF (i < n) THEN
           z = ztop - (i-0.5)*dz
        ELSE
           dz = ztop - zbot - (n-1)*dz
           z  = zbot + 0.5*dz
        ENDIF

        u = uprofile(utop, fc, bee, alpha, z0mg, htop, hbot, z)

        u = max(0._r8, u)
        !uintegral = uintegral + sqrt(u)*dz / (htop-hbot)
! 03/04/2020, yuan: TODO-hard to solve
        ! u开根号后不能解析求解积分，可近似直接对u积分
        ! 如此，最后就不用平方
        uintegral = uintegral + u*dz / (ztop-zbot)
     ENDDO

     !uintegral = uintegral * uintegral

     RETURN

  END FUNCTION uintegral


  !TODO: 计算ztop到zbot之间的effective wind speed
  REAL(r8) FUNCTION ueffect(utop, htop, hbot, &
        ztop, zbot, z0mg, alpha, bee, fc)
     USE precision
     IMPLICIT NONE

     REAL(r8), intent(in) :: utop
     REAL(r8), intent(in) :: htop
     REAL(r8), intent(in) :: hbot
     REAL(r8), intent(in) :: ztop
     REAL(r8), intent(in) :: zbot
     REAL(r8), intent(in) :: z0mg
     REAL(r8), intent(in) :: alpha
     REAL(r8), intent(in) :: bee
     REAL(r8), intent(in) :: fc

     REAL(r8) :: roots(2), uint
     INTEGER  :: rootn

     rootn = 0
     uint  = 0.

     !二分法去找根，满足一定精度，假设最多2个根
     CALL ufindroots(ztop,zbot,(ztop+zbot)/2., &
        utop, htop, hbot, z0mg, alpha, roots, rootn)

! 03/10/2020, yuan: TODO-done, 编写函数fuint
     IF (rootn == 0) THEN ! no root
        uint = uint + fuint(utop, ztop, zbot, &
           htop, hbot, z0mg, alpha, bee, fc)
     ENDIF

     IF (rootn == 1) THEN
        uint = uint + fuint(utop, ztop, roots(1), &
           htop, hbot, z0mg, alpha, bee, fc)
        uint = uint + fuint(utop, roots(1), zbot, &
           htop, hbot, z0mg, alpha, bee, fc)
     ENDIF

     IF (rootn == 2) THEN
        uint = uint + fuint(utop, ztop,     roots(1), &
           htop, hbot, z0mg, alpha, bee, fc)
        uint = uint + fuint(utop, roots(1), roots(2), &
           htop, hbot, z0mg, alpha, bee, fc)
        uint = uint + fuint(utop, roots(2), zbot,     &
           htop, hbot, z0mg, alpha, bee, fc)
     ENDIF

     ueffect = uint / (ztop-zbot)

     RETURN

  END FUNCTION ueffect


  REAL(r8) FUNCTION fuint(utop, ztop, zbot, &
        htop, hbot, z0mg, alpha, bee, fc)

     USE precision
     IMPLICIT NONE

     REAL(r8), intent(in) :: utop, ztop, zbot
     REAL(r8), intent(in) :: htop, hbot
     REAL(r8), intent(in) :: z0mg, alpha
     REAL(r8), intent(in) :: bee, fc

     ! local variables
     REAL(r8) :: fuexpint, fulogint

     fulogint = utop/log(htop/z0mg) *&
        !(ztop*log(ztop/z0mg) - zbot*log(zbot/z0mg) + zbot - ztop) / (ztop-zbot)
        (ztop*log(ztop/z0mg) - zbot*log(zbot/z0mg) + zbot - ztop)

     IF (udif((ztop+zbot)/2.,utop,htop,hbot,z0mg,alpha) <= 0) THEN
        ! uexp is smaller
        fuexpint = utop*(htop-hbot)/alpha*( &
           ! yuan, 12/28/2020:
           exp(-alpha*(htop-ztop)/(htop-hbot))-&
           exp(-alpha*(htop-zbot)/(htop-hbot)) )
           ! yuan, 06/01/2021:
           !exp(-alpha*(ztop-ztop)/(htop-hbot))-&
           !exp(-alpha*(ztop-zbot)/(htop-hbot)) )

        fuint = bee*fc*fuexpint + (1.-bee*fc)*fulogint
     ELSE
        ! ulog is smaller
        fuint = fulogint
     ENDIF

     RETURN

  END FUNCTION fuint


  RECURSIVE SUBROUTINE ufindroots(ztop,zbot,zmid, &
     utop, htop, hbot, z0mg, alpha, roots, rootn)

     USE precision
     IMPLICIT NONE

     REAL(r8), intent(in) :: ztop, zbot, zmid
     REAL(r8), intent(in) :: utop, htop, hbot
     REAL(r8), intent(in) :: z0mg, alpha

     REAL(r8), intent(inout) :: roots(2)
     INTEGER,  intent(inout) :: rootn

     ! local variables
     REAL(r8) :: udif_ub, udif_lb

     udif_ub = udif(ztop,utop,htop,hbot,z0mg,alpha)
     udif_lb = udif(zmid,utop,htop,hbot,z0mg,alpha)

     IF (udif_ub*udif_lb == 0) THEN
        IF (udif_lb == 0) THEN !root found
           rootn = rootn + 1
           IF (rootn > 2) THEN
              print *, "U root number > 2, abort!"
              CALL abort
           ENDIF
           roots(rootn) = zmid
        ENDIF
     ELSE IF (udif_ub*udif_lb < 0) THEN
        IF (ztop-zmid < 0.01) THEN
           rootn = rootn + 1 !root found
           IF (rootn > 2) THEN
              print *, "U root number > 2, abort!"
              CALL abort
           ENDIF
           roots(rootn) = (ztop+zmid)/2.
        ELSE
           CALL ufindroots(ztop,zmid,(ztop+zmid)/2., &
              utop, htop, hbot, z0mg, alpha, roots, rootn)
        ENDIF
     ENDIF

     udif_ub = udif(zmid,utop,htop,hbot,z0mg,alpha)
     udif_lb = udif(zbot,utop,htop,hbot,z0mg,alpha)

     IF (udif_ub*udif_lb == 0) THEN
        IF (udif_ub == 0) THEN !root found
           rootn = rootn + 1
           IF (rootn > 2) THEN
              print *, "U root number > 2, abort!"
              CALL abort
           ENDIF
           roots(rootn) = zmid
        ENDIF
     ELSE IF (udif_ub*udif_lb < 0) THEN
        IF (zmid-zbot < 0.01) THEN
           rootn = rootn + 1 !root found
           IF (rootn > 2) THEN
              print *, "U root number > 2, abort!"
              CALL abort
           ENDIF
           roots(rootn) = (zmid+zbot)/2.
        ELSE
           CALL ufindroots(zmid,zbot,(zmid+zbot)/2., &
              utop, htop, hbot, z0mg, alpha, roots, rootn)
        ENDIF
     ENDIF

  END SUBROUTINE ufindroots


  REAL(r8) FUNCTION udif(z, utop, htop, hbot, z0mg, alpha)

     USE precision
     IMPLICIT NONE

     REAL(r8), intent(in) :: z, utop, htop, hbot
     REAL(r8), intent(in) :: z0mg, alpha

     REAL(r8) :: uexp, ulog

     ! yuan, 12/28/2020:
     !uexp = utop*exp(-alpha*(1-(z-hbot)/(htop-hbot)))
     uexp = utop*exp(-alpha*(htop-z)/(htop-hbot))
     ulog = utop*log(z/z0mg)/log(htop/z0mg)

     udif = uexp - ulog

     RETURN

  END FUNCTION udif


  ! 03/08/2020, yuan: TODO-done, change it to analytical solution
  REAL(r8) FUNCTION kintegral(ktop, fc, bee, alpha, z0mg, &
        displah, htop, hbot, obu, ustar, ztop, zbot)
     USE precision
     IMPLICIT NONE

     REAL(r8), intent(in) :: ktop
     REAL(r8), intent(in) :: fc
     REAL(r8), intent(in) :: bee
     REAL(r8), intent(in) :: alpha
     REAL(r8), intent(in) :: z0mg
     REAL(r8), intent(in) :: displah
     REAL(r8), intent(in) :: htop
     REAL(r8), intent(in) :: hbot
     REAL(r8), intent(in) :: obu
     REAL(r8), intent(in) :: ustar
     REAL(r8), intent(in) :: ztop
     REAL(r8), intent(in) :: zbot

     INTEGER  :: i, n
     REAL(r8) :: dz, z, k

     kintegral = 0.

     IF (ztop <= zbot) THEN
        RETURN
     ENDIF

     ! 09/26/2017: change fixed n -> fixed dz
     ! 10/05/2017: need to improve
     dz = 0.001 ! fordebug
     n  = int( (ztop-zbot) / dz ) + 1

     DO i = 1, n
        IF (i < n) THEN
           z  = ztop - (i-0.5)*dz
        ELSE
           dz = ztop - zbot - (n-1)*dz
           z  = zbot + 0.5*dz
        ENDIF

        k = kprofile(ktop, fc, bee, alpha, &
           displah, htop, hbot, obu, ustar, z)

        kintegral = kintegral + 1./k * dz

     ENDDO

     RETURN

  END FUNCTION kintegral


  REAL(r8) FUNCTION frd(ktop, htop, hbot, &
        ztop, zbot, displah, z0h, obu, ustar, &
        z0mg, alpha, bee, fc)

     USE precision
     IMPLICIT NONE

     REAL(r8), intent(in) :: ktop, htop, hbot
     REAL(r8), intent(in) :: ztop, zbot
     REAL(r8), intent(in) :: displah, z0h, obu, ustar
     REAL(r8), intent(in) :: z0mg, alpha, bee, fc

     ! local parameters
     REAL(r8), parameter :: com1 = 0.4
     REAL(r8), parameter :: com2 = 0.08

     REAL(r8) :: roots(2), fac, kint
     INTEGER  :: rootn

     rootn = 0
     kint  = 0.

     ! calculate fac
     ! yuan, 12/28/2020:
     !fac = 1. / (1.+exp(-(displah-com1)/com2))
! 05/29/2021, yuan: bug. not initialized
     !TODO: 检查fac的设定，为什么设定为0
     fac = 0.
     roots(:) = 0.

     CALL kfindroots(ztop,zbot,(ztop+zbot)/2., &
        ktop, htop, hbot, obu, ustar, fac, alpha, roots, rootn)

     !print *, roots, rootn
     IF (rootn == 0) THEN !no root
        kint = kint + fkint(ktop, ztop, zbot, htop, hbot, &
           z0h, obu, ustar, fac, alpha, bee, fc)
     ENDIF

     IF (rootn == 1) THEN
        kint = kint + fkint(ktop, ztop, roots(1), htop, hbot, &
           z0h, obu, ustar, fac, alpha, bee, fc)
        kint = kint + fkint(ktop, roots(1), zbot, htop, hbot, &
           z0h, obu, ustar, fac, alpha, bee, fc)
     ENDIF

     IF (rootn == 2) THEN
        kint = kint + fkint(ktop, ztop, roots(1), htop, hbot, &
           z0h, obu, ustar, fac, alpha, bee, fc)
        kint = kint + fkint(ktop, roots(1), roots(2), htop, hbot, &
           z0h, obu, ustar, fac, alpha, bee, fc)
        kint = kint + fkint(ktop, roots(2), zbot, htop, hbot, &
           z0h, obu, ustar, fac, alpha, bee, fc)
     ENDIF

     frd = kint

     RETURN

  END FUNCTION frd


  REAL(r8) FUNCTION fkint(ktop, ztop, zbot, htop, hbot, &
        z0h, obu, ustar, fac, alpha, bee, fc)

     USE precision
     USE FRICTION_VELOCITY
     IMPLICIT NONE

     REAL(r8), intent(in) :: ktop, ztop, zbot
     REAL(r8), intent(in) :: htop, hbot
     REAL(r8), intent(in) :: z0h, obu, ustar, fac, alpha
     REAL(r8), intent(in) :: bee, fc

     ! local variables
     REAL(r8) :: fkexpint, fkcobint

     !klin = ktop*z/htop
     !kcob = 1./(fac/klin + (1.-fac)/kmoninobuk(0.,obu,ustar,z))
     fkcobint = fac*htop/ktop*(log(ztop)-log(zbot)) +&
        (1.-fac)*kintmoninobuk(0.,z0h,obu,ustar,ztop,zbot)

     IF (kdif((ztop+zbot)/2.,ktop,htop,hbot,obu,ustar,fac,alpha) <= 0) THEN
        ! kexp is smaller
        IF (alpha > 0) THEN
           fkexpint = -(htop-hbot)/alpha/ktop*( &
              exp(alpha*(htop-ztop)/(htop-hbot))-&
              exp(alpha*(htop-zbot)/(htop-hbot)) )
        ELSE
           fkexpint = (ztop-zbot)/ktop
        ENDIF

        fkint = bee*fc*fkexpint + (1.-bee*fc)*fkcobint
     ELSE
        ! kcob is smaller
        fkint = fkcobint
     ENDIF

     RETURN
  END FUNCTION fkint


  RECURSIVE SUBROUTINE kfindroots(ztop,zbot,zmid, &
     ktop, htop, hbot, obu, ustar, fac, alpha, roots, rootn)

     USE precision
     IMPLICIT NONE

     REAL(r8), intent(in) :: ztop, zbot, zmid
     REAL(r8), intent(in) :: ktop, htop, hbot
     REAL(r8), intent(in) :: obu, ustar, fac, alpha

     REAL(r8), intent(inout) :: roots(2)
     INTEGER,  intent(inout) :: rootn

     ! local variables
     REAL(r8) :: kdif_ub, kdif_lb

     !print *, "*** CALL recursive SUBROUTINE kfindroots!!"
     kdif_ub = kdif(ztop,ktop,htop,hbot,obu,ustar,fac,alpha)
     kdif_lb = kdif(zmid,ktop,htop,hbot,obu,ustar,fac,alpha)

     IF (kdif_ub*kdif_lb == 0) THEN
        IF (kdif_lb == 0) THEN !root found
           rootn = rootn + 1
           IF (rootn > 2) THEN
              print *, "K root number > 2, abort!"
              CALL abort
           ENDIF
           roots(rootn) = zmid
        ENDIF
     ELSE IF (kdif_ub*kdif_lb < 0) THEN
        IF (ztop-zmid < 0.01) THEN
           rootn = rootn + 1 !root found
           IF (rootn > 2) THEN
              print *, "K root number > 2, abort!"
              CALL abort
           ENDIF
           roots(rootn) = (ztop+zmid)/2.
        ELSE
           CALL kfindroots(ztop,zmid,(ztop+zmid)/2., &
              ktop, htop, hbot, obu, ustar, fac, alpha, roots, rootn)
        ENDIF
     ENDIF

     kdif_ub = kdif(zmid,ktop,htop,hbot,obu,ustar,fac,alpha)
     kdif_lb = kdif(zbot,ktop,htop,hbot,obu,ustar,fac,alpha)

     IF (kdif_ub*kdif_lb == 0) THEN
        IF (kdif_ub == 0) THEN !root found
           rootn = rootn + 1
           IF (rootn > 2) THEN
              print *, "K root number > 2, abort!"
              CALL abort
           ENDIF
           roots(rootn) = zmid
        ENDIF
     ELSE IF (kdif_ub*kdif_lb < 0) THEN
        IF (zmid-zbot < 0.01) THEN
           rootn = rootn + 1 !root found
           IF (rootn > 2) THEN
              print *, "K root number > 2, abort!"
              CALL abort
           ENDIF
           roots(rootn) = (zmid+zbot)/2.
        ELSE
           CALL kfindroots(zmid,zbot,(zmid+zbot)/2., &
              ktop, htop, hbot, obu, ustar, fac, alpha, roots, rootn)
        ENDIF
     ENDIF

  END SUBROUTINE kfindroots


  REAL(r8) FUNCTION kdif(z, ktop, htop, hbot, &
        obu, ustar, fac, alpha)

     USE precision
     USE FRICTION_VELOCITY
     IMPLICIT NONE

     REAL(r8), intent(in) :: z, ktop, htop, hbot
     REAL(r8), intent(in) :: obu, ustar, fac, alpha

     REAL(r8) :: kexp, klin, kcob

     ! yuan, 12/28/2020:
     !kexp = ktop*exp(-alpha*(1-(z-hbot)/(htop-hbot)))
     kexp = ktop*exp(-alpha*(htop-z)/(htop-hbot))

     klin = ktop*z/htop
     kcob = 1./(fac/klin + (1.-fac)/kmoninobuk(0.,obu,ustar,z))

     kdif = kexp - kcob

     RETURN

  END FUNCTION kdif


  SUBROUTINE cal_z0_displa (lai, h, fc, z0, displa)

     USE PhysicalConstants, only: vonkar
     IMPLICIT NONE

     REAL(r8), intent(in)  :: lai
     REAL(r8), intent(in)  :: h
     REAL(r8), intent(in)  :: fc
     REAL(r8), intent(out) :: z0
     REAL(r8), intent(out) :: displa

     REAL(r8), parameter :: Cd   = 0.2   !leaf drag coefficient
     REAL(r8), parameter :: cd1  = 7.5   !a free parameter for d/h calculation, Raupach 1992, 1994
     REAL(r8), parameter :: psih = 0.193 !psih = ln(cw) - 1 + cw^-1, cw = 2, Raupach 1994

     ! local variables
     REAL(r8) :: fai, sqrtdragc, temp1, delta , lai0

     ! when assume z0=0.01, displa=0
     ! to calculate lai0, delta displa
     !----------------------------------------------------
     sqrtdragc = -vonkar/(log(0.01/h) - psih)
     sqrtdragc = max(sqrtdragc, 0.0031**0.5)
     IF (sqrtdragc .le. 0.3) THEN
        fai = (sqrtdragc**2-0.003) / 0.3
        fai = min(fai, fc*(1-exp(-20.)))
     ELSE
        fai = 0.29
        print *, "z0m, displa error!"
     ENDIF

     ! calculate delta displa when z0 = 0.01
     lai0  = -log(1.-fai/fc)/0.5
     temp1 = (2.*cd1*fai)**0.5
     delta = -h * ( fc*1.1*log(1. + (Cd*lai0*fc)**0.25) + &
        (1.-fc)*(1.-(1.-exp(-temp1))/temp1) )

     ! calculate z0m, displa
     !----------------------------------------------------
     ! NOTE: potential bug below, ONLY apply for spheric
     ! crowns. For other cases, fc*(...) ==> a*fc*(...)
     fai   = fc*(1. - exp(-0.5*lai))
     sqrtdragc = min( (0.003+0.3*fai)**0.5, 0.3 )
     temp1 = (2.*cd1*fai)**0.5

     IF (lai > lai0) THEN
        displa = delta + h*( &
           (  fc)*1.1*log(1. + (Cd*lai*fc)**0.25) + &
           (1-fc)*(1.-(1.-exp(-temp1))/temp1) )
     ELSE
        displa = h*( &
           (  fc)*1.1*log(1. + (Cd*lai*fc)**0.25) + &
           (1-fc)*(1.-(1.-exp(-temp1))/temp1) )
     ENDIF

     displa = max(displa, 0.)
     z0 = (h-displa) * exp(-vonkar/sqrtdragc + psih)

     IF (z0 < 0.01) THEN
        z0 = 0.01
        displa = 0.
     ENDIF

  END SUBROUTINE cal_z0_displa

END MODULE UrbanFlux
