#include <define.h>

MODULE UrbanBEM

  USE precision
  USE GlobalVars
  USE PhysicalConstants
  USE UrbanShortwave, only: MatrixInverse

  IMPLICIT NONE
  SAVE
  PRIVATE

  ! A simple building energy model to calculate room temperature
  PUBLIC :: SimpleBEM

CONTAINS

  !-------------------------------------------------
  SUBROUTINE SimpleBEM ( deltim, rhoair, fcover, H, troom_max, troom_min, &
                         troof_nl, twsun_nl, twsha_nl, &
                         tkdz_roof, tkdz_wsun, tkdz_wsha, taf, &
                         troom, troof_inner, twsun_inner, twsha_inner, &
                         Fhac, Fwst, Fach )

     IMPLICIT NONE

     REAL(r8), intent(in) :: &
        deltim,     &! seconds in a time step [second]
        rhoair,     &! density air [kg/m3]
        fcover(0:2),&! fractional cover of roof, wall
        H,          &! average building height [m]
        troom_max,  &! maximum temperature of inner building
        troom_min,  &! minimum temperature of inner building
        troof_nl,   &! roof temperature at layer nl_roof
        twsun_nl,   &! sunlit wall temperature at layer nl_wall
        twsha_nl,   &! shaded wall temperature at layer nl_wall
        tkdz_roof,  &! temporal var for heat transfer of roof
        tkdz_wsun,  &! temporal var for heat transfer of sunlit wall
        tkdz_wsha,  &! temporal var for heat transfer of shaded wall
        taf          ! temperature of urban air

     REAL(r8), intent(inout) :: &
        troom        ! temperature of inner building

     REAL(r8), intent(out) :: &
        troof_inner,&! temperature of inner roof
        twsun_inner,&! temperature of inner sunlit wall
        twsha_inner,&! temperature of inner shaded wall
        Fhac,       &! flux from heat or cool AC
        Fwst,       &! waste heat from cool or heat
        Fach         ! flux from air exchange

     ! local variables
     REAL(r8) ::    &
        ACH,        &! air exchange coefficience
        hcv_roof,   &! convective exchange ceofficience for roof<->room
        hcv_wall,   &! convective exchange ceofficience for wall<->room
        waste_cool, &! waste heat for AC cooling
        waste_heat   ! waste heat for AC heating

     REAL(r8) ::    &
        f_wsun,     &! weight factor for sunlit wall
        f_wsha       ! weight factor for shaded wall

     REAL(r8) ::    &
        A(4,4),     &! Heat transfer matrix
        Ainv(4,4),  &! Inverse of Heat transfer matrix
        B(4),       &! B for Ax=B
        X(4)         ! x for Ax=B

  !=================================================================
  !
  ! o 求解以下联立方程组
  ! o 隐式求解troom, troof_inner, twsun_inner, twsha_innter
  !
  !    Hc_roof = Fn_roof        .................................(1)
  !    Hc_wsun = Fn_wsun        .................................(2)
  !    Hc_wsha = Fn_wsha        .................................(3)
  !
  !                   Troom' - Troom
  !    H*rhoair*cpair*-------------- =
  !                         dt
  !     ACH
  !    -----*H*rhoair*cpair*(Taf-Troom') + Hc_roof + Hc_wsun + Hc_wsha
  !    3600? dt?
  !                             .................................(4)
  !=================================================================

     ACH = 0.3          !air exchange coefficience
     hcv_roof   = 4.040 !convective exchange ceofficience for roof<->room
     hcv_wall   = 3.076 !convective exchange ceofficience for wall<->room
     waste_cool = 0.6   !waste heat for AC cooling
     waste_heat = 0.2   !waste heat for AC heating

     f_wsun = fcover(1)/fcover(0) !weight factor for sunlit wall
     f_wsha = fcover(2)/fcover(0) !weight factor for shaded wall

     ! initialization
     Fhac = 0.; Fwst = 0.; Fach = 0.;

     ! Ax = B
     ! set values for heat transfer matrix
     ! 1: roof, 2: sunlit wall, 3: shaded wall, 4: room
     A(:,:) = 0.
     A(1,:) = (/0.5*hcv_roof+0.5*tkdz_roof, 0., 0., -0.5*hcv_roof/)
     A(2,:) = (/0., 0.5*hcv_wall+0.5*tkdz_wsun, 0., -0.5*hcv_wall/)
     A(3,:) = (/0., 0., 0.5*hcv_wall+0.5*tkdz_wsha, -0.5*hcv_wall/)

     A(4,:) = (/-0.5*hcv_roof, -0.5*hcv_wall*f_wsun, -0.5*hcv_wall*f_wsha, &
                 0.5*hcv_roof + 0.5*hcv_wall*f_wsun + 0.5*hcv_wall*f_wsha +&
                 H*rhoair*cpair/deltim + (ACH/3600.)*H*rhoair*cpair /)

     B(1) = -0.5*hcv_roof*(troof_inner-troom) + 0.5*tkdz_roof*(troof_nl-troof_inner) + 0.5*tkdz_roof*troof_nl
     B(2) = -0.5*hcv_wall*(twsun_inner-troom) + 0.5*tkdz_wsun*(twsun_nl-twsun_inner) + 0.5*tkdz_wsun*twsun_nl
     B(3) = -0.5*hcv_wall*(twsha_inner-troom) + 0.5*tkdz_wsha*(twsha_nl-twsha_inner) + 0.5*tkdz_wsha*twsha_nl

     B(4) = H*rhoair*cpair*troom/deltim + (ACH/3600.)*H*rhoair*cpair*taf &
          + 0.5*hcv_roof*(troof_inner-troom) &
          + 0.5*hcv_wall*(twsun_inner-troom)*f_wsun &
          + 0.5*hcv_wall*(twsha_inner-troom)*f_wsha

     ! Inverse of matrix A
     Ainv = MatrixInverse(A)

     ! Matrix computing to revole multiple reflections
     X = matmul(Ainv, B)

     troof_inner = X(1)
     twsun_inner = X(2)
     twsha_inner = X(3)
     troom       = X(4)

     Fach = (ACH/3600.)*H*rhoair*cpair*(troom - taf)

     IF (troom > troom_max) THEN !cooling case
        Fhac  = H*rhoair*cpair*(troom-troom_max)/deltim
        troom = troom_max
        Fwst  = Fhac*waste_cool
     ENDIF

     IF (troom < troom_min) THEN !heating case
        Fhac  = H*rhoair*cpair*(troom-troom_min)/deltim
        troom = troom_min
        Fwst  = abs(Fhac)*waste_heat
        ! nagative value, set it to 0.
        Fhac  = 0.
     ENDIF

     Fach = Fach*fcover(0)
     Fwst = Fwst*fcover(0)
     Fhac = Fhac*fcover(0)

  END SUBROUTINE

END MODULE UrbanBEM
