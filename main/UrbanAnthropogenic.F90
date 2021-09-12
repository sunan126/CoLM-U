#include <define.h>

MODULE UrbanAnthropogenic

  USE precision
  USE GlobalVars
  USE PhysicalConstants
  
  IMPLICIT NONE
  SAVE
  PRIVATE 

  PUBLIC :: SimpleBEM   !a simple building energy model to calculate room temperature

CONTAINS

  !-------------------------------------------------
  SUBROUTINE SimpleBEM ( deltim, rhoair, fcover, H, troom_max, troom_min, &
                         fn_roof, fn_wsun, fn_wsha, taf, troom, Fhac, Fwst, Fach )

     IMPLICIT NONE
     
     REAL(r8), intent(in) :: &
        deltim,     &! 
        rhoair,     &!
        fcover(0:2),&! 
        H,          &! 
        troom_max,  &! 
        troom_min,  &! 
        fn_roof,    &! 
        fn_wsun,    &! 
        fn_wsha,    &! 
        taf          !

     REAL(r8), intent(inout) :: &
        troom,      &! 
        Fhac,       &! 
        Fwst,       &! 
        Fach         ! 

     ! local variables
     REAL(r8) ::    &
        ACH,        &!
        waste_cool, &!
        waste_heat   !
     
     REAL(r8) ::    &
        troom_,     &!
        f_wsun,     &! 
        f_wsha       ! 

  !=================================================================
  !
  !                Troom' - Troom   
  ! H*rhoair*cpair*-------------- =
  !                      dt          
  !    ACH
  !   ------*H*rhoair*cpair*(Taf-Troom) - Fn_roof - Fn_wsun - Fn_wsha
  !    3600? dt? 
  !=================================================================

     ACH = 0.3
     waste_cool = 0.6
     waste_heat = 0.2
     f_wsun = fcover(1)/fcover(0)
     f_wsha = fcover(2)/fcover(0)
     
     Fach   = (ACH/3600.)*H*rhoair*cpair*(taf-troom)
     troom_ = Fach - fn_roof - f_wsun*fn_wsun - f_wsha*fn_wsha
     troom_ = (troom_*deltim)/(H*rhoair*cpair) + troom
            
     IF (troom_ > troom_max) THEN !cooling case
        Fhac  = H*rhoair*cpair*(troom_-troom_max)
        troom = troom_
        Fwst  = Fhac*waste_cool
     ENDIF 

     IF (troom_ < troom_min) THEN !heating case
        Fhac  = H*rhoair*cpair*(troom_-troom_min)
        troom = troom_
        Fwst  = abs(Fhac)*waste_heat
     ENDIF 

     Fach = Fach*fcover(0)
     Fwst = Fwst*fcover(0)
     Fhac = Fhac*fcover(0)

  END SUBROUTINE 

END MODULE UrbanAnthropogenic
