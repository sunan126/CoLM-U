#include <define.h>

MODULE MOD_UrbanTimeVars

! -------------------------------
! Created by Hua Yuan, 12/2020
! -------------------------------

   USE precision
   IMPLICIT NONE
   SAVE
! -----------------------------------------------------------------
! Time-varying state variables which reaquired by restart run

   REAL(r8), allocatable :: fwsun          (:) !sunlit fraction of walls [-]
   REAL(r8), allocatable :: dfwsun         (:) !change of sunlit fraction of walls [-]

   ! shortwave absorption
   REAL(r8), allocatable :: sroof      (:,:,:) !roof aborption [-]
   REAL(r8), allocatable :: swsun      (:,:,:) !sunlit wall absorption [-]
   REAL(r8), allocatable :: swsha      (:,:,:) !shaded wall absorption [-]
   REAL(r8), allocatable :: sgimp      (:,:,:) !impervious absorptioin [-]
   REAL(r8), allocatable :: sgper      (:,:,:) !pervious absorptioin [-]
   REAL(r8), allocatable :: slake      (:,:,:) !urban lake absorptioin [-]

   ! net longwave radiation for last time temperature change
   REAL(r8), allocatable :: lwsun          (:) !net longwave of sunlit wall [W/m2]
   REAL(r8), allocatable :: lwsha          (:) !net longwave of shaded wall [W/m2]
   REAL(r8), allocatable :: lgimp          (:) !net longwave of impervious  [W/m2]
   REAL(r8), allocatable :: lgper          (:) !net longwave of pervious [W/m2]
   REAL(r8), allocatable :: lveg           (:) !net longwave of vegetation [W/m2]

   REAL(r8), allocatable :: z_sno_roof   (:,:) !node depth of roof [m]
   REAL(r8), allocatable :: z_sno_gimp   (:,:) !node depth of impervious [m]
   REAL(r8), allocatable :: z_sno_gper   (:,:) !node depth pervious [m]
   REAL(r8), allocatable :: z_sno_lake   (:,:) !node depth lake [m]

   REAL(r8), allocatable :: dz_sno_roof  (:,:) !interface depth of roof [m]
   REAL(r8), allocatable :: dz_sno_gimp  (:,:) !interface depth of impervious [m]
   REAL(r8), allocatable :: dz_sno_gper  (:,:) !interface depth pervious [m]
   REAL(r8), allocatable :: dz_sno_lake  (:,:) !interface depth lake [m]

   REAL(r8), allocatable :: troof_inner    (:) !temperature of roof [K]
   REAL(r8), allocatable :: twsun_inner    (:) !temperature of sunlit wall [K]
   REAL(r8), allocatable :: twsha_inner    (:) !temperature of shaded wall [K]

   REAL(r8), allocatable :: t_roofsno    (:,:) !temperature of roof [K]
   REAL(r8), allocatable :: t_wallsun    (:,:) !temperature of sunlit wall [K]
   REAL(r8), allocatable :: t_wallsha    (:,:) !temperature of shaded wall [K]
   REAL(r8), allocatable :: t_gimpsno    (:,:) !temperature of impervious [K]
   REAL(r8), allocatable :: t_gpersno    (:,:) !temperature of pervious [K]
   REAL(r8), allocatable :: t_lakesno    (:,:) !temperature of pervious [K]

   REAL(r8), allocatable :: wliq_roofsno (:,:) !liquid water in layers [kg/m2]
   REAL(r8), allocatable :: wliq_gimpsno (:,:) !liquid water in layers [kg/m2]
   REAL(r8), allocatable :: wliq_gpersno (:,:) !liquid water in layers [kg/m2]
   REAL(r8), allocatable :: wliq_lakesno (:,:) !liquid water in layers [kg/m2]
   REAL(r8), allocatable :: wice_roofsno (:,:) !ice lens in layers [kg/m2]
   REAL(r8), allocatable :: wice_gimpsno (:,:) !ice lens in layers [kg/m2]
   REAL(r8), allocatable :: wice_gpersno (:,:) !ice lens in layers [kg/m2]
   REAL(r8), allocatable :: wice_lakesno (:,:) !ice lens in layers [kg/m2]

   REAL(r8), allocatable :: sag_roof       (:) !roof snow age [-]
   REAL(r8), allocatable :: sag_gimp       (:) !impervious ground snow age [-]
   REAL(r8), allocatable :: sag_gper       (:) !pervious ground snow age [-]
   REAL(r8), allocatable :: sag_lake       (:) !urban lake snow age [-]

   REAL(r8), allocatable :: scv_roof       (:) !roof snow cover [-]
   REAL(r8), allocatable :: scv_gimp       (:) !impervious ground snow cover [-]
   REAL(r8), allocatable :: scv_gper       (:) !pervious ground snow cover [-]
   REAL(r8), allocatable :: scv_lake       (:) !urban lake snow cover [-]

   REAL(r8), allocatable :: fsno_roof      (:) !roof snow fraction [-]
   REAL(r8), allocatable :: fsno_gimp      (:) !impervious ground snow fraction [-]
   REAL(r8), allocatable :: fsno_gper      (:) !pervious ground snow fraction [-]
   REAL(r8), allocatable :: fsno_lake      (:) !urban lake snow fraction [-]

   REAL(r8), allocatable :: snowdp_roof    (:) !roof snow depth [m]
   REAL(r8), allocatable :: snowdp_gimp    (:) !impervious ground snow depth [m]
   REAL(r8), allocatable :: snowdp_gper    (:) !pervious ground snow depth [m]
   REAL(r8), allocatable :: snowdp_lake    (:) !urban lake snow depth [m]

   REAL(r8), allocatable :: t_room         (:) !temperature of inner building [K]
   REAL(r8), allocatable :: tafu           (:) !temperature of outer building [K]
   REAL(r8), allocatable :: tu2m           (:) !2 m urban air temperature [K]
   REAL(r8), allocatable :: qu2m           (:) !2 m urban air humidity [kg/kg]
   REAL(r8), allocatable :: Fhac           (:) !sensible flux from heat or cool AC [W/m2]
   REAL(r8), allocatable :: Fwst           (:) !waste heat flux from heat or cool AC [W/m2]
   REAL(r8), allocatable :: Fach           (:) !flux from inner and outter air exchange [W/m2]

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_UrbanTimeVars
   PUBLIC :: deallocate_UrbanTimeVars

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE allocate_UrbanTimeVars ()
! ------------------------------------------------------
! Allocates memory for CLM 1d [numurban] variables
! ------------------------------------------------------
      USE precision
      USE GlobalVars
      IMPLICIT NONE

      allocate (fwsun                         (numurban))
      allocate (dfwsun                        (numurban))

      allocate (sroof                     (2,2,numurban))
      allocate (swsun                     (2,2,numurban))
      allocate (swsha                     (2,2,numurban))
      allocate (sgimp                     (2,2,numurban))
      allocate (sgper                     (2,2,numurban))
      allocate (slake                     (2,2,numurban))

      allocate (lwsun                         (numurban))
      allocate (lwsha                         (numurban))
      allocate (lgimp                         (numurban))
      allocate (lgper                         (numurban))
      allocate (lveg                          (numurban))

      allocate (z_sno_roof         (maxsnl+1:0,numurban))
      allocate (z_sno_gimp         (maxsnl+1:0,numurban))
      allocate (z_sno_gper         (maxsnl+1:0,numurban))
      allocate (z_sno_lake         (maxsnl+1:0,numurban))

      allocate (dz_sno_roof        (maxsnl+1:0,numurban))
      allocate (dz_sno_gimp        (maxsnl+1:0,numurban))
      allocate (dz_sno_gper        (maxsnl+1:0,numurban))
      allocate (dz_sno_lake        (maxsnl+1:0,numurban))

      allocate (t_roofsno    (maxsnl+1:nl_roof,numurban))
      allocate (t_wallsun    (maxsnl+1:nl_wall,numurban))
      allocate (t_wallsha    (maxsnl+1:nl_wall,numurban))
      allocate (t_gimpsno    (maxsnl+1:nl_soil,numurban))
      allocate (t_gpersno    (maxsnl+1:nl_soil,numurban))
      allocate (t_lakesno    (maxsnl+1:nl_soil,numurban))

      allocate (troof_inner                   (numurban))
      allocate (twsun_inner                   (numurban))
      allocate (twsha_inner                   (numurban))

      allocate (wliq_roofsno (maxsnl+1:nl_roof,numurban))
      allocate (wice_roofsno (maxsnl+1:nl_roof,numurban))
      allocate (wliq_gimpsno (maxsnl+1:nl_soil,numurban))
      allocate (wice_gimpsno (maxsnl+1:nl_soil,numurban))
      allocate (wliq_gpersno (maxsnl+1:nl_soil,numurban))
      allocate (wice_gpersno (maxsnl+1:nl_soil,numurban))
      allocate (wliq_lakesno (maxsnl+1:nl_soil,numurban))
      allocate (wice_lakesno (maxsnl+1:nl_soil,numurban))

      allocate (sag_roof                      (numurban))
      allocate (sag_gimp                      (numurban))
      allocate (sag_gper                      (numurban))
      allocate (sag_lake                      (numurban))
      allocate (scv_roof                      (numurban))
      allocate (scv_gimp                      (numurban))
      allocate (scv_gper                      (numurban))
      allocate (scv_lake                      (numurban))
      allocate (fsno_roof                     (numurban))
      allocate (fsno_gimp                     (numurban))
      allocate (fsno_gper                     (numurban))
      allocate (fsno_lake                     (numurban))
      allocate (snowdp_roof                   (numurban))
      allocate (snowdp_gimp                   (numurban))
      allocate (snowdp_gper                   (numurban))
      allocate (snowdp_lake                   (numurban))

      allocate (t_room                        (numurban))
      allocate (tafu                          (numurban))
      allocate (tu2m                          (numurban))
      allocate (qu2m                          (numurban))
      allocate (Fhac                          (numurban))
      allocate (Fwst                          (numurban))
      allocate (Fach                          (numurban))

   END SUBROUTINE allocate_UrbanTimeVars

   SUBROUTINE deallocate_UrbanTimeVars

      deallocate (fwsun        )
      deallocate (dfwsun       )

      deallocate (sroof        )
      deallocate (swsun        )
      deallocate (swsha        )
      deallocate (sgimp        )
      deallocate (sgper        )
      deallocate (slake        )

      deallocate (lwsun        )
      deallocate (lwsha        )
      deallocate (lgimp        )
      deallocate (lgper        )
      deallocate (lveg         )

      deallocate (z_sno_roof   )
      deallocate (z_sno_gimp   )
      deallocate (z_sno_gper   )
      deallocate (z_sno_lake   )

      deallocate (dz_sno_roof  )
      deallocate (dz_sno_gimp  )
      deallocate (dz_sno_gper  )
      deallocate (dz_sno_lake  )

      deallocate (t_roofsno    )
      deallocate (t_wallsun    )
      deallocate (t_wallsha    )
      deallocate (t_gimpsno    )
      deallocate (t_gpersno    )
      deallocate (t_lakesno    )

      deallocate (troof_inner  )
      deallocate (twsun_inner  )
      deallocate (twsha_inner  )

      deallocate (wliq_roofsno )
      deallocate (wice_roofsno )
      deallocate (wliq_gimpsno )
      deallocate (wice_gimpsno )
      deallocate (wliq_gpersno )
      deallocate (wice_gpersno )
      deallocate (wliq_lakesno )
      deallocate (wice_lakesno )

      deallocate (sag_roof     )
      deallocate (sag_gimp     )
      deallocate (sag_gper     )
      deallocate (sag_lake     )
      deallocate (scv_roof     )
      deallocate (scv_gimp     )
      deallocate (scv_gper     )
      deallocate (scv_lake     )
      deallocate (fsno_roof    )
      deallocate (fsno_gimp    )
      deallocate (fsno_gper    )
      deallocate (fsno_lake    )
      deallocate (snowdp_roof  )
      deallocate (snowdp_gimp  )
      deallocate (snowdp_gper  )
      deallocate (snowdp_lake  )

      deallocate (t_room       )
      deallocate (tafu         )
      deallocate (tu2m         )
      deallocate (qu2m         )
      deallocate (Fhac         )
      deallocate (Fwst         )
      deallocate (Fach         )

   END SUBROUTINE deallocate_UrbanTimeVars

END MODULE MOD_UrbanTimeVars
! ---------- EOP ------------
