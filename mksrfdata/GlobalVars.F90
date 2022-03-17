#include <define.h>

MODULE GlobalVars
!-------------------------------------------------------------------------------
! Purpose:
!   Define global variables for mksurdata dir only
!
! ________________
! REVISION HISTORY:
! 08/2019, Hua Yuan: initial code
!-------------------------------------------------------------------------------
   USE precision

   IMPLICIT NONE
   SAVE

   INTEGER, PARAMETER :: nlat30s = 21600 !30sec, 1km
   INTEGER, PARAMETER :: nlon30s = 43200 !30sec, 1km
   INTEGER, PARAMETER :: nlat15s = 43200 !15sec, 500m
   INTEGER, PARAMETER :: nlon15s = 86400 !15sec, 500m

#ifdef USGS_CLASSIFICATION
   INTEGER, PARAMETER :: nlat = nlat30s  !30sec, 1km
   INTEGER, PARAMETER :: nlon = nlon30s  !30sec, 1km
   ! GLCC USGS number of land cover category
   INTEGER, PARAMETER :: N_land_classification = 24

#else
   INTEGER, PARAMETER :: nlat = nlat15s  !15sec, 500m
   INTEGER, PARAMETER :: nlon = nlon15s  !15sec, 500m
   ! MODIS IGBP number of land cover category
   INTEGER, PARAMETER :: N_land_classification = 17
#endif

END MODULE GlobalVars
