#include <define.h>

SUBROUTINE LuLccDRIVER (casename,dir_model_landdata,dir_restart_hist,&
                        idate,greenwich,lon_points,lat_points)

!=======================================================================
! PURPOSE:
!   the main subroutine for Land use and land cover change simulation
!
! Created by Hua Yuan, 04/08/2022
!=======================================================================

   USE precision
   USE MOD_LuLccTimeInvars
   USE MOD_LuLccTimeVars
   USE MOD_LuLccTMatrix

   IMPLICIT NONE

   CHARACTER(LEN=256), intent(in) :: casename           !casename name
   CHARACTER(LEN=256), intent(in) :: dir_model_landdata !surface data directory
   CHARACTER(LEN=256), intent(in) :: dir_restart_hist   !case restart data directory

   LOGICAL, intent(in)    :: greenwich   !true: greenwich time, false: local time
   INTEGER, intent(in)    :: lon_points  !number of longitude points on model grid
   INTEGER, intent(in)    :: lat_points  !number of latitude points on model grid
   INTEGER, intent(inout) :: idate(3)    !year, julian day, seconds of the starting time

   ! allocate LuLcc memory
   CALL allocate_LuLccTimeInvars (lon_points, lat_points)
   CALL allocate_LuLccTimeVars (lon_points, lat_points)

   ! SAVE variables
   CALL SAVE_LuLccTimeInvars()
   CALL SAVE_LuLccTimeVars()

   ! cold start for LuLcc
   print *, ">>> LULCC: initializing..."
   CALL LuLccInitialize (casename,dir_model_landdata,dir_restart_hist,&
                         idate,greenwich,lon_points,lat_points)

   ! simple method for variable recovery
   print *, ">>> LULCC: simple method for variable recovery..."
   CALL REST_LuLccTimeVars (lon_points, lat_points)

   ! conserved method for variable revocery
   print *, ">>> LULCC: Mass&Energy conserve for variable recovery..."
   !CALL READ_LuLccTMatrix()
   !CALL LuLccEnergyConserve()
   !CALL LuLccWaterConserve()

   ! deallocate LuLcc memory
   CALL deallocate_LuLccTimeInvars()
   CALL deallocate_LuLccTimeVars()

END SUBROUTINE LuLccDRIVER
