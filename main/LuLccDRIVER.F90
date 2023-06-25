#include <define.h>

SUBROUTINE LuLccDRIVER (casename,dir_srfdata,dir_restart,&
                        nam_srfdata,nam_urbdata,idate,greenwich,&
                        dir_rawdata,edgen,edgee,edges,edgew)

!=======================================================================
! PURPOSE:
!   the main subroutine for Land use and land cover change simulation
!
! Created by Hua Yuan, 04/08/2022
!=======================================================================

   USE precision
   USE MOD_LuLccTimeInvars
   USE MOD_LuLccTimeVars
   USE MOD_LuLccTransferMatrix

   IMPLICIT NONE

   CHARACTER(LEN=256), intent(in) :: casename    !casename name
   CHARACTER(LEN=256), intent(in) :: dir_srfdata !surface data directory
   CHARACTER(LEN=256), intent(in) :: dir_restart !case restart data directory
   CHARACTER(LEN=256), intent(in) :: nam_srfdata !surface data filename
   CHARACTER(LEN=256), intent(in) :: nam_urbdata !urban data filename

   LOGICAL, intent(in)    :: greenwich           !true: greenwich time, false: local time
   INTEGER, intent(inout) :: idate(3)            !year, julian day, seconds of the starting time

   CHARACTER(len=256), intent(in) :: dir_rawdata
   REAL(r8), intent(in) :: edgen                 !northern edge of grid (degrees)
   REAL(r8), intent(in) :: edgee                 !eastern edge of grid (degrees)
   REAL(r8), intent(in) :: edges                 !southern edge of grid (degrees)
   REAL(r8), intent(in) :: edgew                 !western edge of grid (degrees)

   ! allocate LuLcc memory
   CALL allocate_LuLccTimeInvars
   CALL allocate_LuLccTimeVars

   ! SAVE variables
   CALL SAVE_LuLccTimeInvars
   CALL SAVE_LuLccTimeVars

   ! cold start for LuLcc
   print *, ">>> LULCC: initializing..."
   CALL LuLccInitialize (casename,dir_srfdata,dir_restart,&
                         nam_srfdata,nam_urbdata,idate,greenwich)

   !TODO: to be an option using namelist
   ! ! simple method for variable recovery
   ! print *, ">>> LULCC: simple method for variable recovery..."
   ! CALL REST_LuLccTimeVars

   ! conserved method for variable revocery
   print *, ">>> LULCC: Energy&Mass conserve for variable recovery..."
   CALL MAKE_LuLccTransferMatrix(casename,idate,dir_rawdata,dir_restart,edgen,edgee,edges,edgew)
   CALL READ_LuLccTransferMatrix(casename,idate,dir_restart)
   CALL LuLccEnergyMassConserve()

   ! deallocate LuLcc memory
   CALL deallocate_LuLccTimeInvars()
   CALL deallocate_LuLccTimeVars()

END SUBROUTINE LuLccDRIVER
