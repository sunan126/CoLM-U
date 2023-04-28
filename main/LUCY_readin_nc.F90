#include <define.h>

SUBROUTINE LUCY_readin_nc(dir_srfdata)

      USE precision
      USE GlobalVars
      USE MOD_TimeVariables
      USE MOD_TimeInvariants
      USE MOD_UrbanTimeInvars
      USE ncio
      USE omp_lib

      IMPLICIT NONE

      CHARACTER(LEN=256), INTENT(in) :: dir_srfdata

      CHARACTER(LEN=256) :: lndname

      INTEGER :: ncid
      INTEGER :: reg_vid, pop_vid, veh_vid, wed_vid, weh_vid, wdh_vid, met_vid, hol_vid
      INTEGER :: i, j, t, u, npatch, r

      REAL(r8), allocatable :: lreg_id(:,:)
      REAL(r8), allocatable :: lpopcell(:,:,:)
      REAL(r8):: lweek_holiday(231,7)    ,&  !weekday and weekendday
                 lvehc_prof   (231,24,2) ,&  !diurnal traffic profile
                 lhum_prof    (231,24 )  ,&  !diurnal metabolize profile
                 lfix_holiday (231,365)  ,&  !public holiday
                 lvehicle     (231,3)        !number of cars/mobike/freight


      ! READ in inputdata for LUCY
      lndname = trim("/hard/dongwz/github/AHE/LUCY_AHE.nc")
      !print*,trim(lndname)
      CALL nccheck( nf90_open(trim(lndname), nf90_nowrite, ncid) )

      allocate ( lpopcell(1:lon_points,1:lat_points, 3) )
      allocate ( lreg_id (1:lon_points,1:lat_points   ) )

      CALL nccheck( nf90_inq_varid(ncid, "COUNTRY"    , reg_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "POPCELL"    , pop_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "vehicle"    , veh_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "weekendday" , wed_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "weekendhour", weh_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "weekdayhour", wdh_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "metabolism" , met_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "holiday"    , hol_vid) )


      CALL nccheck( nf90_get_var(ncid, reg_vid, lreg_id ) )
      CALL nccheck( nf90_get_var(ncid, pop_vid, lpopcell) )
      CALL nccheck( nf90_get_var(ncid, veh_vid, lvehicle) )
      CALL nccheck( nf90_get_var(ncid, wed_vid, lweek_holiday) )
      CALL nccheck( nf90_get_var(ncid, weh_vid, lvehc_prof(:,:,2)) )
      CALL nccheck( nf90_get_var(ncid, wdh_vid, lvehc_prof(:,:,1)) )
      CALL nccheck( nf90_get_var(ncid, met_vid, lhum_prof) )
      CALL nccheck( nf90_get_var(ncid, hol_vid, lfix_holiday) )


#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j,t,u,r,npatch)
#endif
      DO u = 1, numurban

         npatch = urb2patch(u)
         i = patch2lon(npatch)
         j = patch2lat(npatch)
         t = urbclass(u)

         reg_id (u)        = 70!lreg_id(i,j)
         r = reg_id(u)
         popcell(u)        = 13000!lpopcell(i,j,t)
         IF (r > 0) THEN
            vehicle     (u,:) = 669!lvehicle(r,:)
            week_holiday(u,:) = lweek_holiday(r,:)
            weh_prof    (u,:) = lvehc_prof   (r,:,2)
            wdh_prof    (u,:) = lvehc_prof   (r,:,1)
            hum_prof    (u,:) = lhum_prof    (r,:)
            fix_holiday (u,:) = lfix_holiday (r,:)
         ENDIF

      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      deallocate ( lreg_id )
      deallocate ( lpopcell)

      CALL nccheck( nf90_close(ncid) )

END SUBROUTINE LUCY_readin_nc
