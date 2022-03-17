#include <define.h>

SUBROUTINE UrbanLAI_readin_nc (lon_points,lat_points,&
           month, dir_model_landdata)

! ===========================================================
! Read in urban LAI, SAI and urban tree cover data
! ===========================================================

      USE precision
      USE GlobalVars
      USE LC_Const
      USE MOD_TimeVariables
      USE MOD_TimeInvariants
      USE MOD_UrbanTimeInvars
      USE ncio
      USE omp_lib

      IMPLICIT NONE

      INTEGER, intent(in) :: lon_points
      INTEGER, intent(in) :: lat_points
      INTEGER, intent(in) :: month
      CHARACTER(LEN=256), intent(in) :: dir_model_landdata

      CHARACTER(LEN=256) :: lndname
      INTEGER :: ncid
      INTEGER :: urbantreelai_vid, urbantreesai_vid
      INTEGER :: i, j, t, u, npatch

      REAL(r8), allocatable :: urbantreelai(:,:,:)
      REAL(r8), allocatable :: urbantreesai(:,:,:)

! READ in Leaf area index and stem area index
      lndname = trim(dir_model_landdata)//'urban_0.5x0.5.MOD2000_v5.nc'
      print*,trim(lndname)
      CALL nccheck( nf90_open(trim(lndname), nf90_nowrite, ncid) )

      allocate ( urbantreelai(1:lon_points,1:lat_points,1:N_URB) )
      allocate ( urbantreesai(1:lon_points,1:lat_points,1:N_URB) )

      CALL nccheck( nf90_inq_varid(ncid, "URBAN_TREE_LAI", urbantreelai_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "URBAN_TREE_SAI", urbantreesai_vid) )

      CALL nccheck( nf90_get_var(ncid, urbantreelai_vid, urbantreelai, &
           start=(/1,1,1,month/), count=(/lon_points,lat_points,N_URB,1/)) )
      CALL nccheck( nf90_get_var(ncid, urbantreesai_vid, urbantreesai, &
           start=(/1,1,1,month/), count=(/lon_points,lat_points,N_URB,1/)) )

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j,t,u,npatch)
#endif
      DO u = 1, numurban

         npatch = urb2patch(u)
         i = patch2lon(npatch)
         j = patch2lat(npatch)
         t = urbclass(u)

         tlai(npatch)  = urbantreelai(i,j,t) !leaf area index
         tsai(npatch)  = urbantreesai(i,j,t) !stem are index
         green(npatch) = 1.                  !fraction of green leaf

      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      deallocate ( urbantreelai )
      deallocate ( urbantreesai )

      CALL nccheck( nf90_close(ncid) )

END SUBROUTINE UrbanLAI_readin_nc
