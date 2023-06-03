#include <define.h>

SUBROUTINE LCC_readin_nc ()

! ===========================================================
! Read in the LAI, the LAI dataset was created by Yuan et al. (2019)
! http://globalchange.bnu.edu.cn
! ===========================================================

      USE precision
      USE GlobalVars
      USE LC_Const
      ! USE PFT_Const
      ! USE MOD_TimeInvariants
      ! USE MOD_TimeVariables
      ! USE MOD_PFTimeVars
      ! USE MOD_PCTimeVars
      ! USE MOD_PFTimeInvars
      ! USE MOD_PCTimeInvars
      USE ncio
      USE omp_lib
      USE MOD_LuLccTimeVars

      IMPLICIT NONE

      ! INTEGER, intent(in) :: year
      ! INTEGER, intent(in) :: month
      ! CHARACTER(LEN=256), intent(in) :: dir_srfdata
      ! CHARACTER(LEN=256), intent(in) :: nam_srfdata

      ! CHARACTER(LEN=256) :: c
      CHARACTER(LEN=256) :: lndname
      ! CHARACTER(len=256) :: cyear
      INTEGER :: ncid
      INTEGER :: lcfr_vid, lcto_vid, dchg_vid
      INTEGER :: i, j, m, npatch

      INTEGER,  allocatable :: lcfrin(:,:,:)
      INTEGER,  allocatable :: lctoin(:,:,:)
      REAL(r8), allocatable :: dchgin(:,:,:)


! READ in Leaf area index and stem area index

      write(cyear,'(i4.4)') year
      ! lndname = trim(dir_srfdata)//trim(cyear)//'/'//trim(nam_srfdata)
      lndname = '/hard/linwy20/mklcc/global_0.5x0.5_2019-2020.MOD.nc'
      print*,trim(lndname)
      CALL nccheck( nf90_open(trim(lndname), nf90_nowrite, ncid) )

#ifdef IGBP_CLASSIFICATION
      allocate ( lcfrin (1:nnum, 1:lon_points, 1:lat_points) )
      allocate ( lctoin (1:nnum, 1:lon_points, 1:lat_points) )
      allocate ( dchgin (1:nnum, 1:lon_points, 1:lat_points) )

      CALL nccheck( nf90_inq_varid(ncid, "LCFR", lcfr_vid ) )
      CALL nccheck( nf90_inq_varid(ncid, "LCTO", lcto_vid ) )
      CALL nccheck( nf90_inq_varid(ncid, "DCHG", dchg_vid) )

      CALL nccheck( nf90_get_var(ncid, lcfr_vid, lcfrin, &
                    start=(/1,1,1/), &
                    count=(/nnum,lon_points,lat_points/)) )
      CALL nccheck( nf90_get_var(ncid, lcto_vid, lctoin, &
                    start=(/1,1,1/), &
                    count=(/nnum,lon_points,lat_points/)) )
      CALL nccheck( nf90_get_var(ncid, dchg_vid, dchgin, &
                    start=(/1,1,1/), &
                    count=(/nnum,lon_points,lat_points/)) )

      allocate ( lcfr (1:nnum, numpatch) )
      allocate ( lcto (1:nnum, numpatch) )
      allocate ( dchg (1:nnum, numpatch) )  

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j,m)
#endif
      DO npatch = 1, numpatch

         i = patch2lon(npatch)
         j = patch2lat(npatch)
         m = patchclass(npatch)

! 12/07/2021, yuan: Urban filter
#ifdef URBAN_MODEL
         IF (m == URBAN) cycle
#endif
         IF ( m == 0 )THEN
            lcfrin(:,npatch)  = -99.
            lctoin(:,npatch)  = -99.
            dchgin(:,npatch)  = -99.
         ELSE
            lcfr(:,npatch)  = lcfrin(:,i,j) 
            lcto(:,npatch)  = lctoin(:,i,j) 
            dchg(:,npatch)  = dchgin(:,i,j)           
         ENDIF
      
      ENDDO

#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      deallocate ( lcfrin )
      deallocate ( lctoin  )
      deallocate ( dchgin )

#endif

      CALL nccheck( nf90_close(ncid) )

END SUBROUTINE LCC_readin_nc