#include <define.h>

SUBROUTINE LAI_readin_nc (year, month, dir_srfdata, nam_srfdata)

! ===========================================================
! Read in the LAI, the LAI dataset was created by Yuan et al. (2019)
! http://globalchange.bnu.edu.cn
! ===========================================================

      USE precision
      USE GlobalVars
      USE LC_Const
      USE PFT_Const
      USE MOD_TimeInvariants
      USE MOD_TimeVariables
      USE MOD_PFTimeVars
      USE MOD_PCTimeVars
      USE MOD_PFTimeInvars
      USE MOD_PCTimeInvars
      USE ncio
      USE omp_lib

      IMPLICIT NONE

      INTEGER, intent(in) :: year
      INTEGER, intent(in) :: month
      CHARACTER(LEN=256), intent(in) :: dir_srfdata
      CHARACTER(LEN=256), intent(in) :: nam_srfdata

      CHARACTER(LEN=256) :: c
      CHARACTER(LEN=256) :: lndname
      CHARACTER(len=256) :: cyear
      INTEGER :: ncid
      INTEGER :: glai_vid,  gsai_vid
      INTEGER :: lclai_vid, lcsai_vid, pftlai_vid, pftsai_vid
      INTEGER :: pclai_vid, pcsai_vid, pctpc_vid
      INTEGER :: i, j, m, n, t, p, ps, pe, pc, npatch

      REAL(r8), allocatable :: glai     (:,:)
      REAL(r8), allocatable :: gsai     (:,:)
      REAL(r8), allocatable :: lclai  (:,:,:)
      REAL(r8), allocatable :: lcsai  (:,:,:)
      REAL(r8), allocatable :: pftlai (:,:,:)
      REAL(r8), allocatable :: pftsai (:,:,:)
      REAL(r8), allocatable :: pclai(:,:,:,:)
      REAL(r8), allocatable :: pcsai(:,:,:,:)
      REAL(r8), allocatable :: pctpc(:,:,:,:)

! READ in Leaf area index and stem area index

      write(cyear,'(i4.4)') year
      lndname = trim(dir_srfdata)//trim(cyear)//'/'//trim(nam_srfdata)
      print*,trim(lndname)
      CALL nccheck( nf90_open(trim(lndname), nf90_nowrite, ncid) )

#ifdef IGBP_CLASSIFICATION
      allocate ( glai (1:lon_points,1:lat_points) )
      allocate ( gsai (1:lon_points,1:lat_points) )
      allocate ( lclai(1:lon_points,1:lat_points,1:N_land_classification) )
      allocate ( lcsai(1:lon_points,1:lat_points,1:N_land_classification) )

      CALL nccheck( nf90_inq_varid(ncid, "MONTHLY_LAI",    glai_vid ) )
      CALL nccheck( nf90_inq_varid(ncid, "MONTHLY_SAI",    gsai_vid ) )
      CALL nccheck( nf90_inq_varid(ncid, "MONTHLY_LC_LAI", lclai_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "MONTHLY_LC_SAI", lcsai_vid) )

      CALL nccheck( nf90_get_var(ncid, glai_vid, glai, &
                    start=(/1,1,month/), &
                    count=(/lon_points,lat_points,1/)) )
      CALL nccheck( nf90_get_var(ncid, gsai_vid, gsai, &
                    start=(/1,1,month/), &
                    count=(/lon_points,lat_points,1/)) )
      CALL nccheck( nf90_get_var(ncid, lclai_vid, lclai, &
                    start=(/1,1,1,month/), &
                    count=(/lon_points,lat_points,N_land_classification,1/)) )
      CALL nccheck( nf90_get_var(ncid, lcsai_vid, lcsai, &
                    start=(/1,1,1,month/), &
                    count=(/lon_points,lat_points,N_land_classification,1/)) )

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
         IF( m == 0 )THEN  !ocean
             fveg(npatch)  = 0.
             tlai(npatch)  = 0.
             tsai(npatch)  = 0.
             green(npatch) = 0.
         ELSE
             fveg(npatch)  = fveg0(m)           !fraction of veg. cover 每种IGBP有植被裸土比例参数
             IF (fveg0(m) > 0) THEN
! 05/19/2022, yuan: in case of land cover TYPE change
#if(defined USE_POINT_DATA)
                tlai(npatch)  = glai(i,j)             !leaf area index 单点模拟会被赋予网格LAI
                tsai(npatch)  = gsai(i,j)             !stem area index
#else
                tlai(npatch)  = lclai(i,j,m)/fveg0(m) !leaf area index  !把IGBP里的植被的裸地分开，LAI完全给植被
                tsai(npatch)  = lcsai(i,j,m)/fveg0(m) !stem area index
#endif
                green(npatch) = 1.                    !fraction of green leaf
             ELSE
                tlai(npatch)  = 0.       !leaf area index
                tsai(npatch)  = 0.       !stem area index
                green(npatch) = 0.       !fraction of green leaf
             ENDIF
         ENDIF
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      deallocate ( glai  )
      deallocate ( gsai  )
      deallocate ( lclai )
      deallocate ( lcsai )

#endif

#ifdef PFT_CLASSIFICATION
      allocate ( pftlai(1:lon_points,1:lat_points,0:N_PFT-1) )
      allocate ( pftsai(1:lon_points,1:lat_points,0:N_PFT-1) )
      allocate ( lclai (1:lon_points,1:lat_points,1:N_land_classification) )
      allocate ( lcsai (1:lon_points,1:lat_points,1:N_land_classification) )

      CALL nccheck( nf90_inq_varid(ncid, "MONTHLY_PFT_LAI", pftlai_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "MONTHLY_PFT_SAI", pftsai_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "MONTHLY_LC_LAI",  lclai_vid ) )
      CALL nccheck( nf90_inq_varid(ncid, "MONTHLY_LC_SAI",  lcsai_vid ) )

      CALL nccheck( nf90_get_var(ncid, pftlai_vid, pftlai, &
                    start=(/1,1,1,month/), &
                    count=(/lon_points,lat_points,N_PFT,1/)) )
      CALL nccheck( nf90_get_var(ncid, pftsai_vid, pftsai, &
                    start=(/1,1,1,month/), &
                    count=(/lon_points,lat_points,N_PFT,1/)) )
      CALL nccheck( nf90_get_var(ncid, lclai_vid, lclai, &
                    start=(/1,1,1,month/), &
                    count=(/lon_points,lat_points,N_land_classification,1/)) )
      CALL nccheck( nf90_get_var(ncid, lcsai_vid, lcsai, &
                    start=(/1,1,1,month/), &
                    count=(/lon_points,lat_points,N_land_classification,1/)) )

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j,m,n,t,p,ps,pe)
#endif
      DO npatch = 1, numpatch

         i = patch2lon(npatch)
         j = patch2lat(npatch)
         t = patchtype(npatch)
         m = patchclass(npatch)

! 12/07/2021, yuan: Urban filter
#ifdef URBAN_MODEL
         IF (m == URBAN) cycle
#endif
         IF (t == 0) THEN
            ps = patch_pft_s(npatch)
            pe = patch_pft_e(npatch)

            DO p = ps, pe
               n = pftclass(p)
               tlai_p(p) = pftlai(i,j,n)
               tsai_p(p) = pftsai(i,j,n)
            ENDDO

            tlai(npatch) = sum(tlai_p(ps:pe) * pftfrac(ps:pe))
            tsai(npatch) = sum(tsai_p(ps:pe) * pftfrac(ps:pe))

         ELSE
! 05/19/2022, yuan: change to lc lai/sai
            tlai(npatch) = lclai(i,j,m)/fveg0(m) !leaf area index
            tsai(npatch) = lcsai(i,j,m)/fveg0(m) !stem area index
         ENDIF

         fveg(npatch)  = fveg0(m)
         green(npatch) = 1.

      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      deallocate ( pftlai )
      deallocate ( pftsai )
      deallocate ( lclai  )
      deallocate ( lcsai  )

#endif

#ifdef PC_CLASSIFICATION
      allocate ( pclai(1:lon_points,1:lat_points, &
                          0:N_PFT-1,1:N_land_classification) )
      allocate ( pcsai(1:lon_points,1:lat_points, &
                          0:N_PFT-1,1:N_land_classification) )
      allocate ( pctpc(1:lon_points,1:lat_points, &
                          0:N_PFT-1,1:N_land_classification) )
      CALL nccheck( nf90_inq_varid(ncid, "MONTHLY_PC_LAI", pclai_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "MONTHLY_PC_SAI", pcsai_vid) )
      CALL nccheck( nf90_inq_varid(ncid, "PCT_PC",         pctpc_vid) )

      CALL nccheck( nf90_get_var(ncid, pclai_vid, pclai, &
                    start=(/1,1,1,1,month/), &
                    count=(/lon_points,lat_points,N_PFT,N_land_classification,1/)) )
      CALL nccheck( nf90_get_var(ncid, pcsai_vid, pcsai, &
                    start=(/1,1,1,1,month/), &
                    count=(/lon_points,lat_points,N_PFT,N_land_classification,1/)) )
      CALL nccheck( nf90_get_var(ncid, pctpc_vid, pctpc ) )

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j,t,m,pc)
#endif
      DO npatch = 1, numpatch

         i = patch2lon(npatch)
         j = patch2lat(npatch)
         t = patchtype(npatch)
         m = patchclass(npatch)

! 12/07/2021, yuan: Urban filter
#ifdef URBAN_MODEL
         IF (m == URBAN) cycle
#endif
         IF (t == 0) THEN
            pc = patch2pc(npatch)
            tlai_c(:,pc) = pclai(i,j,:,m)
            tsai_c(:,pc) = pcsai(i,j,:,m)
! 07/06/2022, yuan: LAICHANGE bug
            tlai(npatch) = sum(tlai_c(:,pc)*pcfrac(:,pc))
            tsai(npatch) = sum(tsai_c(:,pc)*pcfrac(:,pc))
         ELSE
! 12/28/2019, yuan: Bug
            ! pctpc from 1-100% -> 0-1
            tlai(npatch) = sum(pclai(i,j,:,m)*pctpc(i,j,:,m)/100.)
            tsai(npatch) = sum(pcsai(i,j,:,m)*pctpc(i,j,:,m)/100.)
         ENDIF

         fveg(npatch)  = fveg0(m)
         green(npatch) = 1.

      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      deallocate ( pclai )
      deallocate ( pcsai )
      deallocate ( pctpc )

#endif

      CALL nccheck( nf90_close(ncid) )

END SUBROUTINE LAI_readin_nc
