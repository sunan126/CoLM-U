#include <define.h>

MODULE MOD_LuLccTMatrix
! -------------------------------
! Created by Hua Yuan, 04/2022
! May be renamed as "MOD_LuLccTransferPair"
! -------------------------------

  USE precision
  USE GlobalVars
  IMPLICIT NONE
  SAVE
! -----------------------------------------------------------------
  INTEGER :: nlcc = 272 ! maximum pairs of possible lulcc
  INTEGER :: nigbp = 17
  INTEGER :: nnum ! number of lulcc pairs chosen
  INTEGER,  allocatable :: lcfr(:,:)
  INTEGER,  allocatable :: lcto(:,:)
  REAL(r8), allocatable :: dchg(:,:)
  REAL(r8), allocatable :: selfchg(:,:)

  !TODO: need coding below...

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_LuLccTMatrix
  PUBLIC :: deallocate_LuLccTMatrix
  PUBLIC :: READ_LuLccTMatrix

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_LuLccTMatrix
  ! --------------------------------------------------------------------
  ! Allocates memory for LuLcc time invariant variables
  ! --------------------------------------------------------------------

     USE precision
     USE GlobalVars
     IMPLICIT NONE
!TODO: need coding below...

  END SUBROUTINE allocate_LuLccTMatrix


  SUBROUTINE READ_LuLccTMatrix

      USE precision
      USE GlobalVars
      USE MOD_TimeInvariants
      USE MOD_PFTimeInvars
      USE MOD_PCTimeInvars
      USE MOD_UrbanTimeInvars
      USE MOD_LuLccTimeVars
      USE ncio

      IMPLICIT NONE
!TODO: need coding below...

      CHARACTER(LEN=256) :: lndname
      INTEGER :: ncid
      INTEGER :: lcfr_vid, lcto_vid, dchg_vid, selfchg_vid
      INTEGER :: i, j, m, npatch
      INTEGER :: k,numtmp
      REAL(r8):: tmp

      INTEGER,  allocatable :: lcfrin   (:,:,:)
      INTEGER,  allocatable :: lctoin   (:,:,:)
      REAL(r8), allocatable :: dchgin   (:,:,:)
      REAL(r8), allocatable :: selfchgin(:,:,:)

      lndname = '/hard/linwy20/mklcc/global_0.5x0.5_2004-2005.MOD.nc'
      print*,trim(lndname)
      CALL nccheck( nf90_open(trim(lndname), nf90_nowrite, ncid) )

#ifdef IGBP_CLASSIFICATION
      allocate ( lcfrin    (1:nlcc , 1:lon_points, 1:lat_points) )
      allocate ( lctoin    (1:nlcc , 1:lon_points, 1:lat_points) )
      allocate ( dchgin    (1:nlcc , 1:lon_points, 1:lat_points) )
      allocate ( selfchgin (1:nigbp, 1:lon_points, 1:lat_points) )

      CALL nccheck( nf90_inq_varid(ncid, "LCFR",    lcfr_vid    ) )
      CALL nccheck( nf90_inq_varid(ncid, "LCTO",    lcto_vid    ) )
      CALL nccheck( nf90_inq_varid(ncid, "DCHG",    dchg_vid    ) )
      CALL nccheck( nf90_inq_varid(ncid, "SELFCHG", selfchg_vid ) )

      CALL nccheck( nf90_get_var(ncid, lcfr_vid, lcfrin, &
                    start=(/1,1,1/), &
                    count=(/nlcc,lon_points,lat_points/)) )

      CALL nccheck( nf90_get_var(ncid, lcto_vid, lctoin, &
                    start=(/1,1,1/), &
                    count=(/nlcc,lon_points,lat_points/)) )

      CALL nccheck( nf90_get_var(ncid, dchg_vid, dchgin, &
                    start=(/1,1,1/), &
                    count=(/nlcc,lon_points,lat_points/)) )

      CALL nccheck( nf90_get_var(ncid, selfchg_vid, selfchgin, &
                    start=(/1,1,1/), &
                    count=(/nigbp,lon_points,lat_points/)) )


      ! 确定nnum
      nnum = 0
      DO j = 1, lat_points 
        DO i = 1, lon_points
         
          tmp = 0  ! 累积百分比; 先保存自己不变的百分比
          DO k=1,nigbp
            IF (selfchgin(k,i,j) .gt. 0) THEN  
              tmp = tmp + selfchgin(k,i,j)
            ENDIF
          ENDDO

          numtmp = 0  !  记录百分比达到阈值时的索引
          DO k=1,nlcc
            IF (dchgin(k,i,j) .gt. 0) THEN  
              ! print*,'dchg=', dchgin(k,i,j) 
              tmp = tmp + dchgin(k,i,j)
              numtmp = numtmp + 1
            ENDIF
            IF ( (tmp .gt. 99) .and. (nnum .lt. numtmp) ) THEN
              nnum = numtmp
              EXIT
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      print*,'>>>>>>>>>>>>>>>>> nnum=',nnum

      allocate ( lcfr    (1:nnum , numpatch) )
      allocate ( lcto    (1:nnum , numpatch) )
      allocate ( dchg    (1:nnum , numpatch) ) 
      allocate ( selfchg (1:nigbp, numpatch) )  


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
            lcfr    (:,npatch)  = -99
            lcto    (:,npatch)  = -99
            dchg    (:,npatch)  = -99.
            selfchg (:,npatch)  = -99.
         ELSE
            lcfr    (:,npatch)  = lcfrin    (:,i,j) 
            lcto    (:,npatch)  = lctoin    (:,i,j) 
            dchg    (:,npatch)  = dchgin    (:,i,j) 
            selfchg (:,npatch)  = selfchgin (:,i,j)        
         ENDIF
      
      ENDDO

#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      deallocate ( lcfrin )
      deallocate ( lctoin  )
      deallocate ( dchgin )
      deallocate ( selfchgin )

#endif

      CALL nccheck( nf90_close(ncid) )

  END SUBROUTINE READ_LuLccTMatrix


  SUBROUTINE deallocate_LuLccTMatrix
! --------------------------------------------------
! Deallocates memory for LuLcc time invariant variables
! --------------------------------------------------

!TODO: need coding below...


  END SUBROUTINE deallocate_LuLccTMatrix

END MODULE MOD_LuLccTMatrix
! ---------- EOP ------------
