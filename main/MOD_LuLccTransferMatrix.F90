#include <define.h>

MODULE MOD_LuLccTransferMatrix
! -------------------------------
! Created by Wanyi Lin and Hua Yuan, 05/2023
! May be renamed as "MOD_LuLccTransferPair"
! -------------------------------

  USE precision
  USE GlobalVars
  IMPLICIT NONE
  SAVE
! -----------------------------------------------------------------

  INTEGER :: nlc = 18
  !TODO@Wanyi: The below variable is double defined.
  REAL(r8), allocatable :: lccpct(:,:,:)

  CHARACTER(len=*), parameter :: DATASRC = "MOD"  !"ESA"
  CHARACTER(len=256) syear, eyear


! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_LuLccTransferMatrix
  PUBLIC :: deallocate_LuLccTransferMatrix
  PUBLIC :: READ_LuLccTransferMatrix
  PUBLIC :: MAKE_LuLccTransferMatrix

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_LuLccTransferMatrix
  ! --------------------------------------------------------------------
  ! Allocates memory for LuLcc time invariant variables
  ! --------------------------------------------------------------------

     USE precision
     USE GlobalVars
     IMPLICIT NONE
!TODO: need coding below...

  END SUBROUTINE allocate_LuLccTransferMatrix


  SUBROUTINE MAKE_LuLccTransferMatrix(casename,idate,dir_rawdata,dir_restart,edgen,edgee,edges,edgew)

     USE precision
     USE GlobalVars
     USE ncio

     IMPLICIT NONE


     CHARACTER(LEN=256), intent(in) :: casename
     CHARACTER(len=256), intent(in) :: dir_rawdata
     CHARACTER(len=256), intent(in) :: dir_restart

     INTEGER,  intent(in) :: idate(3)
     REAL(r8), intent(in) :: edgen      !northern edge of grid (degrees)
     REAL(r8), intent(in) :: edgee      !eastern edge of grid (degrees)
     REAL(r8), intent(in) :: edges      !southern edge of grid (degrees)
     REAL(r8), intent(in) :: edgew      !western edge of grid (degrees)

     CHARACTER(len=*), parameter :: Title   = "Land surface model input land cover change data"
     CHARACTER(len=*), parameter :: Authors = "Yongjiu Dai's group at Sun Yat-sen University"
     CHARACTER(len=*), parameter :: Address = "School of Atmospheric Sciences, Sun Yat-sen University, Zhuhai, China"
     CHARACTER(len=*), parameter :: Email   = "yuanh25@mail.sysu.edu.cn"

     INTEGER, parameter :: nxy = 1200

     ! define input variables
     REAL(r8), dimension(:,:)    , allocatable :: lcdatas,lcdatae

     ! define output variables
     REAL(r8), dimension(:,:)    , allocatable :: area
     REAL(r8), dimension(:,:,:,:), allocatable :: sumarea,lccpct


     ! define other variables
     ! -----------------------------------------------
     CHARACTER(len=256) lndname,lndnames,lndnamee
     CHARACTER(len=256) :: reg1, reg2, reg3, reg4

     INTEGER :: reg(4)
     INTEGER :: ncid

     ! variable ids
     INTEGER :: lon_vid, lat_vid, lc_vid
     INTEGER :: lat_dimid, lon_dimid, lcs_dimid, lce_dimid

     ! output variables/ids
     INTEGER :: varea,lcs_vid, lce_vid
     INTEGER :: vlccpct
     INTEGER :: lct(22)

     REAL(r8) :: deg2rad, re, dx, dy, wgt, sumpct
     REAL(r8) :: sarea(nxy,nxy)
     REAL(r8) :: lat, lon
     REAL(r8) :: dll, lone(nxy), lonw(nxy), latn(nxy), lats(nxy)
     REAL(r8) :: dllo, latso(lat_points), lonso(lon_points)

     INTEGER :: nlc
     INTEGER :: lcs,lce
     INTEGER :: i, j, io, jo, si, sj, ei, ej
     INTEGER :: reglat,reglon,reglon_,sreglat,sreglon,ereglat,ereglon

     INTEGER :: XY2D(2), GRID3d(3), LC3D(3), PFT3D(3), LC4D(4), PFT4D(4), PC4D(4), PC5D(5)
     LOGICAL :: fileExists

     print*, casename,dir_rawdata,dir_restart,edgen,edgee,edges,edgew


     deg2rad = pi/180.
     re = 6.37122e6 * 0.001

     write(syear,'(i4.4)') idate(1) - 1
     write(eyear,'(i4.4)') idate(1)

     IF (DATASRC == "MOD") THEN
        nlc = 18 !17
     ELSE
        nlc = 22
     ENDIF

     ! allocate memory
     print *, ">>> allocating memory..."
     allocate( lcdatas    (nxy, nxy) )
     allocate( lcdatae    (nxy, nxy) )

     allocate( sumarea    (lon_points,lat_points,nlc,nlc))
     allocate( lccpct     (lon_points,lat_points,nlc,nlc))
     allocate( area       (lon_points,lat_points        ))

     DO i = 1, 22
        lct(i) = i
     ENDDO

     ! initialization
     print *, ">>> initializing..."
     lccpct  (:,:,:,:) = 0.
     area        (:,:) = 0.

     !TODO@Wanyi: the below code for reading data need
     !rewrite for CoLM202X.

     ! calculate output grid size
     dllo = (edgen-edges)/lat_points

     ! calculate the coordinate variable data
     DO i = 1, lon_points
        lonso(i) = edgew + i*dllo - 0.5*dllo
        IF (lonso(i) > 180) THEN
           lonso(i) = lonso(i) - 360
        ENDIF
     ENDDO

     DO i = 1, lat_points
        latso(i) = edgen - i*dllo + 0.5*dllo
     ENDDO

     ! calculate input grid size
     dll = 5./nxy

     ! calculate start region latitude/longitude
     sreglat = 90 - int((90.-edgen)/5.)*5
     ereglat = 90 - int((90.-edges-0.5*dll)/5.)*5
     IF (ereglat > sreglat) ereglat = sreglat

     sreglon = -180 + int((edgew+180)/5.)*5
     ereglon = -180 + int((edgee+180-0.5*dll)/5.)*5
     IF (ereglon < sreglon) ereglon = sreglon

     ! start loop regions
     print *, ">>> start looping regions..."
     DO reglat = sreglat, ereglat, -5
        DO reglon_ = sreglon, ereglon, 5

           reglon = reglon_
           ! IF lon > 180, revise it to nagative value
           IF (reglon_ >= 180) reglon = reglon_ - 360
           reg(1) = reglat
           reg(2) = reglon
           reg(3) = reglat - 5
           reg(4) = reglon + 5

           ! get region file name and open nc file
           write(reg1, "(i4)") reg(1)
           write(reg2, "(i4)") reg(2)
           write(reg3, "(i4)") reg(3)
           write(reg4, "(i4)") reg(4)

           lndnames = trim(dir_rawdata)//'srf_5x5/RG_' &
              //trim(adjustL(reg1))//'_' &
              //trim(adjustL(reg2))//'_' &
              //trim(adjustL(reg3))//'_' &
              //trim(adjustL(reg4))//'.' &
              //DATASRC//trim(syear)//'.nc'

           lndnamee = trim(dir_rawdata)//'srf_5x5/RG_' &
           //trim(adjustL(reg1))//'_' &
           //trim(adjustL(reg2))//'_' &
           //trim(adjustL(reg3))//'_' &
           //trim(adjustL(reg4))//'.' &
           //DATASRC//trim(eyear)//'.nc'

           print *,">>> Processing file ", trim(lndnames), trim(lndnamee), "..."
           inquire (file=lndnames, exist=fileExists)
           IF (fileExists) THEN
              CALL nccheck( nf90_open(trim(lndnames), nf90_nowrite, ncid) )
           ELSE
              print *, "Warning: file ", trim(lndnames), " does not exist!"
              print *, "All zero value assumed! Please Check!"
              CYCLE
           ENDIF

           ! get the raw data
           CALL nccheck( nf90_inq_varid(ncid, "LC" , lcs_vid) )

           CALL nccheck( nf90_get_var(ncid, lcs_vid, lcdatas) )

           ! close file
           CALL nccheck( nf90_close(ncid) )


           inquire (file=lndnamee, exist=fileExists)
           IF (fileExists) THEN
              CALL nccheck( nf90_open(trim(lndnamee), nf90_nowrite, ncid) )
           ELSE
              print *, "Warning: file ", trim(lndnamee), " does not exist!"
              print *, "All zero value assumed! Please Check!"
              CYCLE
           ENDIF

           ! get the raw data
           CALL nccheck( nf90_inq_varid(ncid, "LC" , lce_vid) )

           CALL nccheck( nf90_get_var(ncid, lce_vid, lcdatae) )

           ! close file
           CALL nccheck( nf90_close(ncid) )


           ! calculate the edge of small grids
           DO i = 1, nxy
              lonw(i) = reg(2) + i*dll - dll
              lone(i) = reg(2) + i*dll
              latn(i) = reg(1) - i*dll + dll
              lats(i) = reg(1) - i*dll
           ENDDO

           ! calculate the area size of small grids
           DO i = 1, nxy
              dx = (lone(1)-lonw(1))*deg2rad
              dy = sin(latn(i)*deg2rad) - sin(lats(i)*deg2rad)
              sarea(:,i) = dx*dy*re*re
           ENDDO

           ! default value
           si = 1; ei = nxy
           sj = 1; ej = nxy

           ! calculate start i/j and END i/j
           IF (reglat == sreglat) si = int((reg(1)-edgen)*nxy/5)+1
           IF (reglon == sreglon) sj = int((edgew-reg(2))*nxy/5)+1
           IF (reglat == ereglat) ei = int((reg(1)-edges)*nxy/5)
           IF (reglon == ereglon) ej = int((edgee-reg(2))*nxy/5)
           IF (ei < si) ei = si
           IF (ej < sj) ej = sj

           ! loop for each small grid for aggregation
           DO i = si, ei
              DO j = sj, ej

                 lat = reg(1) - (i-1)*dll - 0.5*dll
                 lon = reg(2) + (j-1)*dll + 0.5*dll

                 IF (dllo > 0) THEN
                    io = int((edgen-lat)/dllo) + 1
                    IF (lon > edgew) THEN
                       jo = int((lon-edgew)/dllo) + 1
                    ELSE
                       jo = int((lon+360-edgew)/dllo) + 1
                    ENDIF
                 ELSE
                    io = 1
                    jo = 1
                 ENDIF

                 ! total area include ocean
                 area(jo,io) = area(jo,io) + sarea(j,i)

                 lcs = lcdatas(j,i)
                 lce = lcdatae(j,i)

                 ! make ocean the 18th index
                 IF (lcs == 0) lcs = 18
                 IF (lce == 0) lce = 18

                 sumarea(jo,io,lcs,lce) = sumarea(jo,io,lcs,lce) + sarea(j,i)

              ENDDO
           ENDDO
        ENDDO
     ENDDO

     print *, ">>> making surface data..."

    !OMP PARALLEL DO NUM_THREADS(92) &
    !OMP PRIVATE(io,jo,sumpct) &

     DO io = 1, lat_points
        DO jo = 1, lon_points

           ! lcc fraction
           lccpct(jo,io,:,:) = sumarea(jo,io,:,:) / area(jo,io) * 100.

           ! nccheck
           sumpct = sum(lccpct(jo,io,:,:))
           IF (abs(sumpct-100) > 1e-3) THEN
              print *, sumpct
              print *, lccpct(jo,io,:,:)
              print *, "Sum of all lcc pairs not equal 1! STOP!"
           ENDIF

        ENDDO
     ENDDO
    !OMP END PARALLEL DO

     print *, ">>> writing out surface data file ..."
     ! create NC file
     lndname = trim(dir_restart)//trim(syear)//'_'//trim(eyear)//'_'//trim(casename)//'.'//DATASRC//'.nc'
     print*, lndname

     CALL nccheck( nf90_create(lndname, NF90_NETCDF4, ncid) )

     ! Define the dimensions.
     CALL nccheck( nf90_def_dim(ncid, "lat",  lat_points , lat_dimid) )
     CALL nccheck( nf90_def_dim(ncid, "lon",  lon_points , lon_dimid) )
     CALL nccheck( nf90_def_dim(ncid, "lcs",  nlc        , lcs_dimid) )
     CALL nccheck( nf90_def_dim(ncid, "lce",  nlc        , lce_dimid) )

     ! Define the coordinate variables.
     CALL nccheck( nf90_def_var(ncid, "lat" , NF90_FLOAT, lat_dimid , lat_vid) )
     CALL nccheck( nf90_def_var(ncid, "lon" , NF90_FLOAT, lon_dimid , lon_vid) )
     CALL nccheck( nf90_def_var(ncid, "lcs" , NF90_INT  , lcs_dimid , lcs_vid) )
     CALL nccheck( nf90_def_var(ncid, "lce" , NF90_INT  , lce_dimid , lce_vid) )


     ! Assign units attributes to coordinate variables.
     CALL nccheck( nf90_put_att(ncid, lat_vid , "long_name", "Latitude"      ) )
     CALL nccheck( nf90_put_att(ncid, lat_vid , "units"    , "degrees_north" ) )
     CALL nccheck( nf90_put_att(ncid, lon_vid , "long_name", "Longitude"     ) )
     CALL nccheck( nf90_put_att(ncid, lon_vid , "units"    , "degrees_east"  ) )
     CALL nccheck( nf90_put_att(ncid, lcs_vid , "long_name", "LC start index") )
     CALL nccheck( nf90_put_att(ncid, lce_vid , "long_name", "LC end index"  ) )


     ! define output variables
     XY2D = (/ lon_dimid, lat_dimid /)
     CALL nccheck( nf90_def_var(ncid, "AREA"         , NF90_FLOAT, XY2D, varea   , deflate_level=6) )

     LC4D  = (/ lon_dimid, lat_dimid, lcs_dimid , lce_dimid /)
     CALL nccheck( nf90_def_var(ncid, "TRANSFER_PCT" , NF90_FLOAT, LC4D, vlccpct , deflate_level=6) )

     ! Assign units attributes to the netCDF variables.
     CALL nccheck( nf90_put_att(ncid, varea  , "units"    , "km^2"                          ) )
     CALL nccheck( nf90_put_att(ncid, varea  , "long_name", "Area of grid"                  ) )

     CALL nccheck( nf90_put_att(ncid, vlccpct, "units"    , "%"                             ) )
     CALL nccheck( nf90_put_att(ncid, vlccpct, "long_name", "Percent land cover type change") )

     CALL nccheck( nf90_put_att(ncid, NF90_GLOBAL, 'Title'  , Title  ) )
     CALL nccheck( nf90_put_att(ncid, NF90_GLOBAL, 'Authors', Authors) )
     CALL nccheck( nf90_put_att(ncid, NF90_GLOBAL, 'Address', Address) )
     CALL nccheck( nf90_put_att(ncid, NF90_GLOBAL, 'Email'  , Email  ) )

     ! End define mode.
     CALL nccheck( nf90_enddef(ncid) )

     CALL nccheck( nf90_put_var(ncid, lat_vid, latso      ) )
     CALL nccheck( nf90_put_var(ncid, lon_vid, lonso      ) )

     CALL nccheck( nf90_put_var(ncid, lcs_vid, lct(1:nlc) ) )
     CALL nccheck( nf90_put_var(ncid, lce_vid, lct(1:nlc) ) )

     ! put variables
     CALL nccheck( nf90_put_var(ncid, varea  , area       ) )
     CALL nccheck( nf90_put_var(ncid, vlccpct, lccpct     ) )

     ! Close the file. This causes netCDF to flush all buffers and make
     ! sure your data are really written to disk.
     CALL nccheck( nf90_close(ncid) )

     print *,"*** SUCCESS writing file ", trim(lndname), "!"

     ! deallocate memory
     deallocate( lcdatas  )
     deallocate( lcdatae  )
     deallocate( area     )
     deallocate( sumarea  )
     deallocate( lccpct   )

  END SUBROUTINE MAKE_LuLccTransferMatrix


  SUBROUTINE READ_LuLccTransferMatrix(casename,idate,dir_restart)

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
      INTEGER :: lccpct_vid !lcfr_vid, lcto_vid, dchg_vid, selfchg_vid
      INTEGER :: i, j, m, npatch
      INTEGER :: k,numtmp
      REAL(r8):: tmp
      INTEGER,  intent(in) :: idate(3)
      CHARACTER(LEN=256), intent(in) :: casename
      CHARACTER(LEN=256), intent(in) :: dir_restart

      REAL(r8), allocatable :: lccpctin(:,:,:,:)

      write(syear,'(i4.4)') idate(1) - 1
      write(eyear,'(i4.4)') idate(1)

      lndname = trim(dir_restart)//trim(syear)//'_'//trim(eyear)//'_'//trim(casename)//'.'//DATASRC//'.nc'
      print*,'read lulccdata',trim(lndname)
      CALL nccheck( nf90_open(trim(lndname), nf90_nowrite, ncid) )

#ifdef IGBP_CLASSIFICATION
      allocate ( lccpctin (1:lon_points, 1:lat_points, 1:nlc, 1:nlc) )

      CALL nccheck( nf90_inq_varid(ncid, "TRANSFER_PCT", lccpct_vid) )

      CALL nccheck( nf90_get_var(ncid, lccpct_vid, lccpctin, &
                    start=(/1,1,1,1/), &
                    count=(/lon_points,lat_points,nlc,nlc/)) )

      allocate ( lccpct  (numpatch,1:nlc, 1:nlc) )


#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j,m)
#endif
      DO npatch = 1, numpatch

         i = patch2lon(npatch)
         j = patch2lat(npatch)
         m = patchclass(npatch)

        lccpct  (npatch,:,:) = lccpctin (i,j,:,:)

      ENDDO

#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      deallocate ( lccpctin )

#endif

      CALL nccheck( nf90_close(ncid) )

  END SUBROUTINE READ_LuLccTransferMatrix

  SUBROUTINE deallocate_LuLccTransferMatrix
! --------------------------------------------------
! Deallocates memory for LuLcc time invariant variables
! --------------------------------------------------

!TODO: need coding below...


  END SUBROUTINE deallocate_LuLccTransferMatrix

END MODULE MOD_LuLccTransferMatrix
! ---------- EOP ------------
