#include <define.h>

  PROGRAM CLM
! ======================================================================
! Reference:
!     [1] Dai et al., 2003: The Common Land Model (CoLM).
!         Bull. of Amer. Meter. Soc., 84: 1013-1023
!     [2] Dai et al., 2004: A two-big-leaf model for canopy temperature,
!         photosynthesis and stomatal conductance. J. Climate, 17: 2281-2299.
!     [3] Dai et al., 2014: The Terrestrial Modeling System (TMS).
!     [4] Dai Yamazaki, 2014: The global river model CaMa-Flood (version 3.6.2)
!
!     Created by Yongjiu Dai, Februay 2004
!     Revised by Yongjiu Dai and Hua Yuan, April 2014
! ======================================================================

      USE precision
      USE PhysicalConstants
      USE GlobalVars
      USE LC_Const
      USE PFT_Const
      USE MOD_TimeInvariants
      USE MOD_TimeVariables
      USE MOD_1D_Forcing
      USE MOD_2D_Forcing
      USE MOD_1D_Fluxes
      USE MOD_2D_Fluxes
      USE MOD_vec2xy
      USE timemanager
      USE GETMETMOD
      USE omp_lib

#if(defined CaMa_Flood)
      USE parkind1   ,only: jpim, jprm
      USE mod_input  ,only: lognam, nxin, nyin
      USE mod_map    ,only: regionall, regionthis
      USE mod_output ,only: csufbin, csufvec, csufcdf
#endif

      IMPLICIT NONE

! ----------------local variables ---------------------------------

      CHARACTER(LEN=256) :: casename  !casename name
      INTEGER :: lc_year          !which year of land cover data used
      INTEGER :: met_year         !which year of meterology data used for fixed_year case
      INTEGER :: idate(3)         !calendar (year, julian day, seconds)
      INTEGER :: edate(3)         !calendar (year, julian day, seconds)
      INTEGER :: pdate(3)         !calendar (year, julian day, seconds)
      INTEGER :: ldate(3)         !calendar (year, julian day, seconds)
      REAL(r8):: deltim           !time step (senconds)
      LOGICAL :: solarin_all_band !downward solar in broad band
      LOGICAL :: greenwich        !greenwich time

      CHARACTER(len=256) :: dir_srfdata   !surface data directory
      CHARACTER(len=256) :: dir_atmdata   !forcing data directory
      CHARACTER(len=256) :: dir_output    !output data directory
      CHARACTER(len=256) :: dir_restart   !restart data directory
      CHARACTER(len=256) :: nam_atmdata   !forcing data filename
      CHARACTER(len=256) :: nam_srfdata   !surface data filename
      CHARACTER(len=256) :: nam_urbdata   !urban data filename
      CHARACTER(len=256) :: cdate         !string date format

      CHARACTER(len=256) :: dir_rawdata, mksrf_file
      real(r8) :: edgen      ! northern edge of grid (degrees)
      real(r8) :: edgee      ! eastern edge of grid (degrees)
      real(r8) :: edges      ! southern edge of grid (degrees)
      real(r8) :: edgew      ! western edge of grid (degrees)

      LOGICAL :: doalb            !true => start up the surface albedo calculation
      LOGICAL :: dolai            !true => start up the time-varying vegetation paramter
      LOGICAL :: dosst            !true => update sst/ice/snow
      LOGICAL :: lwrite           !true: write output  file frequency
      LOGICAL :: rwrite           !true: write restart file frequency

      INTEGER :: istep            !looping step
      INTEGER :: nac              !number of accumulation
      INTEGER :: nac_24           !number of accumulation of days (24h)
      INTEGER,  allocatable :: nac_ln(:,:)      !number of accumulation for local noon time virable
      INTEGER,  allocatable :: nac_dt(:,:)      !number of accumulation for daytime virable
      INTEGER,  allocatable :: nac_nt(:,:)      !number of accumulation for night-time virable
      REAL(r8), allocatable :: oro(:)           !ocean(0)/seaice(2)/ flag
      REAL(r8), allocatable :: a_rnof(:,:)      !total runoff [mm/s]
      INTEGER :: Julian_1day_p, Julian_1day
      INTEGER :: Julian_8day_p, Julian_8day
      INTEGER :: s_year, s_julian, s_seconds
      INTEGER :: e_year, e_julian, e_seconds
      INTEGER :: p_year, p_julian, p_seconds
      INTEGER :: i, j
      INTEGER :: s_month, e_month, p_month
      INTEGER :: s_day, e_day, p_day
      INTEGER :: year, month, mday, month_p, mday_p

      TYPE(timestamp) :: itstamp, etstamp, ptstamp

#if(defined CaMa_Flood)
      INTEGER(kind=jpim) :: iyyyy, imm, idd          !start date
      INTEGER(kind=jpim) :: eyyyy, emm, edd          !end date
      REAL(kind=jprm), allocatable :: r2roffin(:,:)  !input runoff (mm/day)
#ifdef usempi
      include 'mpif.h'
      INTEGER(kind=jpim) :: ierr, nproc, nid
#endif
#endif

      namelist /clmexp/ casename,               &! 1
                        dir_srfdata,            &! 2.1
                        dir_atmdata,            &! 2.2
                        dir_output,             &! 2.3
                        dir_restart,            &! 2.4
                        nam_atmdata,            &! 3.1
                        nam_srfdata,            &! 3.2
                        nam_urbdata,            &! 3.3
                        deltim,                 &! 5
                        solarin_all_band,       &! 6
                        met_year,               &
                        lc_year,                &! 7
                        e_year,                 &! 8.1
                        e_month,                &! 8.2
                        e_day,                  &! 8.3
                        e_seconds,              &! 8.4
                        p_year,                 &! 9.1
                        p_month,                &! 9.2
                        p_day,                  &! 9.3
                        p_seconds,              &! 9.4
                        numpatch,               &!10.1
                        numpft,                 &!10.2
                        numpc,                  &!10.3
                        numurban,               &!10.4
                        greenwich,              &!11
                        s_year,                 &!12.1
                        s_month,                &!12.2
                        s_day,                  &!12.3
                        s_seconds                !12.4

      namelist /mksrfexp/  casename,dir_rawdata,dir_srfdata,&
                           lc_year,edgen,edgee,edges,edgew

! ======================================================================
!     define the run and open files (for off-line use)

      read(5,clmexp)

      mksrf_file = trim(dir_output)//'../'//'mksrf.stdin'
      open(55, status='OLD', file=mksrf_file, form="FORMATTED")
      read(55, nml=mksrfexp)

      CALL Init_GlovalVars
      CALL Init_LC_Const
      CALL Init_PFT_Const
      CALL initimetype(greenwich)

      idate(1) = s_year; idate(3) = s_seconds
      edate(1) = e_year; edate(3) = e_seconds
      pdate(1) = p_year; pdate(3) = p_seconds

      CALL monthday2julian(s_year,s_month,s_day,idate(2))
      CALL monthday2julian(e_year,e_month,e_day,edate(2))
      CALL monthday2julian(p_year,p_month,p_day,pdate(2))

      s_julian = idate(2); e_julian = edate(2); p_julian = pdate(2)

      CALL adj2end(edate)
      CALL adj2end(pdate)

      itstamp = idate
      etstamp = edate
      ptstamp = pdate

      CALL allocate_TimeInvariants
      CALL allocate_TimeVariables
      CALL allocate_1D_Forcing
      CALL allocate_2D_Forcing
      CALL allocate_1D_Fluxes
      CALL allocate_2D_Fluxes
      CALL allocate_vec2xy

      CALL FLUSH_2D_Fluxes

      allocate (oro(numpatch))
      allocate (nac_ln(lon_points,lat_points))
      allocate (nac_dt(lon_points,lat_points))
      allocate (nac_nt(lon_points,lat_points))
! ----------------------------------------------------------------------
    ! Read in the model time invariant constant data
      CALL READ_TimeInvariants(lc_year,dir_restart,casename)

    ! Read in the model time varying data (model state variables)
      CALL READ_TimeVariables (idate,lc_year,dir_restart,casename)


!-----------------------
#if(defined CaMa_Flood)
#if(defined usempi)
      CALL mpi_init(ierr)
      CALL mpi_comm_size(mpi_comm_world, nproc, ierr)
      CALL mpi_comm_rank(mpi_comm_world, nid, ierr)
      regionall =nproc
      regionthis=nid+1
#endif
      !! regional output for mpi run
      !! change suffix of output file for each calculation node
      IF (regionall>=2 )THEN
        write(csufbin,'(a5,i2.2)') '.bin-', regionthis
        write(csufvec,'(a5,i2.2)') '.vec-', regionthis
        write(csufcdf,'(a4,i2.2)') '.nc-',  regionthis
      ENDIF

      iyyyy = s_year
      eyyyy = e_year
      CALL julian2monthday(s_year,s_julian,imm,idd)
      CALL julian2monthday(e_year,e_julian,emm,edd)
      CALL CaMaINI(iyyyy,imm,idd,eyyyy,emm,edd)

      nxin = lon_points
      nyin = lat_points
      allocate (r2roffin(nxin,nyin))  !!input runoff (mm/day)
#endif
      allocate (a_rnof(lon_points,lat_points))
!-----------------------

      doalb = .true.
      dolai = .true.
      dosst = .false.
      oro(:) = 1.

    ! Initialize meteorological forcing data module
      CALL GETMETINI(dir_atmdata, nam_atmdata, deltim)

! ======================================================================
! begin time stepping loop
! ======================================================================

      nac = 0; nac_24 = 0; nac_ln(:,:) = 0
      nac_dt(:,:) = 0; nac_nt(:,:) = 0
      tmax(:) = 0.; tmin(:) = 330.
      istep = 1

      ! date in beginning style
      ldate = idate
      CALL adj2begin(ldate)

      TIMELOOP : DO while (itstamp < etstamp)

         CALL julian2monthday (ldate(1), ldate(2), month_p, mday_p)
         write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') ldate(1),month_p,mday_p,ldate(3)
         print*, 'TIMELOOP: ', istep, "| DATE: ", trim(cdate)

         Julian_1day_p = int(calendarday(ldate)-1)/1*1 + 1
         Julian_8day_p = int(calendarday(ldate)-1)/8*8 + 1

       ! Read in the meteorological forcing
       ! ----------------------------------------------------------------------
         CALL rd_forcing(idate,solarin_all_band,s_year,s_month,s_day,s_seconds,deltim,s_julian,met_year)

       ! Calendar for NEXT time step
       ! ----------------------------------------------------------------------
         CALL TICKTIME (deltim,idate)
         itstamp = itstamp + int(deltim)
         ldate = idate
         CALL adj2begin(ldate)

       ! Call clm driver
       ! ----------------------------------------------------------------------
         CALL CLMDRIVER (idate,deltim,dolai,doalb,dosst,oro)

       ! Get leaf area index
       ! ----------------------------------------------------------------------
#if(!defined DYN_PHENOLOGY)
       ! READ in Leaf area index and stem area index
       ! Update every 8 days (time interval of the MODIS LAI data)
       ! ----------------------------------------------------------------------
#ifdef USGS_CLASSIFICATION
       ! READ in Leaf area index and stem area index
         Julian_8day = int(calendarday(ldate)-1)/8*8 + 1
         IF (Julian_8day /= Julian_8day_p) THEN
            CALL LAI_readin (Julian_8day,numpatch,dir_srfdata)
         ENDIF

#else
! 08/03/2019, yuan: read global LAI/SAI data
         CALL julian2monthday (ldate(1), ldate(2), month, mday)
         year = ldate(1)
         IF (month /= month_p) THEN
#ifdef LAICHANGE
            CALL LAI_readin_nc      (   year,month,dir_srfdata,nam_srfdata)
#ifdef URBAN_MODEL
            CALL UrbanLAI_readin_nc (   year,month,dir_srfdata,nam_urbdata)
#endif
#else
            CALL LAI_readin_nc      (lc_year,month,dir_srfdata,nam_srfdata)
#ifdef URBAN_MODEL
            CALL UrbanLAI_readin_nc (lc_year,month,dir_srfdata,nam_urbdata)
#endif
#endif
         ENDIF
#endif

#else
       ! Update once a day
         dolai = .false.
         Julian_1day = int(calendarday(ldate)-1)/1*1 + 1
         IF (Julian_1day /= Julian_1day_p) THEN
            dolai = .true.
         ENDIF
#endif

       ! Mapping subgrid patch [numpatch] vector of subgrid points to
       !     -> [lon_points]x[lat_points] grid average
       ! ----------------------------------------------------------------------
         CALL vec2xy (istep,deltim,nac,nac_24,nac_ln,nac_dt,nac_nt,a_rnof)

         WHERE (a_rnof < 1.e-10) a_rnof = 0.

#if(defined CaMa_Flood)
       ! Simulating the hydrodynamics in continental-scale rivers
       ! ----------------------------------------------------------------------

         r2roffin(:,:) = a_rnof(:,:)*86400.  !total runoff [mm/s] -> [mm/day]
         CALL CaMaMAIN(iyyyy,imm,idd,istep,r2roffin)

#endif

       ! Logical idenfication for writing output and restart file
         lwrite = .false.
         rwrite = .false.
         CALL lpwrite(idate,deltim,lwrite,rwrite)

         ! 07/10/2017
         ! for the last time step, write the restart file
         IF (itstamp == etstamp) THEN
            rwrite = .true.
         ENDIF

       ! Write out the model variables for restart run and the histroy file
       ! ----------------------------------------------------------------------
         IF ( lwrite ) THEN

            IF ( .not. (itstamp<=ptstamp) ) THEN
               CALL flxwrite (idate,nac,nac_24,nac_ln,nac_dt,nac_nt,dir_output,casename)
            ENDIF

          ! Setting for next output
          ! ----------------------------------------------------------------------
            CALL FLUSH_2D_Fluxes
            nac = 0; nac_24 = 0; nac_ln(:,:) = 0
            nac_dt(:,:) = 0; nac_nt(:,:) = 0
         ENDIF

#ifdef LULCC
         ! DO land use and land cover change simulation
         IF ( isendofyear(idate, deltim) ) THEN
            CALL deallocate_1D_Forcing
            CALL deallocate_1D_Fluxes

            CALL LuLccDRIVER (casename,dir_srfdata,dir_restart,&
                              nam_srfdata,nam_urbdata,idate,greenwich,&
                              dir_rawdata,edgen,edgee,edges,edgew)

            CALL allocate_1D_Forcing
            CALL allocate_1D_Fluxes
         ENDIF
#endif

         IF ( rwrite ) THEN
            ! output restart file for the last timestep of spin-up
            IF ( .not. (itstamp<ptstamp) ) THEN
#ifdef LULCC
               CALL WRITE_TimeVariables (idate,   year,dir_restart,casename)
#else
               CALL WRITE_TimeVariables (idate,lc_year,dir_restart,casename)
#endif
            ENDIF
         ENDIF

         istep = istep + 1

      ENDDO TIMELOOP

      CALL deallocate_TimeInvariants
      CALL deallocate_TimeVariables
      CALL deallocate_1D_Forcing
      CALL deallocate_2D_Forcing
      CALL deallocate_1D_Fluxes
      CALL deallocate_2D_Fluxes
      CALL deallocate_vec2xy

      CALL GETMETFINAL

      deallocate (a_rnof)
      deallocate (oro)
      deallocate (nac_ln)
      deallocate (nac_dt)
      deallocate (nac_nt)
#if(defined CaMa_Flood)
      deallocate (r2roffin)
      IF (regionall>=2 )THEN
        close(lognam)
      ENDIF
#ifdef usempi
      CALL mpi_finalize(ierr)
#endif
#endif

      write(6,*) 'CLM Execution Completed'

  END PROGRAM CLM
! ----------------------------------------------------------------------
! EOP
