#include <define.h>

MODULE MOD_vec2xy
! ----------------------------------------------------------------------
! perfrom the grid average mapping: average a subgrid input 1d vector
! of length numpatch to a output 2d array of length [lon_points,lat_points]
!
! Created by Yongjiu Dai, 03/2014
! Modified by yuan, 06/2022: change to MODULE
! ----------------------------------------------------------------------

   USE precision
   USE GlobalVars

   IMPLICIT NONE
   SAVE

   !---------------------------------------------------------------------
   REAL(r8), allocatable :: a_xy_us  (:,:)    !wind in eastward direction [m/s]
   REAL(r8), allocatable :: a_xy_vs  (:,:)    !wind in northward direction [m/s]
   REAL(r8), allocatable :: a_xy_t   (:,:)    !temperature at reference height [kelvin]
   REAL(r8), allocatable :: a_xy_q   (:,:)    !specific humidity at reference height [kg/kg]
   REAL(r8), allocatable :: a_xy_prc (:,:)    !convective precipitation [mm/s]
   REAL(r8), allocatable :: a_xy_prl (:,:)    !large scale precipitation [mm/s]
   REAL(r8), allocatable :: a_xy_pbot(:,:)    !atmospheric pressure at the surface [pa]
   REAL(r8), allocatable :: a_xy_frl (:,:)    !atmospheric infrared (longwave) radiation [W/m2]
   REAL(r8), allocatable :: a_xy_solarin(:,:) !downward solar radiation at surface [W/m2]
   REAL(r8), allocatable :: a_xy_rain(:,:)    !rain [mm/s]
   REAL(r8), allocatable :: a_xy_snow(:,:)    !snow [mm/s]

   !---------------------------------------------------------------------
   REAL(r8), allocatable :: a_taux   (:,:)    !wind stress: E-W [kg/m/s2]
   REAL(r8), allocatable :: a_tauy   (:,:)    !wind stress: N-S [kg/m/s2]
   REAL(r8), allocatable :: a_fsena  (:,:)    !sensible heat from canopy height to atmosphere [W/m2]
   REAL(r8), allocatable :: a_lfevpa (:,:)    !latent heat flux from canopy height to atmosphere [W/m2]
   REAL(r8), allocatable :: a_fevpa  (:,:)    !evapotranspiration from canopy to atmosphere [mm/s]
   REAL(r8), allocatable :: a_fsenl  (:,:)    !sensible heat from leaves [W/m2]
   REAL(r8), allocatable :: a_fevpl  (:,:)    !evaporation+transpiration from leaves [mm/s]
   REAL(r8), allocatable :: a_etr    (:,:)    !transpiration rate [mm/s]
   REAL(r8), allocatable :: a_fseng  (:,:)    !sensible heat flux from ground [W/m2]
   REAL(r8), allocatable :: a_fevpg  (:,:)    !evaporation heat flux from ground [mm/s]
   REAL(r8), allocatable :: a_fgrnd  (:,:)    !ground heat flux [W/m2]
   REAL(r8), allocatable :: a_sabvsun(:,:)    !solar absorbed by sunlit canopy [W/m2]
   REAL(r8), allocatable :: a_sabvsha(:,:)    !solar absorbed by shaded [W/m2]
   REAL(r8), allocatable :: a_sabg   (:,:)    !solar absorbed by ground  [W/m2]
   REAL(r8), allocatable :: a_olrg   (:,:)    !outgoing long-wave radiation from ground+canopy [W/m2]
   REAL(r8), allocatable :: a_rnet   (:,:)    !net radiation [W/m2]
   REAL(r8), allocatable :: a_xerr   (:,:)    !the error of water banace [mm/s]
   REAL(r8), allocatable :: a_zerr   (:,:)    !the error of energy balance [W/m2]
   REAL(r8), allocatable :: a_rsur   (:,:)    !surface runoff [mm/s]
   REAL(r8), allocatable :: a_qintr  (:,:)    !interception [mm/s]
   REAL(r8), allocatable :: a_qinfl  (:,:)    !inflitration [mm/s]
   REAL(r8), allocatable :: a_qdrip  (:,:)    !throughfall [mm/s]

   REAL(r8), allocatable :: a_assim  (:,:)    !canopy assimilation rate [mol m-2 s-1]
   REAL(r8), allocatable :: a_respc  (:,:)    !respiration (plant+soil) [mol m-2 s-1]
   REAL(r8), allocatable :: a_qcharge(:,:)    !groundwater recharge rate [mm/s]

   !---------------------------------------------------------------------
   REAL(r8), allocatable :: a_t_grnd (:,:)    !ground surface temperature [K]
   REAL(r8), allocatable :: a_tleaf  (:,:)    !sunlit leaf temperature [K]
   REAL(r8), allocatable :: a_ldew   (:,:)    !depth of water on foliage [mm]
   REAL(r8), allocatable :: a_scv    (:,:)    !snow cover, water equivalent [mm]
   REAL(r8), allocatable :: a_snowdp (:,:)    !snow depth [meter]
   REAL(r8), allocatable :: a_fsno   (:,:)    !fraction of snow cover on ground
   REAL(r8), allocatable :: a_sigf   (:,:)    !fraction of veg cover, excluding snow-covered veg [-]
   REAL(r8), allocatable :: a_green  (:,:)    !leaf greenness
   REAL(r8), allocatable :: a_lai    (:,:)    !leaf area index
   REAL(r8), allocatable :: a_laisun (:,:)    !sunlit leaf area index
   REAL(r8), allocatable :: a_laisha (:,:)    !shaded leaf area index
   REAL(r8), allocatable :: a_sai    (:,:)    !stem area index
   REAL(r8), allocatable :: a_alb(:,:,:,:)    !averaged albedo [visible, direct; direct, diffuse]
   REAL(r8), allocatable :: a_emis   (:,:)    !averaged bulk surface emissivity
   REAL(r8), allocatable :: a_z0m    (:,:)    !effective roughness [m]
   REAL(r8), allocatable :: a_trad   (:,:)    !radiative temperature of surface [K]
   REAL(r8), allocatable :: a_tref   (:,:)    !2 m height air temperature [kelvin]
   REAL(r8), allocatable :: a_tmax   (:,:)    !Diurnal Max 2 m height air temperature [kelvin]
   REAL(r8), allocatable :: a_tmin   (:,:)    !Diurnal Min 2 m height air temperature [kelvin]
   REAL(r8), allocatable :: a_tdtr   (:,:)    !DTR of 2 m height air temperature [kelvin]
   REAL(r8), allocatable :: a_qref   (:,:)    !2 m height air specific humidity [kg/kg]

   !---------------------------------------------------------------------
   REAL(r8), allocatable :: a_t_room (:,:)    !temperature of inner building [K]
   REAL(r8), allocatable :: a_tafu   (:,:)    !temperature of outer building [K]
   REAL(r8), allocatable :: a_fhac   (:,:)    !sensible flux from heat or cool AC [W/m2]
   REAL(r8), allocatable :: a_fwst   (:,:)    !waste heat flux from heat or cool AC [W/m2]
   REAL(r8), allocatable :: a_fach   (:,:)    !flux from inner and outter air exchange [W/m2]
   REAL(r8), allocatable :: a_fahe   (:,:)
   REAL(r8), allocatable :: a_fhah   (:,:)
   REAL(r8), allocatable :: a_vehc   (:,:)
   REAL(r8), allocatable :: a_meta   (:,:)

   REAL(r8), allocatable :: a_senroof(:,:)
   REAL(r8), allocatable :: a_senwsun(:,:)
   REAL(r8), allocatable :: a_senwsha(:,:)
   REAL(r8), allocatable :: a_sengimp(:,:)
   REAL(r8), allocatable :: a_sengper(:,:)
   REAL(r8), allocatable :: a_senurl (:,:)

   REAL(r8), allocatable :: a_lfevproof(:,:)
   REAL(r8), allocatable :: a_lfevpgimp(:,:)
   REAL(r8), allocatable :: a_lfevpgper(:,:)
   REAL(r8), allocatable :: a_lfevpurl (:,:)

   REAL(r8), allocatable :: a_troof    (:,:)
   REAL(r8), allocatable :: a_twall    (:,:)

   REAL(r8), allocatable :: a_sabvdt  (:,:)   !solar absorbed by sunlit canopy [W/m2]
   REAL(r8), allocatable :: a_sabgdt  (:,:)   !solar absorbed by ground [W/m2]
   REAL(r8), allocatable :: a_srdt    (:,:)   !total reflected solar radiation (W/m2)
   REAL(r8), allocatable :: a_fsenadt (:,:)   !sensible heat from canopy height to atmosphere [W/m2]
   REAL(r8), allocatable :: a_lfevpadt(:,:)   !latent heat flux from canopy height to atmosphere [W/m2]
   REAL(r8), allocatable :: a_fgrnddt (:,:)   !ground heat flux [W/m2]
   REAL(r8), allocatable :: a_olrgdt  (:,:)   !outgoing long-wave radiation from ground+canopy [W/m2]
   REAL(r8), allocatable :: a_rnetdt  (:,:)   !net radiation [W/m2]
   REAL(r8), allocatable :: a_t_grnddt(:,:)   !ground surface temperature [K]
   REAL(r8), allocatable :: a_traddt  (:,:)   !radiative temperature of surface [K]
   REAL(r8), allocatable :: a_trefdt  (:,:)   !2 m height air temperature [kelvin]
   REAL(r8), allocatable :: a_tafudt  (:,:)   !temperature of outer building [K]

   REAL(r8), allocatable :: a_fsenant (:,:)   !sensible heat from canopy height to atmosphere [W/m2]
   REAL(r8), allocatable :: a_lfevpant(:,:)   !latent heat flux from canopy height to atmosphere [W/m2]
   REAL(r8), allocatable :: a_fgrndnt (:,:)   !ground heat flux [W/m2]
   REAL(r8), allocatable :: a_olrgnt  (:,:)   !outgoing long-wave radiation from ground+canopy [W/m2]
   REAL(r8), allocatable :: a_rnetnt  (:,:)   !net radiation [W/m2]
   REAL(r8), allocatable :: a_t_grndnt(:,:)   !ground surface temperature [K]
   REAL(r8), allocatable :: a_tradnt  (:,:)   !radiative temperature of surface [K]
   REAL(r8), allocatable :: a_trefnt  (:,:)   !2 m height air temperature [kelvin]
   REAL(r8), allocatable :: a_tafunt  (:,:)   !temperature of outer building [K]

   !---------------------------------------------------------------------
   REAL(r8), allocatable :: a_t_soisno   (:,:,:)  !soil temperature [K]
   REAL(r8), allocatable :: a_wliq_soisno(:,:,:)  !liquid water in soil layers [kg/m2]
   REAL(r8), allocatable :: a_wice_soisno(:,:,:)  !ice lens in soil layers [kg/m2]
   REAL(r8), allocatable :: a_h2osoi     (:,:,:)  !volumetric soil water in layers [m3/m3]
   REAL(r8), allocatable :: a_rstfac     (:,:)    !factor of soil water stress
   REAL(r8), allocatable :: a_zwt        (:,:)    !the depth to water table [m]
   REAL(r8), allocatable :: a_wa         (:,:)    !water storage in aquifer [mm]
   REAL(r8), allocatable :: a_wat        (:,:)    !total water storage [mm]

   REAL(r8), allocatable :: a_t_lake      (:,:,:) !lake temperature [K]
   REAL(r8), allocatable :: a_lake_icefrac(:,:,:) !lake ice fraction cover [0-1]

   !---------------------------------------------------------------------
   REAL(r8), allocatable :: a_ustar  (:,:)    !u* in similarity theory [m/s]
   REAL(r8), allocatable :: a_tstar  (:,:)    !t* in similarity theory [kg/kg]
   REAL(r8), allocatable :: a_qstar  (:,:)    !q* in similarity theory [kg/kg]
   REAL(r8), allocatable :: a_zol    (:,:)    !dimensionless height (z/L) used in Monin-Obukhov theory
   REAL(r8), allocatable :: a_rib    (:,:)    !bulk Richardson number in surface layer
   REAL(r8), allocatable :: a_fm     (:,:)    !integral of profile function for momentum
   REAL(r8), allocatable :: a_fh     (:,:)    !integral of profile function for heat
   REAL(r8), allocatable :: a_fq     (:,:)    !integral of profile function for moisture

   REAL(r8), allocatable :: a_us10m  (:,:)    !10m u-velocity [m/s]
   REAL(r8), allocatable :: a_vs10m  (:,:)    !10m v-velocity [m/s]
   REAL(r8), allocatable :: a_fm10m  (:,:)    !integral of profile function for momentum at 10m [-]

   !---------------------------------------------------------------------
   REAL(r8), allocatable :: a_sr     (:,:)    !total reflected solar radiation (W/m2)
   REAL(r8), allocatable :: a_solvd  (:,:)    !incident direct beam vis solar radiation (W/m2)
   REAL(r8), allocatable :: a_solvi  (:,:)    !incident diffuse beam vis solar radiation (W/m2)
   REAL(r8), allocatable :: a_solnd  (:,:)    !incident direct beam nir solar radiation (W/m2)
   REAL(r8), allocatable :: a_solni  (:,:)    !incident diffuse beam nir solar radiation (W/m2)
   REAL(r8), allocatable :: a_srvd   (:,:)    !reflected direct beam vis solar radiation (W/m2)
   REAL(r8), allocatable :: a_srvi   (:,:)    !reflected diffuse beam vis solar radiation (W/m2)
   REAL(r8), allocatable :: a_srnd   (:,:)    !reflected direct beam nir solar radiation (W/m2)
   REAL(r8), allocatable :: a_srni   (:,:)    !reflected diffuse beam nir solar radiation (W/m2)
   REAL(r8), allocatable :: a_solvdln(:,:)    !incident direct beam vis solar radiation at local noon (W/m2)
   REAL(r8), allocatable :: a_solviln(:,:)    !incident diffuse beam vis solar radiation at local noon (W/m2)
   REAL(r8), allocatable :: a_solndln(:,:)    !incident direct beam nir solar radiation at local noon (W/m2)
   REAL(r8), allocatable :: a_solniln(:,:)    !incident diffuse beam nir solar radiation at local noon (W/m2)
   REAL(r8), allocatable :: a_srvdln (:,:)    !reflected direct beam vis solar radiation at local noon (W/m2)
   REAL(r8), allocatable :: a_srviln (:,:)    !reflected diffuse beam vis solar radiation at local noon (W/m2)
   REAL(r8), allocatable :: a_srndln (:,:)    !reflected direct beam nir solar radiation at local noon (W/m2)
   REAL(r8), allocatable :: a_srniln (:,:)    !reflected diffuse beam nir solar radiation at local noon (W/m2)

   PUBLIC  :: vec2xy
   PUBLIC  :: allocate_vec2xy
   PUBLIC  :: deallocate_vec2xy
   PRIVATE :: acc

CONTAINS

   SUBROUTINE vec2xy (istep,deltim,nac,nac_24,nac_ln,nac_dt,nac_nt,a_rnof)

      USE PhysicalConstants, only: vonkar, stefnc, cpair, rgas, grav
      USE MOD_TimeInvariants
      USE MOD_TimeVariables
      USE MOD_1D_Forcing
      USE MOD_2D_Forcing
      USE MOD_1D_Fluxes
      USE MOD_2D_Fluxes
      USE MOD_UrbanTimeVars
      USE MOD_UrbanTimeInvars
      USE FRICTION_VELOCITY
      USE omp_lib

      IMPLICIT NONE

      INTEGER, intent(in)    :: istep
      INTEGER, intent(inout) :: nac
      INTEGER, intent(inout) :: nac_24
      INTEGER, intent(inout) :: nac_ln(lon_points,lat_points)
      INTEGER, intent(inout) :: nac_dt(lon_points,lat_points)
      INTEGER, intent(inout) :: nac_nt(lon_points,lat_points)

      REAL(r8),intent(in ) :: deltim   !seconds in a time-step
      REAL(r8),intent(out) :: a_rnof(lon_points,lat_points)  !total runoff [mm/s]

    !---------------------------------------------------------------------
    ! local variables
      INTEGER  i,j,np,u,l,daystep
      REAL(r8) sumwt(lon_points,lat_points)
      REAL(r8) urbwt(lon_points,lat_points)
      REAL(r8) rhoair,thm,th,thv,ur,displa_av,zldis,hgt_u,hgt_t,hgt_q
      REAL(r8) z0m_av,z0h_av,z0q_av,us,vs,tm,qm,psrf
      REAL(r8) obu,fh2m,fq2m
      REAL(r8) um,thvstar,beta,zii,wc,wc2

      daystep = int(86400/deltim)

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j)
#endif
      DO j = 1, lat_points
         DO i = 1, lon_points

          ! ---------------------------------------------------
          ! Meteorological forcing
          ! ---------------------------------------------------
            a_xy_us     (i,j) = forc_xy_us     (i,j)
            a_xy_vs     (i,j) = forc_xy_vs     (i,j)
            a_xy_t      (i,j) = forc_xy_t      (i,j)
            a_xy_q      (i,j) = forc_xy_q      (i,j)
            a_xy_prc    (i,j) = forc_xy_prc    (i,j)
            a_xy_prl    (i,j) = forc_xy_prl    (i,j)
            a_xy_pbot   (i,j) = forc_xy_pbot   (i,j)
            a_xy_frl    (i,j) = forc_xy_frl    (i,j)

            a_xy_solarin(i,j) = forc_xy_sols (i,j) + forc_xy_soll (i,j) &
                              + forc_xy_solsd(i,j) + forc_xy_solld(i,j)

          ! ------------------------------------------------------------------------------------------
          ! Mapping the fluxes and state variables at patch [numpatch] to grid [lon_points,lat_points]
          ! ------------------------------------------------------------------------------------------
            sumwt     (i,j) = 0.
            urbwt     (i,j) = 0.
            a_taux    (i,j) = 0.
            a_tauy    (i,j) = 0.
            a_fsena   (i,j) = 0.
            a_lfevpa  (i,j) = 0.
            a_fevpa   (i,j) = 0.
            a_fsenl   (i,j) = 0.
            a_fevpl   (i,j) = 0.
            a_etr     (i,j) = 0.
            a_fseng   (i,j) = 0.
            a_fevpg   (i,j) = 0.
            a_fgrnd   (i,j) = 0.
            a_sabvsun (i,j) = 0.
            a_sabvsha (i,j) = 0.
            a_sabg    (i,j) = 0.
            a_olrg    (i,j) = 0.
            a_rnet    (i,j) = 0.
            a_xerr    (i,j) = 0.
            a_zerr    (i,j) = 0.
            a_rsur    (i,j) = 0.
            a_rnof    (i,j) = 0.
            a_qintr   (i,j) = 0.
            a_qinfl   (i,j) = 0.
            a_qdrip   (i,j) = 0.
            a_wat     (i,j) = 0.
            a_assim   (i,j) = 0.
            a_respc   (i,j) = 0.

            a_qcharge (i,j) = 0.
            a_t_grnd  (i,j) = 0.
            a_tleaf   (i,j) = 0.
            a_ldew    (i,j) = 0.
            a_scv     (i,j) = 0.
            a_snowdp  (i,j) = 0.
            a_fsno    (i,j) = 0.
            a_sigf    (i,j) = 0.
            a_green   (i,j) = 0.
            a_lai     (i,j) = 0.
            a_laisun  (i,j) = 0.
            a_laisha  (i,j) = 0.
            a_sai     (i,j) = 0.
            a_alb(:,: ,i,j) = 0.
            a_emis    (i,j) = 0.
            a_z0m     (i,j) = 0.
            a_trad    (i,j) = 0.
            a_tref    (i,j) = 0.
            a_tmax    (i,j) = 0.
            a_tmin    (i,j) = 0.
            a_tdtr    (i,j) = 0.
            a_qref    (i,j) = 0.
            a_xy_rain (i,j) = 0.
            a_xy_snow (i,j) = 0.

            a_t_room  (i,j) = 0.
            a_tafu    (i,j) = 0.
            a_fhac    (i,j) = 0.
            a_fwst    (i,j) = 0.
            a_fach    (i,j) = 0.
            a_fahe    (i,j) = 0.
            a_fhah    (i,j) = 0.
            a_vehc    (i,j) = 0.
            a_meta    (i,j) = 0.

            a_senroof (i,j) = 0.!spval
            a_senwsun (i,j) = 0.!spval
            a_senwsha (i,j) = 0.!spval
            a_sengimp (i,j) = 0.!spval
            a_sengper (i,j) = 0.!spval
            a_senurl  (i,j) = 0.!spval

            a_lfevproof(i,j)= 0.!spval
            a_lfevpgimp(i,j)= 0.!spval
            a_lfevpgper(i,j)= 0.!spval
            a_lfevpurl (i,j)= 0.!spval

            a_troof   (i,j) = 0.!spval
            a_twall   (i,j) = 0.!spval

            a_sabvdt  (i,j) = spval
            a_sabgdt  (i,j) = spval
            a_srdt    (i,j) = spval
            a_fsenadt (i,j) = spval
            a_lfevpadt(i,j) = spval
            a_fgrnddt (i,j) = spval
            a_olrgdt  (i,j) = spval
            a_rnetdt  (i,j) = spval
            a_t_grnddt(i,j) = spval
            a_traddt  (i,j) = spval
            a_trefdt  (i,j) = spval
            a_tafudt  (i,j) = spval

            a_fsenant (i,j) = spval
            a_lfevpant(i,j) = spval
            a_fgrndnt (i,j) = spval
            a_olrgnt  (i,j) = spval
            a_rnetnt  (i,j) = spval
            a_t_grndnt(i,j) = spval
            a_tradnt  (i,j) = spval
            a_trefnt  (i,j) = spval
            a_tafunt  (i,j) = spval

            a_sr      (i,j) = spval
            a_solvd   (i,j) = spval
            a_solvi   (i,j) = spval
            a_solnd   (i,j) = spval
            a_solni   (i,j) = spval
            a_srvd    (i,j) = spval
            a_srvi    (i,j) = spval
            a_srnd    (i,j) = spval
            a_srni    (i,j) = spval
            a_solvdln (i,j) = spval
            a_solviln (i,j) = spval
            a_solndln (i,j) = spval
            a_solniln (i,j) = spval
            a_srvdln  (i,j) = spval
            a_srviln  (i,j) = spval
            a_srndln  (i,j) = spval
            a_srniln  (i,j) = spval
         ENDDO
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,np,u)
#endif
      DO j = 1, lat_points
         DO i = 1, lon_points

            IF (grid_patch_s(i,j) .le. 0) cycle

            DO np = grid_patch_s(i,j), grid_patch_e(i,j)
! 10/05/2021, yuan: only for urban output
#ifdef URBAN_ONLY
               IF (patchclass(np) .ne. URBAN) cycle
#endif
               sumwt(i,j) = sumwt(i,j) + patchfrac(np)

             ! Fluxes
               a_taux   (i,j) = a_taux   (i,j) + patchfrac(np)*taux   (np)
               a_tauy   (i,j) = a_tauy   (i,j) + patchfrac(np)*tauy   (np)
               a_fsena  (i,j) = a_fsena  (i,j) + patchfrac(np)*fsena  (np)
               a_lfevpa (i,j) = a_lfevpa (i,j) + patchfrac(np)*lfevpa (np)
               a_fevpa  (i,j) = a_fevpa  (i,j) + patchfrac(np)*fevpa  (np)
               a_fsenl  (i,j) = a_fsenl  (i,j) + patchfrac(np)*fsenl  (np)
               a_fevpl  (i,j) = a_fevpl  (i,j) + patchfrac(np)*fevpl  (np)
               a_etr    (i,j) = a_etr    (i,j) + patchfrac(np)*etr    (np)
               a_fseng  (i,j) = a_fseng  (i,j) + patchfrac(np)*fseng  (np)
               a_fevpg  (i,j) = a_fevpg  (i,j) + patchfrac(np)*fevpg  (np)
               a_fgrnd  (i,j) = a_fgrnd  (i,j) + patchfrac(np)*fgrnd  (np)
               a_sabvsun(i,j) = a_sabvsun(i,j) + patchfrac(np)*sabvsun(np)
               a_sabvsha(i,j) = a_sabvsha(i,j) + patchfrac(np)*sabvsha(np)
               a_sabg   (i,j) = a_sabg   (i,j) + patchfrac(np)*sabg   (np)
               a_olrg   (i,j) = a_olrg   (i,j) + patchfrac(np)*olrg   (np)
               a_rnet   (i,j) = a_rnet   (i,j) + patchfrac(np)*(sabg(np)+sabvsun(np)+sabvsha(np)-olrg(np))
               a_xerr   (i,j) = a_xerr   (i,j) + patchfrac(np)*xerr   (np)
               a_zerr   (i,j) = a_zerr   (i,j) + patchfrac(np)*zerr   (np)
               a_rsur   (i,j) = a_rsur   (i,j) + patchfrac(np)*rsur   (np)
               a_rnof   (i,j) = a_rnof   (i,j) + patchfrac(np)*rnof   (np)
               a_qintr  (i,j) = a_qintr  (i,j) + patchfrac(np)*qintr  (np)
               a_qinfl  (i,j) = a_qinfl  (i,j) + patchfrac(np)*qinfl  (np)
               a_qdrip  (i,j) = a_qdrip  (i,j) + patchfrac(np)*qdrip  (np)
               a_wat    (i,j) = a_wat    (i,j) + patchfrac(np)*wat    (np)
               a_assim  (i,j) = a_assim  (i,j) + patchfrac(np)*assim  (np)
               a_respc  (i,j) = a_respc  (i,j) + patchfrac(np)*respc  (np)

               a_qcharge(i,j) = a_qcharge(i,j) + patchfrac(np)*qcharge(np)

             ! State and other variables
               a_t_grnd (i,j) = a_t_grnd (i,j) + patchfrac(np)*t_grnd (np)
               a_tleaf  (i,j) = a_tleaf  (i,j) + patchfrac(np)*tleaf  (np)
               a_ldew   (i,j) = a_ldew   (i,j) + patchfrac(np)*ldew   (np)
               a_scv    (i,j) = a_scv    (i,j) + patchfrac(np)*scv    (np)
               a_snowdp (i,j) = a_snowdp (i,j) + patchfrac(np)*snowdp (np)
               a_fsno   (i,j) = a_fsno   (i,j) + patchfrac(np)*fsno   (np)
               a_sigf   (i,j) = a_sigf   (i,j) + patchfrac(np)*sigf   (np)
               a_green  (i,j) = a_green  (i,j) + patchfrac(np)*green  (np)
               a_lai    (i,j) = a_lai    (i,j) + patchfrac(np)*lai    (np)
               a_laisun (i,j) = a_laisun (i,j) + patchfrac(np)*laisun (np)
               a_laisha (i,j) = a_laisha (i,j) + patchfrac(np)*laisha (np)
               a_sai    (i,j) = a_sai    (i,j) + patchfrac(np)*sai    (np)
               a_alb(1,1,i,j) = a_alb(1,1,i,j) + patchfrac(np)*alb(1,1,np)
               a_alb(1,2,i,j) = a_alb(1,2,i,j) + patchfrac(np)*alb(1,2,np)
               a_alb(2,1,i,j) = a_alb(2,1,i,j) + patchfrac(np)*alb(2,1,np)
               a_alb(2,2,i,j) = a_alb(2,2,i,j) + patchfrac(np)*alb(2,2,np)
               a_emis   (i,j) = a_emis   (i,j) + patchfrac(np)*emis   (np)
               a_z0m    (i,j) = a_z0m    (i,j) + patchfrac(np)*z0m    (np)
               a_tref   (i,j) = a_tref   (i,j) + patchfrac(np)*tref   (np)
               a_qref   (i,j) = a_qref   (i,j) + patchfrac(np)*qref   (np)
               a_xy_rain(i,j) = a_xy_rain(i,j) + patchfrac(np)*forc_rain(np)
               a_xy_snow(i,j) = a_xy_snow(i,j) + patchfrac(np)*forc_snow(np)

               IF (mod(istep,daystep) == 0) THEN
                  a_tmax(i,j) = a_tmax   (i,j) + patchfrac(np)*tmax   (np)
                  a_tmin(i,j) = a_tmin   (i,j) + patchfrac(np)*tmin   (np)
                  a_tdtr(i,j) = a_tdtr   (i,j) + patchfrac(np)*(tmax(np)-tmin(np))
                  tmax(np) = 0.; tmin(np) = 330.
               ENDIF

#ifdef URBAN_MODEL
               u = patch2urb(np)
               IF (u > 0) THEN
                  urbwt(i,j) = urbwt(i,j) + patchfrac(np)
                  a_t_room (i,j) = a_t_room (i,j) + patchfrac(np)*t_room (u)
                  a_tafu   (i,j) = a_tafu   (i,j) + patchfrac(np)*tafu   (u)
                  a_fhac   (i,j) = a_fhac   (i,j) + patchfrac(np)*fhac   (u)
                  a_fwst   (i,j) = a_fwst   (i,j) + patchfrac(np)*fwst   (u)
                  a_fach   (i,j) = a_fach   (i,j) + patchfrac(np)*fach   (u)
                  a_fahe   (i,j) = a_fahe   (i,j) + patchfrac(np)*fahe   (u)
                  a_fhah   (i,j) = a_fhah   (i,j) + patchfrac(np)*fhah   (u)
                  a_vehc   (i,j) = a_vehc   (i,j) + patchfrac(np)*vehc   (u)
                  a_meta   (i,j) = a_meta   (i,j) + patchfrac(np)*meta   (u)

                  ! print*, fsen_roof, patchfrac(np)
                  a_senroof(i,j) = a_senroof(i,j) + patchfrac(np)*fsen_roof(u)
                  a_senwsun(i,j) = a_senwsun(i,j) + patchfrac(np)*fsen_wsun(u)
                  a_senwsha(i,j) = a_senwsha(i,j) + patchfrac(np)*fsen_wsha(u)
                  a_sengimp(i,j) = a_sengimp(i,j) + patchfrac(np)*fsen_gimp(u)
                  a_sengper(i,j) = a_sengper(i,j) + patchfrac(np)*fsen_gper(u)
                  a_senurl (i,j) = a_senurl (i,j) + patchfrac(np)*fsen_url (u)

                  a_lfevproof(i,j) = a_lfevproof(i,j) + patchfrac(np)*lfevp_roof(u)
                  a_lfevpgimp(i,j) = a_lfevpgimp(i,j) + patchfrac(np)*lfevp_gimp(u)
                  a_lfevpgper(i,j) = a_lfevpgper(i,j) + patchfrac(np)*lfevp_gper(u)
                  a_lfevpurl (i,j) = a_lfevpurl (i,j) + patchfrac(np)*lfevp_url (u)

                  ! print*, troof
                  a_troof  (i,j) = a_troof  (i,j) + patchfrac(np)*troof  (u)
                  a_twall  (i,j) = a_twall  (i,j) + patchfrac(np)*twall  (u)
               ENDIF
#endif
               ! 根据coszen(np)的正负->daytime or nighttime
               IF (coszen(np) > 0.) THEN
                  CALL acc(sabvsun(np), patchfrac(np), a_sabvdt  (i,j))
                  CALL acc(sabg   (np), patchfrac(np), a_sabgdt  (i,j))
                  CALL acc(sr     (np), patchfrac(np), a_srdt    (i,j))
                  CALL acc(fsena  (np), patchfrac(np), a_fsenadt (i,j))
                  CALL acc(lfevpa (np), patchfrac(np), a_lfevpadt(i,j))
                  CALL acc(fgrnd  (np), patchfrac(np), a_fgrnddt (i,j))
                  CALL acc(olrg   (np), patchfrac(np), a_olrgdt  (i,j))
                  CALL acc(t_grnd (np), patchfrac(np), a_t_grnddt(i,j))
                  CALL acc(tref   (np), patchfrac(np), a_trefdt  (i,j))
#ifdef URBAN_MODEL
                  u = patch2urb(np)
                  IF (u > 0) THEN
                     CALL acc(tafu(u ), patchfrac(np), a_tafudt  (i,j))
                  ENDIF
#endif
                  IF (a_rnetdt(i,j) /= spval) THEN
                     a_rnetdt(i,j) = a_rnetdt(i,j) + patchfrac(np)*(sabg(np)+sabvsun(np)+sabvsha(np)-olrg(np))
                  ELSE
                     a_rnetdt(i,j) = patchfrac(np)*(sabg(np)+sabvsun(np)+sabvsha(np)-olrg(np))
                  ENDIF

               ELSE
                  CALL acc(fsena  (np), patchfrac(np), a_fsenant (i,j))
                  CALL acc(lfevpa (np), patchfrac(np), a_lfevpant(i,j))
                  CALL acc(fgrnd  (np), patchfrac(np), a_fgrndnt (i,j))
                  CALL acc(olrg   (np), patchfrac(np), a_olrgnt  (i,j))
                  CALL acc(t_grnd (np), patchfrac(np), a_t_grndnt(i,j))
                  CALL acc(tref   (np), patchfrac(np), a_trefnt  (i,j))
#ifdef URBAN_MODEL
                  u = patch2urb(np)
                  IF (u > 0) THEN
                     CALL acc(tafu(u ), patchfrac(np), a_tafunt  (i,j))
                  ENDIF
#endif
                  IF (a_rnetnt(i,j) /= spval) THEN
                     a_rnetnt(i,j) = a_rnetnt(i,j) + patchfrac(np)*(sabg(np)+sabvsun(np)+sabvsha(np)-olrg(np))
                  ELSE
                     a_rnetnt(i,j) = patchfrac(np)*(sabg(np)+sabvsun(np)+sabvsha(np)-olrg(np))
                  ENDIF

               ENDIF

               ! radiation fluxes
               CALL acc(sr     (np), patchfrac(np), a_sr     (i,j))
               CALL acc(solvd  (np), patchfrac(np), a_solvd  (i,j))
               CALL acc(solvi  (np), patchfrac(np), a_solvi  (i,j))
               CALL acc(solnd  (np), patchfrac(np), a_solnd  (i,j))
               CALL acc(solni  (np), patchfrac(np), a_solni  (i,j))
               CALL acc(srvd   (np), patchfrac(np), a_srvd   (i,j))
               CALL acc(srvi   (np), patchfrac(np), a_srvi   (i,j))
               CALL acc(srnd   (np), patchfrac(np), a_srnd   (i,j))
               CALL acc(srni   (np), patchfrac(np), a_srni   (i,j))
               ! local noon fluxes
               CALL acc(solvdln(np), patchfrac(np), a_solvdln(i,j))
               CALL acc(solviln(np), patchfrac(np), a_solviln(i,j))
               CALL acc(solndln(np), patchfrac(np), a_solndln(i,j))
               CALL acc(solniln(np), patchfrac(np), a_solniln(i,j))
               CALL acc(srvdln (np), patchfrac(np), a_srvdln (i,j))
               CALL acc(srviln (np), patchfrac(np), a_srviln (i,j))
               CALL acc(srndln (np), patchfrac(np), a_srndln (i,j))
               CALL acc(srniln (np), patchfrac(np), a_srniln (i,j))
            ENDDO
            area(i,j) = gridarea(i,j)
         ENDDO
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j)
#endif
      DO j = 1, lat_points
         DO i = 1, lon_points
            IF(sumwt(i,j).gt.1.00001)THEN
               print*, 'summation of fraction patches (1) = ', sumwt(i,j), i,j
            !  CALL abort
            ENDIF
            IF(sumwt(i,j).gt.0.00001)THEN
               a_taux   (i,j) = a_taux   (i,j) / sumwt(i,j)
               a_tauy   (i,j) = a_tauy   (i,j) / sumwt(i,j)
               a_fsena  (i,j) = a_fsena  (i,j) / sumwt(i,j)
               a_lfevpa (i,j) = a_lfevpa (i,j) / sumwt(i,j)
               a_fevpa  (i,j) = a_fevpa  (i,j) / sumwt(i,j)
               a_fsenl  (i,j) = a_fsenl  (i,j) / sumwt(i,j)
               a_fevpl  (i,j) = a_fevpl  (i,j) / sumwt(i,j)
               a_etr    (i,j) = a_etr    (i,j) / sumwt(i,j)
               a_fseng  (i,j) = a_fseng  (i,j) / sumwt(i,j)
               a_fevpg  (i,j) = a_fevpg  (i,j) / sumwt(i,j)
               a_fgrnd  (i,j) = a_fgrnd  (i,j) / sumwt(i,j)
               a_sabvsun(i,j) = a_sabvsun(i,j) / sumwt(i,j)
               a_sabvsha(i,j) = a_sabvsha(i,j) / sumwt(i,j)
               a_sabg   (i,j) = a_sabg   (i,j) / sumwt(i,j)
               a_olrg   (i,j) = a_olrg   (i,j) / sumwt(i,j)
               a_rnet   (i,j) = a_rnet   (i,j) / sumwt(i,j) + a_xy_frl(i,j)
               a_xerr   (i,j) = a_xerr   (i,j) / sumwt(i,j)
               a_zerr   (i,j) = a_zerr   (i,j) / sumwt(i,j)

               a_rsur   (i,j) = a_rsur   (i,j) / sumwt(i,j)
               a_rnof   (i,j) = a_rnof   (i,j) / sumwt(i,j)
               a_qintr  (i,j) = a_qintr  (i,j) / sumwt(i,j)
               a_qinfl  (i,j) = a_qinfl  (i,j) / sumwt(i,j)
               a_qdrip  (i,j) = a_qdrip  (i,j) / sumwt(i,j)
               a_wat    (i,j) = a_wat    (i,j) / sumwt(i,j)
               a_assim  (i,j) = a_assim  (i,j) / sumwt(i,j)
               a_respc  (i,j) = a_respc  (i,j) / sumwt(i,j)

               a_qcharge(i,j) = a_qcharge(i,j) / sumwt(i,j)

               a_t_grnd (i,j) = a_t_grnd (i,j) / sumwt(i,j)
               a_tleaf  (i,j) = a_tleaf  (i,j) / sumwt(i,j)
               a_ldew   (i,j) = a_ldew   (i,j) / sumwt(i,j)
               a_scv    (i,j) = a_scv    (i,j) / sumwt(i,j)
               a_snowdp (i,j) = a_snowdp (i,j) / sumwt(i,j)
               a_fsno   (i,j) = a_fsno   (i,j) / sumwt(i,j)
               a_sigf   (i,j) = a_sigf   (i,j) / sumwt(i,j)
               a_green  (i,j) = a_green  (i,j) / sumwt(i,j)
               a_lai    (i,j) = a_lai    (i,j) / sumwt(i,j)
               a_laisun (i,j) = a_laisun (i,j) / sumwt(i,j)
               a_laisha (i,j) = a_laisha (i,j) / sumwt(i,j)
               a_sai    (i,j) = a_sai    (i,j) / sumwt(i,j)
               a_alb(1,1,i,j) = a_alb(1,1,i,j) / sumwt(i,j)
               a_alb(1,2,i,j) = a_alb(1,2,i,j) / sumwt(i,j)
               a_alb(2,1,i,j) = a_alb(2,1,i,j) / sumwt(i,j)
               a_alb(2,2,i,j) = a_alb(2,2,i,j) / sumwt(i,j)
               a_emis   (i,j) = a_emis   (i,j) / sumwt(i,j)
               a_z0m    (i,j) = a_z0m    (i,j) / sumwt(i,j)
               a_trad   (i,j) = (a_olrg(i,j)/stefnc)**0.25 ! fordebug
               a_tref   (i,j) = a_tref   (i,j) / sumwt(i,j)
               a_qref   (i,j) = a_qref   (i,j) / sumwt(i,j)
               a_xy_rain(i,j) = a_xy_rain(i,j) / sumwt(i,j)
               a_xy_snow(i,j) = a_xy_snow(i,j) / sumwt(i,j)

               IF (mod(istep,daystep) == 0) THEN
                  a_tmax(i,j) = a_tmax   (i,j) / sumwt(i,j)
                  a_tmin(i,j) = a_tmin   (i,j) / sumwt(i,j)
                  a_tdtr(i,j) = a_tdtr   (i,j) / sumwt(i,j)
               ENDIF

#ifdef URBAN_MODEL
               IF(urbwt(i,j).gt.0.00001)THEN
                  a_t_room (i,j) = a_t_room (i,j) / urbwt(i,j)
                  a_tafu   (i,j) = a_tafu   (i,j) / urbwt(i,j)
                  a_fhac   (i,j) = a_fhac   (i,j) / urbwt(i,j)
                  a_fwst   (i,j) = a_fwst   (i,j) / urbwt(i,j)
                  a_fach   (i,j) = a_fach   (i,j) / urbwt(i,j)
                  a_fahe   (i,j) = a_fahe   (i,j) / urbwt(i,j)
                  a_fhah   (i,j) = a_fhah   (i,j) / urbwt(i,j)
                  a_vehc   (i,j) = a_vehc   (i,j) / urbwt(i,j)
                  a_meta   (i,j) = a_meta   (i,j) / urbwt(i,j)

                  a_senroof(i,j) = a_senroof(i,j) / urbwt(i,j)
                  a_senwsun(i,j) = a_senwsun(i,j) / urbwt(i,j)
                  a_senwsha(i,j) = a_senwsha(i,j) / urbwt(i,j)
                  a_sengimp(i,j) = a_sengimp(i,j) / urbwt(i,j)
                  a_sengper(i,j) = a_sengper(i,j) / urbwt(i,j)
                  a_senurl (i,j) = a_senurl (i,j) / urbwt(i,j)

                  a_lfevproof(i,j) = a_lfevproof(i,j) / urbwt(i,j)
                  a_lfevpgimp(i,j) = a_lfevpgimp(i,j) / urbwt(i,j)
                  a_lfevpgper(i,j) = a_lfevpgper(i,j) / urbwt(i,j)
                  a_lfevpurl (i,j) = a_lfevpurl (i,j) / urbwt(i,j)

                  a_troof  (i,j) = a_troof  (i,j) / urbwt(i,j)
                  a_twall  (i,j) = a_twall  (i,j) / urbwt(i,j)
               ELSE
                  a_t_room (i,j) = spval
                  a_tafu   (i,j) = spval
                  a_fhac   (i,j) = spval
                  a_fwst   (i,j) = spval
                  a_fach   (i,j) = spval
                  a_fahe   (i,j) = spval
                  a_senroof(i,j) = spval
                  a_senwsun(i,j) = spval
                  a_senwsha(i,j) = spval
                  a_sengimp(i,j) = spval
                  a_sengper(i,j) = spval
                  a_senurl (i,j) = spval

                  a_lfevproof(i,j) = spval
                  a_lfevpgimp(i,j) = spval
                  a_lfevpgper(i,j) = spval
                  a_lfevpurl (i,j) = spval

                  a_troof  (i,j) = spval
                  a_twall  (i,j) = spval
               ENDIF
#endif
               IF (a_sabvdt  (i,j) /= spval) a_sabvdt  (i,j) = a_sabvdt  (i,j) / sumwt(i,j)
               IF (a_sabgdt  (i,j) /= spval) a_sabgdt  (i,j) = a_sabgdt  (i,j) / sumwt(i,j)
               IF (a_srdt    (i,j) /= spval) a_srdt    (i,j) = a_srdt    (i,j) / sumwt(i,j)
               IF (a_fsenadt (i,j) /= spval) a_fsenadt (i,j) = a_fsenadt (i,j) / sumwt(i,j)
               IF (a_lfevpadt(i,j) /= spval) a_lfevpadt(i,j) = a_lfevpadt(i,j) / sumwt(i,j)
               IF (a_fgrnddt (i,j) /= spval) a_fgrnddt (i,j) = a_fgrnddt (i,j) / sumwt(i,j)
               IF (a_olrgdt  (i,j) /= spval) a_olrgdt  (i,j) = a_olrgdt  (i,j) / sumwt(i,j)
               IF (a_olrgdt  (i,j) /= spval) a_traddt  (i,j) =(a_olrgdt  (i,j) / stefnc)**0.25
               IF (a_rnetdt  (i,j) /= spval) a_rnetdt  (i,j) = a_rnetdt  (i,j) / sumwt(i,j) + a_xy_frl(i,j)
               IF (a_t_grnddt(i,j) /= spval) a_t_grnddt(i,j) = a_t_grnddt(i,j) / sumwt(i,j)
               IF (a_trefdt  (i,j) /= spval) a_trefdt  (i,j) = a_trefdt  (i,j) / sumwt(i,j)
               IF (a_tafudt  (i,j) /= spval) a_tafudt  (i,j) = a_tafudt  (i,j) / sumwt(i,j)

               IF (a_fsenant (i,j) /= spval) a_fsenant (i,j) = a_fsenant (i,j) / sumwt(i,j)
               IF (a_lfevpant(i,j) /= spval) a_lfevpant(i,j) = a_lfevpant(i,j) / sumwt(i,j)
               IF (a_fgrndnt (i,j) /= spval) a_fgrndnt (i,j) = a_fgrndnt (i,j) / sumwt(i,j)
               IF (a_olrgnt  (i,j) /= spval) a_olrgnt  (i,j) = a_olrgnt  (i,j) / sumwt(i,j)
               IF (a_olrgnt  (i,j) /= spval) a_tradnt  (i,j) =(a_olrgnt  (i,j) / stefnc)**0.25
               IF (a_rnetnt  (i,j) /= spval) a_rnetnt  (i,j) = a_rnetnt  (i,j) / sumwt(i,j) + a_xy_frl(i,j)
               IF (a_t_grndnt(i,j) /= spval) a_t_grndnt(i,j) = a_t_grndnt(i,j) / sumwt(i,j)
               IF (a_trefnt  (i,j) /= spval) a_trefnt  (i,j) = a_trefnt  (i,j) / sumwt(i,j)
               IF (a_tafunt  (i,j) /= spval) a_tafunt  (i,j) = a_tafunt  (i,j) / sumwt(i,j)

               IF (a_sr     (i,j) /= spval) a_sr     (i,j) = a_sr     (i,j) / sumwt(i,j)
               IF (a_solvd  (i,j) /= spval) a_solvd  (i,j) = a_solvd  (i,j) / sumwt(i,j)
               IF (a_solvi  (i,j) /= spval) a_solvi  (i,j) = a_solvi  (i,j) / sumwt(i,j)
               IF (a_solnd  (i,j) /= spval) a_solnd  (i,j) = a_solnd  (i,j) / sumwt(i,j)
               IF (a_solni  (i,j) /= spval) a_solni  (i,j) = a_solni  (i,j) / sumwt(i,j)
               IF (a_srvd   (i,j) /= spval) a_srvd   (i,j) = a_srvd   (i,j) / sumwt(i,j)
               IF (a_srvi   (i,j) /= spval) a_srvi   (i,j) = a_srvi   (i,j) / sumwt(i,j)
               IF (a_srnd   (i,j) /= spval) a_srnd   (i,j) = a_srnd   (i,j) / sumwt(i,j)
               IF (a_srni   (i,j) /= spval) a_srni   (i,j) = a_srni   (i,j) / sumwt(i,j)
               IF (a_solvdln(i,j) /= spval) a_solvdln(i,j) = a_solvdln(i,j) / sumwt(i,j)
               IF (a_solviln(i,j) /= spval) a_solviln(i,j) = a_solviln(i,j) / sumwt(i,j)
               IF (a_solndln(i,j) /= spval) a_solndln(i,j) = a_solndln(i,j) / sumwt(i,j)
               IF (a_solniln(i,j) /= spval) a_solniln(i,j) = a_solniln(i,j) / sumwt(i,j)
               IF (a_srvdln (i,j) /= spval) a_srvdln (i,j) = a_srvdln (i,j) / sumwt(i,j)
               IF (a_srviln (i,j) /= spval) a_srviln (i,j) = a_srviln (i,j) / sumwt(i,j)
               IF (a_srndln (i,j) /= spval) a_srndln (i,j) = a_srndln (i,j) / sumwt(i,j)
               IF (a_srniln (i,j) /= spval) a_srniln (i,j) = a_srniln (i,j) / sumwt(i,j)

               mask(i,j) = 1
               frac(i,j) = sumwt(i,j)

            ELSE
               a_taux   (i,j) = spval
               a_tauy   (i,j) = spval
               a_fsena  (i,j) = spval
               a_lfevpa (i,j) = spval
               a_fevpa  (i,j) = spval
               a_fsenl  (i,j) = spval
               a_fevpl  (i,j) = spval
               a_etr    (i,j) = spval
               a_fseng  (i,j) = spval
               a_fevpg  (i,j) = spval
               a_fgrnd  (i,j) = spval
               a_sabvsun(i,j) = spval
               a_sabvsha(i,j) = spval
               a_sabg   (i,j) = spval
               a_olrg   (i,j) = spval
               a_rnet   (i,j) = spval
               a_xerr   (i,j) = spval
               a_zerr   (i,j) = spval

               a_rsur   (i,j) = spval
               a_rnof   (i,j) = spval
               a_qintr  (i,j) = spval
               a_qinfl  (i,j) = spval
               a_qdrip  (i,j) = spval
               a_wat    (i,j) = spval
               a_assim  (i,j) = spval
               a_respc  (i,j) = spval

               a_qcharge(i,j) = spval

               a_t_grnd (i,j) = spval
               a_tleaf  (i,j) = spval
               a_ldew   (i,j) = spval
               a_scv    (i,j) = spval
               a_snowdp (i,j) = spval
               a_fsno   (i,j) = spval
               a_sigf   (i,j) = spval
               a_green  (i,j) = spval
               a_lai    (i,j) = spval
               a_laisun (i,j) = spval
               a_laisha (i,j) = spval
               a_sai    (i,j) = spval
               a_alb(1,1,i,j) = spval
               a_alb(1,2,i,j) = spval
               a_alb(2,1,i,j) = spval
               a_alb(2,2,i,j) = spval
               a_emis   (i,j) = spval
               a_z0m    (i,j) = spval
               a_trad   (i,j) = spval
               a_tref   (i,j) = spval
               a_tmax   (i,j) = spval
               a_tmin   (i,j) = spval
               a_tdtr   (i,j) = spval
               a_qref   (i,j) = spval
               a_xy_rain(i,j) = spval
               a_xy_snow(i,j) = spval

               a_sabvdt  (i,j) = spval
               a_sabgdt  (i,j) = spval
               a_srdt    (i,j) = spval
               a_fsenadt (i,j) = spval
               a_lfevpadt(i,j) = spval
               a_fgrnddt (i,j) = spval
               a_olrgdt  (i,j) = spval
               a_rnetdt  (i,j) = spval
               a_t_grnddt(i,j) = spval
               a_traddt  (i,j) = spval
               a_trefdt  (i,j) = spval
               a_tafudt  (i,j) = spval

               a_fsenant (i,j) = spval
               a_lfevpant(i,j) = spval
               a_fgrndnt (i,j) = spval
               a_olrgnt  (i,j) = spval
               a_rnetnt  (i,j) = spval
               a_t_grndnt(i,j) = spval
               a_tradnt  (i,j) = spval
               a_trefnt  (i,j) = spval
               a_tafunt  (i,j) = spval

               a_sr     (i,j) = spval
               a_solvd  (i,j) = spval
               a_solvi  (i,j) = spval
               a_solnd  (i,j) = spval
               a_solni  (i,j) = spval
               a_srvd   (i,j) = spval
               a_srvi   (i,j) = spval
               a_srnd   (i,j) = spval
               a_srni   (i,j) = spval
               a_solvdln(i,j) = spval
               a_solviln(i,j) = spval
               a_solndln(i,j) = spval
               a_solniln(i,j) = spval
               a_srvdln (i,j) = spval
               a_srviln (i,j) = spval
               a_srndln (i,j) = spval
               a_srniln (i,j) = spval

               mask(i,j) = 0
               frac(i,j) = 0.
            ENDIF

         ENDDO
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

   ! --------------------------------------------------------------------
   ! Temperature and water (excluding land water bodies and ocean patches)
   ! [soil => 0; urban and built-up => 1; wetland => 2; land ice => 3; land water bodies => 4; ocean => 99]
   ! --------------------------------------------------------------------
      sumwt(:,:) = 0.
      a_t_soisno   (:,:,:) = 0.
      a_wliq_soisno(:,:,:) = 0.
      a_wice_soisno(:,:,:) = 0.

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,np)
#endif
      DO j = 1, lat_points
         DO i = 1, lon_points

            IF (grid_patch_s(i,j) .le. 0) cycle
            DO np = grid_patch_s(i,j), grid_patch_e(i,j)

   ! 10/05/2021, yuan: only for urban output
#ifdef URBAN_ONLY
               IF (patchclass(np) .ne. URBAN) cycle
#endif
               IF(patchtype(np) <= 3)THEN  ! excluded the land water bodies and ocean patches
                  sumwt(i,j) = sumwt(i,j) + patchfrac(np)
                  a_t_soisno   (maxsnl+1:nl_soil,i,j) = a_t_soisno   (maxsnl+1:nl_soil,i,j) &
                     + patchfrac(np)*t_soisno   (maxsnl+1:nl_soil,np)
                  a_wliq_soisno(maxsnl+1:nl_soil,i,j) = a_wliq_soisno(maxsnl+1:nl_soil,i,j) &
                     + patchfrac(np)*wliq_soisno(maxsnl+1:nl_soil,np)
                  a_wice_soisno(maxsnl+1:nl_soil,i,j) = a_wice_soisno(maxsnl+1:nl_soil,i,j) &
                     + patchfrac(np)*wice_soisno(maxsnl+1:nl_soil,np)
               ENDIF
            ENDDO

         ENDDO
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,l)
#endif
      DO j = 1, lat_points
         DO i = 1, lon_points
            IF(sumwt(i,j).gt.1.00001)THEN
               print*, 'summation of fraction patches (2) = ', sumwt(i,j), i,j
            !  CALL abort
            ENDIF
            IF(sumwt(i,j).gt.0.00001)THEN
               DO l = maxsnl+1, nl_soil
                  a_t_soisno   (l,i,j) = a_t_soisno   (l,i,j) / sumwt(i,j)
                  a_wliq_soisno(l,i,j) = a_wliq_soisno(l,i,j) / sumwt(i,j)
                  a_wice_soisno(l,i,j) = a_wice_soisno(l,i,j) / sumwt(i,j)
               ENDDO
            ELSE
               a_t_soisno   (maxsnl+1:nl_soil,i,j) = spval
               a_wliq_soisno(maxsnl+1:nl_soil,i,j) = spval
               a_wice_soisno(maxsnl+1:nl_soil,i,j) = spval
            ENDIF

         ENDDO
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

   ! --------------------------------------------------------------------
   ! additial diagnostic variables for output (vegetated land only <=2)
   ! [soil => 0; urban and built-up => 1; wetland => 2; land ice => 3; land water bodies => 4; ocean => 99]
   ! --------------------------------------------------------------------
      sumwt(:,:) = 0.
      a_h2osoi (:,:,:) = 0.
      a_rstfac (:,:)   = 0.
      a_zwt    (:,:)   = 0.
      a_wa     (:,:)   = 0.

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,np)
#endif
      DO j = 1, lat_points
         DO i = 1, lon_points

            IF (grid_patch_s(i,j) .le. 0) cycle
            DO np = grid_patch_s(i,j), grid_patch_e(i,j)

   ! 10/05/2021, yuan: only for urban output
#ifdef URBAN_ONLY
               IF (patchclass(np) .ne. URBAN) cycle
#endif
               IF(patchtype(np) <= 2)THEN  ! excluded the land water bodies and ocean patches
                  sumwt(i,j) = sumwt(i,j) + patchfrac(np)
                  a_h2osoi(1:nl_soil,i,j) = a_h2osoi(1:nl_soil,i,j) &
                     + patchfrac(np)*h2osoi(1:nl_soil,np)
                  a_rstfac(i,j) = a_rstfac(i,j) + patchfrac(np)*rstfac(np)
                  a_zwt   (i,j) = a_zwt   (i,j) + patchfrac(np)*zwt   (np)
                  a_wa    (i,j) = a_wa    (i,j) + patchfrac(np)*wa    (np)
               ENDIF
            ENDDO

         ENDDO
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,l)
#endif
      DO j = 1, lat_points
         DO i = 1, lon_points
            IF(sumwt(i,j).gt.1.00001)THEN
               print*, 'summation of fraction patches (2) = ', sumwt(i,j), i,j
            !  CALL abort
            ENDIF
            IF(sumwt(i,j).gt.0.00001)THEN
               DO l = 1, nl_soil
                  a_h2osoi (l,i,j) = a_h2osoi (l,i,j) / sumwt(i,j)
               ENDDO
               a_rstfac (i,j) = a_rstfac (i,j) / sumwt(i,j)
               a_zwt    (i,j) = a_zwt    (i,j) / sumwt(i,j)
               a_wa     (i,j) = a_wa     (i,j) / sumwt(i,j)
            ELSE
               a_h2osoi (:,i,j) = spval
               a_rstfac (i,j)   = spval
               a_zwt    (i,j)   = spval
               a_wa     (i,j)   = spval
            ENDIF

         ENDDO
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

   ! -----------------------------------------------
   ! Land water bodies' ice fraction and temperature
   ! -----------------------------------------------
      sumwt(:,:) = 0.
      a_t_lake(:,:,:) = 0.
      a_lake_icefrac(:,:,:) = 0.

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,np,u)
#endif
      DO j = 1, lat_points
         DO i = 1, lon_points

            IF (grid_patch_s(i,j) .le. 0) cycle
            DO np = grid_patch_s(i,j), grid_patch_e(i,j)

   ! 10/05/2021, yuan: only for urban output
#ifdef URBAN_MODEL
               IF (patchclass(np) .ne. URBAN) THEN
                  cycle
               ELSE
                  u = patch2urb(np)
                  IF (u > 0) THEN
                  sumwt(i,j) = sumwt(i,j) + patchfrac(np)*flake(u)
                  a_t_lake(1:nl_lake,i,j) = a_t_lake(1:nl_lake,i,j) + patchfrac(np)*flake(u)*t_lake(1:nl_lake,np)
                  a_lake_icefrac(1:nl_lake,i,j) = a_lake_icefrac(1:nl_lake,i,j) + patchfrac(np)*flake(u)*lake_icefrac(1:nl_lake,np)
                  ENDIF
               ENDIF
#endif
               IF(patchtype(np) == 4)THEN  ! land water bodies only
                  sumwt(i,j) = sumwt(i,j) + patchfrac(np)
                  a_t_lake(1:nl_lake,i,j) = a_t_lake(1:nl_lake,i,j) + patchfrac(np)*t_lake(1:nl_lake,np)
                  a_lake_icefrac(1:nl_lake,i,j) = a_lake_icefrac(1:nl_lake,i,j) + patchfrac(np)*lake_icefrac(1:nl_lake,np)
               ENDIF
            ENDDO

         ENDDO
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,l)
#endif
      DO j = 1, lat_points
         DO i = 1, lon_points
            IF(sumwt(i,j).gt.1.00001)THEN
               print*, 'summation of fraction patches (3) = ', sumwt(i,j), i,j
            !  CALL abort
            ENDIF
            IF(sumwt(i,j).gt.0.00001)THEN
               DO l = 1, nl_lake
                  a_t_lake(l,i,j) = a_t_lake(l,i,j) / sumwt(i,j)
                  a_lake_icefrac(l,i,j) = a_lake_icefrac(l,i,j) / sumwt(i,j)
               ENDDO
            ELSE
               a_t_lake(1:nl_lake,i,j) = spval
               a_lake_icefrac(1:nl_lake,i,j) = spval
            ENDIF
         ENDDO
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif


   ! --------------------------------
   ! Retrieve through averaged fluxes
   ! --------------------------------
      sumwt(:,:) = 0.

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,np)
#endif
      DO j = 1, lat_points
         DO i = 1, lon_points

            IF (grid_patch_s(i,j) .le. 0) cycle
            DO np = grid_patch_s(i,j), grid_patch_e(i,j)

   ! 10/05/2021, yuan: only for urban output
#ifdef URBAN_ONLY
               IF (patchclass(np) .ne. URBAN) cycle
#endif
               sumwt(i,j) = sumwt(i,j) + patchfrac(np)
            ENDDO

         ENDDO
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j) &
!$OMP PRIVATE(rhoair,thm,th,thv,ur,displa_av,zldis,hgt_u,hgt_t,hgt_q) &
!$OMP PRIVATE(z0m_av,z0h_av,z0q_av,us,vs,tm,qm,psrf) &
!$OMP PRIVATE(obu,fh2m,fq2m) &
!$OMP PRIVATE(um,thvstar,beta,zii,wc,wc2)
#endif
      DO j = 1, lat_points
         DO i = 1, lon_points
            IF(sumwt(i,j) > 0.00001)THEN
               z0m_av = a_z0m (i,j)
               z0h_av = a_z0m (i,j)
               z0q_av = a_z0m (i,j)
               displa_av = 2./3.*z0m_av/0.07

               hgt_u = max(forc_xy_hgt_u(i,j),5.+displa_av)
               hgt_t = max(forc_xy_hgt_t(i,j),5.+displa_av)
               hgt_q = max(forc_xy_hgt_q(i,j),5.+displa_av)
               zldis = hgt_u-displa_av

               us = forc_xy_us(i,j)
               vs = forc_xy_vs(i,j)
               tm = forc_xy_t(i,j)
               qm = forc_xy_q(i,j)
               psrf = forc_xy_psrf(i,j)
               rhoair = (psrf - 0.378*qm*psrf/(0.622+0.378*qm)) / (rgas*tm)

               a_ustar(i,j) = sqrt(max(1.e-6,sqrt(a_taux(i,j)**2+a_tauy(i,j)**2))/rhoair)
               a_tstar(i,j) = -a_fsena(i,j)/(rhoair*a_ustar(i,j))/cpair
               a_qstar(i,j) = -a_fevpa(i,j)/(rhoair*a_ustar(i,j))

               thm = tm + 0.0098*hgt_t
               th  = tm*(100000./psrf)**(rgas/cpair)
               thv = th*(1.+0.61*qm)

               a_zol(i,j) = zldis*vonkar*grav &
                          * (a_tstar(i,j)+0.61*th*a_qstar(i,j))/(a_ustar(i,j)**2*thv)

               IF(a_zol(i,j) >= 0.)THEN   !stable
                  a_zol(i,j) = min(2.,max(a_zol(i,j),1.e-6))
               ELSE                       !unstable
                  a_zol(i,j) = max(-100.,min(a_zol(i,j),-1.e-6))
               ENDIF

               beta = 1.
               zii = 1000.
               thvstar=a_tstar(i,j)+0.61*th*a_qstar(i,j)
               ur = sqrt(us*us+vs*vs)
               IF(a_zol(i,j) >= 0.)THEN
                  um = max(ur,0.1)
               ELSE
                  wc = (-grav*a_ustar(i,j)*thvstar*zii/thv)**(1./3.)
                 wc2 = beta*beta*(wc*wc)
                  um = max(0.1,sqrt(ur*ur+wc2))
               ENDIF

               obu = zldis/a_zol(i,j)

   !NOTE: for single point debug [注释下面]
               CALL moninobuk(hgt_u,hgt_t,hgt_q,displa_av,z0m_av,z0h_av,z0q_av,& ! fordebug
                    obu,um,a_ustar(i,j),fh2m,fq2m,&
                    a_fm10m(i,j),a_fm(i,j),a_fh(i,j),a_fq(i,j))

   ! bug found by chen qiying 2013/07/01
               a_rib(i,j) = a_zol(i,j)/vonkar*a_ustar(i,j)**2/(vonkar/a_fh(i,j)*um**2)
               a_rib(i,j) = min(5.,a_rib(i,j))

               a_us10m(i,j) = us/um * a_ustar(i,j)/vonkar * a_fm10m(i,j)
               a_vs10m(i,j) = vs/um * a_ustar(i,j)/vonkar * a_fm10m(i,j)

            ELSE

               a_ustar(i,j) = spval
               a_tstar(i,j) = spval
               a_qstar(i,j) = spval
               a_zol(i,j)   = spval
               a_rib(i,j)   = spval
               a_fm(i,j)    = spval
               a_fh(i,j)    = spval
               a_fq(i,j)    = spval
               a_fm10m(i,j) = spval
               a_us10m(i,j) = spval
               a_vs10m(i,j) = spval

            ENDIF

         ENDDO
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif


   ! ---------------------------------------------------
   ! ACCUMULATION in each time step
   ! ---------------------------------------------------

      nac = nac + 1
      IF (mod(istep,daystep) == 0) nac_24 = nac_24 + 1

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,l)
#endif
      DO j = 1, lat_points
         DO i = 1, lon_points

            CALL acc(a_xy_us     (i,j), 1., f_xy_us     (i,j))
            CALL acc(a_xy_vs     (i,j), 1., f_xy_vs     (i,j))
            CALL acc(a_xy_t      (i,j), 1., f_xy_t      (i,j))
            CALL acc(a_xy_q      (i,j), 1., f_xy_q      (i,j))
            CALL acc(a_xy_prc    (i,j), 1., f_xy_prc    (i,j))
            CALL acc(a_xy_prl    (i,j), 1., f_xy_prl    (i,j))
            CALL acc(a_xy_pbot   (i,j), 1., f_xy_pbot   (i,j))
            CALL acc(a_xy_frl    (i,j), 1., f_xy_frl    (i,j))
            CALL acc(a_xy_solarin(i,j), 1., f_xy_solarin(i,j))

            CALL acc(a_taux      (i,j), 1., f_taux      (i,j))
            CALL acc(a_tauy      (i,j), 1., f_tauy      (i,j))
            CALL acc(a_fsena     (i,j), 1., f_fsena     (i,j))
            CALL acc(a_lfevpa    (i,j), 1., f_lfevpa    (i,j))
            CALL acc(a_fevpa     (i,j), 1., f_fevpa     (i,j))
            CALL acc(a_fsenl     (i,j), 1., f_fsenl     (i,j))
            CALL acc(a_fevpl     (i,j), 1., f_fevpl     (i,j))
            CALL acc(a_etr       (i,j), 1., f_etr       (i,j))
            CALL acc(a_fseng     (i,j), 1., f_fseng     (i,j))
            CALL acc(a_fevpg     (i,j), 1., f_fevpg     (i,j))
            CALL acc(a_fgrnd     (i,j), 1., f_fgrnd     (i,j))
            CALL acc(a_sabvsun   (i,j), 1., f_sabvsun   (i,j))
            CALL acc(a_sabvsha   (i,j), 1., f_sabvsha   (i,j))
            CALL acc(a_sabg      (i,j), 1., f_sabg      (i,j))
            CALL acc(a_olrg      (i,j), 1., f_olrg      (i,j))
            CALL acc(a_rnet      (i,j), 1., f_rnet      (i,j))
            CALL acc(a_xerr      (i,j), 1., f_xerr      (i,j))
            CALL acc(a_zerr      (i,j), 1., f_zerr      (i,j))
            CALL acc(a_rsur      (i,j), 1., f_rsur      (i,j))
            CALL acc(a_rnof      (i,j), 1., f_rnof      (i,j))
            CALL acc(a_qintr     (i,j), 1., f_qintr     (i,j))
            CALL acc(a_qinfl     (i,j), 1., f_qinfl     (i,j))
            CALL acc(a_qdrip     (i,j), 1., f_qdrip     (i,j))
            CALL acc(a_rstfac    (i,j), 1., f_rstfac    (i,j))
            CALL acc(a_zwt       (i,j), 1., f_zwt       (i,j))
            CALL acc(a_wa        (i,j), 1., f_wa        (i,j))
            CALL acc(a_wat       (i,j), 1., f_wat       (i,j))
            CALL acc(a_assim     (i,j), 1., f_assim     (i,j))
            CALL acc(a_respc     (i,j), 1., f_respc     (i,j))

            CALL acc(a_qcharge   (i,j), 1., f_qcharge   (i,j))

            CALL acc(a_t_grnd    (i,j), 1., f_t_grnd    (i,j))
            CALL acc(a_tleaf     (i,j), 1., f_tleaf     (i,j))
            CALL acc(a_ldew      (i,j), 1., f_ldew      (i,j))
            CALL acc(a_scv       (i,j), 1., f_scv       (i,j))
            CALL acc(a_snowdp    (i,j), 1., f_snowdp    (i,j))
            CALL acc(a_fsno      (i,j), 1., f_fsno      (i,j))
            CALL acc(a_sigf      (i,j), 1., f_sigf      (i,j))
            CALL acc(a_green     (i,j), 1., f_green     (i,j))
            CALL acc(a_lai       (i,j), 1., f_lai       (i,j))
            CALL acc(a_laisun    (i,j), 1., f_laisun    (i,j))
            CALL acc(a_laisha    (i,j), 1., f_laisha    (i,j))
            CALL acc(a_sai       (i,j), 1., f_sai       (i,j))
            CALL acc(a_alb   (1,1,i,j), 1., f_alb   (1,1,i,j))
            CALL acc(a_alb   (2,1,i,j), 1., f_alb   (2,1,i,j))
            CALL acc(a_alb   (1,2,i,j), 1., f_alb   (1,2,i,j))
            CALL acc(a_alb   (2,2,i,j), 1., f_alb   (2,2,i,j))
            CALL acc(a_emis      (i,j), 1., f_emis      (i,j))
            CALL acc(a_z0m       (i,j), 1., f_z0m       (i,j))
            CALL acc(a_trad      (i,j), 1., f_trad      (i,j))
            CALL acc(a_tref      (i,j), 1., f_tref      (i,j))
            CALL acc(a_tmax      (i,j), 1., f_tmax      (i,j))
            CALL acc(a_tmin      (i,j), 1., f_tmin      (i,j))
            CALL acc(a_tdtr      (i,j), 1., f_tdtr      (i,j))
            CALL acc(a_qref      (i,j), 1., f_qref      (i,j))
            CALL acc(a_xy_rain   (i,j), 1., f_xy_rain   (i,j))
            CALL acc(a_xy_snow   (i,j), 1., f_xy_snow   (i,j))

            CALL acc(a_t_room    (i,j), 1., f_t_room    (i,j))
            CALL acc(a_tafu      (i,j), 1., f_tafu      (i,j))
            CALL acc(a_fhac      (i,j), 1., f_fhac      (i,j))
            CALL acc(a_fwst      (i,j), 1., f_fwst      (i,j))
            CALL acc(a_fach      (i,j), 1., f_fach      (i,j))
            CALL acc(a_fahe      (i,j), 1., f_fahe      (i,j))
            CALL acc(a_fhah      (i,j), 1., f_fhah      (i,j))
            CALL acc(a_vehc      (i,j), 1., f_vehc      (i,j))
            CALL acc(a_meta      (i,j), 1., f_meta      (i,j))

            CALL acc(a_senroof   (i,j), 1., f_senroof   (i,j))
            CALL acc(a_senwsun   (i,j), 1., f_senwsun   (i,j))
            CALL acc(a_senwsha   (i,j), 1., f_senwsha   (i,j))
            CALL acc(a_sengimp   (i,j), 1., f_sengimp   (i,j))
            CALL acc(a_sengper   (i,j), 1., f_sengper   (i,j))
            CALL acc(a_senurl    (i,j), 1., f_senurl    (i,j))

            CALL acc(a_lfevproof (i,j), 1., f_lfevproof (i,j))
            CALL acc(a_lfevpgimp (i,j), 1., f_lfevpgimp (i,j))
            CALL acc(a_lfevpgper (i,j), 1., f_lfevpgper (i,j))
            CALL acc(a_lfevpurl  (i,j), 1., f_lfevpurl  (i,j))

            CALL acc(a_troof     (i,j), 1., f_troof     (i,j))
            CALL acc(a_twall     (i,j), 1., f_twall     (i,j))

            CALL acc(a_sabvdt    (i,j), 1., f_sabvdt    (i,j))
            CALL acc(a_sabgdt    (i,j), 1., f_sabgdt    (i,j))
            CALL acc(a_srdt      (i,j), 1., f_srdt      (i,j))
            CALL acc(a_fsenadt   (i,j), 1., f_fsenadt   (i,j))
            CALL acc(a_lfevpadt  (i,j), 1., f_lfevpadt  (i,j))
            CALL acc(a_fgrnddt   (i,j), 1., f_fgrnddt   (i,j))
            CALL acc(a_olrgdt    (i,j), 1., f_olrgdt    (i,j))
            CALL acc(a_rnetdt    (i,j), 1., f_rnetdt    (i,j))
            CALL acc(a_t_grnddt  (i,j), 1., f_t_grnddt  (i,j))
            CALL acc(a_traddt    (i,j), 1., f_traddt    (i,j))
            CALL acc(a_trefdt    (i,j), 1., f_trefdt    (i,j))
            CALL acc(a_tafudt    (i,j), 1., f_tafudt    (i,j))

            CALL acc(a_fsenant   (i,j), 1., f_fsenant   (i,j))
            CALL acc(a_lfevpant  (i,j), 1., f_lfevpant  (i,j))
            CALL acc(a_fgrndnt   (i,j), 1., f_fgrndnt   (i,j))
            CALL acc(a_olrgnt    (i,j), 1., f_olrgnt    (i,j))
            CALL acc(a_rnetnt    (i,j), 1., f_rnetnt    (i,j))
            CALL acc(a_t_grndnt  (i,j), 1., f_t_grndnt  (i,j))
            CALL acc(a_tradnt    (i,j), 1., f_tradnt    (i,j))
            CALL acc(a_trefnt    (i,j), 1., f_trefnt    (i,j))
            CALL acc(a_tafunt    (i,j), 1., f_tafunt    (i,j))

            DO l = maxsnl+1, nl_soil
               CALL acc(a_t_soisno    (l,i,j), 1., f_t_soisno    (l,i,j))
               CALL acc(a_wliq_soisno (l,i,j), 1., f_wliq_soisno (l,i,j))
               CALL acc(a_wice_soisno (l,i,j), 1., f_wice_soisno (l,i,j))
            ENDDO

            DO l = 1, nl_soil
               CALL acc(a_h2osoi      (l,i,j), 1., f_h2osoi      (l,i,j))
            ENDDO

            DO l = 1, nl_lake
               CALL acc(a_t_lake      (l,i,j), 1., f_t_lake      (l,i,j))
               CALL acc(a_lake_icefrac(l,i,j), 1., f_lake_icefrac(l,i,j))
            ENDDO

            CALL acc(a_ustar     (i,j), 1., f_ustar     (i,j))
            CALL acc(a_tstar     (i,j), 1., f_tstar     (i,j))
            CALL acc(a_qstar     (i,j), 1., f_qstar     (i,j))
            CALL acc(a_zol       (i,j), 1., f_zol       (i,j))
            CALL acc(a_rib       (i,j), 1., f_rib       (i,j))
            CALL acc(a_fm        (i,j), 1., f_fm        (i,j))
            CALL acc(a_fh        (i,j), 1., f_fh        (i,j))
            CALL acc(a_fq        (i,j), 1., f_fq        (i,j))

            CALL acc(a_us10m     (i,j), 1., f_us10m     (i,j))
            CALL acc(a_vs10m     (i,j), 1., f_vs10m     (i,j))
            CALL acc(a_fm10m     (i,j), 1., f_fm10m     (i,j))

            CALL acc(a_sr        (i,j), 1., f_sr        (i,j))
            CALL acc(a_solvd     (i,j), 1., f_solvd     (i,j))
            CALL acc(a_solvi     (i,j), 1., f_solvi     (i,j))
            CALL acc(a_solnd     (i,j), 1., f_solnd     (i,j))
            CALL acc(a_solni     (i,j), 1., f_solni     (i,j))
            CALL acc(a_srvd      (i,j), 1., f_srvd      (i,j))
            CALL acc(a_srvi      (i,j), 1., f_srvi      (i,j))
            CALL acc(a_srnd      (i,j), 1., f_srnd      (i,j))
            CALL acc(a_srni      (i,j), 1., f_srni      (i,j))
            CALL acc(a_solvdln   (i,j), 1., f_solvdln   (i,j))
            CALL acc(a_solviln   (i,j), 1., f_solviln   (i,j))
            CALL acc(a_solndln   (i,j), 1., f_solndln   (i,j))
            CALL acc(a_solniln   (i,j), 1., f_solniln   (i,j))
            CALL acc(a_srvdln    (i,j), 1., f_srvdln    (i,j))
            CALL acc(a_srviln    (i,j), 1., f_srviln    (i,j))
            CALL acc(a_srndln    (i,j), 1., f_srndln    (i,j))
            CALL acc(a_srniln    (i,j), 1., f_srniln    (i,j))

            IF (a_solvdln(i,j) /= spval) nac_ln(i,j) = nac_ln(i,j) + 1
            IF (a_trefdt (i,j) /= spval) nac_dt(i,j) = nac_dt(i,j) + 1
            IF (a_trefnt (i,j) /= spval) nac_nt(i,j) = nac_nt(i,j) + 1

         ENDDO
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif
   END SUBROUTINE vec2xy


   SUBROUTINE allocate_vec2xy

      IMPLICIT NONE

      !---------------------------------------------------------------------
      allocate (a_xy_us     (lon_points,lat_points))
      allocate (a_xy_vs     (lon_points,lat_points))
      allocate (a_xy_t      (lon_points,lat_points))
      allocate (a_xy_q      (lon_points,lat_points))
      allocate (a_xy_prc    (lon_points,lat_points))
      allocate (a_xy_prl    (lon_points,lat_points))
      allocate (a_xy_pbot   (lon_points,lat_points))
      allocate (a_xy_frl    (lon_points,lat_points))
      allocate (a_xy_solarin(lon_points,lat_points))
      allocate (a_xy_rain   (lon_points,lat_points))
      allocate (a_xy_snow   (lon_points,lat_points))

      !---------------------------------------------------------------------
      allocate (a_taux      (lon_points,lat_points))
      allocate (a_tauy      (lon_points,lat_points))
      allocate (a_fsena     (lon_points,lat_points))
      allocate (a_lfevpa    (lon_points,lat_points))
      allocate (a_fevpa     (lon_points,lat_points))
      allocate (a_fsenl     (lon_points,lat_points))
      allocate (a_fevpl     (lon_points,lat_points))
      allocate (a_etr       (lon_points,lat_points))
      allocate (a_fseng     (lon_points,lat_points))
      allocate (a_fevpg     (lon_points,lat_points))
      allocate (a_fgrnd     (lon_points,lat_points))
      allocate (a_sabvsun   (lon_points,lat_points))
      allocate (a_sabvsha   (lon_points,lat_points))
      allocate (a_sabg      (lon_points,lat_points))
      allocate (a_olrg      (lon_points,lat_points))
      allocate (a_rnet      (lon_points,lat_points))
      allocate (a_xerr      (lon_points,lat_points))
      allocate (a_zerr      (lon_points,lat_points))
      allocate (a_rsur      (lon_points,lat_points))
      allocate (a_qintr     (lon_points,lat_points))
      allocate (a_qinfl     (lon_points,lat_points))
      allocate (a_qdrip     (lon_points,lat_points))

      allocate (a_assim     (lon_points,lat_points))
      allocate (a_respc     (lon_points,lat_points))
      allocate (a_qcharge   (lon_points,lat_points))

      !---------------------------------------------------------------------
      allocate (a_t_grnd    (lon_points,lat_points))
      allocate (a_tleaf     (lon_points,lat_points))
      allocate (a_ldew      (lon_points,lat_points))
      allocate (a_scv       (lon_points,lat_points))
      allocate (a_snowdp    (lon_points,lat_points))
      allocate (a_fsno      (lon_points,lat_points))
      allocate (a_sigf      (lon_points,lat_points))
      allocate (a_green     (lon_points,lat_points))
      allocate (a_lai       (lon_points,lat_points))
      allocate (a_laisun    (lon_points,lat_points))
      allocate (a_laisha    (lon_points,lat_points))
      allocate (a_sai       (lon_points,lat_points))
      allocate (a_alb   (2,2,lon_points,lat_points))
      allocate (a_emis      (lon_points,lat_points))
      allocate (a_z0m       (lon_points,lat_points))
      allocate (a_trad      (lon_points,lat_points))
      allocate (a_tref      (lon_points,lat_points))
      allocate (a_tmax      (lon_points,lat_points))
      allocate (a_tmin      (lon_points,lat_points))
      allocate (a_tdtr      (lon_points,lat_points))
      allocate (a_qref      (lon_points,lat_points))

      !---------------------------------------------------------------------
      allocate (a_t_room    (lon_points,lat_points))
      allocate (a_tafu      (lon_points,lat_points))
      allocate (a_fhac      (lon_points,lat_points))
      allocate (a_fwst      (lon_points,lat_points))
      allocate (a_fach      (lon_points,lat_points))
      allocate (a_fahe      (lon_points,lat_points))
      allocate (a_fhah      (lon_points,lat_points))
      allocate (a_vehc      (lon_points,lat_points))
      allocate (a_meta      (lon_points,lat_points))

      allocate (a_senroof   (lon_points,lat_points))
      allocate (a_senwsun   (lon_points,lat_points))
      allocate (a_senwsha   (lon_points,lat_points))
      allocate (a_sengimp   (lon_points,lat_points))
      allocate (a_sengper   (lon_points,lat_points))
      allocate (a_senurl    (lon_points,lat_points))

      allocate (a_lfevproof (lon_points,lat_points))
      allocate (a_lfevpgimp (lon_points,lat_points))
      allocate (a_lfevpgper (lon_points,lat_points))
      allocate (a_lfevpurl  (lon_points,lat_points))

      allocate (a_troof     (lon_points,lat_points))
      allocate (a_twall     (lon_points,lat_points))

      allocate (a_sabvdt    (lon_points,lat_points))
      allocate (a_sabgdt    (lon_points,lat_points))
      allocate (a_srdt      (lon_points,lat_points))
      allocate (a_fsenadt   (lon_points,lat_points))
      allocate (a_lfevpadt  (lon_points,lat_points))
      allocate (a_fgrnddt   (lon_points,lat_points))
      allocate (a_olrgdt    (lon_points,lat_points))
      allocate (a_rnetdt    (lon_points,lat_points))
      allocate (a_t_grnddt  (lon_points,lat_points))
      allocate (a_traddt    (lon_points,lat_points))
      allocate (a_trefdt    (lon_points,lat_points))
      allocate (a_tafudt    (lon_points,lat_points))

      allocate (a_fsenant   (lon_points,lat_points))
      allocate (a_lfevpant  (lon_points,lat_points))
      allocate (a_fgrndnt   (lon_points,lat_points))
      allocate (a_olrgnt    (lon_points,lat_points))
      allocate (a_rnetnt    (lon_points,lat_points))
      allocate (a_t_grndnt  (lon_points,lat_points))
      allocate (a_tradnt    (lon_points,lat_points))
      allocate (a_trefnt    (lon_points,lat_points))
      allocate (a_tafunt    (lon_points,lat_points))

      !---------------------------------------------------------------------
      allocate (a_t_soisno   (maxsnl+1:nl_soil,lon_points,lat_points))
      allocate (a_wliq_soisno(maxsnl+1:nl_soil,lon_points,lat_points))
      allocate (a_wice_soisno(maxsnl+1:nl_soil,lon_points,lat_points))
      allocate (a_h2osoi            (1:nl_soil,lon_points,lat_points))
      allocate (a_rstfac                      (lon_points,lat_points))
      allocate (a_zwt                         (lon_points,lat_points))
      allocate (a_wa                          (lon_points,lat_points))
      allocate (a_wat                         (lon_points,lat_points))

      allocate (a_t_lake              (nl_lake,lon_points,lat_points))
      allocate (a_lake_icefrac        (nl_lake,lon_points,lat_points))

      !---------------------------------------------------------------------
      allocate (a_ustar     (lon_points,lat_points))
      allocate (a_tstar     (lon_points,lat_points))
      allocate (a_qstar     (lon_points,lat_points))
      allocate (a_zol       (lon_points,lat_points))
      allocate (a_rib       (lon_points,lat_points))
      allocate (a_fm        (lon_points,lat_points))
      allocate (a_fh        (lon_points,lat_points))
      allocate (a_fq        (lon_points,lat_points))

      allocate (a_us10m     (lon_points,lat_points))
      allocate (a_vs10m     (lon_points,lat_points))
      allocate (a_fm10m     (lon_points,lat_points))

      !---------------------------------------------------------------------
      allocate (a_sr        (lon_points,lat_points))
      allocate (a_solvd     (lon_points,lat_points))
      allocate (a_solvi     (lon_points,lat_points))
      allocate (a_solnd     (lon_points,lat_points))
      allocate (a_solni     (lon_points,lat_points))
      allocate (a_srvd      (lon_points,lat_points))
      allocate (a_srvi      (lon_points,lat_points))
      allocate (a_srnd      (lon_points,lat_points))
      allocate (a_srni      (lon_points,lat_points))
      allocate (a_solvdln   (lon_points,lat_points))
      allocate (a_solviln   (lon_points,lat_points))
      allocate (a_solndln   (lon_points,lat_points))
      allocate (a_solniln   (lon_points,lat_points))
      allocate (a_srvdln    (lon_points,lat_points))
      allocate (a_srviln    (lon_points,lat_points))
      allocate (a_srndln    (lon_points,lat_points))
      allocate (a_srniln    (lon_points,lat_points))

   END SUBROUTINE allocate_vec2xy


   SUBROUTINE deallocate_vec2xy

      IMPLICIT NONE

      !---------------------------------------------------------------------
      deallocate ( a_xy_us        )
      deallocate ( a_xy_vs        )
      deallocate ( a_xy_t         )
      deallocate ( a_xy_q         )
      deallocate ( a_xy_prc       )
      deallocate ( a_xy_prl       )
      deallocate ( a_xy_pbot      )
      deallocate ( a_xy_frl       )
      deallocate ( a_xy_solarin   )
      deallocate ( a_xy_rain      )
      deallocate ( a_xy_snow      )

      !---------------------------------------------------------------------
      deallocate ( a_taux         )
      deallocate ( a_tauy         )
      deallocate ( a_fsena        )
      deallocate ( a_lfevpa       )
      deallocate ( a_fevpa        )
      deallocate ( a_fsenl        )
      deallocate ( a_fevpl        )
      deallocate ( a_etr          )
      deallocate ( a_fseng        )
      deallocate ( a_fevpg        )
      deallocate ( a_fgrnd        )
      deallocate ( a_sabvsun      )
      deallocate ( a_sabvsha      )
      deallocate ( a_sabg         )
      deallocate ( a_olrg         )
      deallocate ( a_rnet         )
      deallocate ( a_xerr         )
      deallocate ( a_zerr         )
      deallocate ( a_rsur         )
      deallocate ( a_qintr        )
      deallocate ( a_qinfl        )
      deallocate ( a_qdrip        )

      deallocate ( a_assim        )
      deallocate ( a_respc        )
      deallocate ( a_qcharge      )

      !---------------------------------------------------------------------
      deallocate ( a_t_grnd       )
      deallocate ( a_tleaf        )
      deallocate ( a_ldew         )
      deallocate ( a_scv          )
      deallocate ( a_snowdp       )
      deallocate ( a_fsno         )
      deallocate ( a_sigf         )
      deallocate ( a_green        )
      deallocate ( a_lai          )
      deallocate ( a_laisun       )
      deallocate ( a_laisha       )
      deallocate ( a_sai          )
      deallocate ( a_alb          )
      deallocate ( a_emis         )
      deallocate ( a_z0m          )
      deallocate ( a_trad         )
      deallocate ( a_tref         )
      deallocate ( a_tmax         )
      deallocate ( a_tmin         )
      deallocate ( a_tdtr         )
      deallocate ( a_qref         )

      !---------------------------------------------------------------------
      deallocate ( a_t_room       )
      deallocate ( a_tafu         )
      deallocate ( a_fhac         )
      deallocate ( a_fwst         )
      deallocate ( a_fach         )
      deallocate ( a_fahe         )
      deallocate ( a_fhah         )
      deallocate ( a_vehc         )
      deallocate ( a_meta         )

      deallocate ( a_senroof      )
      deallocate ( a_senwsun      )
      deallocate ( a_senwsha      )
      deallocate ( a_sengimp      )
      deallocate ( a_sengper      )
      deallocate ( a_senurl       )

      deallocate ( a_lfevproof    )
      deallocate ( a_lfevpgimp    )
      deallocate ( a_lfevpgper    )
      deallocate ( a_lfevpurl     )

      deallocate ( a_troof        )
      deallocate ( a_twall        )

      deallocate ( a_sabvdt       )
      deallocate ( a_sabgdt       )
      deallocate ( a_srdt         )
      deallocate ( a_fsenadt      )
      deallocate ( a_lfevpadt     )
      deallocate ( a_fgrnddt      )
      deallocate ( a_olrgdt       )
      deallocate ( a_rnetdt       )
      deallocate ( a_t_grnddt     )
      deallocate ( a_traddt       )
      deallocate ( a_trefdt       )
      deallocate ( a_tafudt       )

      deallocate ( a_fsenant      )
      deallocate ( a_lfevpant     )
      deallocate ( a_fgrndnt      )
      deallocate ( a_olrgnt       )
      deallocate ( a_rnetnt       )
      deallocate ( a_t_grndnt     )
      deallocate ( a_tradnt       )
      deallocate ( a_trefnt       )
      deallocate ( a_tafunt       )

      !---------------------------------------------------------------------
      deallocate ( a_t_soisno     )
      deallocate ( a_wliq_soisno  )
      deallocate ( a_wice_soisno  )
      deallocate ( a_h2osoi       )
      deallocate ( a_rstfac       )
      deallocate ( a_zwt          )
      deallocate ( a_wa           )
      deallocate ( a_wat          )

      deallocate ( a_t_lake       )
      deallocate ( a_lake_icefrac )

      !---------------------------------------------------------------------
      deallocate ( a_ustar        )
      deallocate ( a_tstar        )
      deallocate ( a_qstar        )
      deallocate ( a_zol          )
      deallocate ( a_rib          )
      deallocate ( a_fm           )
      deallocate ( a_fh           )
      deallocate ( a_fq           )

      deallocate ( a_us10m        )
      deallocate ( a_vs10m        )
      deallocate ( a_fm10m        )

      !---------------------------------------------------------------------
      deallocate ( a_sr           )
      deallocate ( a_solvd        )
      deallocate ( a_solvi        )
      deallocate ( a_solnd        )
      deallocate ( a_solni        )
      deallocate ( a_srvd         )
      deallocate ( a_srvi         )
      deallocate ( a_srnd         )
      deallocate ( a_srni         )
      deallocate ( a_solvdln      )
      deallocate ( a_solviln      )
      deallocate ( a_solndln      )
      deallocate ( a_solniln      )
      deallocate ( a_srvdln       )
      deallocate ( a_srviln       )
      deallocate ( a_srndln       )
      deallocate ( a_srniln       )

   END SUBROUTINE deallocate_vec2xy


   SUBROUTINE acc(var, wgt, s)

      USE precision
      USE GlobalVars

      IMPLICIT NONE

      REAL(r8), intent(in)  :: var
      REAL(r8), intent(in)  :: wgt
      REAL(r8), intent(out) :: s

      IF (var /= spval) THEN
         IF (s /= spval) THEN
            s = s + wgt*var
         ELSE
            s = wgt*var
         ENDIF
      ENDIF

   END SUBROUTINE acc

END MODULE MOD_vec2xy
! ----------------------------------------------------------------------
! EOP
