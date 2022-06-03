#include <define.h>

!TODO: 改成module,动态分配内存,SAVE
SUBROUTINE vec2xy (nac,nac_ln,nac_dt,nac_nt,a_rnof)
! ----------------------------------------------------------------------
! perfrom the grid average mapping: average a subgrid input 1d vector
! of length numpatch to a output 2d array of length [lon_points,lat_points]
!
! Created by Yongjiu Dai, 03/2014
!---------------------------------------------------------------------

USE precision
USE GlobalVars
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

INTEGER, intent(inout) :: nac
INTEGER, intent(inout) :: nac_ln(lon_points,lat_points)
INTEGER, intent(inout) :: nac_dt(lon_points,lat_points)
INTEGER, intent(inout) :: nac_nt(lon_points,lat_points)

!---------------------------------------------------------------------
REAL(r8) a_xy_us  (lon_points,lat_points)  ! wind in eastward direction [m/s]
REAL(r8) a_xy_vs  (lon_points,lat_points)  ! wind in northward direction [m/s]
REAL(r8) a_xy_t   (lon_points,lat_points)  ! temperature at reference height [kelvin]
REAL(r8) a_xy_q   (lon_points,lat_points)  ! specific humidity at reference height [kg/kg]
REAL(r8) a_xy_prc (lon_points,lat_points)  ! convective precipitation [mm/s]
REAL(r8) a_xy_prl (lon_points,lat_points)  ! large scale precipitation [mm/s]
REAL(r8) a_xy_pbot(lon_points,lat_points)  ! atmospheric pressure at the surface [pa]
REAL(r8) a_xy_frl (lon_points,lat_points)  ! atmospheric infrared (longwave) radiation [W/m2]
REAL(r8) a_xy_solarin(lon_points,lat_points)! downward solar radiation at surface [W/m2]
REAL(r8) a_xy_rain(lon_points,lat_points)  ! rain [mm/s]
REAL(r8) a_xy_snow(lon_points,lat_points)  ! snow [mm/s]

!---------------------------------------------------------------------
REAL(r8) a_taux   (lon_points,lat_points)  ! wind stress: E-W [kg/m/s2]
REAL(r8) a_tauy   (lon_points,lat_points)  ! wind stress: N-S [kg/m/s2]
REAL(r8) a_fsena  (lon_points,lat_points)  ! sensible heat from canopy height to atmosphere [W/m2]
REAL(r8) a_lfevpa (lon_points,lat_points)  ! latent heat flux from canopy height to atmosphere [W/m2]
REAL(r8) a_fevpa  (lon_points,lat_points)  ! evapotranspiration from canopy to atmosphere [mm/s]
REAL(r8) a_fsenl  (lon_points,lat_points)  ! sensible heat from leaves [W/m2]
REAL(r8) a_fevpl  (lon_points,lat_points)  ! evaporation+transpiration from leaves [mm/s]
REAL(r8) a_etr    (lon_points,lat_points)  ! transpiration rate [mm/s]
REAL(r8) a_fseng  (lon_points,lat_points)  ! sensible heat flux from ground [W/m2]
REAL(r8) a_fevpg  (lon_points,lat_points)  ! evaporation heat flux from ground [mm/s]
REAL(r8) a_fgrnd  (lon_points,lat_points)  ! ground heat flux [W/m2]
REAL(r8) a_sabvsun(lon_points,lat_points)  ! solar absorbed by sunlit canopy [W/m2]
REAL(r8) a_sabvsha(lon_points,lat_points)  ! solar absorbed by shaded [W/m2]
REAL(r8) a_sabg   (lon_points,lat_points)  ! solar absorbed by ground  [W/m2]
REAL(r8) a_olrg   (lon_points,lat_points)  ! outgoing long-wave radiation from ground+canopy [W/m2]
REAL(r8) a_rnet   (lon_points,lat_points)  ! net radiation [W/m2]
REAL(r8) a_xerr   (lon_points,lat_points)  ! the error of water banace [mm/s]
REAL(r8) a_zerr   (lon_points,lat_points)  ! the error of energy balance [W/m2]
REAL(r8) a_rsur   (lon_points,lat_points)  ! surface runoff [mm/s]
REAL(r8) a_rnof   (lon_points,lat_points)  ! total runoff [mm/s]
REAL(r8) a_qintr  (lon_points,lat_points)  ! interception [mm/s]
REAL(r8) a_qinfl  (lon_points,lat_points)  ! inflitration [mm/s]
REAL(r8) a_qdrip  (lon_points,lat_points)  ! throughfall [mm/s]

REAL(r8) a_assim  (lon_points,lat_points)  ! canopy assimilation rate [mol m-2 s-1]
REAL(r8) a_respc  (lon_points,lat_points)  ! respiration (plant+soil) [mol m-2 s-1]
REAL(r8) a_qcharge(lon_points,lat_points)  ! groundwater recharge rate [mm/s]

!---------------------------------------------------------------------
REAL(r8) a_t_grnd (lon_points,lat_points)  ! ground surface temperature [K]
REAL(r8) a_tleaf  (lon_points,lat_points)  ! sunlit leaf temperature [K]
REAL(r8) a_ldew   (lon_points,lat_points)  ! depth of water on foliage [mm]
REAL(r8) a_scv    (lon_points,lat_points)  ! snow cover, water equivalent [mm]
REAL(r8) a_snowdp (lon_points,lat_points)  ! snow depth [meter]
REAL(r8) a_fsno   (lon_points,lat_points)  ! fraction of snow cover on ground
REAL(r8) a_sigf   (lon_points,lat_points)  ! fraction of veg cover, excluding snow-covered veg [-]
REAL(r8) a_green  (lon_points,lat_points)  ! leaf greenness
REAL(r8) a_lai    (lon_points,lat_points)  ! leaf area index
REAL(r8) a_laisun (lon_points,lat_points)  ! sunlit leaf area index
REAL(r8) a_laisha (lon_points,lat_points)  ! shaded leaf area index
REAL(r8) a_sai    (lon_points,lat_points)  ! stem area index
REAL(r8) a_alb(2,2,lon_points,lat_points)  ! averaged albedo [visible, direct; direct, diffuse]
REAL(r8) a_emis   (lon_points,lat_points)  ! averaged bulk surface emissivity
REAL(r8) a_z0m    (lon_points,lat_points)  ! effective roughness [m]
REAL(r8) a_trad   (lon_points,lat_points)  ! radiative temperature of surface [K]
REAL(r8) a_tref   (lon_points,lat_points)  ! 2 m height air temperature [kelvin]
REAL(r8) a_qref   (lon_points,lat_points)  ! 2 m height air specific humidity [kg/kg]

!---------------------------------------------------------------------
REAL(r8) a_t_room (lon_points,lat_points)  ! temperature of inner building [K]
REAL(r8) a_tafu   (lon_points,lat_points)  ! temperature of outer building [K]
REAL(r8) a_fhac   (lon_points,lat_points)  ! sensible flux from heat or cool AC [W/m2]
REAL(r8) a_fwst   (lon_points,lat_points)  ! waste heat flux from heat or cool AC [W/m2]
REAL(r8) a_fach   (lon_points,lat_points)  ! flux from inner and outter air exchange [W/m2]

REAL(r8) a_sabvdt  (lon_points,lat_points) ! solar absorbed by sunlit canopy [W/m2]
REAL(r8) a_sabgdt  (lon_points,lat_points) ! solar absorbed by ground [W/m2]
REAL(r8) a_srdt    (lon_points,lat_points) ! total reflected solar radiation (W/m2)
REAL(r8) a_fsenadt (lon_points,lat_points) ! sensible heat from canopy height to atmosphere [W/m2]
REAL(r8) a_lfevpadt(lon_points,lat_points) ! latent heat flux from canopy height to atmosphere [W/m2]
REAL(r8) a_fgrnddt (lon_points,lat_points) ! ground heat flux [W/m2]
REAL(r8) a_olrgdt  (lon_points,lat_points) ! outgoing long-wave radiation from ground+canopy [W/m2]
REAL(r8) a_rnetdt  (lon_points,lat_points) ! net radiation [W/m2]
REAL(r8) a_t_grnddt(lon_points,lat_points) ! ground surface temperature [K]
REAL(r8) a_traddt  (lon_points,lat_points) ! radiative temperature of surface [K]
REAL(r8) a_trefdt  (lon_points,lat_points) ! 2 m height air temperature [kelvin]
REAL(r8) a_tafudt  (lon_points,lat_points) ! temperature of outer building [K]

REAL(r8) a_fsenant (lon_points,lat_points) ! sensible heat from canopy height to atmosphere [W/m2]
REAL(r8) a_lfevpant(lon_points,lat_points) ! latent heat flux from canopy height to atmosphere [W/m2]
REAL(r8) a_fgrndnt (lon_points,lat_points) ! ground heat flux [W/m2]
REAL(r8) a_olrgnt  (lon_points,lat_points) ! outgoing long-wave radiation from ground+canopy [W/m2]
REAL(r8) a_rnetnt  (lon_points,lat_points) ! net radiation [W/m2]
REAL(r8) a_t_grndnt(lon_points,lat_points) ! ground surface temperature [K]
REAL(r8) a_tradnt  (lon_points,lat_points) ! radiative temperature of surface [K]
REAL(r8) a_trefnt  (lon_points,lat_points) ! 2 m height air temperature [kelvin]
REAL(r8) a_tafunt  (lon_points,lat_points) ! temperature of outer building [K]

!---------------------------------------------------------------------
REAL(r8) a_t_soisno   (maxsnl+1:nl_soil,lon_points,lat_points)  ! soil temperature [K]
REAL(r8) a_wliq_soisno(maxsnl+1:nl_soil,lon_points,lat_points)  ! liquid water in soil layers [kg/m2]
REAL(r8) a_wice_soisno(maxsnl+1:nl_soil,lon_points,lat_points)  ! ice lens in soil layers [kg/m2]
REAL(r8) a_h2osoi            (1:nl_soil,lon_points,lat_points)  ! volumetric soil water in layers [m3/m3]
REAL(r8) a_rstfac                      (lon_points,lat_points)  ! factor of soil water stress
REAL(r8) a_zwt                         (lon_points,lat_points)  ! the depth to water table [m]
REAL(r8) a_wa                          (lon_points,lat_points)  ! water storage in aquifer [mm]
REAL(r8) a_wat                         (lon_points,lat_points)  ! total water storage [mm]

REAL(r8) a_t_lake      (nl_lake,lon_points,lat_points) ! lake temperature [K]
REAL(r8) a_lake_icefrac(nl_lake,lon_points,lat_points) ! lake ice fraction cover [0-1]

!---------------------------------------------------------------------
REAL(r8) a_ustar  (lon_points,lat_points)  ! u* in similarity theory [m/s]
REAL(r8) a_tstar  (lon_points,lat_points)  ! t* in similarity theory [kg/kg]
REAL(r8) a_qstar  (lon_points,lat_points)  ! q* in similarity theory [kg/kg]
REAL(r8) a_zol    (lon_points,lat_points)  ! dimensionless height (z/L) used in Monin-Obukhov theory
REAL(r8) a_rib    (lon_points,lat_points)  ! bulk Richardson number in surface layer
REAL(r8) a_fm     (lon_points,lat_points)  ! integral of profile function for momentum
REAL(r8) a_fh     (lon_points,lat_points)  ! integral of profile function for heat
REAL(r8) a_fq     (lon_points,lat_points)  ! integral of profile function for moisture

REAL(r8) a_us10m  (lon_points,lat_points)  ! 10m u-velocity [m/s]
REAL(r8) a_vs10m  (lon_points,lat_points)  ! 10m v-velocity [m/s]
REAL(r8) a_fm10m  (lon_points,lat_points)  ! integral of profile function for momentum at 10m [-]

!---------------------------------------------------------------------
REAL(r8) a_sr     (lon_points,lat_points)  ! total reflected solar radiation (W/m2)
REAL(r8) a_solvd  (lon_points,lat_points)  ! incident direct beam vis solar radiation (W/m2)
REAL(r8) a_solvi  (lon_points,lat_points)  ! incident diffuse beam vis solar radiation (W/m2)
REAL(r8) a_solnd  (lon_points,lat_points)  ! incident direct beam nir solar radiation (W/m2)
REAL(r8) a_solni  (lon_points,lat_points)  ! incident diffuse beam nir solar radiation (W/m2)
REAL(r8) a_srvd   (lon_points,lat_points)  ! reflected direct beam vis solar radiation (W/m2)
REAL(r8) a_srvi   (lon_points,lat_points)  ! reflected diffuse beam vis solar radiation (W/m2)
REAL(r8) a_srnd   (lon_points,lat_points)  ! reflected direct beam nir solar radiation (W/m2)
REAL(r8) a_srni   (lon_points,lat_points)  ! reflected diffuse beam nir solar radiation (W/m2)
REAL(r8) a_solvdln(lon_points,lat_points)  ! incident direct beam vis solar radiation at local noon (W/m2)
REAL(r8) a_solviln(lon_points,lat_points)  ! incident diffuse beam vis solar radiation at local noon (W/m2)
REAL(r8) a_solndln(lon_points,lat_points)  ! incident direct beam nir solar radiation at local noon (W/m2)
REAL(r8) a_solniln(lon_points,lat_points)  ! incident diffuse beam nir solar radiation at local noon (W/m2)
REAL(r8) a_srvdln (lon_points,lat_points)  ! reflected direct beam vis solar radiation at local noon (W/m2)
REAL(r8) a_srviln (lon_points,lat_points)  ! reflected diffuse beam vis solar radiation at local noon (W/m2)
REAL(r8) a_srndln (lon_points,lat_points)  ! reflected direct beam nir solar radiation at local noon (W/m2)
REAL(r8) a_srniln (lon_points,lat_points)  ! reflected diffuse beam nir solar radiation at local noon (W/m2)

!---------------------------------------------------------------------
! local variables

      INTEGER  i,j,np,u,l
      REAL(r8) sumwt(lon_points,lat_points)
      REAL(r8) urbwt(lon_points,lat_points)
      REAL(r8) rhoair,thm,th,thv,ur,displa_av,zldis,hgt_u,hgt_t,hgt_q
      REAL(r8) z0m_av,z0h_av,z0q_av,us,vs,tm,qm,psrf
      REAL(r8) obu,fh2m,fq2m
      REAL(r8) um,thvstar,beta,zii,wc,wc2

! ---------------------------------------------------
! Meteorological forcing
! ---------------------------------------------------
      a_xy_us     (:,:) = forc_xy_us     (:,:)
      a_xy_vs     (:,:) = forc_xy_vs     (:,:)
      a_xy_t      (:,:) = forc_xy_t      (:,:)
      a_xy_q      (:,:) = forc_xy_q      (:,:)
      a_xy_prc    (:,:) = forc_xy_prc    (:,:)
      a_xy_prl    (:,:) = forc_xy_prl    (:,:)
      a_xy_pbot   (:,:) = forc_xy_pbot   (:,:)
      a_xy_frl    (:,:) = forc_xy_frl    (:,:)

      a_xy_solarin(:,:) = forc_xy_sols (:,:) + forc_xy_soll (:,:) &
                        + forc_xy_solsd(:,:) + forc_xy_solld(:,:)

! ------------------------------------------------------------------------------------------
! Mapping the fluxes and state variables at patch [numpatch] to grid [lon_points,lat_points]
! ------------------------------------------------------------------------------------------
      sumwt    (:,:) = 0.
      urbwt    (:,:) = 0.
      a_taux   (:,:) = 0.
      a_tauy   (:,:) = 0.
      a_fsena  (:,:) = 0.
      a_lfevpa (:,:) = 0.
      a_fevpa  (:,:) = 0.
      a_fsenl  (:,:) = 0.
      a_fevpl  (:,:) = 0.
      a_etr    (:,:) = 0.
      a_fseng  (:,:) = 0.
      a_fevpg  (:,:) = 0.
      a_fgrnd  (:,:) = 0.
      a_sabvsun(:,:) = 0.
      a_sabvsha(:,:) = 0.
      a_sabg   (:,:) = 0.
      a_olrg   (:,:) = 0.
      a_rnet   (:,:) = 0.
      a_xerr   (:,:) = 0.
      a_zerr   (:,:) = 0.
      a_rsur   (:,:) = 0.
      a_rnof   (:,:) = 0.
      a_qintr  (:,:) = 0.
      a_qinfl  (:,:) = 0.
      a_qdrip  (:,:) = 0.
      a_wat    (:,:) = 0.
      a_assim  (:,:) = 0.
      a_respc  (:,:) = 0.

      a_qcharge(:,:) = 0.
      a_t_grnd (:,:) = 0.
      a_tleaf  (:,:) = 0.
      a_ldew   (:,:) = 0.
      a_scv    (:,:) = 0.
      a_snowdp (:,:) = 0.
      a_fsno   (:,:) = 0.
      a_sigf   (:,:) = 0.
      a_green  (:,:) = 0.
      a_lai    (:,:) = 0.
      a_laisun (:,:) = 0.
      a_laisha (:,:) = 0.
      a_sai    (:,:) = 0.
      a_alb(:,:,:,:) = 0.
      a_emis   (:,:) = 0.
      a_z0m    (:,:) = 0.
      a_trad   (:,:) = 0.
      a_tref   (:,:) = 0.
      a_qref   (:,:) = 0.
      a_xy_rain(:,:) = 0.
      a_xy_snow(:,:) = 0.

      a_t_room (:,:) = 0.
      a_tafu   (:,:) = 0.
      a_fhac   (:,:) = 0.
      a_fwst   (:,:) = 0.
      a_fach   (:,:) = 0.

!TODO: 添加a_tmax, a_tmin, a_tavg, a_tdtr

      a_sabvdt  (:,:) = spval
      a_sabgdt  (:,:) = spval
      a_srdt    (:,:) = spval
      a_fsenadt (:,:) = spval
      a_lfevpadt(:,:) = spval
      a_fgrnddt (:,:) = spval
      a_olrgdt  (:,:) = spval
      a_rnetdt  (:,:) = spval
      a_t_grnddt(:,:) = spval
      a_traddt  (:,:) = spval
      a_trefdt  (:,:) = spval
      a_tafudt  (:,:) = spval

      a_fsenant (:,:) = spval
      a_lfevpant(:,:) = spval
      a_fgrndnt (:,:) = spval
      a_olrgnt  (:,:) = spval
      a_rnetnt  (:,:) = spval
      a_t_grndnt(:,:) = spval
      a_tradnt  (:,:) = spval
      a_trefnt  (:,:) = spval
      a_tafunt  (:,:) = spval

      a_sr     (:,:) = spval
      a_solvd  (:,:) = spval
      a_solvi  (:,:) = spval
      a_solnd  (:,:) = spval
      a_solni  (:,:) = spval
      a_srvd   (:,:) = spval
      a_srvi   (:,:) = spval
      a_srnd   (:,:) = spval
      a_srni   (:,:) = spval
      a_solvdln(:,:) = spval
      a_solviln(:,:) = spval
      a_solndln(:,:) = spval
      a_solniln(:,:) = spval
      a_srvdln (:,:) = spval
      a_srviln (:,:) = spval
      a_srndln (:,:) = spval
      a_srniln (:,:) = spval

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

#ifdef URBAN_MODEL
               u = patch2urb(np)
               IF (u > 0) THEN
                  urbwt(i,j) = urbwt(i,j) + patchfrac(np)
                  a_t_room (i,j) = a_t_room (i,j) + patchfrac(np)*t_room (u)
                  a_tafu   (i,j) = a_tafu   (i,j) + patchfrac(np)*tafu   (u)
                  a_fhac   (i,j) = a_fhac   (i,j) + patchfrac(np)*fhac   (u)
                  a_fwst   (i,j) = a_fwst   (i,j) + patchfrac(np)*fwst   (u)
                  a_fach   (i,j) = a_fach   (i,j) + patchfrac(np)*fach   (u)
               ENDIF
#endif
               !TODO: 根据coszen(np)的正负->daytime or nighttime
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

#ifdef URBAN_MODEL
               IF(urbwt(i,j).gt.0.00001)THEN
                  a_t_room (i,j) = a_t_room (i,j) / urbwt(i,j)
                  a_tafu   (i,j) = a_tafu   (i,j) / urbwt(i,j)
                  a_fhac   (i,j) = a_fhac   (i,j) / urbwt(i,j)
                  a_fwst   (i,j) = a_fwst   (i,j) / urbwt(i,j)
                  a_fach   (i,j) = a_fach   (i,j) / urbwt(i,j)
               ELSE
                  a_t_room (i,j) = spval
                  a_tafu   (i,j) = spval
                  a_fhac   (i,j) = spval
                  a_fwst   (i,j) = spval
                  a_fach   (i,j) = spval
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

            CALL acc(a_taux   (i,j), 1., f_taux   (i,j))
            CALL acc(a_tauy   (i,j), 1., f_tauy   (i,j))
            CALL acc(a_fsena  (i,j), 1., f_fsena  (i,j))
            CALL acc(a_lfevpa (i,j), 1., f_lfevpa (i,j))
            CALL acc(a_fevpa  (i,j), 1., f_fevpa  (i,j))
            CALL acc(a_fsenl  (i,j), 1., f_fsenl  (i,j))
            CALL acc(a_fevpl  (i,j), 1., f_fevpl  (i,j))
            CALL acc(a_etr    (i,j), 1., f_etr    (i,j))
            CALL acc(a_fseng  (i,j), 1., f_fseng  (i,j))
            CALL acc(a_fevpg  (i,j), 1., f_fevpg  (i,j))
            CALL acc(a_fgrnd  (i,j), 1., f_fgrnd  (i,j))
            CALL acc(a_sabvsun(i,j), 1., f_sabvsun(i,j))
            CALL acc(a_sabvsha(i,j), 1., f_sabvsha(i,j))
            CALL acc(a_sabg   (i,j), 1., f_sabg   (i,j))
            CALL acc(a_olrg   (i,j), 1., f_olrg   (i,j))
            CALL acc(a_rnet   (i,j), 1., f_rnet   (i,j))
            CALL acc(a_xerr   (i,j), 1., f_xerr   (i,j))
            CALL acc(a_zerr   (i,j), 1., f_zerr   (i,j))
            CALL acc(a_rsur   (i,j), 1., f_rsur   (i,j))
            CALL acc(a_rnof   (i,j), 1., f_rnof   (i,j))
            CALL acc(a_qintr  (i,j), 1., f_qintr  (i,j))
            CALL acc(a_qinfl  (i,j), 1., f_qinfl  (i,j))
            CALL acc(a_qdrip  (i,j), 1., f_qdrip  (i,j))
            CALL acc(a_rstfac (i,j), 1., f_rstfac (i,j))
            CALL acc(a_zwt    (i,j), 1., f_zwt    (i,j))
            CALL acc(a_wa     (i,j), 1., f_wa     (i,j))
            CALL acc(a_wat    (i,j), 1., f_wat    (i,j))
            CALL acc(a_assim  (i,j), 1., f_assim  (i,j))
            CALL acc(a_respc  (i,j), 1., f_respc  (i,j))

            CALL acc(a_qcharge(i,j), 1., f_qcharge(i,j))

            CALL acc(a_t_grnd (i,j), 1., f_t_grnd (i,j))
            CALL acc(a_tleaf  (i,j), 1., f_tleaf  (i,j))
            CALL acc(a_ldew   (i,j), 1., f_ldew   (i,j))
            CALL acc(a_scv    (i,j), 1., f_scv    (i,j))
            CALL acc(a_snowdp (i,j), 1., f_snowdp (i,j))
            CALL acc(a_fsno   (i,j), 1., f_fsno   (i,j))
            CALL acc(a_sigf   (i,j), 1., f_sigf   (i,j))
            CALL acc(a_green  (i,j), 1., f_green  (i,j))
            CALL acc(a_lai    (i,j), 1., f_lai    (i,j))
            CALL acc(a_laisun (i,j), 1., f_laisun (i,j))
            CALL acc(a_laisha (i,j), 1., f_laisha (i,j))
            CALL acc(a_sai    (i,j), 1., f_sai    (i,j))
            CALL acc(a_alb(1,1,i,j), 1., f_alb(1,1,i,j))
            CALL acc(a_alb(2,1,i,j), 1., f_alb(2,1,i,j))
            CALL acc(a_alb(1,2,i,j), 1., f_alb(1,2,i,j))
            CALL acc(a_alb(2,2,i,j), 1., f_alb(2,2,i,j))
            CALL acc(a_emis   (i,j), 1., f_emis   (i,j))
            CALL acc(a_z0m    (i,j), 1., f_z0m    (i,j))
            CALL acc(a_trad   (i,j), 1., f_trad   (i,j))
            CALL acc(a_tref   (i,j), 1., f_tref   (i,j))
            CALL acc(a_qref   (i,j), 1., f_qref   (i,j))
            CALL acc(a_xy_rain(i,j), 1., f_xy_rain(i,j))
            CALL acc(a_xy_snow(i,j), 1., f_xy_snow(i,j))

            CALL acc(a_t_room (i,j), 1., f_t_room (i,j))
            CALL acc(a_tafu   (i,j), 1., f_tafu   (i,j))
            CALL acc(a_fhac   (i,j), 1., f_fhac   (i,j))
            CALL acc(a_fwst   (i,j), 1., f_fwst   (i,j))
            CALL acc(a_fach   (i,j), 1., f_fach   (i,j))

            CALL acc(a_sabvdt  (i,j), 1., f_sabvdt  (i,j))
            CALL acc(a_sabgdt  (i,j), 1., f_sabgdt  (i,j))
            CALL acc(a_srdt    (i,j), 1., f_srdt    (i,j))
            CALL acc(a_fsenadt (i,j), 1., f_fsenadt (i,j))
            CALL acc(a_lfevpadt(i,j), 1., f_lfevpadt(i,j))
            CALL acc(a_fgrnddt (i,j), 1., f_fgrnddt (i,j))
            CALL acc(a_olrgdt  (i,j), 1., f_olrgdt  (i,j))
            CALL acc(a_rnetdt  (i,j), 1., f_rnetdt  (i,j))
            CALL acc(a_t_grnddt(i,j), 1., f_t_grnddt(i,j))
            CALL acc(a_traddt  (i,j), 1., f_traddt  (i,j))
            CALL acc(a_trefdt  (i,j), 1., f_trefdt  (i,j))
            CALL acc(a_tafudt  (i,j), 1., f_tafudt  (i,j))

            CALL acc(a_fsenant (i,j), 1., f_fsenant (i,j))
            CALL acc(a_lfevpant(i,j), 1., f_lfevpant(i,j))
            CALL acc(a_fgrndnt (i,j), 1., f_fgrndnt (i,j))
            CALL acc(a_olrgnt  (i,j), 1., f_olrgnt  (i,j))
            CALL acc(a_rnetnt  (i,j), 1., f_rnetnt  (i,j))
            CALL acc(a_t_grndnt(i,j), 1., f_t_grndnt(i,j))
            CALL acc(a_tradnt  (i,j), 1., f_tradnt  (i,j))
            CALL acc(a_trefnt  (i,j), 1., f_trefnt  (i,j))
            CALL acc(a_tafunt  (i,j), 1., f_tafunt  (i,j))

            DO l = maxsnl+1, nl_soil
               CALL acc(a_t_soisno   (l,i,j), 1., f_t_soisno   (l,i,j))
               CALL acc(a_wliq_soisno(l,i,j), 1., f_wliq_soisno(l,i,j))
               CALL acc(a_wice_soisno(l,i,j), 1., f_wice_soisno(l,i,j))
            ENDDO

            DO l = 1, nl_soil
               CALL acc(a_h2osoi     (l,i,j), 1., f_h2osoi     (l,i,j))
            ENDDO

            DO l = 1, nl_lake
               CALL acc(a_t_lake(l,i,j), 1., f_t_lake(l,i,j))
               CALL acc(a_lake_icefrac(l,i,j), 1., f_lake_icefrac(l,i,j))
            ENDDO

            CALL acc(a_ustar(i,j), 1., f_ustar(i,j))
            CALL acc(a_tstar(i,j), 1., f_tstar(i,j))
            CALL acc(a_qstar(i,j), 1., f_qstar(i,j))
            CALL acc(a_zol  (i,j), 1., f_zol  (i,j))
            CALL acc(a_rib  (i,j), 1., f_rib  (i,j))
            CALL acc(a_fm   (i,j), 1., f_fm   (i,j))
            CALL acc(a_fh   (i,j), 1., f_fh   (i,j))
            CALL acc(a_fq   (i,j), 1., f_fq   (i,j))

            CALL acc(a_us10m(i,j), 1., f_us10m(i,j))
            CALL acc(a_vs10m(i,j), 1., f_vs10m(i,j))
            CALL acc(a_fm10m(i,j), 1., f_fm10m(i,j))

            CALL acc(a_sr     (i,j), 1., f_sr     (i,j))
            CALL acc(a_solvd  (i,j), 1., f_solvd  (i,j))
            CALL acc(a_solvi  (i,j), 1., f_solvi  (i,j))
            CALL acc(a_solnd  (i,j), 1., f_solnd  (i,j))
            CALL acc(a_solni  (i,j), 1., f_solni  (i,j))
            CALL acc(a_srvd   (i,j), 1., f_srvd   (i,j))
            CALL acc(a_srvi   (i,j), 1., f_srvi   (i,j))
            CALL acc(a_srnd   (i,j), 1., f_srnd   (i,j))
            CALL acc(a_srni   (i,j), 1., f_srni   (i,j))
            CALL acc(a_solvdln(i,j), 1., f_solvdln(i,j))
            CALL acc(a_solviln(i,j), 1., f_solviln(i,j))
            CALL acc(a_solndln(i,j), 1., f_solndln(i,j))
            CALL acc(a_solniln(i,j), 1., f_solniln(i,j))
            CALL acc(a_srvdln (i,j), 1., f_srvdln (i,j))
            CALL acc(a_srviln (i,j), 1., f_srviln (i,j))
            CALL acc(a_srndln (i,j), 1., f_srndln (i,j))
            CALL acc(a_srniln (i,j), 1., f_srniln (i,j))

            IF (a_solvdln(i,j) /= spval) nac_ln(i,j) = nac_ln(i,j) + 1
            IF (a_trefdt (i,j) /= spval) nac_dt(i,j) = nac_dt(i,j) + 1
            IF (a_trefnt (i,j) /= spval) nac_nt(i,j) = nac_nt(i,j) + 1

         ENDDO
      ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

END SUBROUTINE vec2xy

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

! ----- EOP ---------
