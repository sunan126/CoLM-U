#include <define.h>

 SUBROUTINE flxwrite (idate,nac,nac_24,nac_ln,nac_dt,nac_nt,dir_output,casename)

!=======================================================================
! Original version: Yongjiu Dai, September 15, 1999, 03/2014
!=======================================================================

  USE precision
  USE GlobalVars
  USE MOD_2D_Fluxes
  USE MOD_TimeInvariants, only: gridlond, gridlatd
  USE timemanager
  USE netcdf
  USE ncio

  IMPLICIT NONE

  INTEGER, intent(in) :: idate(3)
  INTEGER, intent(in) :: nac
  INTEGER, intent(inout) :: nac_24
  INTEGER, intent(in) :: nac_ln(lon_points,lat_points)
  INTEGER, intent(in) :: nac_dt(lon_points,lat_points)
  INTEGER, intent(in) :: nac_nt(lon_points,lat_points)

  CHARACTER(LEN=256) :: dir_output
  CHARACTER(LEN=256) :: casename

  INTEGER luout, month, day, i, j, l
  CHARACTER(LEN=256) fout
  CHARACTER(LEN=256) cdate
  REAL(r8) a

  INTEGER  :: ix,iy,ilev
  INTEGER  :: ncid, xid, yid, varid, sslevid, lakelevid, bandid
  INTEGER  :: lake_lev(nl_lake), sslev(nl_soil-maxsnl), band(2)
  REAL(r4) :: lons(lon_points), lats(lat_points)
  REAL(r8) :: tmp(lon_points, lat_points, 2)
  REAL(r8) :: tmp1(lon_points, lat_points, maxsnl+1:nl_soil)
  REAL(r8) :: tmp2(lon_points, lat_points, nl_lake)

! ----------------------------------------------------------------------
! Open for model time varying data (model state variables) and history filed

     luout = 100
#if(defined WO_MONTHLY)
     CALL julian2monthday(idate(1), idate(2), month, day)
     write(cdate,'(i4.4,"-",i2.2)') idate(1), month
#else
     write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1),idate(2),idate(3)
#endif
     fout = trim(dir_output)//trim(casename)//'_'//'2D_Fluxes'//'_'//trim(cdate)
     print*,trim(fout)
     open(unit=luout,file=fout,access='sequential',form='unformatted',&
                     status='unknown',action='write')

     a = float(nac)

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) PRIVATE(i,j,l)
#endif
     DO j = 1, lat_points
        DO i = 1, lon_points

           IF (f_taux   (i,j) /= spval) f_taux   (i,j) = f_taux   (i,j) / a  ! wind stress: E-W [kg/m/s2]
           IF (f_tauy   (i,j) /= spval) f_tauy   (i,j) = f_tauy   (i,j) / a  ! wind stress: N-S [kg/m/s2]
           IF (f_fsena  (i,j) /= spval) f_fsena  (i,j) = f_fsena  (i,j) / a  ! sensible heat from canopy height to atmosphere [W/m2]
           IF (f_lfevpa (i,j) /= spval) f_lfevpa (i,j) = f_lfevpa (i,j) / a  ! latent heat flux from canopy height to atmosphere [W/m2]
           IF (f_fevpa  (i,j) /= spval) f_fevpa  (i,j) = f_fevpa  (i,j) / a  ! evapotranspiration from canopy to atmosphere [mm/s]
           IF (f_fsenl  (i,j) /= spval) f_fsenl  (i,j) = f_fsenl  (i,j) / a  ! sensible heat from leaves [W/m2]
           IF (f_fevpl  (i,j) /= spval) f_fevpl  (i,j) = f_fevpl  (i,j) / a  ! evaporation+transpiration from leaves [mm/s]
           IF (f_etr    (i,j) /= spval) f_etr    (i,j) = f_etr    (i,j) / a  ! transpiration rate [mm/s]
           IF (f_fseng  (i,j) /= spval) f_fseng  (i,j) = f_fseng  (i,j) / a  ! sensible heat flux from ground [W/m2]
           IF (f_fevpg  (i,j) /= spval) f_fevpg  (i,j) = f_fevpg  (i,j) / a  ! evaporation heat flux from ground [mm/s]
           IF (f_fgrnd  (i,j) /= spval) f_fgrnd  (i,j) = f_fgrnd  (i,j) / a  ! ground heat flux [W/m2]
           IF (f_sabvsun(i,j) /= spval) f_sabvsun(i,j) = f_sabvsun(i,j) / a  ! solar absorbed by sunlit canopy [W/m2]
           IF (f_sabvsha(i,j) /= spval) f_sabvsha(i,j) = f_sabvsha(i,j) / a  ! solar absorbed by shaded [W/m2]
           IF (f_sabg   (i,j) /= spval) f_sabg   (i,j) = f_sabg   (i,j) / a  ! solar absorbed by ground  [W/m2]
           IF (f_olrg   (i,j) /= spval) f_olrg   (i,j) = f_olrg   (i,j) / a  ! outgoing long-wave radiation from ground+canopy [W/m2]
           IF (f_rnet   (i,j) /= spval) f_rnet   (i,j) = f_rnet   (i,j) / a  ! net radiation [W/m2]
           IF (f_xerr   (i,j) /= spval) f_xerr   (i,j) = f_xerr   (i,j) / a  ! the error of water banace [mm/s]
           IF (f_zerr   (i,j) /= spval) f_zerr   (i,j) = f_zerr   (i,j) / a  ! the error of energy balance [W/m2]
           IF (f_rsur   (i,j) /= spval) f_rsur   (i,j) = f_rsur   (i,j) / a  ! surface runoff [mm/s]
           IF (f_rnof   (i,j) /= spval) f_rnof   (i,j) = f_rnof   (i,j) / a  ! total runoff [mm/s]
           IF (f_qintr  (i,j) /= spval) f_qintr  (i,j) = f_qintr  (i,j) / a  ! interception [mm/s]
           IF (f_qinfl  (i,j) /= spval) f_qinfl  (i,j) = f_qinfl  (i,j) / a  ! inflitraton [mm/s]
           IF (f_qdrip  (i,j) /= spval) f_qdrip  (i,j) = f_qdrip  (i,j) / a  ! throughfall [mm/s]
           IF (f_rstfac (i,j) /= spval) f_rstfac (i,j) = f_rstfac (i,j) / a  ! factor of soil water stress
           IF (f_zwt    (i,j) /= spval) f_zwt    (i,j) = f_zwt    (i,j) / a  ! water depth [m]
           IF (f_wa     (i,j) /= spval) f_wa     (i,j) = f_wa     (i,j) / a  ! water storage in aquifer [mm]
           IF (f_wat    (i,j) /= spval) f_wat    (i,j) = f_wat    (i,j) / a  ! total water storage [mm]
           IF (f_assim  (i,j) /= spval) f_assim  (i,j) = f_assim  (i,j) / a  ! canopy assimilation rate [mol m-2 s-1]
           IF (f_respc  (i,j) /= spval) f_respc  (i,j) = f_respc  (i,j) / a  ! respiration (plant+soil) [mol m-2 s-1]
           IF (f_qcharge(i,j) /= spval) f_qcharge(i,j) = f_qcharge(i,j) / a  ! groundwater recharge rate [mm/s]

!---------------------------------------------------------------------
           IF (nac_24 == 0) nac_24 = 1
           IF (f_t_grnd (i,j) /= spval) f_t_grnd (i,j) = f_t_grnd (i,j) / a  ! ground surface temperature [K]
           IF (f_tleaf  (i,j) /= spval) f_tleaf  (i,j) = f_tleaf  (i,j) / a  ! sunlit leaf temperature [K]
           IF (f_ldew   (i,j) /= spval) f_ldew   (i,j) = f_ldew   (i,j) / a  ! depth of water on foliage [mm]
           IF (f_scv    (i,j) /= spval) f_scv    (i,j) = f_scv    (i,j) / a  ! snow cover, water equivalent [mm]
           IF (f_snowdp (i,j) /= spval) f_snowdp (i,j) = f_snowdp (i,j) / a  ! snow depth [meter]
           IF (f_fsno   (i,j) /= spval) f_fsno   (i,j) = f_fsno   (i,j) / a  ! fraction of snow cover on ground
           IF (f_sigf   (i,j) /= spval) f_sigf   (i,j) = f_sigf   (i,j) / a  ! fraction of veg cover, excluding snow-covered veg [-]
           IF (f_green  (i,j) /= spval) f_green  (i,j) = f_green  (i,j) / a  ! leaf greenness
           IF (f_lai    (i,j) /= spval) f_lai    (i,j) = f_lai    (i,j) / a  ! leaf area index
           IF (f_laisun (i,j) /= spval) f_laisun (i,j) = f_laisun (i,j) / a  ! leaf area index
           IF (f_laisha (i,j) /= spval) f_laisha (i,j) = f_laisha (i,j) / a  ! leaf area index
           IF (f_sai    (i,j) /= spval) f_sai    (i,j) = f_sai    (i,j) / a  ! stem area index
           IF (f_alb(1,1,i,j) /= spval) f_alb(1,1,i,j) = f_alb(1,1,i,j) / a  ! averaged albedo [visible, direct; direct, diffuse]
           IF (f_alb(2,1,i,j) /= spval) f_alb(2,1,i,j) = f_alb(2,1,i,j) / a  ! averaged albedo [visible, direct; direct, diffuse]
           IF (f_alb(1,2,i,j) /= spval) f_alb(1,2,i,j) = f_alb(1,2,i,j) / a  ! averaged albedo [visible, direct; direct, diffuse]
           IF (f_alb(2,2,i,j) /= spval) f_alb(2,2,i,j) = f_alb(2,2,i,j) / a  ! averaged albedo [visible, direct; direct, diffuse]
           IF (f_emis   (i,j) /= spval) f_emis   (i,j) = f_emis   (i,j) / a  ! averaged bulk surface emissivity
           IF (f_z0m    (i,j) /= spval) f_z0m    (i,j) = f_z0m    (i,j) / a  ! effective roughness [m]
           IF (f_trad   (i,j) /= spval) f_trad   (i,j) = f_trad   (i,j) / a  ! radiative temperature of surface [K]
           IF (f_tref   (i,j) /= spval) f_tref   (i,j) = f_tref   (i,j) / a  ! 2 m height air temperature [kelvin]
           IF (f_tmax   (i,j) /= spval) f_tmax   (i,j) = f_tmax   (i,j) / nac_24  ! Diurnal Max 2 m height air temperature [kelvin]
           IF (f_tmin   (i,j) /= spval) f_tmin   (i,j) = f_tmin   (i,j) / nac_24  ! Diurnal Min 2 m height air temperature [kelvin]
           IF (f_tdtr   (i,j) /= spval) f_tdtr   (i,j) = f_tdtr   (i,j) / nac_24  ! DTR of 2 m height air temperature [kelvin]
           IF (f_qref   (i,j) /= spval) f_qref   (i,j) = f_qref   (i,j) / a  ! 2 m height air specific humidity [kg/kg]
           IF (f_t_room (i,j) /= spval) f_t_room (i,j) = f_t_room (i,j) / a  ! temperature of inner building [K]
           IF (f_tafu   (i,j) /= spval) f_tafu   (i,j) = f_tafu   (i,j) / a  ! temperature of outer building [K]
           IF (f_fhac   (i,j) /= spval) f_fhac   (i,j) = f_fhac   (i,j) / a  ! sensible flux from heat or cool AC [W/m2]
           IF (f_fwst   (i,j) /= spval) f_fwst   (i,j) = f_fwst   (i,j) / a  ! waste heat flux from heat or cool AC [W/m2]
           IF (f_fach   (i,j) /= spval) f_fach   (i,j) = f_fach   (i,j) / a  ! flux from inner and outter air exchange [W/m2]
           IF (f_fahe   (i,j) /= spval) f_fahe   (i,j) = f_fahe   (i,j) / a  ! flux from metabolic and vehicle [W/m2]
           IF (f_fhah   (i,j) /= spval) f_fhah   (i,j) = f_fhah   (i,j) / a  ! flux from heatiog [W/m2]  
           IF (f_vehc   (i,j) /= spval) f_vehc   (i,j) = f_vehc   (i,j) / a  ! flux from vehicle [W/m2]
           IF (f_meta   (i,j) /= spval) f_meta   (i,j) = f_meta   (i,j) / a  ! flux from metabolic [W/m2]
           IF (f_xy_rain(i,j) /= spval) f_xy_rain(i,j) = f_xy_rain(i,j) / a  ! rain [mm/s]
           IF (f_xy_snow(i,j) /= spval) f_xy_snow(i,j) = f_xy_snow(i,j) / a  ! snow [mm/s]

           IF (f_senroof(i,j) /= spval) f_senroof   (i,j) = f_senroof   (i,j) / a  ! sensible heat flux from roof [W/m2]
           IF (f_senwsun(i,j) /= spval) f_senwsun   (i,j) = f_senwsun   (i,j) / a  ! sensible heat flux from sunlit wall [W/m2]
           IF (f_senwsha(i,j) /= spval) f_senwsha   (i,j) = f_senwsha   (i,j) / a  ! sensible heat flux from shaded wall [W/m2]
           IF (f_sengimp(i,j) /= spval) f_sengimp   (i,j) = f_sengimp   (i,j) / a  ! sensible heat flux from impervious road [W/m2]
           IF (f_sengper(i,j) /= spval) f_sengper   (i,j) = f_sengper   (i,j) / a  ! sensible heat flux from pervious road [W/m2]
           IF (f_senurbl(i,j) /= spval) f_senurbl   (i,j) = f_senurbl   (i,j) / a  ! sensible heat flux from urban vegetation [W/m2]

           IF (f_lfevproof(i,j) /= spval) f_lfevproof   (i,j) = f_lfevproof   (i,j) / a  ! latent heat flux from roof [W/m2]
           IF (f_lfevpgimp(i,j) /= spval) f_lfevpgimp   (i,j) = f_lfevpgimp   (i,j) / a  ! latent heat flux from impervious road [W/m2]
           IF (f_lfevpgper(i,j) /= spval) f_lfevpgper   (i,j) = f_lfevpgper   (i,j) / a  ! latent heat flux from pervious road [W/m2]
           IF (f_lfevpurbl(i,j) /= spval) f_lfevpurbl   (i,j) = f_lfevpurbl   (i,j) / a  ! latent heat flux from urban vegetation [W/m2]

           IF (f_troof   (i,j) /= spval) f_troof   (i,j) = f_troof   (i,j) / a  ! temperature of roof [K]
           IF (f_twall   (i,j) /= spval) f_twall   (i,j) = f_twall   (i,j) / a  ! temperature of wall [K]

           IF (f_sabvdt  (i,j) /= spval) f_sabvdt  (i,j) = f_sabvdt  (i,j) / nac_dt(i,j)  ! solar absorbed by sunlit canopy [w/m2]
           IF (f_sabgdt  (i,j) /= spval) f_sabgdt  (i,j) = f_sabgdt  (i,j) / nac_dt(i,j)  ! solar absorbed by ground [w/m2]
           IF (f_srdt    (i,j) /= spval) f_srdt    (i,j) = f_srdt    (i,j) / nac_dt(i,j)  ! total reflected solar radiation (w/m2)
           IF (f_fsenadt (i,j) /= spval) f_fsenadt (i,j) = f_fsenadt (i,j) / nac_dt(i,j)  ! sensible heat from canopy height to atmosphere [w/m2]
           IF (f_lfevpadt(i,j) /= spval) f_lfevpadt(i,j) = f_lfevpadt(i,j) / nac_dt(i,j)  ! latent heat flux from canopy height to atmosphere [w/m2]
           IF (f_fgrnddt (i,j) /= spval) f_fgrnddt (i,j) = f_fgrnddt (i,j) / nac_dt(i,j)  ! ground heat flux [w/m2]
           IF (f_olrgdt  (i,j) /= spval) f_olrgdt  (i,j) = f_olrgdt  (i,j) / nac_dt(i,j)  ! outgoing long-wave radiation from ground+canopy [w/m2]
           IF (f_rnetdt  (i,j) /= spval) f_rnetdt  (i,j) = f_rnetdt  (i,j) / nac_dt(i,j)  ! net radiation [w/m2]
           IF (f_t_grnddt(i,j) /= spval) f_t_grnddt(i,j) = f_t_grnddt(i,j) / nac_dt(i,j)  ! ground surface temperature [k]
           IF (f_traddt  (i,j) /= spval) f_traddt  (i,j) = f_traddt  (i,j) / nac_dt(i,j)  ! radiative temperature of surface [k]
           IF (f_trefdt  (i,j) /= spval) f_trefdt  (i,j) = f_trefdt  (i,j) / nac_dt(i,j)  ! 2 m height air temperature [kelvin]
           IF (f_tafudt  (i,j) /= spval) f_tafudt  (i,j) = f_tafudt  (i,j) / nac_dt(i,j)  ! temperature of outer building [K]

           IF (f_fsenant (i,j) /= spval) f_fsenant (i,j) = f_fsenant (i,j) / nac_nt(i,j)  ! sensible heat from canopy height to atmosphere [w/m2]
           IF (f_lfevpant(i,j) /= spval) f_lfevpant(i,j) = f_lfevpant(i,j) / nac_nt(i,j)  ! latent heat flux from canopy height to atmosphere [w/m2]
           IF (f_fgrndnt (i,j) /= spval) f_fgrndnt (i,j) = f_fgrndnt (i,j) / nac_nt(i,j)  ! ground heat flux [w/m2]
           IF (f_olrgnt  (i,j) /= spval) f_olrgnt  (i,j) = f_olrgnt  (i,j) / nac_nt(i,j)  ! outgoing long-wave radiation from ground+canopy [w/m2]
           IF (f_rnetnt  (i,j) /= spval) f_rnetnt  (i,j) = f_rnetnt  (i,j) / nac_nt(i,j)  ! net radiation [w/m2]
           IF (f_t_grndnt(i,j) /= spval) f_t_grndnt(i,j) = f_t_grndnt(i,j) / nac_nt(i,j)  ! ground surface temperature [k]
           IF (f_tradnt  (i,j) /= spval) f_tradnt  (i,j) = f_tradnt  (i,j) / nac_nt(i,j)  ! radiative temperature of surface [k]
           IF (f_trefnt  (i,j) /= spval) f_trefnt  (i,j) = f_trefnt  (i,j) / nac_nt(i,j)  ! 2 m height air temperature [kelvin]
           IF (f_tafunt  (i,j) /= spval) f_tafunt  (i,j) = f_tafunt  (i,j) / nac_nt(i,j)  ! temperature of outer building [K]

!---------------------------------------------------------------------
           DO l = maxsnl+1, nl_soil
              IF (f_t_soisno   (l,i,j) /= spval) f_t_soisno   (l,i,j) = f_t_soisno   (l,i,j) / a  ! soil temperature [K]
              IF (f_wliq_soisno(l,i,j) /= spval) f_wliq_soisno(l,i,j) = f_wliq_soisno(l,i,j) / a  ! liquid water in soil layers [kg/m2]
              IF (f_wice_soisno(l,i,j) /= spval) f_wice_soisno(l,i,j) = f_wice_soisno(l,i,j) / a  ! ice lens in soil layers [kg/m2]
           ENDDO

           DO l = 1, nl_soil
              IF (f_h2osoi     (l,i,j) /= spval) f_h2osoi     (l,i,j) = f_h2osoi     (l,i,j) / a  ! volumetric soil water in layers [m3/m3]
           ENDDO

           DO l = 1, nl_lake
              IF (f_t_lake      (l,i,j) /= spval) f_t_lake      (l,i,j) = f_t_lake      (l,i,j) / a  ! lake temperature [K]
              IF (f_lake_icefrac(l,i,j) /= spval) f_lake_icefrac(l,i,j) = f_lake_icefrac(l,i,j) / a  ! lake ice fraction cover [0-1]
           ENDDO

           IF (f_ustar  (i,j) /= spval) f_ustar  (i,j) = f_ustar  (i,j) / a  ! u* in similarity theory [m/s]
           IF (f_tstar  (i,j) /= spval) f_tstar  (i,j) = f_tstar  (i,j) / a  ! t* in similarity theory [kg/kg]
           IF (f_qstar  (i,j) /= spval) f_qstar  (i,j) = f_qstar  (i,j) / a  ! q* in similarity theory [kg/kg]
           IF (f_zol    (i,j) /= spval) f_zol    (i,j) = f_zol    (i,j) / a  ! dimensionless height (z/L) used in Monin-Obukhov theory
           IF (f_rib    (i,j) /= spval) f_rib    (i,j) = f_rib    (i,j) / a  ! bulk Richardson number in surface layer
           IF (f_fm     (i,j) /= spval) f_fm     (i,j) = f_fm     (i,j) / a  ! integral of profile function for momentum
           IF (f_fh     (i,j) /= spval) f_fh     (i,j) = f_fh     (i,j) / a  ! integral of profile function for heat
           IF (f_fq     (i,j) /= spval) f_fq     (i,j) = f_fq     (i,j) / a  ! integral of profile function for moisture
           IF (f_us10m  (i,j) /= spval) f_us10m  (i,j) = f_us10m  (i,j) / a  ! 10m u-velocity [m/s]
           IF (f_vs10m  (i,j) /= spval) f_vs10m  (i,j) = f_vs10m  (i,j) / a  ! 10m v-velocity [m/s]
           IF (f_fm10m  (i,j) /= spval) f_fm10m  (i,j) = f_fm10m  (i,j) / a  ! integral of profile function for momentum at 10m [-]

           IF (f_xy_us  (i,j) /= spval) f_xy_us  (i,j) = f_xy_us  (i,j) / a  ! wind in eastward direction [m/s]
           IF (f_xy_vs  (i,j) /= spval) f_xy_vs  (i,j) = f_xy_vs  (i,j) / a  ! wind in northward direction [m/s]
           IF (f_xy_t   (i,j) /= spval) f_xy_t   (i,j) = f_xy_t   (i,j) / a  ! temperature at reference height [kelvin]
           IF (f_xy_q   (i,j) /= spval) f_xy_q   (i,j) = f_xy_q   (i,j) / a  ! specific humidity at reference height [kg/kg]
           IF (f_xy_prc (i,j) /= spval) f_xy_prc (i,j) = f_xy_prc (i,j) / a  ! convective precipitation [mm/s]
           IF (f_xy_prl (i,j) /= spval) f_xy_prl (i,j) = f_xy_prl (i,j) / a  ! large scale precipitation [mm/s]
           IF (f_xy_pbot(i,j) /= spval) f_xy_pbot(i,j) = f_xy_pbot(i,j) / a  ! atmospheric pressure at the surface [pa]
           IF (f_xy_frl (i,j) /= spval) f_xy_frl (i,j) = f_xy_frl (i,j) / a  ! atmospheric infrared (longwave) radiation [W/m2]
           IF (f_xy_solarin(i,j) /= spval) f_xy_solarin(i,j) = f_xy_solarin(i,j) / a  ! downward solar radiation at surface [W/m2]

           IF (f_sr     (i,j) /= spval) f_sr     (i,j) = f_sr     (i,j) / a  ! total reflected solar radiation (W/m2)
           IF (f_solvd  (i,j) /= spval) f_solvd  (i,j) = f_solvd  (i,j) / a  ! incident direct beam vis solar radiation (W/m2)
           IF (f_solvi  (i,j) /= spval) f_solvi  (i,j) = f_solvi  (i,j) / a  ! incident diffuse beam vis solar radiation (W/m2)
           IF (f_solnd  (i,j) /= spval) f_solnd  (i,j) = f_solnd  (i,j) / a  ! incident direct beam nir solar radiation (W/m2)
           IF (f_solni  (i,j) /= spval) f_solni  (i,j) = f_solni  (i,j) / a  ! incident diffuse beam nir solar radiation (W/m2)
           IF (f_srvd   (i,j) /= spval) f_srvd   (i,j) = f_srvd   (i,j) / a  ! reflected direct beam vis solar radiation (W/m2)
           IF (f_srvi   (i,j) /= spval) f_srvi   (i,j) = f_srvi   (i,j) / a  ! reflected diffuse beam vis solar radiation (W/m2)
           IF (f_srnd   (i,j) /= spval) f_srnd   (i,j) = f_srnd   (i,j) / a  ! reflected direct beam nir solar radiation (W/m2)
           IF (f_srni   (i,j) /= spval) f_srni   (i,j) = f_srni   (i,j) / a  ! reflected diffuse beam nir solar radiation (W/m2)
           IF (f_solvdln(i,j) /= spval) f_solvdln(i,j) = f_solvdln(i,j) / nac_ln(i,j) ! incident direct beam vis solar radiation at local noon (W/m2)
           IF (f_solviln(i,j) /= spval) f_solviln(i,j) = f_solviln(i,j) / nac_ln(i,j) ! incident diffuse beam vis solar radiation at local noon (W/m2)
           IF (f_solndln(i,j) /= spval) f_solndln(i,j) = f_solndln(i,j) / nac_ln(i,j) ! incident direct beam nir solar radiation at local noon (W/m2)
           IF (f_solniln(i,j) /= spval) f_solniln(i,j) = f_solniln(i,j) / nac_ln(i,j) ! incident diffuse beam nir solar radiation at local noon (W/m2)
           IF (f_srvdln (i,j) /= spval) f_srvdln (i,j) = f_srvdln (i,j) / nac_ln(i,j) ! reflected direct beam vis solar radiation at local noon (W/m2)
           IF (f_srviln (i,j) /= spval) f_srviln (i,j) = f_srviln (i,j) / nac_ln(i,j) ! reflected diffuse beam vis solar radiation at local noon (W/m2)
           IF (f_srndln (i,j) /= spval) f_srndln (i,j) = f_srndln (i,j) / nac_ln(i,j) ! reflected direct beam nir solar radiation at local noon (W/m2)
           IF (f_srniln (i,j) /= spval) f_srniln (i,j) = f_srniln (i,j) / nac_ln(i,j) ! reflected diffuse beam nir solar radiation at local noon(W/m2)

        ENDDO
     ENDDO
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

! 单点、时间步长输出，减少后处理时间，增加部分变量，便于Urban-PLUMBER2比较
#ifdef NC_OUTPUT
    ! create netcdf and define dimensions
      CALL nccheck( nf90_create(trim(fout)//'.nc', nf90_64bit_offset, ncid) )

    ! dimensions
      CALL nccheck( nf90_def_dim(ncid, 'lon',           lon_points,     xid) )
      CALL nccheck( nf90_def_dim(ncid, 'lat',           lat_points,     yid) )
      CALL nccheck( nf90_def_dim(ncid, 'soil_snow_lev', nl_soil-maxsnl, sslevid) )
      CALL nccheck( nf90_def_dim(ncid, 'lake_lev',      nl_lake,        lakelevid) )
      CALL nccheck( nf90_def_dim(ncid, 'band',          2,              bandid) )

    ! global attr
      CALL nccheck( nf90_put_att(ncid, NF90_GLOBAL, 'title','CLM SIMULATED SURFACE FLUXES') )

    ! variables
    ! dimension variables
      CALL nccheck( nf90_def_var(ncid, 'lon', nf90_float, (/xid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','longitude') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','degrees_east') )

      CALL nccheck( nf90_def_var(ncid, 'lat', nf90_float, (/yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','latitude') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','degrees_north') )

      CALL nccheck( nf90_def_var(ncid, 'soil_snow_lev', nf90_int, (/sslevid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name',"soil plus snow level") )

      CALL nccheck( nf90_def_var(ncid, 'lake_lev', nf90_int, (/lakelevid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name',"lake level") )

      CALL nccheck( nf90_def_var(ncid, 'band', nf90_int, (/bandid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name',"band (vis/nir)") )

      ! grid mask
      CALL nccheck( nf90_def_var(ncid, 'landmask', nf90_int, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','grid mask') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','none') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', -1) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', -1) )

      ! grid total fraction
      CALL nccheck( nf90_def_var(ncid, 'landfrac', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','grid total fraction') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','fraction') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! grid cell area [km2]
      CALL nccheck( nf90_def_var(ncid, 'area', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','grid cell area [km2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','fraction') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! wind stress: E-W [kg/m/s2]
      CALL nccheck( nf90_def_var(ncid, 'f_taux', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','wind stress: E-W [kg/m/s2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','kg/m/s2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! wind stress: N-S [kg/m/s2]
      CALL nccheck( nf90_def_var(ncid, 'f_tauy', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','wind stress: N-S [kg/m/s2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','kg/m/s2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! sensible heat from canopy height to atmosphere [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_fsena', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','sensible heat from canopy height to atmosphere [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      CALL nccheck( nf90_def_var(ncid, 'f_fsenroof', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','sensible heat from urban roof [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      CALL nccheck( nf90_def_var(ncid, 'f_fsenwsun', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','sensible heat from urban sunwall [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      CALL nccheck( nf90_def_var(ncid, 'f_fsenwsha', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','sensible heat from urban shadedwall [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      CALL nccheck( nf90_def_var(ncid, 'f_fsengimp', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','sensible heat from urban imperivous ground [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      CALL nccheck( nf90_def_var(ncid, 'f_fsengper', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','sensible heat from urban perivous ground [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      CALL nccheck( nf90_def_var(ncid, 'f_fsenurl', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','sensible heat from urban tree [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! latent heat flux from canopy height to atmosphere [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_lfevpa', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','latent heat flux from canopy height to atmosphere [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      CALL nccheck( nf90_def_var(ncid, 'f_lfevproof', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','latent heat flux from urban roof [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      CALL nccheck( nf90_def_var(ncid, 'f_lfevpgimp', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','latent heat flux from urban imperivous ground [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      CALL nccheck( nf90_def_var(ncid, 'f_lfevpgper', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','latent heat flux from urban perivous ground [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      CALL nccheck( nf90_def_var(ncid, 'f_lfevpurl', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','latent heat flux from urban tree [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! evapotranspiration from canopy to atmosphere [mm/s]
      CALL nccheck( nf90_def_var(ncid, 'f_fevpa', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','evapotranspiration from canopy to atmosphere [mm/s]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','mm/s') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! sensible heat from leaves [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_fsenl', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','sensible heat from leaves [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! evaporation+transpiration from leaves [mm/s]
      CALL nccheck( nf90_def_var(ncid, 'f_fevpl', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','evaporation+transpiration from leaves [mm/s]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','mm/s') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! transpiration rate [mm/s]
      CALL nccheck( nf90_def_var(ncid, 'f_etr', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','transpiration rate [mm/s]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','mm/s') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! sensible heat flux from ground [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_fseng', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','sensible heat flux from ground [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! evaporation heat flux from ground [mm/s]
      CALL nccheck( nf90_def_var(ncid, 'f_fevpg', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','evaporation heat flux from ground [mm/s]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','mm/s') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! ground heat flux [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_fgrnd', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','ground heat flux [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! solar absorbed by sunlit canopy [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_sabvsun', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','solar absorbed by sunlit canopy [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! solar absorbed by shaded [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_sabvsha', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','solar absorbed by shaded [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! solar absorbed by ground  [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_sabg', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','solar absorbed by ground [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! outgoing long-wave radiation from ground+canopy [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_olrg', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','outgoing long-wave radiation from ground+canopy [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! net radiation [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_rnet', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','net radiation [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! the error of water banace [mm/s]
      CALL nccheck( nf90_def_var(ncid, 'f_xerr', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','the error of water banace [mm/s]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','mm/s') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! the error of energy balance [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_zerr', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','the error of energy balance [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! surface runoff [mm/s]
      CALL nccheck( nf90_def_var(ncid, 'f_rsur', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','surface runoff [mm/s]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','mm/s') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! total runoff [mm/s]
      CALL nccheck( nf90_def_var(ncid, 'f_rnof', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','total runoff [mm/s]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','mm/s') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! canopy assimilation rate [mol m-2 s-1]
      CALL nccheck( nf90_def_var(ncid, 'f_assim', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','canopy assimilation rate [mol m-2 s-1]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','mol m-2 s-1') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! respiration (plant+soil) [mol m-2 s-1]
      CALL nccheck( nf90_def_var(ncid, 'f_respc', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','respiration (plant+soil) [mol m-2 s-1]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','mol m-2 s-1') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! groundwater recharge rate [mm/s]
      CALL nccheck( nf90_def_var(ncid, 'f_qcharge', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','groundwater recharge rate [mm/s] ') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','mm/s') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      !---------------------------------------------------------------------
      ! ground surface temperature [K]
      CALL nccheck( nf90_def_var(ncid, 'f_t_grnd', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','ground surface temperature [K]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','K') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! sunlit leaf temperature [K]
      CALL nccheck( nf90_def_var(ncid, 'f_tleaf', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','sunlit leaf temperature [K]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','K') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! depth of water on foliage [mm]
      CALL nccheck( nf90_def_var(ncid, 'f_ldew', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','depth of water on foliage [mm]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','mm') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! snow cover, water equivalent [mm]
      CALL nccheck( nf90_def_var(ncid, 'f_scv', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','snow cover, water equivalent [mm]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','mm') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! snow depth [meter]
      CALL nccheck( nf90_def_var(ncid, 'f_snowdp', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','snow depth [meter]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','meter') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! fraction of snow cover on ground [-]
      CALL nccheck( nf90_def_var(ncid, 'f_fsno', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','fraction of snow cover on ground [-]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','-') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! fraction of veg cover, excluding snow-covered veg [-]
      CALL nccheck( nf90_def_var(ncid, 'f_sigf', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','fraction of veg cover, excluding snow-covered veg [-]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','-') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! leaf greenness [fraction]
      CALL nccheck( nf90_def_var(ncid, 'f_green', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','leaf greenness [fraction]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','fraction') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! leaf area index [m2/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_lai', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','leaf area index [m2/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','m2/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! stem area index [m2/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_sai', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','stem area index [m2/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','m2/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! averaged albedo direct [%]
      CALL nccheck( nf90_def_var(ncid, 'f_albd', nf90_double, (/xid,yid,bandid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','averaged albedo direct [-]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','-') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! averaged albedo diffuse [%]
      CALL nccheck( nf90_def_var(ncid, 'f_albi', nf90_double, (/xid,yid,bandid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','averaged albedo diffuse [-]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','-') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! averaged bulk surface emissivity [-]
      CALL nccheck( nf90_def_var(ncid, 'f_emis', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','averaged bulk surface emissivity [-]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','-') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! effective roughness [m]
      CALL nccheck( nf90_def_var(ncid, 'f_z0m', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','effective roughness [m]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','m') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! radiative temperature of surface [K]
      CALL nccheck( nf90_def_var(ncid, 'f_trad', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','radiative temperature of surface [K]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','K') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! 2 m height air temperature [kelvin]
      CALL nccheck( nf90_def_var(ncid, 'f_tref', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','2 m height air temperature [kelvin]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','kelvin') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! Diurnal Max 2 m height air temperature [kelvin]
      CALL nccheck( nf90_def_var(ncid, 'f_tmax', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','Diurnal Max 2 m height air temperature [kelvin]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','kelvin') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! Diurnal Min 2 m height air temperature [kelvin]
      CALL nccheck( nf90_def_var(ncid, 'f_tmin', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','Diurnal Min 2 m height air temperature [kelvin]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','kelvin') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! DTR of 2 m height air temperature [kelvin]
      CALL nccheck( nf90_def_var(ncid, 'f_tdtr', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','DTR of 2 m height air temperature [kelvin]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','kelvin') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! 2 m height air specific humidity [kg/kg]
      CALL nccheck( nf90_def_var(ncid, 'f_qref', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','2 m height air specific humidity [kg/kg]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','kg/kg') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! rain [mm/s]
      CALL nccheck( nf90_def_var(ncid, 'f_xy_rain', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','rain [mm/s]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','mm/s') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! snow [mm/s]
      CALL nccheck( nf90_def_var(ncid, 'f_xy_snow', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','snow [mm/s]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','mm/s') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! temperature of inner building [K]
      CALL nccheck( nf90_def_var(ncid, 'f_t_room', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','temperature of inner building [K]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','K') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! temperature of outer building [K]
      CALL nccheck( nf90_def_var(ncid, 'f_tafu', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','temperature of outer building [K]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','K') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! sensible flux from heat or cool AC [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_fhac', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','sensible flux from heat or cool AC [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! waste heat flux from heat or cool AC [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_fwst', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','waste heat flux from heat or cool AC [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! flux from inner and outter air exchange [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_fach', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','flux from inner and outter air exchange [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! flux from metabolism [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_fahe', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','flux from people and cars [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! flux from heatiog [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_fhah', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','flux from heating/cooling [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! flux from vehicle [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_vehc', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','flux from cars [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! flux from metabolism [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_meta', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','flux from people [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! solar absorbed by sunlit canopy in daytime [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_sabvdt', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','solar absorbed by sunlit canopy in daytime [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! solar absorbed by ground in daytime [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_sabgdt', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','solar absorbed by ground in daytime [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! reflected solar radiation at surface in daytime [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_srdt', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','reflected solar radiation at surface in daytime [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! sensible heat from canopy height to atmosphere in daytime [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_fsenadt', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','sensible heat from canopy height to atmosphere in daytime [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! latent heat flux from canopy height to atmosphere in daytime [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_lfevpadt', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','latent heat flux from canopy height to atmosphere in daytime [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! ground heat flux in daytime [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_fgrnddt', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','ground heat flux in daytime [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! outgoing long-wave radiation from ground+canopy in daytime [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_olrgdt', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','outgoing long-wave radiation from ground+canopy in daytime [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! net radiation in daytime [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_rnetdt', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','net radiation in daytime [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! ground surface temperature in daytime [K]
      CALL nccheck( nf90_def_var(ncid, 'f_t_grnddt', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','ground surface temperature in daytime [K]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','K') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! radiative temperature of surface in daytime [K]
      CALL nccheck( nf90_def_var(ncid, 'f_traddt', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','radiative temperature of surface in daytime [K]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','K') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! 2 m height air temperature in daytime [kelvin]
      CALL nccheck( nf90_def_var(ncid, 'f_trefdt', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','2 m height air temperature in daytime [kelvin]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','kelvin') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! temperature of outer building in daytime [K]
      CALL nccheck( nf90_def_var(ncid, 'f_tafudt', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','temperature of outer building in daytime [K]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','K') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! temperature of roof [K]
      CALL nccheck( nf90_def_var(ncid, 'f_t_roof', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','temperature of roof [K]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','K') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! temperature of wall [K]
      CALL nccheck( nf90_def_var(ncid, 'f_t_wall', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','temperature of wall [K]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','K') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! sensible heat from canopy height to atmosphere in night-time [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_fsenant', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','sensible heat from canopy height to atmosphere in night-time [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! latent heat flux from canopy height to atmosphere in night-time [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_lfevpant', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','latent heat flux from canopy height to &
                                 atmosphere in night-time [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! ground heat flux in night-time [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_fgrndnt', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','ground heat flux in night-time [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! outgoing long-wave radiation from ground+canopy in night-time [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_olrgnt', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','outgoing long-wave radiation from ground+canopy in night-time [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! net radiation in night-time [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_rnetnt', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','net radiation in night-time [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! ground surface temperature in night-time [K]
      CALL nccheck( nf90_def_var(ncid, 'f_t_grndnt', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','ground surface temperature in night-time [K]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','K') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! radiative temperature of surface in night-time [K]
      CALL nccheck( nf90_def_var(ncid, 'f_tradnt', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','radiative temperature of surface in night-time [K]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','K') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! 2 m height air temperature in night-time [kelvin]
      CALL nccheck( nf90_def_var(ncid, 'f_trefnt', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','2 m height air temperature in night-time [kelvin]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','kelvin') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! temperature of outer building in nighttime [K]
      CALL nccheck( nf90_def_var(ncid, 'f_tafunt', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','temperature of outer building in nighttime [K]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','K') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      !---------------------------------------------------------------------
      ! soil temperature [K]
      CALL nccheck( nf90_def_var(ncid, 'f_t_soisno', nf90_double, (/xid,yid,sslevid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','soil temperature [K]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','K') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! liquid water in soil layers [kg/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_wliq_soisno', nf90_double, (/xid,yid,sslevid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','liquid water in soil layers [kg/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','kg/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! ice lens in soil layers [kg/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_wice_soisno', nf90_double, (/xid,yid,sslevid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','ice lens in soil layers [kg/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','kg/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! lake temperature [K]
      CALL nccheck( nf90_def_var(ncid, 'f_t_lake', nf90_double, (/xid,yid,lakelevid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','lake temperature [K]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','K') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! lake ice fraction cover [0-1]
      CALL nccheck( nf90_def_var(ncid, 'f_lake_icefrac', nf90_double, (/xid,yid,lakelevid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','lake ice fraction cover [0-1]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','0-1') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! u* in similarity theory [m/s]
      CALL nccheck( nf90_def_var(ncid, 'f_ustar', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','u* in similarity theory [m/s]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','m/s') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! t* in similarity theory [kg/kg]
      CALL nccheck( nf90_def_var(ncid, 'f_tstar', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','t* in similarity theory [kg/kg]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','kg/kg') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! q* in similarity theory [kg/kg]
      CALL nccheck( nf90_def_var(ncid, 'f_qstar', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','q* in similarity theory [kg/kg]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','kg/kg') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! dimensionless height (z/L) used in Monin-Obukhov theory [-]
      CALL nccheck( nf90_def_var(ncid, 'f_zol', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','dimensionless height (z/L) used in Monin-Obukhov theory [-]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','-') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! bulk Richardson number in surface layer [-]
      CALL nccheck( nf90_def_var(ncid, 'f_rib', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','bulk Richardson number in surface layer [-]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','-') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! integral of profile function for momentum [-]
      CALL nccheck( nf90_def_var(ncid, 'f_fm', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','integral of profile function for momentum [-]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','-') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! integral of profile function for heat [-]
      CALL nccheck( nf90_def_var(ncid, 'f_fh', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','integral of profile function for heat [-]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','-') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! integral of profile function for moisture [-]
      CALL nccheck( nf90_def_var(ncid, 'f_fq', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','integral of profile function for moisture [-]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','-') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! 10m u-velocity [m/s]
      CALL nccheck( nf90_def_var(ncid, 'f_us10m', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','10m u-velocity [m/s]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','m/s') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! 10m v-velocity [m/s]
      CALL nccheck( nf90_def_var(ncid, 'f_vs10m', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','10m v-velocity [m/s]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','m/s') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! integral of profile function for momentum at 10m [-]
      CALL nccheck( nf90_def_var(ncid, 'f_fm10m', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','integral of profile function for momentum at 10m [-]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','-') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! wind in eastward direction [m/s]
      CALL nccheck( nf90_def_var(ncid, 'f_xy_us', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','wind in eastward direction [m/s]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','m/s') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! wind in northward direction [m/s]
      CALL nccheck( nf90_def_var(ncid, 'f_xy_vs', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','wind in northward direction [m/s]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','m/s') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! temperature at reference height [kelvin]
      CALL nccheck( nf90_def_var(ncid, 'f_xy_t', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','temperature at reference height [kelvin]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','kelvin') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! specific humidity at reference height [kg/kg]
      CALL nccheck( nf90_def_var(ncid, 'f_xy_q', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','specific humidity at reference height [kg/kg]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','kg/kg') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! convective precipitation [mm/s]
      CALL nccheck( nf90_def_var(ncid, 'f_xy_prc', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','convective precipitation [mm/s]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','mm/s') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! large scale precipitation [mm/s]
      CALL nccheck( nf90_def_var(ncid, 'f_xy_prl', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','large scale precipitation [mm/s]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','mm/s') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! atmospheric pressure at the surface [pa]
      CALL nccheck( nf90_def_var(ncid, 'f_xy_pbot', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','atmospheric pressure at the surface [pa]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','pa') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! atmospheric infrared (longwave) radiation [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_xy_frl', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','atmospheric infrared (longwave) radiation [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )

      ! downward solar radiation at surface [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_xy_solarin', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','downward solar radiation at surface [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! reflected solar radiation at surface [W/m2]
      CALL nccheck( nf90_def_var(ncid, 'f_sr', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','reflected solar radiation at surface [W/m2]') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! incident direct beam vis solar radiation (W/m2)
      CALL nccheck( nf90_def_var(ncid, 'f_solvd', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','incident direct beam vis solar radiation (W/m2)') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! incident diffuse beam vis solar radiation (W/m2)
      CALL nccheck( nf90_def_var(ncid, 'f_solvi', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','incident diffuse beam vis solar radiation (W/m2)') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! incident direct beam nir solar radiation (W/m2)
      CALL nccheck( nf90_def_var(ncid, 'f_solnd', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','incident direct beam nir solar radiation (W/m2)') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! incident diffuse beam nir solar radiation (W/m2)
      CALL nccheck( nf90_def_var(ncid, 'f_solni', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','incident diffuse beam nir solar radiation (W/m2)') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! reflected direct beam vis solar radiation (W/m2)
      CALL nccheck( nf90_def_var(ncid, 'f_srvd', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','reflected direct beam vis solar radiation (W/m2)') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! reflected diffuse beam vis solar radiation (W/m2)
      CALL nccheck( nf90_def_var(ncid, 'f_srvi', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','reflected diffuse beam vis solar radiation (W/m2)') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! reflected direct beam nir solar radiation (W/m2)
      CALL nccheck( nf90_def_var(ncid, 'f_srnd', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','reflected direct beam nir solar radiation (W/m2)') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! reflected diffuse beam nir solar radiation (W/m2)
      CALL nccheck( nf90_def_var(ncid, 'f_srni', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','reflected diffuse beam nir solar radiation (W/m2)') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! incident direct beam vis solar radiation at local noon(W/m2)
      CALL nccheck( nf90_def_var(ncid, 'f_solvdln', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','incident direct beam vis solar radiation at local noon(W/m2)') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! incident diffuse beam vis solar radiation at local noon(W/m2)
      CALL nccheck( nf90_def_var(ncid, 'f_solviln', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','incident diffuse beam vis solar radiation at local noon(W/m2)') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! incident direct beam nir solar radiation at local noon(W/m2)
      CALL nccheck( nf90_def_var(ncid, 'f_solndln', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','incident direct beam nir solar radiation at local noon(W/m2)') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! incident diffuse beam nir solar radiation at local noon(W/m2)
      CALL nccheck( nf90_def_var(ncid, 'f_solniln', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','incident diffuse beam nir solar radiation at local noon(W/m2)') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! reflected direct beam vis solar radiation at local noon(W/m2)
      CALL nccheck( nf90_def_var(ncid, 'f_srvdln', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','reflected direct beam vis solar radiation at local noon(W/m2)') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! reflected diffuse beam vis solar radiation at local noon(W/m2)
      CALL nccheck( nf90_def_var(ncid, 'f_srviln', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','reflected diffuse beam vis solar radiation at local noon(W/m2)') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! reflected direct beam nir solar radiation at local noon(W/m2)
      CALL nccheck( nf90_def_var(ncid, 'f_srndln', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','reflected direct beam nir solar radiation at local noon(W/m2)') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )

      ! reflected diffuse beam nir solar radiation at local noon(W/m2)
      CALL nccheck( nf90_def_var(ncid, 'f_srniln', nf90_double, (/xid,yid/), varid) )
      CALL nccheck( nf90_put_att(ncid, varid, 'long_name','reflected diffuse beam nir solar radiation at local noon(W/m2)') )
      CALL nccheck( nf90_put_att(ncid, varid, 'units','W/m2') )
      CALL nccheck( nf90_put_att(ncid, varid, 'missing_value', spval) )
      CALL nccheck( nf90_put_att(ncid, varid, '_FillValue', spval) )


    ! end defination
      CALL nccheck( nf90_enddef(ncid) )

    ! write data
    ! ------------------------------------------------------------

    ! dimension data
      lons = gridlond
      CALL nccheck( nf90_inq_varid(ncid,'lon',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,lons) )

      lats = gridlatd
      CALL nccheck( nf90_inq_varid(ncid,'lat',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,lats) )

      DO ilev = 1, nl_soil-maxsnl
         sslev(ilev) = maxsnl+ilev
      ENDDO
      CALL nccheck( nf90_inq_varid(ncid,'soil_snow_lev',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,sslev) )

      DO ilev = 1, nl_lake
         lake_lev(ilev) = ilev
      ENDDO
      CALL nccheck( nf90_inq_varid(ncid,'lake_lev',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,lake_lev) )

      DO ilev = 1, 2
         band(ilev) = ilev
      ENDDO
      CALL nccheck( nf90_inq_varid(ncid,'band',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,band) )

      ! grid mask
      CALL nccheck( nf90_inq_varid(ncid,'landmask',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,mask) )

      ! grid total fraction
      CALL nccheck( nf90_inq_varid(ncid,'landfrac',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,frac) )

      ! grid cell area [km2]
      CALL nccheck( nf90_inq_varid(ncid,'area',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,area) )

      ! wind stress: E-W [kg/m/s2]
      CALL nccheck( nf90_inq_varid(ncid,'f_taux',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_taux) )

      ! wind stress: N-S [kg/m/s2]
      CALL nccheck( nf90_inq_varid(ncid,'f_tauy',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_tauy) )

      ! sensible heat from canopy height to atmosphere [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_fsena',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fsena) )

      CALL nccheck( nf90_inq_varid(ncid,'f_fsenroof',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_senroof) )

      CALL nccheck( nf90_inq_varid(ncid,'f_fsenwsun',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_senwsun) )

      CALL nccheck( nf90_inq_varid(ncid,'f_fsenwsha',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_senwsha) )

      CALL nccheck( nf90_inq_varid(ncid,'f_fsengimp',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_sengimp) )

      CALL nccheck( nf90_inq_varid(ncid,'f_fsengper',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_sengper) )

      CALL nccheck( nf90_inq_varid(ncid,'f_fsenurbl',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_senurbl) )

      ! latent heat flux from canopy height to atmosphere [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_lfevpa',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_lfevpa) )

      CALL nccheck( nf90_inq_varid(ncid,'f_lfevproof',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_lfevproof) )

      CALL nccheck( nf90_inq_varid(ncid,'f_lfevpgimp',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_lfevpgimp) )

      CALL nccheck( nf90_inq_varid(ncid,'f_lfevpgper',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_lfevpgper) )

      CALL nccheck( nf90_inq_varid(ncid,'f_lfevpurbl',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_lfevpurbl) )

      ! evapotranspiration from canopy to atmosphere [mm/s]
      CALL nccheck( nf90_inq_varid(ncid,'f_fevpa',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fevpa) )

      ! sensible heat from leaves [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_fsenl',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fsenl) )

      ! evaporation+transpiration from leaves [mm/s]
      CALL nccheck( nf90_inq_varid(ncid,'f_fevpl',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fevpl) )

      ! transpiration rate [mm/s]
      CALL nccheck( nf90_inq_varid(ncid,'f_etr',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_etr) )

      ! sensible heat flux from ground [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_fseng',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fseng) )

      ! evaporation heat flux from ground [mm/s]
      CALL nccheck( nf90_inq_varid(ncid,'f_fevpg',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fevpg) )

      ! ground heat flux [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_fgrnd',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fgrnd) )

      ! solar absorbed by sunlit canopy [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_sabvsun',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_sabvsun) )

      ! solar absorbed by shaded [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_sabvsha',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_sabvsha) )

      ! solar absorbed by ground  [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_sabg',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_sabg) )

      ! outgoing long-wave radiation from ground+canopy [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_olrg',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_olrg) )

      ! net radiation [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_rnet',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_rnet) )

      ! the error of water banace [mm/s]
      CALL nccheck( nf90_inq_varid(ncid,'f_xerr',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_xerr) )

      ! the error of energy balance [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_zerr',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_zerr) )

      ! surface runoff [mm/s]
      CALL nccheck( nf90_inq_varid(ncid,'f_rsur',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_rsur) )

      ! total runoff [mm/s]
      CALL nccheck( nf90_inq_varid(ncid,'f_rnof',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_rnof) )

      ! canopy assimilation rate [mol m-2 s-1]
      CALL nccheck( nf90_inq_varid(ncid,'f_assim',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_assim) )

      ! respiration (plant+soil) [mol m-2 s-1]
      CALL nccheck( nf90_inq_varid(ncid,'f_respc',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_respc) )

      ! groundwater recharge rate [mm/s]
      CALL nccheck( nf90_inq_varid(ncid,'f_qcharge',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_qcharge) )

!---------------------------------------------------------------------
      ! ground surface temperature [K]
      CALL nccheck( nf90_inq_varid(ncid,'f_t_grnd',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_t_grnd) )

      ! sunlit leaf temperature [K]
      CALL nccheck( nf90_inq_varid(ncid,'f_tleaf',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_tleaf) )

      ! shaded leaf temperature [K]
      !CALL nccheck( nf90_inq_varid(ncid,'f_tlsha',varid) )
      !CALL nccheck( nf90_put_var(ncid,varid,f_tlsha) )

      ! depth of water on foliage [mm]
      CALL nccheck( nf90_inq_varid(ncid,'f_ldew',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_ldew) )

      ! snow cover, water equivalent [mm]
      CALL nccheck( nf90_inq_varid(ncid,'f_scv',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_scv) )

      ! snow depth [meter]
      CALL nccheck( nf90_inq_varid(ncid,'f_snowdp',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_snowdp) )

      ! fraction of snow cover on ground
      CALL nccheck( nf90_inq_varid(ncid,'f_fsno',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fsno) )

      ! fraction of veg cover, excluding snow-covered veg [-]
      CALL nccheck( nf90_inq_varid(ncid,'f_sigf',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_sigf) )

      ! leaf greenness
      CALL nccheck( nf90_inq_varid(ncid,'f_green',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_green) )

      ! leaf area index
      CALL nccheck( nf90_inq_varid(ncid,'f_lai',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_lai) )

      ! stem area index
      CALL nccheck( nf90_inq_varid(ncid,'f_sai',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_sai) )

      ! averaged albedo direct
      tmp(:,:,1) = f_alb(1,1,:,:)
      tmp(:,:,2) = f_alb(2,1,:,:)
      CALL nccheck( nf90_inq_varid(ncid,'f_albd',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,tmp) )

      ! averaged albedo diffuse
      tmp(:,:,1) = f_alb(1,2,:,:)
      tmp(:,:,2) = f_alb(2,2,:,:)
      CALL nccheck( nf90_inq_varid(ncid,'f_albi',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,tmp) )

      ! averaged bulk surface emissivity
      CALL nccheck( nf90_inq_varid(ncid,'f_emis',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_emis) )

      ! effective roughness [m]
      CALL nccheck( nf90_inq_varid(ncid,'f_z0m',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_z0m) )

      ! radiative temperature of surface [K]
      CALL nccheck( nf90_inq_varid(ncid,'f_trad',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_trad) )

      ! 2 m height air temperature [kelvin]
      CALL nccheck( nf90_inq_varid(ncid,'f_tref',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_tref) )

      ! Diurnal Max 2 m height air temperature [kelvin]
      CALL nccheck( nf90_inq_varid(ncid,'f_tmax',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_tmax) )

      ! Dirnal 2 m height air temperature [kelvin]
      CALL nccheck( nf90_inq_varid(ncid,'f_tmin',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_tmin) )

      ! DTR of 2 m height air temperature [kelvin]
      CALL nccheck( nf90_inq_varid(ncid,'f_tdtr',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_tdtr) )

      ! 2 m height air specific humidity [kg/kg]
      CALL nccheck( nf90_inq_varid(ncid,'f_qref',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_qref) )

      ! rain [mm/s]
      CALL nccheck( nf90_inq_varid(ncid,'f_xy_rain',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_xy_rain) )

      ! snow [mm/s]
      CALL nccheck( nf90_inq_varid(ncid,'f_xy_snow',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_xy_snow) )

      ! temperature of inner building [K]
      CALL nccheck( nf90_inq_varid(ncid,'f_t_room',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_t_room) )

      ! temperature of outer building [K]
      CALL nccheck( nf90_inq_varid(ncid,'f_tafu',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_tafu) )

      ! temperature of roof [K]
      CALL nccheck( nf90_inq_varid(ncid,'f_t_roof',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_troof) )

      ! temperature of wall [K]
      CALL nccheck( nf90_inq_varid(ncid,'f_t_wall',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_twall) )

      ! sensible flux from heat or cool AC [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_fhac',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fhac) )

      ! waste heat flux from heat or cool AC [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_fwst',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fwst) )

      ! flux from inner and outter air exchange [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_fach',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fach) )

      ! flux from metabolism and vehicle [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_fahe',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fahe) )

      ! flux from heating [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_fhah',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fhah) )

      ! flux from vehicle [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_vehc',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_vehc) )

      ! flux from metabolism [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_meta',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_meta) )

      ! solar absorbed by sunlit canopy in daytime [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_sabvdt',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_sabvdt) )

      ! solar absorbed by ground in daytime [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_sabgdt',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_sabgdt) )

      ! reflected solar radiation at surface in daytime [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_srdt',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_srdt) )

      ! sensible heat from canopy height to atmosphere in daytime [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_fsenadt',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fsenadt) )

      ! latent heat flux from canopy height to atmosphere in daytime [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_lfevpadt',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_lfevpadt) )

      ! ground heat flux in daytime [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_fgrnddt',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fgrnddt) )

      ! outgoing long-wave radiation from ground+canopy in daytime [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_olrgdt',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_olrgdt) )

      ! net radiation in daytime [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_rnetdt',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_rnetdt) )

      ! ground surface temperature in daytime [K]
      CALL nccheck( nf90_inq_varid(ncid,'f_t_grnddt',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_t_grnddt) )

      ! radiative temperature of surface in daytime [K]
      CALL nccheck( nf90_inq_varid(ncid,'f_traddt',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_traddt) )

      ! 2 m height air temperature in daytime [kelvin]
      CALL nccheck( nf90_inq_varid(ncid,'f_trefdt',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_trefdt) )

      ! temperature of outer building in daytime [K]
      CALL nccheck( nf90_inq_varid(ncid,'f_tafudt',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_tafudt) )

      ! sensible heat from canopy height to atmosphere in night-time [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_fsenant',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fsenant) )

      ! latent heat flux from canopy height to atmosphere in night-time [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_lfevpant',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_lfevpant) )

      ! ground heat flux in night-time [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_fgrndnt',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fgrndnt) )

      ! outgoing long-wave radiation from ground+canopy in night-time [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_olrgnt',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_olrgnt) )

      ! net radiation in night-time [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_rnetnt',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_rnetnt) )

      ! ground surface temperature in night-time [K]
      CALL nccheck( nf90_inq_varid(ncid,'f_t_grndnt',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_t_grndnt) )

      ! radiative temperature of surface in night-time [K]
      CALL nccheck( nf90_inq_varid(ncid,'f_tradnt',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_tradnt) )

      ! 2 m height air temperature in night-time [kelvin]
      CALL nccheck( nf90_inq_varid(ncid,'f_trefnt',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_trefnt) )

      ! temperature of outer building in nighttime [K]
      CALL nccheck( nf90_inq_varid(ncid,'f_tafunt',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_tafunt) )

!---------------------------------------------------------------------
      ! soil temperature [K]
      DO i = maxsnl+1, nl_soil
         tmp1(:,:,i) = f_t_soisno(i,:,:)
      ENDDO
      CALL nccheck( nf90_inq_varid(ncid,'f_t_soisno',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,tmp1) )

      ! liquid water in soil layers [kg/m2]
      DO i = maxsnl+1, nl_soil
         tmp1(:,:,i) = f_wliq_soisno(i,:,:)
      ENDDO
      CALL nccheck( nf90_inq_varid(ncid,'f_wliq_soisno',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,tmp1) )

      ! ice lens in soil layers [kg/m2]
      DO i = maxsnl+1, nl_soil
         tmp1(:,:,i) = f_wice_soisno(i,:,:)
      ENDDO
      CALL nccheck( nf90_inq_varid(ncid,'f_wice_soisno',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,tmp1) )

      ! lake temperature [K]
      DO i = 1, nl_lake
         tmp2(:,:,i) = f_t_lake(i,:,:)
      ENDDO
      CALL nccheck( nf90_inq_varid(ncid,'f_t_lake',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,tmp2) )

      ! lake ice fraction cover [0-1]
      DO i = 1, nl_lake
         tmp2(:,:,i) = f_lake_icefrac(i,:,:)
      ENDDO
      CALL nccheck( nf90_inq_varid(ncid,'f_lake_icefrac',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,tmp2) )

      ! u* in similarity theory [m/s]
      CALL nccheck( nf90_inq_varid(ncid,'f_ustar',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_ustar) )

      ! t* in similarity theory [kg/kg]
      CALL nccheck( nf90_inq_varid(ncid,'f_tstar',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_tstar) )

      ! q* in similarity theory [kg/kg]
      CALL nccheck( nf90_inq_varid(ncid,'f_qstar',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_qstar) )

      ! dimensionless height (z/L) used in Monin-Obukhov theory
      CALL nccheck( nf90_inq_varid(ncid,'f_zol',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_zol) )

      ! bulk Richardson number in surface layer
      CALL nccheck( nf90_inq_varid(ncid,'f_rib',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_rib) )

      ! integral of profile function for momentum
      CALL nccheck( nf90_inq_varid(ncid,'f_fm',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fm) )

      ! integral of profile function for heat
      CALL nccheck( nf90_inq_varid(ncid,'f_fh',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fh) )

      ! integral of profile function for moisture
      CALL nccheck( nf90_inq_varid(ncid,'f_fq',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fq) )

      ! 10m u-velocity [m/s]
      CALL nccheck( nf90_inq_varid(ncid,'f_us10m',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_us10m) )

      ! 10m v-velocity [m/s]
      CALL nccheck( nf90_inq_varid(ncid,'f_vs10m',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_vs10m) )

      ! integral of profile function for momentum at 10m [-]
      CALL nccheck( nf90_inq_varid(ncid,'f_fm10m',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_fm10m) )

      ! wind in eastward direction [m/s]
      CALL nccheck( nf90_inq_varid(ncid,'f_xy_us',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_xy_us) )

      ! wind in northward direction [m/s]
      CALL nccheck( nf90_inq_varid(ncid,'f_xy_vs',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_xy_vs) )

      ! temperature at reference height [kelvin]
      CALL nccheck( nf90_inq_varid(ncid,'f_xy_t',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_xy_t) )

      ! specific humidity at reference height [kg/kg]
      CALL nccheck( nf90_inq_varid(ncid,'f_xy_q',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_xy_q) )

      ! convective precipitation [mm/s]
      CALL nccheck( nf90_inq_varid(ncid,'f_xy_prc',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_xy_prc) )

      ! large scale precipitation [mm/s]
      CALL nccheck( nf90_inq_varid(ncid,'f_xy_prl',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_xy_prl) )

      ! atmospheric pressure at the surface [pa]
      CALL nccheck( nf90_inq_varid(ncid,'f_xy_pbot',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_xy_pbot) )

      ! atmospheric infrared (longwave) radiation [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_xy_frl',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_xy_frl) )

      ! downward solar radiation at surface [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_xy_solarin',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_xy_solarin) )

      ! total reflected solar radiation at surface [W/m2]
      CALL nccheck( nf90_inq_varid(ncid,'f_sr',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_sr) )

      ! incident direct beam vis solar radiation (W/m2)
      CALL nccheck( nf90_inq_varid(ncid,'f_solvd',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_solvd) )

      ! incident diffuse beam vis solar radiation (W/m2)
      CALL nccheck( nf90_inq_varid(ncid,'f_solvi',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_solvi) )

      ! incident direct beam nir solar radiation (W/m2)
      CALL nccheck( nf90_inq_varid(ncid,'f_solnd',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_solnd) )

      ! incident diffuse beam nir solar radiation (W/m2)
      CALL nccheck( nf90_inq_varid(ncid,'f_solni',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_solni) )

      ! reflected direct beam vis solar radiation (W/m2)
      CALL nccheck( nf90_inq_varid(ncid,'f_srvd',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_srvd) )

      ! reflected diffuse beam vis solar radiation (W/m2)
      CALL nccheck( nf90_inq_varid(ncid,'f_srvi',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_srvi) )

      ! reflected direct beam nir solar radiation (W/m2)
      CALL nccheck( nf90_inq_varid(ncid,'f_srnd',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_srnd) )

      ! reflected diffuse beam nir solar radiation (W/m2)
      CALL nccheck( nf90_inq_varid(ncid,'f_srni',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_srni) )

      ! incident direct beam vis solar radiation at local noon(W/m2)
      CALL nccheck( nf90_inq_varid(ncid,'f_solvdln',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_solvdln) )

      ! incident diffuse beam vis solar radiation at local noon(W/m2)
      CALL nccheck( nf90_inq_varid(ncid,'f_solviln',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_solviln) )

      ! incident direct beam nir solar radiation at local noon(W/m2)
      CALL nccheck( nf90_inq_varid(ncid,'f_solndln',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_solndln) )

      ! incident diffuse beam nir solar radiation at local noon(W/m2)
      CALL nccheck( nf90_inq_varid(ncid,'f_solniln',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_solniln) )

      ! reflected direct beam vis solar radiation at local noon(W/m2)
      CALL nccheck( nf90_inq_varid(ncid,'f_srvdln',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_srvdln) )

      ! reflected diffuse beam vis solar radiation at local noon(W/m2)
      CALL nccheck( nf90_inq_varid(ncid,'f_srviln',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_srviln) )

      ! reflected direct beam nir solar radiation at local noon(W/m2)
      CALL nccheck( nf90_inq_varid(ncid,'f_srndln',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_srndln) )

      ! reflected diffuse beam nir solar radiation at local noon(W/m2)
      CALL nccheck( nf90_inq_varid(ncid,'f_srniln',varid) )
      CALL nccheck( nf90_put_var(ncid,varid,f_srniln) )

      CALL nccheck( nf90_close(ncid) )
#else
     write(luout) gridlond (:)    ! longitude in degree
     write(luout) gridlatd (:)    ! latitude in degree
     write(luout) mask     (:,:)  ! grid mask
     write(luout) frac     (:,:)  ! grid total fraction
     write(luout) area     (:,:)  ! grid cell area

     write(luout) f_taux   (:,:)  ! wind stress: E-W [kg/m/s2]
     write(luout) f_tauy   (:,:)  ! wind stress: N-S [kg/m/s2]
     write(luout) f_fsena  (:,:)  ! sensible heat from canopy height to atmosphere [W/m2]
     write(luout) f_lfevpa (:,:)  ! latent heat flux from canopy height to atmosphere [W/m2]
     write(luout) f_fevpa  (:,:)  ! evapotranspiration from canopy to atmosphere [mm/s]
     write(luout) f_fsenl  (:,:)  ! sensible heat from leaves [W/m2]
     write(luout) f_fevpl  (:,:)  ! evaporation+transpiration from leaves [mm/s]
     write(luout) f_etr    (:,:)  ! transpiration rate [mm/s]
     write(luout) f_fseng  (:,:)  ! sensible heat flux from ground [W/m2]
     write(luout) f_fevpg  (:,:)  ! evaporation heat flux from ground [mm/s]
     write(luout) f_fgrnd  (:,:)  ! ground heat flux [W/m2]
     write(luout) f_sabvsun(:,:)  ! solar absorbed by sunlit canopy [W/m2]
     write(luout) f_sabvsha(:,:)  ! solar absorbed by shaded [W/m2]
     write(luout) f_sabg   (:,:)  ! solar absorbed by ground  [W/m2]
     write(luout) f_olrg   (:,:)  ! outgoing long-wave radiation from ground+canopy [W/m2]
     write(luout) f_rnet   (:,:)  ! net radiation [W/m2]
     write(luout) f_xerr   (:,:)  ! the error of water banace [mm/s]
     write(luout) f_zerr   (:,:)  ! the error of energy balance [W/m2]
     write(luout) f_rsur   (:,:)  ! surface runoff [mm/s]
     write(luout) f_rnof   (:,:)  ! total runoff [mm/s]
     write(luout) f_qintr  (:,:)  ! interception [mm/s]
     write(luout) f_qinfl  (:,:)  ! inflitration [mm/s]
     write(luout) f_qdrip  (:,:)  ! throughfall [mm/s]
     write(luout) f_assim  (:,:)  ! canopy assimilation rate [mol m-2 s-1]
     write(luout) f_respc  (:,:)  ! respiration (plant+soil) [mol m-2 s-1]
     write(luout) f_qcharge(:,:)  ! groundwater recharge rate [mm/s]

!---------------------------------------------------------------------
     write(luout) f_t_grnd (:,:)  ! ground surface temperature [K]
     write(luout) f_tleaf  (:,:)  ! sunlit leaf temperature [K]
     write(luout) f_ldew   (:,:)  ! depth of water on foliage [mm]
     write(luout) f_scv    (:,:)  ! snow cover, water equivalent [mm]
     write(luout) f_snowdp (:,:)  ! snow depth [meter]
     write(luout) f_fsno   (:,:)  ! fraction of snow cover on ground
     write(luout) f_sigf   (:,:)  ! fraction of veg cover, excluding snow-covered veg [-]
     write(luout) f_green  (:,:)  ! leaf greenness
     write(luout) f_lai    (:,:)  ! leaf area index
     write(luout) f_laisun (:,:)  ! sunlit leaf area index
     write(luout) f_laisha (:,:)  ! shaded leaf area index
     write(luout) f_sai    (:,:)  ! stem area index
     write(luout) f_alb(:,:,:,:)  ! averaged albedo [visible, direct; direct, diffuse]
     write(luout) f_emis   (:,:)  ! averaged bulk surface emissivity
     write(luout) f_z0m    (:,:)  ! effective roughness [m]
     write(luout) f_trad   (:,:)  ! radiative temperature of surface [K]
     write(luout) f_tref   (:,:)  ! 2 m height air temperature [kelvin]
     write(luout) f_tmax   (:,:)  ! Diurnal Max 2 m height air temperature [kelvin]
     write(luout) f_tmin   (:,:)  ! Diurnal Min 2 m height air temperature [kelvin]
     write(luout) f_tdtr   (:,:)  ! DTR of 2 m height air temperature [kelvin]
     write(luout) f_qref   (:,:)  ! 2 m height air specific humidity [kg/kg]
     write(luout) f_xy_rain(:,:)  ! rain [mm/s]
     write(luout) f_xy_snow(:,:)  ! snow [mm/s]
     write(luout) f_t_room (:,:)  ! temperature of inner building [K]
     write(luout) f_tafu   (:,:)  ! temperature of outer building [K]
     write(luout) f_fhah   (:,:)  ! flux from heation [W/m2]
     write(luout) f_fhac   (:,:)  ! sensible flux from heat or cool AC [W/m2]
     write(luout) f_fwst   (:,:)  ! waste heat flux from heat or cool AC [W/m2]
     write(luout) f_fach   (:,:)  ! flux from inner and outter air exchange [W/m2]
     write(luout) f_fahe   (:,:)  ! flux from metabolism and vehicle [W/m2]

     write(luout) f_sabvdt  (:,:) ! solar absorbed by sunlit canopy [w/m2]
     write(luout) f_sabgdt  (:,:) ! solar absorbed by ground [w/m2]
     write(luout) f_srdt    (:,:) ! total reflected solar radiation (w/m2)
     write(luout) f_fsenadt (:,:) ! sensible heat from canopy height to atmosphere [w/m2]
     write(luout) f_lfevpadt(:,:) ! latent heat flux from canopy height to atmosphere [w/m2]
     write(luout) f_fgrnddt (:,:) ! ground heat flux [w/m2]
     write(luout) f_olrgdt  (:,:) ! outgoing long-wave radiation from ground+canopy [w/m2]
     write(luout) f_rnetdt  (:,:) ! net radiation [w/m2]
     write(luout) f_t_grnddt(:,:) ! ground surface temperature [k]
     write(luout) f_traddt  (:,:) ! radiative temperature of surface [k]
     write(luout) f_trefdt  (:,:) ! 2 m height air temperature [kelvin]
     write(luout) f_tafudt  (:,:) ! temperature of outer building [K]

     write(luout) f_fsenant (:,:) ! sensible heat from canopy height to atmosphere [w/m2]
     write(luout) f_lfevpant(:,:) ! latent heat flux from canopy height to atmosphere [w/m2]
     write(luout) f_fgrndnt (:,:) ! ground heat flux [w/m2]
     write(luout) f_olrgnt  (:,:) ! outgoing long-wave radiation from ground+canopy [w/m2]
     write(luout) f_rnetnt  (:,:) ! net radiation [w/m2]
     write(luout) f_t_grndnt(:,:) ! ground surface temperature [k]
     write(luout) f_tradnt  (:,:) ! radiative temperature of surface [k]
     write(luout) f_trefnt  (:,:) ! 2 m height air temperature [kelvin]
     write(luout) f_tafunt  (:,:) ! temperature of outer building [K]

!---------------------------------------------------------------------
     write(luout) f_t_soisno   (:,:,:)  ! soil temperature [K]
     write(luout) f_wliq_soisno(:,:,:)  ! liquid water in soil layers [kg/m2]
     write(luout) f_wice_soisno(:,:,:)  ! ice lens in soil layers [kg/m2]
     write(luout) f_h2osoi     (:,:,:)  ! volumetric soil water in layers [m3/m3]
     write(luout) f_rstfac     (:,:)    ! factor of soil water stress
     write(luout) f_zwt        (:,:)    ! the depth to water table [m]
     write(luout) f_wa         (:,:)    ! water storage in aquifer [mm]
     write(luout) f_wat        (:,:)    ! total water storage [mm]

     write(luout) f_t_lake      (:,:,:) ! lake temperature [K]
     write(luout) f_lake_icefrac(:,:,:) ! lake ice fraction cover [0-1]

     write(luout) f_ustar  (:,:)   ! u* in similarity theory [m/s]
     write(luout) f_tstar  (:,:)   ! t* in similarity theory [kg/kg]
     write(luout) f_qstar  (:,:)   ! q* in similarity theory [kg/kg]
     write(luout) f_zol    (:,:)   ! dimensionless height (z/L) used in Monin-Obukhov theory
     write(luout) f_rib    (:,:)   ! bulk Richardson number in surface layer
     write(luout) f_fm     (:,:)   ! integral of profile function for momentum
     write(luout) f_fh     (:,:)   ! integral of profile function for heat
     write(luout) f_fq     (:,:)   ! integral of profile function for moisture
     write(luout) f_us10m  (:,:)   ! 10m u-velocity [m/s]
     write(luout) f_vs10m  (:,:)   ! 10m v-velocity [m/s]
     write(luout) f_fm10m  (:,:)   ! integral of profile function for momentum at 10m [-]

     write(luout) f_xy_us  (:,:)   ! wind in eastward direction [m/s]
     write(luout) f_xy_vs  (:,:)   ! wind in northward direction [m/s]
     write(luout) f_xy_t   (:,:)   ! temperature at reference height [kelvin]
     write(luout) f_xy_q   (:,:)   ! specific humidity at reference height [kg/kg]
     write(luout) f_xy_prc (:,:)   ! convective precipitation [mm/s]
     write(luout) f_xy_prl (:,:)   ! large scale precipitation [mm/s]
     write(luout) f_xy_pbot(:,:)   ! atmospheric pressure at the surface [pa]
     write(luout) f_xy_frl (:,:)   ! atmospheric infrared (longwave) radiation [W/m2]
     write(luout) f_xy_solarin(:,:)   ! downward solar radiation at surface [W/m2]

     write(luout) f_sr     (:,:)  ! total reflected solar radiation (W/m2)
     write(luout) f_solvd  (:,:)  ! incident direct beam vis solar radiation (W/m2)
     write(luout) f_solvi  (:,:)  ! incident diffuse beam vis solar radiation (W/m2)
     write(luout) f_solnd  (:,:)  ! incident direct beam nir solar radiation (W/m2)
     write(luout) f_solni  (:,:)  ! incident diffuse beam nir solar radiation (W/m2)
     write(luout) f_srvd   (:,:)  ! reflected direct beam vis solar radiation (W/m2)
     write(luout) f_srvi   (:,:)  ! reflected diffuse beam vis solar radiation (W/m2)
     write(luout) f_srnd   (:,:)  ! reflected direct beam nir solar radiation (W/m2)
     write(luout) f_srni   (:,:)  ! reflected diffuse beam nir solar radiation (W/m2)
     write(luout) f_solvdln(:,:)  ! incident direct beam vis solar radiation at local noon(W/m2)
     write(luout) f_solviln(:,:)  ! incident diffuse beam vis solar radiation at local noon(W/m2)
     write(luout) f_solndln(:,:)  ! incident direct beam nir solar radiation at local noon(W/m2)
     write(luout) f_solniln(:,:)  ! incident diffuse beam nir solar radiation at local noon(W/m2)
     write(luout) f_srvdln (:,:)  ! reflected direct beam vis solar radiation at local noon(W/m2)
     write(luout) f_srviln (:,:)  ! reflected diffuse beam vis solar radiation at local noon(W/m2)
     write(luout) f_srndln (:,:)  ! reflected direct beam nir solar radiation at local noon(W/m2)
     write(luout) f_srniln (:,:)  ! reflected diffuse beam nir solar radiation at local noon(W/m2)

     CLOSE (luout)
#endif
  END SUBROUTINE flxwrite

! ----------------------------------------------------------------------
! EOP
