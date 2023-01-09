
  !-----------------------------------------------------------------------
  subroutine SnowAlbedo( use_snicar_frc,use_snicar_ad ,coszen_col    ,&
                         albsod        ,albsoi        ,snl           ,frac_sno      ,&
                         h2osno        ,h2osno_liq    ,h2osno_ice    ,snw_rds       ,& 

                         mss_cnc_bcphi ,mss_cnc_bcpho ,mss_cnc_ocphi ,mss_cnc_ocpho ,&
                         mss_cnc_dst1  ,mss_cnc_dst2  ,mss_cnc_dst3  ,mss_cnc_dst4  ,&

                         albgrd        ,albgri        ,albgrd_pur    ,albgri_pur    ,&
                         albgrd_bc     ,albgri_bc     ,albgrd_oc     ,albgri_oc     ,&
                         albgrd_dst    ,albgri_dst    ,flx_absdv     ,flx_absdn     ,&
                         flx_absiv     ,flx_absin      )

  ! !DESCRIPTION:
  ! The calling sequence is:
  ! -> SNICAR_RT:   snow albedos: direct beam (SNICAR)
  !    or
  !    SNICAR_AD_RT: snow albedos: direct beam (SNICAR-AD)
  ! -> SNICAR_RT:   snow albedos: diffuse (SNICAR)
  !    or
  !    SNICAR_AD_RT:   snow albedos: diffuse (SNICAR-AD)
  !
  !-----------------------------------------------------------------------
  ! !USES:
    use SnowSnicarMod , only : SNICAR_RT, SNICAR_AD_RT

  ! and the evolution of snow effective radius
  !
! DAI, Dec. 28, 2022

    IMPLICIT NONE

    integer, parameter :: r8 = selected_real_kind(12) ! 8 byte real

!-------------------------------------------------------------------------
! temporay setting
    integer, parameter :: begc = 1               !  beginning column index
    integer, parameter :: endc = 1               !  beginning and ending column index
!-------------------------------------------------------------------------

    integer, parameter :: numrad  = 2            !  number of solar radiation bands: vis, nir
    integer, parameter :: nlevsno = 5            !  maximum number of snow layers

    integer, parameter :: sno_nbr_aer = 8        !  number of aerosol species in snowpack
    logical, parameter :: DO_SNO_OC   = .true.   !  parameter to include organic carbon (OC)
    logical, parameter :: DO_SNO_AER  = .true.   !  parameter to include aerosols in snowpack radiative calculations
    integer, parameter :: subgridflag = 1        !  = 0 use subgrid fluxes, = 1 not use subgrid fluxes
    !
    ! !ARGUMENTS:
    !
    logical , INTENT(in) :: use_snicar_frc       !  true : if radiative forcing is being calculated, first estimate clean-snow albedo
    logical , INTENT(in) :: use_snicar_ad        !  true: use SNICAR_AD_RT, false: use SNICAR_RT

    real(r8), INTENT(in) :: coszen_col    ( begc:endc )  !  cosine of solar zenith angle
    real(r8), INTENT(in) :: albsod        ( begc:endc,numrad )  !  direct-beam soil albedo (col,bnd) [frc]
    real(r8), INTENT(in) :: albsoi        ( begc:endc,numrad )  !  diffuse soil albedo (col,bnd) [frc]

    integer , INTENT(in) :: snl           ( begc:endc )  !  negative number of snow layers (col) [nbr]
    real(r8), INTENT(in) :: frac_sno      ( begc:endc )  ! fraction of ground covered by snow (0 to 1)
    real(r8), INTENT(in) :: h2osno        ( begc:endc )  ! snow water equivalent (mm H2O)
    real(r8), INTENT(in) :: h2osno_liq    ( begc:endc,-nlevsno+1:0 )  ! liquid water content (col,lyr) [kg/m2]
    real(r8), INTENT(in) :: h2osno_ice    ( begc:endc,-nlevsno+1:0 )  ! ice lens content (col,lyr) [kg/m2]
    real(r8), INTENT(in) :: snw_rds       ( begc:endc,-nlevsno+1:0 )  ! snow grain radius (col,lyr) [microns]

    real(r8), INTENT(in) :: mss_cnc_bcphi ( begc:endc,-nlevsno+1:0 )  !  mass concentration of hydrophilic BC (col,lyr) [kg/kg]
    real(r8), INTENT(in) :: mss_cnc_bcpho ( begc:endc,-nlevsno+1:0 )  !  mass concentration of hydrophobic BC (col,lyr) [kg/kg]
    real(r8), INTENT(in) :: mss_cnc_ocphi ( begc:endc,-nlevsno+1:0 )  !  mass concentration of hydrophilic OC (col,lyr) [kg/kg]
    real(r8), INTENT(in) :: mss_cnc_ocpho ( begc:endc,-nlevsno+1:0 )  !  mass concentration of hydrophobic OC (col,lyr) [kg/kg]
    real(r8), INTENT(in) :: mss_cnc_dst1  ( begc:endc,-nlevsno+1:0 )  !  mass concentration of dust aerosol species 1 (col,lyr) [kg/kg]
    real(r8), INTENT(in) :: mss_cnc_dst2  ( begc:endc,-nlevsno+1:0 )  !  mass concentration of dust aerosol species 2 (col,lyr) [kg/kg]
    real(r8), INTENT(in) :: mss_cnc_dst3  ( begc:endc,-nlevsno+1:0 )  !  mass concentration of dust aerosol species 3 (col,lyr) [kg/kg]
    real(r8), INTENT(in) :: mss_cnc_dst4  ( begc:endc,-nlevsno+1:0 )  !  mass concentration of dust aerosol species 4 (col,lyr) [kg/kg]

    real(r8), INTENT(out) :: albgrd       ( begc:endc,numrad )  !  ground albedo (direct)
    real(r8), INTENT(out) :: albgri       ( begc:endc,numrad )  !  ground albedo (diffuse)
    real(r8), INTENT(out) :: albgrd_pur   ( begc:endc,numrad )  !  pure snow ground albedo (direct)
    real(r8), INTENT(out) :: albgri_pur   ( begc:endc,numrad )  !  pure snow ground albedo (diffuse)
    real(r8), INTENT(out) :: albgrd_bc    ( begc:endc,numrad )  !  ground albedo without BC (direct)
    real(r8), INTENT(out) :: albgri_bc    ( begc:endc,numrad )  !  ground albedo without BC (diffuse)
    real(r8), INTENT(out) :: albgrd_oc    ( begc:endc,numrad )  !  ground albedo without OC (direct)
    real(r8), INTENT(out) :: albgri_oc    ( begc:endc,numrad )  !  ground albedo without OC (diffuse)
    real(r8), INTENT(out) :: albgrd_dst   ( begc:endc,numrad )  !  ground albedo without dust (direct)
    real(r8), INTENT(out) :: albgri_dst   ( begc:endc,numrad )  !  ground albedo without dust (diffuse)
    real(r8), INTENT(out) :: flx_absdv    ( begc:endc,-nlevsno+1:1 )  !  direct flux absorption factor (col,lyr): VIS [frc]
    real(r8), INTENT(out) :: flx_absdn    ( begc:endc,-nlevsno+1:1 )  !  direct flux absorption factor (col,lyr): NIR [frc]
    real(r8), INTENT(out) :: flx_absiv    ( begc:endc,-nlevsno+1:1 )  !  diffuse flux absorption factor (col,lyr): VIS [frc]
    real(r8), INTENT(out) :: flx_absin    ( begc:endc,-nlevsno+1:1 )  !  diffuse flux absorption factor (col,lyr): NIR [frc]

  !-----------------------------------------------------------------------
    !
    ! !LOCAL VARIABLES:
    integer  :: i            ! index for layers [idx]
    integer  :: aer          ! index for sno_nbr_aer
    integer  :: c            ! column indices
    integer  :: ib           ! band index
    integer  :: ic           ! 0=unit incoming direct; 1=unit incoming diffuse
    integer  :: flg_slr      ! flag for SNICAR (=1 if direct, =2 if diffuse)
    integer  :: flg_snw_ice  ! flag for SNICAR (=1 when called from ELM, =2 when called from sea-ice)

    real(r8) :: mss_cnc_aer_in_frc_pur (begc:endc,-nlevsno+1:0,sno_nbr_aer) ! mass concentration of aerosol species for forcing calculation (zero) (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_frc_bc  (begc:endc,-nlevsno+1:0,sno_nbr_aer) ! mass concentration of aerosol species for BC forcing (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_frc_oc  (begc:endc,-nlevsno+1:0,sno_nbr_aer) ! mass concentration of aerosol species for OC forcing (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_frc_dst (begc:endc,-nlevsno+1:0,sno_nbr_aer) ! mass concentration of aerosol species for dust forcing (col,lyr,aer) [kg kg-1]
    real(r8) :: mss_cnc_aer_in_fdb     (begc:endc,-nlevsno+1:0,sno_nbr_aer) ! mass concentration of all aerosol species for feedback calculation (col,lyr,aer) [kg kg-1]

    real(r8) :: albsfc       (begc:endc,numrad)               ! albedo of surface underneath snow (col,bnd)
    real(r8) :: albsnd       (begc:endc,numrad)               ! snow albedo (direct)
    real(r8) :: albsni       (begc:endc,numrad)               ! snow albedo (diffuse)
    real(r8) :: albsnd_pur   (begc:endc,numrad)               ! direct pure snow albedo (radiative forcing)
    real(r8) :: albsni_pur   (begc:endc,numrad)               ! diffuse pure snow albedo (radiative forcing)
    real(r8) :: albsnd_bc    (begc:endc,numrad)               ! direct snow albedo without BC (radiative forcing)
    real(r8) :: albsni_bc    (begc:endc,numrad)               ! diffuse snow albedo without BC (radiative forcing)
    real(r8) :: albsnd_oc    (begc:endc,numrad)               ! direct snow albedo without OC (radiative forcing)
    real(r8) :: albsni_oc    (begc:endc,numrad)               ! diffuse snow albedo without OC (radiative forcing)
    real(r8) :: albsnd_dst   (begc:endc,numrad)               ! direct snow albedo without dust (radiative forcing)
    real(r8) :: albsni_dst   (begc:endc,numrad)               ! diffuse snow albedo without dust (radiative forcing)
    real(r8) :: flx_absd_snw (begc:endc,-nlevsno+1:1,numrad)  ! flux absorption factor for just snow (direct) [frc]
    real(r8) :: flx_absi_snw (begc:endc,-nlevsno+1:1,numrad)  ! flux absorption factor for just snow (diffuse) [frc]
    real(r8) :: foo_snw      (begc:endc,-nlevsno+1:1,numrad)  ! dummy array for forcing calls

    integer  :: snw_rds_in   (begc:endc,-nlevsno+1:0)         ! snow grain size sent to SNICAR (col,lyr) [microns]

    integer , parameter :: nband =numrad   ! number of solar radiation waveband classes

  !-----------------------------------------------------------------------

    ! Initialize output because solar radiation only done if coszen > 0

    do ib = 1, numrad
       do c = begc, endc
          albgrd(c,ib)     = 0._r8
          albgri(c,ib)     = 0._r8
          albgrd_pur(c,ib) = 0._r8
          albgri_pur(c,ib) = 0._r8
          albgrd_bc(c,ib)  = 0._r8
          albgri_bc(c,ib)  = 0._r8
          albgrd_oc(c,ib)  = 0._r8
          albgri_oc(c,ib)  = 0._r8
          albgrd_dst(c,ib) = 0._r8
          albgri_dst(c,ib) = 0._r8
          do i=-nlevsno+1,1,1
             flx_absdv(c,i) = 0._r8
             flx_absdn(c,i) = 0._r8
             flx_absiv(c,i) = 0._r8
             flx_absin(c,i) = 0._r8
          enddo
       end do
    end do  ! end of numrad loop

    ! set variables to pass to SNICAR.

    flg_snw_ice = 1   
    do c = begc, endc
       albsfc(c,:)     = albsoi(c,:)
       snw_rds_in(c,:) = nint(snw_rds(c,:))
    end do

    ! zero aerosol input arrays
    do aer = 1, sno_nbr_aer
       do i = -nlevsno+1, 0
          do c = begc, endc
             mss_cnc_aer_in_frc_pur(c,i,aer) = 0._r8
             mss_cnc_aer_in_frc_bc(c,i,aer)  = 0._r8
             mss_cnc_aer_in_frc_oc(c,i,aer)  = 0._r8
             mss_cnc_aer_in_frc_dst(c,i,aer) = 0._r8
             mss_cnc_aer_in_fdb(c,i,aer)     = 0._r8
          end do
       end do
    end do

    ! If radiative forcing is being calculated, first estimate clean-snow albedo

    if (use_snicar_frc) then

       ! 1. PURE SNOW ALBEDO CALCULATIONS
          flg_slr = 1  ! direct-beam
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col(begc:endc), &
                             snl(begc:endc), &
                             h2osno(begc:endc), &
                             frac_sno(begc:endc), &
                             h2osno_liq(begc:endc, :), &
                             h2osno_ice(begc:endc, :), &
                             snw_rds_in(begc:endc, :), &
                             mss_cnc_aer_in_frc_pur(begc:endc, :, :), &
                             albsfc(begc:endc, :), &
                             albsnd_pur(begc:endc, :), &
                             foo_snw(begc:endc, :, :) )
          else
              call SNICAR_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col(begc:endc), &
                             snl(begc:endc), &
                             h2osno(begc:endc), &
                             frac_sno(begc:endc), &
                             h2osno_liq(begc:endc, :), &
                             h2osno_ice(begc:endc, :), &
                             snw_rds_in(begc:endc, :), &
                             mss_cnc_aer_in_frc_pur(begc:endc, :, :), &
                             albsfc(begc:endc, :), &
                             albsnd_pur(begc:endc, :), &
                             foo_snw(begc:endc, :, :) )
          endif ! end if use_snicar_ad

          flg_slr = 2  ! diffuse
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col(begc:endc), &
                             snl(begc:endc), &
                             h2osno(begc:endc), &
                             frac_sno(begc:endc), &
                             h2osno_liq(begc:endc, :), &
                             h2osno_ice(begc:endc, :), &
                             snw_rds_in(begc:endc, :), &
                             mss_cnc_aer_in_frc_pur(begc:endc, :, :), &
                             albsfc(begc:endc, :), &
                             albsni_pur(begc:endc, :), &
                             foo_snw(begc:endc, :, :) )
          else
              call SNICAR_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col(begc:endc), &
                             snl(begc:endc), &
                             h2osno(begc:endc), &
                             frac_sno(begc:endc), &
                             h2osno_liq(begc:endc, :), &
                             h2osno_ice(begc:endc, :), &
                             snw_rds_in(begc:endc, :), &
                             mss_cnc_aer_in_frc_pur(begc:endc, :, :), &
                             albsfc(begc:endc, :), &
                             albsni_pur(begc:endc, :), &
                             foo_snw(begc:endc, :, :) )
          endif ! end if use_snicar_ad

       ! 2. BC input array:
       !  set dust and (optionally) OC concentrations, so BC_FRC=[(BC+OC+dust)-(OC+dust)]
       if (DO_SNO_OC) then
          mss_cnc_aer_in_frc_bc(begc:endc,:,3) = mss_cnc_ocphi(begc:endc,:)
          mss_cnc_aer_in_frc_bc(begc:endc,:,4) = mss_cnc_ocpho(begc:endc,:)
       endif
       mss_cnc_aer_in_frc_bc(begc:endc,:,5) = mss_cnc_dst1(begc:endc,:)
       mss_cnc_aer_in_frc_bc(begc:endc,:,6) = mss_cnc_dst2(begc:endc,:)
       mss_cnc_aer_in_frc_bc(begc:endc,:,7) = mss_cnc_dst3(begc:endc,:)
       mss_cnc_aer_in_frc_bc(begc:endc,:,8) = mss_cnc_dst4(begc:endc,:)

       ! BC FORCING CALCULATIONS
       flg_slr = 1  ! direct-beam
       if (use_snicar_ad) then
           call SNICAR_AD_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col(begc:endc), &
                             snl(begc:endc), &
                             h2osno(begc:endc), &
                             frac_sno(begc:endc), &
                             h2osno_liq(begc:endc, :), &
                             h2osno_ice(begc:endc, :), &
                             snw_rds_in(begc:endc, :), &
                             mss_cnc_aer_in_frc_bc(begc:endc, :, :), &
                             albsfc(begc:endc, :), &
                             albsnd_bc(begc:endc, :), &
                             foo_snw(begc:endc, :, :) )
       else
           call SNICAR_RT   (flg_snw_ice, &
                             flg_slr, &
                             coszen_col(begc:endc), &
                             snl(begc:endc), &
                             h2osno(begc:endc), &
                             frac_sno(begc:endc), &
                             h2osno_liq(begc:endc, :), &
                             h2osno_ice(begc:endc, :), &
                             snw_rds_in(begc:endc, :), &
                             mss_cnc_aer_in_frc_bc(begc:endc, :, :), &
                             albsfc(begc:endc, :), &
                             albsnd_bc(begc:endc, :), &
                             foo_snw(begc:endc, :, :) )
       endif ! end if use_snicar_ad

       flg_slr = 2  ! diffuse
       if (use_snicar_ad) then
           call SNICAR_AD_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col(begc:endc), &
                             snl(begc:endc), &
                             h2osno(begc:endc), &
                             frac_sno(begc:endc), &
                             h2osno_liq(begc:endc, :), &
                             h2osno_ice(begc:endc, :), &
                             snw_rds_in(begc:endc, :), &
                             mss_cnc_aer_in_frc_bc(begc:endc, :, :), &
                             albsfc(begc:endc, :), &
                             albsni_bc(begc:endc, :), &
                             foo_snw(begc:endc, :, :) )
       else
           call SNICAR_RT   (flg_snw_ice, &
                             flg_slr, &
                             coszen_col(begc:endc), &
                             snl(begc:endc), &
                             h2osno(begc:endc), &
                             frac_sno(begc:endc), &
                             h2osno_liq(begc:endc, :), &
                             h2osno_ice(begc:endc, :), &
                             snw_rds_in(begc:endc, :), &
                             mss_cnc_aer_in_frc_bc(begc:endc, :, :), &
                             albsfc(begc:endc, :), &
                             albsni_bc(begc:endc, :), &
                             foo_snw(begc:endc, :, :) )
       endif ! end if use_snicar_ad

       ! 3. OC input array:
       !  set BC and dust concentrations, so OC_FRC=[(BC+OC+dust)-(BC+dust)]
       if (DO_SNO_OC) then
          mss_cnc_aer_in_frc_oc(begc:endc,:,1) = mss_cnc_bcphi(begc:endc,:)
          mss_cnc_aer_in_frc_oc(begc:endc,:,2) = mss_cnc_bcpho(begc:endc,:)

          mss_cnc_aer_in_frc_oc(begc:endc,:,5) = mss_cnc_dst1(begc:endc,:)
          mss_cnc_aer_in_frc_oc(begc:endc,:,6) = mss_cnc_dst2(begc:endc,:)
          mss_cnc_aer_in_frc_oc(begc:endc,:,7) = mss_cnc_dst3(begc:endc,:)
          mss_cnc_aer_in_frc_oc(begc:endc,:,8) = mss_cnc_dst4(begc:endc,:)

       ! OC FORCING CALCULATIONS
          flg_slr = 1  ! direct-beam
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col(begc:endc), &
                             snl(begc:endc), &
                             h2osno(begc:endc), &
                             frac_sno(begc:endc), &
                             h2osno_liq(begc:endc, :), &
                             h2osno_ice(begc:endc, :), &
                             snw_rds_in(begc:endc, :), &
                             mss_cnc_aer_in_frc_oc(begc:endc, :, :), &
                             albsfc(begc:endc, :), &
                             albsnd_oc(begc:endc, :), &
                             foo_snw(begc:endc, :, :) )
          else
              call SNICAR_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col(begc:endc), &
                             snl(begc:endc), &
                             h2osno(begc:endc), &
                             frac_sno(begc:endc), &
                             h2osno_liq(begc:endc, :), &
                             h2osno_ice(begc:endc, :), &
                             snw_rds_in(begc:endc, :), &
                             mss_cnc_aer_in_frc_oc(begc:endc, :, :), &
                             albsfc(begc:endc, :), &
                             albsnd_oc(begc:endc, :), &
                             foo_snw(begc:endc, :, :) )
          endif ! end if use_snicar_ad

          flg_slr = 2  ! diffuse
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col(begc:endc), &
                             snl(begc:endc), &
                             h2osno(begc:endc), &
                             frac_sno(begc:endc), &
                             h2osno_liq(begc:endc, :), &
                             h2osno_ice(begc:endc, :), &
                             snw_rds_in(begc:endc, :), &
                             mss_cnc_aer_in_frc_oc(begc:endc, :, :), &
                             albsfc(begc:endc, :), &
                             albsni_oc(begc:endc, :), &
                             foo_snw(begc:endc, :, :) )
          else
              call SNICAR_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col(begc:endc), &
                             snl(begc:endc), &
                             h2osno(begc:endc), &
                             frac_sno(begc:endc), &
                             h2osno_liq(begc:endc, :), &
                             h2osno_ice(begc:endc, :), &
                             snw_rds_in(begc:endc, :), &
                             mss_cnc_aer_in_frc_oc(begc:endc, :, :), &
                             albsfc(begc:endc, :), &
                             albsni_oc(begc:endc, :), &
                             foo_snw(begc:endc, :, :) )
          endif ! end if use_snicar_ad
       endif  ! end if (DO_SNO_OC) 

       ! 4. DUST FORCING CALCULATIONS
          ! DUST input array:
          ! set BC and OC concentrations, so DST_FRC=[(BC+OC+dust)-(BC+OC)]
          mss_cnc_aer_in_frc_dst(begc:endc,:,1) = mss_cnc_bcphi(begc:endc,:)
          mss_cnc_aer_in_frc_dst(begc:endc,:,2) = mss_cnc_bcpho(begc:endc,:)

          if (DO_SNO_OC) then
              mss_cnc_aer_in_frc_dst(begc:endc,:,3) = mss_cnc_ocphi(begc:endc,:)
              mss_cnc_aer_in_frc_dst(begc:endc,:,4) = mss_cnc_ocpho(begc:endc,:)
          endif

          flg_slr = 1  ! direct-beam
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col(begc:endc), &
                             snl(begc:endc), &
                             h2osno(begc:endc), &
                             frac_sno(begc:endc), &
                             h2osno_liq(begc:endc, :), &
                             h2osno_ice(begc:endc, :), &
                             snw_rds_in(begc:endc, :), &
                             mss_cnc_aer_in_frc_dst(begc:endc, :, :), &
                             albsfc(begc:endc, :), &
                             albsnd_dst(begc:endc, :), &
                             foo_snw(begc:endc, :, :) )
          else
              call SNICAR_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col(begc:endc), &
                             snl(begc:endc), &
                             h2osno(begc:endc), &
                             frac_sno(begc:endc), &
                             h2osno_liq(begc:endc, :), &
                             h2osno_ice(begc:endc, :), &
                             snw_rds_in(begc:endc, :), &
                             mss_cnc_aer_in_frc_dst(begc:endc, :, :), &
                             albsfc(begc:endc, :), &
                             albsnd_dst(begc:endc, :), &
                             foo_snw(begc:endc, :, :) )
          endif ! end if use_snicar_ad

          flg_slr = 2  ! diffuse
          if (use_snicar_ad) then
              call SNICAR_AD_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col(begc:endc), &
                             snl(begc:endc), &
                             h2osno(begc:endc), &
                             frac_sno(begc:endc), &
                             h2osno_liq(begc:endc, :), &
                             h2osno_ice(begc:endc, :), &
                             snw_rds_in(begc:endc, :), &
                             mss_cnc_aer_in_frc_dst(begc:endc, :, :), &
                             albsfc(begc:endc, :), &
                             albsni_dst(begc:endc, :), &
                             foo_snw(begc:endc, :, :) )
          else
              call SNICAR_RT(flg_snw_ice, &
                             flg_slr, &
                             coszen_col(begc:endc), &
                             snl(begc:endc), &
                             h2osno(begc:endc), &
                             frac_sno(begc:endc), &
                             h2osno_liq(begc:endc, :), &
                             h2osno_ice(begc:endc, :), &
                             snw_rds_in(begc:endc, :), &
                             mss_cnc_aer_in_frc_dst(begc:endc, :, :), &
                             albsfc(begc:endc, :), &
                             albsni_dst(begc:endc, :), &
                             foo_snw(begc:endc, :, :)  )
          endif ! end if use_snicar_ad

    end if !end if use_snicar_frc


    ! --------------------------------------------
    ! CLIMATE FEEDBACK CALCULATIONS, ALL AEROSOLS:
    ! --------------------------------------------
    ! Set aerosol input arrays
    ! feedback input arrays have been zeroed
    ! set soot and dust aerosol concentrations:
    if (DO_SNO_AER) then
        mss_cnc_aer_in_fdb(begc:endc,:,1) = mss_cnc_bcphi(begc:endc,:)
        mss_cnc_aer_in_fdb(begc:endc,:,2) = mss_cnc_bcpho(begc:endc,:)

        ! DO_SNO_OC is set in SNICAR_varpar. Default case is to ignore OC concentrations because:
        !  1) Knowledge of their optical properties is primitive
        !  2) When 'water-soluble' OPAC optical properties are applied to OC in snow,
        !     it has a negligible darkening effect.
        if (DO_SNO_OC) then
           mss_cnc_aer_in_fdb(begc:endc,:,3) = mss_cnc_ocphi(begc:endc,:)
           mss_cnc_aer_in_fdb(begc:endc,:,4) = mss_cnc_ocpho(begc:endc,:)
        endif

        mss_cnc_aer_in_fdb(begc:endc,:,5) = mss_cnc_dst1(begc:endc,:)
        mss_cnc_aer_in_fdb(begc:endc,:,6) = mss_cnc_dst2(begc:endc,:)
        mss_cnc_aer_in_fdb(begc:endc,:,7) = mss_cnc_dst3(begc:endc,:)
        mss_cnc_aer_in_fdb(begc:endc,:,8) = mss_cnc_dst4(begc:endc,:)
    endif

    flg_slr = 1  ! direct-beam
    if (use_snicar_ad) then
        call SNICAR_AD_RT(flg_snw_ice, &
                          flg_slr, &
                          coszen_col(begc:endc), &
                          snl(begc:endc), &
                          h2osno(begc:endc), &
                          frac_sno(begc:endc), &
                          h2osno_liq(begc:endc, :), &
                          h2osno_ice(begc:endc, :), &
                          snw_rds_in(begc:endc, :), &
                          mss_cnc_aer_in_fdb(begc:endc, :, :), &
                          albsfc(begc:endc, :), &
                          albsnd(begc:endc, :), &
                          flx_absd_snw(begc:endc, :, :) )
    else
        call SNICAR_RT   (flg_snw_ice, &
                          flg_slr, &
                          coszen_col(begc:endc), &
                          snl(begc:endc), &
                          h2osno(begc:endc), &
                          frac_sno(begc:endc), &
                          h2osno_liq(begc:endc, :), &
                          h2osno_ice(begc:endc, :), &
                          snw_rds_in(begc:endc, :), &
                          mss_cnc_aer_in_fdb(begc:endc, :, :), &
                          albsfc(begc:endc, :), &
                          albsnd(begc:endc, :), &
                          flx_absd_snw(begc:endc, :, :) )
    endif ! end if use_snicar_ad

    flg_slr = 2  ! diffuse
    if (use_snicar_ad) then
        call SNICAR_AD_RT(flg_snw_ice, &
                          flg_slr, &
                          coszen_col(begc:endc), &
                          snl(begc:endc), &
                          h2osno(begc:endc), &
                          frac_sno(begc:endc), &
                          h2osno_liq(begc:endc, :), &
                          h2osno_ice(begc:endc, :), &
                          snw_rds_in(begc:endc, :), &
                          mss_cnc_aer_in_fdb(begc:endc, :, :), &
                          albsfc(begc:endc, :), &
                          albsni(begc:endc, :), &
                          flx_absi_snw(begc:endc, :, :) )
    else
        call SNICAR_RT   (flg_snw_ice, &
                          flg_slr, &
                          coszen_col(begc:endc), &
                          snl(begc:endc), &
                          h2osno(begc:endc), &
                          frac_sno(begc:endc), &
                          h2osno_liq(begc:endc, :), &
                          h2osno_ice(begc:endc, :), &
                          snw_rds_in(begc:endc, :), &
                          mss_cnc_aer_in_fdb(begc:endc, :, :), &
                          albsfc(begc:endc, :), &
                          albsni(begc:endc, :), &
                          flx_absi_snw(begc:endc, :, :) )
    endif ! end if use_snicar_ad


    ! ground albedos and snow-fraction weighting of snow absorption factors
    do ib = 1, nband
       do c = begc, endc
          if (coszen_col(c) > 0._r8) then
             ! ground albedo was originally computed in SoilAlbedo, but is now computed here
             ! because the order of SoilAlbedo and SNICAR_RT/SNICAR_AD_RT was switched for SNICAR/SNICAR_AD_RT.
             albgrd(c,ib) = albsod(c,ib)*(1._r8-frac_sno(c)) + albsnd(c,ib)*frac_sno(c)
             albgri(c,ib) = albsoi(c,ib)*(1._r8-frac_sno(c)) + albsni(c,ib)*frac_sno(c)

             ! albedos for radiative forcing calculations:
             if (use_snicar_frc) then
                ! pure snow albedo for all-aerosol radiative forcing
                albgrd_pur(c,ib) = albsod(c,ib)*(1.-frac_sno(c)) + albsnd_pur(c,ib)*frac_sno(c)
                albgri_pur(c,ib) = albsoi(c,ib)*(1.-frac_sno(c)) + albsni_pur(c,ib)*frac_sno(c)

                ! BC forcing albedo
                albgrd_bc(c,ib) = albsod(c,ib)*(1.-frac_sno(c)) + albsnd_bc(c,ib)*frac_sno(c)
                albgri_bc(c,ib) = albsoi(c,ib)*(1.-frac_sno(c)) + albsni_bc(c,ib)*frac_sno(c)

                if (DO_SNO_OC) then
                   ! OC forcing albedo
                   albgrd_oc(c,ib) = albsod(c,ib)*(1.-frac_sno(c)) + albsnd_oc(c,ib)*frac_sno(c)
                   albgri_oc(c,ib) = albsoi(c,ib)*(1.-frac_sno(c)) + albsni_oc(c,ib)*frac_sno(c)
                endif

                ! dust forcing albedo
                albgrd_dst(c,ib) = albsod(c,ib)*(1.-frac_sno(c)) + albsnd_dst(c,ib)*frac_sno(c)
                albgri_dst(c,ib) = albsoi(c,ib)*(1.-frac_sno(c)) + albsni_dst(c,ib)*frac_sno(c)
             end if

             ! also in this loop (but optionally in a different loop for vectorized code)
             !  weight snow layer radiative absorption factors based on snow fraction and soil albedo
             !  (NEEDED FOR ENERGY CONSERVATION)
             do i = -nlevsno+1,1,1
                if (subgridflag == 0 ) then
                   if (ib == 1) then
                      flx_absdv(c,i) = flx_absd_snw(c,i,ib)*frac_sno(c) + &
                           ((1.-frac_sno(c))*(1-albsod(c,ib))*(flx_absd_snw(c,i,ib)/(1.-albsnd(c,ib))))
                      flx_absiv(c,i) = flx_absi_snw(c,i,ib)*frac_sno(c) + &
                           ((1.-frac_sno(c))*(1-albsoi(c,ib))*(flx_absi_snw(c,i,ib)/(1.-albsni(c,ib))))
                   elseif (ib == 2) then
                      flx_absdn(c,i) = flx_absd_snw(c,i,ib)*frac_sno(c) + &
                           ((1.-frac_sno(c))*(1-albsod(c,ib))*(flx_absd_snw(c,i,ib)/(1.-albsnd(c,ib))))
                      flx_absin(c,i) = flx_absi_snw(c,i,ib)*frac_sno(c) + &
                           ((1.-frac_sno(c))*(1-albsoi(c,ib))*(flx_absi_snw(c,i,ib)/(1.-albsni(c,ib))))
                   endif
                else
                   if (ib == 1) then
                      flx_absdv(c,i) = flx_absd_snw(c,i,ib)*(1.-albsnd(c,ib))
                      flx_absiv(c,i) = flx_absi_snw(c,i,ib)*(1.-albsni(c,ib))
                   elseif (ib == 2) then
                      flx_absdn(c,i) = flx_absd_snw(c,i,ib)*(1.-albsnd(c,ib))
                      flx_absin(c,i) = flx_absi_snw(c,i,ib)*(1.-albsni(c,ib))
                   endif
                endif
             enddo
          endif
       enddo
    enddo

  end subroutine SnowAlbedo
