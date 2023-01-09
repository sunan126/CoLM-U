
MODULE AerosolMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  !
  IMPLICIT NONE
  SAVE
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: AerosolMasses
  public :: AerosolFluxes
  !
  ! !PUBLIC DATA MEMBERS:
  !-----------------------------------------------------------------------
    integer, parameter :: r8 = selected_real_kind(12) ! 8 byte real

    integer, parameter :: begc = 1      !  beginning column index
    integer, parameter :: endc = 1      !  beginning and ending column index
    real(r8) :: dtime = 1800.0_r8       !  land model time step (sec)

    logical, parameter :: use_extrasnowlayers = .false.
    integer, parameter :: nlevsno = 5   !  maximum number of snow layers

    real(r8), public, parameter :: snw_rds_min = 54.526_r8  ! minimum allowed snow effective radius (also "fresh snow" value) [microns

contains

  !-----------------------------------------------------------------------
  subroutine AerosolMasses( snl           ,do_capsnow     ,&
             h2osno_ice    ,h2osno_liq    ,qflx_snwcp_ice ,snw_rds       ,&

             mss_bcpho     ,mss_bcphi     ,mss_ocpho      ,mss_ocphi     ,&
             mss_dst1      ,mss_dst2      ,mss_dst3       ,mss_dst4      ,&

             mss_cnc_bcphi ,mss_cnc_bcpho ,mss_cnc_ocphi  ,mss_cnc_ocpho ,& 
             mss_cnc_dst1  ,mss_cnc_dst2  ,mss_cnc_dst3   ,mss_cnc_dst4  ) 

    !
    ! !DESCRIPTION:
    ! Calculate column-integrated aerosol masses, and
    ! mass concentrations for radiative calculations and output
    ! (based on new snow level state, after SnowFilter is rebuilt.
    ! NEEDS TO BE AFTER SnowFiler is rebuilt in Hydrology2, otherwise there
    ! can be zero snow layers but an active column in filter)

    IMPLICIT NONE

    ! !ARGUMENTS:
    !
    integer, intent(in)     ::   snl ( begc:endc )            !  number of snow layers

    logical,  intent(in)    ::  do_capsnow    ( begc:endc )   !  true => do snow capping
    real(r8), intent(in)    ::  h2osno_ice    ( begc:endc, -nlevsno+1:0 ) !  ice lens (kg/m2)
    real(r8), intent(in)    ::  h2osno_liq    ( begc:endc, -nlevsno+1:0 ) !  liquid water (kg/m2)
    real(r8), intent(in)    ::  qflx_snwcp_ice( begc:endc )   !  excess snowfall due to snow capping (mm H2O /s) [+]

    real(r8), intent(inout) ::  snw_rds       ( begc:endc, -nlevsno+1:0 ) !  effective snow grain radius (col,lyr) [microns, m^-6]

    real(r8), intent(inout) ::  mss_bcpho     ( begc:endc, -nlevsno+1:0 ) !  mass of hydrophobic BC in snow (col,lyr) [kg]
    real(r8), intent(inout) ::  mss_bcphi     ( begc:endc, -nlevsno+1:0 ) !  mass of hydrophillic BC in snow (col,lyr) [kg]
    real(r8), intent(inout) ::  mss_ocpho     ( begc:endc, -nlevsno+1:0 ) !  mass of hydrophobic OC in snow (col,lyr) [kg]
    real(r8), intent(inout) ::  mss_ocphi     ( begc:endc, -nlevsno+1:0 ) !  mass of hydrophillic OC in snow (col,lyr) [kg]
    real(r8), intent(inout) ::  mss_dst1      ( begc:endc, -nlevsno+1:0 ) !  mass of dust species 1 in snow (col,lyr) [kg]
    real(r8), intent(inout) ::  mss_dst2      ( begc:endc, -nlevsno+1:0 ) !  mass of dust species 2 in snow (col,lyr) [kg]
    real(r8), intent(inout) ::  mss_dst3      ( begc:endc, -nlevsno+1:0 ) !  mass of dust species 3 in snow (col,lyr) [kg]
    real(r8), intent(inout) ::  mss_dst4      ( begc:endc, -nlevsno+1:0 ) !  mass of dust species 4 in snow (col,lyr) [kg]

    real(r8), intent(out)   ::  mss_cnc_bcphi ( begc:endc, -nlevsno+1:0 ) !  mass concentration of BC species 1 (col,lyr) [kg/kg]
    real(r8), intent(out)   ::  mss_cnc_bcpho ( begc:endc, -nlevsno+1:0 ) !  mass concentration of BC species 2 (col,lyr) [kg/kg]
    real(r8), intent(out)   ::  mss_cnc_ocphi ( begc:endc, -nlevsno+1:0 ) !  mass concentration of OC species 1 (col,lyr) [kg/kg]
    real(r8), intent(out)   ::  mss_cnc_ocpho ( begc:endc, -nlevsno+1:0 ) !  mass concentration of OC species 2 (col,lyr) [kg/kg]
    real(r8), intent(out)   ::  mss_cnc_dst1  ( begc:endc, -nlevsno+1:0 ) !  mass concentration of dust species 1 (col,lyr) [kg/kg]
    real(r8), intent(out)   ::  mss_cnc_dst2  ( begc:endc, -nlevsno+1:0 ) !  mass concentration of dust species 2 (col,lyr) [kg/kg]
    real(r8), intent(out)   ::  mss_cnc_dst3  ( begc:endc, -nlevsno+1:0 ) !  mass concentration of dust species 3 (col,lyr) [kg/kg]
    real(r8), intent(out)   ::  mss_cnc_dst4  ( begc:endc, -nlevsno+1:0 ) !  mass concentration of dust species 4 (col,lyr) [kg/kg]

    ! !LOCAL VARIABLES:
    integer  :: c,j             ! indices
    real(r8) :: snowmass        ! liquid+ice snow mass in a layer [kg/m2]
    real(r8) :: snowcap_scl_fct ! temporary factor used to correct for snow capping

    !-----------------------------------------------------------------------

      do c = begc, endc
         do j = -nlevsno+1, 0

            ! layer mass of snow:
            snowmass = h2osno_ice(c,j) + h2osno_liq(c,j)

            if (.not. use_extrasnowlayers) then
               ! Correct the top layer aerosol mass to account for snow capping. 
               ! This approach conserves the aerosol mass concentration
               ! (but not the aerosol amss) when snow-capping is invoked

               if (j == snl(c)+1) then
                  if (do_capsnow(c)) then 

                     snowcap_scl_fct = snowmass / (snowmass + (qflx_snwcp_ice(c)*dtime))

                     mss_bcpho(c,j) = mss_bcpho(c,j)*snowcap_scl_fct
                     mss_bcphi(c,j) = mss_bcphi(c,j)*snowcap_scl_fct
                     mss_ocpho(c,j) = mss_ocpho(c,j)*snowcap_scl_fct
                     mss_ocphi(c,j) = mss_ocphi(c,j)*snowcap_scl_fct

                     mss_dst1(c,j)  = mss_dst1(c,j)*snowcap_scl_fct
                     mss_dst2(c,j)  = mss_dst2(c,j)*snowcap_scl_fct
                     mss_dst3(c,j)  = mss_dst3(c,j)*snowcap_scl_fct
                     mss_dst4(c,j)  = mss_dst4(c,j)*snowcap_scl_fct
                  endif
               endif
            endif

            if (j >= snl(c)+1) then

               mss_cnc_bcphi(c,j) = mss_bcphi(c,j) / snowmass
               mss_cnc_bcpho(c,j) = mss_bcpho(c,j) / snowmass

               mss_cnc_ocphi(c,j) = mss_ocphi(c,j) / snowmass
               mss_cnc_ocpho(c,j) = mss_ocpho(c,j) / snowmass

               mss_cnc_dst1(c,j)  = mss_dst1(c,j)  / snowmass
               mss_cnc_dst2(c,j)  = mss_dst2(c,j)  / snowmass
               mss_cnc_dst3(c,j)  = mss_dst3(c,j)  / snowmass
               mss_cnc_dst4(c,j)  = mss_dst4(c,j)  / snowmass

            else
               !set variables of empty snow layers to zero
               snw_rds(c,j)       = 0._r8

               mss_bcpho(c,j)     = 0._r8
               mss_bcphi(c,j)     = 0._r8
               mss_cnc_bcphi(c,j) = 0._r8
               mss_cnc_bcpho(c,j) = 0._r8

               mss_ocpho(c,j)     = 0._r8
               mss_ocphi(c,j)     = 0._r8
               mss_cnc_ocphi(c,j) = 0._r8
               mss_cnc_ocpho(c,j) = 0._r8

               mss_dst1(c,j)      = 0._r8
               mss_dst2(c,j)      = 0._r8
               mss_dst3(c,j)      = 0._r8
               mss_dst4(c,j)      = 0._r8
               mss_cnc_dst1(c,j)  = 0._r8
               mss_cnc_dst2(c,j)  = 0._r8
               mss_cnc_dst3(c,j)  = 0._r8
               mss_cnc_dst4(c,j)  = 0._r8
            endif
         enddo

      enddo

  end subroutine AerosolMasses



  !-----------------------------------------------------------------------
  subroutine AerosolFluxes( snl, forc_aer, &
                            mss_bcphi  ,mss_bcpho  ,mss_ocphi  ,mss_ocpho ,& 
                            mss_dst1   ,mss_dst2   ,mss_dst3   ,mss_dst4   ) 
    !
    ! !DESCRIPTION:
    ! Compute aerosol fluxes through snowpack and aerosol deposition fluxes into top layere
    !
    IMPLICIT NONE
    !
    !-----------------------------------------------------------------------
    ! !ARGUMENTS:
    integer, intent(in) :: snl ( begc:endc )   ! number of snow layers

    real(r8), intent(in) :: forc_aer ( begc:endc,14 )  ! aerosol deposition from atmosphere model (grd,aer) [kg m-1 s-1]

    real(r8), intent(inout) :: mss_bcphi  ( begc:endc,-nlevsno+1:0 )  ! hydrophillic BC mass in snow (col,lyr) [kg]
    real(r8), intent(inout) :: mss_bcpho  ( begc:endc,-nlevsno+1:0 )  ! hydrophobic  BC mass in snow (col,lyr) [kg]
    real(r8), intent(inout) :: mss_ocphi  ( begc:endc,-nlevsno+1:0 )  ! hydrophillic OC mass in snow (col,lyr) [kg]
    real(r8), intent(inout) :: mss_ocpho  ( begc:endc,-nlevsno+1:0 )  ! hydrophobic  OC mass in snow (col,lyr) [kg]
    real(r8), intent(inout) :: mss_dst1   ( begc:endc,-nlevsno+1:0 )  ! mass of dust species 1 in snow (col,lyr) [kg]
    real(r8), intent(inout) :: mss_dst2   ( begc:endc,-nlevsno+1:0 )  ! mass of dust species 2 in snow (col,lyr) [kg]
    real(r8), intent(inout) :: mss_dst3   ( begc:endc,-nlevsno+1:0 )  ! mass of dust species 3 in snow (col,lyr) [kg]
    real(r8), intent(inout) :: mss_dst4   ( begc:endc,-nlevsno+1:0 )  ! mass of dust species 4 in snow (col,lyr) [kg]

    ! !LOCAL VARIABLES:
    real(r8) :: flx_bc_dep       ( begc:endc )   ! total BC deposition (col) [kg m-2 s-1]
    real(r8) :: flx_bc_dep_phi   ( begc:endc )   ! hydrophillic BC deposition (col) [kg m-1 s-1]
    real(r8) :: flx_bc_dep_pho   ( begc:endc )   ! hydrophobic BC deposition (col) [kg m-1 s-1]
    real(r8) :: flx_oc_dep       ( begc:endc )   ! total OC deposition (col) [kg m-2 s-1]
    real(r8) :: flx_oc_dep_phi   ( begc:endc )   ! hydrophillic OC deposition (col) [kg m-1 s-1]
    real(r8) :: flx_oc_dep_pho   ( begc:endc )   ! hydrophobic OC deposition (col) [kg m-1 s-1]
    real(r8) :: flx_dst_dep      ( begc:endc )   ! total dust deposition (col) [kg m-2 s-1]

    real(r8) :: flx_dst_dep_wet1 ( begc:endc )   ! wet dust (species 1) deposition (col) [kg m-2 s-1]
    real(r8) :: flx_dst_dep_dry1 ( begc:endc )   ! dry dust (species 1) deposition (col) [kg m-2 s-1]
    real(r8) :: flx_dst_dep_wet2 ( begc:endc )   ! wet dust (species 2) deposition (col) [kg m-2 s-1]
    real(r8) :: flx_dst_dep_dry2 ( begc:endc )   ! dry dust (species 2) deposition (col) [kg m-2 s-1]
    real(r8) :: flx_dst_dep_wet3 ( begc:endc )   ! wet dust (species 3) deposition (col) [kg m-2 s-1]
    real(r8) :: flx_dst_dep_dry3 ( begc:endc )   ! dry dust (species 3) deposition (col) [kg m-2 s-1]
    real(r8) :: flx_dst_dep_wet4 ( begc:endc )   ! wet dust (species 4) deposition (col) [kg m-2 s-1]
    real(r8) :: flx_dst_dep_dry4 ( begc:endc )   ! dry dust (species 4) deposition (col) [kg m-2 s-1]

    integer  :: c

    !-----------------------------------------------------------------------
    ! set aerosol deposition fluxes from forcing array
    ! The forcing array is either set from an external file
    ! or from fluxes received from the atmosphere model
#ifdef MODAL_AER
    ! Mapping for modal aerosol scheme where within-hydrometeor and
    ! interstitial aerosol fluxes are differentiated. Here, "phi"
    ! flavors of BC and OC correspond to within-hydrometeor
    ! (cloud-borne) aerosol, and "pho" flavors are interstitial
    ! aerosol. "wet" and "dry" fluxes of BC and OC specified here are
    ! purely diagnostic
    do c = begc, endc

       flx_bc_dep_phi(c)   = forc_aer(c,3)
       flx_bc_dep_pho(c)   = forc_aer(c,1) + forc_aer(c,2)
       flx_bc_dep(c)       = forc_aer(c,1) + forc_aer(c,2) + forc_aer(c,3)

       flx_oc_dep_phi(c)   = forc_aer(c,6)
       flx_oc_dep_pho(c)   = forc_aer(c,4) + forc_aer(c,5)
       flx_oc_dep(c)       = forc_aer(c,4) + forc_aer(c,5) + forc_aer(c,6)

       flx_dst_dep_wet1(c) = forc_aer(c,7)
       flx_dst_dep_dry1(c) = forc_aer(c,8)
       flx_dst_dep_wet2(c) = forc_aer(c,9)
       flx_dst_dep_dry2(c) = forc_aer(c,10)
       flx_dst_dep_wet3(c) = forc_aer(c,11)
       flx_dst_dep_dry3(c) = forc_aer(c,12)
       flx_dst_dep_wet4(c) = forc_aer(c,13)
       flx_dst_dep_dry4(c) = forc_aer(c,14)
       flx_dst_dep(c)      = forc_aer(c,7) + forc_aer(c,8) + forc_aer(c,9) + &
                             forc_aer(c,10) + forc_aer(c,11) + forc_aer(c,12) + &
                             forc_aer(c,13) + forc_aer(c,14)
    end do

#else

    ! Original mapping for bulk aerosol deposition. phi and pho BC/OC
    ! species are distinguished in model, other fluxes (e.g., dry and
    ! wet BC/OC) are purely diagnostic.

      do c = begc,endc

         flx_bc_dep_phi(c)   = forc_aer(c,1) + forc_aer(c,3)
         flx_bc_dep_pho(c)   = forc_aer(c,2)
         flx_bc_dep(c)       = forc_aer(c,1) + forc_aer(c,2) + forc_aer(c,3)

         flx_oc_dep_phi(c)   = forc_aer(c,4) + forc_aer(c,6)
         flx_oc_dep_pho(c)   = forc_aer(c,5)
         flx_oc_dep(c)       = forc_aer(c,4) + forc_aer(c,5) + forc_aer(c,6)

         flx_dst_dep_wet1(c) = forc_aer(c,7)
         flx_dst_dep_dry1(c) = forc_aer(c,8)
         flx_dst_dep_wet2(c) = forc_aer(c,9)
         flx_dst_dep_dry2(c) = forc_aer(c,10)
         flx_dst_dep_wet3(c) = forc_aer(c,11)
         flx_dst_dep_dry3(c) = forc_aer(c,12)
         flx_dst_dep_wet4(c) = forc_aer(c,13)
         flx_dst_dep_dry4(c) = forc_aer(c,14)
         flx_dst_dep(c)      = forc_aer(c,7) + forc_aer(c,8) + forc_aer(c,9) + &
                               forc_aer(c,10) + forc_aer(c,11) + forc_aer(c,12) + &
                               forc_aer(c,13) + forc_aer(c,14)
      end do
#endif

      ! aerosol deposition fluxes into top layer
      ! This is done after the inter-layer fluxes so that some aerosol
      ! is in the top layer after deposition, and is not immediately
      ! washed out before radiative calculations are done


      do c = begc,endc
         mss_bcphi(c,snl(c)+1) = mss_bcphi(c,snl(c)+1) + (flx_bc_dep_phi(c)*dtime)
         mss_bcpho(c,snl(c)+1) = mss_bcpho(c,snl(c)+1) + (flx_bc_dep_pho(c)*dtime)
         mss_ocphi(c,snl(c)+1) = mss_ocphi(c,snl(c)+1) + (flx_oc_dep_phi(c)*dtime)
         mss_ocpho(c,snl(c)+1) = mss_ocpho(c,snl(c)+1) + (flx_oc_dep_pho(c)*dtime)

         mss_dst1(c,snl(c)+1) = mss_dst1(c,snl(c)+1) + (flx_dst_dep_dry1(c) + flx_dst_dep_wet1(c))*dtime
         mss_dst2(c,snl(c)+1) = mss_dst2(c,snl(c)+1) + (flx_dst_dep_dry2(c) + flx_dst_dep_wet2(c))*dtime
         mss_dst3(c,snl(c)+1) = mss_dst3(c,snl(c)+1) + (flx_dst_dep_dry3(c) + flx_dst_dep_wet3(c))*dtime
         mss_dst4(c,snl(c)+1) = mss_dst4(c,snl(c)+1) + (flx_dst_dep_dry4(c) + flx_dst_dep_wet4(c))*dtime
      end do

  end subroutine AerosolFluxes


END MODULE AerosolMod

