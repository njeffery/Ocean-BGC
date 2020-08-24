! copyright (c) 2013,  los alamos national security, llc (lans)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_benthos
!
!> \brief benthos_parms
!> \author Nicole Jeffery
!> \date   07/21/2020
!> \details
!> This module contains parameters and type definitions for the
!> ocean benthos biogeochemical model.
!
!-----------------------------------------------------------------------
MODULE benthos_parms

  !-----------------------------------------------------------------------------

  implicit none
  save
  private

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: benthos_parms_init

  !-----------------------------------------------------------------------------
  !   kinds definitions
  !-----------------------------------------------------------------------------
   integer, parameter, public ::               &
      benthos_char           = 256                    ,&
      benthos_log            = kind(.true.)           ,&
      benthos_i4             = selected_int_kind(6)   ,&
      benthos_i8             = selected_int_kind(13)  ,&
      benthos_r4             = selected_real_kind(6)  ,&
      benthos_r8             = selected_real_kind(13)

  real (KIND=benthos_r8), public :: T0_Kelvin_benthos
  !-----------------------------------------------------------------------------
  !   floating point constants used across ecosystem module
  !-----------------------------------------------------------------------------

  real (KIND=benthos_r8), parameter, public :: &
       c0_benthos  =   0.0_benthos_r8,    &
       c1_benthos  =   1.0_benthos_r8,    &
       c2_benthos  =   2.0_benthos_r8,    &
       c10_benthos = 10.0_benthos_r8,     &
       p5_benthos  =  0.5_benthos_r8,     &
       p1_benthos  =  0.1_benthos_r8

  real(KIND=benthos_r8), parameter, public :: &
       puny          = 1.0e-15_benthos_r8, & ! small number
       sec_per_day   = 86400.0_benthos_r8,  & ! number of seconds in a day
       sec_per_year  = 31536000_benthos_r8, & ! number pf seconds per year
       days_per_sec  = c1_benthos / sec_per_day,            & ! number of days in a second
       years_per_sec = c1_benthos / sec_per_year, & ! number of years in a second
       m2percm2      = 0.0001_benthos_r8, &      ! m2 per cm2
       mM_umolperL   = 1000.0_benthos_r8, &  ! mM to umol/L
       fluxSed_o     = 1.22_benthos_r8 / sec_per_year, & ! Prescribed sediment flux (kg/m2/s)
       deepFlux_h2s  = -50.0_benthos_r8 / sec_per_year, & ! upward flux for h2s (mmol/m2/s)
       pocFluxRate_o = 6000.0_benthos_r8 / sec_per_year, &   ! POC flux rate (mmol/m2/s)
       bottom_water_density = 1027.0_benthos_r8, & ! kg/m3
       sediment_density = 2650.0_benthos_r8, & ! (kg/m3)
       mass_to_vol   = 1.0e3_benthos_r8 * bottom_water_density, & ! mol/kg -> mmol/m3
!       mass_to_vol   = 1.0e6_benthos_68 * sediment_density, &
       vol_to_mass   = c1_benthos / mass_to_vol   ! mmol/m3 -> mol/kg

  !-----------------------------------------------------------------------------
  !   Fixed biogeochemical parameters
  !-----------------------------------------------------------------------------

  integer (KIND=benthos_i4), parameter, public :: &
       nPrimaryReactions   = 6,  & ! number of primary reactions
       nSecondaryReactions = 16, & ! number of secondary reactions
       nCarbonateReactions = 3     ! number of carbonate dissolution reactions

  real(KIND=benthos_r8), dimension(:,:), allocatable, public :: &
       secondarySourceStoich, & ! contains mole ratios of secondary reaction sources
       secondarySinkStoich,   & ! contains mole ratios of secondary reaction sinks
       primarySourceStoich,   & ! contains mole ratios of primary reaction sources
       primarySinkStoich,     & ! contains mole ratios of primary reaction sinks
       carbonateSourceStoich, & ! contains mole ratios of carbonate reaction sources
       carbonateSinkStoich,    & ! contains mole ratios of carbonate reaction sinks
       elementRatios                  ! converts mol tracer to mol elements

  real(KIND=benthos_r8), parameter, public :: &
       benthosDepth     = 0.3_benthos_r8,   & ! (m) depth of the benthos active layer
       benthos_salt_min = 0.1_benthos_r8,   & ! minimum salinity (ppt)
       benthos_dic_min  = benthos_salt_min/35.0_benthos_r8 * 1944.0_benthos_r8, & ! minium dic (mmol/m3)
       benthos_alk_min  = benthos_salt_min/35.0_benthos_r8 * 2225.0_benthos_r8, & ! minium alkalinity (mmolEq/m3)
       oceanMagnesium = 53.0_benthos_r8/sediment_density, & ! (mol/kg)  concentration of Mg in the ocean
       benthos_phlo_3d_init = 6.0_benthos_r8, & ! lower bound for pH
       benthos_phhi_3d_init = 9.0_benthos_r8, & ! higher bound for pH
       CtoP = 106.0_benthos_r8, & ! carbon to phosphorus ratio of sinking organic matter  (117 in BGC)
       NtoP = 16.0_benthos_r8,  & ! carbon to nitrogen ratio of organic sinking matter)
       ironBoundPFraction = 0.175 ! fraction of bound p to iron

  !-----------------------------------------------------------------------------
  !   functional group
  !-----------------------------------------------------------------------------

  type, public :: benthos_indices_type
   integer (KIND=benthos_i4) :: &
      poca_ind,          & ! particulate labile organic carbon
      pocb_ind,          & ! particulate semi-labile organic carbon
      pocc_ind,          & ! particulate refractory organic carbon
      pona_ind,          & ! particulate labile organic nitrogen
      ponb_ind,          & ! particulate semi-labile organic nitrogen
      ponc_ind,          & ! particulate refractory organic nitrogen
      popa_ind,          & ! particulate labile organic phophorus
      popb_ind,          & ! particulate semi-labile organic phophorus
      popc_ind,          & ! particulate refractory organic phophorus
      o2_ind,            & ! dissolved oxygen
      nh4_ind,           & ! dissolved ammonia
      h2po4_ind,         & ! dissolved inorganic phophate
      co2_ind,           & ! dissolved carbon dioxide
      no3_ind,           & ! dissolved nitrate
      mno2a_ind,         & ! particulate amorphous manganese oxyhydroxide
      mno2b_ind,         & ! particulate crystalline manganese oxyhydroxide
      mn_ind,            & ! dissolved manganese(II)
      feoh3a_ind,        & ! particulate amorphous iron oxyhydroxide
      feoh3b_ind,        & ! particulate crystalline iron oxyhydroxide
      fe_ind,            & ! dissolved iron(II)
      fepa_ind,          & ! particulate iron-bound phosphate amorphous
      fepb_ind,          & ! particulate iron-bound phosphate crystalline
      so4_ind,           & ! dissolved sulfate
      h2s_ind,           & ! dissolved hydrogen sulfide
      ch4_ind,           & ! dissolved methane
      hco3_ind,          & ! dissolved bicarbonate
      fes_ind,           & ! particulate iron monosulfide
      fes2_ind,          & ! particulate pyrite
      s_ind,             & ! particulate elemental sulfur
      caco3a_ind,        & ! particulate calcite
      caco3b_ind,        & ! particulate aragonite
      co3_ind,           & ! dissolved carbonate
      camgco3_ind,       & ! partiuclate 15% magnesium calcite
      alk_ind,           & ! alkalinity
      dic_ind              ! dissolved inorganic carbon

   character (benthos_char), allocatable, dimension(:) ::  &
      short_name,      & ! short name of variable
      long_name,       & ! long name of variable
      units              ! units of variable

  end type benthos_indices_type

  type, public :: benthos_element_indices_type
   integer (KIND=benthos_i4) :: &
      carbon_ind,        & ! total carbon
      oxygen_ind,        & ! total oxygen
      nitrogen_ind,      & ! total nitrogen
      phosphorus_ind,    & ! total phosphorus
      sulfur_ind,        & ! total sulfur
      manganese_ind,     & ! total manganese
      iron_ind             ! total iron

   character (benthos_char), allocatable, dimension(:) ::  &
      short_name,      & ! short name of variable
      long_name,       & ! long name of variable
      units              ! units of variable

  end type benthos_element_indices_type

  type, public :: benthos_input_type

   real (KIND=benthos_r8), allocatable, dimension(:,:,:) :: &
      tracerConc, &
      benthosTracerBulk
   real (KIND=benthos_r8), allocatable, dimension(:,:) :: &
      PH_PREV_3D, &
      deepStorage
    integer (KIND=benthos_i4), allocatable, dimension(:) :: &
      tracerType, &
      bottom_level
    real (KIND=benthos_r8), allocatable, dimension(:) :: &
      oceanBottomDepth, &
      oceanBottomTemperature, &
      oceanBottomSalinity, &
      oceanBottomSilicate   ! and any other bottom values that are needed!!!

  end type benthos_input_type

  type, public :: benthos_forcing_type   ! not sure I need this?
   real (KIND=benthos_r8), allocatable, dimension(:,:) :: &
      deepStoragePrescribedFlux,        &
      benthosOceanSedimentFlux

   real (KIND=benthos_r8), allocatable, dimension(:) :: &
      sediment_FLUX_IN,          &
      iceFraction,               &
      iceKeelDepth,              &
      fastIceFraction

   real (KIND=benthos_r8), allocatable, dimension(:,:) :: &
      depositionFlux, &
      riverFlux, &
      gasFlux, &
      seaIceFlux, &
      netFlux

  end type benthos_forcing_type

  type, public :: benthos_output_type
   real (KIND=benthos_r8), allocatable, dimension(:,:,:) :: &
      benthosTendencies, &
      benthosReactionTendencies, &
      benthosTransportTendencies, &
      elementTendencies, &
      elementReactionTendencies, &
      elementTransportTendencies

   real (KIND=benthos_r8), allocatable, dimension(:,:) :: &
      elementDeepStorage

  end type benthos_output_type

  type, public :: benthos_flux_diagnostics_type
   real (KIND=benthos_r8), allocatable, dimension(:,:) :: &
      deepStorageDiffVelFlux,           &
      deepStorageFlux,                  &
      benthosOceanDiffVelFlux,          &
      benthosOceanFlux,                 &
      elementSedimentation,             &
      elementBurialExchange,            &
      elementOceanExchange

  end type benthos_flux_diagnostics_type

  type, public :: benthos_diagnostics_type
! 2D stuff
    real (KIND=benthos_r8), allocatable, dimension(:,:) :: &
      diag_primarySourceTend,   &
      diag_primarySinkTend,     &
      diag_secondarySourceTend, &
      diag_secondarySinkTend,   &
      diag_carbonateSourceTend, &
      diag_carbonateSinkTend

! more 2D stuff
    real (KIND=benthos_r8), allocatable, dimension(:,:) :: &
      diag_netElements,             &
      diag_netBenthosTracersCell,   &
      diag_netBenthosReactionTend,  &
      diag_netElementsReactionTend,  &
      diag_netBenthosTransportTend, &
      diag_netElementsTransportTend, &
      diag_netBenthosTend,          &
      diag_netElementsTend

  end type benthos_diagnostics_type

  !----------------------------------------------------------------------------
  !   ecosystem parameters accessible via namelist input
  !----------------------------------------------------------------------------

   integer(KIND=benthos_i4), public :: &
       benthosTestCase

   logical(KIND=benthos_log), public :: &
       useCarbonateSaturation,     & ! Turns on carbonate calculations
       useSecondaryReactions,      & ! Compute secondary reactionset
       useDepthDependentPorosity,  & ! fix Porosity to 0.5 or use  data estimate of porosity distribution with depth
       useNonZeroDiffusivity,      & ! set diffusivity to zero or use biodiffusion and depth dependent molecular diffusion
       useSedimentation,           & ! Turn on or off sedimentation rates
       usePositiveSedimentation,   & ! Set sign of sedimentation rates
       useFastSedimentation,       & ! use fast sedimentation rates (for testing advection)
       useBGCSinkingFlux,          & ! Turn on and off  biochemical sediment fluxes
       useStepInitialProfiles,     & ! use Step function Initial profiles
       useBenthicReactions,        & ! Turn on or off reaction terms
       useFluxCorrection,          & ! correct flux in check_conservation_FCT.m
       useOceanConc,               & ! Turn on dissolved biogeochemical tracers in the ocean
       useDeepSource,              & ! Turn on deep source from below the active benthos layers
       useConstantDiffusivity        ! Turn off the depth dependence of diffusivty except no flux
                                     ! bottom boundary and no flux top boundary for solids
  !---------------------------------------------------------------------
  !     Vertical transport parameters
  !---------------------------------------------------------------------

   real(KIND=benthos_r8), allocatable, dimension(:), public :: &
        vertBenthosGridThickI,   & ! normalized thickness of interface vertical grid levels
        benthosMidPointI,           & ! (m) depth of mid-point of intefrace vertical Grid
        benthosGridPointsI,      & ! location of interface grid values (positive down) n+1
        benthosGridPoints,       & ! location of grid values n
        benthosPorosity,         & ! liquid fraction of benthos grid levels
        benthosPorosityI,        & ! liquid fraction at benthos interface levels
        benthosSolidPorosity,    & ! solid fraction of benthos grid levels
        benthosSolidPorosityI,   & ! solid fraction at benthos interface levels
        benthosDiffusivity,      & ! solute diffusivity at benthos grid levels
        benthosDiffusivityI,     & ! solute diffusivity at benthos interface grid levels
        benthosSolidDiffusivity, & ! solid diffusivity at benthos grid levels
        benthosSolidDiffusivityI   ! solid diffusivity at benthos interface grid levels

  real(KIND=benthos_r8), parameter, public :: &
        porosity_inf   = 0.877_benthos_r8, & ! porosity at the bottom boundary
        porosity_o     = 0.943_benthos_r8, & ! porosity at the benthos surface
        porosity_gamma = 0.05301_benthos_r8, & ! (m) Porosity e-folding distance
        fauna_biomass_max = 1000.0_benthos_r8, & ! (g/m2) saturation biomass
        fauna_biomass     = 1000.0_benthos_r8, & ! (g/m2) surface biomass use in biodiffusion
        molecular_diff    = 0.035_benthos_r8, & ! (m2/y) molecular diffusivity (Hensen et al 1998 for nitrate)
        max_bio_diff      = 5.411_benthos_r8    ! (cm2/s) maximum biodiffusivity

  !---------------------------------------------------------------------
  !     Temperature parameters
  !---------------------------------------------------------------------

!  real(KIND=benthos_r8), parameter :: &
!       Tref = 30.0_BGC_r8, & ! reference temperature (C)
!       Q_10 = 1.5_BGC_r8     ! factor for temperature dependence (non-dim)

  !*****************************************************************************

CONTAINS

  !*****************************************************************************

  subroutine benthos_parms_init

! !INPUT/OUTPUT parameterS:

    !---------------------------------------------------------------------------
    !   default namelist settings:  NONE YET
    !---------------------------------------------------------------------------

    ! turn of carbonate saturation subroutines
    useCarbonateSaturation = .true.

    ! turn of carbonate saturation subroutines
    useSecondaryReactions = .true.

    ! fix Porosity to 0.5 or use  data estimate of porosity distribution with depth
    useDepthDependentPorosity = .true.

    ! set diffusivity to zero or use biodiffusion and depth dependent molecular diffusion
    useNonZeroDiffusivity = .true.

    ! Turn on or off sedimentation rates
    useSedimentation = .true.

    ! Set sign of sedimentation rates
    usePositiveSedimentation = .true.

    ! use fast sedimentation rates (for testing advection)
    useFastSedimentation = .false.

    ! Turn on and off  biochemical sediment fluxes
    useBGCSinkingFlux = .true.

    ! use Step function Initial profiles
    useStepInitialProfiles = .false.

    ! Turn on or off reaction terms
    useBenthicReactions = .true.

    ! correct flux in check_conservation_FCT.m
    useFluxCorrection = .true.

    ! Turn on dissolved biogeochemical tracers in the ocean
    useOceanConc = .true.

    ! Turn on deep source from below the active benthos layers
    useDeepSource = .true.

    ! Turn off the depth dependence of diffusivty except no flux bottom boundary and no flux top boundary for solids
    useConstantDiffusivity = .false.

    ! Set to default configuration
    benthosTestCase = 0

  END SUBROUTINE benthos_parms_init

  !*****************************************************************************

END MODULE benthos_parms
