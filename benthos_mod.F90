!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module benthos_mod

!BOP
!  MODULE: benthos_mod
!
!> \brief MPAS ocean benthos
!> \author Nicole Jeffery
!> \date   07/15/2020
!> \details
!------------------------------------------------------------------------------
!> This subroutine computes the biogeochemistry of a ~30 cm sea-floor sediment
!> layer resolved at the cm scale.  The 35 component coupled-biogeochemistry
!> is based on Reed et al 2011 and Krumins et al 2013.  The vertical transport scheme
!> is similar to Jeffery et al. 2011 for sea ice and modified for a  fixed depth
!> porous sediment layer.
!------------------------------------------------------------------------------

!  REVISION HISTORY:

!  SVN:$Id:  $

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!  USES:


!  INPUT PARAMETERS:
!-----------------------------------------------------------------------------
!   include benthos parameters
!   all variables from these modules have a parm_ prefix
 !-----------------------------------------------------------------------------

   use benthos_parms

   implicit none
   save
   private

!-----------------------------------------------------------------------
!  public/private declarations
!-----------------------------------------------------------------------

   public :: &
      benthos_init,           &
      benthos_flags_init,     &
      benthos_SurfaceFluxes,  &
      benthos_SourceSink

!-----------------------------------------------------------------------
!  module variables
!-----------------------------------------------------------------------

   integer (KIND=benthos_i4), public, parameter :: &
      benthos_tracer_cnt = 35, &
      element_tracer_cnt = 7

   integer (KIND=benthos_i4), public :: &
      nBenthicTracers, &
      nElements

   integer (KIND=benthos_i4), public :: &
      poca_ind, &
      pocb_ind, &
      pocc_ind, &
      pona_ind, &
      ponb_ind, &
      ponc_ind, &
      popa_ind, &
      popb_ind, &
      popc_ind, &
      o2_ind, &
      nh4_ind, &
      h2po4_ind, &
      co2_ind, &
      no3_ind, &
      mno2a_ind, &
      mno2b_ind, &
      mn_ind, &
      feoh3a_ind, &
      feoh3b_ind, &
      fe_ind, &
      fepa_ind, &
      fepb_ind, &
      so4_ind, &
      h2s_ind, &
      ch4_ind, &
      dic_ind, &
      alk_ind, &
      hco3_ind, &
      fes_ind, &
      fes2_ind, &
      s_ind, &
      caco3a_ind, &
      caco3b_ind, &
      co3_ind, &
      camgco3_ind

!*****************************************************************************

contains

!*****************************************************************************
!BOP
! !IROUTINE: benthos_init
! !INTERFACE:

 subroutine benthos_init(benthos_indices, benthos_element_indices, &
    nBenthicVertLevels, nBenthicTracers_in, nElements_in, benthosInterfaceLayerThickness)

! !DESCRIPTION:
!  Initialize tracegas tracer module. This involves setting metadata, reading
!  the module namelist, setting initial conditions, setting up forcing,
!  and defining additional tavg variables.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer(KIND=benthos_i4), intent(in) :: &
     nBenthicVertLevels, &
     nBenthicTracers_in, &
     nElements_in

! !INPUT/OUTPUT PARAMETERS:

  type(benthos_indices_type), intent(inout) :: benthos_indices
  type(benthos_element_indices_type), intent(inout) :: benthos_element_indices

! !OUTPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(:), intent(out) :: &
     benthosInterfaceLayerThickness ! (m)

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

  !NJ-TEST
  logical(KIND=benthos_log) :: test_flag
  integer(KIND=benthos_i4) :: iBenthicTracers, iSecondaryReactions, iElements, iBenthicVertLevels
  !NJ-END
!-----------------------------------------------------------------------
!  Local copy of indices
  !-----------------------------------------------------------------------

   poca_ind    = benthos_indices%poca_ind
   pocb_ind    = benthos_indices%pocb_ind
   pocc_ind    = benthos_indices%pocc_ind
   pona_ind    = benthos_indices%pona_ind
   ponb_ind    = benthos_indices%ponb_ind
   ponc_ind    = benthos_indices%ponc_ind
   popa_ind    = benthos_indices%popa_ind
   popb_ind    = benthos_indices%popb_ind
   popc_ind    = benthos_indices%popc_ind
   o2_ind      = benthos_indices%o2_ind
   nh4_ind     = benthos_indices%nh4_ind
   h2po4_ind   = benthos_indices%h2po4_ind
   co2_ind     = benthos_indices%co2_ind
   no3_ind     = benthos_indices%no3_ind
   mno2a_ind   = benthos_indices%mno2a_ind
   mno2b_ind   = benthos_indices%mno2b_ind
   mn_ind      = benthos_indices%mn_ind
   feoh3a_ind  = benthos_indices%feoh3a_ind
   feoh3b_ind  = benthos_indices%feoh3b_ind
   fe_ind      = benthos_indices%fe_ind
   fepa_ind    = benthos_indices%fepa_ind
   fepb_ind    = benthos_indices%fepb_ind
   so4_ind     = benthos_indices%so4_ind
   h2s_ind     = benthos_indices%h2s_ind
   ch4_ind     = benthos_indices%ch4_ind
   dic_ind     = benthos_indices%dic_ind
   alk_ind     = benthos_indices%alk_ind
   hco3_ind    = benthos_indices%hco3_ind
   fes_ind     = benthos_indices%fes_ind
   fes2_ind    = benthos_indices%fes2_ind
   s_ind       = benthos_indices%s_ind
   caco3a_ind  = benthos_indices%caco3a_ind
   caco3b_ind  = benthos_indices%caco3b_ind
   co3_ind     = benthos_indices%co3_ind
   camgco3_ind = benthos_indices%camgco3_ind

!-------------------------------------------------
! initialize tracer and element numbers
!-------------------------------------------------

   nElements = nElements_in
   nBenthicTracers = nBenthicTracers_in

   if (nBenthicTracers.ne.benthos_tracer_cnt) then
     !  write(*,*) 'error: nBenthicTracers is inconsistent with assumed value'
     ! err = 1
     ! write(*,*) 'nBenthicTracers:', nBenthicTracers
     ! write(*,*) 'benthos_tracer_cnt:', benthos_tracer_cnt
   end if
!-----------------------------------------------------------------------
!  Define benthos tracers
!-----------------------------------------------------------------------
   !NJ-TEST
   !write(*,*) 'start benthos_init'
   !NJ-END
   
   benthos_indices%units(:) = 'mmol/m3'
   benthos_element_indices%units(:) = 'mmol/m3'

   benthos_indices%short_name(poca_ind)='POCa'
   benthos_indices%long_name(poca_ind)='Particulate labile organic carbon'
   benthos_indices%units(poca_ind)='mmol/kg'

   benthos_indices%short_name(pocb_ind)='POCb'
   benthos_indices%long_name(pocb_ind)='Particulate semi-labile organic carbon'
   benthos_indices%units(pocb_ind)='mmol/kg'
   benthos_indices%short_name(pocc_ind)='POCc'
   benthos_indices%long_name(pocc_ind)='Particulate refractory organic carbon'
   benthos_indices%units(pocc_ind)='mmol/kg'

   benthos_indices%short_name(pona_ind)='PONa'
   benthos_indices%long_name(pona_ind)='Particulate labile organic nitrogen'
   benthos_indices%units(pona_ind)='mmol/kg'

   benthos_indices%short_name(ponb_ind)='PONb'
   benthos_indices%long_name(ponb_ind)='Particulate semi-labile organic nitrogen'
   benthos_indices%units(ponb_ind)='mmol/kg'

   benthos_indices%short_name(ponc_ind)='PONc'
   benthos_indices%long_name(ponc_ind)='Particulate refractory organic nitrogen'
   benthos_indices%units(ponc_ind)='mmol/kg'

   benthos_indices%short_name(popa_ind)='POPa'
   benthos_indices%long_name(popa_ind)='Particulate labile organic phosphorus'
   benthos_indices%units(popa_ind)='mmol/kg'

   benthos_indices%short_name(popb_ind)='POPb'
   benthos_indices%long_name(popb_ind)='Particulate semi-labile organic phosphorus'
   benthos_indices%units(popb_ind)='mmol/kg'

   benthos_indices%short_name(popc_ind)='POPc'
   benthos_indices%long_name(popc_ind)='Particulate refractory organic phosphorus'
   benthos_indices%units(popc_ind)='mmol/kg'

   benthos_indices%short_name(o2_ind)='O2'
   benthos_indices%long_name(o2_ind)='Dissolved oxygen'

   benthos_indices%short_name(nh4_ind)='NH4'
   benthos_indices%long_name(nh4_ind)='Dissolved ammonium'

   benthos_indices%short_name(h2po4_ind)='H2PO4'
   benthos_indices%long_name(h2po4_ind)='Dissolved phosphate'

   benthos_indices%short_name(co2_ind)='CO2'
   benthos_indices%long_name(co2_ind)='Dissolved carbon dioxide'

   benthos_indices%short_name(no3_ind)='NO3'
   benthos_indices%long_name(no3_ind)='Dissolved nitrate'

   benthos_indices%short_name(mno2a_ind)='MnO2a'
   benthos_indices%long_name(mno2a_ind)='Particulate amorphous manganese oxyhydroxide'
   benthos_indices%units(mno2a_ind)='mmol/kg'

   benthos_indices%short_name(mno2b_ind)='MnO2b'
   benthos_indices%long_name(mno2b_ind)='Particulate crystalline manganese oxyhydroxide'
   benthos_indices%units(mno2b_ind)='mmol/kg'

   benthos_indices%short_name(mn_ind)='Mn'
   benthos_indices%long_name(mn_ind)='Dissolved manganese'

   benthos_indices%short_name(feoh3a_ind)='Fe(OH)3a'
   benthos_indices%long_name(feoh3a_ind)='Particulate amorphous iron oxyhydroxide'
   benthos_indices%units(feoh3a_ind)='mmol/kg'

   benthos_indices%short_name(feoh3b_ind)='Fe(OH)3b'
   benthos_indices%long_name(feoh3b_ind)='Particulate crystalline iron oxyhydroxide'
   benthos_indices%units(feoh3b_ind)='mmol/kg'

   benthos_indices%short_name(fe_ind)='Fe(II)'
   benthos_indices%long_name(fe_ind)='Dissolved iron'

   benthos_indices%short_name(fepa_ind)='[Fe-P]a'
   benthos_indices%long_name(fepa_ind)='particulate amorphous iron bound phosphate'
   benthos_indices%units(fepa_ind)='mmol/kg'

   benthos_indices%short_name(fepb_ind)='[Fe-P]b'
   benthos_indices%long_name(fepb_ind)='particulate crystalline iron bound phosphate'
   benthos_indices%units(fepb_ind)='mmol/kg'

   benthos_indices%short_name(so4_ind)='SO4'
   benthos_indices%long_name(so4_ind)='Dissolved sulfate'

   benthos_indices%short_name(h2s_ind)='H2S'
   benthos_indices%long_name(h2s_ind)='Dissolved hydrogen sulfide'

   benthos_indices%short_name(ch4_ind)='CH4'
   benthos_indices%long_name(ch4_ind)='Dissolved methane'

   benthos_indices%short_name(dic_ind)='DIC'
   benthos_indices%long_name(dic_ind)='Dissolved inorganic carbon'

   benthos_indices%short_name(alk_ind)='ALK'
   benthos_indices%long_name(alk_ind)='Total alkalinity'
   benthos_indices%units(fepa_ind)='meq/m3'

   benthos_indices%short_name(hco3_ind)='HCO3'
   benthos_indices%long_name(hco3_ind)='Dissolved bicarbonate'

   benthos_indices%short_name(fes_ind)='FeS'
   benthos_indices%long_name(fes_ind)='Partculate iron monosulfide'
   benthos_indices%units(fes_ind)='mmol/kg'

   benthos_indices%short_name(fes2_ind)='FeS2'
   benthos_indices%long_name(fes2_ind)='Particulate pyrite'
   benthos_indices%units(fes2_ind)='mmol/kg'

   benthos_indices%short_name(s_ind)='S'
   benthos_indices%long_name(s_ind)='Particulate elemental sulfur'
   benthos_indices%units(s_ind)='mmol/kg'

   benthos_indices%short_name(caco3a_ind)='CaCO3a'
   benthos_indices%long_name(caco3a_ind)='Particulate calcite'
   benthos_indices%units(caco3a_ind)='mmol/kg'

   benthos_indices%short_name(caco3b_ind)='CaCO3b'
   benthos_indices%long_name(caco3b_ind)='Particulate aragonite'
   benthos_indices%units(caco3b_ind)='mmol/kg'

   benthos_indices%short_name(co3_ind)='CO3'
   benthos_indices%long_name(co3_ind)='Dissolved carbonate'

   benthos_indices%short_name(camgco3_ind)='Ca.85Mg.15CO3'
   benthos_indices%long_name(camgco3_ind)='Particulate 15% magnesium calcite'
   benthos_indices%units(camgco3_ind)='mmol/kg'

   benthos_element_indices%short_name(benthos_element_indices%carbon_ind)='Carbon'
   benthos_element_indices%long_name(benthos_element_indices%carbon_ind)='Total carbon'

   benthos_element_indices%short_name(benthos_element_indices%oxygen_ind)='Oxygen'
   benthos_element_indices%long_name(benthos_element_indices%oxygen_ind)='Total oxygen'

   benthos_element_indices%short_name(benthos_element_indices%nitrogen_ind)='Nitrogen'
   benthos_element_indices%long_name(benthos_element_indices%nitrogen_ind)='Total nitrogen'

   benthos_element_indices%short_name(benthos_element_indices%phosphorus_ind)='Phosphorus'
   benthos_element_indices%long_name(benthos_element_indices%phosphorus_ind)='Total phosphorus'

   benthos_element_indices%short_name(benthos_element_indices%sulfur_ind)='Sulfur'
   benthos_element_indices%long_name(benthos_element_indices%sulfur_ind)='Total sulfur'

   benthos_element_indices%short_name(benthos_element_indices%manganese_ind)='Manganese'
   benthos_element_indices%long_name(benthos_element_indices%manganese_ind)='Total manganese'

   benthos_element_indices%short_name(benthos_element_indices%iron_ind)='Iron'
   benthos_element_indices%long_name(benthos_element_indices%iron_ind)='Total iron'

   ! initialize transport parameters

   !NJ-TEST
   !write(*,*) 'before benthos_parameters_init'
   !NJ-END
   call benthos_parameters_init (nBenthicVertLevels)
   
   do iBenthicVertLevels = 1, nBenthicVertLevels+1
      benthosInterfaceLayerThickness(iBenthicVertLevels) =  vertBenthosGridThickI(iBenthicVertLevels)*benthosDepth
   end do
   
   !NJ-TEST
   !write(*,*) 'after benthos_parameters_init'
   !NJ-END

   !NJ-TEST: VERIFIED
   test_flag = .false.
   if (test_flag) then
      write(*,*) 'diagnostics: benthos_parameters_init'
      do iBenthicVertLevels = 1, nBenthicVertLevels+2
         if (iBenthicVertLevels .lt. nBenthicVertLevels+2) then
            write(*,*) ' Interface variables '
            write(*,*) 'vertBenthosGridThickI(iBenthicVertLevels):',vertBenthosGridThickI(iBenthicVertLevels)
            write(*,*) 'benthosMidPointI(iBenthicVertLevels):',benthosMidPointI(iBenthicVertLevels)
            write(*,*) 'benthosGridPointsI(iBenthicVertLevels):',benthosGridPointsI(iBenthicVertLevels)
            write(*,*) 'benthosPorosityI(iBenthicVertLevels):',benthosPorosityI(iBenthicVertLevels)
            write(*,*) 'benthosSolidPorosityI(iBenthicVertLevels):',benthosSolidPorosityI(iBenthicVertLevels)
            write(*,*) 'benthosDiffusivityI(iBenthicVertLevels):',benthosDiffusivityI(iBenthicVertLevels)
            write(*,*) 'benthosSolidDiffusivityI(iBenthicVertLevels):',benthosSolidDiffusivityI(iBenthicVertLevels)

            write(*,*) ' Grid point variables '
            write(*,*) 'benthosGridPoints(iBenthicVertLevels):',benthosGridPoints(iBenthicVertLevels)
            write(*,*) 'benthosPorosity(iBenthicVertLevels):',benthosPorosity(iBenthicVertLevels)
            write(*,*) 'benthosSolidPorosity(iBenthicVertLevels):',benthosSolidPorosity(iBenthicVertLevels)
            write(*,*) 'benthosDiffusivity(iBenthicVertLevels):',benthosDiffusivity(iBenthicVertLevels)
            write(*,*) 'benthosSolidDiffusivity(iBenthicVertLevels):',benthosSolidDiffusivity(iBenthicVertLevels)
         else
            write(*,*) 'benthosGridPoints(iBenthicVertLevels):',benthosGridPoints(iBenthicVertLevels)
            write(*,*) 'benthosPorosity(iBenthicVertLevels):',benthosPorosity(iBenthicVertLevels)
            write(*,*) 'benthosSolidPorosity(iBenthicVertLevels):',benthosSolidPorosity(iBenthicVertLevels)
            write(*,*) 'benthosDiffusivity(iBenthicVertLevels):',benthosDiffusivity(iBenthicVertLevels)
            write(*,*) 'benthosSolidDiffusivity(iBenthicVertLevels):',benthosSolidDiffusivity(iBenthicVertLevels)
         end if
      end do

      do iElements = 1,nElements
         do iBenthicTracers = 1,nBenthicTracers
            write(*,*) 'elementRatio:', elementRatios(iBenthicTracers, iElements)
         end do
      end do
      
   end if
   !NJ-END
   
      
   ! define stoichiometry of secondary reactions

   call secondaryStoichMatrix

   !NJ-TEST : VERIFIED
   test_flag = .false.
   if (test_flag) then
      write(*,*) 'diagnostics: secondaryStoichMatrix'
      do iBenthicTracers = 1,nBenthicTracers
         do iSecondaryReactions = 1,nSecondaryReactions
            write(*,*) 'iBenthicTracers, iSecondaryReactions):',iBenthicTracers,iSecondaryReactions
            write(*,*) 'secondarySourceStoich(iBenthicTracers,iSecondaryReactions):',secondarySourceStoich(iBenthicTracers,iSecondaryReactions)
            write(*,*) 'secondarySinkStoich(iBenthicTracers,iSecondaryReactions):',secondarySinkStoich(iBenthicTracers,iSecondaryReactions)
         end do
      end do
   end if
   !NJ-END
     
   ! define stoichiometry of carbonate reactions

   call carbonateStoichMatrix

   !NJ-TEST : VERIFIED
   test_flag = .false.
   if (test_flag) then
      write(*,*) 'diagnostics: carbonatesStoichMatrix'
      do iBenthicTracers = 1,nBenthicTracers
         do iSecondaryReactions = 1,nCarbonateReactions
            write(*,*) 'iBenthicTracers, iSecondaryReactions:',iBenthicTracers,iSecondaryReactions
            write(*,*) 'carbonateSourceStoich(iBenthicTracers,iSecondaryReactions):',carbonateSourceStoich(iBenthicTracers,iSecondaryReactions)
            write(*,*) 'carbonateSinkStoich(iBenthicTracers,iSecondaryReactions):',carbonateSinkStoich(iBenthicTracers,iSecondaryReactions)
         end do
      end do
      write(*,*) 'end benthos_init'
   end if
   write(*,*) 'end benthos_init'
   !NJ-END
     
!-----------------------------------------------------------------------
!EOC

 end subroutine benthos_init

!*****************************************************************************
!BOP
! !IROUTINE: benthos_parameters_init
! !INTERFACE:

 subroutine benthos_parameters_init (nBenthicVertLevels)


! !DESCRIPTION:
!  Specify parameter profiles for diffusivity, porosity and the vertical grid
!
! !REVISION HISTORY:
!  same as module

   implicit none

! !INPUT PARAMETERS:

   integer (KIND=benthos_i4), intent(in) :: nBenthicVertLevels

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
   !-----------------------------------------------------------------------

   integer(KIND=benthos_i4) :: &
        k, &
        iElements

   real(KIND=benthos_r8) :: &
     dz    , &
     Dm    , & ! molecular diffusion in m2/s
     iDin_o    ! maximum biodiffusion in m2/s

   real(KIND=benthos_r8), dimension(nBenthicVertLevels+1) :: &
        tortuosity

   logical(KIND=benthos_log) :: &
        test_flag

   ! allocate parameter functions

   allocate(vertBenthosGridThickI(nBenthicVertLevels+1))
   allocate(benthosMidPointI(nBenthicVertLevels+1))
   allocate(benthosGridPointsI(nBenthicVertLevels+1))
   allocate(benthosGridPoints(nBenthicVertLevels+2))
   allocate(benthosPorosity(nBenthicVertLevels+2))
   allocate(benthosPorosityI(nBenthicVertLevels+1))
   allocate(benthosSolidPorosity(nBenthicVertLevels+2))
   allocate(benthosSolidPorosityI(nBenthicVertLevels+1))
   allocate(benthosDiffusivity(nBenthicVertLevels+2))
   allocate(benthosDiffusivityI(nBenthicVertLevels+1))
   allocate(benthosSolidDiffusivity(nBenthicVertLevels+2))
   allocate(benthosSolidDiffusivityI(nBenthicVertLevels+1))

   ! allocate primary reaction matrices

   allocate(primarySourceStoich(nBenthicTracers,nPrimaryReactions))
   allocate(primarySinkStoich(nBenthicTracers,nPrimaryReactions))

   ! allocate element ratios for tracer to net element conversion

   allocate(elementRatios(nBenthicTracers,nElements))

   ! carbon (mmol/ m^3), solid units  (mmol/kg)
   elementRatios(poca_ind,1) = c1_benthos*sediment_density
   elementRatios(pocb_ind,1) = c1_benthos*sediment_density
   elementRatios(pocc_ind,1) = c1_benthos*sediment_density
   elementRatios(dic_ind,1) = c1_benthos
   !elementRatios(co3_ind,1) = c1_benthos
   !elementRatios(co2_ind,1) = c1_benthos
   !elementRatios(hco3_ind,1) = c1_benthos
   elementRatios(caco3a_ind,1) = c1_benthos*sediment_density
   elementRatios(caco3b_ind,1) = c1_benthos*sediment_density
   elementRatios(camgco3_ind,1) = c1_benthos*sediment_density
   elementRatios(ch4_ind,1) = c1_benthos

  ! oxygen (mmol/m^3)
   elementRatios(poca_ind,2) = c1_benthos*sediment_density
   elementRatios(pocb_ind,2) = c1_benthos*sediment_density
   elementRatios(pocc_ind,2) = c1_benthos*sediment_density
   elementRatios(popa_ind,2) = 4.0_benthos_r8*sediment_density
   elementRatios(popb_ind,2) = 4.0_benthos_r8*sediment_density
   elementRatios(popc_ind,2) = 4.0_benthos_r8*sediment_density
   elementRatios(o2_ind,2) = c2_benthos
   elementRatios(h2po4_ind,2) = 4.0_benthos_r8
   elementRatios(co2_ind,2) = c2_benthos
   elementRatios(no3_ind,2) = 3.0_benthos_r8
   elementRatios(mno2a_ind,2) = c2_benthos*sediment_density
   elementRatios(mno2b_ind,2) = c2_benthos*sediment_density
   elementRatios(feoh3a_ind,2) = 3.0_benthos_r8*sediment_density
   elementRatios(feoh3b_ind,2) = 3.0_benthos_r8*sediment_density
   elementRatios(so4_ind,2) = 4.0_benthos_r8
   elementRatios(hco3_ind,2) = 3.0_benthos_r8
   elementRatios(caco3a_ind,2) = 3.0*sediment_density
   elementRatios(caco3b_ind,2) = 3.0*sediment_density
   elementRatios(co3_ind,2) = 3.0_benthos_r8

  ! nitrogen (mmol/m^3)
    elementRatios(nh4_ind,3) =c1_benthos
   elementRatios(pona_ind,3) = c1_benthos*sediment_density
   elementRatios(ponb_ind,3) = c1_benthos*sediment_density
   elementRatios(ponc_ind,3) = c1_benthos*sediment_density
     elementRatios(no3_ind,3) = c1_benthos

! phosphorus (mmol/m^3)
   elementRatios(popa_ind,4) = c1_benthos*sediment_density
   elementRatios(popb_ind,4) = c1_benthos*sediment_density
   elementRatios(popc_ind,4) = c1_benthos*sediment_density
   elementRatios(h2po4_ind,4) = c1_benthos
   elementRatios(fepa_ind,4) = c1_benthos*sediment_density
   elementRatios(fepb_ind,4) = c1_benthos*sediment_density

! sulfur (mmol/m3)
   elementRatios(so4_ind,5) = c1_benthos
   elementRatios(h2s_ind,5) = c1_benthos
   elementRatios(fes_ind,5) = c1_benthos*sediment_density
   elementRatios(fes2_ind,5) = c2_benthos*sediment_density
   elementRatios(s_ind,5) = c1_benthos*sediment_density

! manganese (mmol/m^3)
   elementRatios(mno2a_ind,6) = c1_benthos*sediment_density
   elementRatios(mno2b_ind,6) = c1_benthos*sediment_density
   elementRatios(mn_ind,6) = c1_benthos

! iron (mmol/m3)
   elementRatios(feoh3a_ind,7) = c1_benthos*sediment_density
   elementRatios(feoh3b_ind,7) = c1_benthos*sediment_density
   elementRatios(fe_ind,7) = c1_benthos
   elementRatios(fes_ind,7) = c1_benthos*sediment_density
   elementRatios(fes2_ind,7) = c1_benthos*sediment_density


   ! Define the benthos vertical grid

   dz = c1_benthos/real(nBenthicVertLevels,Kind=benthos_r8)
   vertBenthosGridThickI(:)      = c0_benthos
   benthosDiffusivity(:)         = c0_benthos
   benthosDiffusivityI(:)        = c0_benthos
   benthosSolidDiffusivity(:)    = c0_benthos
   benthosSolidDiffusivityI(:)   = c0_benthos
   !benthosStorageConc = linspace(1,nBenthicTracers,nBenthicTracers)*0   !!!  (mmol/m3)

   benthosGridPointsI(1) = c0_benthos
   benthosGridPoints(1) = c0_benthos
   benthosGridPoints(nBenthicVertLevels+2) = c1_benthos
   vertBenthosGridThickI(1) = dz/c2_benthos

   benthosMidPointI(1) = benthosDepth*(benthosGridPointsI(1) + vertBenthosGridThickI(1)*p5_benthos)
   do k = 2, nBenthicVertLevels+1
      vertBenthosGridThickI(k) = dz
      benthosGridPointsI(k) = benthosGridPointsI(k-1) + dz
      benthosGridPoints(k) = benthosGridPointsI(k) - dz/c2_benthos
      benthosMidPointI(k) = benthosDepth*(benthosGridPointsI(k) + vertBenthosGridThickI(k)*p5_benthos)
   enddo
   vertBenthosGridThickI(nBenthicVertLevels+1) = dz/c2_benthos
   benthosMidPointI(nBenthicVertLevels+1) = benthosDepth*(benthosGridPointsI(nBenthicVertLevels+1) + vertBenthosGridThickI(nBenthicVertLevels+1)*p5_benthos)

   ! Initialize porosity profile (Reed et al 2012)

   if (useDepthDependentPorosity) then

      benthosPorosityI(1) = porosity_o
      benthosPorosityI(nBenthicVertLevels+1) = porosity_inf
      do k = 2, nBenthicVertLevels
         benthosPorosityI(k) = porosity_inf + (porosity_o-porosity_inf)*exp(-(k-1)*dz*benthosDepth/porosity_gamma)
      enddo

      benthosPorosity(1) = c1_benthos
      benthosPorosity(nBenthicVertLevels+2) = porosity_inf
      benthosPorosity(2) = benthosPorosityI(1)
      benthosPorosity(nBenthicVertLevels+1) = benthosPorosityI(nBenthicVertLevels+1) 

      do k = 3, nBenthicVertLevels
         benthosPorosity(k) = (benthosPorosityI(k-1) + benthosPorosityI(k))/c2_benthos
      enddo

   else
      benthosPorosityI(:) = p5_benthos
      benthosPorosity(:)  = p5_benthos
   endif
 
   ! Define  solid porosity

   do k = 1, nBenthicVertLevels+2
      benthosSolidPorosityI(k) = c1_benthos - benthosPorosityI(k)
      benthosSolidPorosity(k) = c1_benthos - benthosPorosity(k)

      ! Tortuosity

      tortuosity(k) = c1_benthos - log(benthosPorosityI(k))
   enddo

   benthosSolidPorosity(1) = benthosSolidPorosityI(1)
   benthosSolidPorosity(nBenthicVertLevels+2) = benthosSolidPorosity(nBenthicVertLevels+1)

   ! Define diffusivity (m2/s)

   Dm = molecular_diff/sec_per_year
   iDin_o = max_bio_diff*m2percm2/sec_per_year

   if (useNonZeroDiffusivity) then
      if (useConstantDiffusivity) then
         benthosDiffusivityI(:)      = iDin_o
         benthosSolidDiffusivityI(:) = iDin_o
         benthosDiffusivity(:)       = iDin_o
         benthosSolidDiffusivity(:)  = iDin_o
      else
         benthosDiffusivityI(1) = benthosPorosityI(1)*(Dm/tortuosity(1)+iDin_o)/benthosDepth**2  
         benthosSolidDiffusivityI(1) = benthosSolidPorosityI(1)*iDin_o/benthosDepth**2

         do k = 2,nBenthicVertLevels+1
            benthosDiffusivityI(k) = benthosPorosityI(k)*(Dm/tortuosity(k) + iDin_o *exp(-benthosGridPointsI(k)*benthosDepth/x_biodiffusion))/benthosDepth**2
            benthosSolidDiffusivityI(k) =benthosSolidPorosityI(k)*(iDin_o *exp(-benthosGridPointsI(k)*benthosDepth/x_biodiffusion))/benthosDepth**2
         enddo

         benthosDiffusivity(1) = (Dm+iDin_o)/benthosDepth**2
         benthosDiffusivity(2) = benthosDiffusivityI(1)
         benthosDiffusivityi(nBenthicVertLevels+2) = benthosDiffusivityI(nBenthicVertLevels+1)
         benthosSolidDiffusivity(nBenthicVertLevels+1) = benthosDiffusivityI(nBenthicVertLevels+1)

         benthosSolidDiffusivity(1) = benthosSolidDiffusivityI(1)
         benthosSolidDiffusivity(2) = benthosSolidDiffusivityI(1)
         benthosSolidDiffusivity(nBenthicVertLevels+2) = benthosSolidDiffusivityI(nBenthicVertLevels+1)
         benthosSolidDiffusivity(nBenthicVertLevels+1) = benthosSolidDiffusivityI(nBenthicVertLevels+1)

         do k = 3,nBenthicVertLevels
            benthosDiffusivity(k) = (benthosDiffusivityI(k-1) + benthosDiffusivityI(k))/c2_benthos
            benthosSolidDiffusivity(k) = (benthosSolidDiffusivityI(k-1) + benthosSolidDiffusivityI(k))/c2_benthos
         enddo

     endif   ! constant diffusivity
  endif
  test_flag = .false.
  if (test_flag) then
     do k = 1,nBenthicVertLevels+1
        write(*,*) 'benthosDiffusivity(k):', benthosDiffusivity(k)
        write(*,*) 'benthosSolidDiffusivity(k):', benthosSolidDiffusivity(k)
        write(*,*) 'benthosPorosity(k):', benthosPorosity(k)
        write(*,*) 'benthosSolidPorosity(k):', benthosSolidPorosity(k)
     end do
  end if

  
!
!  Initialize the deep benthic storage concentrations

! benthosStorageConc(:) = 0


!EOC

 end subroutine benthos_parameters_init

!*****************************************************************************
!BOP
! !IROUTINE: benthos_parameters_init
! !INTERFACE:

 subroutine benthos_parameters_dt (nBenthicVertLevels,oceanBottomTemperature)

! !DESCRIPTION:
!  Update time-dependent parameters
! 
!
! !REVISION HISTORY:
!  same as module

   implicit none

! !INPUT PARAMETERS:

   integer (KIND=benthos_i4), intent(in) :: nBenthicVertLevels

   real(KIND=benthos_r8), intent(in) :: oceanBottomTemperature

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer(KIND=benthos_i4) :: &
        k, &
        iElements

   real(KIND=benthos_r8) :: &
     dz    , &
     Dm    , & ! molecular diffusion in m2/s
     iDin_o    ! maximum biodiffusion in m2/s

   real(KIND=benthos_r8), dimension(nBenthicVertLevels+1) :: &
        tortuosity

   logical(KIND=benthos_log) :: &
        test_flag

   benthosDiffusivity(:)         = c0_benthos
   benthosDiffusivityI(:)        = c0_benthos
   benthosSolidDiffusivity(:)    = c0_benthos
   benthosSolidDiffusivityI(:)   = c0_benthos
   
   do k = 1, nBenthicVertLevels+2
      ! Tortuosity
      tortuosity(k) = c1_benthos - log(benthosPorosityI(k))
   enddo

   ! Define diffusivity (m2/s)

   Dm = molecular_diff/sec_per_year
   iDin_o = max_bio_diff*m2percm2/sec_per_year*abs(oceanBottomTemperature/T_max_biodiffusion)

   if (useNonZeroDiffusivity) then
      if (useConstantDiffusivity) then
         benthosDiffusivityI(:)      = iDin_o
         benthosSolidDiffusivityI(:) = iDin_o
         benthosDiffusivity(:)       = iDin_o
         benthosSolidDiffusivity(:)  = iDin_o
      else
         benthosDiffusivityI(1) = benthosPorosityI(1)*(Dm/tortuosity(1)+iDin_o)/benthosDepth**2  
         benthosSolidDiffusivityI(1) = benthosSolidPorosityI(1)*iDin_o/benthosDepth**2

         do k = 2,nBenthicVertLevels+1
            benthosDiffusivityI(k) = benthosPorosityI(k)*(Dm/tortuosity(k) + iDin_o *exp(-benthosGridPointsI(k)*benthosDepth/x_biodiffusion))/benthosDepth**2
            benthosSolidDiffusivityI(k) =benthosSolidPorosityI(k)*(iDin_o *exp(-benthosGridPointsI(k)*benthosDepth/x_biodiffusion))/benthosDepth**2
         enddo

         benthosDiffusivity(1) = (Dm+iDin_o)/benthosDepth**2
         benthosDiffusivity(2) = benthosDiffusivityI(1)
         benthosDiffusivityi(nBenthicVertLevels+2) = benthosDiffusivityI(nBenthicVertLevels+1)
         benthosSolidDiffusivity(nBenthicVertLevels+1) = benthosDiffusivityI(nBenthicVertLevels+1)

         benthosSolidDiffusivity(1) = benthosSolidDiffusivityI(1)
         benthosSolidDiffusivity(2) = benthosSolidDiffusivityI(1)
         benthosSolidDiffusivity(nBenthicVertLevels+2) = benthosSolidDiffusivityI(nBenthicVertLevels+1)
         benthosSolidDiffusivity(nBenthicVertLevels+1) = benthosSolidDiffusivityI(nBenthicVertLevels+1)

         do k = 3,nBenthicVertLevels
            benthosDiffusivity(k) = (benthosDiffusivityI(k-1) + benthosDiffusivityI(k))/c2_benthos
            benthosSolidDiffusivity(k) = (benthosSolidDiffusivityI(k-1) + benthosSolidDiffusivityI(k))/c2_benthos
         enddo

     endif   ! constant diffusivity
  endif

!EOC

 end subroutine benthos_parameters_dt

!*****************************************************************************
!BOP
! !IROUTINE: carbonateStoichMatrix
! !INTERFACE:

 subroutine carbonateStoichMatrix


! !DESCRIPTION:
!  Specify mole fractions for sources and sinks in the following carbonate reactions
!
!  carbonate reactions
! 3 total carbonate equations
!
!
! p1. CaCO3 (calcite) --> Ca + CO3
! p2. CaCO3 (aragonite) --> Ca + CO3
! p3. Ca.85Mg.15CO3 --> 0.85 Ca + 0.15 Mg + CO3
!
! !REVISION HISTORY:
!  same as module

   implicit none

! !INPUT PARAMETERS:
!
!EOP
!BOC

   ! allocate parameter functions

   allocate(carbonateSourceStoich(nBenthicTracers,nCarbonateReactions))
   allocate(carbonateSinkStoich(nBenthicTracers,nCarbonateReactions))

   carbonateSourceStoich(:,:) = c0_benthos
   carbonateSinkStoich(:,:) = c0_benthos

   ! SINK TERMS: Rates are in mmol/m3/s
   carbonateSinkStoich(caco3a_ind,1)  = c1_benthos/sediment_density
   carbonateSinkStoich(caco3b_ind,2)  = c1_benthos/sediment_density
   carbonateSinkStoich(camgco3_ind,3) = c1_benthos/sediment_density

   ! SOURCE TERMS
   carbonateSourceStoich(co3_ind,:) = c1_benthos
   carbonateSourceStoich(dic_ind,:) = c1_benthos
   carbonateSourceStoich(alk_ind,:) = c2_benthos

 end subroutine carbonateStoichMatrix

!*****************************************************************************
!BOP
! !IROUTINE: secondaryStoichMatrix
! !INTERFACE:

 subroutine secondaryStoichMatrix


! !DESCRIPTION:
!  Specify mole fractions for sources and sinks in the following secondary reactions

!  secondary reactions
!
! s1. O2 + NH4/2 + HCO3 --> NO3/2 + CO2 + 3/2H2O
! s2. O2 + 2Mn + 4HCO3 --> 2MnO2^a + 4CO2 + 2H2O
! s3. O2 + 4Fe + 8HCO3 + 2 H20 + 4 gamma H2PO4 + 24 gamma H^+ -->
!          4Fe(OH)3^a + 4 gamma [Fe-P]^a + 8CO2 + 16 gamma H2O
! s4. O2 + FeS/2 --> SO4/2 + Fe/2
! s5. O2+ 2/7FeS2 +  2/7H2O --> 4/7SO4 +  2/7Fe+ 4/7H
! s6. O2 + H2S/2 + HCO3 -->  SO4/2 + CO2 + H2O
! s7. O2 + CH4 --> CO2 + 4H^+
! s8. MnO2a + 2Fe + 2 gamma H2PO4 + 2H2O + 2HCO3 + 12gamma H -->
!             2Fe(OH)3^a +  2 gamma [Fe- P]^a + Mn + 2CO2 + 8gamma H2O
! (s9)s8b. MnO2b + 2Fe + 2 gamma H2PO4 + 2H2O + 2HCO3 + 12gammaH -->
!             2Fe(OH)3^a +  2 gamma [Fe- P]^a + Mn + 2CO2 + 8gamma H2O
! (s10)s9. Mn02a + H2S +  2CO2 -->  Mn +  S +  2HCO3
! (s11)s9b. Mn02b + H2S +  2CO2 -->  Mn +  S +  2HCO3
! (s12)s10. Fe(OH)3a+   gamma[Fe-P]a + H2S/2+ 2CO2 + 4 gamma H2O-->
!             Fe +   gamma H2PO4 + S/2 + 2HCO3 + (2+6gamma)H
! (s13)s10b. Fe(OH)3b+   gamma[Fe-P]b + H2S/2+ 2CO2 + 4 gamma H2O-->
!             Fe +  gamma H2PO4+ S/2 + 2HCO3 + (2+6 gamma)H
! (s14)s11. Fe + H2S -->  FeS + 2H
! (s15)s12. SO4 +  CH4 +  CO2 -->  2HCO3 + H2S
! (s16)s13. S +  H2O -->  3/4H2S +  SO4/4 +  H/2
! (s17)s14. FeS +  S --> FeS2
! (s18)s15. Fe(OH)3a + gamma [Fe-P]a --> Fe(OH)3b + gamma [Fe-P]b
! (s19)s16 MnO2a --> MnO2b
!
! 19 total secondary equations with 3 repetitions
! gamma = "ironBoundPFraction" is the fraction of  co-precipitated P with iron
!
!
! !REVISION HISTORY:
!  same as module

   implicit none

! !INPUT PARAMETERS:
!
!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

   ! allocate parameter functions

   allocate(secondarySourceStoich(nBenthicTracers,nSecondaryReactions))
   allocate(secondarySinkStoich(nBenthicTracers,nSecondaryReactions))

   secondarySourceStoich(:,:) = c0_benthos
   secondarySinkStoich(:,:) = c0_benthos
   
   ! add units conversion for secondary rates which are in mmol/m3/s except for the last
   ! three reactions (s17-s19) which only involve solid phase material. These are in mmol/kg/s

   ! SINK TERMS
   ! s1
   secondarySinkStoich(o2_ind,1) = c1_benthos
   !secondarySinkStoich(hco3_ind,1) = c1_benthos
   secondarySinkStoich(dic_ind,1) = c1_benthos
   secondarySinkStoich(nh4_ind,1) = p5_benthos
   secondarySinkStoich(alk_ind,1) = c1_benthos

   !s2

   secondarySinkStoich(mn_ind,2) = c2_benthos
   !secondarySinkStoich(hco3_ind,2) = 4.0_benthos_r8
   secondarySinkStoich(dic_ind,2) = 4.0_benthos_r8
   secondarySinkStoich(o2_ind,2) = c1_benthos
   secondarySinkStoich(alk_ind,2) = 4.0_benthos_r8

   ! s3
   secondarySinkStoich(o2_ind,3) = c1_benthos
   !secondarySinkStoich(hco3_ind,3) = 8.0_benthos_r8
   secondarySinkStoich(fe_ind,3) = 4.0_benthos_r8
   secondarySinkStoich(h2po4_ind,3) = 4.0_benthos_r8*ironBoundPFraction
   secondarySinkStoich(alk_ind,3) = 8.0_benthos_r8

   !s4
   secondarySinkStoich(o2_ind,4) = c1_benthos
   secondarySinkStoich(fes_ind,4) = p5_benthos/sediment_density

   !s5
   secondarySinkStoich(o2_ind,5) = c1_benthos
   secondarySinkStoich(fes2_ind,5) = c2_benthos/7.0_benthos_r8/sediment_density
   secondarySinkStoich(alk_ind,5) = 4.0_benthos_r8/7.0_benthos_r8

   !s6
   secondarySinkStoich(o2_ind,6) = c1_benthos
   !secondarySinkStoich(hco3_ind,6) = c1_benthos
   secondarySinkStoich(dic_ind,6) = c1_benthos
   secondarySinkStoich(h2s_ind,6) = p5_benthos
   secondarySinkStoich(alk_ind,6) = c1_benthos

   !s7
   secondarySinkStoich(o2_ind,7) = c1_benthos
   secondarySinkStoich(ch4_ind,7) = c1_benthos

   !s8 and s9
   secondarySinkStoich(mno2a_ind,8) = c1_benthos/sediment_density
   secondarySinkStoich(fe_ind,8) = c2_benthos
   secondarySinkStoich(fe_ind,9) = c2_benthos
   secondarySinkStoich(h2po4_ind,8) = c2_benthos*ironBoundPFraction
   secondarySinkStoich(h2po4_ind,9) = c2_benthos*ironBoundPFraction
   !secondarySinkStoich(hco3_ind,8:9) = c2_benthos
   secondarySinkStoich(dic_ind,8) = c2_benthos
   secondarySinkStoich(dic_ind,9) = c2_benthos
   secondarySinkStoich(mno2b_ind,9) = c1_benthos/sediment_density
   secondarySinkStoich(alk_ind,8) = c2_benthos
   secondarySinkStoich(alk_ind,9) = c2_benthos

   !s10 and s11
   secondarySinkStoich(h2s_ind,10) = c1_benthos
   secondarySinkStoich(h2s_ind,11) = c1_benthos
   secondarySinkStoich(mno2a_ind,10) = c1_benthos/sediment_density
   secondarySinkStoich(mno2b_ind,11) = c1_benthos/sediment_density
   !secondarySinkStoich(co2_ind,10:11) = c2_benthos
   secondarySinkStoich(dic_ind,10) = c2_benthos
   secondarySinkStoich(dic_ind,11) = c2_benthos

   !s12 and s13
   secondarySinkStoich(h2s_ind,12) = p5_benthos
   secondarySinkStoich(h2s_ind,13) = p5_benthos
   !econdarySinkStoich(co2_ind,12:13) = c2_benthos
   secondarySinkStoich(dic_ind,12) = c2_benthos
   secondarySinkStoich(dic_ind,13) = c2_benthos
   secondarySinkStoich(feoh3a_ind,12) = c1_benthos/sediment_density
   secondarySinkStoich(feoh3b_ind,13) = c1_benthos/sediment_density
   secondarySinkStoich(fepa_ind,12) = c1_benthos*ironBoundPFraction/sediment_density
   secondarySinkStoich(fepb_ind,13) = c1_benthos*ironBoundPFraction/sediment_density
   secondarySinkStoich(alk_ind,12) = 6.0_benthos_r8*ironBoundPFraction
   secondarySinkStoich(alk_ind,13) = 6.0_benthos_r8*ironBoundPFraction

   !s14
   secondarySinkStoich(h2s_ind,14) = c1_benthos
   secondarySinkStoich(fe_ind,14) = c1_benthos
   secondarySinkStoich(alk_ind,14) = c2_benthos

   !s15
   secondarySinkStoich(ch4_ind,15) = c1_benthos
   !secondarySinkStoich(co2_ind,15) = c1_benthos
   secondarySinkStoich(so4_ind,15) = c1_benthos
   secondarySinkStoich(dic_ind,15) = c1_benthos

   !s16
   secondarySinkStoich(s_ind,16) = c1_benthos/sediment_density

   !s17
   secondarySinkStoich(fes_ind,17) = c1_benthos
   secondarySinkStoich(s_ind,17) = c1_benthos

   !s18
   secondarySinkStoich(feoh3a_ind,18) = c1_benthos
   secondarySinkStoich(fepa_ind,18) = c1_benthos*ironBoundPFraction

   !s19
   secondarySinkStoich(mno2a_ind,19) = c1_benthos

   ! SOURCE TERMS
   !s1
   secondarySourceStoich(no3_ind,1) = p5_benthos
   !secondarySourceStoich(co2_ind,1) = c1_benthos
   secondarySourceStoich(dic_ind,1) = c1_benthos

   !s2
   secondarySourceStoich(mno2a_ind,2) = c2_benthos/sediment_density
   !secondarySourceStoich(co2_ind,2) = 4.0_benthos_r8
   secondarySourceStoich(dic_ind,2) = 4.0_benthos_r8

   !s3
   secondarySourceStoich(feoh3a_ind,3) = 4.0_benthos_r8/sediment_density
   secondarySourceStoich(fepa_ind,3) = 4.0_benthos_r8/sediment_density*ironBoundPFraction
   !secondarySourceStoich(co2_ind,3) = 8.0_benthos_r8
   secondarySourceStoich(dic_ind,3) = 8.0_benthos_r8
   secondarySourceStoich(alk_ind,3) = 24.0_benthos_r8*ironBoundPFraction

   !s4
   secondarySourceStoich(so4_ind,4) = p5_benthos
   secondarySourceStoich(fe_ind,4) = p5_benthos

   !s5
   secondarySourceStoich(so4_ind,5) = 4.0_benthos_r8/7.0_benthos_r8
   secondarySourceStoich(fe_ind,5) = c2_benthos/7.0_benthos_r8

   !s6
   secondarySourceStoich(so4_ind,6) = p5_benthos
   !secondarySourceStoich(co2_ind,6) = c1_benthos
   secondarySourceStoich(dic_ind,6) = c1_benthos

   !s7
   !secondarySourceStoich(co2_ind,7) = c1_benthos
   secondarySourceStoich(dic_ind,7) = c1_benthos
   secondarySourceStoich(alk_ind,7) = 4.0_benthos_r8

   !s8
   secondarySourceStoich(feoh3a_ind,8) = c2_benthos/sediment_density
   secondarySourceStoich(fepa_ind,8) = c2_benthos/sediment_density*ironBoundPFraction
   secondarySourceStoich(mn_ind,8) = c1_benthos
   !secondarySourceStoich(co2_ind,8) = c2_benthos
   secondarySourceStoich(dic_ind,8) = c2_benthos
   secondarySourceStoich(alk_ind,8) = 12.0_benthos_r8*ironBoundPFraction

   !s9
   secondarySourceStoich(feoh3a_ind,9) = c2_benthos/sediment_density
   secondarySourceStoich(fepa_ind,9) = c2_benthos/sediment_density*ironBoundPFraction
   secondarySourceStoich(mn_ind,9) = c1_benthos
   !secondarySourceStoich(co2_ind,9) = c2_benthos
   secondarySourceStoich(dic_ind,9) = c2_benthos
   secondarySourceStoich(alk_ind,9) = 12.0_benthos_r8*ironBoundPFraction

   !s10
   secondarySourceStoich(mn_ind,10) = c1_benthos
   secondarySourceStoich(s_ind,10) = c1_benthos/sediment_density
   !secondarySourceStoich(hco3_ind,10) = c2_benthos
   secondarySourceStoich(dic_ind,10) = c2_benthos
   secondarySourceStoich(alk_ind,10) = c2_benthos

   !s11
   secondarySourceStoich(mn_ind,11) = c1_benthos
   secondarySourceStoich(s_ind,11) = c1_benthos/sediment_density
   !secondarySourceStoich(hco3_ind,11) = c2_benthos
   secondarySourceStoich(dic_ind,11) = c2_benthos
   secondarySourceStoich(alk_ind,11) = c2_benthos

   !s12
   secondarySourceStoich(fe_ind,12) = c1_benthos
   secondarySourceStoich(h2po4_ind,12) = ironBoundPFraction
   secondarySourceStoich(s_ind,12) = p5_benthos/sediment_density
   !secondarySourceStoich(hco3_ind,12) = c2_benthos
   secondarySourceStoich(dic_ind,12) = c2_benthos

   !s13
   secondarySourceStoich(fe_ind,13) = c1_benthos
   secondarySourceStoich(h2po4_ind,13) = ironBoundPFraction
   secondarySourceStoich(s_ind,13) = p5_benthos/sediment_density
   !secondarySourceStoich(hco3_ind,13) = c2_benthos
   secondarySourceStoich(dic_ind,13) = c2_benthos

   !s14
   secondarySourceStoich(fes_ind,14) = c1_benthos/sediment_density

   !s15
   !secondarySourceStoich(hco3_ind,15) = c2_benthos
   secondarySourceStoich(h2s_ind,15) = c1_benthos
   secondarySourceStoich(dic_ind,15) = c2_benthos

   !s16
   secondarySourceStoich(h2s_ind,16) = 3.0_benthos_r8/4.0_benthos_r8
   secondarySourceStoich(so4_ind,16) = 0.25_benthos_r8

   !s17
   secondarySourceStoich(fes2_ind,17) = c1_benthos

   !s18
   secondarySourceStoich(feoh3b_ind,18) = c1_benthos
   secondarySourceStoich(fepb_ind,18) = c1_benthos*ironBoundPFraction

   !s19
   secondarySourceStoich(mno2b_ind,19) = c1_benthos

 end subroutine secondaryStoichMatrix
 
!*****************************************************************************
!BOP
! !IROUTINE: computeNetElements
! !INTERFACE:

subroutine computeNetElements (netElements, nBenthicVertLevels, benthosTracerBulk)

! !DESCRIPTION:
!  Initialize tracegas tracer module.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer(KIND=benthos_i4), intent(in) :: nBenthicVertLevels
  real(KIND=benthos_r8), dimension(:,:), intent(in) :: benthosTracerBulk

! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(nElements), intent(out) :: netElements

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------
!   real(KIND=benthos_r8), dimension(:,:), allocatable :: &
!     elementRatios

   real(KIND=benthos_r8), dimension(nBenthicTracers) :: &
     columnBenthosTracerBulk

   integer (KIND=benthos_i4) :: &
     iTracers, &
     iElements

   columnBenthosTracerBulk(:) = c0_benthos

   do iTracers = 1,nBenthicTracers
      call sum_benthos_column(columnBenthosTracerBulk(iTracers), benthosTracerBulk(iTracers,:),nBenthicVertLevels)
   end do
   do iElements = 1,nElements
      do iTracers = 1,nBenthicTracers
           netElements(iElements) = netElements(iElements) + elementRatios(iTracers,iElements)*columnBenthosTracerBulk(iTracers)
      end do
   end do

 end subroutine computeNetElements
 
!*****************************************************************************
!BOP
! !IROUTINE: sum_benthos_column
! !INTERFACE:

  subroutine sum_benthos_column(columnBenthosTracerBulk, benthosTracerBulk,nBenthicVertLevels)

! !DESCRIPTION:
!  Initialize tracegas tracer module.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer(KIND=benthos_i4), intent(in) :: nBenthicVertLevels

  real(KIND=benthos_r8), dimension(:), intent(in):: benthosTracerBulk

! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:

  real(KIND=benthos_r8), intent(out) :: columnBenthosTracerBulk

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------
   integer(KIND=benthos_i4) :: &
     iLevels

   columnBenthosTracerBulk = c0_benthos

   do iLevels = 1,nBenthicVertLevels+1
      columnBenthosTracerBulk = columnBenthosTracerBulk + benthosTracerBulk(iLevels)*vertBenthosGridThickI(iLevels)*benthosDepth
   end do

 end subroutine sum_benthos_column

!*****************************************************************************
!BOP
! !IROUTINE: updateBenthosBoundaryConditions
! !INTERFACE:

 subroutine updateBenthosBoundaryConditions  (oceanSedFlux, oceanPrecipFlux, totalOceanSolidFlux, &
      oceanTracerConc, deepBenthosBottomFlux, oceanBottomDensity, nBenthicVertLevels,FluxSed, dhtop, vel, dt, &
      oceanBottomPhosphate, oceanBottomAmmonium, oceanBottomNitrate, oceanBottomOxygen, &
      oceanBottomIron, oceanBottomDIC, oceanBottomAlkalinity, oceanPOCFlux, oceanPONFlux, &
      oceanPOPFlux, oceanParticulateIronFlux, oceanCalciteFlux)

! !DESCRIPTION:
!  compute benthos boundary fluxes and concentrations at the start of each timestep
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer(KIND=benthos_i4), intent(in) :: &
     nBenthicVertLevels  ! number of vertical interfaces in the benthos

  real(KIND=benthos_r8), intent(in) :: &
       dt, &
       oceanBottomDensity, &
       oceanBottomPhosphate, &
       oceanBottomAmmonium, &
       oceanBottomNitrate, &
       oceanBottomOxygen, &
       oceanBottomIron, &
       oceanBottomDIC, &
       oceanBottomAlkalinity, &
       oceanPOCFlux, &
       oceanPOPFlux, &
       oceanPONFlux, &
       oceanParticulateIronFlux, &
       oceanCalciteFlux

! !INPUT/OUTPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(:), intent(inout) :: &
     oceanSedFlux,    & ! particulate flux of sinking ocean tracers (mmol/kg*m/s)
     oceanTracerConc, & ! bottom ocean concentration of solutes (mmol/m3)
     oceanPrecipFlux, & ! flux of precipitated material from ocean bottom (mmol/kg*m/s)
     totalOceanSolidFlux, & ! total solid flux divided by benthosDepth (mmol/kg/s)
     deepBenthosBottomFlux  ! prescribed deep storage flux (mmol/m2/s)

! !OUTPUT PARAMETERS:

  real(KIND=benthos_r8), intent(out) :: &
     FluxSed, & ! kg/m2/s
     dhtop,   & ! (m) increase in benthos surface over the timestep
     vel        ! (1/s) normalized velocity of the surface and bottom boundaries

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

  real (KIND=benthos_r8) :: &
       pocFluxRate, &  ! mmol (m/kg)/s
       ponFluxRate, &
       popFluxRate

!-----------------------------------------------------------------------
!EOC

   ! pocFluxRate = pocFluxRate_o/sediment_density
   pocFluxRate = oceanPOCFlux/sediment_density
   ponFluxRate = oceanPONFlux/sediment_density
   popFluxRate = oceanPOPFlux/sediment_density

   ! sediment sinking flux will eventually come from the ocean model
   !  sea bed accumulation rates (ob, Yenisey  0.2-0.7 cm/y

   if (useSedimentation) then
      if (useFastSedimentation) then
         fluxSed = 1000.0_benthos_r8*fluxSed_o
      else
         fluxSed = fluxSed_o
      end if
      if (.not. usePositiveSedimentation) fluxSed = -fluxSed
   else
      fluxSed = c0_benthos
   end if

   dhtop = fluxSed/sediment_density * dt
   vel = dhtop/benthosDepth/dt

   if (useDeepSource) then
   ! only a couple have a source term  (H2S has an upward flux -50.0  mmol/m2/y

      deepBenthosBottomFlux(h2s_ind) = deepFlux_h2s
   else
      deepBenthosBottomFlux(h2s_ind) = c0_benthos
   end if

   if (useBGCSinkingFlux) then

      oceanSedFlux(poca_ind) = pocFluxRate*fracHighlyReactive   ! mmol (m/kg)/s  (49% highy reactive, Westrich and Berner 1984)
      oceanSedFlux(pocb_ind) = pocFluxRate*fracLessReactive ! 15% moderately reactive
      oceanSedFlux(pocc_ind) = pocFluxRate*fracRefractory ! 34% nonreactive
      
      oceanSedFlux(popa_ind) = popFluxRate*fracHighlyReactive   ! mmol (m/kg)/s  (49% highy reactive, Westrich and Berner 1984)
      oceanSedFlux(popb_ind) = popFluxRate*fracLessReactive ! 15% moderately reactive
      oceanSedFlux(popc_ind) = popFluxRate*fracRefractory ! 34% nonreactive
      
      oceanSedFlux(pona_ind) = ponFluxRate*fracHighlyReactive   ! mmol (m/kg)/s  (49% highy reactive, Westrich and Berner 1984)
      oceanSedFlux(ponb_ind) = ponFluxRate*fracLessReactive ! 15% moderately reactive
      oceanSedFlux(ponc_ind) = ponFluxRate*fracRefractory ! 34% nonreactive
      oceanSedFlux(mno2a_ind) = 0.03_benthos_r8/sec_per_year/sediment_density !mmol/m2/y to mmol (m/kg)/s
      oceanSedFlux(feoh3a_ind) = oceanParticulateIronFlux/sediment_density  !57.5_benthos_r8/sec_per_year/sediment_density! !mmol/m2/y to mmol (m/kg)/s
      oceanSedFlux(fepa_ind) = oceanSedFlux(feoh3a_ind)*ironBoundPFraction !
      oceanSedFlux(caco3a_ind) = oceanCalciteFlux * fracCalcite  !300.0_benthos_r8/sec_per_year/sediment_density! mmol/m2/y to mmol (m/kg)/s Krumins
      oceanSedFlux(caco3b_ind) = oceanCalciteFlux * fracAragonite !Calcite150.0_benthos_r8/sec_per_year/sediment_density! mmol/m2/y  Krumins to mmol (m/kg)/s
      oceanSedFlux(camgco3_ind) = oceanCalciteFlux * fracMgCalcite !100.0_benthos_r8/sec_per_year/sediment_density! to mmol (m/kg)/s

   else
      oceanSedFlux(:) = c0_benthos
   end if

   ! initial Ocean tracer concentrations  : zero for solids
   ! non-zero for solutes  (some from ocean model)
   !
   ! concentrations in mmol/m3 (umol/L)
   !
   oceanTracerConc(:) = c0_benthos

   if (useOceanConc) then
      oceanTracerConc(o2_ind) = oceanBottomOxygen !350.0_benthos_r8  !umol/L (mmol/m3) Reed et al, 2011
      oceanTracerConc(nh4_ind) = oceanBottomAmmonium !c2_benthos  !umol/L  (Reed)
      oceanTracerConc(h2po4_ind) = oceanBottomPhosphate !c1_benthos  !umol/L Reed et al. 2011 
      oceanTracerConc(co2_ind) = 0.125_benthos_r8*mM_umolperL  !mM Krumins et al, 2013
      oceanTracerConc(no3_ind) = oceanBottomNitrate !6.0_benthos_r8 ! umol/L Reed.
      oceanTracerConc(mn_ind) = c0_benthos ! umol/L
      oceanTracerConc(fe_ind) = oceanBottomIron !c1_benthos ! umol/L Reed
      oceanTracerConc(so4_ind) = 12.0_benthos_r8*mM_umolperL ! mmol/L Reed  to mmol/m3
      oceanTracerConc(h2s_ind) = c0_benthos ! Reed
      oceanTracerConc(ch4_ind) = c0_benthos !
      oceanTracerConc(hco3_ind) = c2_benthos*mM_umolperL  !mM  Krumins et al  (essentially *DIC*) and TA?
      oceanTracerConc(co3_ind) = 70.0_benthos_r8*oceanBottomDensity !bottom_water_density  ! mmol/m3
      !umol/kg*density of water 1.027 g/m3  (Sulpis et al 2018)
      oceanTracerConc(dic_ind) = oceanBottomDIC !oceanTracerConc(co3_ind) + oceanTracerConc(co2_ind) + oceanTracerConc(hco3_ind)
      oceanTracerConc(alk_ind) = oceanBottomAlkalinity !2306.0_benthos_r8  ! mmol/m3  Krumins et al.
   else
      oceanTracerConc(:) = c0_benthos
   end if

   !  make sure dic is consistent
   !initialBenthicTracerBulk(dic_ind,:) = initialBenthicTracerBulk(hco3_ind,:) + initialBenthicTracerBulk(co2_ind,:)  +  initialBenthicTracerBulk(co3_ind,:)

   ! define ocean precipitation flux: ignore this right now
   oceanPrecipFlux(:) = 0
   oceanPrecipFlux(mno2a_ind) = 0.09_benthos_r8*oceanTracerConc(mn_ind)*oceanTracerConc(o2_ind) * k_s2 *  c2_benthos/sediment_density*benthosDepth
   oceanPrecipFlux(feoh3a_ind) = 0.04_benthos_r8*oceanTracerConc(fe_ind) * oceanTracerConc(o2_ind)*k_s3 * 4.0_benthos_r8/sediment_density*benthosDepth
   oceanPrecipFlux(fepa_ind) = oceanPrecipFlux(feoh3a_ind)*ironBoundPFraction

   totalOceanSolidFlux(:) = (oceanPrecipFlux(:) + oceanSedFlux(:))/benthosDepth

 end subroutine updateBenthosBoundaryConditions

!*****************************************************************************
!BOP
! !IROUTINE: benthos_flags_init
! !INTERFACE:

 subroutine benthos_flags_init (setBenthosFlags, useCarbonateSaturation, useSecondaryReactions, &
            useDepthDependentPorosity, useNonZeroDiffusivity, useSedimentation, &
            useFastSedimentation, usePositiveSedimentation, useBgcSinkingFlux, useStepInitialProfiles, &
            useBenthicReactions, useFluxCorrection, useOceanConc, useDeepSource, &
            useConstantDiffusivity)

! !DESCRIPTION:
!  Initialize flags for testcases
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer (KIND=benthos_i4), intent(in) :: setBenthosFlags

! !INPUT/OUTPUT PARAMETERS:

  logical (KIND=benthos_log), intent(inout) :: &
       useCarbonateSaturation,    & ! true turns on carbonate chemistry
       useSecondaryReactions,     & ! true computes secondary reaction set
       useDepthDependentPorosity, & ! true uses a depth dependent porosity
       useNonZeroDiffusivity,     & ! false sets the diffusivity to zero
       useSedimentation,          & ! true allows sediment accumulation or loss at the benthos surface
       useFastSedimentation,      & ! true enhances sedimentation rates to test advection
       usePositiveSedimentation,  & ! false switches from a sediment flux to resuspension at the benthic surface
       useBgcSinkingFlux,         & ! true uses ocean sinking flux of particulate benthos tracers
       useStepInitialProfiles,    & ! true sets the initial benthos profiles as step functions
       useBenthicReactions,       & ! turns on and off reaction bio-chemistry
       useFluxCorrection,         & ! true tracks transport errors in the storage flux
       useOceanConc,              & ! false sets ocean concentrations of solutes to zero
       useDeepSource,             & ! false sets to zero all deep storage fluxes
       useConstantDiffusivity       ! true removes the depth dependence of the diffusivity

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
  
  !NJ-TEST
  logical(KIND=benthos_log) :: test_flag
  integer(KIND=benthos_i4) :: iBenthicTracers, iSecondaryReactions, iElements, iBenthicVertLevels
  !NJ-END
!-----------------------------------------------------------------------
!EOC
  !NJ-TEST
!  write(*,*) 'test benthos_flags_init'
  !NJ-END
  
 SELECT CASE (setBenthosFlags)
    CASE(:-1)

       ! 'User Specified Test Case: Use Flags in Namelist'

       useCarbonateSaturation = useCarbonateSaturation
       useSecondaryReactions = useSecondaryReactions
       useDepthDependentPorosity = useDepthDependentPorosity
       useNonZeroDiffusivity = useNonZeroDiffusivity
       useSedimentation = useSedimentation
       usePositiveSedimentation = usePositiveSedimentation
       useFastSedimentation = useFastSedimentation
       useBgcSinkingFlux = useBgcSinkingFlux
       useStepInitialProfiles = useStepInitialProfiles
       useBenthicReactions = useBenthicReactions
       useFluxCorrection = useFluxCorrection
       useOceanConc = useOceanConc
       useDeepSource = useDeepSource
       useConstantDiffusivity = useConstantDiffusivity

    CASE(1)

       !  'Positive fast advection test case, fixed porosity'

       useDepthDependentPorosity = .false.
       useNonZeroDiffusivity = .false.
       useSedimentation = .true.
       usePositiveSedimentation = .true.
       useFastSedimentation = .true.
       useBgcSinkingFlux = .false.
       useStepInitialProfiles = .true.
       useBenthicReactions = .false.
       useFluxCorrection = .false.
       useOceanConc = .false.
       useDeepSource = .false.
       useConstantDiffusivity = .false.
       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.

    CASE(2)

       ! 'Negative (resuspension) fast advection test case'

       useDepthDependentPorosity = .false.
       useNonZeroDiffusivity = .false.
       useSedimentation = .true.
       usePositiveSedimentation = .false.
       useFastSedimentation = .true.
       useBgcSinkingFlux = .false.
       useStepInitialProfiles = .true.
       useBenthicReactions = .false.
       useFluxCorrection = .false.
       useOceanConc = .false.
       useDeepSource = .false.
       useConstantDiffusivity = .false.
       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.

    CASE(3)

       ! 'Positive fast advection with variable porosity test case'

       useDepthDependentPorosity = .true.
       useNonZeroDiffusivity = .false.
       useSedimentation = .true.
       usePositiveSedimentation = .true.
       useFastSedimentation = .true.
       useBgcSinkingFlux = .false.
       useStepInitialProfiles = .true.
       useBenthicReactions = .false.
       useFluxCorrection = .false.
       useOceanConc = .false.
       useDeepSource = .false.
       useConstantDiffusivity = .false.
       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.

    CASE(4)

       ! 'Positive fast advection with surface bio source and  variable porosity test case'

       useDepthDependentPorosity = .true.
       useNonZeroDiffusivity = .false.
       useSedimentation = .true.
       usePositiveSedimentation = .true.
       useFastSedimentation = .true.
       useBgcSinkingFlux = .true.
       useStepInitialProfiles = .true.
       useBenthicReactions = .false.
       useFluxCorrection = .false.
       useOceanConc = .false.
       useDeepSource = .false.
       useConstantDiffusivity = .false.
       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.

    CASE(5)

       ! 'Positive fast advection with surface Bio source + reaction and variable porosity test case'

       useDepthDependentPorosity = .true.
       useNonZeroDiffusivity = .false.
       useSedimentation = .true.
       usePositiveSedimentation = .true.
       useFastSedimentation = .true.
       useBgcSinkingFlux = .true.
       useStepInitialProfiles = .true.
       useBenthicReactions = .true.
       useFluxCorrection = .false.
       useOceanConc = .false.
       useDeepSource = .false.
       useConstantDiffusivity = .false.
       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.

    CASE(6)

        ! 'Constant Diffusion (no bottom fluxes and no solid surface fluxes balanced by Q_top and Q_bot)  test case'

       useDepthDependentPorosity = .false.
       useNonZeroDiffusivity = .true.
       useSedimentation = .false.
       usePositiveSedimentation = .true.
       useFastSedimentation = .true.
       useBgcSinkingFlux = .false.
       useStepInitialProfiles = .true.
       useBenthicReactions = .false.
       useFluxCorrection = .false.
       useOceanConc = .true.
       useDeepSource = .false.
       useConstantDiffusivity = .true.
       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.

   CASE(7)

      ! 'Constant Diffusion with variable porosity non-zero ocean Conc for solutes  test case'

       useDepthDependentPorosity = .true.
       useNonZeroDiffusivity = .true.
       useSedimentation = .false.
       usePositiveSedimentation = .true.
       useFastSedimentation = .true.
       useBgcSinkingFlux = .false.
       useStepInitialProfiles = .true.
       useBenthicReactions = .false.
       useFluxCorrection = .false.
       useOceanConc = .true.
       useDeepSource = .false.
       useConstantDiffusivity = .true.
       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.

    CASE(8)

       ! 'Variable Diffusion and variable porosity non-zero ocean Conc for solutes  test case'

       useDepthDependentPorosity = .true.
       useNonZeroDiffusivity = .true.
       useConstantDiffusivity = .false.
       useSedimentation = .false.
       usePositiveSedimentation = .true.
       useFastSedimentation = .true.
       useBgcSinkingFlux = .false.
       useStepInitialProfiles = .true.
       useBenthicReactions = .false.
       useFluxCorrection = .true.
       useOceanConc = .true.
       useDeepSource = .false.
       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.

    CASE(9)

       ! 'Variable Diffusion, variable porosity + velocities'

       useDepthDependentPorosity = .true.
       useNonZeroDiffusivity = .true.
       useConstantDiffusivity = .false.
       useSedimentation = .true.
       usePositiveSedimentation = .true.
       useFastSedimentation = .true.
       useBgcSinkingFlux = .true.
       useStepInitialProfiles = .true.
       useBenthicReactions = .false.
       useFluxCorrection = .true.
       useOceanConc = .true.
       useDeepSource = .false.
       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.

    CASE(10)

       ! 'Reactions only with variable porosity'

       useDepthDependentPorosity = .true.
       useNonZeroDiffusivity = .false.
       useConstantDiffusivity = .false.
       useSedimentation = .false.
       usePositiveSedimentation = .true.
       useFastSedimentation = .false.
       useBgcSinkingFlux = .false.
       useStepInitialProfiles = .false.
       useBenthicReactions = .true.
       useFluxCorrection = .false.
       useOceanConc = .false.
       useDeepSource = .false.
       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.

    CASE(11)

       ! 'Reactions only with variable porosity'

       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.
       useDepthDependentPorosity = .true.
       useNonZeroDiffusivity = .false.
       useConstantDiffusivity = .false.
       useSedimentation = .false.
       usePositiveSedimentation = .true.
       useFastSedimentation = .false.
       useBgcSinkingFlux = .false.
       useStepInitialProfiles = .false.
       useBenthicReactions = .true.
       useFluxCorrection = .false.
       useOceanConc = .false.
       useDeepSource = .false.

    CASE DEFAULT

       ! 'Default test case'

       useCarbonateSaturation = .true.
       useSecondaryReactions = .true.
       useDepthDependentPorosity = .true.
       useNonZeroDiffusivity = .true.
       useSedimentation = .true.
       usePositiveSedimentation = .true.
       useFastSedimentation = .false.
       useBgcSinkingFlux = .true.
       useStepInitialProfiles = .false.
       useBenthicReactions = .true.
       useFluxCorrection = .true.
       useOceanConc = .true.
       useDeepSource = .true.
       useConstantDiffusivity = .false.

 END SELECT

 end subroutine benthos_flags_init

!***********************************************************************
!BOP
! !IROUTINE: benthos_SourceSink
! !INTERFACE:

 subroutine benthos_SourceSink(benthos_input, benthos_forcing, &
                           benthos_output, benthos_diagnostic_fields, &
                           nBenthicVertLevels,numBenthicColumnsMax, &
                           dt,err)

! !DESCRIPTION:
!  Compute time derivatives for benthos tracers
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  type(benthos_forcing_type),     intent(in ) :: benthos_forcing

  integer (KIND=benthos_i4), intent(in) :: nBenthicVertLevels, numBenthicColumnsMax

  real (KIND=benthos_r8), intent(in) :: &
       dt

! !INPUT/OUTPUT PARAMETERS:

  integer (KIND=benthos_i4), intent(inout) :: err
  type(benthos_input_type),        intent(inout ) :: benthos_input   ! may want to change this to in only
  type(benthos_output_type), intent(inout) :: benthos_output
  type(benthos_diagnostics_type), intent(inout) :: benthos_diagnostic_fields

! !!OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

  integer (KIND=benthos_i4), parameter :: &
       nOMtype = 3, &
       column = 1

  logical(KIND=benthos_log) :: &
       l_stop, &
       tooLarge, &
       lcomp_co3_coeff_benthos, &
       lcalc_co2_terms

  real(KIND=benthos_r8) :: &
       dts,              & ! timestep during subcycling
       oceanBottomTemperature, & ! oC
       oceanBottomSalinity, & !  (g salt/kg seawater)
       oceanBottomSilicate, & !
       oceanBottomDepth, &
       oceanBottomDensity, &
       oceanBottomPhosphate, &
       oceanBottomAmmonium, &
       oceanBottomNitrate, &
       oceanBottomOxygen, &
       oceanBottomIron, &
       oceanBottomDIC, &
       oceanBottomAlkalinity, &
       oceanPOCFlux, &
       oceanPONFlux, &
       oceanPOPFlux, &
       oceanParticulateIronFlux, &
       oceanCalciteFlux

  real(KIND=benthos_r8) :: &
       work1, work2, work3, work4, work5, & !
       FluxSed,      & ! sediment sinking flux (kg/m2/s)
       dhtop,          & ! change in benthos upper surface (m)
       vel,              &  ! velocity of upper surface (m/s)
       C_topm,      & ! surface boundary concentration for solute or solid
       C_botm,      & ! bottom boundary concentration for solute or solid
       C_topSolid, & !
       sum_old, &
       sum_new, &
       sum_tot, &
       sum_init, &
       totalOceanSolidFluxm, & !
       benthosBottomFluxm, &
       subN, &
       dt_old, &
       del_ph, &
       phlo_3d_init, &
       phhi_3d_init, &
       reactionTend, &
       reactionTendprimary, &
       reactionTendsecondary, &
       reactionTendcarbonate, &
       trueTend, &
       tempR, &
       tempSol, &
       tempS, &
       rel_error_carbon,&
       benthosTracerBeforeReaction

  real(KIND=benthos_r8), dimension(:), allocatable :: &
       oceanSedFluxdt, & ! particulate flux * dt of sinking ocean tracers (mmol/kg m/s dt)
       oceanDiffVelFluxdt, & ! diffusive/vel surface boundary flux (mmol/kg or mmol/m3) * m/s dt
       benthosSedFluxdt, & ! particulate flux to deep storage (mmol/kg or mmol/m3) * m/s dt
       benthosDiffVelFluxdt, & ! diffusive/vel bottom boundary flux (mmol/kg or mmol/m3) * m/s * dt
       diffError, & ! mmol/m3 or mmol/kg
       totFluxesdt, & ! sum of tracer boundary fluxes * dt
       oceanSedFlux,    & ! particulate flux of sinking ocean tracers (mmol/kg*m/s)
       oceanTracerConc, & ! bottom ocean concentration of solutes (mmol/m3)
       oceanPrecipFlux, & ! flux of precipitated material from ocean bottom (mmol/kg*m/s)
       totalOceanSolidFlux, & ! total solid flux divided by benthosDepth (mmol/kg/s)
       deepBenthosBottomFlux, &  ! prescribed deep storage flux (mmol/m2/s)
       netElementsStart, & ! initial column element concentration (mmol/m2)
       totalColumnTracersInitial, &  ! initial column integrated tracer (mmol/m2)
       totalColumnTracersFinal, &  ! final column integrated tracer (mmol/m2)
       bDiff, & ! fluid or solid Diffusivity on the benthos grid
       bpor, & ! fluid or solid porosity on the benthos grid
       iDiff, & ! fluid or solid Diffusivity on the interfaces
       ipor, &  ! fluid or solid porosities on the interfaces
       C_top_vel, & ! (mmol/m3) surface concentration for velocity boundary condition
       C_bot_vel, & ! (mmol/m3) bottom concentration for velocity boundary condition
       C_top, & ! (mmol/m3) surface concntration for diffusive boundary condition
       C_bot, & ! (mmol/m3) bottom concentration for diffusive boundary condition
       source, & !
       bottomsource, &
       C_tot !, &
!       netTransport_tend

  real(KIND=benthos_r8), dimension(nBenthicTracers) :: &
       Source_top, &
       Source_bot, &
       Sink_bot_d, &
       Sink_bot_s, &
       Sink_top_d, &
       Sink_top_s, &
       flux_bio, &
       benthosStorageConc

  real(KIND=benthos_r8), dimension(nElements) :: &
       oceanElementExchange, &
       burialElementExchange, &
       netElements_end_Transport, &
       oceanElementSedimentation, &
       deepBurial, &
       netElements_start_reactions, &
       netElements_start_reactions_tmp, &
       netElements_end_reactions_tmp, &
       netElements_end_reactions

  real(KIND=benthos_r8), dimension(nCarbonateReactions) :: &
       carbonateRates

  real(KIND=benthos_r8), dimension(nSecondaryReactions) :: &
       secondaryRates

  real(KIND=benthos_r8), dimension(nOMtype,nPrimaryReactions) :: &
       primaryRates

  real(KIND=benthos_r8), dimension(nBenthicVertLevels+1) :: &
       initcons, & ! array for initial tracer concentration
       biocons, & ! array for tracer concentration
       biomat_low, &
       D_spdiag, D_sbdiag, ML_diag, &
       rhs, spdiag, diag, sbdiag

  real(KIND=benthos_r8), dimension(:,:), allocatable :: &
       benthosTracerBulk, &    ! local tracer concentration
       initBenthosTracerBulk, &
       primarySourceTend, &
       primarySinkTend, &
       carbonateSourceTend, &
       carbonateSinkTend, &
       secondarySourceTend, &
       secondarySinkTend, &
       netReactionTend

  integer (KIND=benthos_i4) :: &
       iLevels, &
       iTracers, &
       subtt

  real (KIND=benthos_r8) :: &
       co3_mol, &   ! comp_CO3terms
       hco3_mol, &
       co2_mol, &
       sat_calc, sat_arag, sat_mgcalc, &
       K_calc, K_arag, K_mgcalc, &
       phlo, &
       phhi, &
       ph, &
       kw, &
       kb, &
       ks, &
       kf, &
       k1p, &
       k2p, &
       k3p, &
       ksi, &
       bt, &
       st, &
       ft, &
       dic, &
       ta, &
       pt, &
       sit

!-----------------------------------------------------------------------
!EOC
!-----------------------------------------------------------------------
  write(*,*)'begin benthos_SourceSink'
  
  ! initialize arrays

   allocate(deepBenthosBottomFlux(nBenthicTracers))
   allocate(oceanSedFlux(nBenthicTracers))
   allocate(oceanPrecipFlux(nBenthicTracers))
   allocate(totalOceanSolidFlux(nBenthicTracers))
   allocate(oceanTracerConc(nBenthicTracers))
   allocate(netElementsStart(nElements))
   allocate(benthosTracerBulk(nBenthicTracers, nBenthicVertLevels+1))
   allocate(initBenthosTracerBulk(nBenthicTracers, nBenthicVertLevels+1))
   allocate(totalColumnTracersInitial(nBenthicTracers))
   allocate(totalColumnTracersFinal(nBenthicTracers))
   allocate(bDiff(nBenthicVertLevels+2))
   allocate(iDiff(nBenthicVertLevels+1))
   allocate(bpor(nBenthicVertLevels+2))
   allocate(ipor(nBenthicVertLevels+1))
   allocate(C_top_vel(nBenthicTracers))
   allocate(C_bot_vel(nBenthicTracers))
   allocate(C_top(nBenthicTracers))
   allocate(C_bot(nBenthicTracers))
   allocate(source(nBenthicTracers))
   allocate(bottomsource(nBenthicTracers))
   allocate(C_tot(nBenthicTracers))
!   allocate(netTransport_tend(nBenthicTracers))    
   allocate(oceanSedFluxdt(nBenthicTracers))
   allocate(oceanDiffVelFluxdt(nBenthicTracers))
   allocate(benthosSedFluxdt(nBenthicTracers))
   allocate(benthosDiffVelFluxdt(nBenthicTracers))
   allocate(diffError(nBenthicTracers))
   allocate(totFluxesdt(nBenthicTracers))
   allocate(primarySourceTend(nBenthicTracers,nBenthicVertLevels+1))
   allocate(primarySinkTend(nBenthicTracers,nBenthicVertLevels+1))
   allocate(carbonateSourceTend(nBenthicTracers,nBenthicVertLevels+1))
   allocate(carbonateSinkTend(nBenthicTracers,nBenthicVertLevels+1))
   allocate(secondarySourceTend(nBenthicTracers,nBenthicVertLevels+1))
   allocate(secondarySinkTend(nBenthicTracers,nBenthicVertLevels+1))
   allocate(netReactionTend(nBenthicTracers,nBenthicVertLevels+1))

   err = 0

   deepBenthosBottomFlux(:) = c0_benthos !   Bottom boundary condition  (mmol/m2/s)
   oceanSedFlux(:)          = c0_benthos
   oceanPrecipFlux (:)      = c0_benthos
   totalOceanSolidFlux(:)   = c0_benthos
   oceanTracerConc(:)       = c0_benthos

   !Initialize local variables
   oceanBottomDepth = benthos_input%oceanBottomDepth(column)
   oceanBottomTemperature = benthos_input%oceanBottomTemperature(column)
   oceanBottomSalinity = benthos_input%oceanBottomSalinity(column)
   oceanBottomSilicate = benthos_input%oceanBottomSilicate(column)
   oceanBottomDensity = benthos_input%oceanBottomDensity(column)
   oceanBottomPhosphate  = benthos_input%oceanBottomPhosphate(column)
   oceanBottomAmmonium  = benthos_input%oceanBottomAmmonium(column)
   oceanBottomNitrate  = benthos_input%oceanBottomNitrate(column)
   oceanBottomOxygen  = benthos_input%oceanBottomOxygen(column)
   oceanBottomIron  = benthos_input%oceanBottomIron(column)
   oceanBottomDIC  = benthos_input%oceanBottomDIC(column)
   oceanBottomAlkalinity  = benthos_input%oceanBottomAlkalinity(column)
   oceanPOCFlux  = benthos_input%oceanPOCFlux(column)
   oceanPONFlux  = benthos_input%oceanPONFlux(column)
   oceanPOPFlux  = benthos_input%oceanPOPFlux(column)
   oceanParticulateIronFlux  = benthos_input%oceanParticulateIronFlux(column)
   oceanCalciteFlux  = benthos_input%oceanCalciteFlux(column)

   write(*,*) 'column:',column
   write(*,*) 'oceanBottomDepth:', oceanBottomDepth
   write(*,*) 'oceanBottomTemperature:', oceanBottomTemperature
   write(*,*) 'oceanBottomSalinity:', oceanBottomSalinity
   write(*,*) 'oceanBottomSilicate:', oceanBottomSilicate
   write(*,*) 'oceanBottomDensity:', oceanBottomDensity
   write(*,*) 'oceanBottomPhosphate:', oceanBottomPhosphate
   write(*,*) 'oceanBottomAmmonium:', oceanBottomAmmonium
   write(*,*) 'oceanBottomNitrate:', oceanBottomNitrate
   write(*,*) 'oceanBottomOxygen:', oceanBottomOxygen
   write(*,*) 'oceanBottomIron:', oceanBottomIron
   write(*,*) 'oceanBottomDIC:', oceanBottomDIC
   write(*,*) 'oceanBottomAlkalinity:', oceanBottomAlkalinity
   write(*,*) 'oceanPOCFlux:', oceanPOCFlux
   write(*,*) 'oceanPONFlux:', oceanPONFlux
   write(*,*) 'oceanPOPFlux:', oceanPOPFlux
   write(*,*) 'oceanParticulateIronFlux:', oceanParticulateIronFlux
   write(*,*) 'oceanCalciteFlux:', oceanCalciteFlux

   do iTracers = 1, nBenthicTracers
      benthosStorageConc(iTracers) = benthos_input%deepStorage(column,iTracers)
      do iLevels = 1, nBenthicVertLevels+1
        benthosTracerBulk(iTracers,iLevels) = benthos_input%benthosTracerBulk(iLevels,column,iTracers)
      end do
   end do

   ! update time-dependent parameters (temperature dependent diffusion)
   call benthos_parameters_dt (nBenthicVertLevels, oceanBottomTemperature)
   
   ! Initialize the bottom flux boundary conditions
   call updateBenthosBoundaryConditions (oceanSedFlux, oceanPrecipFlux, totalOceanSolidFlux, &
        oceanTracerConc, deepBenthosBottomFlux, oceanBottomDensity, nBenthicVertLevels, &
        FluxSed, dhtop, vel, dt, oceanBottomPhosphate, oceanBottomAmmonium, oceanBottomNitrate, &
        oceanBottomOxygen, oceanBottomIron, oceanBottomDIC, oceanBottomAlkalinity, oceanPOCFlux, &
        oceanPONFlux, oceanPOPFlux, oceanParticulateIronFlux, oceanCalciteFlux)

!   write(*,*) 'Before computeNetElements'
   call computeNetElements (netElementsStart, nBenthicVertLevels, benthosTracerBulk)

   do iTracers = 1,nBenthicTracers
      call sum_benthos_column (totalColumnTracersInitial(iTracers),benthosTracerBulk(iTracers,:),nBenthicVertLevels)

      bDiff(nBenthicVertLevels+2) = benthosDiffusivity(nBenthicVertLevels+2)*&
           real(benthos_input%tracerType(iTracers),KIND=benthos_r8) +  &
           benthosSolidDiffusivity(nBenthicVertLevels+2)*(c1_benthos-&
           real(benthos_input%tracerType(iTracers),KIND=benthos_r8))

      bpor(nBenthicVertLevels+2) = benthosPorosity(nBenthicVertLevels+2)*&
           real(benthos_input%tracerType(iTracers),KIND=benthos_r8) +  &
           benthosSolidPorosity(nBenthicVertLevels+2)*(c1_benthos-&
           real(benthos_input%tracerType(iTracers),KIND=benthos_r8))

      do iLevels = 1,nBenthicVertLevels+1

         iDiff(iLevels) = benthosDiffusivityI(iLevels)*real(benthos_input%tracerType(iTracers),KIND=benthos_r8) + &
              benthosSolidDiffusivityI(iLevels)*(c1_benthos - real(benthos_input%tracerType(iTracers),KIND=benthos_r8))
         bDiff(iLevels) = benthosDiffusivity(iLevels)*real(benthos_input%tracerType(iTracers),KIND=benthos_r8) + &
              benthosSolidDiffusivity(iLevels)*(c1_benthos-real(benthos_input%tracerType(iTracers),KIND=benthos_r8))

         ipor(iLevels) = benthosPorosityI(iLevels)*real(benthos_input%tracerType(iTracers),KIND=benthos_r8) + &
              benthosSolidPorosityI(iLevels)*(c1_benthos-real(benthos_input%tracerType(iTracers),KIND=benthos_r8))
         bpor(iLevels) = benthosSolidPorosity(iLevels)*real(benthos_input%tracerType(iTracers),KIND=benthos_r8) + &
              benthosSolidPorosity(iLevels)*(c1_benthos-real(benthos_input%tracerType(iTracers),KIND=benthos_r8))

!          benthicTracerBulk(iTracers,iLevels) = initialBenthicTracerBuliLevels(iTracers,iLevels)
!          benthicTracerBrine(iTracers,iLevels) =   benthicTracerBulk(iTracers,iLevels)/ipor(iLevels)
!          initailBenthicTracerBrine(iTracers,iLevels) = benthicTracerBulk(iTracers,iLevels)/ipor(iLevels)

          initcons(iLevels) = benthosTracerBulk(iTracers,iLevels)
          biocons(iLevels) = benthosTracerBulk(iTracers,iLevels)
       end do

     C_top_vel(iTracers) = benthosTracerBulk(iTracers,1)
     C_bot_vel(iTracers) = benthosTracerBulk(iTracers,nBenthicVertLevels+1)
     C_top(iTracers) = oceanTracerConc(iTracers)
     C_bot(iTracers) = c0_benthos
     !if (vel > 0)
       C_bot(iTracers) = C_bot_vel(iTracers)
     !end

      C_topSolid = c0_benthos
      !C_topm = C_top(iTracers)

      C_botm = C_bot(iTracers)
      C_topm = C_top(iTracers)*real(benthos_input%tracerType(iTracers),KIND=benthos_r8) + C_top_vel(iTracers)*(c1_benthos-real(benthos_input%tracerType(iTracers),KIND=benthos_r8))

      totalOceanSolidFluxm = totalOceanSolidFlux(iTracers)
      benthosBottomFluxm = deepBenthosBottomFlux(iTracers)/benthosDepth  ! normalized   positive down
      source(iTracers) = totalOceanSolidFlux(iTracers)    ! from ocean model
      bottomsource(iTracers) = deepBenthosBottomFlux(iTracers)/benthosDepth   ! positive down  *Prescribed*

      write(*,*) 'iTracers:',iTracers
      write(*,*) 'Before compute_FCT_matrix, totalOceanSolidFluxm:', totalOceanSolidFluxm
      write(*,*) 'benthosBottomFluxm:', totalOceanSolidFluxm
      write(*,*) 'source(iTracers):', source(iTracers)
      write(*,*) 'bottomsource(iTracers):', bottomsource(iTracers)
      write(*,*) 'biocons:', biocons
      
      call compute_FCT_matrix_CN (sbdiag, diag, spdiag, rhs, ML_diag, D_sbdiag, D_spdiag, &
           Source_top(iTracers), Source_bot(iTracers), Sink_bot_d(iTracers), Sink_bot_s(iTracers), &
           Sink_top_d(iTracers), Sink_top_s(iTracers),  &
           biocons, dt, nBenthicVertLevels,  dhtop,  ipor,  iDiff, benthosDepth,    &
           totalOceanSolidFluxm, bpor,   C_topm, C_botm,   C_top_vel(iTracers), &
           C_bot_vel(iTracers),  benthosBottomFluxm,vel)

      call tridiag_solverz (biocons, nBenthicVertLevels+1, sbdiag,  diag, spdiag, rhs)

      ! goes negative at boundaries right away !!!
      write(*,*) 'after tridiag_solverz, biocons:', biocons

      call check_conservation_FCT (flux_bio(iTracers), l_stop, benthosStorageConc(iTracers), &
           oceanSedFluxdt(iTracers),oceanDiffVelFluxdt(iTracers),benthosSedFluxdt(iTracers),     &
           benthosDiffVelFluxdt(iTracers),diffError(iTracers),totFluxesdt(iTracers), initcons,             &
           biocons, Source_top(iTracers), Source_bot(iTracers), Sink_bot_d(iTracers),                   &
           Sink_bot_s(iTracers), Sink_top_d(iTracers), Sink_top_s(iTracers), dt,                              &
           nBenthicVertLevels, bottomsource(iTracers), vel, benthosDepth,        &
           iTracers,benthosBottomFluxm,totalOceanSolidFluxm)

      if (l_stop) then
         err = 1
         return
      end if
      do iLevels = 1, nBenthicVertLevels+1
         biomat_low(iLevels) = biocons(iLevels)
      end do

      write(*,*) 'before compute_FCT_corr'

      call compute_FCT_corr (initcons,  &
                                 biocons, dt, nBenthicVertLevels, &
                                 D_sbdiag, D_spdiag, ML_diag)

      do iLevels = 1, nBenthicVertLevels+1
         benthosTracerBulk(iTracers,iLevels) = biocons(iLevels)
      end do

      call sum_benthos_column &
          (sum_old,biomat_low,nBenthicVertLevels)

      call sum_benthos_column &
          (sum_new,biocons,nBenthicVertLevels)

      call sum_benthos_column &
          (sum_tot,benthosTracerBulk(iTracers,:),nBenthicVertLevels)

       call sum_benthos_column &
          (sum_init,initcons,nBenthicVertLevels)

       if (abs(sum_new-sum_old) .gt. puny*sum_old .or. &
                minval(biocons) .lt. c0_benthos   .or.  l_stop) then
                write(*,*) 'Benthic FCT tracer solution failed,iTracers', iTracers
                write(*,*) 'sum_new,sum_old:',sum_new,sum_old
                write(*,*) 'iTracers,biocons(:):',iTracers,biocons(:)
                write(*,*) 'biomat_low:',biomat_low
                write(*,*) 'iDiff(:):',iDiff(:)
                write(*,*)  'initcons(:):',initcons(:)
                write(*,*) 'benthosTracerBulk(iTracers,:):',benthosTracerBulk(iTracers,:)
                write(*,*) 'dhtop',dhtop
                write(*,*)'source(iTracers):', source(iTracers)
                write(*,*)'bottomsource(iTracers):', bottomsource(iTracers)
                !stop_label = 'zbgc FCT tracer solution failed'
            end if
            if (sum_new .gt. 1e12_benthos_r8 ) then
                write(*,*)'Benthic tracer blowing up after FCT,iTracers', iTracers
                write(*,*)'sum_new,sum_old:',sum_new,sum_old
                write(*,*)'iTracers,biocons(:):',iTracers,biocons(:)
                write(*,*)'biomat_low:',biomat_low
                write(*,*)'iDiff(:):',iDiff(:)
                write(*,*) 'initcons(:):',initcons(:)
                write(*,*) 'benthosTracerBulk(iTracers,:):',benthosTracerBulk(iTracers,:)
                write(*,*)'dhtop:',dhtop
                write(*,*)'source(iTracers):', source(iTracers)
                write(*,*)'bottomsource(iTracers):', bottomsource(iTracers)
                !stop_label = 'Benthic tracer blowing up'
            end if

            do iLevels = 1,nBenthicVertLevels+1

                if (benthosTracerBulk(iTracers,iLevels) .lt. puny .and. &
                   benthosTracerBulk(iTracers,iLevels) .ge. -puny*1000.0_benthos_r8) then
                    benthosTracerBulk(iTracers,iLevels) = c0_benthos
                end if

                benthos_output%benthosTransportTendencies(iLevels,column,iTracers) = &
                    benthosTracerBulk(iTracers,iLevels) - initcons(iLevels)

                !tracer_bulk(iTracers,k,tt+1) = benthosTracerBulk(iTracers,iLevels)
                !tracer_brine(iTracers,k,tt+1) = benthicTracerBrine(iTracers,iLevels)

            end do  ! iLevels

      call sum_benthos_column (C_tot(iTracers),benthosTracerBulk(iTracers,:),nBenthicVertLevels)
      call sum_benthos_column (benthos_diagnostic_fields % diag_netBenthosTransportTend(column,iTracers), &
         benthos_output%benthosTransportTendencies(:,column,iTracers),nBenthicVertLevels)

   end  do !iTracers)

   call computeElementFluxes (oceanElementExchange, burialElementExchange, oceanSedFluxdt, &
        oceanDiffVelFluxdt, benthosSedFluxdt, benthosDiffVelFluxdt,benthosDepth)

   call computeNetElements (netElements_end_Transport, nBenthicVertLevels, benthosTracerBulk)

   call netElementExchange (oceanElementSedimentation, deepBurial, totalOceanSolidFlux, &
        benthosStorageConc,dt,benthosDepth)

   !
   ! Reaction terms : all the tracers have been updated for
   !  advection and transport.  Now define the reaction tendencies
   !

   ! initialize
   dts = dt
   dt_old = dt
   subN = p1_benthos
   tooLarge = .true.

   !  make sure there is consistency to start:
   do iLevels = 1, nBenthicVertLevels+1

      benthosTracerBulk(dic_ind,iLevels) = benthosTracerBulk(hco3_ind,iLevels) + benthosTracerBulk(co2_ind,iLevels)  +  benthosTracerBulk(co3_ind,iLevels)
      benthosTracerBulk(fepa_ind,iLevels) = benthosTracerBulk(feoh3a_ind,iLevels)*ironBoundPFraction
      benthosTracerBulk(fepb_ind,iLevels) = benthosTracerBulk(feoh3b_ind,iLevels)*ironBoundPFraction

      initBenthosTracerBulk(iTracers,iLevels) = benthosTracerBulk(iTracers,iLevels)
   end do

   !
   !  Save initial mmols of each ion
   ! Count the total mols of the following ions:
   ! In the order: C, O, N, P, S, Mn, Fe
   !

   call computeNetElements (netElements_start_reactions,nBenthicVertLevels,benthosTracerBulk)

   if (useBenthicReactions) then

      do while (tooLarge)
         subN = NINT(subN*10.0_benthos_r8)
         dts = dts/subN
         tooLarge = .false.
         do iTracers = 1, nBenthicTracers
            do iLevels = 1, nBenthicVertLevels+1
               benthosTracerBulk(iTracers,iLevels) = initBenthosTracerBulk(iTracers,iLevels)
            end do
         end do

         do subtt = 1,subN
            do iLevels = 1,nBenthicVertLevels+1

               !if (iLevels == 1)
               !netCarbon_carbonate_err(tt) = 0
               !netco3_carbonate_tend(tt) = 0
               !netco2_carbonate_tend(tt) = 0
               !nethco3_carbonate_tend(tt) = 0
               !avgpH(tt) = 0
               !end

               call computeNetElements (netElements_start_reactions_tmp, nBenthicVertLevels, benthosTracerBulk)

               !  Primary reactions

               primarySourceTend(:,iLevels) = c0_benthos
               primarySinkTend(:,iLevels) = c0_benthos

               call primaryStoichMatrix (primarySourceStoich, primarySinkStoich, &
                    benthosTracerBulk(:,iLevels))

               call primaryRateConstants (primaryRates,benthosTracerBulk(:,iLevels),nOMtype)

               call primaryBenthosReactions (primarySourceTend(:,iLevels),primarySinkTend(:,iLevels),&
                    benthosTracerBulk(:,iLevels),primarySourceStoich,&
                    primarySinkStoich,primaryRates,dts,nOMtype)

               !-------------------------------
               !  Carbonate chemistry
               !-------------------------------

               carbonateSourceTend(:,iLevels) = c0_benthos
               carbonateSinkTend(:,iLevels) = c0_benthos

               ! equilibrium concentrations
               lcomp_co3_coeff_benthos = .true.
               lcalc_co2_terms = .true.
               del_ph = 0.2_benthos_r8
               phlo_3d_init = 6.0_benthos_r8
               phhi_3d_init = 9.0_benthos_r8

               !sumC_initial = benthicTracerBulk(dic_ind,iLevels) %benthicTracerBulk(co3_ind,iLevels) + benthicTracerBulk(co2_ind,iLevels) + benthicTracerBulk(hco3_ind,iLevels)
               !sumco3_initial = benthicTracerBulk(co3_ind,iLevels)
               !sumco2_initial = benthicTracerBulk(co2_ind,iLevels)
               !sumhco3_initial = benthicTracerBulk(hco3_ind,iLevels)

               if (useCarbonateSaturation) then

                  work3 = c0_benthos
                  work4 = c0_benthos

                  if (lcalc_co2_terms) then
                     if (benthos_input%PH_PREV_3D(iLevels,column) .ne. c0_benthos) then
                        work1 = benthos_input%PH_PREV_3D(iLevels,column) - del_ph
                        work2 = benthos_input%PH_PREV_3D(iLevels,column) + del_ph
                     else
                        work1 = phlo_3d_init
                        work2 = phhi_3d_init
                     end if

                     work5 = oceanBottomDepth + benthosMidPointI(iLevels)   ! m

                     call comp_CO3terms_benthos (co3_mol,hco3_mol,co2_mol,ph,kw,kb, &
                          ks,kf,k1p,k2p,k3p,ksi,bt,st,ft,dic,ta,pt,sit, oceanBottomTemperature,  &
                          oceanBottomSalinity,work5, &
                          lcomp_co3_coeff_benthos, oceanBottomSilicate, work1, work2, &
                          benthosTracerBulk(dic_ind,iLevels), benthosTracerBulk(alk_ind,iLevels), &
                          benthosTracerBulk(h2po4_ind,iLevels), iLevels, oceanBottomDensity)

                     work3 = ph   ! do i need this?
                     benthos_input % PH_PREV_3D(iLevels,column) = ph
                  else
                     co2_mol  = c0_benthos !H2CO3 = 0
                     hco3_mol = c0_benthos !HCO3  = 0
                     co3_mol  = c0_benthos !CO3   = 0
                     benthos_input % PH_PREV_3D(iLevels,column) = 8.0_benthos_r8
                  end if

                  benthosTracerBulk(alk_ind,iLevels) = ta
                  benthosTracerBulk(h2po4_ind,iLevels) = pt
                  oceanBottomSilicate = sit

                  !sumC_final = dic  %co3_mol + co2_mol + hco3_mol
                  !netCarbon_carbonate_err(tt) = netCarbon_carbonate_err(tt) + (sumC_final-sumC_initial)*zspace(iLevels)*benthosDepth
                  !sumco3_final = benthicTracerBulk(co3_ind,iLevels)
                  !sumco2_final = benthicTracerBulk(co2_ind,iLevels)
                  !sumhco3_final = benthicTracerBulk(hco3_ind,iLevels)


                  !netco3_carbonate_tend(tt) = netco3_carbonate_tend(tt) + (sumco3_final - sumco3_initial)*zspace(iLevels)*benthosDepth
                  !netco2_carbonate_tend(tt) = netco2_carbonate_tend(tt) + (sumco2_final - sumco2_initial)*zspace(iLevels)*benthosDepth
                  !nethco3_carbonate_tend(tt) = nethco3_carbonate_tend(tt) + (sumhco3_final - sumhco3_initial)*zspace(iLevels)*benthosDepth
                  !avgpH(tt) = avgpH(tt) + ph*zspace(iLevels)


                  ! check for conserved carbon
                  !if (abs(sumC_initial - sumC_final) > accuracy*max(sumC_final,sumC_initial))
                  !  'carbon not conserved in carbonate chemistry'
                  !  hco3_mol = max(sumC_initial-co3_mol-co2_mol,0.)
                  !end

                  benthosTracerBulk(co3_ind,iLevels) = co3_mol
                  benthosTracerBulk(co2_ind,iLevels) = co2_mol
                  benthosTracerBulk(hco3_ind,iLevels) = hco3_mol
                  benthosTracerBulk(dic_ind,iLevels) = dic

                  ! saturation state

                  call  comp_co3_sat_benthos(sat_calc, sat_arag, sat_mgcalc,&
                       K_arag,K_calc,K_mgcalc,oceanBottomTemperature,oceanBottomSalinity,&
                       oceanBottomDepth,oceanMagnesium,co3_mol,oceanBottomDensity)  ! oceanBottomDepth to work5?

                  call carbonateRateConstants (carbonateRates, &
                       benthosTracerBulk(:,iLevels),dts,sat_calc, sat_arag,&
                       sat_mgcalc,oceanBottomDensity)

                  call carbonateBenthosReactions (carbonateSourceTend(:,iLevels),carbonateSinkTend(:,iLevels),&
                      carbonateSourceStoich,carbonateSinkStoich,carbonateRates,dts)

               end if ! useCarbonateSaturation

               ! Secondary Reactions

               secondarySourceTend(:,iLevels) = c0_benthos
               secondarySinkTend(:,iLevels) = c0_benthos


               if (useSecondaryReactions) then

                  call secondaryRateConstants(secondaryRates, benthosTracerBulk(:,iLevels))

                  call secondaryBenthosReactions(secondarySourceTend(:,iLevels),secondarySinkTend(:,iLevels), &
                       secondarySourceStoich,secondarySinkStoich,secondaryRates,dts)

               end if

               do iTracers = 1,nBenthicTracers

                  reactionTend = c0_benthos
                  trueTend = c0_benthos
                  tempR = c0_benthos
                  tempSol = c0_benthos
                  tempS = c0_benthos

                  !ipor(iLevels) = iphi(iLevels)*real(benthos_input%tracerType(iTracers),KIND=benthos_r8) +  &
                  !      iphis(iLevels)*(c1_benthos-real(benthos_input%tracerType(iTracers),KIND=benthos_r8))

                  if (benthosTracerBulk(iTracers,iLevels) .ge. puny) then
                     tempR = (primarySinkTend(iTracers,iLevels) + secondarySinkTend(iTracers,iLevels)+&
                          carbonateSinkTend(iTracers,iLevels))/(benthosTracerBulk(iTracers,iLevels) )
                  else
                     tempR = c0_benthos
                  end if


                  tempS = (primarySourceTend(iTracers,iLevels) + secondarySourceTend(iTracers,iLevels)+carbonateSourceTend(iTracers,iLevels))


                  if (tempR .lt. puny) then
                     tempSol = benthosTracerBulk(iTracers,iLevels)  + tempS* (c1_benthos - tempR/c2_benthos)
                  else
                     tempSol = benthosTracerBulk(iTracers,iLevels)*exp(-tempR) + tempS*(c1_benthos-exp(-tempR))/tempR
                  end if

                  trueTend =  (tempSol-benthosTracerBulk(iTracers,iLevels))

                  reactionTendprimary = (primarySourceTend(iTracers,iLevels)-primarySinkTend(iTracers,iLevels))
                  reactionTendsecondary = (secondarySourceTend(iTracers,iLevels)-secondarySinkTend(iTracers,iLevels))
                  reactionTendcarbonate =  (carbonateSourceTend(iTracers,iLevels)-carbonateSinkTend(iTracers,iLevels))

                  benthos_output % benthosReactionTendencies(iLevels, column, iTracers) = trueTend

                  ! For testing  elemental conservation of reaction set at each iLevels

                  netReactionTend(iTracers,iLevels) = reactionTendsecondary + reactionTendprimary+ reactionTendcarbonate  ! netTend?


                  if (iTracers .eq. nBenthicTracers .and. iLevels .eq. nBenthicVertLevels+1) then
                     call computeNetElements (netElements_end_reactions_tmp,nBenthicVertLevels,benthosTracerBulk)
                     rel_error_carbon = (netElements_end_reactions_tmp(1) - netElements_start_reactions_tmp(1))/(netElements_end_reactions_tmp(1)+c1_benthos)
                  else
                     rel_error_carbon = c0_benthos
                  end if


                  if (iTracers .eq. nBenthicTracers .and. abs(rel_error_carbon) .gt. 1.0e-10_benthos_r8 .and. &
                       .not. tooLarge .and.  subN .lt. 100.0_benthos_r8) then
                     tooLarge = .true. ! subcycle
                     exit
                  end if

                  benthosTracerBeforeReaction = benthosTracerBulk(iTracers,iLevels)
                  benthosTracerBulk(iTracers,iLevels) = benthosTracerBulk(iTracers,iLevels) + trueTend

                  !  benthosTracerBrine(iTracers,iLevels)= benthosTracerBulk(iTracers,iLevels)/ipor(k)


                  if (benthosTracerBulk(iTracers,iLevels) .lt. c0_benthos .and. benthosTracerBulk(iTracers,iLevels) .ge. -puny) then

                     !   'Values Slightly Negative'
                     benthosTracerBulk(iTracers,iLevels) = c0_benthos
                  end if

                  if (benthosTracerBulk(iTracers,iLevels) .lt. -puny) then


                     !'Negative benthic tracers after  reactions'
                     !'tempR,tempS:',tempR,tempS
                     !       'tracer index (iTracers),vertical index (k), timestep(tt),subtt:',iTracers,k,tt,subtt
                     !      'benthosTracerBulk(iTracers,iLevels):',benthosTracerBulk(iTracers,iLevels)
                     !       'benthosTracerBulk(iTracers,iLevels) before reactions:',benthosTracerBeforeReaction
                     !       'reactionTendBenthicTracerBulk(iTracers,iLevels),trueTend,reactionTend:',reactionTendBenthicTracerBulk(iTracers,iLevels),trueTend,reactionTend
                     !      'reactionTendprimary, reactionTendsecondary:',reactionTendprimary,reactionTendsecondary
                     !      'ipor(k):',ipor(k)
                     err = 1

                  else if (benthosTracerBulk(iTracers,iLevels) .gt. 1.0e20_benthos_r8) then

                     !      'benthic tracers blowing up after  reactions'
                     !'tracer index (iTracers),vertical index (k), timestep(tt),subtt:',iTracers,k,tt,subtt
                     ! 'benthosTracerBulk(iTracers,iLevels):',benthosTracerBulk(iTracers,iLevels)
                     ! 'benthosTracerBulk(iTracers,iLevels) before reaction:',benthosTracerBeforeReaction
                     !'reactionTendBenthicTracerBulk(iTracers,iLevels),trueTend,reactionTend:',reactionTendBenthicTracerBulk(iTracers,iLevels),trueTend,reactionTend
                     !'ipor(k):',ipor(k)
                     err = 1

                  end if

                  if (subtt .eq. subN .and. iLevels .eq. nBenthicVertLevels+1) then

                     ! tracer_bulk(iTracers,k,tt+1) = benthosTracerBulk(iTracers,iLevels)
                     call sum_benthos_column(benthos_diagnostic_fields % diag_netBenthosReactionTend(column,iTracers), &
                          benthos_output % benthosReactionTendencies(:, column, iTracers),nBenthicVertLevels)

                  end if ! if
               end do ! iTracers
            end do  ! iLevels
         end do !subtt  subgrid timestep
      end do !  while

   end if ! if useBenthicReactions

   call computeNetElements (netElements_end_reactions, nBenthicVertLevels, benthosTracerBulk)
   
   benthos_diagnostic_fields % diag_netBenthosTend(column,:) = c0_benthos
   
   ! New Solution, but Maybe needs to be done in integration routine?
   ! define diagnostic fields

   do iTracers = 1, nElements
      benthos_diagnostic_fields % diag_netElements(column,iTracers) = &
           netElements_end_reactions(iTracers)
   end do
   
   do iTracers = 1, nBenthicTracers

      benthos_diagnostic_fields % diag_netBenthosTend(column,iTracers) = &
         benthos_diagnostic_fields % diag_netBenthosTransportTend(column,iTracers) + &
         benthos_diagnostic_fields % diag_netBenthosReactionTend(column,iTracers)
      
      benthos_input%deepStorage(column,iTracers) = benthosStorageConc(iTracers)
      do iLevels = 1, nBenthicVertLevels+1
         benthos_input%benthosTracerBulk(iLevels,column,iTracers) = benthosTracerBulk(iTracers,iLevels)
      end do

   end do

   !
   ! call compute conservation and error diagnostics
   !

   deallocate(deepBenthosBottomFlux)
   deallocate(oceanSedFlux)
   deallocate(oceanPrecipFlux)
   deallocate(totalOceanSolidFlux)
   deallocate(oceanTracerConc)
   deallocate(netElementsStart)
   deallocate(benthosTracerBulk)
   deallocate(initBenthosTracerBulk)
   deallocate(totalColumnTracersInitial)
   deallocate(totalColumnTracersFinal)
   deallocate(bDiff)
   deallocate(iDiff)
   deallocate(bpor)
   deallocate(ipor)
   deallocate(C_top_vel)
   deallocate(C_bot_vel)
   deallocate(C_top)
   deallocate(C_bot)
   deallocate(source)
   deallocate(bottomsource)
   deallocate(C_tot)
!   deallocate(netTransport_tend)
   deallocate(oceanSedFluxdt)
   deallocate(oceanDiffVelFluxdt)
   deallocate(benthosSedFluxdt)
   deallocate(benthosDiffVelFluxdt)
   deallocate(diffError)
   deallocate(totFluxesdt)
   deallocate(primarySourceTend)
   deallocate(primarySinkTend)
   deallocate(carbonateSourceTend)
   deallocate(carbonateSinkTend)
   deallocate(secondarySourceTend)
   deallocate(secondarySinkTend)

 end subroutine benthos_SourceSink

!*****************************************************************************
!BOP
! !IROUTINE: compute_FCT_matrix_CN
! !INTERFACE:

 subroutine compute_FCT_matrix_CN (sbdiag, diag, spdiag, rhs, ML_diag, D_sbdiag, D_spdiag, &
      Qtop_out, Qbot_out, Sink_bot_diag,Sink_bot_sbdiag, Sink_top_diag,Sink_top_spdiag, &
      bioncons,  dt,  nBenthicVertLevels,   dhtop, &
      ipor, iDiff, benthosDepth,  totalOceanSolidFluxm, bpor,  &
      C_top, C_bot, C_top_vel, C_bot_vel,benthosBottomFluxm,vel)

! !DESCRIPTION:
!  compute the matrix elements for the flux-corrected conservative transport scheme
! Crank-nicholson version:
! [Mjk - dt * gamma *(Kjk + Sjk)]C'k = Mjk Ck + dt qk + (1-gamma)*(Kjk + Sjk)Ck
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer(KIND=benthos_i4), intent(in) :: nBenthicVertLevels

   real(KIND=benthos_r8), dimension(:), intent(in) :: &
        bpor,  & ! grid fluid or solid porosity
        ipor,  & ! interface fluid or solid porosity
        iDiff, & !  (m/s) normalized interface diffusivity
        bioncons ! mmol/m3 or mmol/kg tracer concentration

   real(KIND=benthos_r8), intent(in) :: &
        vel,                               & ! per s (normalized rate of surface growth or loss)
        benthosBottomFluxm, & ! (mmol/m3/s)  normalized positive down prescribed bottom flux
        C_bot_vel,                   & ! bottom concentration for velocity fluxes (mmol/m3 or mmol/kg)
        C_top_vel,                   & ! top concentration for velocity fluxes (mmol/m3 or mmol/kg)
        C_bot,                          & ! bottom concentration for diffusion(mmol/m3)
        C_top,                          & ! top concentration for diffusion (mmol/m3)
        totalOceanSolidFluxm, & ! precipitation and sinking flux of solids
        benthosDepth,              & ! (m) depth of the active benthos layer
        dhtop,                            & ! (m) change in the surface boundary in dt
        dt                                      ! (s) time-step

  ! !INPUT/OUTPUT PARAMETERS:

  ! !OUTPUT PARAMETERS:

  real(KIND=benthos_r8), intent(out) :: &
       Sink_top_diag,   & ! (m/s * dt) diagonal matrix contribution to surface flux
       Sink_top_spdiag, & ! (m/s * dt) off-diagonal matrix contribution to surface flux
       Sink_bot_sbdiag, & ! (m/s * dt) off-diagonal matrix contribution to bottom flux
       Sink_bot_diag,   & ! (m/s * dt) diagonal matrix contribution to bottom flux
       Qbot_out,        & ! (mmol/m2/s * dt) rhs contribution to bottom flux
       Qtop_out           ! (mmol/m2/s * dt) ths contribution to surface flux

  real(KIND=benthos_r8), dimension(nBenthicVertLevels+1), intent(out) :: &
       D_spdiag,             & !  matrix elements
       D_sbdiag,             & !
       ML_diag,              & !
       rhs,                  & ! right hand side
       spdiag,               & !
       diag,                 &
       sbdiag

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
  !-----------------------------------------------------------------------------
  real(KIND=benthos_r8) :: &
       zspace, &
       Sink_bot, &
       Sink_top, &
       dphi_dx, &
       dphi_dx1, &
       dphi_dxN, &
       gamma

  real (KIND=benthos_r8), dimension(nBenthicVertLevels+1) :: &
       ML, & ! lumped mass matrix
       K_diag, & !
       D_diag, & !
       K_spdiag, & !
       K_sbdiag, & !
       S_diag, &
       S_spdiag, &
       S_sbdiag, &
       Q_bot, &
       Q_top, &
       diagTemp, &
       ML_diag_CN

  integer (KIND=benthos_i4) :: &
       iLevels
!---------------------------------------------------------------------
!  Diag (jj) solve for j = 1:nBenthicVertLevels+1
!  spdiag(j) == (j,j+1) solve for j = 1:nBenthicVertLevels otherwise 0
!  sbdiag(j) == (j,j-1) solve for j = 2:nBenthicVertLevels+1 otherwise 0
!---------------------------------------------------------------------
     zspace = c1_benthos/real(nBenthicVertLevels, benthos_r8)
     Qbot_out = c0_benthos
     Qtop_out = c0_benthos
     Sink_bot = c0_benthos
     Sink_top = c0_benthos

    ! compute the lumped mass matrix

     ML_diag(1) = zspace/c2_benthos
     ML_diag(nBenthicVertLevels+1) = zspace/c2_benthos

    ! compute matrix K: K_diag , K_sbdiag, K_spdiag
    ! compute matrix S: S_diag, S_sbdiag, S_spdiag

    do iLevels = 1 , nBenthicVertLevels+1

       if (iLevels == 1) then

          dphi_dx = (ipor(iLevels+1)-bpor(iLevels))/(c2_benthos*zspace)
          dphi_dx = c0_benthos
          dphi_dx1 = dphi_dx

          K_sbdiag(iLevels)= c0_benthos
          K_spdiag(iLevels)= -p5_benthos*(min(c0_benthos,vel) +  &
               iDiff(iLevels)/(ipor(iLevels)*ipor(iLevels+1))*dphi_dx)
          !K_diag(k) =   -max(0,vel)- (iDiff(k))/(zspace*ipor(k)) ...
          ! + 0.5*(vel + iDiff(k)/ipor(k)^2*dphi_dx) !

          ! original
          K_diag(iLevels) =   -max(c0_benthos,vel)- (iDiff(iLevels))/(zspace*ipor(iLevels)) &
              + p5_benthos*(vel + iDiff(iLevels)/ipor(iLevels)**2*dphi_dx)

          S_sbdiag(iLevels) = c0_benthos
          S_diag(iLevels) = -(iDiff(iLevels))/(ipor(iLevels)*zspace)
          S_spdiag(iLevels) =  (iDiff(iLevels)+iDiff(iLevels+1))*p5_benthos/(ipor(iLevels+1)*zspace)

       else if (iLevels == nBenthicVertLevels+1) then

          dphi_dx = (bpor(iLevels+1)-ipor(iLevels-1))/(c2_benthos*zspace)
          dphi_dx = c0_benthos
          dphi_dxN = dphi_dx

          K_spdiag(iLevels) = c0_benthos

          K_sbdiag(iLevels)= p5_benthos*(max(c0_benthos,vel) +  &
               iDiff(iLevels)/(ipor(iLevels)*ipor(iLevels-1))*dphi_dx)
          K_diag(iLevels) =   - (iDiff(iLevels)) / (zspace*ipor(iLevels)) - &
               p5_benthos*(vel + iDiff(iLevels)/ipor(iLevels)**2*dphi_dx )+min(c0_benthos,vel)

          S_spdiag(iLevels) = c0_benthos
          S_sbdiag(iLevels) =  (iDiff(iLevels)+iDiff(iLevels-1))*p5_benthos/(ipor(iLevels-1)*zspace)
          S_diag(iLevels) = -(iDiff(iLevels))/(ipor(iLevels)*zspace)

       else

          dphi_dx = (ipor(iLevels+1)-ipor(iLevels-1))/(c2_benthos*zspace)
          dphi_dx = c0_benthos

          K_sbdiag(iLevels)= p5_benthos*(max(c0_benthos,vel))  + &
               p5_benthos*( iDiff(iLevels)/(ipor(iLevels)*ipor(iLevels-1))*dphi_dx)
          K_spdiag(iLevels)= -p5_benthos*(min(c0_benthos,vel)) - &
               p5_benthos*( iDiff(iLevels)/(ipor(iLevels)*ipor(iLevels+1))*dphi_dx)
          K_diag(iLevels) = -K_sbdiag(iLevels)-K_spdiag(iLevels)

          S_diag(iLevels) =    -(c2_benthos*iDiff(iLevels))/(ipor(iLevels)*zspace)
          S_sbdiag(iLevels)   =(iDiff(iLevels)+iDiff(iLevels-1))*p5_benthos/(ipor(iLevels-1)*zspace)
          S_spdiag(iLevels) = (iDiff(iLevels)+iDiff(iLevels+1))*p5_benthos/(ipor(iLevels+1)*zspace)
       end if
    end do

    ! compute matrix artificial D: D_spdiag, D_diag  (D_spdiag(iLevels) = D_sbdiag(iLevels+1))

     do iLevels = 1,nBenthicVertLevels
        D_spdiag(iLevels)    = max(c0_benthos,max(-K_spdiag(iLevels), -K_sbdiag(iLevels+1)))
        D_sbdiag(iLevels+1)  = D_spdiag(iLevels)
      end do
      do  iLevels = 1,nBenthicVertLevels+1
         D_diag(iLevels) =  D_diag(iLevels) - D_spdiag(iLevels) - D_sbdiag(iLevels)
      end do

      ! compute Q_top and Q_bot: top and bottom sources

      Q_top(:) = c0_benthos
      Q_top(1) =  (iDiff(1))*C_top/(zspace*bpor(1))+ &
           p5_benthos*max(c0_benthos,totalOceanSolidFluxm)

      Qtop_out = Q_top(1)

      Q_bot(:) = 0
      Q_bot(nBenthicVertLevels+1) =  (iDiff(nBenthicVertLevels+1))*C_bot &
           /(zspace*bpor(nBenthicVertLevels+2)) - &
           min(c0_benthos,benthosBottomFluxm)*p5_benthos- min(c0_benthos,vel)*C_bot

       Qbot_out = Q_bot(nBenthicVertLevels+1)

        if (vel .LE. c0_benthos) then

           Sink_top_diag = K_diag(1)
           Sink_top_spdiag = -p5_benthos*(iDiff(1)/(ipor(1)*ipor(2))*dphi_dx1)

           Sink_bot_diag =   - (iDiff(nBenthicVertLevels+1)) / &
                (zspace*ipor(nBenthicVertLevels+1)) - &
                p5_benthos*(iDiff(nBenthicVertLevels+1)/ipor(nBenthicVertLevels+1)**2*dphi_dxN )
           Sink_bot_sbdiag = K_sbdiag(nBenthicVertLevels+1)

        else   ! velocity advects downward away from the surface boundary

           Sink_top_diag = -(iDiff(1))/(zspace*ipor(1))  + p5_benthos*(iDiff(1)/ipor(1)**2*dphi_dx1)
           Sink_top_spdiag = K_spdiag(1)

           Sink_bot_diag  = K_diag(nBenthicVertLevels+1)

           Sink_bot_sbdiag =  p5_benthos*( iDiff(nBenthicVertLevels+1)/(ipor(nBenthicVertLevels+1)&
                *ipor(nBenthicVertLevels))*dphi_dxN)

        end if
        do iLevels = 1, nBenthicVertLevels+1
           spdiag(iLevels) = -dt *(D_spdiag(iLevels) + K_spdiag(iLevels) + S_spdiag(iLevels))
           sbdiag(iLevels) = -dt * (D_sbdiag(iLevels) + K_sbdiag(iLevels) + S_sbdiag(iLevels))
           diagTemp(iLevels) = -dt *(D_diag(iLevels) + K_diag(iLevels) + S_diag(iLevels))
           diag = ML_diag(iLevels) - dt * (D_diag(iLevels) + K_diag(iLevels) + S_diag(iLevels))
           ML_diag_CN((iLevels)) = c0_benthos
        end do

        ML_diag_CN(1) = spdiag(1) * bioncons(2) +diagTemp(1)* bioncons(1)
        ML_diag_CN(nBenthicVertLevels+1) = sbdiag(nBenthicVertLevels+1) * &
             bioncons(nBenthicVertLevels)  + diagTemp(nBenthicVertLevels+1) * &
             bioncons(nBenthicVertLevels+1)

        do  iLevels = 2,nBenthicVertLevels
           ML_diag_CN(iLevels) = -spdiag(iLevels) * bioncons(iLevels+1) - &
                sbdiag(iLevels) * bioncons(iLevels-1) - diagTemp(iLevels)*bioncons(iLevels)
        end do

       !
       ! Crank-Nicholson Formulation
        !

       gamma = c1_benthos
       do iLevels = 1, nBenthicVertLevels+1
          spdiag(iLevels) = -dt *gamma*(D_spdiag(iLevels) + K_spdiag(iLevels) + S_spdiag(iLevels))
          sbdiag(iLevels) = -dt * gamma*(D_sbdiag(iLevels) + K_sbdiag(iLevels) + S_sbdiag(iLevels))
          diag(iLevels) = ML_diag(iLevels) - dt *gamma* (D_diag(iLevels) + K_diag(iLevels) + &
               S_diag(iLevels))

          rhs(iLevels) = ML_diag(iLevels) * bioncons(iLevels) + dt * Q_top(iLevels) + &
               dt* Q_bot(iLevels) + (c1_benthos-gamma)*ML_diag_CN(iLevels)*dt
       end do

     end subroutine compute_FCT_matrix_CN

!*****************************************************************************
!BOP
! !IROUTINE: tridiag_solverz
! !INTERFACE:

subroutine tridiag_solverz(xout, nmat, sbdiagz, diagz, spdiagz, rhsz)

! !DESCRIPTION:
!  tridiagonal matrix solver
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer(KIND=benthos_i4), intent(in) :: nmat

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       diagz,     &
       spdiagz, &
       sbdiagz, &
       rhsz

! !INPUT/OUTPUT PARAMETERS:

  real(KIND=benthos_r8),dimension(:), intent(inout) :: xout

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------
  integer(KIND=benthos_i4) :: k

  real(KIND=benthos_r8) :: wbeta

  real(KIND=benthos_r8), dimension(nmat) :: &
       wgamma

   wbeta = diagz(1)
   xout(1) = rhsz(1) / wbeta

   do k = 2, nmat
      wgamma(k) = spdiagz(k-1) / wbeta
      wbeta = diagz(k) - sbdiagz(k)*wgamma(k)
      xout(k) = (rhsz(k) - sbdiagz(k)*xout(k-1)) / wbeta
   end do

    do k = nmat-1, 1, -1
         xout(k) = xout(k) - wgamma(k+1)*xout(k+1)
    end do

  end subroutine tridiag_solverz

!*****************************************************************************
!BOP
! !IROUTINE: check_conservation_FCT
! !INTERFACE:

  subroutine check_conservation_FCT(fluxbio, l_stop,deep_storage,oceanSedFluxdt, &
       oceanDiffVelFluxdt, benthosSedFluxdt, benthosDiffVelFluxdt, diffError, sources, &
       C_init, C_new, S_top,  S_bot, L_bot_d, L_bot_s,L_top_d,L_top_s, dt,     &
       nBenthicVertLevels, bottomsource, vel,benthosDepth, iTracers, &
       benthosBottomFluxm, totalOceanSolidFluxm)

! !DESCRIPTION:
!  Initialize tracegas tracer module.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer(KIND=benthos_i4), intent(in) :: &
       nBenthicVertLevels, &
       iTracers

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       C_init, &
       C_new

  real(KIND=benthos_r8), intent(in) :: &
       S_top, &
       S_bot, &
       L_bot_d, &
       L_bot_s, &
       L_top_s, &
       L_top_d, &
       dt, &
       bottomsource, &
       vel, &
       benthosDepth, &
       benthosBottomFluxm, &
       totalOceanSolidFluxm

  ! !INPUT/OUTPUT PARAMETERS:

  real(KIND=benthos_r8), intent(inout) :: &
       deep_storage, & !
       oceanSedFluxdt, & !
       benthosSedFluxdt, &
       sources

  ! !OUTPUT PARAMETERS:

  logical(KIND=benthos_log), intent(out) :: l_stop

  real(KIND=benthos_r8), intent(out) :: &
       fluxbio, & !
       oceanDiffVelFluxdt, & !
       benthosDiffVelFluxdt, &
       diffError

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------
  real(KIND=benthos_r8), dimension(nBenthicVertLevels+1) :: &
       dC

  real(KIND=benthos_r8) :: &
       C_init_tot, &
       C_new_tot, &
       dC_tot, &
       sources_o, &
       sources_top, &
       sources_bot, &
       deep_storage_tmp, &
       diffCnewCinit, &
       diffError_o, &
       sources_top_corr, &
       sources_bot_corr, &
       sources_corr, &
       diffError_corr, &
       accuracy

  integer(KIND=benthos_i4) :: &
       iLevels

  l_stop = .false.

!-------------------------------------
!  Ocean flux: positive into the ocean
!-------------------------------------

   dC(:) = c0_benthos
   do iLevels = 1, nBenthicVertLevels + 1
      dC(iLevels) = C_new(iLevels) - C_init(iLevels)
   end do

   call sum_benthos_column (C_init_tot, C_init, nBenthicVertLevels)
   call sum_benthos_column (C_new_tot, C_new,nBenthicVertLevels)
   call sum_benthos_column (dC_tot, dC,nBenthicVertLevels)

    ! This is the estimate of the source and source from the top and bottom
    !       sources = (S_top+S_bot+L_bot_d*C_new(nBenthicVertLevels+1)+
    !  L_top_d*C_new(1)+L_bot_s*C_new(nBenthicVertLevels)+L_top_s*C_new(2))*dt*benthosDepth

   sources_o = (S_top+S_bot+L_bot_d*C_init(nBenthicVertLevels+1)+L_top_d*C_init(1)+ &
        L_bot_s*C_init(nBenthicVertLevels)+L_top_s*C_init(2))*dt*benthosDepth

   sources_top = (S_top+L_top_d*C_new(1)+L_top_s*C_new(2))*dt*benthosDepth

   sources_bot = (S_bot+L_bot_d*C_new(nBenthicVertLevels+1)+ &
        L_bot_s*C_new(nBenthicVertLevels))*dt*benthosDepth   ! mmol/kg * m or mmol/m2
              ! positive implies positive increase in benthic  so flux into the Benthic from below

   sources = sources_top + sources_bot

   deep_storage_tmp = deep_storage-sources_bot
   oceanSedFluxdt = totalOceanSolidFluxm*dt*benthosDepth
   benthosSedFluxdt = benthosBottomFluxm*dt*benthosDepth
   oceanDiffVelFluxdt = sources_top - oceanSedFluxdt
   benthosDiffVelFluxdt = -(sources_bot + benthosSedFluxdt)


   diffCnewCinit = (C_new_tot-C_init_tot)      ! mmol/m2
   diffError = (diffCnewCinit - sources)
   diffError_o = (diffCnewCinit-sources_o)    ! estimate of sources from old C_init


   if (useFluxCorrection) then
      benthosDiffVelFluxdt = benthosDiffVelFluxdt - diffError
      deep_storage = deep_storage_tmp - diffError
   else
      deep_storage = deep_storage_tmp
   end if

   sources_top_corr = oceanDiffVelFluxdt + oceanSedFluxdt
   sources_bot_corr = -benthosDiffVelFluxdt - benthosSedFluxdt
   sources_corr = sources_top_corr + sources_bot_corr
   diffError_corr = (diffCnewCinit - sources_corr)

   fluxbio = oceanDiffVelFluxdt/dt    !passed to the ocean bottom layers

! implicit time integration,  errors are linear in dt  and decrease with the grid spacing zspace*h. 


   accuracy = max(puny*max(maxval(C_new),maxval(C_init)),1000.0_benthos_r8*puny)

    if (minval(C_new) < c0_benthos) then
         write(*,*)'Positivity of zbgc low order solution failed: C_new:',C_new
         l_stop = .true.
     end if

     if (abs(diffError_corr) > accuracy) then ! & iTracers > 9)
            l_stop = .true.
            write(*,*) 'Conservation of zbgc low order solution failed: diffError (mmol/m2):', diffError
            write(*,*) 'diffError_o', diffError_o
            write(*,*) 'diffError_corr', diffError_corr
            write(*,*)'accuracy:',accuracy
            write(*,*)'iTracers',iTracers
            write(*,*) 'sources', sources
            write(*,*) 'sources_o',sources_o
            write(*,*)  'sources_corr', sources_corr
            write(*,*) 'sources_bot:',sources_bot
            write(*,*) 'sources_bot_corr:',sources_bot_corr
            write(*,*) 'sources_bot:',sources_top
            write(*,*) 'sources_bot_corr:',sources_top_corr
           write(*,*) 'Start deep_storage',deep_storage
           write(*,*) 'intermediate deep_storage',deep_storage_tmp
           write(*,*) 'Final deep_storage',deep_storage
           write(*,*)   'diffCnewCinit (mmol/m2):',diffCnewCinit
           write(*,*)   'oceanSedFluxdt',oceanSedFluxdt
           write(*,*)   'oceanDiffVelFluxdt',oceanDiffVelFluxdt
           write(*,*)  'totalOceanSolidFluxm,',totalOceanSolidFluxm
           write(*,*)  'benthosSedFluxdt',benthosSedFluxdt
           write(*,*)   'S_bot', S_bot
           write(*,*)   'L_bot_d', L_bot_d
           write(*,*)    'L_bot_s',L_bot_s
           write(*,*)  'fluxbio:',fluxbio
           write(*,*) 'bottom final tracer', C_new(nBenthicVertLevels+1)
           write(*,*) 'top final tracer', C_new(1)
           write(*,*) 'bottom init tracer', C_init(nBenthicVertLevels+1)
           write(*,*) 'top init tracer', C_init(1)
           write(*,*)   'vel:',vel
           write(*,*)   'S_top:',S_top
           write(*,*) 'Total initial tracer', C_init_tot
           write(*,*) 'Total final1  tracer', C_new_tot
     end if

   end subroutine check_conservation_FCT

!*****************************************************************************
!BOP
! !IROUTINE: computeElementFluxes
! !INTERFACE:

   subroutine computeElementFluxes (oceanElementExchange, burialElementExchange, oceanSedFluxdt, &
        oceanDiffVelFluxdt, benthosSedFluxdt, benthosDiffVelFluxdt,benthosDepth)

! !DESCRIPTION:
!  calculates net boundary exchange fluxes for each of 7 elements: C, O, N, P, S, Mn, Fe
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  real(KIND=benthos_r8), intent(in) :: &
       benthosDepth

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       benthosDiffVelFluxdt, &
       benthosSedFluxdt, &
       oceanDiffVelFluxdt, &
       oceanSedFluxdt
  
! !INPUT/OUTPUT PARAMETERS:

! !OUTPUT PARAMETERS:
  real(KIND=benthos_r8), dimension(nElements), intent(out) :: &
       burialElementExchange, &
       oceanElementExchange

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------
   real (KIND=benthos_r8), dimension(nBenthicTracers) :: &
       topExchange, & !
       botExchange

   integer (KIND=benthos_i4) :: &
        iTracers, &
        iElements

   do iTracers = 1,nBenthicTracers
      topExchange(iTracers) = oceanDiffVelFluxdt(iTracers) + oceanSedFluxdt(iTracers)      ! positive down
      botExchange(iTracers) = benthosSedFluxdt(iTracers) + benthosDiffVelFluxdt(iTracers)   ! positive down
   end do

   do iElements = 1,nElements
      burialElementExchange(iElements) = c0_benthos
      oceanElementExchange(iElements) = c0_benthos
      do iTracers = 1,nBenthicTracers
         oceanElementExchange(iElements) = oceanElementExchange(iElements) + &
            elementRatios(iTracers,iElements)*topExchange(iTracers)
         burialElementExchange(iElements) = burialElementExchange(iElements) + &
            elementRatios(iTracers,iElements)*botExchange(iTracers)
      end do
   end do

   end subroutine computeElementFluxes

!*****************************************************************************
!BOP
! !IROUTINE: netElementExchange
! !INTERFACE:

   subroutine netElementExchange (oceanElementSedimentation, deepBurial, &
        totalOceanSolidFlux,benthosStorageConc,dt,benthosDepth)

! !DESCRIPTION:
!  calculates net boundary exchange fluxes for each of 7 elements: C, O, N, P, S, Mn, Fe
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  real(KIND=benthos_r8), intent(in) :: &
       benthosDepth, &
       dt

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       benthosStorageConc, &
       totalOceanSolidFlux
  
! !INPUT/OUTPUT PARAMETERS:

! !OUTPUT PARAMETERS:
  real(KIND=benthos_r8), dimension(nElements), intent(out) :: &
       oceanElementSedimentation, &
       deepBurial

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

   integer (KIND=benthos_i4) :: &
        iTracers, &
        iElements

   do iElements = 1,nElements
      do iTracers = 1,nBenthicTracers
         oceanElementSedimentation(iElements) = oceanElementSedimentation(iElements) + &
              elementRatios(iTracers,iElements)*totalOceanSolidFlux(iTracers)*benthosDepth*dt
         deepBurial(iElements) = deepBurial(iElements) + elementRatios(iTracers,iElements)*benthosStorageConc(iTracers)
      end do
   end do

   end subroutine netElementExchange

!*****************************************************************************
!BOP
! !IROUTINE: compute_FCT_corr
! !INTERFACE:

subroutine compute_FCT_corr (C_in,  &
                                 C_low, dt, nBenthicVertLevels, &
                                 D_sbdiag, D_spdiag, ML)

! !DESCRIPTION:
!  computes the correction to the initial solver
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer(KIND=benthos_i4), intent(in) :: nBenthicVertLevels

  real(KIND=benthos_r8), intent(in) ::  dt

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       C_in, & !
       D_sbdiag, & !
       D_spdiag, &
       ML

! !INPUT/OUTPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(:), intent(inout) :: &
       C_low

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

  real(KIND=benthos_r8) :: &
       zspace

  real(KIND=benthos_r8), dimension(nBenthicVertLevels+1) :: &
       M_spdiag, &
       M_sbdiag, &
       F_diag, &
       F_spdiag, &
       F_sbdiag, &
       a_spdiag, &
       a_sbdiag, &
       Pplus, &
       Pminus, &
       Qplus, &
       Qminus, &
       Rplus, &
       Rminus

  integer(KIND=benthos_i4) :: &
       iLevels

   zspace = 1/real(nBenthicVertLevels,KIND=benthos_r8)
   !C_low_new = C_low
   ! compute the mass matrix

   M_spdiag(:) = zspace/6.0_benthos_r8
   M_spdiag(nBenthicVertLevels+1) = c0_benthos
   M_sbdiag(:) = zspace/6.0_benthos_r8
   M_sbdiag(1) = c0_benthos

    ! compute off matrix F

    F_diag(:) = c0_benthos
    F_spdiag(:) =  c0_benthos
    F_sbdiag(:) =  c0_benthos

    do iLevels = 1,nBenthicVertLevels
       F_spdiag(iLevels) = M_spdiag(iLevels)*(C_low(iLevels)-C_in(iLevels) - C_low(iLevels+1)+ &
            C_in(iLevels+1))/dt + D_spdiag(iLevels)*(C_low(iLevels)-C_low(iLevels+1))
       F_sbdiag(iLevels+1) =  M_sbdiag(iLevels+1)*(C_low(iLevels+1)-C_in(iLevels+1) - C_low(iLevels)+ &
            C_in(iLevels))/dt + D_sbdiag(iLevels+1)*(C_low(iLevels+1)-C_low(iLevels))

        if (F_spdiag(iLevels)*(C_low(iLevels) - C_low(iLevels+1)) > c0_benthos) F_spdiag(iLevels) = c0_benthos
        if (F_sbdiag(iLevels+1)*(C_low(iLevels+1) - C_low(iLevels)) > c0_benthos) F_sbdiag(iLevels+1) = c0_benthos
     end do

    if (maxval(abs(F_spdiag)) > c0_benthos) then

     ! compute the weighting factors: a_spdiag, a_sbdiag

        a_spdiag(:) = c0_benthos
        a_sbdiag(:) = c0_benthos

        Pplus(1)  = max(c0_benthos, F_spdiag(1))
        Pminus(1) = min(c0_benthos, F_spdiag(1))
        Pplus(nBenthicVertLevels+1)  = max(c0_benthos, F_sbdiag(nBenthicVertLevels+1))
        Pminus(nBenthicVertLevels+1) = min(c0_benthos, F_sbdiag(nBenthicVertLevels+1))
        Qplus(1) = max(c0_benthos,C_low(2)-C_low(1))
        Qminus(1)= min(c0_benthos,C_low(2)-C_low(1))
        Qplus(nBenthicVertLevels+1) = max(c0_benthos,C_low(nBenthicVertLevels)-C_low(nBenthicVertLevels+1))
        Qminus(nBenthicVertLevels+1)= min(c0_benthos,C_low(nBenthicVertLevels)-C_low(nBenthicVertLevels+1))
        Rplus(1)  = min(c1_benthos, ML(1)*Qplus(1)/dt/(Pplus(1)+puny))
        Rminus(1) = min(c1_benthos, ML(1)*Qminus(1)/dt/(Pminus(1)-puny))
        Rplus(nBenthicVertLevels+1)  = min(c1_benthos, ML(nBenthicVertLevels+1)*Qplus(nBenthicVertLevels+1)&
             /dt/(Pplus(nBenthicVertLevels+1)+puny))
        Rminus(nBenthicVertLevels+1) = min(c1_benthos, ML(nBenthicVertLevels+1)*Qminus(nBenthicVertLevels+1) &
             /dt/(Pminus(nBenthicVertLevels+1)-puny))

        do iLevels = 2,nBenthicVertLevels
           Pplus(iLevels)  = max(c0_benthos,F_spdiag(iLevels)) + max(c0_benthos,F_sbdiag(iLevels))
           Pminus(iLevels) = min(c0_benthos,F_spdiag(iLevels)) + min(c0_benthos,F_sbdiag(iLevels))
           Qplus(iLevels)  = max(c0_benthos, max(C_low(iLevels+1)-C_low(iLevels),C_low(iLevels-1)-C_low(iLevels)))
           Qminus(iLevels) = min(c0_benthos, min(C_low(iLevels+1)-C_low(iLevels),C_low(iLevels-1)-C_low(iLevels)))
           Rplus(iLevels)  = min(c1_benthos, ML(iLevels)*Qplus(iLevels)/dt/(Pplus(iLevels)+puny))
           Rminus(iLevels) = min(c1_benthos, ML(iLevels)*Qminus(iLevels)/dt/(Pminus(iLevels)-puny))
         end do

         do  iLevels = 1,nBenthicVertLevels

            a_spdiag(iLevels) = min(Rminus(iLevels),Rplus(iLevels+1))
            if (F_spdiag(iLevels) > c0_benthos) a_spdiag(iLevels) = min(Rplus(iLevels),Rminus(iLevels+1))
            a_sbdiag(iLevels+1) = min(Rminus(iLevels+1),Rplus(iLevels))
            if (F_sbdiag(iLevels+1) > c0_benthos) a_sbdiag(iLevels+1) = min(Rplus(iLevels+1),Rminus(iLevels))

         end do

         ! compute F_diag

         F_diag(1) = a_spdiag(1)*F_spdiag(1)
         F_diag(nBenthicVertLevels+1) = a_sbdiag(nBenthicVertLevels+1)* F_sbdiag(nBenthicVertLevels+1)
         C_low(1) = C_low(1) + dt*F_diag(1)/ML(1)
         C_low(nBenthicVertLevels+1) = C_low(nBenthicVertLevels+1) + dt*F_diag(nBenthicVertLevels+1)/ &
           ML(nBenthicVertLevels+1)

         do iLevels = 2,nBenthicVertLevels

            F_diag(iLevels) = a_spdiag(iLevels)*F_spdiag(iLevels) + a_sbdiag(iLevels)*F_sbdiag(iLevels)
            C_low(iLevels) = C_low(iLevels) + dt*F_diag(iLevels)/ML(iLevels)

         end do

     end if  !F_spdiag is nonzero

   end subroutine compute_FCT_corr

!*****************************************************************************
!BOP
! !IROUTINE: primaryStoichMatrix
! !INTERFACE:

subroutine primaryStoichMatrix (primarySourceStoich, primarySinkStoich, &
                    benthosTracerBulkLevel)

! !DESCRIPTION:
!
!  computes the stoichiometry matrix for the primary reactions
!
!  primary reactions OM has the form (CH2O)a(NH3)b(H3PO4)c or POC_a PON_b POP_c
!
! P1.  OM/a+ O2 + (b-c)/aH^+ --> CO2 + b/aNH4 + c/aH2PO4 + H2O
! P2.  OM/a + 4/5 NO3 + (4/5 + (b-c)/a)H^+ ---> CO2 + 2/5 N2 + b/aNH4 + c/aH2PO4 + 7/5H2O
! P3.  OM/a + 2MnO2 + 4H + (b-c)/aH --> CO2 + b/aNH4 + c/aH2PO4 + 2Mn + 3H2O
! P4.  OM/a + 4Fe(OH)3 + 4(gamma)[Fe-P] + (12 + (b-c)/a-4gamma)H^+ --> CO2 + b/aNH4 +  (c/a+4gamma) H2PO4 + 4Fe + 13H2O
! P5.  OM/a +  1/2 SO4 + (1+(b-c)/a)H^+ --> CO2 + b/aNH4 + c/aH2PO4 + 1/2H2S + H2O
! P6.  OM/a + (b-c)/aH^+ --> 1/2CO2 + 1/2 CH4 + b/aNH4 + c/aH2PO4 + H2O
!  gamma is the fraction of iron-bound phosphorus
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       benthosTracerBulkLevel

! !INPUT/OUTPUT PARAMETERS:
! !OUTPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(:,:), intent(out) :: &
       primarySourceStoich, &
       primarySinkStoich

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
  !-----------------------------------------------------------------------------
  real(KIND=benthos_r8) :: &
       a_ratio, &
       b_ratio

   a_ratio = CtoP
   b_ratio = NtoP

   primarySourceStoich(:,:) = c0_benthos
   primarySinkStoich(:,:) = c0_benthos

   !define ratios based on concentrations
   !
   ! if (benthosTracerBulkLevel(popa_ind) .gt.  puny)
   !    a_ratio = benthosTracerBulkLevel(poca_ind)/benthosTracerBulkLevel(popa_ind)
   !    b_ratio = benthosTracerBulkLevel(pona_ind)/benthosTracerBulkLevel(popa_ind)
   ! end
   ! Rates are in mmol/m3/s.  Add conversion for solids to mmol/kg/s

   primarySinkStoich(poca_ind,:) = c1_benthos/sediment_density
   primarySinkStoich(pocb_ind,:) = c1_benthos/sediment_density
   primarySinkStoich(pocc_ind,:) = c1_benthos/sediment_density
   primarySinkStoich(pona_ind,:) = b_ratio/a_ratio/sediment_density
   primarySinkStoich(ponb_ind,:) = b_ratio/a_ratio/sediment_density
   primarySinkStoich(ponc_ind,:) = b_ratio/a_ratio/sediment_density
   primarySinkStoich(popa_ind,:) = c1_benthos/a_ratio/sediment_density
   primarySinkStoich(popb_ind,:) = c1_benthos/a_ratio/sediment_density
   primarySinkStoich(popc_ind,:) = c1_benthos/a_ratio/sediment_density
   primarySinkStoich(o2_ind,1) = c1_benthos
   primarySinkStoich(alk_ind,:) = c1_benthos/a_ratio
   primarySinkStoich(alk_ind,4) = (c1_benthos/a_ratio+4.0*ironBoundPFraction)
   primarySinkStoich(no3_ind,2) = 4.0_benthos_r8/5.0_benthos_r8
   primarySinkStoich(mno2a_ind,3) = c2_benthos/sediment_density
   primarySinkStoich(feoh3a_ind,4) = 4.0_benthos_r8/sediment_density
   primarySinkStoich(fepa_ind,4) = 4.0_benthos_r8*ironBoundPFraction/sediment_density
   primarySinkStoich(so4_ind,5) = p5_benthos

   primarySourceStoich(dic_ind,:) = c1_benthos
   primarySourceStoich(dic_ind,6) = p5_benthos
   primarySourceStoich(co2_ind,:) = c1_benthos
   primarySourceStoich(co2_ind,6) = p5_benthos
   primarySourceStoich(alk_ind,1) = b_ratio/a_ratio
   primarySourceStoich(alk_ind,2) = (4.0_benthos_r8/5.0_benthos_r8 + b_ratio/a_ratio)
   primarySourceStoich(alk_ind,3) = (4.0_benthos_r8+b_ratio/a_ratio)
   primarySourceStoich(alk_ind,4) = (12.0_benthos_r8+b_ratio/a_ratio)
   primarySourceStoich(alk_ind,5) = (c1_benthos+b_ratio/a_ratio)
   primarySourceStoich(alk_ind,6) = b_ratio/a_ratio
   primarySourceStoich(nh4_ind,:) = b_ratio/a_ratio
   primarySourceStoich(h2po4_ind,:) = c1_benthos/a_ratio
   primarySourceStoich(h2po4_ind,4) =  (c1_benthos/a_ratio + 4.0_benthos_r8*ironBoundPFraction)
   primarySourceStoich(mn_ind,3) = c2_benthos
   primarySourceStoich(fe_ind,4) = 4.0_benthos_r8
   primarySourceStoich(h2s_ind,5) = p5_benthos
   primarySourceStoich(ch4_ind,6) = p5_benthos

   end subroutine primaryStoichMatrix

!*****************************************************************************
!BOP
! !IROUTINE: primaryRateConstants
! !INTERFACE:

   subroutine primaryRateConstants (primaryRates,benthosTracerBulkLevel,nOMtype)

! !DESCRIPTION:
! compute the rates for the primary redox reactions
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer(KIND=benthos_i4), intent(in) :: &
       nOMtype

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       benthosTracerBulkLevel

! !INPUT/OUTPUT PARAMETERS:
! !OUTPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(nOMtype,nPrimaryReactions), intent(out) :: &
       primaryRates   ! mmol/m3/s

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------
!  Constants from Reed et al 2011
!
  real (KIND=benthos_r8), parameter :: &
       k_alpha = 1.63_benthos_r8/sec_per_year, &   ! per year  Labile rate constant for primary redox of OM^a
       k_beta = 0.0086_benthos_r8/sec_per_year, & ! per year  Semi-labile rate constant for primary redox of OM^b
       k_ref = c0_benthos, &   ! refractory
       k_o2 = 20.0_benthos_r8, & ! mmol/m3  half sat/inhibition for O2
       k_no3 = 4.0_benthos_r8, & ! mmol/m3 half sat/inhibition for NO3
       k_mno2a = 4.0_benthos_r8, & ! mmol/kg half sat/inhibition for MnO2 (umol/g)
       k_feoh3a = 65.0_benthos_r8, & ! mmol/kg half sat/inhibition for fe(OH)3
       k_so4 = 1600.0_benthos_r8, & ! mmol/m3 half sat/inhibition for SO4
       reductionFractionSO4 = 0.075_benthos_r8 ! reduces redox rates for P5-P6

  real(KIND=benthos_r8), dimension(nOMtype) :: &
       k_OM

  real(KIND=benthos_r8), dimension(nOMtype) :: &
       OM

  integer (KIND=benthos_i4) :: &
       iType

   k_OM = (/k_alpha, k_beta, k_ref/)
   OM = (/ benthosTracerBulkLevel(poca_ind), benthosTracerBulkLevel(pocb_ind), &
        benthosTracerBulkLevel(pocc_ind) /)*sediment_density

   primaryRates(:,:) = c0_benthos

   do iType = 1,nOMtype

      primaryRates(iType,1) = &
           k_OM(iType)*OM(iType)*(benthosTracerBulkLevel(o2_ind)/(k_o2 +  benthosTracerBulkLevel(o2_ind)))

      primaryRates(iType,2) = &
           k_OM(iType)*OM(iType)*(benthosTracerBulkLevel(no3_ind)/(k_no3 +  benthosTracerBulkLevel(no3_ind)))* &
           (k_o2/(benthosTracerBulkLevel(o2_ind)+k_o2))

      primaryRates(iType,3) = &
           k_OM(iType)*OM(iType)*(benthosTracerBulkLevel(mno2a_ind)/(k_mno2a + &
           benthosTracerBulkLevel(mno2a_ind)))*(k_o2/(benthosTracerBulkLevel(o2_ind)+k_o2))*&
           (k_no3/(benthosTracerBulkLevel(no3_ind)+k_no3))

      primaryRates(iType,4) = &
           k_OM(iType)*OM(iType)*(benthosTracerBulkLevel(feoh3a_ind)/(k_feoh3a + &
           benthosTracerBulkLevel(feoh3a_ind)))*(k_o2/(benthosTracerBulkLevel(o2_ind)+k_o2))*&
           (k_no3/(benthosTracerBulkLevel(no3_ind)+k_no3))*&
           (k_mno2a/(benthosTracerBulkLevel(mno2a_ind)+k_mno2a))

      primaryRates(iType,5) = &
           reductionFractionSO4*k_OM(iType)*OM(iType)*(benthosTracerBulkLevel(so4_ind)/(k_so4 + &
           benthosTracerBulkLevel(so4_ind)))*(k_o2/(benthosTracerBulkLevel(o2_ind)+k_o2))* &
           (k_no3/(benthosTracerBulkLevel(no3_ind)+k_no3))*&
           (k_mno2a/(benthosTracerBulkLevel(mno2a_ind)+k_mno2a))* &
           (k_feoh3a/(benthosTracerBulkLevel(feoh3a_ind)+k_feoh3a))

      primaryRates(iType,6) = &
           reductionFractionSO4*k_OM(iType)*OM(iType)*(k_so4/(k_so4 + &
           benthosTracerBulkLevel(so4_ind)))*(k_o2/(benthosTracerBulkLevel(o2_ind)+k_o2))*&
           (k_no3/(benthosTracerBulkLevel(no3_ind)+k_no3))*&
           (k_mno2a/(benthosTracerBulkLevel(mno2a_ind)+k_mno2a))* &
           (k_feoh3a/(benthosTracerBulkLevel(feoh3a_ind)+k_feoh3a))

   end do

 end subroutine primaryRateConstants

!*****************************************************************************
!BOP
! !IROUTINE: primaryBenthosReactions
! !INTERFACE:

 subroutine primaryBenthosReactions (primarySourceTendk,primarySinkTendk,&
                    benthosTracerBulkLevels,primarySourceStoich,&
                    primarySinkStoich,primaryRates,dt,nOMtype)

! !DESCRIPTION:
!  compute tendencies for the primary redox reactions
!
!  primary reactions OM has the form (CH2O)a(NH3)b(H3PO4)c or POC_a PON_b POP_c
!
! P1.  OM/a+ O2 + (b-c)/aH^+ --> CO2 + b/aNH4 + c/aH2PO4 + H2O
! P2.  OM/a + 4/5 NO3 + (4/5 + (b-c)/a)H^+ ---> CO2 + 2/5 N2 + b/aNH4 + c/aH2PO4 + 7/5H2O
! P3.  OM/a + 2MnO2 + 4H + (b-c)/aH --> CO2 + b/aNH4 + c/aH2PO4 + 2Mn + 3H2O
! P4.  OM/a + 4Fe(OH)3 + 5(gamma)[Fe-P] + (12 + (b-c)/a-4gamma)H^+ --> CO2 + b/aNH4 +  (c/a+4gamma) H2PO4 + 4Fe + 13H2O
! P5.  OM/a +  1/2 SO4 + (1+(b-c)/a)H^+ --> CO2 + b/aNH4 + c/aH2PO4 + 1/2H2S + H2O
! P6.  OM/a + (b-c)/aH^+ --> 1/2CO2 + 1/2 CH4 + b/aNH4 + c/aH2PO4 + H2O
!  gamma is the fraction of iron-bound phosphorus

!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  integer(KIND=benthos_i4), intent(in) :: &
       nOMtype

  real(KIND=benthos_r8), intent(in) :: &
       dt

  real(KIND=benthos_r8),dimension(:,:), intent(in) :: &
       primarySourceStoich, &
       primarySinkStoich, &
       primaryRates

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       benthosTracerBulklevels

! !INPUT/OUTPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(:),  intent(inout) :: &
       primarySourceTendk, &
       primarySinkTendk

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------
   integer (KIND=benthos_i4) :: &
        iReactions, &
        iType, &
        iTracers

   do iReactions  = 1,nPrimaryReactions

      primarySourceTendk(poca_ind) = &
           primarySourceTendk(poca_ind) + primarySourceStoich(poca_ind,iReactions)*&
           primaryRates(1,iReactions)*dt

      primarySourceTendk(pocb_ind) = &
           primarySourceTendk(pocb_ind) + primarySourceStoich(pocb_ind,iReactions)* &
           primaryRates(2,iReactions)*dt

      primarySourceTendk(pocc_ind) = &
           primarySourceTendk(pocc_ind) + primarySourceStoich(pocc_ind,iReactions)* &
           primaryRates(3,iReactions)*dt

      primarySourceTendk(pona_ind) = &
           primarySourceTendk(pona_ind) + primarySourceStoich(pona_ind,iReactions)* &
           primaryRates(1,iReactions)*dt

      primarySourceTendk(ponb_ind) = &
           primarySourceTendk(ponb_ind) + primarySourceStoich(ponb_ind,iReactions)* &
           primaryRates(2,iReactions)*dt

      primarySourceTendk(ponc_ind) = &
           primarySourceTendk(ponc_ind) + primarySourceStoich(ponc_ind,iReactions)*&
           primaryRates(3,iReactions)*dt

      primarySourceTendk(popa_ind) = &
           primarySourceTendk(popa_ind) + primarySourceStoich(popa_ind,iReactions)* &
           primaryRates(1,iReactions)*dt

      primarySourceTendk(popb_ind) = &
           primarySourceTendk(popb_ind) + primarySourceStoich(popb_ind,iReactions)* &
           primaryRates(2,iReactions)*dt

      primarySourceTendk(popc_ind) = &
           primarySourceTendk(popc_ind) + primarySourceStoich(popc_ind,iReactions)* &
           primaryRates(3,iReactions)*dt

      primarySinkTendk(poca_ind) = &
           primarySinkTendk(poca_ind) + primarySinkStoich(poca_ind,iReactions)* &
           primaryRates(1,iReactions)*dt

      primarySinkTendk(pocb_ind) = &
           primarySinkTendk(pocb_ind) +    primarySinkStoich(pocb_ind,iReactions)* &
           primaryRates(2,iReactions)*dt

      primarySinkTendk(pocc_ind) = &
           primarySinkTendk(pocc_ind) +    primarySinkStoich(pocc_ind,iReactions)* &
           primaryRates(3,iReactions)*dt

      primarySinkTendk(pona_ind) = &
           primarySinkTendk(pona_ind) +    primarySinkStoich(pona_ind,iReactions)* &
           primaryRates(1,iReactions)*dt

      primarySinkTendk(ponb_ind) = &
           primarySinkTendk(ponb_ind) +    primarySinkStoich(ponb_ind,iReactions)* &
           primaryRates(2,iReactions)*dt

      primarySinkTendk(ponc_ind) = &
           primarySinkTendk(ponc_ind) +    primarySinkStoich(ponc_ind,iReactions)* &
           primaryRates(3,iReactions)*dt

      primarySinkTendk(popa_ind) = &
           primarySinkTendk(popa_ind) +    primarySinkStoich(popa_ind,iReactions)* &
           primaryRates(1,iReactions)*dt

      primarySinkTendk(popb_ind) = &
           primarySinkTendk(popb_ind) +    primarySinkStoich(popb_ind,iReactions)* &
           primaryRates(2,iReactions)*dt

      primarySinkTendk(popc_ind) = &
           primarySinkTendk(popc_ind) +    primarySinkStoich(popc_ind,iReactions)* &
           primaryRates(3,iReactions)*dt


      do iTracers = o2_ind,nBenthicTracers
         do iType = 1,3
            primarySourceTendk(iTracers) = primarySourceTendk(iTracers) + &
                 primarySourceStoich(iTracers,iReactions)*primaryRates(iType,iReactions)*dt

            primarySinkTendk(iTracers) = primarySinkTendk(iTracers) + &
                 primarySinkStoich(iTracers,iReactions)*primaryRates(iType,iReactions)*dt
         end do ! iType
      end do ! iTracers
   end do ! iReactions

 end subroutine primaryBenthosReactions

!*****************************************************************************
!BOP
! !IROUTINE: comp_CO3terms
! !INTERFACE:

subroutine comp_CO3terms_benthos (CO3,HCO3,H2CO3,ph,kw,kb, &
                       ks,kf,k1p,k2p,k3p,ksi,bt,st,ft,dic,ta,pt,sit, temp, salt,depth, &
                       lcomp_co3_coeffs, sit_in, phlo, phhi, dic_in,ta_in,pt_in,k, oceanBottomDensity)

! !DESCRIPTION:
!
!---------------------------------------------------------------------------
!   after SUBROUTINE comp_CO3terms
!
!   PURPOSE : Calculate H2CO3, HCO3, CO3 from
!             total alkalinity, total CO2, temp, salinity (s), etc.
!---------------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

  ! !INPUT PARAMETERS:
  integer (KIND=benthos_i4), intent(in) :: k

  real (KIND=benthos_r8), intent(in) :: &
       temp, &
       salt, &
       depth, &
       sit_in, &
       dic_in, &
       ta_in, &
       pt_in, &
       oceanBottomDensity

   logical (KIND=benthos_log), intent(in) :: &
       lcomp_co3_coeffs

! !INPUT/OUTPUT PARAMETERS:
  real(KIND=benthos_r8), intent(inout):: &
       phlo, &
       phhi

! !OUTPUT PARAMETERS:

  real(KIND=benthos_r8), intent(out) :: &
       CO3, &
       HCO3, &
       H2CO3, &
       ph, &
       kw, &
       kb, &
       ks, &
       kf, &
       k1p, &
       k2p, &
       k3p, &
       ksi, &
       bt, &
       st, &
       ft, &
       dic, &
       ta, &
       pt, &
       sit

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

   logical (KIND=benthos_log) :: k1_k2_pH_tot

   integer(KIND=benthos_i4) :: &
        i

   real(KIND=benthos_r8) :: &
        htotal2, denom, &
        massToVol, &
        volToMass

   real(KIND=benthos_r8) :: &
        htotal,       & ! free concentration of H ion
        k0,k1,k2,     & ! equilibrium constants for CO2 species
        ff              ! fugacity of CO2

  !  Initialize

   CO3 = c0_benthos
   HCO3 = c0_benthos
   H2CO3 = c0_benthos
   !phlo_out = phlo
   !phhi_out = phhi
   ph = phlo
   kw = c0_benthos
   kb = c0_benthos
   ks = c0_benthos
   kf = c0_benthos
   k1p = c0_benthos
   k2p = c0_benthos
   k3p = c0_benthos
   ksi = c0_benthos
   bt = c0_benthos
   st = c0_benthos
   ft = c0_benthos
   dic = dic_in
   ta = ta_in
   pt = pt_in
   sit = sit_in

    !------------------------------------------------------------------------
    !   compute thermodynamic CO3 coefficients
    !------------------------------------------------------------------------

    !    IF (lcomp_co3_coeffs) THEN

    k1_k2_pH_tot=.true.

    call comp_co3_coeffs_benthos(k0,k1,k2,ff,kw,kb,ks,kf,k1p,k2p,k3p,ksi,bt,st,ft,k, depth, temp, salt, k1_k2_pH_tot)

    !    END IF

    !------------------------------------------------------------------------
    !   compute htotal
    !------------------------------------------------------------------------

    call  comp_htotal_benthos(htotal,dic,ta,pt,sit, &
         k, temp, dic_in, ta_in, pt_in, sit_in, k1, &
         k2,phlo,phhi,kw,kb,ks,kf,k1p,k2p,k3p,ksi,bt,st,ft, &
         oceanBottomDensity)

    !------------------------------------------------------------------------
    !   Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2,
    !   ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49-51)
    !------------------------------------------------------------------------

    htotal2  = htotal**2
    denom    = 1 / (htotal2 + k1 * htotal + k1 * k2)

    !  this insures H2CO3 + HCO3 + CO3 = 1

    H2CO3 = dic* htotal2 * denom
    HCO3  = dic * k1 * htotal * denom
    CO3   = dic * k1 * k2 * denom

    ph    = -LOG10(htotal)

    !------------------------------------------------------------------
    !   Convert units of output arguments
    !   use (massToVol rather than mass_to_vol)
    !------------------------------------------------------------------
    massToVol = oceanBottomDensity * 1.0e3_benthos_r8
    
    H2CO3 = H2CO3 * massToVol
    HCO3  = HCO3 * massToVol
    CO3   = CO3 * massToVol
    dic = dic * massToVol
    ta = ta * massToVol
    pt = pt * massToVol
    sit = sit * massToVol

   end subroutine comp_CO3terms_benthos

 !*****************************************************************************

   SUBROUTINE comp_co3_coeffs_benthos(sk0,sk1,sk2,sff,kw,kb,ks,kf,k1p,k2p,k3p,ksi,&
        bt,st,ft,k,depth,temp,salt,k1_k2_pH_tot)

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    integer(KIND=benthos_i4), intent(in) :: k
    real(KIND=benthos_r8), intent(in) :: &
         depth,    & ! depth (meters)
         temp,     & ! temperature (degrees C)
         salt        ! salinity (PSU)
    LOGICAL(KIND=benthos_log), intent(in) :: k1_k2_pH_tot

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

!maltrud these are scalar versions--need to copy from array(1) due to shr_vmath
    real(KIND=benthos_r8), intent(out) :: &
         sk0,sk1,sk2,     & ! equilibrium constants for CO2 species
         sff, &              ! fugacity of CO2
         kw, kb, ks, kf, &
         k1p, k2p, k3p, &
         ksi, bt, st, ft
    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    real(KIND=benthos_r8), dimension(1) :: &  ! need to be arrays for shr_vmath
         k0,k1,k2,     & ! equilibrium constants for CO2 species
         ff              ! fugacity of CO2

    integer(KIND=benthos_i4) :: i

    real(KIND=benthos_r8) :: &
         press_bar       ! pressure at level k [bars]

    real(KIND=benthos_r8), dimension(1) :: &  !  need to be arrays to use shr_vmath
         salt_lim,     & ! bounded salt
         tk,           & ! temperature (K)
         is,           & ! ionic strength
         scl,          & ! chlorinity
         tk100, tk1002, invtk, dlogtk, is2, sqrtis, &
         s2, sqrts, s15, invRtk, arg, &
         deltaV,Kappa,lnKfac,Kfac, & ! pressure correction terms
         log_1_m_1p005em3_s, &
         log_1_p_tot_sulfate_div_ks

    !---------------------------------------------------------------------------
    !initialize output

    sk0 = c0_benthos
    sk1 = c0_benthos
    sk2 = c0_benthos
    sff = c0_benthos
    kw = c0_benthos
    kb = c0_benthos
    ks = c0_benthos
    kf = c0_benthos
    k1p = c0_benthos
    k2p = c0_benthos
    k3p = c0_benthos
    ksi = c0_benthos
    bt = c0_benthos
    st = c0_benthos
    ft = c0_benthos

!   press_bar = ref_pressure(k)
!  below is from POP ref_pressure
    press_bar = 0.059808_benthos_r8*(exp(-0.025_benthos_r8*depth) - c1_benthos)     &
            + 0.100766_benthos_r8*depth + 2.28405e-7_benthos_r8*depth**2

    !---------------------------------------------------------------------------
    !   Calculate all constants needed to convert between various
    !   measured carbon species. References for each equation are
    !   noted in the code.  Once calculated, the constants are stored
    !   and passed in the common block "const". The original version
    !   of this code was based on the code by Dickson in Version 2 of
    !   "Handbook of Methods for the Analysis of the Various Parameters
    !   of the Carbon Dioxide System in Seawater", DOE, 1994 (SOP No. 3,
    !   p25-26).
    !   Derive simple terms used more than once
    !---------------------------------------------------------------------------

    salt_lim = max(salt,benthos_salt_min)
    tk       = T0_Kelvin_benthos + temp
    tk100    = tk * 1e-2_benthos_r8
    tk1002   = tk100 * tk100
    invtk    = c1_benthos / tk
#ifdef CCSMCOUPLED
    CALL shr_vmath_log(tk, dlogtk, 1)
#else
    dlogtk   = LOG(tk)
#endif
    invRtk   = (c1_benthos / 83.1451_benthos_r8) * invtk

    is       = 19.924_benthos_r8 * salt_lim / (1000.0_benthos_r8 - 1.005_benthos_r8 * salt_lim)
    is2      = is * is
#ifdef CCSMCOUPLED
    CALL shr_vmath_sqrt(is, sqrtis, 1)
    CALL shr_vmath_sqrt(salt_lim, sqrts, 1)
#else
    sqrtis   = SQRT(is)
    sqrts    = SQRT(salt_lim)
#endif
    s2       = salt_lim * salt_lim
    scl      = salt_lim / 1.80655_benthos_r8

    arg = c1_benthos - 0.001005_benthos_r8 * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_log(arg, log_1_m_1p005em3_s, 1)
#else
    log_1_m_1p005em3_s = LOG(arg)
#endif

    !---------------------------------------------------------------------------
    !   f = k0(1-pH2O)*correction term for non-ideality
    !   Weiss & Price (1980, Mar. Chem., 8, 347-359
    !                 Eq 13 with table 6 values)
    !---------------------------------------------------------------------------

    arg = -162.8301_benthos_r8 + 218.2968_benthos_r8 / tk100 + &
          90.9241_benthos_r8 * (dlogtk + LOG(1e-2_benthos_r8)) - 1.47696_benthos_r8 * tk1002 + &
          salt_lim * (.025695_benthos_r8 - .025225_benthos_r8 * tk100 + 0.0049867_benthos_r8 * tk1002)
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, ff, 1)
#else
    ff = EXP(arg)
#endif
    sff = ff(1)

    !---------------------------------------------------------------------------
    !   K0 from Weiss 1974
    !---------------------------------------------------------------------------

    arg = 93.4517_benthos_r8 / tk100 - 60.2409_benthos_r8 + 23.3585_benthos_r8 * (dlogtk + LOG(1e-2_benthos_r8)) + &
          salt_lim * (.023517_benthos_r8 - 0.023656_benthos_r8 * tk100 + 0.0047036_benthos_r8 * tk1002)
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k0, 1)
#else
    k0 = EXP(arg)
#endif
    sk0 = k0(1)

    !---------------------------------------------------------------------------
    !   k1 = [H][HCO3]/[H2CO3]
    !   k2 = [H][CO3]/[HCO3]
    !   if k1_k2_pH_tot == .true., then use
    !      Lueker, Dickson, Keeling (2000) using Mehrbach et al. data on total scale
    !   otherwise, use
    !      Millero p.664 (1995) using Mehrbach et al. data on seawater scale
    !      this is only present to be consistent w/ OCMIP2 code
    !      it should not be used for new runs
    !      the only reason to use it is to be compatible with prior
    !      long spun up runs that had used it
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    IF (k1_k2_pH_tot) THEN
       ! total pH scale
       arg = 3633.86_benthos_r8 * invtk - 61.2172_benthos_r8 + &
             9.67770_benthos_r8 * dlogtk - 0.011555_benthos_r8 * salt_lim + &
             0.0001152_benthos_r8 * s2
    ELSE
       ! seawater pH scale, see comment above
       arg = 3670.7_benthos_r8 * invtk - 62.008_benthos_r8 + &
             9.7944_benthos_r8 * dlogtk - 0.0118_benthos_r8 * salt_lim + &
             0.000116_benthos_r8 * s2
    END IF
    arg = -LOG(c10_benthos) * arg
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k1, 1)
#else
    k1 = EXP(arg)
#endif
    sk1 = k1(1)

    IF (k > 1) THEN
       deltaV = -25.5_benthos_r8 + 0.1271_benthos_r8 * temp
       Kappa  = (-3.08_benthos_r8 + 0.0877_benthos_r8 * temp) * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       k1 = k1 * Kfac
    END IF

    IF (k1_k2_pH_tot) THEN
       ! total pH scale
       arg = 471.78_benthos_r8 * invtk + 25.9290_benthos_r8 - &
             3.16967_benthos_r8 * dlogtk - 0.01781_benthos_r8 * salt_lim + 0.0001122_benthos_r8 * s2
    ELSE
       ! seawater pH scale, see comment above
       arg = 1394.7_benthos_r8 * invtk + 4.777_benthos_r8 - &
             0.0184_benthos_r8 * salt_lim + 0.000118_benthos_r8 * s2
    END IF
    arg = -LOG(c10_benthos) * arg
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k2, 1)
#else
    k2 = EXP(arg)
#endif
    sk2 = k2(1)

    IF (k > 1) THEN
       deltaV = -15.82_benthos_r8 - 0.0219_benthos_r8 * temp
       Kappa  = (1.13_benthos_r8 - 0.1475_benthos_r8 * temp) * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       k2 = k2 * Kfac
    END IF

    !---------------------------------------------------------------------------
    !   kb = [H][BO2]/[HBO2]
    !   Millero p.669 (1995) using data from Dickson (1990)
    !   CO2SYS states that this in on total pH scale
    !   pressure correction from Millero 1979, p. 1657
    !      omitting salinity contribution
    !---------------------------------------------------------------------------

    arg = (-8966.90_benthos_r8 - 2890.53_benthos_r8 * sqrts - &
           77.942_benthos_r8 * salt_lim + 1.728_benthos_r8 * salt_lim * sqrts - &
           0.0996_benthos_r8 * s2) * invtk + &
          (148.0248_benthos_r8 + 137.1942_benthos_r8 * sqrts + 1.62142_benthos_r8 * salt_lim) + &
          (-24.4344_benthos_r8 - 25.085_benthos_r8 * sqrts - 0.2474_benthos_r8 * salt_lim) * dlogtk + &
          0.053105_benthos_r8 * sqrts * tk
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, kb, 1)
#else
    kb = exp(arg(1))
#endif

    IF (k > 1) THEN
       deltaV = -29.48_benthos_r8 + (0.1622_benthos_r8 - 0.002608_benthos_r8 * temp) * temp
       Kappa  = -2.84_benthos_r8 * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       kb = kb * Kfac(1)
    END IF

    !---------------------------------------------------------------------------
    !   k1p = [H][H2PO4]/[H3PO4]
    !   DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    arg = -4576.752_benthos_r8 * invtk + 115.525_benthos_r8 - &
          18.453_benthos_r8 * dlogtk + &
          (-106.736_benthos_r8 * invtk + 0.69171_benthos_r8) * sqrts + &
          (-0.65643_benthos_r8 * invtk - 0.01844_benthos_r8) * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k1p, 1)
#else
    k1p = EXP(arg(1))
#endif

    IF (k > 1) THEN
       deltaV = -14.51_benthos_r8 + (0.1211_benthos_r8 - 0.000321_benthos_r8 * temp) * temp
       Kappa  = (-2.67_benthos_r8 + 0.0427_benthos_r8 * temp) * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       k1p = k1p * Kfac(1)
    END IF

    !---------------------------------------------------------------------------
    !   k2p = [H][HPO4]/[H2PO4]
    !   DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    arg = -8814.715_benthos_r8 * invtk + 172.0883_benthos_r8 - &
          27.927_benthos_r8 * dlogtk + &
          (-160.340_benthos_r8 * invtk + 1.3566_benthos_r8) * sqrts + &
          (0.37335_benthos_r8 * invtk - 0.05778_benthos_r8) * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k2p, 1)
#else
    k2p = EXP(arg(1))
#endif

    IF (k > 1) THEN
       deltaV = -23.12_benthos_r8 + (0.1758_benthos_r8 - 0.002647_benthos_r8 * temp) * temp
       Kappa  = (-5.15_benthos_r8 + 0.09_benthos_r8 * temp) * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       k2p = k2p * Kfac(1)
    END IF

    !---------------------------------------------------------------------------
    !   k3p = [H][PO4]/[HPO4]
    !   DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    arg = -3070.75_benthos_r8 * invtk - 18.141_benthos_r8 + &
          (17.27039_benthos_r8 * invtk + 2.81197_benthos_r8) * sqrts + &
          (-44.99486_benthos_r8 * invtk - 0.09984_benthos_r8) * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k3p, 1)
#else
    k3p = EXP(arg(1))
#endif

    IF (k > 1) THEN
       deltaV = -26.57_benthos_r8 + (0.202_benthos_r8 - 0.003042_benthos_r8 * temp) * temp
       Kappa  = (-4.08_benthos_r8 + 0.0714_benthos_r8 * temp) * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       k3p = k3p * Kfac(1)
    END IF

    !---------------------------------------------------------------------------
    !   ksi = [H][SiO(OH)3]/[Si(OH)4]
    !   Millero p.671 (1995) using data from Yao and Millero (1995)
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !      apply boric acid values
    !---------------------------------------------------------------------------

    arg = -8904.2_benthos_r8 * invtk + 117.385_benthos_r8 - &
          19.334_benthos_r8 * dlogtk + &
          (-458.79_benthos_r8 * invtk + 3.5913_benthos_r8) * sqrtis + &
          (188.74_benthos_r8 * invtk - 1.5998_benthos_r8) * is + &
          (-12.1652_benthos_r8 * invtk + 0.07871_benthos_r8) * is2 + &
          log_1_m_1p005em3_s
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, ksi, 1)
#else
    ksi = EXP(arg(1))
#endif

    IF (k > 1) THEN
       deltaV = -29.48_benthos_r8 + (0.1622_benthos_r8 - 0.002608_benthos_r8 * temp) * temp
       Kappa  = -2.84_benthos_r8 * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       ksi = ksi * Kfac(1)
    END IF

    !---------------------------------------------------------------------------
    !   kw = [H][OH]
    !   Millero p.670 (1995) using composite data
    !   following DOE Handbook, 0.015 substracted from constant to
    !   approximately convert from SWS pH scale to total pH scale
    !   pressure correction from Millero 1983
    !      note that deltaV coeffs in Millero 1995 are those actually
    !      freshwater deltaV coeffs from Millero 1983
    !---------------------------------------------------------------------------

    arg = -13847.26_benthos_r8 * invtk + 148.9652_benthos_r8 - 23.6521_benthos_r8 * dlogtk + &
          (118.67_benthos_r8 * invtk - 5.977_benthos_r8 + 1.0495_benthos_r8 * dlogtk) * sqrts - &
          0.01615_benthos_r8 * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, kw, 1)
#else
    kw = EXP(arg(1))
#endif

    IF (k > 1) THEN
       deltaV = -20.02_benthos_r8 + (0.1119_benthos_r8 - 0.001409_benthos_r8 * temp) * temp
       Kappa  = (-5.13_benthos_r8 + 0.0794_benthos_r8 * temp) * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       kw = kw * Kfac(1)
    END IF

    !---------------------------------------------------------------------------
    !   ks = [H][SO4]/[HSO4], free pH scale
    !   Dickson (1990, J. chem. Thermodynamics 22, 113)
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    arg = -4276.1_benthos_r8 * invtk + 141.328_benthos_r8 - 23.093_benthos_r8 * dlogtk + &
          (-13856.0_benthos_r8 * invtk + 324.57_benthos_r8 - 47.986_benthos_r8 * dlogtk) * sqrtis + &
          (35474.0_benthos_r8 * invtk - 771.54_benthos_r8 + 114.723_benthos_r8 * dlogtk) * is - &
          2698.0_benthos_r8 * invtk * is * sqrtis + &
          1776.0_benthos_r8 * invtk * is2 + &
          log_1_m_1p005em3_s
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, ks, 1)
#else
    ks = EXP(arg(1))
#endif

    IF (k > 1) THEN
       deltaV = -18.03_benthos_r8 + (0.0466_benthos_r8 + 0.000316_benthos_r8 * temp) * temp
       Kappa  = (-4.53_benthos_r8 + 0.09_benthos_r8 * temp) * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       ks = ks * Kfac(1)
    END IF

    !---------------------------------------------------------------------
    !   kf = [H][F]/[HF]
    !   Dickson and Riley (1979) -- change pH scale to total
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------

    arg = c1_benthos + (0.1400_benthos_r8 / 96.062_benthos_r8) * (scl) / ks
#ifdef CCSMCOUPLED
       CALL shr_vmath_log(arg, log_1_p_tot_sulfate_div_ks, 1)
#else
    log_1_p_tot_sulfate_div_ks = LOG(arg)
#endif
    arg = 1590.2_benthos_r8 * invtk - 12.641_benthos_r8 + 1.525_benthos_r8 * sqrtis + &
          log_1_m_1p005em3_s + log_1_p_tot_sulfate_div_ks
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, kf, 1)
#else
    kf = EXP(arg(1))
#endif

    IF (k > 1) THEN
       deltaV = -9.78_benthos_r8 - (0.009_benthos_r8 + 0.000942_benthos_r8 * temp) * temp
       Kappa  = (-3.91_benthos_r8 + 0.054_benthos_r8 * temp) * 0.001_benthos_r8
       lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       kf = kf * Kfac(1)
    END IF

    !---------------------------------------------------------------------
    !   Calculate concentrations for borate, sulfate, and fluoride
    !   bt : Uppstrom (1974)
    !   st : Morris & Riley (1966)
    !   ft : Riley (1965)
    !---------------------------------------------------------------------

    bt = 0.000232_benthos_r8 / 10.811_benthos_r8 * scl(1)
    st = 0.14_benthos_r8 / 96.062_benthos_r8 * scl(1)
    ft = 0.000067_benthos_r8 / 18.9984_benthos_r8 * scl(1)

  END SUBROUTINE comp_co3_coeffs_benthos

!*****************************************************************************
!BOP
! !IROUTINE: comp_htotal_benthos
! !INTERFACE:

  subroutine comp_htotal_benthos(htotal,dic,ta,pt,sit, &
       k, temp, dic_in, ta_in, pt_in, sit_in, &
       k1, k2,phlo,phhi,kw,kb,ks,kf,k1p,k2p,k3p,ksi,bt,st,ft, oceanBottomDensity)

! !DESCRIPTION:
!
!---------------------------------------------------------------------------
!   SUBROUTINE comp_htotal
!
!   PURPOSE : Calculate htotal from total alkalinity, total CO2,
!             temp, salinity (s), etc.
!---------------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

    integer(KIND=benthos_i4), intent(in) :: k
    real(KIND=benthos_r8), intent(in) :: &
         temp,     & ! temperature (degrees C)
         dic_in,   & ! total inorganic carbon (nmol/cm^3)
         ta_in,    & ! total alkalinity (neq/cm^3)
         pt_in,    & ! inorganic phosphate (nmol/cm^3)
         sit_in,   & ! inorganic silicate (nmol/cm^3)
         k1,k2       ! equilibrium constants for CO2 species

    real(KIND=benthos_r8), intent(in) :: &
         kw, kb, ks, &
         kf, k1p, k2p, &
         k3p, ksi, bt, st, ft, &
         oceanBottomDensity

! !INPUT/OUTPUT PARAMETERS:

    real(KIND=benthos_r8), intent(inout) :: &
         phlo,     & ! lower limit of pH range
         phhi        ! upper limit of pH range

! !OUTPUT PARAMETERS:

    real(KIND=benthos_r8), intent(out) :: &
         htotal, &      ! free concentration of H ion
         dic, &
         ta, &
         pt, &
         sit

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

    integer(KIND=benthos_i4) :: i

    real(KIND=benthos_r8) :: &
         x1, x2, &         ! bounds on htotal for solver
         xacc, &
         volToMass

    !---------------------------------------------------------------------------
    !   convert tracer units to per mass
    !  (switch vol_to_mass to volToMass using prognostic density)
    !---------------------------------------------------------------------------

    volToMass = c1_benthos/1.0e3_benthos_r8/oceanBottomDensity
    
    dic  = max(dic_in,benthos_dic_min) * volToMass
    ta   = max(ta_in,benthos_alk_min)  * volToMass
    pt   = max(pt_in,c0_benthos)       * volToMass
    sit  = max(sit_in,c0_benthos)      * volToMass

    x1 = c10_benthos**(-phhi)
    x2 = c10_benthos**(-phlo)

    !---------------------------------------------------------------------------
    !   If DIC and TA are known then either a root finding or iterative
    !   method must be used to calculate htotal. In this case we use
    !   the Newton-Raphson "safe" method taken from "Numerical Recipes"
    !   (function "rtsafe.f" with error trapping removed).
    !
    !   As currently set, this procedure iterates about 12 times. The
    !   x1 and x2 values set below will accomodate ANY oceanographic
    !   values. If an initial guess of the pH is known, then the
    !   number of iterations can be reduced to about 5 by narrowing
    !   the gap between x1 and x2. It is recommended that the first
    !   few time steps be run with x1 and x2 set as below. After that,
    !   set x1 and x2 to the previous value of the pH +/- ~0.5.
    !---------------------------------------------------------------------------
!      x1_in = x1
!      x2_in  = x2
      xacc = c10_benthos **(-10)

      call  drtsafe_row_benthos(htotal, k, k1, k2, x1, x2, xacc, dic,kw,kb,ks,kf,k1p,k2p,k3p,ksi,bt,st,ft,ta,pt,sit)

    end subroutine comp_htotal_benthos

!*****************************************************************************
!BOP
! !IROUTINE: drtsafe_row_benthos
! !INTERFACE:

   subroutine drtsafe_row_benthos(soln, k, k1, k2, x1, x2, xacc, dic,kw,kb,ks,kf,k1p,&
        k2p,k3p,ksi,bt,st,ft,ta,pt,sit)

! !DESCRIPTION:
!
!---------------------------------------------------------------------------
! from drtsafe_row:
!   Vectorized version of drtsafe, which was a modified version of
!      Numerical Recipes algorithm.
!   Keith Lindsay, Oct 1999
!
!   Algorithm comment :
!      Iteration from Newtons method is used unless it leaves
!      bracketing interval or the dx is > 0.5 the previous dx.
!      In that case, bisection method is used.
!---------------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

    integer(KIND=benthos_i4), intent(in) :: k
    real(KIND=benthos_r8), intent(in) :: &
         k1, k2,&
         xacc, &
         dic, &
         kw, kb, ks, kf, &
         k1p, k2p, k3p, &
         ksi, bt, st, ft, &
         ta, pt, sit

!  !INPUT/OUTPUT PARAMETERS:

    real(KIND=benthos_r8), intent(INOUT) :: x1, x2

!  !OUTPUT PARAMETERS:

    real(KIND=benthos_r8), intent(OUT) :: soln

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

    logical(KIND=benthos_log) :: leave_bracket, dx_decrease, mask
    integer(KIND=benthos_i4) ::  i, it, maxit, max_bracket_grow_it
    real(KIND=benthos_r8) :: temp
    real(KIND=benthos_r8) :: xlo, xhi, flo, fhi, f, df, dxold, dx

!---------------------------------------------------------------------------
!   bracket root at each location and set up first iteration
!---------------------------------------------------------------------------
    max_bracket_grow_it = 3
    maxit = 100
    it = 0

    DO
!    mask =. true.
!    while (mask)


       call talk_row_benthos(flo,df,k1, k2, x1, dic,kw,kb,ks,kf,k1p,k2p,k3p,ksi,bt,st,ft,ta,pt,sit)
       call talk_row_benthos(fhi,df,k1, k2, x2, dic,kw,kb,ks,kf,k1p,k2p,k3p,ksi,bt,st,ft,ta,pt,sit)
       mask = (flo .gt.  c0_benthos .and.  fhi .gt. c0_benthos) .or. &
            (flo .lt. c0_benthos .and. fhi .lt. c0_benthos)

       if (.not. mask) EXIT

       it = it + 1

    !    if (it .gt. max_bracket_grow_it) then
    !         'bounding bracket for pH solution not found'
    !        'drtsafe_row'
    !         return
    !         CALL shr_sys_abort('bounding bracket for pH solution not found')
    !       end if


       dx = sqrt(x2 / x1)
       x2 = x2 * dx
       x1 = x1 / dx
    end do

    if (flo .lt. c0_benthos) then
       xlo = x1
       xhi = x2
    else
       xlo = x2
       xhi = x1
       temp = flo
       flo = fhi
       fhi = temp
    end if
    soln = p5_benthos* (xlo+xhi)
    dxold = abs(xlo - xhi)
    dx = dxold

    call talk_row_benthos(f,df,k1, k2, soln, dic,kw,kb,ks,kf,k1p,k2p,k3p,ksi,bt,st,ft,ta,pt,sit)

    !---------------------------------------------------------------------------
    !   perform iterations, zeroing mask when a location has converged
    !---------------------------------------------------------------------------

    mask = .true.
    do it = 1,maxit
       leave_bracket = ((soln - xhi) * df - f) * &
       ((soln - xlo) * df - f) .ge. c0_benthos
       dx_decrease = abs(c2_benthos * f) .le. abs(dxold * df)
       if (leave_bracket .or. .not. dx_decrease) then
          dxold = dx
          dx = p5_benthos * (xhi - xlo)
          soln = xlo + dx
          if (xlo .eq. soln) mask = .false.
       else
          dxold = dx
          dx = -f / df
          temp = soln
          soln = soln + dx
          if (temp .eq. soln) mask = .false.
       end if
       if (abs(dx) .lt. xacc) mask = .false.

       if (.not. mask) return

       call talk_row_benthos(f,df,k1, k2, soln, dic,kw,kb,ks,kf,&
            k1p,k2p,k3p,ksi,bt,st,ft,ta,pt,sit)

       if (f .lt. c0_benthos) then
          xlo = soln
          flo = f
       else
          xhi = soln
          fhi = f
       end  if
    end do

#ifdef CCSMCOUPLED
!   CALL shr_sys_abort('lack of convergence in drtsafe_row')
#endif

  end subroutine drtsafe_row_benthos

!*****************************************************************************
!BOP
! !IROUTINE: talk_row_benthos
! !INTERFACE:

   subroutine talk_row_benthos(fn,df, &
        k1, k2, x,dic,kw,kb,ks,kf,k1p,k2p,k3p,ksi,bt,st,ft,ta,pt,sit)

! !DESCRIPTION:
! from talk_row:
!---------------------------------------------------------------------------
!   This routine computes TA as a function of DIC, htotal and constants.
!   It also calculates the derivative of this function with respect to
!   htotal. It is used in the iterative solution for htotal. In the call
!   "x" is the input value for htotal, "fn" is the calculated value for
!   TA and "df" is the value for dTA/dhtotal.
!---------------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
     real(KIND=benthos_r8), intent(in) :: &
          k1, k2, x, &
          dic, kw, kb, ks, kf, &
          k1p, k2p, k3p, ksi, &
          bt, st, ft, ta, pt, sit

! !INPUT/OUTPUT PARAMETERS:
! !OUTPUT PARAMETERS:

     real(KIND=benthos_r8), intent(out) :: &
          fn, &
          df

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

     integer(KIND=benthos_i4) :: i
     real(KIND=benthos_r8) :: &
         x1, x1_r, x2, x2_r, x3, k12, k12p, k123p, &
         a, a_r, a2_r, da, b, b_r, b2_r, db, c, c_r, &
         kb_p_x1_r, ksi_p_x1_r, c1_p_c_ks_x1_r_r, c1_p_kf_x1_r_r

     x1 = x
     x1_r = c1_benthos / x1
     x2 = x1 * x1
     x2_r = x1_r * x1_r
     x3 = x2 * x1
     k12 = k1 * k2
     k12p = k1p * k2p
     k123p = k12p * k3p
     a = x3 + k1p * x2 + k12p * x1 + k123p
     a_r = c1_benthos/ a
     a2_r = a_r * a_r
     da = 3.0_benthos_r8 * x2 + c2_benthos * k1p * x1 + k12p
     b = x2 + k1 * x1 + k12
     b_r = c1_benthos / b
     b2_r = b_r * b_r
     db = c2_benthos * x1 + k1
     c = c1_benthos + st/ks
     c_r = c1_benthos/ c
     kb_p_x1_r = c1_benthos / (kb + x1)
     ksi_p_x1_r = c1_benthos/ (ksi + x1)
     c1_p_c_ks_x1_r_r = c1_benthos / (c1_benthos + c * ks * x1_r)
     c1_p_kf_x1_r_r = c1_benthos/ (c1_benthos+ kf * x1_r)

    !---------------------------------------------------------------------
    !   fn = hco3+co3+borate+oh+hpo4+2*po4+silicate-hfree-hso4-hf-h3po4-ta
    !---------------------------------------------------------------------

     fn = k1 * dic * x1 * b_r &
          + c2_benthos * dic * k12 * b_r &
          + bt * kb * kb_p_x1_r &
          + kw * x1_r &
          + pt * k12p * x1 * a_r &
          + c2_benthos * pt * k123p * a_r &
          + sit * ksi * ksi_p_x1_r &
          - x1 * c_r &
          - st * c1_p_c_ks_x1_r_r &
          - ft * c1_p_kf_x1_r_r &
          - pt * x3 * a_r &
          - ta

     !---------------------------------------------------------------------
     !   df = d(fn)/dx
     !---------------------------------------------------------------------

     df = k1 * dic * (b - x1 * db) * b2_r &
          - c2_benthos * dic * k12 * db * b2_r &
          - bt * kb * kb_p_x1_r * kb_p_x1_r &
          - kw * x2_r &
          + (pt * k12p * (a - x1 * da)) * a2_r &
          - c2_benthos * pt * k123p * da * a2_r &
          - sit * ksi * ksi_p_x1_r * ksi_p_x1_r &
          - c1_benthos * c_r &
          - st * c1_p_c_ks_x1_r_r * c1_p_c_ks_x1_r_r * (c * ks * x2_r) &
          - ft * c1_p_kf_x1_r_r * c1_p_kf_x1_r_r * kf * x2_r &
          - pt * x2 * (3.0_benthos_r8 * a - x1 * da) * a2_r

   end subroutine talk_row_benthos

!*****************************************************************************
!BOP
! !IROUTINE: comp_co3_sat_benthos
! !INTERFACE:

subroutine comp_co3_sat_benthos(sat_calc, sat_arag, sat_mgcalc,&
                       K_arag,K_calc,K_mgcalc,temp,salt, &
                       depth,mg_ions,co3_mol,oceanBottomDensity)

! !DESCRIPTION:
!
!---------------------------------------------------------------------------
!   Based on SUBROUTINE after comp_co3_sat_vals
!
!   PURPOSE : Calculate co3 concentration at calcite,  aragonite and 15! Mg-calcite saturation 
!             from temp, salinity (s), press
!---------------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

    real(KIND=benthos_r8), intent(in) :: &
         depth,      & ! depth (meters)
         temp,       & ! temperature (degrees C)
         salt,       & ! salinity (PSU)
         mg_ions,    & ! magnesium ion concentration (mmol/m3)
         co3_mol,    & ! bottom ocean density
         oceanBottomDensity

! !INPUT/OUTPUT PARAMETERS:
! !OUTPUT PARAMETERS:
    real(KIND=benthos_r8), intent(out) :: &
         sat_calc,  & ! co3 concentration at calcite saturation : change to sat_calc  for saturation state
         sat_arag,  & ! co3 concentration at aragonite saturation
         sat_mgcalc ! co3 concentration at 15% magnesium calcite

    real(KIND=benthos_r8), intent(out):: &
         K_calc,       & ! thermodynamic constant for calcite
         K_arag,       & ! thermodynamic constant for aragonite
         K_mgcalc       ! thermodynamic constant for mg calcite

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

    integer(KIND=benthos_i4) :: i, j

    real(KIND=benthos_r8) :: &
         press_bar, &          ! pressure at level k [bars]
         volToMass

    real(KIND=benthos_r8), dimension(1) :: &  !  need to be arrays to use shr_vmath
         salt_lim,     & ! bounded salt
         tk,             & ! temperature (K)
         log10tk,invtk,sqrts,s15,invRtk,arg,&
         deltaV,Kappa, & ! pressure correction terms
         lnKfac,Kfac,    & ! pressure correction terms
         inv_K_calc,     & ! inverse K_calc
         inv_K_arag,     & ! inverse K_arag
         inv_K_mgcalc,   & ! inverse K_mgcalc
         inv_Ca, &         ! inverse of Calcium concentration
         ca_ions, &
         ca85, &
         mg15

    !---------------------------------------------------------------------------
    !   set unit conversion factors
    !  mmol/kg to mmol/m^3
    !---------------------------------------------------------------------------

    volToMass = c1_benthos/1.0e3_benthos_r8/oceanBottomDensity

    press_bar = 0.059808 * (exp(-0.025_benthos_r8 * depth) - c1_benthos) + &
         0.100766_benthos_r8 * depth + 2.28405_benthos_r8*c10_benthos**(-7) * depth**2


    salt_lim = max(salt,benthos_salt_min)
    tk       = T0_Kelvin_benthos + temp

#ifdef CCSMCOUPLED
       CALL shr_vmath_log(tk, log10tk, 1)
#else
       log10tk  = log(tk)
#endif

       !maltrud dividing by non-shared log(10) here.  i guess since it is a scalar?

       log10tk  = log10tk/log(c10_benthos)
       invtk    = c1_benthos/tk
       invRtk   = (c1_benthos/83.1451_benthos_r8)*invtk

#ifdef CCSMCOUPLED
       CALL shr_vmath_sqrt(salt_lim, sqrts, 1)
#else
       sqrts    =sqrt(salt_lim)
#endif
       s15      = sqrts * salt_lim

       !------------------------------------------------------------------------
       !   Compute K_calc, K_arag, apply pressure factor
       !   Mucci, Amer. J. of Science 283:781-799, 1983 & Millero 1979
       !------------------------------------------------------------------------

       arg = -171.9065_benthos_r8 - 0.077993_benthos_r8 * tk + 2839.319_benthos_r8 * invtk + &
            71.595_benthos_r8 * log10tk + &
             (-0.77712_benthos_r8 + 0.0028426_benthos_r8 * tk + 178.34_benthos_r8 * invtk) * sqrts - &
             0.07711_benthos_r8 * salt_lim + 0.0041249_benthos_r8 * s15

       arg = log(c10_benthos) * arg
#ifdef CCSMCOUPLED
!       CALL shr_vmath_exp(arg, K_calc, 1)
#else
       K_calc = exp(arg(1))
#endif

!       if (k .ge. 1) then
          deltaV = -48.76_benthos_r8 + 0.5304_benthos_r8 * temp
          Kappa  = (-11.76_benthos_r8 + 0.3692_benthos_r8 * temp) * 0.001_benthos_r8
          lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk

#ifdef CCSMCOUPLED
          CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
          Kfac = exp(lnKfac)
#endif

          K_calc = K_calc * Kfac(1)    ! mol/kg
!       end if

       arg = -171.945_benthos_r8 - 0.077993_benthos_r8 * tk + 2903.293_benthos_r8 * invtk + &
            71.595_benthos_r8 * log10tk +  &
           (-0.068393_benthos_r8 + 0.0017276_benthos_r8 * tk + 88.135_benthos_r8 * invtk) * sqrts -  &
            0.10018_benthos_r8 * salt_lim + 0.0059415_benthos_r8 * s15

       arg = log(c10_benthos) * arg
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(arg, K_arag, 1)
#else
      K_arag = exp(arg(1))
#endif

!       if  (k .ge. 1) then
          deltaV = deltaV + 2.8_benthos_r8
          lnKfac = (-deltaV + p5_benthos * Kappa * press_bar) * press_bar * invRtk

#ifdef CCSMCOUPLED
          CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else

          Kfac = exp(lnKfac) !EXP(lnKfac)
#endif
          K_arag = K_arag * Kfac(1)    !mol/kg
!       end if

       !
       ! Mg -Calcite
       !  The solubility factor is assumed to be 21! greater than that of  aragonite
       !  Morse et al (2006).

       K_mgcalc = K_arag * 1.21_benthos_r8    !mol/kg

    !------------------------------------------------------------------
    !   Compute CO3 concentration at calcite & aragonite saturation
    !  but, we are now keeping track of calcium ions...
    !------------------------------------------------------------------

    !   inv_Ca(1) = (35.0_benthos_r8 / 0.01028_benthos_r8) / salt_lim(1)
   !   co3_sat_calc = K_calc(1) * inv_Ca(1)   ! only need to be arrays if using shared math
   !   co3_sat_arag = K_arag(1) * inv_Ca(1)
   !   co3_sat_mgcalc = K_mgcalc(1) * inv_Ca(1)

     inv_K_calc = c1_benthos/K_calc  !(kg/mol)^2
     ca_ions = 0.01028_benthos_r8/35.0_benthos_r8*salt_lim  !(mol/kg)
     sat_calc = co3_mol * volToMass * ca_ions(1) * inv_K_calc(1)  ! (co3_mol from mmol/m3 to mol/kg)  now unitless

     inv_K_arag = c1_benthos/K_arag
     sat_arag = co3_mol* volToMass * ca_ions(1) * inv_K_arag(1)

     inv_K_mgcalc = c1_benthos/K_mgcalc
     ca85 = (ca_ions)**(0.85_benthos_r8)
     mg15 = (mg_ions)**(0.15_benthos_r8)
     sat_mgcalc = ca85(1) * mg15(1) * inv_K_mgcalc(1)

   end subroutine comp_co3_sat_benthos

!*****************************************************************************
!BOP
! !IROUTINE: carbonateRateConstants
! !INTERFACE:

subroutine carbonateRateConstants(carbonateRates,benthosTracerBulkLevels,&
     dt,sat_calc, sat_arag,sat_mgcalc, oceanBottomDensity)

! !DESCRIPTION:
!  compute the rate constants for calcium carbonate dissociation
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  real(KIND=benthos_r8), intent(in) :: &
       dt, &
       sat_calc, &
       sat_arag, &
       sat_mgcalc, &
       oceanBottomDensity

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       benthosTracerBulkLevels

! !INPUT/OUTPUT PARAMETERS:
! !OUTPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(nCarbonateReactions), intent(out) :: &
       carbonateRates

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
  !-----------------------------------------------------------------------------
  real (KIND=benthos_r8) :: &
       k_p1, k_p2, k_p3, & ! mmol/kg/hr (Walter and Morse 1985);  (per mM per y) Krumins et al.
       n1, n2, n3, &   ! power law
       inhibition, &   ! inhibition factor of 0.1 to reduce rates relative to experimental result (Krumins et al)
       sec_per_hour    ! seconds per hour

   !  Constants from Krumins et al.  Something odd about the units

   ! k_p1 = 40.0_benthos_r8/sec_per_year/mM_umolperL  ! /(mmol/m3)/s   from  40.0 per mM per yr
   ! k_p2 = 110.0_benthos_r8/sec_per_year/mM_umolperL ! /(mmol/m3)/s from 110.0 per mM per yr
   ! k_p3 = 58.0_benthos_r8/sec_per_year/mM_umolperL ! /(mmol/m3)/s from 58 per mM per yr

   ! Constants from Walter and Morse (1985) reduced by inhibition factor of 0.1 as in Krumins et al
   inhibition = 0.1_benthos_r8
   sec_per_hour = 3600.0_benthos_r8

   k_p1 = 457.1_benthos_r8 * inhibition * sec_per_hour * sediment_density  ! barnacle Balanus
   k_p2 = 1258.9_benthos_r8 * inhibition * sec_per_hour * sediment_density ! green algae Halimeda
   k_p3 = 660.69_benthos_r8 * inhibition * sec_per_hour * sediment_density ! foraminifier Peneroplis
   n1 = 2.74_benthos_r8  ! power law for p1
   n2 = 2.43_benthos_r8
   n3 = 3.61_benthos_r8

   carbonateRates(:) = c0_benthos

   carbonateRates(1) = k_p1 * (c1_benthos - min(c1_benthos,sat_calc))**n1   ! mmol/m3/s
   carbonateRates(2) = k_p2 * (c1_benthos - min(c1_benthos,sat_arag))**n2   ! mmol/m3/s
   carbonateRates(3) = k_p3 * (c1_benthos - min(c1_benthos,sat_mgcalc))**n3 ! mmol/m3/s

   if (benthosTracerBulkLevels(caco3a_ind)*0.9_benthos_r8 < carbonateRates(1) * dt/sediment_density) then
      carbonateRates(1) = benthosTracerBulkLevels(caco3a_ind)*0.9_benthos_r8/dt*sediment_density
   end if

   if (benthosTracerBulkLevels(caco3b_ind)*0.9_benthos_r8 < carbonateRates(2) * dt/sediment_density) then
      carbonateRates(2) = benthosTracerBulkLevels(caco3b_ind)*0.9_benthos_r8/dt*sediment_density
   end if

   if (benthosTracerBulkLevels(camgco3_ind)*0.9_benthos_r8 < carbonateRates(3) * dt/sediment_density) then
      carbonateRates(3) = benthosTracerBulkLevels(camgco3_ind)*0.9_benthos_r8/dt*sediment_density
   end if

 end subroutine carbonateRateConstants

!*****************************************************************************
!BOP
! !IROUTINE: carbonateBenthosReactions
! !INTERFACE:

subroutine carbonateBenthosReactions (carbonateSourceTendk,carbonateSinkTendk,&
                    carbonateSourceStoich,carbonateSinkStoich,carbonateRates,dt)

! !DESCRIPTION:
!  Initialize tracegas tracer module.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  real(KIND=benthos_r8), intent(in) :: &
       dt

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       carbonateRates  ! 1/s

  real(KIND=benthos_r8), dimension(:,:), intent(in) :: &
       carbonateSourceStoich, &
       carbonateSinkStoich

! !INPUT/OUTPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(nBenthicTracers), intent(inout) :: &
       carbonateSourceTendk, &
       carbonateSinkTendk

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------
  integer (KIND=benthos_i4) :: &
       iTracers, &
       iReactions

   carbonateSourceTendk(:) = c0_benthos
   carbonateSinkTendk(:) = c0_benthos

   do iTracers = 1,nBenthicTracers
      do iReactions = 1,nCarbonateReactions

         if (carbonateSourceStoich(iTracers,iReactions) .gt. c0_benthos .and. &
              carbonateRates(iReactions) .gt. puny) then

            carbonateSourceTendk(iTracers) = carbonateSourceTendk(iTracers) +&
                 carbonateSourceStoich(iTracers,iReactions)*carbonateRates(iReactions)*dt

         end if

         if (carbonateSinkStoich(iTracers,iReactions) .gt. c0_benthos .and. &
              carbonateRates(iReactions) .gt. puny) then
            carbonateSinkTendk(iTracers) = carbonateSinkTendk(iTracers) +&
                 carbonateSinkStoich(iTracers,iReactions)*carbonateRates(iReactions)*dt

         end if
      end do ! iReactions
   end do    ! iTracers

 end subroutine carbonateBenthosReactions

!*****************************************************************************
!BOP
! !IROUTINE: secondaryRateConstants
! !INTERFACE:

 subroutine secondaryRateConstants(secondaryRates, benthosTracerBulkLevels)

! !DESCRIPTION:
!  Initialize tracegas tracer module.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real(KIND=benthos_r8), dimension(:), intent(in) ::&
        benthosTracerBulkLevels

! !INPUT/OUTPUT PARAMETERS:
! !OUTPUT PARAMETERS:

   real(KIND=benthos_r8), dimension(nSecondaryReactions), intent(out) :: &
        secondaryRates
!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

   !  Constants from Reed et al 2011 in benthos_parms.F90

   secondaryRates(:) = c0_benthos

   ! all rates changed to mmol/m3/s

   secondaryRates(1) = max(c0_benthos,k_s1*benthosTracerBulkLevels(nh4_ind)  * benthosTracerBulkLevels(o2_ind))  !mmol/m3/s

   secondaryRates(2) = max(c0_benthos,k_s2 * benthosTracerBulkLevels(mn_ind)* benthosTracerBulkLevels(o2_ind))  ! mmol/m3/s

   secondaryRates(3) =  max(c0_benthos,k_s3*benthosTracerBulkLevels(fe_ind) * benthosTracerBulkLevels(o2_ind) )!mmol/m3/s

   secondaryRates(4) = max(c0_benthos,k_s4 * benthosTracerBulkLevels(fes_ind)* benthosTracerBulkLevels(o2_ind)*sediment_density)  ! mmol/m3/s

   secondaryRates(5) = max(c0_benthos,k_s5 * benthosTracerBulkLevels(fes2_ind)* benthosTracerBulkLevels(o2_ind)*sediment_density)  !  umol/g/s

   secondaryRates(6) = max(c0_benthos,k_s6*benthosTracerBulkLevels(h2s_ind)  * benthosTracerBulkLevels(o2_ind) )   !mmol/m3/s

   secondaryRates(7) =max(c0_benthos,k_s7*benthosTracerBulkLevels(ch4_ind)  * benthosTracerBulkLevels(o2_ind) )   !mmol/m3/s

   secondaryRates(8) = max(c0_benthos,k_s8 * benthosTracerBulkLevels(mno2a_ind)* benthosTracerBulkLevels(fe_ind)*sediment_density )  ! umol/g/s

   secondaryRates(9) =max(c0_benthos,k_s8 * benthosTracerBulkLevels(mno2b_ind)* benthosTracerBulkLevels(fe_ind)*sediment_density)  !  umol/g/s

   secondaryRates(10) = max(c0_benthos,k_s9 * benthosTracerBulkLevels(mno2a_ind)* benthosTracerBulkLevels(h2s_ind)*sediment_density)  ! umol/g/s

   secondaryRates(11) = max(c0_benthos,k_s9 * benthosTracerBulkLevels(mno2b_ind)*benthosTracerBulkLevels(h2s_ind)*sediment_density)  !  umol/g/s

   secondaryRates(12) =max(c0_benthos,k_s10 * benthosTracerBulkLevels(feoh3a_ind)*benthosTracerBulkLevels(h2s_ind)*sediment_density) ! umol/g/s

   secondaryRates(13) = max(c0_benthos,k_s10 * benthosTracerBulkLevels(feoh3b_ind)* benthosTracerBulkLevels(h2s_ind)*sediment_density)   ! umol/g/s

   secondaryRates(14) = max(c0_benthos,k_s11*benthosTracerBulkLevels(h2s_ind)* benthosTracerBulkLevels(fe_ind))    !mmol/m3/s

   secondaryRates(15) = max(c0_benthos,k_s12*benthosTracerBulkLevels(ch4_ind)* benthosTracerBulkLevels(so4_ind)) ! mmol/m3/s

   secondaryRates(16) = max(c0_benthos,k_s13 * benthosTracerBulkLevels(s_ind)*sediment_density)   ! mmol/m3/s

   secondaryRates(17) = max(c0_benthos,k_s14*benthosTracerBulkLevels(s_ind) * benthosTracerBulkLevels(fes_ind)*sediment_density) !  umol/g/s

   secondaryRates(18) = max(c0_benthos,k_s15 * benthosTracerBulkLevels(feoh3a_ind)) !umol/g/s

   secondaryRates(19) = max(c0_benthos,k_s16 * benthosTracerBulkLevels(mno2a_ind)) !umol/g/s

 end subroutine secondaryRateConstants

!*****************************************************************************
!BOP
! !IROUTINE: secondaryBenthosReactions
! !INTERFACE:

subroutine secondaryBenthosReactions(secondarySourceTendk,secondarySinkTendk, &
     secondarySourceStoich,secondarySinkStoich,secondaryRates,dt)

! !DESCRIPTION:
!
!  compute tendencies for secondary reactions
!
! s1. O2 + NH4/2 + HCO3 --> NO3/2 + CO2 + 3/2H2O
! ** s2. O2 + 2Mn + 4HCO3 --> 2MnO2^a + 4CO2 + 2H2O
! ** s3. O2 + 4Fe + 8HCO3 + 2 H20 + 4 gamma H2PO4 + 24 gamma H^+ --> 4Fe(OH)3^a + 4 gamma [Fe-P]^a + 8CO2 + 16 gamma H2O
! s4. O2 + FeS/2 --> SO4/2 + Fe/2
! s5. O2+ 2/7FeS2 +  2/7H2O --> 4/7SO4 +  2/7Fe+ 4/7H
! s6. O2 + H2S/2 + HCO3 -->  SO4/2 + CO2 + H2O
! s7. O2 + CH4 --> CO2 + 4H^+
! s8. MnO2a + 2Fe + 2 gamma H2PO4 + 2H2O + 2HCO3 + 12gamma H --> 2Fe(OH)3^a +  2 gamma [Fe- P]^a + Mn + 2CO2 + 8gamma H2O
! (s9)s8b. MnO2b + 2Fe + 2 gamma H2PO4 + 2H2O + 2HCO3 + 12gammaH --> 2Fe(OH)3^a +  2 gamma [Fe- P]^a + Mn + 2CO2 + 8gamma H2O
! (s10)s9. Mn02a + H2S +  2CO2 -->  Mn +  S +  2HCO3
! (s11)s9b. Mn02b + H2S +  2CO2 -->  Mn +  S +  2HCO3
! (s12)s10. Fe(OH)3a+   gamma[Fe-P]a + H2S/2+ 2CO2 + 4 gamma H2O--> Fe +   gamma H2PO4 + S/2 + 2HCO3 + (2+6gamma)H
! (s13)s10b. Fe(OH)3b+   gamma[Fe-P]b + H2S/2+ 2CO2 + 4 gamma H2O--> Fe +  gamma H2PO4+ S/2 + 2HCO3 + (2+6 gamma)H
! (s14)s11. Fe + H2S -->  FeS + 2H
! (s15)s12. SO4 +  CH4 +  CO2 -->  2HCO3 + H2S
! (s16)s13. S +  H2O -->  3/4H2S +  SO4/4 +  H/2
! (s17)s14. FeS +  S --> FeS2
! (s18)s15. Fe(OH)3a + gamma [Fe-P]a --> Fe(OH)3b + gamma [Fe-P]b
! (s19)s16 MnO2a --> MnO2b
!
! 19 total secondary equations with 3 repetitions
! gamma is the fraction of  co-precipitated P with iron

!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  real(KIND=benthos_r8), dimension(:), intent(in) :: &
       secondaryRates

  real(KIND=benthos_r8), dimension(:,:),intent(in) :: &
       secondarySourceStoich, &
       secondarySinkStoich

  real(KIND=benthos_r8), intent(in) :: &
       dt

  ! !INPUT/OUTPUT PARAMETERS:
  real(KIND=benthos_r8), dimension(nBenthicTracers), intent(inout) :: &
       secondarySourceTendk, &
       secondarySinkTendk

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------------

  integer (KIND=benthos_i4) :: &
       iTracers, &
       iReactions

   secondarySourceTendk(:) = c0_benthos
   secondarySinkTendk(:) = c0_benthos

   do iTracers = 1,nBenthicTracers
      do iReactions = 1,nSecondaryReactions

         if (secondarySourceStoich(iTracers,iReactions) .gt. c0_benthos) then
            secondarySourceTendk(iTracers) = secondarySourceTendk(iTracers) + &

                 secondarySourceStoich(iTracers,iReactions)*secondaryRates(iReactions)*dt

         end if

         if (secondarySinkStoich(iTracers,iReactions) .gt. c0_benthos) then

            secondarySinkTendk(iTracers) = secondarySinkTendk(iTracers) + &
                 secondarySinkStoich(iTracers,iReactions)*secondaryRates(iReactions)*dt
         end if

      end do ! iReactions
   end do ! iTracers

   end subroutine secondaryBenthosReactions

!***********************************************************************
!BOP
! !IROUTINE: benthos_SurfaceFluxes
! !INTERFACE:

 subroutine benthos_SurfaceFluxes(benthos_input, benthos_forcing,   &
                              benthos_flux_diagnostic_fields, &
                              numColumnsMax, numColumns)

! !DESCRIPTION:
!  Compute surface fluxes for tracegas tracer module.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  type(benthos_input_type),   intent(in )   :: benthos_input

  integer (KIND=benthos_i4), intent(in)  :: numColumnsMax, numColumns

! !INPUT/OUTPUT PARAMETERS:

  type(benthos_forcing_type), intent(inout) :: benthos_forcing
  type(benthos_flux_diagnostics_type), intent(inout) :: benthos_flux_diagnostic_fields

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (KIND=benthos_r8) :: &
      seaSurfaceTemp
!-----------------------------------------------------------------------
!EOC

 end subroutine benthos_SurfaceFluxes

  !*****************************************************************************

end module benthos_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
