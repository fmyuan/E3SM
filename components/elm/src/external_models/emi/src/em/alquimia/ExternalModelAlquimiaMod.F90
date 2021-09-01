module ExternalModelAlquimiaMod

  use abortutils                   , only : endrun
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use shr_log_mod                  , only : errMsg => shr_log_errMsg
  use decompMod                    , only : bounds_type
  use EMI_DataMod, only : emi_data_list, emi_data
  use elm_varctl                   , only : iulog

  use ExternalModelBaseType        , only : em_base_type
  use ExternalModelConstants
  use EMI_Atm2LndType_Constants
  use EMI_CanopyStateType_Constants
  use EMI_ColumnType_Constants
  use EMI_EnergyFluxType_Constants
  use EMI_Filter_Constants
  use EMI_Landunit_Constants
  use EMI_SoilHydrologyType_Constants
  use EMI_SoilStateType_Constants
  use EMI_WaterFluxType_Constants
  use EMI_WaterStateType_Constants
  use EMI_CNCarbonStateType_Constants
  use EMI_CNNitrogenStateType_Constants
  use EMI_CNNitrogenFluxType_Constants
  use EMI_CNCarbonFluxType_Constants
  use EMI_ChemStateType_Constants
  use EMI_ColumnWaterStateType_Constants
  use EMI_ColumnWaterFluxType_Constants
  use EMI_ColumnEnergyStateType_Constants

#ifdef USE_ALQUIMIA_LIB
   use AlquimiaContainers_module, only : AlquimiaSizes,AlquimiaProblemMetaData,AlquimiaProperties,&
            AlquimiaState,AlquimiaAuxiliaryData,AlquimiaAuxiliaryOutputData, AlquimiaEngineStatus, &
            AlquimiaEngineFunctionality,AlquimiaGeochemicalCondition
   use AlquimiaContainers_module, only : kAlquimiaMaxStringLength
   use alquimia_fortran_interface_mod, only : AlquimiaFortranInterface
   use iso_c_binding, only : c_ptr
   use c_f_interface_module, only : c_f_string_ptr, f_c_string_ptr
#endif

  use, intrinsic :: iso_c_binding, only : C_CHAR, c_double, c_int, c_bool, c_f_pointer

  implicit none

  type, public, extends(em_base_type) :: em_alquimia_type
    ! Initialization data needed
    integer :: index_l2e_init_filter_soilc
    integer :: index_l2e_init_filter_num_soilc
    integer :: index_l2e_init_state_temperature_soil
    integer :: index_l2e_init_state_h2osoi_liq
    integer :: index_l2e_init_state_h2osoi_ice

    integer :: index_l2e_col_dz
    
    ! Solve data needed
    integer :: index_l2e_state_watsatc ! Porosity
    integer :: index_l2e_filter_soilc
    integer :: index_l2e_filter_num_soilc
    integer :: index_l2e_state_h2osoi_liqvol
    integer :: index_l2e_state_decomp_cpools
    integer :: index_l2e_state_decomp_npools
    integer :: index_l2e_state_temperature_soil
    integer :: index_l2e_soil_pool_decomp_k
    integer :: index_l2e_state_nh4
    integer :: index_l2e_state_no3
    integer :: index_l2e_flux_plantNdemand
    integer :: index_l2e_flux_qflx_adv
    integer :: index_l2e_flux_qflx_lat_aqu_layer
    
    ! Solve data returned to land model
    integer :: index_e2l_state_decomp_cpools
    integer :: index_e2l_state_decomp_npools
    integer :: index_e2l_flux_hr
    integer :: index_e2l_state_nh4
    integer :: index_e2l_state_no3
    integer :: index_e2l_state_DOC
    integer :: index_e2l_state_DON
    integer :: index_e2l_state_DIC

    integer :: index_e2l_state_ph
    integer :: index_e2l_state_salinity
    integer :: index_e2l_state_sulfate
    integer :: index_e2l_state_O2
    integer :: index_e2l_state_Fe2
    integer :: index_e2l_state_FeOxide
    integer :: index_e2l_state_carbonate

    integer :: index_e2l_flux_Nimm
    integer :: index_e2l_flux_Nimp
    integer :: index_e2l_flux_Nmin

    integer :: index_e2l_flux_plantNO3uptake
    integer :: index_e2l_flux_plantNH4uptake

    integer :: index_e2l_flux_NO3runoff
    integer :: index_e2l_flux_DONrunoff
    integer :: index_e2l_flux_DICrunoff
    integer :: index_e2l_flux_DOCrunoff

    ! Alquimia state data gets passed back and forth
    integer :: index_e2l_water_density
    integer :: index_l2e_water_density
    integer :: index_e2l_aqueous_pressure
    integer :: index_l2e_aqueous_pressure
    integer :: index_e2l_total_mobile
    integer :: index_l2e_total_mobile
    integer :: index_e2l_total_immobile
    integer :: index_l2e_total_immobile
    integer :: index_e2l_mineral_volume_fraction
    integer :: index_l2e_mineral_volume_fraction
    integer :: index_e2l_mineral_specific_surface_area
    integer :: index_l2e_mineral_specific_surface_area
    integer :: index_e2l_surface_site_density
    integer :: index_l2e_surface_site_density
    integer :: index_e2l_cation_exchange_capacity
    integer :: index_l2e_cation_exchange_capacity
    integer :: index_e2l_aux_doubles
    integer :: index_l2e_aux_doubles
    integer :: index_e2l_aux_ints
    integer :: index_l2e_aux_ints

    integer :: index_e2l_chem_dt
    
#ifdef USE_ALQUIMIA_LIB
    ! Chemistry engine: Should be one per thread
    type(AlquimiaFortranInterface)       :: chem
    type(AlquimiaEngineStatus)    :: chem_status
    type(c_ptr)                   :: chem_engine
    
    ! Chemistry metadata
    type(AlquimiaSizes)           :: chem_sizes
    type(AlquimiaProblemMetaData) :: chem_metadata
    
    ! Chemical properties and state
    type(AlquimiaProperties) :: chem_properties  ! One copy per processor
    type(AlquimiaState)      :: chem_state       ! Contains a list of species in the structure
    type(AlquimiaAuxiliaryData)   :: chem_aux_data
    type(AlquimiaAuxiliaryOutputData)   :: chem_aux_output
    
    ! Initial condition. Maybe this can just be created and destroyed in a subroutine?
    type(AlquimiaGeochemicalCondition) :: chem_ic
    
#endif

    
    ! Mapping between ELM and alquimia decomp pools
    integer, pointer, dimension(:)       :: carbon_pool_mapping
    integer, pointer, dimension(:)       :: nitrogen_pool_mapping
    integer, pointer, dimension(:)       :: pool_reaction_mapping
    integer                              :: CO2_pool_number
    integer                              :: NH4_pool_number,NO3_pool_number
    integer                              :: Nimm_pool_number,Nmin_pool_number,Nimp_pool_number
    integer                              :: plantNO3uptake_pool_number,plantNH4uptake_pool_number
    integer                              :: plantNO3demand_pool_number,plantNH4demand_pool_number
    integer                              :: plantNO3uptake_reaction_number,plantNH4uptake_reaction_number
    integer                              :: Hplus_pool_number,sulfate_pool_number,O2_pool_number,chloride_pool_number,Fe2_pool_number,FeOH3_pool_number
    logical, pointer, dimension(:)       :: is_dissolved_gas
    real(r8),pointer,dimension(:)        :: DOC_content,DIC_content,DON_content,carbonate_C_content ! Also add extra SOM content tracker for pools beyond ELM's litter and SOM?
    real(r8),pointer,dimension(:)        :: bc ! Boundary condition (len of chem_sizes%num_primary)
    
   contains
     procedure, public :: Populate_L2E_Init_List  => EMAlquimia_Populate_L2E_Init_List
     procedure, public :: Populate_E2L_Init_List  => EMAlquimia_Populate_E2L_Init_List
     procedure, public :: Populate_L2E_List       => EMAlquimia_Populate_L2E_List
     procedure, public :: Populate_E2L_List       => EMAlquimia_Populate_E2L_List
     procedure, public :: Init                    => EMAlquimia_Init
     procedure, public :: Solve                   => EMAlquimia_Solve
#ifdef USE_ALQUIMIA_LIB
     procedure, private :: Copy_Alquimia_To_ELM
     procedure, private :: Copy_ELM_To_Alquimia
     procedure, private :: map_alquimia_pools
#endif
  end type em_alquimia_type


  real(r8),parameter :: min_dt = 1.0 ! Minimum time step length(s) before crashing model on non-convergence in ReactionStepOperatorSplit
#ifndef USE_ALQUIMIA_LIB
  integer, parameter :: kAlquimiaMaxStringLength = 512
#endif

contains

  !------------------------------------------------------------------------
  subroutine EMAlquimia_Populate_L2E_Init_List(this, l2e_init_list)
    !
    ! !DESCRIPTION:
    ! Initialze an emi_list for exchanging data from land model to external
    ! model during initialization stage
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_Alquimia_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_init_list
    
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index
    
    ! number_em_stages = 1
    ! allocate(em_stages(number_em_stages))
    ! em_stages(1) = EM_INITIALIZATION_STAGE


    ! deallocate(em_stages)
    
    write(iulog,*)'L2EInit List:'
    call l2e_init_list%PrintInfo()

  end subroutine EMAlquimia_Populate_L2E_Init_List

  !------------------------------------------------------------------------
  subroutine EMAlquimia_Populate_E2L_Init_List(this, e2l_init_list)
    !
    ! !DESCRIPTION:
    ! Initialze an emi_list for exchanging data from external model to land
    ! model during initialization stage
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_Alquimia_type)                 :: this
    class(emi_data_list), intent(inout) :: e2l_init_list

    ! write(iulog,*)'EMAlquimia_Populate_E2L_Init_List must be extended by a child class.'
    ! call endrun(msg=errMsg(__FILE__, __LINE__))
    ! write(iulog,*)'EMAlquimia_Populate_E2L_Init_List is empty.'
    write(iulog,*)'E2LInit List:'
    call e2l_init_list%PrintInfo()

  end subroutine EMAlquimia_Populate_E2L_Init_List

  !------------------------------------------------------------------------
  subroutine EMAlquimia_Populate_L2E_List(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Initialze an emi_list for exchanging data from land model to external
    ! model during time integration stage
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_Alquimia_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_list
    
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_Alquimia_SOLVE_STAGE



    ! Liquid water
    id                                   = L2E_STATE_SOIL_LIQ_VOL_COL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_liqvol      = index
    
    ! Carbon pools
    id                                             = L2E_STATE_CARBON_POOLS_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_decomp_cpools              = index
    
    ! Nitrogen pools
    id                                             = L2E_STATE_NITROGEN_POOLS_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_decomp_npools              = index

    ! Soil temperature
    id                                             = L2E_STATE_TSOIL_NLEVSOI_COL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_temperature_soil              = index
    
    ! Decomposition rate constants
    id                                             = L2E_FLUX_SOIL_POOL_DECOMP_K
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_soil_pool_decomp_k              = index
    
    id                                             = L2E_STATE_NH4_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_nh4              = index
    
    id                                             = L2E_STATE_NO3_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_no3              = index

    id                                             = L2E_FLUX_PLANT_NDEMAND_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_plantNdemand              = index


    ! Alquimia data is sent from ELM to Alquimia only at solve stage (not set yet at cold start stage)

    id                                             = L2E_STATE_WATER_DENSITY
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_water_density      = index

    id                                             = L2E_STATE_AQUEOUS_PRESSURE
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_aqueous_pressure      = index

    id                                             = L2E_STATE_TOTAL_IMMOBILE
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_total_immobile      = index

    id                                             = L2E_STATE_TOTAL_MOBILE
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_total_mobile      = index

    id                                             = L2E_STATE_MINERAL_VOLUME_FRACTION
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_mineral_volume_fraction      = index

    id                                             = L2E_STATE_MINERAL_SPECIFIC_SURFACE_AREA
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_mineral_specific_surface_area      = index

    id                                             = L2E_STATE_SURFACE_SITE_DENSITY
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_surface_site_density      = index

    id                                             = L2E_STATE_CATION_EXCHANGE_CAPACITY
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_cation_exchange_capacity     = index

    id                                             = L2E_STATE_AUX_DOUBLES
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_aux_doubles      = index

    id                                             = L2E_STATE_AUX_INTS
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_aux_ints      = index

    ! Water flow
    id                                   = L2E_FLUX_SOIL_QFLX_ADV_COL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_qflx_adv      = index

    id                                   = L2E_FLUX_SOIL_QFLX_LAT_COL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_qflx_lat_aqu_layer      = index



    ! Needed for both stages
    deallocate(em_stages)
    number_em_stages = 2
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_ALQUIMIA_SOLVE_STAGE
    em_stages(2) = EM_ALQUIMIA_COLDSTART_STAGE

    id                                             = L2E_PARAMETER_WATSATC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_watsatc              = index

    id                                             = L2E_COLUMN_DZ
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_dz              = index

    id                                   = L2E_FILTER_SOILC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_soilc     = index

    id                                   = L2E_FILTER_NUM_SOILC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_num_soilc = index


    deallocate(em_stages)
    
    write(iulog,*)'L2E List:'
    call l2e_list%PrintInfo()

  end subroutine EMAlquimia_Populate_L2E_List

  !------------------------------------------------------------------------
  subroutine EMAlquimia_Populate_E2L_List(this, e2l_list)
    !
    ! !DESCRIPTION:
    ! Initialze an emi_list for exchanging data from external model to land
    ! model during time integration stage
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_Alquimia_type)                 :: this
    class(emi_data_list), intent(inout) :: e2l_list
    
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    ! Updated Carbon pools
    ! May want to change this to rates of change instead?
    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_ALQUIMIA_SOLVE_STAGE
    
    id                                             = E2L_STATE_CARBON_POOLS_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_decomp_cpools              = index
    
    ! Nitrogen pools
    id                                             = E2L_STATE_NITROGEN_POOLS_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_decomp_npools              = index

    ! Heterotrophic respiration flux
    id                                             = E2L_FLUX_HETEROTROPHIC_RESP!_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_hr              = index
    
    id                                             = E2L_STATE_NH4_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_nh4              = index
    
    id                                             = E2L_STATE_NO3_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_no3              = index

    id                                             = E2L_STATE_DOC_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_DOC              = index

    id                                             = E2L_STATE_DON_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_DON              = index

    id                                             = E2L_STATE_DIC_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_DIC              = index

    id                                             = E2L_FLUX_NIMM_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_Nimm              = index
    
    id                                             = E2L_FLUX_NIMP_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_Nimp              = index

    id                                             = E2L_FLUX_NMIN_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_Nmin              = index

    id                                             = E2L_FLUX_SMIN_NO3_TO_PLANT_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_plantNO3uptake      = index

    id                                             = E2L_FLUX_SMIN_NH4_TO_PLANT_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_plantNH4uptake      = index

    id                                             = E2L_FLUX_NO3_RUNOFF
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_NO3runoff      = index

    id                                             = E2L_FLUX_DON_RUNOFF
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_DONrunoff      = index

    id                                             = E2L_FLUX_DIC_RUNOFF
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_DICrunoff      = index

    id                                             = E2L_FLUX_DOC_RUNOFF
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_DOCrunoff      = index

    id                                             = E2L_STATE_SOIL_PH
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_ph              = index

    id                                             = E2L_STATE_SOIL_SALINITY
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_salinity             = index

    id                                             = E2L_STATE_SOIL_O2
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_O2              = index

    id                                             = E2L_STATE_SOIL_SULFATE
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_sulfate              = index

    id                                             = E2L_STATE_SOIL_FE2
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_Fe2              = index

    id                                             = E2L_STATE_SOIL_FE_OXIDE
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_FeOxide              = index

    id                                             = E2L_STATE_SOIL_CARBONATE
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_carbonate              = index

    id                                             = E2L_STATE_CHEM_DT
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_chem_dt      = index

    ! These need to be exchanged in both stages
    deallocate(em_stages)
    number_em_stages = 2
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_ALQUIMIA_SOLVE_STAGE
    em_stages(2) = EM_ALQUIMIA_COLDSTART_STAGE

    id                                             = E2L_STATE_WATER_DENSITY
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_water_density      = index

    id                                             = E2L_STATE_AQUEOUS_PRESSURE
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_aqueous_pressure      = index

    id                                             = E2L_STATE_TOTAL_IMMOBILE
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_total_immobile      = index

    id                                             = E2L_STATE_TOTAL_MOBILE
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_total_mobile      = index

    id                                             = E2L_STATE_MINERAL_VOLUME_FRACTION
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_mineral_volume_fraction      = index

    id                                             = E2L_STATE_MINERAL_SPECIFIC_SURFACE_AREA
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_mineral_specific_surface_area      = index

    id                                             = E2L_STATE_SURFACE_SITE_DENSITY
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_surface_site_density      = index

    id                                             = E2L_STATE_CATION_EXCHANGE_CAPACITY
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_cation_exchange_capacity     = index

    id                                             = E2L_STATE_AUX_DOUBLES
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_aux_doubles      = index

    id                                             = E2L_STATE_AUX_INTS
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_aux_ints      = index

    deallocate(em_stages)

    write(iulog,*)'E2L List:'
    call e2l_list%PrintInfo()

  end subroutine EMAlquimia_Populate_E2L_List

  !------------------------------------------------------------------------
  subroutine EMAlquimia_Init(this, l2e_init_list, e2l_init_list, iam, bounds_clump)
    !
    ! !DESCRIPTION:
    ! Initialize an emi_list for exchanging data from land model to external
    ! model during time integration stage
    !
#ifdef USE_ALQUIMIA_LIB
    use alquimia_fortran_interface_mod, only : AllocateAlquimiaEngineStatus, &
                                            AllocateAlquimiaProblemMetaData,&
                                            AllocateAlquimiaState,&
                                            AllocateAlquimiaProperties,&
                                            AllocateAlquimiaAuxiliaryData,&
                                            AllocateAlquimiaAuxiliaryOutputData, &
                                            AllocateAlquimiaGeochemicalCondition
                                            
    use elm_varctl, only : alquimia_inputfile,alquimia_engine_name,alquimia_IC_name,alquimia_handsoff

    use elm_varpar            , only : alquimia_num_primary, alquimia_num_minerals,&
                                       alquimia_num_surface_sites, alquimia_num_ion_exchange_sites, &
                                       alquimia_num_aux_doubles, alquimia_num_aux_ints
    use landunit_varcon, only : istcrop,istsoil

    use PFloTranAlquimiaInterface_module, only : PrintSizes,PrintProblemMetaData, ProcessCondition,PrintState

    implicit none
    !
    ! !ARGUMENTS
    class(em_Alquimia_type)                  :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    integer              , intent(in)    :: iam
    type(bounds_type)    , intent (in)   :: bounds_clump
    
    
    ! Local variables

    
    
    ! Should read this from a namelist
    character(kind=C_CHAR,len=kAlquimiaMaxStringLength) :: inputfile
    character(kind=C_CHAR,len=kAlquimiaMaxStringLength) :: engine_name
    character(kind=C_CHAR,len=kAlquimiaMaxStringLength) :: IC_name
    
    
    logical(C_BOOL) :: hands_off
    character(kind=C_CHAR,len=kAlquimiaMaxStringLength) :: status_message
    type(AlquimiaEngineFunctionality) :: chem_engine_functionality

    
    write(iulog,*), 'Entering Alquimia setup'
    
    inputfile   = alquimia_inputfile
    engine_name = alquimia_engine_name
    IC_name     = alquimia_IC_name  ! Name of initial condition
    hands_off = alquimia_handsoff  ! hands_off = .false. allows/requires rate constants, mineral rate const, CEC, complexation site density, and isotherms to be passed through alquimia 
    
    ! Make sure these are not defined until explicitly set
    this%carbon_pool_mapping   => NULL()
    this%nitrogen_pool_mapping => NULL()
    this%pool_reaction_mapping => NULL()
    this%is_dissolved_gas      => NULL()
    this%bc                    => NULL()

    ! Allocate memory for status container
    call AllocateAlquimiaEngineStatus(this%chem_status)
    ! Point Alquimia interface to correct subroutines (based on engine that was specified in engine_name)
    call this%chem%CreateInterface(engine_name, this%chem_status)
    
    ! Print out the result of the interface creation call
    call c_f_string_ptr(this%chem_status%message,status_message)
    if(this%chem_status%error /= 0) then
      call endrun(msg='Alquimia error: '//status_message)
    endif
    
    ! Set up the engine and get the storage requirements
    ! Should this only happen on one processor and then broadcast?
    call this%chem%Setup(inputfile, hands_off, this%chem_engine, this%chem_sizes, chem_engine_functionality, this%chem_status)
    ! Print out the result of the interface creation call
    call c_f_string_ptr(this%chem_status%message,status_message)
    if(this%chem_status%error /= 0) then
      call endrun(msg='Alquimia error: '//status_message)
    endif
    
    ! Copy array sizes over to clm_varpar
    ! EMI xml system doesn't seem to allow single integers to transferred very easily so we are writing directly to clm_varpar
    alquimia_num_primary            = this%chem_sizes%num_primary
    alquimia_num_minerals           = this%chem_sizes%num_minerals
    alquimia_num_surface_sites      = this%chem_sizes%num_surface_sites
    alquimia_num_ion_exchange_sites = this%chem_sizes%num_ion_exchange_sites
    alquimia_num_aux_doubles        = this%chem_sizes%num_aux_doubles
    alquimia_num_aux_ints           = this%chem_sizes%num_aux_integers
    
    ! Allocate memory for chemistry data
    call AllocateAlquimiaProblemMetaData(this%chem_sizes, this%chem_metadata)
    
    call this%chem%GetProblemMetaData(this%chem_engine, this%chem_metadata, this%chem_status)
    if(this%chem_status%error /= 0) then
      call c_f_string_ptr(this%chem_status%message,status_message)
      call endrun(msg='Alquimia error: '//status_message)
    endif
    
    ! Transfer metadata back to ELM? Does EMI allow character data transfers?
    call printproblemmetadata(this%chem_metadata)
    

    ! Initial condition. The zero length for constraints suggest that it must be read in from input file
    ! In principle the input deck could also include constraints for upper boundary condition and lateral boundary conditions (saltwater, freshwater?)
    ! but that could get tricky if those conditions are not constant over time. Would we have to reprocess the BC every time step?
    ! I think we do need some kind of alquimia condition for boundaries because ELM/MOSART/etc won't necessarily have the same chemicals as the alquimia reaction network
    call AllocateAlquimiaGeochemicalCondition(len_trim(ic_name,C_INT),0,0,this%chem_ic)
    call f_c_string_ptr(ic_name,this%chem_ic%name,len_trim(ic_name)+1)
    
    ! Allocate alquimia's data structures. One copy per processor which will be written into as needed
    call AllocateAlquimiaState(this%chem_sizes, this%chem_state)
    call AllocateAlquimiaProperties(this%chem_sizes, this%chem_properties)
    call AllocateAlquimiaAuxiliaryData(this%chem_sizes, this%chem_aux_data)
    call AllocateAlquimiaAuxiliaryOutputData(this%chem_sizes, this%chem_aux_output)

    allocate(this%bc(this%chem_sizes%num_primary))
    

#else
  implicit none
  !
  ! !ARGUMENTS
  class(em_Alquimia_type)                  :: this
  class(emi_data_list) , intent(in)    :: l2e_init_list
  class(emi_data_list) , intent(inout) :: e2l_init_list
  integer              , intent(in)    :: iam
  type(bounds_type)    , intent (in)   :: bounds_clump

  call endrun(msg='ERROR: Attempting to run with alquimia when model not compiled with USE_ALQUIMIA_LIB')
#endif

  end subroutine EMAlquimia_Init


      !------------------------------------------------------------------------
  subroutine EMAlquimia_Solve(this, em_stage, dt, nstep, clump_rank, l2e_list, e2l_list, &
    bounds_clump)
!
! !DESCRIPTION:
! 
!
! !USES:
use shr_kind_mod              , only : r8 => shr_kind_r8
use abortutils                , only : endrun
use shr_log_mod               , only : errMsg => shr_log_errMsg
use elm_varctl                , only : iulog
use ExternalModelConstants    , only : EM_ALQUIMIA_SOLVE_STAGE,EM_ALQUIMIA_COLDSTART_STAGE
                                      
!
implicit none
!
! !ARGUMENTS:
class(em_alquimia_type)              :: this
integer              , intent(in)    :: em_stage
real(r8)             , intent(in)    :: dt
integer              , intent(in)    :: nstep
integer              , intent(in)    :: clump_rank
class(emi_data_list) , intent(in)    :: l2e_list
class(emi_data_list) , intent(inout) :: e2l_list
type(bounds_type)    , intent (in)   :: bounds_clump

select case(em_stage)

case (EM_ALQUIMIA_SOLVE_STAGE)
   call EMAlquimia_Solve_BGC(this, dt, nstep, clump_rank, l2e_list, e2l_list, &
        bounds_clump)


case (EM_ALQUIMIA_COLDSTART_STAGE)
  call EMAlquimia_Coldstart(this, clump_rank, l2e_list, e2l_list, bounds_clump)

case default
   write(iulog,*)'EM_Alquimia_Solve: Unknown em_stage.'
   call endrun(msg=errMsg(__FILE__, __LINE__))
end select

end subroutine EMAlquimia_Solve


subroutine EMAlquimia_Coldstart(this, clump_rank, l2e_list, e2l_list, bounds_clump)

  use elm_varpar, only : nlevdecomp

  class(em_alquimia_type)              :: this
  integer              , intent(in)    :: clump_rank
  class(emi_data_list) , intent(in)    :: l2e_list
  class(emi_data_list) , intent(inout) :: e2l_list
  type(bounds_type)    , intent (in)   :: bounds_clump

#ifdef USE_ALQUIMIA_LIB

  real(r8) , pointer, dimension(:,:)   ::  porosity_l2e, dz, h2o_liqvol
  real(r8) , pointer, dimension(:,:)   ::  water_density_e2l,aqueous_pressure_e2l
  real(r8) , pointer, dimension(:,:,:) ::  total_mobile_e2l
  real(r8) , pointer, dimension(:,:,:) ::  total_immobile_e2l
  real(r8) , pointer, dimension(:,:,:) ::  mineral_volume_fraction_e2l
  real(r8) , pointer, dimension(:,:,:) ::  mineral_specific_surface_area_e2l
  real(r8) , pointer, dimension(:,:,:) ::  surface_site_density_e2l
  real(r8) , pointer, dimension(:,:,:) ::  cation_exchange_capacity_e2l
  real(r8) , pointer, dimension(:,:,:) ::  aux_doubles_e2l
  real(r8) , dimension(nlevdecomp, this%chem_sizes%num_primary) :: free_mobile
  integer  , pointer, dimension(:,:,:)   ::  aux_ints_e2l
  integer   , pointer                  :: filter_soilc(:)

  integer :: c, fc, j, num_soilc
  character(kind=C_CHAR,len=kAlquimiaMaxStringLength) :: status_message

  call l2e_list%GetPointerToInt1D(this%index_l2e_filter_soilc , filter_soilc   )
  call l2e_list%GetIntValue(this%index_l2e_filter_num_soilc          , num_soilc   )
  call l2e_list%GetPointerToReal2D(this%index_l2e_col_dz, dz)  

  call l2e_list%GetPointerToReal2D(this%index_l2e_state_watsatc       , porosity_l2e     )
  call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_liqvol, h2o_liqvol) ! m3/m3

  ! Alquimia state data to set on ELM side
  call e2l_list%GetPointerToReal2D(this%index_e2l_water_density, water_density_e2l)
  call e2l_list%GetPointerToReal2D(this%index_e2l_aqueous_pressure, aqueous_pressure_e2l)
  call e2l_list%GetPointerToReal3D(this%index_e2l_total_mobile, total_mobile_e2l)
  call e2l_list%GetPointerToReal3D(this%index_e2l_total_immobile, total_immobile_e2l)
  call e2l_list%GetPointerToReal3D(this%index_e2l_mineral_volume_fraction, mineral_volume_fraction_e2l)
  call e2l_list%GetPointerToReal3D(this%index_e2l_mineral_specific_surface_area, mineral_specific_surface_area_e2l)
  call e2l_list%GetPointerToReal3D(this%index_e2l_surface_site_density, surface_site_density_e2l)
  call e2l_list%GetPointerToReal3D(this%index_e2l_cation_exchange_capacity, cation_exchange_capacity_e2l)
  call e2l_list%GetPointerToReal3D(this%index_e2l_aux_doubles, aux_doubles_e2l)
  call e2l_list%GetPointerToInt3D(this%index_e2l_aux_ints, aux_ints_e2l)

  do fc = 1, num_soilc
    c = filter_soilc(fc)

      do j = 1, nlevdecomp

          
          ! Initialize the state for the cell
          this%chem_properties%volume = dz(c,j)
          this%chem_properties%saturation = 0.5_r8 ! h2o_liqvol(c,j)/porosity_l2e(c,j)
          this%chem_state%water_density = 1.0e3_r8
          this%chem_state%porosity = porosity_l2e(c,j)
          this%chem_state%aqueous_pressure = 101325.0
          this%chem_state%temperature = 20.0_r8 ! Temperature may not have been initialized yet

          call this%chem%ProcessCondition(this%chem_engine, this%chem_ic, this%chem_properties, this%chem_state, &
                                         this%chem_aux_data, this%chem_status)
          if(this%chem_status%error /= 0) then
            call c_f_string_ptr(this%chem_status%message,status_message)
            call endrun(msg='Alquimia error in ProcessCondition: '//status_message)
          endif

          this%chem_state%porosity = porosity_l2e(c,j)
          ! But this can only happen after ELM allocation step, so this whole thing might need to move somewhere else
          call this%copy_Alquimia_to_ELM(c,j,water_density_e2l,&
                                        aqueous_pressure_e2l,&
                                        total_mobile_e2l,free_mobile,&
                                        total_immobile_e2l,&
                                        mineral_volume_fraction_e2l,&
                                        mineral_specific_surface_area_e2l,&
                                        surface_site_density_e2l,&
                                        cation_exchange_capacity_e2l,&
                                        aux_doubles_e2l,&
                                        aux_ints_e2l) 
          
      enddo
  enddo
  ! Save condition to use as surface boundary condition. Units here are converted back to mol/m3 H2O
  ! Note: Boundary condition also needs to be set (or saved/read) when initializing from restart
  this%bc(1:this%chem_sizes%num_primary) = total_mobile_e2l(c,1,1:this%chem_sizes%num_primary)/(porosity_l2e(c,1)*0.5_r8)
#endif
end subroutine EMAlquimia_Coldstart

  !------------------------------------------------------------------------
  subroutine EMAlquimia_Solve_BGC(this, dt, nstep, clump_rank, l2e_list, e2l_list, &
       bounds_clump)

    
#ifdef USE_ALQUIMIA_LIB

    use elm_varpar, only : nlevdecomp,ndecomp_pools
    use landunit_varcon, only : istcrop,istsoil
    ! use clm_varcon, only : catomw,natomw ! Replacing these with constants that are the same as PFLOTRAN defs
    use AlquimiaContainers_module, only : AlquimiaEngineStatus
    use alquimia_fortran_interface_mod, only :  ReactionStepOperatorSplit, GetAuxiliaryOutput
    use PFloTranAlquimiaInterface_module, only : printState
    
    use CNDecompCascadeConType, only : decomp_cascade_con

    implicit none
    !
    ! !ARGUMENTS
    class(em_alquimia_type)              :: this
    real(r8)             , intent(in)    :: dt ! s
    integer              , intent(in)    :: nstep
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent (in)   :: bounds_clump
    
    
    ! Local variables
    integer                              :: c,fc,j,k,poolnum
    integer   , pointer                  :: filter_soilc(:)
    integer                              :: num_soilc
    integer                              :: max_cuts
    real(r8) , pointer, dimension(:,:,:)    :: soilcarbon_l2e,soilcarbon_e2l 
    real(r8) , pointer, dimension(:,:,:)    :: soilnitrogen_l2e,soilnitrogen_e2l 
    real(r8) , pointer, dimension(:,:,:)    :: decomp_k
    real(r8) , pointer, dimension(:,:)    :: temperature, h2o_liqvol
    real(r8) , pointer, dimension(:)     :: hr_e2l ! 1D total surface emission
    real(r8) , pointer, dimension(:)     :: NO3runoff_e2l,DONrunoff_e2l ! 1D total column runoff (gN/m2/s)
    real(r8) , pointer, dimension(:)     :: DICrunoff_e2l,DOCrunoff_e2l ! 1D total column runoff (gN/m2/s)
    real(r8) , pointer, dimension(:,:)  :: no3_e2l,no3_l2e,nh4_e2l,nh4_l2e
    real(r8) , pointer, dimension(:,:)  :: Nimm_e2l, Nimp_e2l, Nmin_e2l
    real(r8) , pointer, dimension(:,:)  :: plantNO3uptake_e2l,plantNH4uptake_e2l, plantNdemand_l2e
    real(r8) , pointer, dimension(:,:)  :: water_density_l2e,water_density_e2l,aqueous_pressure_l2e,aqueous_pressure_e2l,porosity_l2e,dz
    real(r8) , pointer, dimension(:,:,:) :: total_mobile_l2e , total_mobile_e2l
    real(r8) , pointer, dimension(:,:,:) :: total_immobile_l2e , total_immobile_e2l
    real(r8) , pointer, dimension(:,:,:) :: mineral_volume_fraction_l2e , mineral_volume_fraction_e2l
    real(r8) , pointer, dimension(:,:,:) :: mineral_specific_surface_area_l2e , mineral_specific_surface_area_e2l
    real(r8) , pointer, dimension(:,:,:) :: surface_site_density_l2e , surface_site_density_e2l
    real(r8) , pointer, dimension(:,:,:) :: cation_exchange_capacity_l2e , cation_exchange_capacity_e2l
    real(r8) , pointer, dimension(:,:,:) :: aux_doubles_l2e , aux_doubles_e2l
    integer  , pointer, dimension(:,:,:)   :: aux_ints_l2e, aux_ints_e2l
    real(r8) , pointer, dimension(:,:)    :: qflx_adv_l2e, qflx_lat_aqu_l2e
    real(r8) , pointer, dimension(:,:)    :: DOC_e2l, DON_e2l, DIC_e2l
    real(r8) , pointer, dimension(:,:)    :: pH_e2l, O2_e2l, salinity_e2l, sulfate_e2l, Fe2_e2l, FeOxide_e2l, carbonate_e2l
    real(r8) , pointer, dimension(:)     :: actual_dt_e2l
    real(r8)                            :: CO2_before, molperL_to_molperm3
    real(r8), parameter                 :: minval = 1.e-30_r8 ! Minimum value to pass to PFLOTRAN to avoid numerical errors with concentrations of 0

    ! Setting these to the values in PFLOTRAN clm_rspfuncs.F90
    real(r8), parameter :: natomw = 14.0067d0 ! Value in clmvarcon is 14.007
    real(r8), parameter :: catomw = 12.0110d0 ! Value in clmvarcon is 12.011
    real(r8),dimension(this%chem_sizes%num_primary)  :: surf_flux, surf_bc, lat_flux, lat_bc
    real(r8),dimension(nlevdecomp,this%chem_sizes%num_primary) :: free_mobile
    
    character(kind=C_CHAR,len=kAlquimiaMaxStringLength) :: status_message
    procedure(ReactionStepOperatorSplit), pointer :: engine_ReactionStepOperatorSplit
    procedure(GetAuxiliaryOutput), pointer   :: engine_getAuxiliaryOutput
    ! real (c_double), pointer :: alquimia_mobile_data(:), alquimia_immobile_data(:), alquimia_rates_data(:)

    ! write(iulog,*) 'Alquimia solving step!'
    
    ! Pass data from ELM
    
    ! Column filters
    call l2e_list%GetPointerToInt1D(this%index_l2e_filter_soilc , filter_soilc   )
    call l2e_list%GetIntValue(this%index_l2e_filter_num_soilc          , num_soilc   )

    call l2e_list%GetPointerToReal2D(this%index_l2e_col_dz, dz)
    
    ! C and N pools. Units: gC/m2, gN/m2
    call l2e_list%GetPointerToReal3D(this%index_l2e_state_decomp_cpools , soilcarbon_l2e)
    call l2e_list%GetPointerToReal3D(this%index_l2e_state_decomp_npools , soilnitrogen_l2e)
    
    ! (gN/m3)
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_no3 , no3_l2e)
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_nh4 , nh4_l2e)

    ! Abiotic factors
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_temperature_soil , temperature  ) ! K
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_liqvol, h2o_liqvol) ! m3/m3
    ! call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_ice, h2o_ice) ! kg/m2
    
    ! Pool turnover rate constants calculated in ELM, incorporating T and moisture effects (1/s)
    call l2e_list%GetPointerToReal3D(this%index_l2e_soil_pool_decomp_k, decomp_k)

    call l2e_list%GetPointerToReal2D(this%index_l2e_flux_plantNdemand, plantNdemand_l2e)

    ! C and N pools
    call e2l_list%GetPointerToReal3D(this%index_e2l_state_decomp_cpools , soilcarbon_e2l) ! gC/m2
    call e2l_list%GetPointerToReal3D(this%index_e2l_state_decomp_npools , soilnitrogen_e2l) ! gN/m2
    ! call e2l_list%GetPointerToReal2D(this%index_e2l_flux_hr , hr_e2l) ! (gC/m3/s)
    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_hr , hr_e2l) ! (gC/m2/s)
    
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_no3 , no3_e2l) ! gN/m3
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_nh4 , nh4_e2l) ! gN/m3

    call e2l_list%GetPointerToReal2D(this%index_e2l_flux_Nimm , Nimm_e2l) ! gN/m3/s
    call e2l_list%GetPointerToReal2D(this%index_e2l_flux_Nimp , Nimp_e2l) ! gN/m3/s
    call e2l_list%GetPointerToReal2D(this%index_e2l_flux_Nmin , Nmin_e2l) ! gN/m3/s

    call e2l_list%GetPointerToReal2D(this%index_e2l_flux_plantNO3uptake , plantNO3uptake_e2l) ! gN/m3/s
    call e2l_list%GetPointerToReal2D(this%index_e2l_flux_plantNH4uptake , plantNH4uptake_e2l) ! gN/m3/s

    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_NO3runoff , NO3runoff_e2l) ! gN/m2/s
    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_DONrunoff , DONrunoff_e2l) ! gN/m2/s
    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_DOCrunoff , DOCrunoff_e2l) ! gC/m2/s
    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_DICrunoff , DICrunoff_e2l) ! gC/m2/s

    ! Alquimia state data on ELM side
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_watsatc       , porosity_l2e     )
    call e2l_list%GetPointerToReal2D(this%index_e2l_water_density, water_density_e2l)
    call l2e_list%GetPointerToReal2D(this%index_l2e_water_density, water_density_l2e)
    call e2l_list%GetPointerToReal2D(this%index_e2l_aqueous_pressure, aqueous_pressure_e2l)
    call l2e_list%GetPointerToReal2D(this%index_l2e_aqueous_pressure, aqueous_pressure_l2e)
    call e2l_list%GetPointerToReal3D(this%index_e2l_total_mobile, total_mobile_e2l) ! Note total mobile is stored as mol/m3 bulk and only converted to mol/L water when passed to/from alquimia
    call l2e_list%GetPointerToReal3D(this%index_l2e_total_mobile, total_mobile_l2e)
    call e2l_list%GetPointerToReal3D(this%index_e2l_total_immobile, total_immobile_e2l)
    call l2e_list%GetPointerToReal3D(this%index_l2e_total_immobile, total_immobile_l2e)
    call e2l_list%GetPointerToReal3D(this%index_e2l_mineral_volume_fraction, mineral_volume_fraction_e2l)
    call l2e_list%GetPointerToReal3D(this%index_l2e_mineral_volume_fraction, mineral_volume_fraction_l2e)
    call e2l_list%GetPointerToReal3D(this%index_e2l_mineral_specific_surface_area, mineral_specific_surface_area_e2l)
    call l2e_list%GetPointerToReal3D(this%index_l2e_mineral_specific_surface_area, mineral_specific_surface_area_l2e)
    call e2l_list%GetPointerToReal3D(this%index_e2l_surface_site_density, surface_site_density_e2l)
    call l2e_list%GetPointerToReal3D(this%index_l2e_surface_site_density, surface_site_density_l2e)
    call e2l_list%GetPointerToReal3D(this%index_e2l_cation_exchange_capacity, cation_exchange_capacity_e2l)
    call l2e_list%GetPointerToReal3D(this%index_l2e_cation_exchange_capacity, cation_exchange_capacity_l2e)
    call e2l_list%GetPointerToReal3D(this%index_e2l_aux_doubles, aux_doubles_e2l)
    call l2e_list%GetPointerToReal3D(this%index_l2e_aux_doubles, aux_doubles_l2e)
    call e2l_list%GetPointerToInt3D(this%index_e2l_aux_ints, aux_ints_e2l)
    call l2e_list%GetPointerToInt3D(this%index_l2e_aux_ints, aux_ints_l2e)

    call l2e_list%GetPointerToReal2D(this%index_l2e_flux_qflx_adv       , qflx_adv_l2e     )
    call l2e_list%GetPointerToReal2D(this%index_l2e_flux_qflx_lat_aqu_layer       , qflx_lat_aqu_l2e     )

    call e2l_list%GetPointerToReal2D(this%index_e2l_state_DIC , DIC_e2l)
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_DOC , DOC_e2l)
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_DON , DON_e2l)

    call e2l_list%GetPointerToReal2D(this%index_e2l_state_pH , pH_e2l)
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_salinity , salinity_e2l)
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_O2 , O2_e2l)
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_sulfate , sulfate_e2l)
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_Fe2 , Fe2_e2l)
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_FeOxide , FeOxide_e2l)
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_carbonate , carbonate_e2l)

    call e2l_list%GetPointerToReal1D(this%index_e2l_chem_dt , actual_dt_e2l)

    ! First check if pools have been mapped between ELM and Alquimia
    if(.not. associated(this%carbon_pool_mapping)) then
      call this%map_alquimia_pools()

      c = filter_soilc(1)

      ! At this point, also make sure boundary condition is set
      ! Initialize the state for the cell
      this%chem_properties%volume = dz(c,1)
      this%chem_properties%saturation = 0.5_r8 ! h2o_liqvol(c,j)/porosity_l2e(c,j)
      this%chem_state%water_density = 1.0e3_r8
      this%chem_state%porosity = porosity_l2e(c,1)
      this%chem_state%aqueous_pressure = 101325.0
      this%chem_state%temperature = 20.0_r8 ! Temperature may not have been initialized yet

      call this%chem%ProcessCondition(this%chem_engine, this%chem_ic, this%chem_properties, this%chem_state, &
                                     this%chem_aux_data, this%chem_status)
      if(this%chem_status%error /= 0) then
        call c_f_string_ptr(this%chem_status%message,status_message)
        call endrun(msg='Alquimia error in ProcessCondition: '//status_message)
      endif

      this%chem_state%porosity = porosity_l2e(c,1)
      if(isnan(porosity_l2e(c,1)) .or. porosity_l2e(c,1)<0.01) then
        write(iulog,*),'Porosity = ',porosity_l2e(c,1)
        write(iulog,*),'Boundary condition: ',this%bc(1:this%chem_sizes%num_primary)
        call endrun(msg='Alquimia: Problem with porosity in boundary condition')
      endif

      ! Does this need to loop over all columns?
      call this%copy_Alquimia_to_ELM(c,1,water_density_e2l,&
                aqueous_pressure_e2l,&
                total_mobile_e2l,free_mobile,&
                total_immobile_e2l,&
                mineral_volume_fraction_e2l,&
                mineral_specific_surface_area_e2l,&
                surface_site_density_e2l,&
                cation_exchange_capacity_e2l,&
                aux_doubles_e2l,&
                aux_ints_e2l) 

      this%bc(1:this%chem_sizes%num_primary) = total_mobile_e2l(c,1,1:this%chem_sizes%num_primary)/(porosity_l2e(c,1)*0.5_r8)
      write(iulog,*),'Initialized boundary condition: ',this%bc(1:this%chem_sizes%num_primary)

    endif
    
     ! Run the reactions engine for a step. Alquimia works on one cell at a time
     ! TODO: Transport needs to be integrated somehow. 
    do fc = 1, num_soilc
      c = filter_soilc(fc)

         do j = 1, nlevdecomp  

             ! Set soil carbon and nitrogen from land model
             ! Convert soil C,N from g/m3 to mol/m3. Assumes pool is defined as immobile, not aqueous
             ! May need to deal with case that pools are all zero (initial condition) which PFLOTRAN will not be able to solve.
             
            !  write(iulog,*),'Before solve'
             do poolnum=1,ndecomp_pools
               if(this%carbon_pool_mapping(poolnum)>0) &  
                 total_immobile_l2e(c,j,this%carbon_pool_mapping(poolnum)) = max(soilcarbon_l2e(c,j,poolnum)/catomw,minval)
               ! Separate N pool only exists if floating CN ratio
                !  write(iulog,*),poolnum,soilnitrogen_l2e(c,j,poolnum)
               if(decomp_cascade_con%floating_cn_ratio_decomp_pools(poolnum) .and. this%nitrogen_pool_mapping(poolnum)>0) &
                  total_immobile_l2e(c,j,this%nitrogen_pool_mapping(poolnum)) = max(soilnitrogen_l2e(c,j,poolnum)/natomw,minval/20)
             enddo
             
             CO2_before = total_immobile_l2e(c,j,this%CO2_pool_number)*catomw + &
                          total_mobile_l2e(c,j,this%CO2_pool_number)*catomw

             ! Copy dissolved nitrogen species. Units need to be converted from gN/m3 to M/L. Currently assuming saturated porosity
             
             if(this%NO3_pool_number>0) total_mobile_l2e(c,j,this%NO3_pool_number) = max(no3_l2e(c,j)/natomw,minval)
             if(this%NH4_pool_number>0) total_mobile_l2e(c,j,this%NH4_pool_number) = max(nh4_l2e(c,j)/natomw,minval)

             ! Set rate constant based on plant N demand. Convert from gN/m3/s to mol/L/s
             ! Also scale rates by relative concentrations of NO3 and NH4 so total uptake doesn't exceed demand
             ! Assumes alquimia is running in hands-off mode: Biomass term of N uptake microbial reaction is set to plant NO3 or NH4 demand
             ! This assumes the rate constant of the reaction is set to 1 in the input deck!
            if(this%plantNH4demand_pool_number>0) then
              ! Limits demand to not be too much higher than N availability to avoid cutting to tiny time step at low available N
              total_immobile_l2e(c,j,this%plantNH4demand_pool_number) = min(plantNdemand_l2e(c,j),(nh4_l2e(c,j)+no3_l2e(c,j))/dt*2)/natomw/(1000.0*porosity_l2e(c,j)*max(h2o_liqvol(c,j)/porosity_l2e(c,j),0.01))
              if(this%NO3_pool_number>0 .and. this%NH4_pool_number>0 .and. (no3_l2e(c,j)+nh4_l2e(c,j)>0)) &
                  total_immobile_l2e(c,j,this%plantNH4demand_pool_number) = total_immobile_l2e(c,j,this%plantNH4demand_pool_number)*nh4_l2e(c,j)/(nh4_l2e(c,j)+no3_l2e(c,j))
              total_immobile_l2e(c,j,this%plantNH4demand_pool_number) = max(total_immobile_l2e(c,j,this%plantNH4demand_pool_number),minval)
            endif
            if(this%plantNO3demand_pool_number>0) then
              total_immobile_l2e(c,j,this%plantNO3demand_pool_number) = min(plantNdemand_l2e(c,j),(nh4_l2e(c,j)+no3_l2e(c,j))/dt*2)/natomw/(1000.0*porosity_l2e(c,j)*max(h2o_liqvol(c,j)/porosity_l2e(c,j),0.01))
              if(this%NO3_pool_number>0 .and. this%NH4_pool_number>0 .and. (no3_l2e(c,j)+nh4_l2e(c,j)>0)) &
                total_immobile_l2e(c,j,this%plantNO3demand_pool_number) = total_immobile_l2e(c,j,this%plantNO3demand_pool_number)*no3_l2e(c,j)/(nh4_l2e(c,j)+no3_l2e(c,j))
              total_immobile_l2e(c,j,this%plantNO3demand_pool_number) = max(total_immobile_l2e(c,j,this%plantNO3demand_pool_number),minval)
            endif

           
            ! Reset diagnostic N immobilization, mineralization
             if(this%Nimm_pool_number>0) total_immobile_l2e(c,j,this%Nimm_pool_number) = minval
             if(this%Nimp_pool_number>0) total_immobile_l2e(c,j,this%Nimp_pool_number) = minval
             if(this%Nmin_pool_number>0) total_immobile_l2e(c,j,this%Nmin_pool_number) = minval

             if(this%plantNO3uptake_pool_number>0) total_immobile_l2e(c,j,this%plantNO3uptake_pool_number) = minval
             if(this%plantNO3uptake_pool_number>0) total_mobile_l2e(c,j,this%plantNO3uptake_pool_number) = minval
             if(this%plantNH4uptake_pool_number>0) total_immobile_l2e(c,j,this%plantNH4uptake_pool_number) = minval
             if(this%plantNH4uptake_pool_number>0) total_mobile_l2e(c,j,this%plantNH4uptake_pool_number) = minval

          enddo ! End of layer loop setting things up

              ! Step the chemistry solver, including advection/diffusion and timestep cutting capability for whole column
              ! Need to set surface and lateral boundary condition concentrations
          ! Surface boundary condition should be atmosphere unless there is surface water?
          ! Lateral boundary condition in MARSH mode would be saltwater if we are in the marsh column
          ! If we're in the tidal column and we want to keep track, it's concentrations in water flowing out of the marsh... Makes it trickier
          ! Need to save lateral flow for C balance
              surf_flux(:) = 0.0_r8 ! Positive means into soil
              lat_flux(:)  = 0.0_r8
              lat_bc(:) = this%bc(:) ! Currently setting to initial condition. Should update so it tracks saline/fresh
              surf_bc(:) = this%bc(:) ! Currently setting to initial condition. Should update so it tracks atmospheric O2, CO2, CH4 concentrations
              ! Assume surface water has no dissolved N. At some point should track N content of surface water though
              if(this%NO3_pool_number>0) surf_bc(this%NO3_pool_number) = 0.0_r8
              if(this%NH4_pool_number>0) surf_bc(this%NH4_pool_number) = 0.0_r8
              ! write(iulog,*),'Boundary condition',this%bc
              ! write(iulog,*),__LINE__,'adv_flow',qflx_adv_l2e(c,:)
              ! This changes total_mobile_l2e so we need to make sure we aren't using that for conservation checks

              call run_column_onestep(this, c, dt,0,max_cuts,&
                  water_density_l2e,&
                  aqueous_pressure_l2e,&
                  total_mobile_l2e,free_mobile,&
                  total_immobile_l2e,&
                  mineral_volume_fraction_l2e,&
                  mineral_specific_surface_area_l2e,&
                  surface_site_density_l2e,&
                  cation_exchange_capacity_l2e,&
                  aux_doubles_l2e,&
                  aux_ints_l2e,&
                  porosity_l2e,temperature,dz,h2o_liqvol/porosity_l2e,-qflx_adv_l2e(:,0:nlevdecomp),qflx_lat_aqu_l2e,lat_bc,lat_flux,surf_bc,surf_flux)

              ! if(max_cuts>3) write(iulog,'(a,i2,a,2i3)'),"Alquimia converged after",max_cuts," cuts. Column",c
              actual_dt_e2l(c)=dt/2**max_cuts
              ! write(iulog,*), 'lat_flux (mol/m2) = ',lat_flux
              ! write(iulog,*), 'surf_flux (mol/m2) = ',surf_flux
              ! write(iulog,*), 'bc',this%bc

              ! Save back to ELM
              water_density_e2l                 = water_density_l2e
              aqueous_pressure_e2l              = aqueous_pressure_l2e
              total_mobile_e2l                  = total_mobile_l2e
              total_immobile_e2l                = total_immobile_l2e
              mineral_volume_fraction_e2l       = mineral_volume_fraction_l2e
              mineral_specific_surface_area_e2l = mineral_specific_surface_area_l2e
              surface_site_density_e2l          = surface_site_density_l2e
              cation_exchange_capacity_e2l      = cation_exchange_capacity_l2e
              aux_doubles_e2l                   = aux_doubles_l2e
              aux_ints_e2l                      = aux_ints_l2e

              if(this%CO2_pool_number>0) then
                hr_e2l(c) = -surf_flux(this%CO2_pool_number)*catomw/dt ! Is this an issue if there is surface water?
              else
                hr_e2l(c) = 0.0_r8
              endif
              ! Surface flow of dissolved NO3 and NH4 need to be accounted for either by adding to runoff/leaching or tracking content in h2osfc
              ! Infiltration is a potential issue currently since we should really be tracking dissolved N stock in surface water as part of the column
              ! We will need to add DOC and DON runoff to ELM balance calculations eventually as well
              if(this%NO3_pool_number>0) then
                NO3runoff_e2l(c) = -surf_flux(this%NO3_pool_number)*natomw/dt - lat_flux(this%NO3_pool_number)*natomw/dt
              else
                NO3runoff_e2l(c) = 0.0_r8
              endif
              if(this%NH4_pool_number>0) then
                ! For now, including NO3 and NH4 in NO3 runoff since ELM does not include any NH4 runoff
                ! This also allows runoff to be negative if nitrogen is being carried in laterally or through infiltration
                NO3runoff_e2l(c) = NO3runoff_e2l(c) - surf_flux(this%NH4_pool_number)*natomw/dt - lat_flux(this%NH4_pool_number)*natomw/dt
              endif

              DONrunoff_e2l(c) = 0.0_r8
              DOCrunoff_e2l(c) = 0.0_r8
              DICrunoff_e2l(c) = -hr_e2l(c) ! Subtract HR to avoid double counting surface flux
              do k=1, this%chem_sizes%num_primary
                DONrunoff_e2l(c) = DONrunoff_e2l(c) - (surf_flux(k)+lat_flux(k))*this%DON_content(k)*natomw/dt
                DOCrunoff_e2l(c) = DOCrunoff_e2l(c) - (surf_flux(k)+lat_flux(k))*this%DOC_content(k)*catomw/dt
                DICrunoff_e2l(c) = DICrunoff_e2l(c) - (surf_flux(k)+lat_flux(k))*this%DIC_content(k)*catomw/dt
              enddo

          ! Loop through layers after solve and update ELM values
          do j=1,nlevdecomp

              ! Set updated land model values. Should this be moved into copy subroutine?
              ! Convert from mol/m3 to gC/m2
              do poolnum=1,ndecomp_pools
                if(this%carbon_pool_mapping(poolnum)>0) &
                  soilcarbon_e2l(c,j,poolnum) = total_immobile_e2l(c,j,this%carbon_pool_mapping(poolnum))*catomw
                ! Separate N pool only exists if floating CN ratio
                if(decomp_cascade_con%floating_cn_ratio_decomp_pools(poolnum) .and. this%nitrogen_pool_mapping(poolnum)>0) then
                   soilnitrogen_e2l(c,j,poolnum) = total_immobile_e2l(c,j,this%nitrogen_pool_mapping(poolnum))*natomw
                 elseif (this%carbon_pool_mapping(poolnum)>0) then
                   ! Calculate from CN ratio and C pool
                   soilnitrogen_e2l(c,j,poolnum) = soilcarbon_e2l(c,j,poolnum)/decomp_cascade_con%initial_cn_ratio(poolnum)
                endif
                
                ! write(iulog,*),poolnum,soilnitrogen_e2l(c,j,poolnum)
              enddo
              ! Sum together mobile and immobile pools
              ! hr_e2l goes to hr_vr (gC/m3/s)
              ! With vertical transport, comparing CO2 before/after is no longer accurate and also ignores surface exchange
              ! Best bet may be to update total HR instead of vertically resolved HR
              ! Need to add soil DIC and DOC fields to balance C
              ! if(this%CO2_pool_number>0) then 
              !   hr_e2l(c,j) = - CO2_before
              !   ! Immobile: Convert from mol/m3 to gC/m3/s
              !   hr_e2l(c,j) = hr_e2l(c,j) + total_immobile_e2l(c,j,this%CO2_pool_number)*catomw
              !   ! Mobile: convert from mol/L to gC/m3/s. mol/L*gC/mol*1000L/m3*porosity
              !   hr_e2l(c,j) = hr_e2l(c,j) + total_mobile_e2l(c,j,this%CO2_pool_number)*catomw
              !   hr_e2l(c,j) = hr_e2l(c,j)/dt
              ! endif

              DOC_e2l(c,j) = 0.0_r8
              DON_e2l(c,j) = 0.0_r8
              DIC_e2l(c,j) = 0.0_r8
              do k=1, this%chem_sizes%num_primary
                DOC_e2l(c,j) = DOC_e2l(c,j) + total_mobile_e2l(c,j,k)*catomw*this%DOC_content(k)
                DON_e2l(c,j) = DON_e2l(c,j) + total_mobile_e2l(c,j,k)*natomw*this%DON_content(k)
                DIC_e2l(c,j) = DIC_e2l(c,j) + total_mobile_e2l(c,j,k)*catomw*this%DIC_content(k)
              enddo

              carbonate_e2l(c,j) = 0.0_r8
              do k=1, this%chem_sizes%num_minerals
                ! volume fraction (m3 mineral/m3 bulk) * C_content (mol C/m3 mineral) * catomw (gC/mol) = gC/m3 bulk
                carbonate_e2l(c,j) = carbonate_e2l(c,j) + mineral_volume_fraction_e2l(c,j,k)*this%carbonate_C_content(k)*catomw
              enddo

              ! This should probably use free ion concentration or pH (in aux_output) instead of total concentration
              molperL_to_molperm3 = 1000.0*h2o_liqvol(c,j)
              if(this%Hplus_pool_number>0) then
                  pH_e2l(c,j) = -log10(free_mobile(j,this%Hplus_pool_number)/molperL_to_molperm3)
              else
                  pH_e2l(c,j) = 0.0_r8
              endif

              ! Maybe these should also use free ion concentration?
              if(this%sulfate_pool_number>0) then
                  sulfate_e2l(c,j) = total_mobile_e2l(c,j,this%sulfate_pool_number)
              else
                  sulfate_e2l(c,j) = 0.0_r8
              endif

              if(this%O2_pool_number>0) then
                  O2_e2l(c,j) = total_mobile_e2l(c,j,this%O2_pool_number)
              else
                  O2_e2l(c,j) = 0.0_r8
              endif

              if(this%chloride_pool_number>0) then
                  ! Chloride concentration needs to be converted to ppt (by mass) in water = mg/L. mol/L Cl- * 35.453 g/mol * 1.8066 g salt/g Cl * 1000 mg/g
                  salinity_e2l(c,j) = total_mobile_e2l(c,j,this%chloride_pool_number)/(1000.0*porosity_l2e(c,j)*max(h2o_liqvol(c,j)/porosity_l2e(c,j),0.01))*35.453*1.80655*1000.0
              else
                  salinity_e2l(c,j) = 0.0_r8
              endif

              if(this%Fe2_pool_number>0) then
                  Fe2_e2l(c,j) = total_mobile_e2l(c,j,this%Fe2_pool_number)
              else
                  Fe2_e2l(c,j) = 0.0_r8
              endif

              if(this%FeOH3_pool_number>0) then
                  ! Minerals need to be divided by molar volume (m3/mol) since alquimia units are m3/m3
                  ! Molar volume of Fe(OH)3 is 34.3600 cm3/mol from hanford.dat
                  FeOxide_e2l(c,j) = mineral_volume_fraction_e2l(c,j,this%FeOH3_pool_number)/34.36e-6
              else
                  FeOxide_e2l(c,j) = 0.0_r8
              endif

              if(this%NO3_pool_number>0) no3_e2l(c,j) = total_mobile_e2l(c,j,this%NO3_pool_number)*natomw
              if(this%NH4_pool_number>0) nh4_e2l(c,j) = total_mobile_e2l(c,j,this%NH4_pool_number)*natomw

              if(this%Nimm_pool_number>0) Nimm_e2l(c,j) = total_immobile_e2l(c,j,this%Nimm_pool_number)*natomw/dt
              if(this%Nimp_pool_number>0) Nimp_e2l(c,j) = total_immobile_e2l(c,j,this%Nimp_pool_number)*natomw/dt
              ! Nmin will be added to the NH4 pool elsewhere in ELM so skip that for now
              ! if(this%Nmin_pool_number>0) Nmin_e2l(c,j) = alquimia_immobile_data(this%Nmin_pool_number)*natomw/dt

              ! PFLOTRAN may use an aqueous tracer to model plant N uptake if defining using Microbial reaction
              if(this%plantNO3uptake_pool_number>0) plantNO3uptake_e2l(c,j) = (total_immobile_e2l(c,j,this%plantNO3uptake_pool_number)-minval)*natomw/dt + &
                                                (total_mobile_e2l(c,j,this%plantNO3uptake_pool_number)-minval)*natomw/dt
                if(this%plantNH4uptake_pool_number>0) plantNH4uptake_e2l(c,j) = (total_immobile_e2l(c,j,this%plantNH4uptake_pool_number)-minval)*natomw/dt + &
                                                (total_mobile_e2l(c,j,this%plantNH4uptake_pool_number)-minval)*natomw/dt



              ! Todo: Add C check
              ! Note: Generates errors if not multiplied by layer volume (imbalance on the order of 1e-8 gN/m3)
              ! Note: Generates error after restart at precision of 1e-9. But doesn't set off N conservation errors in model when precision here is relaxed.
              ! if(abs(sum(soilnitrogen_l2e(c,j,:))+no3_l2e(c,j)+nh4_l2e(c,j)-&
              !           (sum(soilnitrogen_e2l(c,j,:))+no3_e2l(c,j)+nh4_e2l(c,j)+plantNO3uptake_e2l(c,j)*dt+plantNH4uptake_e2l(c,j)*dt))*dz(c,j)>1e-5) then
              !   write(iulog,'(a,1x,i3,a,i5)'),'Nitrogen imbalance after alquimia solve step in layer',j,' Column ',c,__FILE__,__LINE__
              !   call print_alquimia_state(this,c,j)
                
              !   write(iulog,'(a25,3e20.8)'),'Total N: ', sum(soilnitrogen_l2e(c,j,:))+no3_l2e(c,j)+nh4_l2e(c,j),&
              !                               sum(soilnitrogen_e2l(c,j,:))+no3_e2l(c,j)+nh4_e2l(c,j)+plantNH4uptake_e2l(c,j)*dt+plantNO3uptake_e2l(c,j)*dt,&
              !                               sum(soilnitrogen_e2l(c,j,:))+no3_e2l(c,j)+nh4_e2l(c,j)+plantNO3uptake_e2l(c,j)*dt+plantNH4uptake_e2l(c,j)*dt-(sum(soilnitrogen_l2e(c,j,:))+no3_l2e(c,j)+nh4_l2e(c,j))
              !   write(iulog,'(a25,3e20.8)'),'SON pools: ' ,sum(soilnitrogen_l2e(c,j,:)),sum(soilnitrogen_e2l(c,j,:)),sum(soilnitrogen_e2l(c,j,:)-soilnitrogen_l2e(c,j,:))
              !   write(iulog,'(a25,3e20.8)'),'NO3: ',no3_l2e(c,j),no3_e2l(c,j),no3_e2l(c,j)-no3_l2e(c,j)
              !   write(iulog,'(a25,3e20.8)'),'NH4: ',nh4_l2e(c,j),nh4_e2l(c,j),nh4_e2l(c,j)-nh4_l2e(c,j)
              !   write(iulog,'(a25,3e20.8)'),'Plant NO3, NH4 uptake: ',plantNO3uptake_e2l(c,j)*dt,plantNH4uptake_e2l(c,j)*dt,plantNO3uptake_e2l(c,j)*dt+plantNH4uptake_e2l(c,j)*dt
              !   call endrun(msg='N imbalance after alquimia solve')
              ! endif
        enddo
     enddo
     

     ! Alquimia here calls GetAuxiliaryOutput which copies data back to interface arrays. We should do that here for EMI arrays
     ! Again, need to convert units back to ELM style, keeping track of what kind of species we are using so units are correct

#else
  implicit none
  !
  ! !ARGUMENTS
  class(em_alquimia_type)              :: this
  real(r8)             , intent(in)    :: dt ! s
  integer              , intent(in)    :: nstep
  integer              , intent(in)    :: clump_rank
  class(emi_data_list) , intent(in)    :: l2e_list
  class(emi_data_list) , intent(inout) :: e2l_list
  type(bounds_type)    , intent (in)   :: bounds_clump
  
  call endrun(msg='ERROR: Attempting to run with alquimia when model not compiled with USE_ALQUIMIA_LIB')
#endif

  end subroutine EMAlquimia_Solve_BGC
  


  
#ifdef USE_ALQUIMIA_LIB

  subroutine copy_Alquimia_to_ELM(this,c,j,water_density,&
    aqueous_pressure,&
    total_mobile,free_mobile,&
    total_immobile,&
    mineral_volume_fraction,&
    mineral_specific_surface_area,&
    surface_site_density,&
    cation_exchange_capacity,&
    aux_doubles,&
    aux_ints)

    implicit None

    ! !ARGUMENTS
    class(em_alquimia_type)              :: this
    integer                              :: c,j ! Column, layer
    ! Pointer arrays that were previously mapped using EMI
    real(r8) :: water_density(:,:), aqueous_pressure(:,:)
    real(r8) :: total_mobile(:,:,:), total_immobile(:,:,:),free_mobile(:,:)
    real(r8) :: mineral_volume_fraction(:,:,:), mineral_specific_surface_area(:,:,:)
    real(r8) :: surface_site_density(:,:,:), cation_exchange_capacity(:,:,:), aux_doubles(:,:,:)
    integer  :: aux_ints(:,:,:)

    real (c_double), pointer :: alquimia_data(:)
    integer (c_int)   , pointer :: alquimia_int_data(:)
    real(r8) :: molperL_to_molperm3

    water_density(c,j) = this%chem_state%water_density
    aqueous_pressure(c,j) = this%chem_state%aqueous_pressure

    ! We will store mobile concentrations as  mol/m3 bulk on ELM side and mol/L on alquimia side
    ! This is so changes in layer water content across time steps are properly reflected in concentrations
    molperL_to_molperm3 = 1000.0*this%chem_state%porosity*this%chem_properties%saturation
    ! write(iulog,*),'molperL_to_molperm3',molperL_to_molperm3

    ! c_f_pointer just points an array to the right data, so it needs to be actually copied
    call c_f_pointer(this%chem_state%total_mobile%data, alquimia_data, (/this%chem_sizes%num_primary/))
    ! total_mobile is converted to mol/m3 units for ELM side
    total_mobile(c,j,1:this%chem_sizes%num_primary)   = alquimia_data(1:this%chem_sizes%num_primary)*molperL_to_molperm3
    call c_f_pointer(this%chem_state%total_immobile%data, alquimia_data, (/this%chem_sizes%num_primary/))
    total_immobile(c,j,1:this%chem_sizes%num_primary)   = alquimia_data(1:this%chem_sizes%num_primary)
    call c_f_pointer(this%chem_aux_output%primary_free_ion_concentration%data, alquimia_data, (/this%chem_sizes%num_primary/))
    ! free_mobile coming out of alquimia is in molal units (mol/kg H2O). Convert to mol/m3 to match totals
    free_mobile(j,1:this%chem_sizes%num_primary)   = alquimia_data(1:this%chem_sizes%num_primary)*water_density(c,j)/1000.0_r8
    call c_f_pointer(this%chem_state%mineral_volume_fraction%data, alquimia_data, (/this%chem_sizes%num_minerals/))
    mineral_volume_fraction(c,j,1:this%chem_sizes%num_minerals)   = alquimia_data(1:this%chem_sizes%num_minerals)
    call c_f_pointer(this%chem_state%mineral_specific_surface_area%data, alquimia_data, (/this%chem_sizes%num_minerals/))
    mineral_specific_surface_area(c,j,1:this%chem_sizes%num_minerals)   = alquimia_data(1:this%chem_sizes%num_minerals)
    call c_f_pointer(this%chem_state%surface_site_density%data, alquimia_data, (/this%chem_sizes%num_surface_sites/))
    surface_site_density(c,j,1:this%chem_sizes%num_surface_sites)   = alquimia_data(1:this%chem_sizes%num_surface_sites)
    call c_f_pointer(this%chem_state%cation_exchange_capacity%data, alquimia_data, (/this%chem_sizes%num_ion_exchange_sites/))
    cation_exchange_capacity(c,j,1:this%chem_sizes%num_ion_exchange_sites)   = alquimia_data(1:this%chem_sizes%num_ion_exchange_sites)
    call c_f_pointer(this%chem_aux_data%aux_doubles%data, alquimia_data, (/this%chem_sizes%num_aux_doubles/))
    aux_doubles(c,j,1:this%chem_sizes%num_aux_doubles)   = alquimia_data(1:this%chem_sizes%num_aux_doubles)
    call c_f_pointer(this%chem_aux_data%aux_ints%data, alquimia_int_data, (/this%chem_sizes%num_aux_integers/))
    aux_ints(c,j,1:this%chem_sizes%num_aux_integers)   = alquimia_int_data(1:this%chem_sizes%num_aux_integers)

  end subroutine copy_Alquimia_to_ELM


  subroutine Copy_ELM_To_Alquimia(this,c,j,water_density,&
    aqueous_pressure,&
    total_mobile,&
    total_immobile,&
    mineral_volume_fraction,&
    mineral_specific_surface_area,&
    surface_site_density,&
    cation_exchange_capacity,&
    aux_doubles,&
    aux_ints)


    implicit None

    ! !ARGUMENTS
    class(em_alquimia_type)              :: this
    integer                              :: c,j,k ! Column, layer
    ! Pointer arrays that were previously mapped using EMI
    real(r8) :: water_density(:,:), aqueous_pressure(:,:)
    real(r8) :: total_mobile(:,:,:), total_immobile(:,:,:)
    real(r8) :: mineral_volume_fraction(:,:,:), mineral_specific_surface_area(:,:,:)
    real(r8) :: surface_site_density(:,:,:), cation_exchange_capacity(:,:,:), aux_doubles(:,:,:)
    integer  :: aux_ints(:,:,:)

    real (c_double), pointer :: alquimia_data(:)
    integer (c_int)   , pointer :: alquimia_int_data(:)

    real(r8) :: molperL_to_molperm3
    real(r8), parameter   :: minval = 1.e-35_r8

    this%chem_state%water_density = water_density(c,j)
    this%chem_state%aqueous_pressure = aqueous_pressure(c,j)

    ! We will store mobile concentrations as  mol/m3 bulk on ELM side and mol/L on alquimia side
    ! This is so changes in layer water content across time steps are properly reflected in concentrations
    molperL_to_molperm3 = 1000.0*this%chem_state%porosity*this%chem_properties%saturation

    ! c_f_pointer just points an array to the right data, so it needs to be actually copied
    call c_f_pointer(this%chem_state%total_mobile%data, alquimia_data, (/this%chem_sizes%num_primary/))
    alquimia_data(1:this%chem_sizes%num_primary) = total_mobile(c,j,1:this%chem_sizes%num_primary)/molperL_to_molperm3
    ! Don't let concentrations get below a minimal value to prevent crashes
    do k=1,this%chem_sizes%num_primary
      if(abs(alquimia_data(k))>0.0 .and. abs(alquimia_data(k))<minval) alquimia_data(k)=minval
    enddo
    call c_f_pointer(this%chem_state%total_immobile%data, alquimia_data, (/this%chem_sizes%num_primary/))
    alquimia_data(1:this%chem_sizes%num_primary) = total_immobile(c,j,1:this%chem_sizes%num_primary)
    call c_f_pointer(this%chem_state%mineral_volume_fraction%data, alquimia_data, (/this%chem_sizes%num_minerals/))
    alquimia_data(1:this%chem_sizes%num_minerals) = mineral_volume_fraction(c,j,1:this%chem_sizes%num_minerals)
    call c_f_pointer(this%chem_state%mineral_specific_surface_area%data, alquimia_data, (/this%chem_sizes%num_minerals/))
    alquimia_data(1:this%chem_sizes%num_minerals) = mineral_specific_surface_area(c,j,1:this%chem_sizes%num_minerals)
    call c_f_pointer(this%chem_state%surface_site_density%data, alquimia_data, (/this%chem_sizes%num_surface_sites/))
    alquimia_data(1:this%chem_sizes%num_surface_sites) = surface_site_density(c,j,1:this%chem_sizes%num_surface_sites)
    call c_f_pointer(this%chem_state%cation_exchange_capacity%data, alquimia_data, (/this%chem_sizes%num_ion_exchange_sites/))
    alquimia_data(1:this%chem_sizes%num_ion_exchange_sites) = cation_exchange_capacity(c,j,1:this%chem_sizes%num_ion_exchange_sites) 
    call c_f_pointer(this%chem_aux_data%aux_doubles%data, alquimia_data, (/this%chem_sizes%num_aux_doubles/))
    alquimia_data(1:this%chem_sizes%num_aux_doubles) = aux_doubles(c,j,1:this%chem_sizes%num_aux_doubles) 
    ! Messing with aux data is probably frowned upon, but very low values of free ion concentrations (<1e-200) in aux_doubles were causing crashes
    do k=1,this%chem_sizes%num_primary
      if(abs(alquimia_data(k))>0.0 .and. abs(alquimia_data(k))<minval) alquimia_data(k)=minval
    enddo
    call c_f_pointer(this%chem_aux_data%aux_ints%data, alquimia_int_data, (/this%chem_sizes%num_aux_integers/))
    alquimia_int_data(1:this%chem_sizes%num_aux_integers) = aux_ints(c,j,1:this%chem_sizes%num_aux_integers) 

  end subroutine Copy_ELM_To_Alquimia

  integer function find_alquimia_pool(pool_name,name_list,n_names) result(pool_number)
    use c_f_interface_module, only : c_f_string_ptr

    implicit none
    
    character(*),intent(in) :: pool_name
    type (c_ptr), pointer,intent(in) :: name_list(:)
    integer, intent(in) :: n_names
    
    integer :: jj
    character(len=kAlquimiaMaxStringLength) :: alq_poolname
    
    
    pool_number=-1
    
    do jj=1, n_names
      call c_f_string_ptr(name_list(jj),alq_poolname)
      if(trim(alq_poolname) == trim(pool_name)) then
        pool_number=jj
        exit
      endif
    enddo
    
  end function find_alquimia_pool

  subroutine map_alquimia_pools(this)


    use elm_varpar, only : ndecomp_pools
    use CNDecompCascadeConType, only : decomp_cascade_con
    use elm_varctl, only : alquimia_IC_name,alquimia_CO2_name,&
        alquimia_NO3_name,alquimia_NH4_name,alquimia_Nimp_name,alquimia_Nmin_name,alquimia_Nimm_name,&
        alquimia_plantNO3uptake_name,alquimia_plantNH4uptake_name,alquimia_plantNO3demand_name,alquimia_plantNH4demand_name
    use elm_varpar, only : nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions

    class(em_alquimia_type)              :: this

    integer :: ii
    character(len=kAlquimiaMaxStringLength) :: alq_poolname,donor_poolname,receiver_poolname
    type (c_ptr), pointer :: name_list(:)
    logical :: found_pool
    integer :: pool_num

    ! Map out the location of pertinent pools in Alquimia data structure
    ! Assumes that organic matter pools in PFLOTRAN are named the same as decomp_pool_name_history
    ! Currently we are not mapping any non-CTC pools.
    ! This could be a problem if chemstate_vars is initialized before the decomp cascade pool structure in ELM
    ! CN pool names in ELM are assigned in init_decompcascade_cn which is called after chemstatemod initialization and restart reading that require alquimia sizes to be set
    ! Maybe best to move this to solve step but only do it if it hasn't been done previously?
    write(iulog,*),'Alquimia carbon pool mapping:'
    allocate(this%carbon_pool_mapping(ndecomp_pools))
    call c_f_pointer(this%chem_metadata%primary_names%data, name_list, (/this%chem_sizes%num_primary/))
    do ii=1, ndecomp_pools
      if(decomp_cascade_con%floating_cn_ratio_decomp_pools(ii)) then
        alq_poolname = trim(decomp_cascade_con%decomp_pool_name_history(ii))//'C'
      else
        alq_poolname = trim(decomp_cascade_con%decomp_pool_name_history(ii))
      endif
      pool_num = find_alquimia_pool(alq_poolname,name_list,this%chem_sizes%num_primary)
      if(pool_num>0) then 
        write(iulog, '(a, i3, 1X,a7, a, i3, 1X, a)'),'ELM pool',ii,trim(decomp_cascade_con%decomp_pool_name_history(ii)),' <-> Alquimia pool',pool_num,trim(alq_poolname)
      else
        write(iulog,*),'WARNING: No match for pool',ii,trim(decomp_cascade_con%decomp_pool_name_history(ii))
      endif
      this%carbon_pool_mapping(ii)=pool_num
    enddo
    
    pool_num = find_alquimia_pool(alquimia_CO2_name,name_list,this%chem_sizes%num_primary)
    if (pool_num>0) then
      write(iulog, '(a,6x,a,i3,1x,a)'),'CO2 production', '<-> Alquimia pool',pool_num,trim(alquimia_CO2_name)
    else
      write(iulog, '(a,i3,1X,a)'),'WARNING: No match for pool',ii,trim(alquimia_CO2_name)
    endif
    this%CO2_pool_number = pool_num
    
    
    write(iulog,*),'Alquimia nitrogen pool mapping:'
    allocate(this%nitrogen_pool_mapping(ndecomp_pools))
    do ii=1, ndecomp_pools
      alq_poolname = trim(decomp_cascade_con%decomp_pool_name_history(ii))//'N'
      pool_num = find_alquimia_pool(alq_poolname,name_list,this%chem_sizes%num_primary)
      if(pool_num>0) then 
        write(iulog, '(a, i3, 1X,a7, a, i3, 1X, a)'),'ELM pool',ii,trim(decomp_cascade_con%decomp_pool_name_history(ii)),' <-> Alquimia pool',pool_num,trim(alq_poolname)
      elseif  (decomp_cascade_con%floating_cn_ratio_decomp_pools(ii)) then
        write(iulog, '(a,i3,1X,a)'),'WARNING: No match for pool',ii,trim(decomp_cascade_con%decomp_pool_name_history(ii))
      endif
      this%nitrogen_pool_mapping(ii)=pool_num
    enddo
    
    pool_num = find_alquimia_pool(alquimia_NH4_name,name_list,this%chem_sizes%num_primary)
    if (pool_num>0) then
      write(iulog, '(a,6x,a,i3,1x,a)'),'NH4', '<-> Alquimia pool',pool_num,trim(alquimia_NH4_name)
    else
      write(iulog, '(a,i3,1X,a)'),'WARNING: No match for pool',ii,trim(alquimia_NH4_name)
    endif
    this%NH4_pool_number = pool_num
    
    pool_num = find_alquimia_pool(alquimia_NO3_name,name_list,this%chem_sizes%num_primary)
    if (pool_num>0) then
      write(iulog,'(a,6x,a,i3,1x,a)'),'NO3', '<-> Alquimia pool',pool_num,trim(alquimia_NO3_name)
    else
      write(iulog,'(a,i3,1X,a)'),'WARNING: No match for pool',ii,trim(alquimia_NO3_name)
    endif
    this%NO3_pool_number = pool_num
    ! write(iulog,*),this%carbon_pool_mapping
    ! write(iulog,*),this%nitrogen_pool_mapping
    pool_num = find_alquimia_pool(alquimia_Nimm_name,name_list,this%chem_sizes%num_primary)
    if (pool_num>0) then
      write(iulog,'(a,6x,a,i3,1x,a)'),'N immobilization', '<-> Alquimia pool',pool_num,trim(alquimia_Nimm_name)
    else
      write(iulog,'(a,i3,1X,a)'),'WARNING: No match for pool',ii,trim(alquimia_Nimm_name)
    endif
    this%Nimm_pool_number = pool_num
    
    pool_num = find_alquimia_pool(alquimia_Nimp_name,name_list,this%chem_sizes%num_primary)
    if (pool_num>0) then
      write(iulog,'(a,6x,a,i3,1x,a)'),'N potential immobilization', '<-> Alquimia pool',pool_num,trim(alquimia_Nimp_name)
    else
      write(iulog,'(a,i3,1X,a)'),'WARNING: No match for pool',ii,trim(alquimia_Nimp_name)
    endif
    this%Nimp_pool_number = pool_num

    pool_num = find_alquimia_pool(alquimia_Nmin_name,name_list,this%chem_sizes%num_primary)
    if (pool_num>0) then
      write(iulog,'(a,6x,a,i3,1x,a)'),'N mineralization', '<-> Alquimia pool',pool_num,trim(alquimia_Nmin_name)
    else
      write(iulog,'(a,i3,1X,a)'),'WARNING: No match for pool',ii,trim(alquimia_Nmin_name)
    endif
    this%Nmin_pool_number = pool_num

    pool_num = find_alquimia_pool(alquimia_plantNH4uptake_name,name_list,this%chem_sizes%num_primary)
    if (pool_num>0) then
      write(iulog,'(a,6x,a,i3,1x,a)'),'Plant NH4 uptake', '<-> Alquimia pool',pool_num,trim(alquimia_plantNH4uptake_name)
    else
      write(iulog,'(a,i3,1X,a)'),'WARNING: No match for pool',ii,trim(alquimia_plantNH4uptake_name)
    endif
    this%plantNH4uptake_pool_number = pool_num

    pool_num = find_alquimia_pool(alquimia_plantNO3uptake_name,name_list,this%chem_sizes%num_primary)
    if (pool_num>0) then
      write(iulog,'(a,6x,a,i3,1x,a)'),'Plant NO3 uptake', '<-> Alquimia pool',pool_num,trim(alquimia_plantNO3uptake_name)
    else
      write(iulog,'(a,i3,1X,a)'),'WARNING: No match for pool',ii,trim(alquimia_plantNO3uptake_name)
    endif
    this%plantNO3uptake_pool_number = pool_num

    pool_num = find_alquimia_pool(alquimia_plantNH4demand_name,name_list,this%chem_sizes%num_primary)
    if (pool_num>0) then
      write(iulog,'(a,6x,a,i3,1x,a)'),'Plant NH4 demand', '<-> Alquimia pool',pool_num,trim(alquimia_plantNH4demand_name)
    else
      write(iulog,'(a,i3,1X,a)'),'WARNING: No match for pool',ii,trim(alquimia_plantNH4demand_name)
    endif
    this%plantNH4demand_pool_number = pool_num

    pool_num = find_alquimia_pool(alquimia_plantNO3demand_name,name_list,this%chem_sizes%num_primary)
    if (pool_num>0) then
      write(iulog,'(a,6x,a,i3,1x,a)'),'Plant NO3 demand', '<-> Alquimia pool',pool_num,trim(alquimia_plantNO3demand_name)
    else
      write(iulog,'(a,i3,1X,a)'),'WARNING: No match for pool',ii,trim(alquimia_plantNO3demand_name)
    endif
    this%plantNO3demand_pool_number = pool_num

    ! Need to map out reactions as well
    allocate(this%pool_reaction_mapping(ndecomp_pools))
    call c_f_pointer(this%chem_metadata%aqueous_kinetic_names%data, name_list, (/this%chem_metadata%aqueous_kinetic_names%size/))
    write(iulog,*),'Alquimia reactions:'
    do ii=1,this%chem_metadata%aqueous_kinetic_names%size
      call c_f_string_ptr(name_list(ii),alq_poolname)
      write(iulog,*),trim(alq_poolname)
    enddo
    ! cascade_receiver_pool goes to ndecomp_cascade_transitions, not ndecomp_pools
    ! But decomp_k_pools is actually by pool not by transition. So we should map based on donor pool
    do ii=1, ndecomp_cascade_transitions
      donor_poolname = decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(ii))
      if(decomp_cascade_con%cascade_receiver_pool(ii)>0) then
        receiver_poolname = decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(ii))
      else
        receiver_poolname = 'CO2'
      endif
      ! This depends on a particular PFLOTRAN/alquimia naming convention and so is not very flexible
      ! Would it be better to provide full names as inputs in varctl or something?
      alq_poolname = trim(donor_poolname)//' decay to '// trim(receiver_poolname)//' (SOMDEC sandbox)'
      pool_num = find_alquimia_pool(alq_poolname,name_list,this%chem_metadata%aqueous_kinetic_names%size)
      if(pool_num>0) then 
        write(iulog,'(a, i3, 1X,a7,a, i3, 1X, a)'),'ELM reaction',ii,trim(decomp_cascade_con%cascade_step_name(ii)),' <-> Alquimia reaction',pool_num,trim(alq_poolname)
      else
        write(iulog,'(a,i3,1x,a,1x,a)'),'WARNING: No match for reaction',ii,trim(decomp_cascade_con%cascade_step_name(ii)),':'//trim(alq_poolname)
      endif
      ! Here the index of the mapping needs to be the index of the donor pool, not the index of the transition
      this%pool_reaction_mapping(decomp_cascade_con%cascade_donor_pool(ii))=pool_num
    enddo

    ! Find plant NO3 and NH4 uptake reactions to rate constants can be set
    ! This is trickier for Microbial reactions because they are named by stoichiometry and representation depends on precision in input deck
    ! Best long-term solution is probably running in hands-off mode to avoid this entirely
    alq_poolname = '1.00000000e+00 NH4+  -> 1.00000000e+00 Tracer2' ! Todo: Fix this! !
    this%plantNH4uptake_reaction_number = find_alquimia_pool(alq_poolname,name_list,this%chem_metadata%aqueous_kinetic_names%size)
    if(this%plantNH4uptake_reaction_number>0) then 
      write(iulog,'(a, i3, 1X, a)'),'ELM plant NH4+ uptake <-> Alquimia reaction',this%plantNH4uptake_reaction_number,trim(alq_poolname)
    else
      write(iulog,'(a,1x,a)'),'WARNING: No match for plant NH4+ uptake reaction',trim(alq_poolname)
    endif
    alq_poolname = '1.00000000e+00 NO3-  -> 1.00000000e+00 Tracer' ! Todo: Fix this! !
    this%plantNO3uptake_reaction_number = find_alquimia_pool(alq_poolname,name_list,this%chem_metadata%aqueous_kinetic_names%size)
    if(this%plantNO3uptake_reaction_number>0) then 
      write(iulog,'(a, i3, 1X, a)'),'ELM plant NO3- uptake <-> Alquimia reaction',this%plantNO3uptake_reaction_number,trim(alq_poolname)
    else
      write(iulog,'(a,1x,a)'),'WARNING: No match for plant NO3- uptake reaction',trim(alq_poolname)
    endif

    ! Find aqueous gas pools
    allocate(this%is_dissolved_gas(this%chem_sizes%num_primary))
    this%is_dissolved_gas(:) = .FALSE.
    call c_f_pointer(this%chem_metadata%primary_names%data, name_list, (/this%chem_sizes%num_primary/))
    do ii=1, this%chem_sizes%num_primary
      call c_f_string_ptr(name_list(ii),alq_poolname)
      if((trim(alq_poolname) == 'CO2(aq)') .or. &
         (trim(alq_poolname) == 'HCO3-') .or. & ! This one might be tricky because of pH balance?
         (trim(alq_poolname) == 'CH4(aq)') .or. &
         (trim(alq_poolname) == 'O2(aq)')  .or. &
         (trim(alq_poolname) == 'H2S(aq)')  .or. &
         (trim(alq_poolname) == 'N2(aq)')  .or. &
         (trim(alq_poolname) == 'N2O(aq)')  .or. &
         (trim(alq_poolname) == 'H2(aq)')  ) then
        this%is_dissolved_gas(ii) = .TRUE.
        write(iulog,*),'Found alquimia dissolved gas pool: ',trim(alq_poolname)

      endif
    enddo

    ! Map DOC and DIC pools. Think about better approaches than hard coding names here
    allocate(this%DOC_content(this%chem_sizes%num_primary))
    allocate(this%DON_content(this%chem_sizes%num_primary))
    allocate(this%DIC_content(this%chem_sizes%num_primary))
    this%DOC_content(:) = 0.0_r8
    this%DIC_content(:) = 0.0_r8
    this%DON_content(:) = 0.0_r8
    call c_f_pointer(this%chem_metadata%primary_names%data, name_list, (/this%chem_sizes%num_primary/))
    do ii=1, this%chem_sizes%num_primary
      call c_f_string_ptr(name_list(ii),alq_poolname)
      if((trim(alq_poolname) == 'CO2(aq)') .or. &
         (trim(alq_poolname) == 'HCO3-') .or. &
         (trim(alq_poolname) == 'CH4(aq)') ) then
        this%DIC_content(ii) = 1.0_r8
      endif
      ! Not sure if there's a good way to pass C:N ratios from PFLOTRAN to here, but this is really clunky
      if(trim(alq_poolname) == 'DOM1') then
        this%DOC_content(ii) = 1.0_r8
        this%DON_content(ii) = 1.0_r8/100_r8*12.0110_r8/14.0067_r8
      endif
      if(trim(alq_poolname) == 'DOM2') then
        this%DOC_content(ii) = 1.0_r8
        this%DON_content(ii) = 1.0_r8/12.0_r8*12.0110_r8/14.0067_r8
      endif
      if(trim(alq_poolname) == 'DOM3') then
        this%DOC_content(ii) = 1.0_r8
        this%DON_content(ii) = 1.0_r8/16.0_r8*12.0110_r8/14.0067_r8
      endif
      if(trim(alq_poolname) == 'Acetate-') this%DOC_content(ii) = 2.0_r8
      if(this%DIC_content(ii)>0) write(iulog,'(a,1x,a,f5.2)'),'Found alquimia DIC pool: ',trim(alq_poolname),this%DIC_content(ii)
      if(this%DOC_content(ii)>0) write(iulog,'(a,1x,a,f5.2)'),'Found alquimia DOC pool: ',trim(alq_poolname),this%DOC_content(ii)
      if(this%DON_content(ii)>0) write(iulog,'(a,1x,a,f7.4)'),'Found alquimia DON pool: ',trim(alq_poolname),this%DON_content(ii)
    enddo


    ! Find other important aqueous pools to pass back to ELM
    call c_f_pointer(this%chem_metadata%primary_names%data, name_list, (/this%chem_sizes%num_primary/))

    this%Hplus_pool_number = find_alquimia_pool('H+',name_list,this%chem_sizes%num_primary)
    this%sulfate_pool_number = find_alquimia_pool('SO4--',name_list,this%chem_sizes%num_primary)
    this%O2_pool_number = find_alquimia_pool('O2(aq)',name_list,this%chem_sizes%num_primary)
    this%chloride_pool_number = find_alquimia_pool('Cl-',name_list,this%chem_sizes%num_primary)
    this%Fe2_pool_number = find_alquimia_pool('Fe++',name_list,this%chem_sizes%num_primary)

    if(this%Hplus_pool_number>0) write(iulog,'(a,6x,a,i3,1x)'),'H+', '<-> Alquimia pool',this%Hplus_pool_number
    if(this%sulfate_pool_number>0) write(iulog,'(a,6x,a,i3,1x)'),'SO4--', '<-> Alquimia pool',this%sulfate_pool_number
    if(this%O2_pool_number>0) write(iulog,'(a,6x,a,i3,1x)'),'O2(aq)', '<-> Alquimia pool',this%O2_pool_number
    if(this%chloride_pool_number>0) write(iulog,'(a,6x,a,i3,1x)'),'Cl-', '<-> Alquimia pool',this%chloride_pool_number
    if(this%Fe2_pool_number>0) write(iulog,'(a,6x,a,i3,1x)'),'Fe++', '<-> Alquimia pool',this%Fe2_pool_number
    
    ! Minerals might be trickier because they could have different stoichiometries and molar volumes
    ! Might be better to do this similar to DIC_content to allow for different Fe oxide minerals
    call c_f_pointer(this%chem_metadata%mineral_names%data, name_list, (/this%chem_sizes%num_minerals/))
    this%FeOH3_pool_number = find_alquimia_pool('Fe(OH)3',name_list,this%chem_sizes%num_minerals)
    if(this%FeOH3_pool_number>0) write(iulog,'(a,6x,a,i3,1x)'),'Fe(OH)3', '<-> Alquimia mineral',this%FeOH3_pool_number
    
    allocate(this%carbonate_C_content(this%chem_sizes%num_minerals))
    this%carbonate_C_content(:) = 0.0_r8
    do ii=1, this%chem_sizes%num_minerals
      call c_f_string_ptr(name_list(ii),alq_poolname)
      ! One C per mol of calcite divided by molar volume (m3/mol) to give units of mol C/m3
      if(trim(alq_poolname) == 'Calcite') this%carbonate_C_content(ii)=1.0_r8/36.9340e-6_r8
      if(this%carbonate_C_content(ii)>0) write(iulog,'(a,1x,a,x,f10.4,x,a)'),'Found alquimia carbonate mineral pool: ',trim(alq_poolname),this%carbonate_C_content(ii),'mol C/m^3'
    enddo


  end subroutine map_alquimia_pools

  
  subroutine print_alquimia_state(this,c,j)

    use elm_varpar, only : ndecomp_pools,ndecomp_cascade_transitions
    use iso_c_binding, only : c_f_pointer, c_double
    use c_f_interface_module, only : c_f_string_ptr
    use CNDecompCascadeConType, only : decomp_cascade_con

    implicit none

    class(em_alquimia_type)              :: this
    integer, intent(in) :: c,j

    integer :: poolnum
    character(len=256) :: poolname
    real (c_double), pointer :: alquimia_mobile_data(:), alquimia_free_data(:), alquimia_immobile_data(:), &
                                alquimia_mineral_data(:),alquimia_mineral_SSA(:),alquimia_aux_data(:)
    type (c_ptr), pointer :: name_list(:)

    call c_f_pointer(this%chem_state%total_immobile%data, alquimia_immobile_data, (/this%chem_sizes%num_primary/))
    call c_f_pointer(this%chem_state%total_mobile%data, alquimia_mobile_data, (/this%chem_sizes%num_primary/))
    call c_f_pointer(this%chem_aux_output%primary_free_ion_concentration%data, alquimia_free_data, (/this%chem_sizes%num_primary/))
    call c_f_pointer(this%chem_state%mineral_volume_fraction%data, alquimia_mineral_data, (/this%chem_sizes%num_minerals/))
    call c_f_pointer(this%chem_state%mineral_specific_surface_area%data, alquimia_mineral_SSA, (/this%chem_sizes%num_minerals/))
    ! call c_f_pointer(this%chem_properties%aqueous_kinetic_rate_cnst%data, alquimia_rates_data, (/this%chem_properties%aqueous_kinetic_rate_cnst%size/))
    
    call c_f_pointer(this%chem_aux_data%aux_doubles%data, alquimia_aux_data, (/this%chem_sizes%num_aux_doubles/))
    write(iulog,'(a)'),'Alquimia aux doubles:'
    do poolnum=1,this%chem_sizes%num_aux_doubles
      write(iulog,'(a,i3,e18.5)'),'Aux double ',poolnum,alquimia_aux_data(poolnum)
    enddo

    call c_f_pointer(this%chem_metadata%primary_names%data, name_list, (/this%chem_sizes%num_primary/))
    write(iulog,'(a)'),'Alquimia primary species (mol/m3 bulk):'
    write(iulog,'(23x,a22,5x,a22,5x,a22)'),'Immobile (mol/m3 bulk)','Tot Mobile (mol/L H2O)','Free (mol/L H2O)'
    do poolnum=1,this%chem_sizes%num_primary
      call c_f_string_ptr(name_list(poolnum),poolname)
      write(iulog,'(i3,a20,e22.5,5x,e22.5,5x,e22.5)'),poolnum,trim(poolname),alquimia_immobile_data(poolnum),alquimia_mobile_data(poolnum),alquimia_free_data(poolnum)
    enddo

    call c_f_pointer(this%chem_metadata%mineral_names%data, name_list, (/this%chem_sizes%num_minerals/))
    write(iulog,'(a)'),'Alquimia minerals:'
    write(iulog,'(23x,a22,5x,a22)'),'Vol. frac. (m3/m3)', 'SSA (m2/m3)'
    do poolnum=1,this%chem_sizes%num_minerals
      call c_f_string_ptr(name_list(poolnum),poolname)
      write(iulog,'(i3,a20,e22.5,5x,e22.5)'),poolnum,trim(poolname),alquimia_mineral_data(poolnum),alquimia_mineral_SSA(poolnum)
    enddo

    write(iulog,'(a,f10.3)'),'Porosity      = ',this%chem_state%porosity
    write(iulog,'(a,f10.3)'),'Saturation    = ',this%chem_properties%saturation
    write(iulog,'(a,f10.3)'),'Temperature   = ',this%chem_state%temperature
    write(iulog,'(a,f10.3)'),'Pressure      = ',this%chem_state%aqueous_pressure
    write(iulog,'(a,f10.3)'),'Water density = ',this%chem_state%water_density
    write(iulog,'(a,f10.3)'),'Volume        = ',this%chem_properties%volume


  end subroutine print_alquimia_state

  recursive subroutine run_onestep(this,c,j,dt,num_cuts,max_cuts)
    
    use c_f_interface_module, only : c_f_string_ptr
    
    implicit none
    
    class(em_alquimia_type)              :: this
    integer,intent(out)                  :: max_cuts
    integer,intent(in)                   :: num_cuts,c,j
    real(r8),intent(in)                  :: dt
    
    real(r8) :: actual_dt,porosity
    character(512) :: msg
    character(kind=C_CHAR,len=kAlquimiaMaxStringLength) :: status_message
    integer :: ncuts2,ncuts,ii
    
    max_cuts = num_cuts
    actual_dt = dt/(2**num_cuts)
    
    ncuts=0
    ncuts2=0
    
    porosity=this%chem_state%porosity
    call this%chem%ReactionStepOperatorSplit(this%chem_engine, actual_dt, this%chem_properties, this%chem_state, &
                                           this%chem_aux_data, this%chem_status)
    ! Reset porosity because Pflotran tends to mess it up
    this%chem_state%porosity=porosity
    ! write(iulog,*),'Converged =',this%chem_status%converged,"ncuts =",num_cuts
    if (this%chem_status%converged) then
      ! Success. Can get aux output and finish execution of the subroutine
      ! Get auxiliary output
      call this%chem%getAuxiliaryOutput(this%chem_engine, this%chem_properties, this%chem_state, &
                                  this%chem_aux_data, this%chem_aux_output, this%chem_status)
      if(this%chem_status%error /= 0) then
        call c_f_string_ptr(this%chem_status%message,status_message)
        call endrun(msg='Alquimia error in ReactionStepOperatorSplit: '//status_message)
      endif
      
    else ! Solve did not converge. Cut timestep, and bail out if too short
      if(actual_dt/2 < min_dt) then
        call c_f_string_ptr(this%chem_status%message,status_message)
        write(msg,'(a,i3,a,f5.2,a,i4,a,i3,a,i5)') "Error: Alquimia ReactionStepOperatorSplit failed to converge after ",num_cuts," cuts to dt = ",actual_dt,' s. Newton iterations = ',this%chem_status%num_newton_iterations,' Layer = ',j," Col = ",c
        call print_alquimia_state(this,c,j)
        call endrun(msg=msg)
      else
        ! If we are not at minimum timestep yet, cut and keep going
        ! Need to run the step two times because we have cut the timestep in half
        call run_onestep(this, c,j, dt,num_cuts+1,ncuts)
        if(ncuts>max_cuts) max_cuts=ncuts
        ! write(iulog,*),'Converged =',this%chem_status%converged,"ncuts =",ncuts,'(Substep 1)'
        
        ! The second one starts from the maximum number of cuts from the first one so it doesn't waste time retrying a bunch of failed timestep lengths
         do ii=1,2**(max_cuts-(num_cuts+1))
           call run_onestep(this, c,j, dt,ncuts,ncuts2)
           if(ncuts2>max_cuts) max_cuts=ncuts2
        !   write(iulog,*),'Converged =',this%chem_status%converged,"ncuts =",ncuts2,'. Substep 2 +',ii
         enddo
        ! call run_onestep(this, c,j, dt,num_cuts+1,ncuts)
        ! if(ncuts>max_cuts) max_cuts=ncuts
        ! write(iulog,*),'Converged =',this%chem_status%converged,"ncuts =",ncuts,'(Substep 2)'
      endif
    endif
      
      

  end subroutine run_onestep

  recursive subroutine run_column_onestep(this,c,dt,num_cuts,max_cuts, &
          water_density,&
          aqueous_pressure,&
          total_mobile,free_mobile,&
          total_immobile,&
          mineral_volume_fraction,&
          mineral_specific_surface_area,&
          surface_site_density,&
          cation_exchange_capacity,&
          aux_doubles,&
          aux_ints,&
          porosity,temperature,volume,saturation,adv_flux,lat_flow,lat_bc,lat_flux,surf_bc,surf_flux)
    
  use c_f_interface_module, only : c_f_string_ptr
  use elm_varpar       , only : nlevdecomp
  use elm_varcon, only : dzsoi_decomp
  use shr_infnan_mod         , only : isnan => shr_infnan_isnan
  
  implicit none
  
  class(em_alquimia_type)              :: this
  integer,intent(out)                  :: max_cuts
  integer,intent(in)                   :: num_cuts,c
  real(r8),intent(in)                  :: dt
  real(r8),intent(inout),pointer       :: water_density(:,:),&
                                          aqueous_pressure(:,:),&
                                          total_mobile(:,:,:),&
                                          total_immobile(:,:,:),&
                                          mineral_volume_fraction(:,:,:),&
                                          mineral_specific_surface_area(:,:,:),&
                                          surface_site_density(:,:,:),&
                                          cation_exchange_capacity(:,:,:),&
                                          aux_doubles(:,:,:)
  integer,intent(inout)   ,pointer     :: aux_ints(:,:,:)
  real(r8),intent(in),dimension(:,:)  :: porosity,temperature,volume,saturation,lat_flow
  real(r8),intent(in),dimension(:,:)   :: adv_flux
  real(r8),intent(in),dimension(:)   :: lat_bc, surf_bc
  real(r8),intent(inout)             :: surf_flux(:), lat_flux(:),free_mobile(:,:) ! Total (cumulative) surface flux in time step. Units of mol/time step

    real(r8)             :: water_density_tmp(1,nlevdecomp),&
                            aqueous_pressure_tmp(1,nlevdecomp),&
                            total_mobile_tmp(1,nlevdecomp,this%chem_sizes%num_primary),&
                            free_mobile_tmp(nlevdecomp,this%chem_sizes%num_primary),&
                            total_immobile_tmp(1,nlevdecomp,this%chem_sizes%num_primary),&
                            mineral_volume_fraction_tmp(1,nlevdecomp,this%chem_sizes%num_minerals),&
                            mineral_specific_surface_area_tmp(1,nlevdecomp,this%chem_sizes%num_minerals),&
                            surface_site_density_tmp(1,nlevdecomp,this%chem_sizes%num_surface_sites),&
                            cation_exchange_capacity_tmp(1,nlevdecomp,this%chem_sizes%num_ion_exchange_sites),&
                            aux_doubles_tmp(1,nlevdecomp,this%chem_sizes%num_aux_doubles)
    integer            ::   aux_ints_tmp(1,nlevdecomp,this%chem_sizes%num_aux_integers)
    real(r8) :: diffus(nlevdecomp), sat(nlevdecomp)
    real(r8) :: transport_change_rate(nlevdecomp,this%chem_sizes%num_primary),source_term(nlevdecomp,this%chem_sizes%num_primary)
    real(r8) :: surf_adv_step(this%chem_sizes%num_primary),surf_equil_step(this%chem_sizes%num_primary), lat_flux_step(this%chem_sizes%num_primary)
    ! real(r8) :: bot_adv_step(this%chem_sizes%num_primary)
  
  real(r8) :: actual_dt,porosity_tmp
  character(512) :: msg
  character(kind=C_CHAR,len=kAlquimiaMaxStringLength) :: status_message
  integer :: ncuts2,ncuts,ii,j,k
  
  max_cuts = num_cuts
  actual_dt = dt/(2**num_cuts)
  
  ncuts=0
  ncuts2=0

  do j=1,nlevdecomp
    sat(j) = min(max(saturation(c,j),0.01),1.0)
  enddo

  do k=1,this%chem_sizes%num_primary
    diffus(:) = 0.0_r8
    surf_equil_step(k) = 0.0_r8
    lat_flux_step(k) = 0.0_r8
    surf_adv_step(k) = 0.0_r8
    ! Set diffusion coefficient depending on saturation and whether species is aqueous gas or not
    ! Need to set boundary condition concentrations for adv flux (top layer infiltration) and lateral flux (source)

    ! Skip species that are not actually mobile
    if(k == this%plantNH4uptake_pool_number .or. k == this%plantNO3uptake_pool_number) cycle

    if(this%is_dissolved_gas(k)) then
      ! For gases, diffusion rates are set using gas diffusive transport (Meslin et al., SSSAJ, 2010. doi:10.2136/sssaj2009.0474)
      ! Estimating gas diffusion coefficient of 0.2 cm2/s and dry soil diffusion coefficient of 30% of gas (Moldrup et al 2004, SSSAJ)
      do j=1,nlevdecomp
        diffus(j) = 2.0e-5_r8*0.3_r8*(1.0_r8 - sat(j))**2.5
      enddo

      ! Equilibrate top layer of dissolved gases w.r.t. upper BC. BC is in mol/m3 H2O units and total_mobile is in mol/m3 units
      ! Unless this should be treated as a source term in advection-diffusion?
      ! Possible issue in low-moisture conditions: Total layer stock of O2 might be really low because not that much fits in the small amount of water
      ! In reality water should be in equilibrium with soil pore air space
      ! If we don't multiply by sat here, I guess we shove all the layer gas into the water...
      surf_equil_step(k) = ( surf_bc(k)*porosity(c,1)*max(sat(1),0.3) - total_mobile(c,1,k) )*dzsoi_decomp(1)
      ! write(iulog,*),'Dissolved gas',k,'BC',surf_bc(k)*porosity(c,1)*sat(1),'Surf conc',total_mobile(c,1,k),'(mol m-3 equivalent)','porosity',porosity(c,1),'saturation',sat(1),'flux',surf_equil_step(k)
      total_mobile(c,1,k) = surf_bc(k)*porosity(c,1)*max(sat(1),0.3)
    endif
    
    do j=1,nlevdecomp
      if(isnan(total_mobile(c,j,k))) then
        write(iulog,*),__LINE__,'Chem spec',k,total_mobile(c,:,k)
        call endrun(msg="Mobile species is NaN")
      endif
      ! Assume diffusion through water according to Wright (1990)
      ! In that paper diffus_water = 0.000025 cm2/s
      diffus(j) = diffus(j) + 2.5e-9_r8*0.005_r8*exp(10.0_r8*sat(j)*porosity(c,j))

      ! Source term is lateral flow. For inflow, use lateral boundary condition. For outflow, use local concentration
      ! lat_flux units are mm H2O/s = 1e-3 m3 h2o/m2/s
      ! lat_bc in units of mol/m3 H2O
      ! source_term in mol/m3 bulk/s
      if(lat_flow(c,j) > 0) then
        source_term(j,k) = lat_flow(c,j)*1e-3_r8 * lat_bc(k)*porosity(c,j)*sat(j) ! mol/m3 bulk/s
      else
        source_term(j,k) = lat_flow(c,j)*1e-3_r8 * total_mobile(c,j,k)
      endif
      lat_flux_step(k) = lat_flux_step(k) + source_term(j,k)*dzsoi_decomp(j)

    enddo
    
      ! adv_flux units are mm H2O/s
    if(adv_flux(c,1)<0.0_r8) then ! Downward flow uses surface boundary condition
      surf_adv_step(k) = - adv_flux(c,1)*1e-3_r8*surf_bc(k)*actual_dt/2 
    else ! Upward flow uses surface layer concentration. Should this concentration be per bulk volume or per water volume?
      surf_adv_step(k) = - adv_flux(c,1)*1e-3_r8*total_mobile(c,1,k)*actual_dt/2
    endif
    ! if(adv_flux(c,nlevdecomp+1)<0.0_r8) then
      ! bot_adv_step(k) = -adv_flux(c,nlevdecomp+1)*1e-3_r8*total_mobile(c,nlevdecomp,k)*actual_dt/2
      ! write(iulog,*) 'Flow at bottom',adv_flux(c,nlevdecomp+1),-adv_flux(c,nlevdecomp+1)*1e-3_r8*total_mobile(c,nlevdecomp,k)*actual_dt/2
    ! endif

    ! At this point, total_mobile is stored as mol/m3 bulk (ELM side). Dividing by porosity*saturation converts to mol/m3 water
    ! Note adv_flux is defined in advection_diffusion as <0 being downward
    ! write(iulog,*) 'Before adv_diff. ncuts = ',num_cuts
    ! write(iulog,*) 'diffus',diffus
    ! write(iulog,*) 'adv_flux',adv_flux(c,:)
    ! write(iulog,*) 'source',source_term(:,k)/(porosity(c,:)*sat(:))
    ! write(iulog,*) 'total_mobile',total_mobile(c,:,k)/(porosity(c,:)*sat(:))
    ! write(iulog,*) 'lat_flow',lat_flow(c,:)
    ! write(iulog,*) 'porosity',porosity(c,:)
    ! write(iulog,*) 'saturation',sat(:) ! Need to account for when saturation is 0
    ! write(iulog,*),'Mobile spec',k,'Before: ',total_mobile(c,1:nlevdecomp,k)
    ! write(iulog,*),__LINE__,'adv_flux',adv_flux(c,1:nlevdecomp+1)
    call advection_diffusion(total_mobile(c,1:nlevdecomp,k),adv_flux(c,1:nlevdecomp+1)*1e-3,diffus(1:nlevdecomp),& 
      source_term(1:nlevdecomp,k),&
      surf_bc(k),actual_dt/2,transport_change_rate(1:nlevdecomp,k))
    ! At this point perhaps we should go through and re-equilibrate dissolved gases in top layer if unsaturated?
    ! write(iulog,*) 'change rate',transport_change_rate(:,k)

    total_mobile(c,1:nlevdecomp,k) = total_mobile(c,1:nlevdecomp,k) + transport_change_rate(1:nlevdecomp,k)*actual_dt/2
    ! write(iulog,*),'Mobile spec',k,'After: ',total_mobile(c,1:nlevdecomp,k)
    ! write(iulog,*),'Diff rate',transport_change_rate(1:nlevdecomp,k)*dzsoi_decomp(1:nlevdecomp)
    ! write(iulog,*),k,'Total diff',sum(transport_change_rate(1:nlevdecomp,k)*dzsoi_decomp(1:nlevdecomp))*actual_dt/2,'Surf adv',surf_adv_step(k),'Surf equil',surf_equil_step(k)
  enddo



  do j=1,nlevdecomp

    ! Update properties from ELM
    this%chem_state%porosity =    porosity(c,j)
    this%chem_state%temperature = temperature(c,j) - 273.15
    this%chem_properties%volume = volume(c,j)
    this%chem_properties%saturation = sat(j) ! Set minimum saturation to stop concentrations from blowing up at low soil moisture
    
    call this%copy_ELM_to_Alquimia(c,j,water_density,&
          aqueous_pressure,&
          total_mobile,&
          total_immobile,&
          mineral_volume_fraction,&
          mineral_specific_surface_area,&
          surface_site_density,&
          cation_exchange_capacity,&
          aux_doubles,&
          aux_ints) 

    porosity_tmp=this%chem_state%porosity
    call this%chem%ReactionStepOperatorSplit(this%chem_engine, actual_dt, this%chem_properties, this%chem_state, &
                                         this%chem_aux_data, this%chem_status)
    ! Reset porosity because Pflotran tends to mess it up
    this%chem_state%porosity=porosity_tmp

    ! Get auxiliary output
    call this%chem%getAuxiliaryOutput(this%chem_engine, this%chem_properties, this%chem_state, &
    this%chem_aux_data, this%chem_aux_output, this%chem_status)

    ! Check for error
    if(this%chem_status%error /= 0) then
      call c_f_string_ptr(this%chem_status%message,status_message)
      call endrun(msg='Alquimia error in getAuxiliaryOutput: '//status_message)
    endif

    ! In top layer, cut time step if gas species being absorbed very fast because it should be close to equilibrium
    if(j==1) then
      do k=1,this%chem_sizes%num_primary
        if(this%is_dissolved_gas(k) .and. (surf_bc(k) > 0.0) .and. &
            ((surf_bc(k)*porosity(c,1)*max(sat(1),0.3) - total_mobile(c,1,k) )/(surf_bc(k)*porosity(c,1)*max(sat(1),0.3)) > 0.25)) then
              this%chem_status%converged = .FALSE.
              if(num_cuts>6) write(iulog,'(a,f5.2,x,a,i4,a)'),'Cutting time step to dt = ',actual_dt,' because species',k,'reduced too fast in layer 1'
        endif
      enddo
    endif


    if (this%chem_status%converged) then
      ! Success. Can finish execution of the subroutine

      ! Copy back to column structure
      call this%copy_Alquimia_to_ELM(1,j,water_density_tmp,&
        aqueous_pressure_tmp,&
        total_mobile_tmp,free_mobile_tmp,&
        total_immobile_tmp,&
        mineral_volume_fraction_tmp,&
        mineral_specific_surface_area_tmp,&
        surface_site_density_tmp,&
        cation_exchange_capacity_tmp,&
        aux_doubles_tmp,&
        aux_ints_tmp) 
      
    else ! Solve did not converge. Cut timestep, and bail out if too short
      if(actual_dt/2 < min_dt) then
        call c_f_string_ptr(this%chem_status%message,status_message)
        ! Sometimes solve fails because time step is too short (not sure why)
        ! I wonder if at bailout we should first try pausing transport and solving each layer separately?
        write(msg,'(a,i3,a,f5.3,a,i4,a,i3,a,i5)') "Error: Alquimia ReactionStepOperatorSplit failed to converge after ",num_cuts," cuts to dt = ",actual_dt,' s. Newton iterations = ',this%chem_status%num_newton_iterations,' Layer = ',j," Col = ",c
        call print_alquimia_state(this,c,j)
        call endrun(msg=msg)
      else
        exit ! Drop out of the layer loop to start over at shorter time step
      endif
    endif
    enddo ! Layer loop

    if(.not. this%chem_status%converged) then
        ! If we are not at minimum timestep yet, cut and keep going

      if(actual_dt<=-30.0_r8) then
        write(iulog,*),'Alquimia: Time step cut to 30 s. Attempting to solve by pausing transport and solving layer by layer'
        do j=1,nlevdecomp
              ! Update properties from ELM
          this%chem_state%porosity =    porosity(c,j)
          this%chem_state%temperature = temperature(c,j) - 273.15
          this%chem_properties%volume = volume(c,j)
          this%chem_properties%saturation = sat(j) ! Set minimum saturation to stop concentrations from blowing up at low soil moisture
          
          call this%copy_ELM_to_Alquimia(c,j,water_density,&
                                            aqueous_pressure,&
                                            total_mobile,&
                                            total_immobile,&
                                            mineral_volume_fraction,&
                                            mineral_specific_surface_area,&
                                            surface_site_density,&
                                            cation_exchange_capacity,&
                                            aux_doubles,&
                                            aux_ints) 
          call run_onestep(this,c,j,dt,num_cuts,ncuts)
          call this%copy_Alquimia_to_ELM(1,j,water_density_tmp,&
                                          aqueous_pressure_tmp,&
                                          total_mobile_tmp,free_mobile_tmp,&
                                          total_immobile_tmp,&
                                          mineral_volume_fraction_tmp,&
                                          mineral_specific_surface_area_tmp,&
                                          surface_site_density_tmp,&
                                          cation_exchange_capacity_tmp,&
                                          aux_doubles_tmp,&
                                          aux_ints_tmp)
        enddo
      else
    
        ! Here we are basically throwing out all the _tmp array values and starting over with the originals

        ! Also need to undo transport because we are starting this time step over
      ! Unless we change transport to act on temp arrays
       do k=1,this%chem_sizes%num_primary
        if(k == this%plantNH4uptake_pool_number .or. k == this%plantNO3uptake_pool_number) cycle
        total_mobile(c,1:nlevdecomp,k) = total_mobile(c,1:nlevdecomp,k) &
                        - transport_change_rate(1:nlevdecomp,k)*actual_dt/2
        total_mobile(c,1,k) = total_mobile(c,1,k) - surf_equil_step(k)/dzsoi_decomp(1)
      enddo
        ! write(iulog,'(a,f8.1,a,i3,a,i3)'),'Cutting time step to',actual_dt/2,' s, layer',j,', column',c
        ! Need to run the step two times because we have cut the timestep in half
        call run_column_onestep(this, c, dt,num_cuts+1,ncuts,&
          water_density,&
          aqueous_pressure,&
          total_mobile,free_mobile,&
          total_immobile,&
          mineral_volume_fraction,&
          mineral_specific_surface_area,&
          surface_site_density,&
          cation_exchange_capacity,&
          aux_doubles,&
          aux_ints,porosity,temperature,volume,saturation,adv_flux,lat_flow,lat_bc,lat_flux,surf_bc,surf_flux)

        if(ncuts>max_cuts) max_cuts=ncuts
        ! write(iulog,*),'Converged =',this%chem_status%converged,"ncuts =",ncuts,'(Substep 1)'
        
        ! The second one starts from the maximum number of cuts from the first one so it doesn't waste time retrying a bunch of failed timestep lengths
        do ii=1,2**(max_cuts-(num_cuts+1))
          call run_column_onestep(this, c, dt,ncuts,ncuts2,&
          water_density,&
          aqueous_pressure,&
          total_mobile,free_mobile,&
          total_immobile,&
          mineral_volume_fraction,&
          mineral_specific_surface_area,&
          surface_site_density,&
          cation_exchange_capacity,&
          aux_doubles,&
          aux_ints,porosity,temperature,volume,saturation,adv_flux,lat_flow,lat_bc,lat_flux,surf_bc,surf_flux)
          if(ncuts2>max_cuts) max_cuts=ncuts2
        !   write(iulog,*),'Converged =',this%chem_status%converged,"ncuts =",ncuts2,'. Substep 2 +',ii
        enddo
      endif ! Attempt at running layer by layer
        ! call run_onestep(this, c,j, dt,num_cuts+1,ncuts)
        ! if(ncuts>max_cuts) max_cuts=ncuts
        ! write(iulog,*),'Converged =',this%chem_status%converged,"ncuts =",ncuts,'(Substep 2)'
      else ! It did converge and we made it to the end of the layer loop

    ! At this point we've successfully updated the column chemistry for all layers. Copy back to inout arrays
    water_density(c,:) = water_density_tmp(1,:)
    aqueous_pressure(c,:) = aqueous_pressure_tmp(1,:)
    total_mobile(c,:,:) = total_mobile_tmp(1,:,:)
    free_mobile(:,:) = free_mobile_tmp(:,:)
    total_immobile(c,:,:) = total_immobile_tmp(1,:,:)
    mineral_volume_fraction(c,:,:) = mineral_volume_fraction_tmp(1,:,:)
    mineral_specific_surface_area(c,:,:) = mineral_specific_surface_area_tmp(1,:,:)
    surface_site_density(c,:,:) = surface_site_density_tmp(1,:,:)
    cation_exchange_capacity(c,:,:) = cation_exchange_capacity_tmp(1,:,:)
    aux_doubles(c,:,:) = aux_doubles_tmp(1,:,:)
    aux_ints(c,:,:) = aux_ints_tmp(1,:,:)

    surf_flux = surf_flux + surf_adv_step + surf_equil_step
    lat_flux  = lat_flux  + lat_flux_step
      
    ! Second half of transport (Strang splitting)
    ! This is only done if we converged at this time step for all layers

    do k=1,this%chem_sizes%num_primary
      ! Set diffusion coefficient depending on saturation and whether species is aqueous gas or not
      ! Need to set boundary condition concentrations for adv flux (top layer infiltration) and lateral flux (source)
      diffus(:) = 0.0_r8
      surf_equil_step(k) = 0.0_r8
      lat_flux_step(k) = 0.0_r8
      surf_adv_step(k) = 0.0_r8
  
      ! Skip species that are not actually mobile
      if(k == this%plantNH4uptake_pool_number .or. k == this%plantNO3uptake_pool_number) cycle

      
      if(this%is_dissolved_gas(k)) then
        ! For gases, diffusion rates are set using gas diffusive transport (Meslin et al., SSSAJ, 2010. doi:10.2136/sssaj2009.0474)
        ! Estimating gas diffusion coefficient of 0.2 cm2/s and dry soil diffusion coefficient of 30% of gas (Moldrup et al 2004, SSSAJ)
        ! Consider some different diffusion exponents? Previous CORPSE work ended up with an exponent of 0.6 rather than 2.5. But Cusack project gave vals around 1.5-2.5
        do j=1,nlevdecomp
          diffus(j) = 2.0e-5_r8*0.3_r8*(1.0_r8 - sat(j))**2.5
        enddo
  
        ! Equilibrate top layer of dissolved gases w.r.t. upper BC. BC is in mol/m3 H2O units and total_mobile is in mol/m3 units
        ! write(iulog,*),'Dissolved gas',k,'BC',surf_bc(k)*porosity(c,1)*sat(1),'Surf conc',total_mobile(c,1,k),'(mol m-3 equivalent)','porosity',porosity(c,1),'saturation',sat(1)
        surf_equil_step(k) = ( surf_bc(k)*porosity(c,1)*max(sat(1),0.3) - total_mobile(c,1,k) )*dzsoi_decomp(1)
        total_mobile(c,1,k) = surf_bc(k)*porosity(c,1)*max(sat(1),0.3)
      
      endif
      
      do j=1,nlevdecomp
        if(isnan(total_mobile(c,j,k))) then
          write(iulog,*),__LINE__,'Chem ',k,total_mobile(c,:,k)
          call endrun(msg="Mobile species is NaN")
        endif
        ! Assume diffusion through water according to Wright (1990)
        ! In that paper diffus_water = 0.000025 cm2/s
        diffus(j) = diffus(j) + 2.5e-9_r8*0.005_r8*exp(10.0_r8*sat(j)*porosity(c,j))
  
        ! Source term is lateral flow. For inflow, use lateral boundary condition. For outflow, use local concentration
        ! lat_flux units are mm H2O/s = 1e-3 m3 h2o/m2/s
        ! lat_bc in units of mol/m3 H2O
        ! source_term in mol/m3 bulk/s
        if(lat_flow(c,j) > 0) then
          source_term(j,k) = lat_flow(c,j)*1e-3_r8 * lat_bc(k)*porosity(c,j)*sat(j) ! mol/m3 bulk/s
        else
          source_term(j,k) = lat_flow(c,j)*1e-3_r8 * total_mobile(c,j,k)
        endif
        lat_flux_step(k) = lat_flux_step(k) + source_term(j,k)*dzsoi_decomp(j)
  
      enddo

      ! adv_flux units are mm H2O/s
      if(adv_flux(c,1)<0.0_r8) then ! Downward flow uses surface boundary condition
        surf_adv_step(k) = - adv_flux(c,1)*1e-3_r8*surf_bc(k)*actual_dt/2 
      else ! Upward flow uses surface layer concentration. Should this concentration be per bulk volume or per water volume?
        surf_adv_step(k) = - adv_flux(c,1)*1e-3_r8*total_mobile(c,1,k)*actual_dt/2
      endif    
      ! if(adv_flux(c,nlevdecomp+1)>0.0_r8) then
        ! bot_adv_step(k) = -adv_flux(c,nlevdecomp+1)*1e-3_r8*total_mobile(c,nlevdecomp,k)*actual_dt/2
        ! write(iulog,*) 'Flow at bottom',adv_flux(c,nlevdecomp+1),-adv_flux(c,nlevdecomp+1)*1e-3_r8*total_mobile(c,nlevdecomp,k)*actual_dt/2
      ! endif

      ! At this point, total_mobile is stored as mol/m3 bulk (ELM side). Dividing by porosity*saturation converts to mol/m3 water
      ! Note adv_flux is defined in advection_diffusion as <0 being downward
      ! write(iulog,*) 'Before adv_diff. ncuts = ',num_cuts
      ! write(iulog,*) 'diffus',diffus
      ! write(iulog,*) 'adv_flux',adv_flux(c,:)
      ! write(iulog,*) 'source',source_term(:,k)/(porosity(c,:)*sat(:))
      ! write(iulog,*) 'total_mobile',total_mobile(c,:,k)/(porosity(c,:)*sat(:))
      ! write(iulog,*) 'lat_flow',lat_flow(c,:)
      ! write(iulog,*) 'porosity',porosity(c,:)
      ! write(iulog,*) 'saturation',sat(:) ! Need to account for when saturation is 0
      ! write(iulog,*),'Mobile spec',k,'Before: ',total_mobile(c,1:nlevdecomp,k)
      ! write(iulog,*),__LINE__,'adv_flux',adv_flux(c,1:nlevdecomp+1)
      call advection_diffusion(total_mobile(c,1:nlevdecomp,k),adv_flux(c,1:nlevdecomp+1)*1e-3,diffus(1:nlevdecomp),& 
        source_term(1:nlevdecomp,k),&
        surf_bc(k),actual_dt/2,transport_change_rate(1:nlevdecomp,k))
      ! At this point perhaps we should go through and re-equilibrate dissolved gases in top layer if unsaturated?
      ! write(iulog,*) 'change rate',transport_change_rate(:,k)
  
    ! Here need to convert back from mol/m3 water to mol/m3 bulk
      total_mobile(c,1:nlevdecomp,k) = total_mobile(c,1:nlevdecomp,k) + transport_change_rate(1:nlevdecomp,k)*actual_dt/2
      ! write(iulog,*),'Mobile spec',k,'After: ',total_mobile(c,1:nlevdecomp,k)
      ! write(iulog,*),'Diff rate',transport_change_rate(1:nlevdecomp,k)*dzsoi_decomp(1:nlevdecomp)
      ! write(iulog,*),k,'Total diff',sum(transport_change_rate(1:nlevdecomp,k)*dzsoi_decomp(1:nlevdecomp))*actual_dt/2,'Surf adv',surf_adv_step(k),'Surf equil',surf_equil_step(k)
    enddo
  
  
    surf_flux = surf_flux + surf_equil_step + surf_adv_step
    lat_flux  = lat_flux  + lat_flux_step
  endif ! if converged

end subroutine run_column_onestep
  
#endif

! Should make sure this is available when alquimia is turned off/not compiled in case we want to track e.g. salinity without BGC
subroutine advection_diffusion(conc_trcr,adv_flux,diffus,source,surf_bc,dtime,conc_change_rate)
  ! Advection and diffusion for a single tracer in one column given diffusion coefficient, flow, and source-sink terms
  ! Based on SoilLittVertTranspMod, which implements S. V. Patankar, Numerical Heat Transfer and Fluid Flow, Series in Computational Methods in Mechanics and Thermal Sciences, Hemisphere Publishing Corp., 1980. Chapter 5
  ! Not sure if this belongs here or somewhere else. Is it bad to do this in the EMI subroutine?

  use elm_varpar       , only : nlevdecomp
  use elm_varcon       , only : zsoi, zisoi, dzsoi_decomp

  real(r8), intent(in) :: conc_trcr(1:nlevdecomp) ! Bulk concentration (e.g. mol/m3). Or should it be concentration in water??
  real(r8), intent(in) :: adv_flux(1:nlevdecomp+1)    ! (m/s), vertical into layer (down is negative)
  real(r8), intent(in) :: diffus(1:nlevdecomp)  ! diffusivity (m2/s)
  real(r8), intent(in) :: source(1:nlevdecomp)  ! Source term (mol/m3/s)
 
  real(r8), intent(in) :: surf_bc                 ! Surface boundary layer concentration (for infiltration)
  real(r8), intent(in) :: dtime                   ! Time step (s)
  real(r8), intent(out):: conc_change_rate(1:nlevdecomp) ! Bulk concentration (e.g. mol/m3/s). Or should it be concentration in water??

  ! Local variables
  real(r8) :: aaa                                                ! "A" function in Patankar
  real(r8) :: pe                                                 ! Pe for "A" function in Patankar
  real(r8) :: w_m1, w_p1                                         ! Weights for calculating harmonic mean of diffusivity
  real(r8) :: d_m1, d_p1                                         ! Harmonic mean of diffusivity
  real(r8) :: a_tri(0:nlevdecomp+1)      ! "a" vector for tridiagonal matrix
  real(r8) :: b_tri(0:nlevdecomp+1)      ! "b" vector for tridiagonal matrix
  real(r8) :: c_tri(0:nlevdecomp+1)      ! "c" vector for tridiagonal matrix
  real(r8) :: r_tri(0:nlevdecomp+1)      ! "r" vector for tridiagonal solution
  real(r8) :: d_p1_zp1(1:nlevdecomp+1)   ! diffusivity/delta_z for next j  (set to zero for no diffusion)
  real(r8) :: d_m1_zm1(1:nlevdecomp+1)   ! diffusivity/delta_z for previous j (set to zero for no diffusion)
  real(r8) :: f_p1(1:nlevdecomp+1)       ! water flux for next j
  real(r8) :: f_m1(1:nlevdecomp+1)       ! water flux for previous j
  real(r8) :: pe_p1(1:nlevdecomp+1)      ! Peclet # for next j
  real(r8) :: pe_m1(1:nlevdecomp+1)      ! Peclet # for previous j
  real(r8) :: dz_node(1:nlevdecomp+1)                            ! difference between nodes
  real(r8) :: a_p_0
  real(r8) :: conc_after(0:nlevdecomp+1)
  real(r8) :: rho(1:nlevdecomp)     ! Water density (bulk) in layer
  
  integer :: j, info
  
  ! Statement function
  aaa (pe) = max (0._r8, (1._r8 - 0.1_r8 * abs(pe))**5)  ! "A" function from Patankar, Table 5.2, pg 95

  rho(1:nlevdecomp) = 1.0_r8 ! Placeholder in case we want to account for varying water content

  ! Set the distance between the node and the one ABOVE it   
  dz_node(1) = zsoi(1)
  do j = 2,nlevdecomp+1
     dz_node(j)= zsoi(j) - zsoi(j-1)
  enddo

  ! write(iulog,*) 'adv_flux',adv_flux(1:nlevdecomp+1)
  ! write(iulog,*) 'diffus',diffus(1:nlevdecomp)
  ! write(iulog,*) 'source',source(1:nlevdecomp)

  ! Calculate the D and F terms in the Patankar algorithm
  ! d: diffusivity
  ! f: flow
  ! m: layer above
  ! p: layer below
  ! pe: Peclet number (ratio of convection to diffusion)
  do j = 1,nlevdecomp
    if (j == 1) then
      d_m1_zm1(j) = 0._r8
      w_p1 = (zsoi(j+1) - zisoi(j)) / dz_node(j+1)
      if ( diffus(j+1) > 0._r8 .and. diffus(j) > 0._r8) then
        d_p1 = 1._r8 / ((1._r8 - w_p1) / diffus(j) + w_p1 / diffus(j+1)) ! Harmonic mean of diffus
      else
        d_p1 = 0._r8
      endif
      d_p1_zp1(j) = d_p1 / dz_node(j+1)
      f_m1(j) = adv_flux(j)  ! Include infiltration here
      f_p1(j) = adv_flux(j+1)
      pe_m1(j) = 0._r8
      pe_p1(j) = f_p1(j) / d_p1_zp1(j) ! Peclet #
    elseif (j == nlevdecomp) then
        ! At the bottom, assume no gradient in d_z (i.e., they're the same)
        w_m1 = (zisoi(j-1) - zsoi(j-1)) / dz_node(j)
        if ( diffus(j) > 0._r8 .and. diffus(j-1) > 0._r8) then
          d_m1 = 1._r8 / ((1._r8 - w_m1) / diffus(j) + w_m1 / diffus(j-1)) ! Harmonic mean of diffus
        else
          d_m1 = 0._r8
        endif
        d_m1_zm1(j) = d_m1 / dz_node(j)
        d_p1_zp1(j) = d_m1_zm1(j) ! Set to be the same
        f_m1(j) = adv_flux(j)
        !f_p1(j) = adv_flux(j+1)
        f_p1(j) = 0._r8
        pe_m1(j) = f_m1(j) / d_m1_zm1(j) ! Peclet #
        pe_p1(j) = f_p1(j) / d_p1_zp1(j) ! Peclet #
    else
        ! Use distance from j-1 node to interface with j divided by distance between nodes
        w_m1 = (zisoi(j-1) - zsoi(j-1)) / dz_node(j)
        if ( diffus(j-1) > 0._r8 .and. diffus(j) > 0._r8) then
          d_m1 = 1._r8 / ((1._r8 - w_m1) / diffus(j) + w_m1 / diffus(j-1)) ! Harmonic mean of diffus
        else
          d_m1 = 0._r8
        endif
        w_p1 = (zsoi(j+1) - zisoi(j)) / dz_node(j+1)
        if ( diffus(j+1) > 0._r8 .and. diffus(j) > 0._r8) then
          d_p1 = 1._r8 / ((1._r8 - w_p1) / diffus(j) + w_p1 / diffus(j+1)) ! Harmonic mean of diffus
        else
          d_p1 = (1._r8 - w_p1) * diffus(j) + w_p1 * diffus(j+1) ! Arithmetic mean of diffus
        endif
        d_m1_zm1(j) = d_m1 / dz_node(j)
        d_p1_zp1(j) = d_p1 / dz_node(j+1)
        f_m1(j) = adv_flux(j)
        f_p1(j) = adv_flux(j+1)
        pe_m1(j) = f_m1(j) / d_m1_zm1(j) ! Peclet #
        pe_p1(j) = f_p1(j) / d_p1_zp1(j) ! Peclet #
    end if
  enddo ! j; nlevdecomp


  ! Calculate the tridiagonal coefficients
  ! Coefficients of tridiagonal problem: a_i*x_(i-1) + b_i*(x_i) + c_i*x_(i+1) = r_i
  ! Here, this is equivalent to Patankar equation 5.56 and 5.57 (but in one dimension):
  ! a_P*phi_P = a_E*phi_E + a_W*phi_W + b [phi is concentration, = x in tridiagonal]. Converting East/West to above/below
  ! -> -a_E*phi_E + a_P*phi_P - a_W+phi_W = b
  ! -a_tri = a_above = D_above*A(Pe)+max(-F_above,0); D_above=diffus_above/dz
  ! b_tri = a_above+a_below+rho*dz/dt
  ! -c_tri = D_below*A(Pe)+max(F_below,0); D_below = diffus_below/dz
  ! r_tri = b = source_const*dz + conc*rho*dz/dt
  do j = 0,nlevdecomp +1

    if (j > 0 .and. j < nlevdecomp+1) then
        a_p_0 =  dzsoi_decomp(j) / dtime * rho(j) ! Should this be multiplied by layer water content (for rho)?
    endif

    if (j == 0) then ! top layer (atmosphere)
        a_tri(j) = 0._r8
        b_tri(j) = 1._r8
        c_tri(j) = -1._r8
        r_tri(j) = 0._r8
    elseif (j == 1) then
        a_tri(j) = -(d_m1_zm1(j) * aaa(pe_m1(j)) + max( f_m1(j), 0._r8)) ! Eqn 5.47 Patankar
        c_tri(j) = -(d_p1_zp1(j) * aaa(pe_p1(j)) + max(-f_p1(j), 0._r8))
        b_tri(j) = -a_tri(j) - c_tri(j) + a_p_0
        ! r_tri includes infiltration assuming same concentration as top layer. May want to change to either provide upper boundary condition or include in source term
        ! r_tri(j) = source(j) * dzsoi_decomp(j) + (a_p_0 - adv_flux(j)) * conc_trcr(j)
        r_tri(j) = source(j) * dzsoi_decomp(j) + a_p_0 * conc_trcr(j)
        if(adv_flux(j)<0) then ! downward flow (infiltration)
           r_tri(j) = r_tri(j) - adv_flux(j)*surf_bc
          !  write(iulog,*),__LINE__,adv_flux(j),surf_bc,adv_flux(j)*surf_bc
        else ! upward flow to the surface
          r_tri(j) = r_tri(j) - adv_flux(j)*conc_trcr(j)
          ! write(iulog,*),__LINE__,adv_flux(j),conc_trcr(j),adv_flux(j)*conc_trcr(j)
        endif
        
    elseif (j < nlevdecomp+1) then
        a_tri(j) = -(d_m1_zm1(j) * aaa(pe_m1(j)) + max( f_m1(j), 0._r8)) ! Eqn 5.47 Patankar
        c_tri(j) = -(d_p1_zp1(j) * aaa(pe_p1(j)) + max(-f_p1(j), 0._r8))
        b_tri(j) = -a_tri(j) - c_tri(j) + a_p_0
        r_tri(j) = source(j) * dzsoi_decomp(j) + a_p_0 * conc_trcr(j) ! Eq. 5.57
    else ! j==nlevdecomp+1; 0 concentration gradient at bottom
        a_tri(j) = -1._r8
        b_tri(j) = 1._r8
        c_tri(j) = 0._r8 
        r_tri(j) = 0._r8
    endif
  enddo ! j; nlevdecomp

  ! write(iulog,'(11a18)'),'a','b','c','r','ap0','pe_m','pe_p','f_m','f_p','d_m','d_p'
  ! j=0
  ! write(iulog,'(i3,4e18.9)'),j,a_tri(j),b_tri(j),c_tri(j),r_tri(j)
  ! do j=1,nlevdecomp
  !   write(iulog,'(i3,11e18.9)'),j,a_tri(j),b_tri(j),c_tri(j),r_tri(j),dzsoi_decomp(j) / dtime * rho(j) ,pe_m1(j),pe_p1(j),f_m1(j),f_p1(j),d_m1_zm1(j)*dz_node(j),d_p1_zp1(j)*dz_node(j+1)
  ! enddo
  ! j=nlevdecomp+1
  ! write(iulog,'(i3,4e18.9)'),j,a_tri(j),b_tri(j),c_tri(j),r_tri(j)

  ! Solve for the concentration profile for this time step
  ! call Tridiagonal(0, nlevdecomp+1, 0, a_tri, b_tri, c_tri, r_tri, conc_after)
  ! This is the LAPACK tridiagonal solver which gave more accurate results in my testing
  call dgtsv( nlevdecomp+2, 1, c_tri(0:nlevdecomp), b_tri, a_tri(1:nlevdecomp+1),  & 
              r_tri, nlevdecomp+2, info )

  if(info < 0) call endrun(msg='dgtsv error in adv_diff line __LINE__: illegal argument')
  if(info > 0) call endrun(msg='dgtsv error in adv_diff line __LINE__: singular matrix')
  conc_after = r_tri

  ! write(iulog,*),'conc_before',conc_trcr
  ! write(iulog,*),'conc_after',conc_after
  ! write(iulog,*),'Diff=',sum((conc_after(1:nlevdecomp)-conc_trcr)*dzsoi_decomp)
  ! write(iulog,*),'Flow',adv_flux(1:nlevdecomp+1)
  ! write(iulog,*),'Diffus',diffus
  ! write(iulog,*),'dz',dzsoi_decomp
  ! write(iulog,*),'dznode',dz_node

  conc_change_rate = (conc_after(1:nlevdecomp)-conc_trcr)/dtime

end subroutine advection_diffusion


  !-----------------------------------------------------------------------
! Modified to operate on a single column instead of passing bounds and filters
subroutine Tridiagonal (lbj, ubj, jtop, a, b, c, r, u)
  !
  ! !DESCRIPTION:
  ! Tridiagonal matrix solution
  ! A x = r
  ! where x and r are vectors

  !
  ! !ARGUMENTS:
  implicit none

  integer           , intent(in)    :: lbj, ubj                                 ! lbinning and ubing level indices
  integer           , intent(in)    :: jtop         ! top level for each column [col]
  real(r8)          , intent(in)    :: a(lbj:ubj)    ! "a" left off diagonal of tridiagonal matrix [col , j]
  real(r8)          , intent(in)    :: b(lbj:ubj)    ! "b" diagonal column for tridiagonal matrix [col  , j]
  real(r8)          , intent(in)    :: c(lbj:ubj)    ! "c" right off diagonal tridiagonal matrix [col   , j]
  real(r8)          , intent(in)    :: r(lbj:ubj)    ! "r" forcing term of tridiagonal matrix [col      , j]
  real(r8)          , intent(inout) :: u(lbj:ubj)    ! solution [col                                    , j]
                                                                                !
  integer                           :: j                                 ! indices

  real(r8)                          :: gam(lbj:ubj)     ! temporary
  real(r8)                          :: bet            ! temporary

  !-----------------------------------------------------------------------

  bet = b(jtop)

  do j = lbj, ubj
    if (j >= jtop) then
      if (j == jtop) then
        u(j) = r(j) / bet
      else
        gam(j) = c(j-1) / bet
        bet = b(j) - a(j) * gam(j)
        u(j) = (r(j) - a(j)*u(j-1)) / bet
      end if
    end if
  end do

  do j = ubj-1,lbj,-1
    if (j >= jtop) then
      u(j) = u(j) - gam(j+1) * u(j+1)
    end if
  end do

end subroutine Tridiagonal

end module ExternalModelAlquimiaMod
