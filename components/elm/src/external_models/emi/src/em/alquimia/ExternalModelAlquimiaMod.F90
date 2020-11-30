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
  use EMI_ColumnDataType_Constants
  use EMI_WaterFluxType_Constants
  use EMI_WaterStateType_Constants
  use EMI_CNCarbonStateType_Constants
  use EMI_CNNitrogenStateType_Constants
  use EMI_CNNitrogenFluxType_Constants
  use EMI_CNCarbonFluxType_Constants
  use EMI_ChemStateType_Constants
  use EMI_ColumnWaterStateType_Constants

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
    
    ! Solve data returned to land model
    integer :: index_e2l_state_decomp_cpools
    integer :: index_e2l_state_decomp_npools
    integer :: index_e2l_flux_hr
    integer :: index_e2l_state_nh4
    integer :: index_e2l_state_no3

    integer :: index_e2l_flux_Nimm
    integer :: index_e2l_flux_Nimp
    integer :: index_e2l_flux_Nmin

    integer :: index_e2l_flux_plantNO3uptake
    integer :: index_e2l_flux_plantNH4uptake

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
    id                                             = E2L_FLUX_HETEROTROPHIC_RESP_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_hr              = index
    
    id                                             = E2L_STATE_NH4_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_nh4              = index
    
    id                                             = E2L_STATE_NO3_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_no3              = index

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
                                            
    use clm_varctl, only : alquimia_inputfile,alquimia_engine_name,alquimia_IC_name,alquimia_handsoff

    use clm_varpar            , only : alquimia_num_primary, alquimia_num_minerals,&
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
    

    ! Todo: Keep track of other (non-SOMDEC) reactions too somehow
    
    ! Initial condition. The zero length for constraints suggest that it will be read in from input file?
    call AllocateAlquimiaGeochemicalCondition(len_trim(ic_name,C_INT),0,0,this%chem_ic)
    call f_c_string_ptr(ic_name,this%chem_ic%name,len_trim(ic_name)+1)
    
    ! Allocate alquimia's data structures. One copy per processor which will be written into as needed
    call AllocateAlquimiaState(this%chem_sizes, this%chem_state)
    call AllocateAlquimiaProperties(this%chem_sizes, this%chem_properties)
    call AllocateAlquimiaAuxiliaryData(this%chem_sizes, this%chem_aux_data)
    call AllocateAlquimiaAuxiliaryOutputData(this%chem_sizes, this%chem_aux_output)
    

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
use clm_varctl                , only : iulog
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

  use clm_varpar, only : nlevdecomp

  class(em_alquimia_type)              :: this
  integer              , intent(in)    :: clump_rank
  class(emi_data_list) , intent(in)    :: l2e_list
  class(emi_data_list) , intent(inout) :: e2l_list
  type(bounds_type)    , intent (in)   :: bounds_clump

  real(r8) , pointer, dimension(:,:)   ::  porosity_l2e, dz, h2o_liqvol
  real(r8) , pointer, dimension(:,:)   ::  water_density_e2l,aqueous_pressure_e2l
  real(r8) , pointer, dimension(:,:,:) ::  total_mobile_e2l
  real(r8) , pointer, dimension(:,:,:) ::  total_immobile_e2l
  real(r8) , pointer, dimension(:,:,:) ::  mineral_volume_fraction_e2l
  real(r8) , pointer, dimension(:,:,:) ::  mineral_specific_surface_area_e2l
  real(r8) , pointer, dimension(:,:,:) ::  surface_site_density_e2l
  real(r8) , pointer, dimension(:,:,:) ::  cation_exchange_capacity_e2l
  real(r8) , pointer, dimension(:,:,:) ::  aux_doubles_e2l
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

#ifdef USE_ALQUIMIA_LIB
  do fc = 1, num_soilc
    c = filter_soilc(fc)

      do j = 1, nlevdecomp

          
          ! Initialize the state for the cell
          this%chem_properties%volume = dz(c,j)
          this%chem_properties%saturation = h2o_liqvol(c,j)/porosity_l2e(c,j)
          this%chem_state%water_density = 1.0e3_r8
          this%chem_state%porosity = porosity_l2e(c,j)
          this%chem_state%aqueous_pressure = 101325.0
          this%chem_state%temperature = 250.0_r8 - 273.15 ! Temperature may not have been initialized yet

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
                                        total_mobile_e2l,&
                                        total_immobile_e2l,&
                                        mineral_volume_fraction_e2l,&
                                        mineral_specific_surface_area_e2l,&
                                        surface_site_density_e2l,&
                                        cation_exchange_capacity_e2l,&
                                        aux_doubles_e2l,&
                                        aux_ints_e2l) 
          
      enddo
  enddo
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
    integer                              :: c,fc,j,l,poolnum
    integer   , pointer                  :: filter_soilc(:)
    integer                              :: num_soilc
    integer                              :: max_cuts
    real(r8) , pointer, dimension(:,:,:)    :: soilcarbon_l2e,soilcarbon_e2l 
    real(r8) , pointer, dimension(:,:,:)    :: soilnitrogen_l2e,soilnitrogen_e2l 
    real(r8) , pointer, dimension(:,:,:)    :: decomp_k
    real(r8) , pointer, dimension(:,:)    :: hr_e2l , temperature, h2o_liqvol
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
    real(r8)                            :: CO2_before
    real(r8)                            :: molperL_to_molperm3
    real(r8), parameter                 :: minval = 1.e-30_r8 ! Minimum value to pass to PFLOTRAN to avoid numerical errors with concentrations of 0
    
    ! Setting these to the values in PFLOTRAN clm_rspfuncs.F90
    real(r8), parameter :: natomw = 14.0067d0 ! Value in clmvarcon is 14.007
    real(r8), parameter :: catomw = 12.0110d0 ! Value in clmvarcon is 12.011
    
    character(kind=C_CHAR,len=kAlquimiaMaxStringLength) :: status_message
    procedure(ReactionStepOperatorSplit), pointer :: engine_ReactionStepOperatorSplit
    procedure(GetAuxiliaryOutput), pointer   :: engine_getAuxiliaryOutput
    real (c_double), pointer :: alquimia_mobile_data(:), alquimia_immobile_data(:), alquimia_rates_data(:)

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
    call e2l_list%GetPointerToReal2D(this%index_e2l_flux_hr , hr_e2l) ! (gC/m3/s)
    
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_no3 , no3_e2l) ! gN/m3
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_nh4 , nh4_e2l) ! gN/m3

    call e2l_list%GetPointerToReal2D(this%index_e2l_flux_Nimm , Nimm_e2l) ! gN/m3/s
    call e2l_list%GetPointerToReal2D(this%index_e2l_flux_Nimp , Nimp_e2l) ! gN/m3/s
    call e2l_list%GetPointerToReal2D(this%index_e2l_flux_Nmin , Nmin_e2l) ! gN/m3/s

    call e2l_list%GetPointerToReal2D(this%index_e2l_flux_plantNO3uptake , plantNO3uptake_e2l) ! gN/m3/s
    call e2l_list%GetPointerToReal2D(this%index_e2l_flux_plantNH4uptake , plantNH4uptake_e2l) ! gN/m3/s

    ! Alquimia state data on ELM side
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_watsatc       , porosity_l2e     )
    call e2l_list%GetPointerToReal2D(this%index_e2l_water_density, water_density_e2l)
    call l2e_list%GetPointerToReal2D(this%index_l2e_water_density, water_density_l2e)
    call e2l_list%GetPointerToReal2D(this%index_e2l_aqueous_pressure, aqueous_pressure_e2l)
    call l2e_list%GetPointerToReal2D(this%index_l2e_aqueous_pressure, aqueous_pressure_l2e)
    call e2l_list%GetPointerToReal3D(this%index_e2l_total_mobile, total_mobile_e2l)
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


    ! First check if pools have been mapped between ELM and Alquimia
    if(.not. associated(this%carbon_pool_mapping)) then
      call this%map_alquimia_pools()
    endif
    
     ! Run the reactions engine for a step. Alquimia works on one cell at a time
     ! TODO: Transport needs to be integrated somehow. 
    do fc = 1, num_soilc
      c = filter_soilc(fc)
         do j = 1, nlevdecomp  
             ! First copy Alquimia data from ELM to alquimia data structure
              call this%copy_ELM_to_Alquimia(c,j,water_density_l2e,&
                                            aqueous_pressure_l2e,&
                                            total_mobile_l2e,&
                                            total_immobile_l2e,&
                                            mineral_volume_fraction_l2e,&
                                            mineral_specific_surface_area_l2e,&
                                            surface_site_density_l2e,&
                                            cation_exchange_capacity_l2e,&
                                            aux_doubles_l2e,&
                                            aux_ints_l2e) 
             
            ! Update properties from ELM. Does it make more sense to put these in above subroutine call?
             this%chem_state%porosity =    porosity_l2e(c,j)
             this%chem_state%temperature = temperature(c,j) - 273.15
             this%chem_properties%volume = dz(c,j)
             this%chem_properties%saturation = max(h2o_liqvol(c,j)/porosity_l2e(c,j),0.01) ! Set minimum saturation to stop concentrations from blowing up at low soil moisture
             
             ! Conversion from aqueous Alquimia units (mol/L water) to ELM units (mol/m3 bulk)
             molperL_to_molperm3 = 1000.0*this%chem_state%porosity*this%chem_properties%saturation
             ! Should we read aqueous pressure in from ELM?
             ! Set volume and saturation?
             ! Saturation determines volume of water and concentration of aqueous species
             ! Not sure if volume actually matters

             ! Set rate constants in properties.aqueous_kinetic_rate_cnst if in hands_on mode
             ! Should we check if we are in hands_on mode?
             ! Discussions with Peter, Fengming, Scott, Ethan suggested we should avoid setting rate constants this way and let PFLOTRAN handle it
             ! to make linking identical BGC with ATS or other models easier
             call c_f_pointer(this%chem_properties%aqueous_kinetic_rate_cnst%data, alquimia_rates_data, (/this%chem_properties%aqueous_kinetic_rate_cnst%size/))
             call c_f_pointer(this%chem_state%total_mobile%data, alquimia_mobile_data, (/this%chem_sizes%num_primary/))
             call c_f_pointer(this%chem_state%total_immobile%data, alquimia_immobile_data, (/this%chem_sizes%num_primary/))

             do poolnum=1,ndecomp_pools
               if(this%pool_reaction_mapping(poolnum)>0) then
                !  write(iulog, '(a,i3,a,i3,1x,a,e12.3)'), 'Setting pool',poolnum,' -> alquimia reaction ',this%pool_reaction_mapping(poolnum),'kinetic rate const to',decomp_k(c,j,poolnum)
                 alquimia_rates_data(this%pool_reaction_mapping(poolnum)) = decomp_k(c,j,poolnum)
               endif
             enddo
             
             
             ! Set soil carbon and nitrogen from land model
             ! Convert soil C,N from g/m3 to mol/m3. Assumes pool is defined as immobile, not aqueous
             ! May need to deal with case that pools are all zero (initial condition) which PFLOTRAN will not be able to solve.
             
            !  write(iulog,*),'Before solve'
             do poolnum=1,ndecomp_pools
               if(this%carbon_pool_mapping(poolnum)>0) &  
                 alquimia_immobile_data(this%carbon_pool_mapping(poolnum)) = max(soilcarbon_l2e(c,j,poolnum)/catomw,minval)
               ! Separate N pool only exists if floating CN ratio
                !  write(iulog,*),poolnum,soilnitrogen_l2e(c,j,poolnum)
               if(decomp_cascade_con%floating_cn_ratio_decomp_pools(poolnum) .and. this%nitrogen_pool_mapping(poolnum)>0) &
                  alquimia_immobile_data(this%nitrogen_pool_mapping(poolnum)) = max(soilnitrogen_l2e(c,j,poolnum)/natomw,minval/20)
             enddo
             
             CO2_before = alquimia_immobile_data(this%CO2_pool_number)*catomw + &
                          alquimia_mobile_data(this%CO2_pool_number)*catomw*(molperL_to_molperm3)

             ! Copy dissolved nitrogen species. Units need to be converted from gN/m3 to M/L. Currently assuming saturated porosity
             
             if(this%NO3_pool_number>0) alquimia_mobile_data(this%NO3_pool_number) = max(no3_l2e(c,j)/natomw/(molperL_to_molperm3),minval)
             if(this%NH4_pool_number>0) alquimia_mobile_data(this%NH4_pool_number) = max(nh4_l2e(c,j)/natomw/(molperL_to_molperm3),minval)

             ! Set rate constant based on plant N demand. Convert from gN/m3/s to mol/L/s
             ! Also scale rates by relative concentrations of NO3 and NH4 so total uptake doesn't exceed demand
             ! Assumes alquimia is running in hands-off mode: Biomass term of N uptake microbial reaction is set to plant NO3 or NH4 demand
             ! This assumes the rate constant of the reaction is set to 1 in the input deck!
            if(this%plantNH4demand_pool_number>0) then
              alquimia_immobile_data(this%plantNH4demand_pool_number) = plantNdemand_l2e(c,j)/natomw/(molperL_to_molperm3)
              if(this%NO3_pool_number>0 .and. this%NH4_pool_number>0 .and. (no3_l2e(c,j)+nh4_l2e(c,j)>0)) &
                  alquimia_immobile_data(this%plantNH4demand_pool_number) = alquimia_immobile_data(this%plantNH4demand_pool_number)*nh4_l2e(c,j)/(nh4_l2e(c,j)+no3_l2e(c,j))
              alquimia_immobile_data(this%plantNH4demand_pool_number) = max(alquimia_immobile_data(this%plantNH4demand_pool_number),minval)
            endif
            if(this%plantNO3demand_pool_number>0) then
              alquimia_immobile_data(this%plantNO3demand_pool_number) = plantNdemand_l2e(c,j)/natomw/(molperL_to_molperm3)
              if(this%NO3_pool_number>0 .and. this%NH4_pool_number>0 .and. (no3_l2e(c,j)+nh4_l2e(c,j)>0)) &
                alquimia_immobile_data(this%plantNO3demand_pool_number) = alquimia_immobile_data(this%plantNO3demand_pool_number)*no3_l2e(c,j)/(nh4_l2e(c,j)+no3_l2e(c,j))
              alquimia_immobile_data(this%plantNO3demand_pool_number) = max(alquimia_immobile_data(this%plantNO3demand_pool_number),minval)
            endif

           
            ! Reset diagnostic N immobilization, mineralization
             if(this%Nimm_pool_number>0) alquimia_immobile_data(this%Nimm_pool_number) = minval
             if(this%Nimp_pool_number>0) alquimia_immobile_data(this%Nimp_pool_number) = minval
             if(this%Nmin_pool_number>0) alquimia_immobile_data(this%Nmin_pool_number) = minval

             if(this%plantNO3uptake_pool_number>0) alquimia_immobile_data(this%plantNO3uptake_pool_number) = minval
             if(this%plantNO3uptake_pool_number>0) alquimia_mobile_data(this%plantNO3uptake_pool_number) = minval
             if(this%plantNH4uptake_pool_number>0) alquimia_immobile_data(this%plantNH4uptake_pool_number) = minval
             if(this%plantNH4uptake_pool_number>0) alquimia_mobile_data(this%plantNH4uptake_pool_number) = minval

              ! Step the chemistry solver, including timestep cutting capability
              call run_onestep(this, c,j, dt,0,max_cuts)
              if(max_cuts>3) write(iulog,'(a,i2,a,2i3)'),"Alquimia converged after",max_cuts," cuts",c,j

              ! Save back to ELM
              call this%copy_Alquimia_to_ELM(c,j,water_density_e2l,&
                                              aqueous_pressure_e2l,&
                                              total_mobile_e2l,&
                                              total_immobile_e2l,&
                                              mineral_volume_fraction_e2l,&
                                              mineral_specific_surface_area_e2l,&
                                              surface_site_density_e2l,&
                                              cation_exchange_capacity_e2l,&
                                              aux_doubles_e2l,&
                                              aux_ints_e2l) 

              ! Set updated land model values. Should this be moved into copy subroutine?
              ! Convert from mol/m3 to gC/m2
              do poolnum=1,ndecomp_pools
                if(this%carbon_pool_mapping(poolnum)>0) &
                  soilcarbon_e2l(c,j,poolnum) = alquimia_immobile_data(this%carbon_pool_mapping(poolnum))*catomw
                ! Separate N pool only exists if floating CN ratio
                if(decomp_cascade_con%floating_cn_ratio_decomp_pools(poolnum) .and. this%nitrogen_pool_mapping(poolnum)>0) then
                   soilnitrogen_e2l(c,j,poolnum) = alquimia_immobile_data(this%nitrogen_pool_mapping(poolnum))*natomw
                 elseif (this%carbon_pool_mapping(poolnum)>0) then
                   ! Calculate from CN ratio and C pool
                   soilnitrogen_e2l(c,j,poolnum) = soilcarbon_e2l(c,j,poolnum)/decomp_cascade_con%initial_cn_ratio(poolnum)
                endif
                
                ! write(iulog,*),poolnum,soilnitrogen_e2l(c,j,poolnum)
              enddo
              ! Sum together mobile and immobile pools
              ! hr_e2l goes to hr_vr (gC/m3/s)
              if(this%CO2_pool_number>0) then 
                hr_e2l(c,j) = - CO2_before
                ! Immobile: Convert from mol/m3 to gC/m3/s
                hr_e2l(c,j) = hr_e2l(c,j) + alquimia_immobile_data(this%CO2_pool_number)*catomw
                ! Mobile: convert from mol/L to gC/m3/s. mol/L*gC/mol*1000L/m3*porosity
                hr_e2l(c,j) = hr_e2l(c,j) + alquimia_mobile_data(this%CO2_pool_number)*catomw*(molperL_to_molperm3)
                hr_e2l(c,j) = hr_e2l(c,j)/dt
              endif
              
              if(this%NO3_pool_number>0) no3_e2l(c,j) = alquimia_mobile_data(this%NO3_pool_number)*natomw*(molperL_to_molperm3)
              if(this%NH4_pool_number>0) nh4_e2l(c,j) = alquimia_mobile_data(this%NH4_pool_number)*natomw*(molperL_to_molperm3)

              if(this%Nimm_pool_number>0) Nimm_e2l(c,j) = alquimia_immobile_data(this%Nimm_pool_number)*natomw/dt
              if(this%Nimp_pool_number>0) Nimp_e2l(c,j) = alquimia_immobile_data(this%Nimp_pool_number)*natomw/dt
              ! Nmin will be added to the NH4 pool elsewhere in ELM so skip that for now
              ! if(this%Nmin_pool_number>0) Nmin_e2l(c,j) = alquimia_immobile_data(this%Nmin_pool_number)*natomw/dt

              ! PFLOTRAN may use an aqueous tracer to model plant N uptake if defining using Microbial reaction
              if(this%plantNO3uptake_pool_number>0) plantNO3uptake_e2l(c,j) = (alquimia_immobile_data(this%plantNO3uptake_pool_number)-minval)*natomw/dt + &
                                                (alquimia_mobile_data(this%plantNO3uptake_pool_number)-minval)*natomw*(molperL_to_molperm3)/dt
                if(this%plantNH4uptake_pool_number>0) plantNH4uptake_e2l(c,j) = (alquimia_immobile_data(this%plantNH4uptake_pool_number)-minval)*natomw/dt + &
                                                (alquimia_mobile_data(this%plantNH4uptake_pool_number)-minval)*natomw*(molperL_to_molperm3)/dt



              ! Todo: Add C check
              ! Note: Generates errors if not multiplied by layer volume (imbalance on the order of 1e-8 gN/m3)
              ! Note: Generates error after restart at precision of 1e-9. But doesn't set off N conservation errors in model when precision here is relaxed.
              if(abs(sum(soilnitrogen_l2e(c,j,:))+no3_l2e(c,j)+nh4_l2e(c,j)-&
                        (sum(soilnitrogen_e2l(c,j,:))+no3_e2l(c,j)+nh4_e2l(c,j)+plantNO3uptake_e2l(c,j)*dt+plantNH4uptake_e2l(c,j)*dt))*this%chem_properties%volume>1e-5) then
                write(iulog,'(a,1x,i3,a,i5)'),'Nitrogen imbalance after alquimia solve step in layer',j,' Column ',c,__FILE__,__LINE__
                call print_pools(this,c,j)
                
                write(iulog,'(a25,3e20.8)'),'Total N: ', sum(soilnitrogen_l2e(c,j,:))+no3_l2e(c,j)+nh4_l2e(c,j),&
                                            sum(soilnitrogen_e2l(c,j,:))+no3_e2l(c,j)+nh4_e2l(c,j)+plantNH4uptake_e2l(c,j)*dt+plantNO3uptake_e2l(c,j)*dt,&
                                            sum(soilnitrogen_e2l(c,j,:))+no3_e2l(c,j)+nh4_e2l(c,j)+plantNO3uptake_e2l(c,j)*dt+plantNH4uptake_e2l(c,j)*dt-(sum(soilnitrogen_l2e(c,j,:))+no3_l2e(c,j)+nh4_l2e(c,j))
                write(iulog,'(a25,3e20.8)'),'SON pools: ' ,sum(soilnitrogen_l2e(c,j,:)),sum(soilnitrogen_e2l(c,j,:)),sum(soilnitrogen_e2l(c,j,:)-soilnitrogen_l2e(c,j,:))
                write(iulog,'(a25,3e20.8)'),'NO3: ',no3_l2e(c,j),no3_e2l(c,j),no3_e2l(c,j)-no3_l2e(c,j)
                write(iulog,'(a25,3e20.8)'),'NH4: ',nh4_l2e(c,j),nh4_e2l(c,j),nh4_e2l(c,j)-nh4_l2e(c,j)
                write(iulog,'(a25,3e20.8)'),'Plant NO3, NH4 uptake: ',plantNO3uptake_e2l(c,j)*dt,plantNH4uptake_e2l(c,j)*dt,plantNO3uptake_e2l(c,j)*dt+plantNH4uptake_e2l(c,j)*dt
                call endrun(msg='N imbalance after alquimia solve')
              endif
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
    integer                              :: c,j ! Column, layer
    ! Pointer arrays that were previously mapped using EMI
    real(r8), pointer :: water_density(:,:), aqueous_pressure(:,:)
    real(r8), pointer :: total_mobile(:,:,:), total_immobile(:,:,:)
    real(r8), pointer :: mineral_volume_fraction(:,:,:), mineral_specific_surface_area(:,:,:)
    real(r8), pointer :: surface_site_density(:,:,:), cation_exchange_capacity(:,:,:), aux_doubles(:,:,:)
    integer , pointer :: aux_ints(:,:,:)

    real (c_double), pointer :: alquimia_data(:)
    integer (c_int)   , pointer :: alquimia_int_data(:)

    water_density(c,j) = this%chem_state%water_density
    aqueous_pressure(c,j) = this%chem_state%aqueous_pressure

    ! c_f_pointer just points an array to the right data, so it needs to be actually copied
    call c_f_pointer(this%chem_state%total_mobile%data, alquimia_data, (/this%chem_sizes%num_primary/))
    total_mobile(c,j,1:this%chem_sizes%num_primary)   = alquimia_data(1:this%chem_sizes%num_primary)
    call c_f_pointer(this%chem_state%total_immobile%data, alquimia_data, (/this%chem_sizes%num_primary/))
    total_immobile(c,j,1:this%chem_sizes%num_primary)   = alquimia_data(1:this%chem_sizes%num_primary)
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
    integer                              :: c,j ! Column, layer
    ! Pointer arrays that were previously mapped using EMI
    real(r8), pointer :: water_density(:,:), aqueous_pressure(:,:)
    real(r8), pointer :: total_mobile(:,:,:), total_immobile(:,:,:)
    real(r8), pointer :: mineral_volume_fraction(:,:,:), mineral_specific_surface_area(:,:,:)
    real(r8), pointer :: surface_site_density(:,:,:), cation_exchange_capacity(:,:,:), aux_doubles(:,:,:)
    integer , pointer :: aux_ints(:,:,:)

    real (c_double), pointer :: alquimia_data(:)
    integer (c_int)   , pointer :: alquimia_int_data(:)

    this%chem_state%water_density = water_density(c,j)
    this%chem_state%aqueous_pressure = aqueous_pressure(c,j)

    ! c_f_pointer just points an array to the right data, so it needs to be actually copied
    call c_f_pointer(this%chem_state%total_mobile%data, alquimia_data, (/this%chem_sizes%num_primary/))
    alquimia_data(1:this%chem_sizes%num_primary) = total_mobile(c,j,1:this%chem_sizes%num_primary)
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


    use clm_varpar, only : ndecomp_pools
    use CNDecompCascadeConType, only : decomp_cascade_con
    use clm_varctl, only : alquimia_IC_name,alquimia_CO2_name,&
        alquimia_NO3_name,alquimia_NH4_name,alquimia_Nimp_name,alquimia_Nmin_name,alquimia_Nimm_name,&
        alquimia_plantNO3uptake_name,alquimia_plantNH4uptake_name,alquimia_plantNO3demand_name,alquimia_plantNH4demand_name
    use clm_varpar, only : nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions

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
    alq_poolname = '1.0000e+00 NH4+  -> 1.0000e+00 Tracer2' ! Todo: Fix this! !
    this%plantNH4uptake_reaction_number = find_alquimia_pool(alq_poolname,name_list,this%chem_metadata%aqueous_kinetic_names%size)
    if(this%plantNH4uptake_reaction_number>0) then 
      write(iulog,'(a, i3, 1X, a)'),'ELM plant NH4+ uptake <-> Alquimia reaction',this%plantNH4uptake_reaction_number,trim(alq_poolname)
    else
      write(iulog,'(a,1x,a)'),'WARNING: No match for plant NH4+ uptake reaction',trim(alq_poolname)
    endif
    alq_poolname = '1.0000e+00 NO3-  -> 1.0000e+00 Tracer' ! Todo: Fix this! !
    this%plantNO3uptake_reaction_number = find_alquimia_pool(alq_poolname,name_list,this%chem_metadata%aqueous_kinetic_names%size)
    if(this%plantNO3uptake_reaction_number>0) then 
      write(iulog,'(a, i3, 1X, a)'),'ELM plant NO3- uptake <-> Alquimia reaction',this%plantNO3uptake_reaction_number,trim(alq_poolname)
    else
      write(iulog,'(a,1x,a)'),'WARNING: No match for plant NO3- uptake reaction',trim(alq_poolname)
    endif
  end subroutine map_alquimia_pools

  
  subroutine print_pools(this,c,j)

    use clm_varpar, only : ndecomp_pools,ndecomp_cascade_transitions
    use iso_c_binding, only : c_f_pointer, c_double
    use CNDecompCascadeConType, only : decomp_cascade_con

    implicit none

    class(em_alquimia_type)              :: this
    integer, intent(in) :: c,j

    integer :: poolnum
    character(len=256) :: poolname
    real (c_double), pointer :: alquimia_mobile_data(:), alquimia_immobile_data(:), alquimia_rates_data(:)

    call c_f_pointer(this%chem_state%total_immobile%data, alquimia_immobile_data, (/this%chem_sizes%num_primary/))
    call c_f_pointer(this%chem_state%total_mobile%data, alquimia_mobile_data, (/this%chem_sizes%num_primary/))
    call c_f_pointer(this%chem_properties%aqueous_kinetic_rate_cnst%data, alquimia_rates_data, (/this%chem_properties%aqueous_kinetic_rate_cnst%size/))

    write(iulog,*), "Carbon pool values: Immobile pools"
    do poolnum=1,ndecomp_pools
      poolname = trim(decomp_cascade_con%decomp_pool_name_history(poolnum))

      if(this%carbon_pool_mapping(poolnum)>0) then  
        write(iulog,'(a8,i4,e12.5)'), trim(poolname),this%carbon_pool_mapping(poolnum),alquimia_immobile_data(this%carbon_pool_mapping(poolnum))
      else
        write(iulog,*), trim(poolname),'Not in alquimia structure'
      endif

    enddo

    if(this%CO2_pool_number>0) then
      write(iulog,'(a8,i4,e12.5)'),'CO2',this%CO2_pool_number,alquimia_immobile_data(this%CO2_pool_number)
    endif

    write(iulog,*) "Nitrogen pool values: Immobile pools"
    do poolnum=1,ndecomp_pools
      if(decomp_cascade_con%floating_cn_ratio_decomp_pools(poolnum) .and. this%nitrogen_pool_mapping(poolnum)>0) then
          poolname = trim(decomp_cascade_con%decomp_pool_name_history(poolnum))
          write(iulog,'(a8,i4,e12.5)'),trim(poolname),this%nitrogen_pool_mapping(poolnum),alquimia_immobile_data(this%nitrogen_pool_mapping(poolnum))
      endif
    enddo

    if(this%NH4_pool_number>0) then
      write(iulog,'(a8,i4,e12.5)'),'NH4',this%NH4_pool_number,alquimia_immobile_data(this%NH4_pool_number)
    else
      write(iulog,*),'NH4 pool not in alquimia'
    endif
    if(this%NO3_pool_number>0) then
      write(iulog,'(a8,i4,e12.5)'),'NO3',this%NO3_pool_number,alquimia_immobile_data(this%NO3_pool_number)
    else
      write(iulog,*),'NO3 pool not in alquimia'
    endif
    


    write(iulog,*) "Carbon pool values: Aqueous pools"
    
    do poolnum=1,ndecomp_pools
      poolname = trim(decomp_cascade_con%decomp_pool_name_history(poolnum))

      if(this%carbon_pool_mapping(poolnum)>0) then  
        write(iulog,'(a8,i4,e12.5)') trim(poolname),this%carbon_pool_mapping(poolnum),alquimia_mobile_data(this%carbon_pool_mapping(poolnum))
      else
        write(iulog,*) trim(poolname),'Not in alquimia structure'
      endif

    enddo

    if(this%CO2_pool_number>0) then
      write(iulog,'(a8,i4,e12.5)'),'CO2',this%CO2_pool_number,alquimia_mobile_data(this%CO2_pool_number)
    endif

    write(iulog,*) "Nitrogen pool values: Aqueous pools"
    do poolnum=1,ndecomp_pools
      if(decomp_cascade_con%floating_cn_ratio_decomp_pools(poolnum) .and. this%nitrogen_pool_mapping(poolnum)>0) then
          poolname = trim(decomp_cascade_con%decomp_pool_name_history(poolnum))
          write(iulog,'(a8,i4,e12.5)'),trim(poolname),this%nitrogen_pool_mapping(poolnum),alquimia_mobile_data(this%nitrogen_pool_mapping(poolnum))
      endif
    enddo

    if(this%NH4_pool_number>0) then
      write(iulog,'(a8,i4,e12.5)'),'NH4',this%NH4_pool_number,alquimia_mobile_data(this%NH4_pool_number)
    else
      write(iulog,*),'NH4 pool not in alquimia'
    endif
    if(this%NO3_pool_number>0) then
      write(iulog,'(a8,i4,e12.5)'),'NO3',this%NO3_pool_number,alquimia_mobile_data(this%NO3_pool_number)
    else
      write(iulog,*),'NO3 pool not in alquimia'
    endif

    if(this%plantNO3demand_pool_number>0) then
      write(iulog,'(a,i4,e12.5)'),'Plant NO3 demand',this%plantNO3demand_pool_number,alquimia_immobile_data(this%plantNO3demand_pool_number)
    else
      write(iulog,*),'Plant NO3 demand pool not in alquimia'
    endif
    if(this%plantNH4demand_pool_number>0) then
      write(iulog,'(a,i4,e12.5)'),'Plant NH4 demand',this%plantNH4demand_pool_number,alquimia_immobile_data(this%plantNH4demand_pool_number)
    else
      write(iulog,*),'Plant NH4 demand pool not in alquimia'
    endif
    
    write(iulog,*) "Rate constants"
    do poolnum=1, ndecomp_cascade_transitions
      if(this%pool_reaction_mapping(decomp_cascade_con%cascade_donor_pool(poolnum))>0) then
        write(iulog,*),trim(decomp_cascade_con%cascade_step_name(poolnum)),alquimia_rates_data(this%pool_reaction_mapping(decomp_cascade_con%cascade_donor_pool(poolnum)))
      endif
    enddo

    if(this%plantNH4uptake_reaction_number>0) write(iulog,*),'Plant NH4 uptake',alquimia_rates_data(this%plantNH4uptake_reaction_number)
    if(this%plantNO3uptake_reaction_number>0) write(iulog,*),'Plant NO3 uptake',alquimia_rates_data(this%plantNO3uptake_reaction_number)

    write(iulog,*),'Porosity =',this%chem_state%porosity

  end subroutine print_pools

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
        call print_pools(this,c,j)
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
  
#endif

end module ExternalModelAlquimiaMod
