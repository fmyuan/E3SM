module ExternalModelATSMod

#ifdef USE_ATS_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module provides
  !
  use abortutils                   , only : endrun
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use shr_log_mod                  , only : errMsg => shr_log_errMsg
  use EMI_DataMod, only : emi_data_list, emi_data
  use mpp_varctl                   , only : iulog
  use ExternalModelBaseType        , only : em_base_type
  use decompMod                    , only : bounds_type
  use ExternalModelConstants
  use EMI_Atm2LndType_Constants
  use EMI_CanopyStateType_Constants
  use EMI_ColumnType_Constants
  use EMI_EnergyFluxType_Constants
  use EMI_Filter_Constants
  use EMI_Landunit_Constants
  use EMI_SoilHydrologyType_Constants
  use EMI_SoilStateType_Constants
  use EMI_TemperatureType_Constants
  use EMI_WaterFluxType_Constants
  use EMI_WaterStateType_Constants
  !
  !use ats_elm_interface            , only : ats_elm_type
  !
  implicit none
  !

  type, public, extends(em_base_type) :: em_ats_type
     ! ----------------------------------------------------------------------
     ! Indicies required during the initialization
     ! ----------------------------------------------------------------------
     integer :: index_l2e_init_col_active
     integer :: index_l2e_init_col_type
     integer :: index_l2e_init_col_landunit_index
     integer :: index_l2e_init_col_zi
     integer :: index_l2e_init_col_dz
     integer :: index_l2e_init_col_z
     integer :: index_l2e_init_col_area

     integer :: index_l2e_init_state_wtd
     integer :: index_l2e_init_state_soilp

     integer :: index_l2e_init_h2osoi_liq
     integer :: index_l2e_init_h2osoi_ice

     integer :: index_e2l_init_state_h2osoi_liq
     integer :: index_e2l_init_state_h2osoi_ice
     integer :: index_e2l_init_state_smp
     integer :: index_e2l_init_state_wtd

     integer :: index_e2l_init_flux_mflx_snowlyr_col
     integer :: index_l2e_init_flux_mflx_snowlyr_col

     integer :: index_l2e_init_landunit_type
     integer :: index_l2e_init_landunit_lakepoint
     integer :: index_l2e_init_landunit_urbanpoint

     integer :: index_l2e_init_parameter_watsatc
     integer :: index_l2e_init_parameter_hksatc
     integer :: index_l2e_init_parameter_bswc
     integer :: index_l2e_init_parameter_sucsatc
     integer :: index_l2e_init_parameter_effporosityc

     ! ----------------------------------------------------------------------
     ! Indicies required during timestepping
     ! ----------------------------------------------------------------------

     integer :: index_l2e_state_tsoil
     integer :: index_l2e_state_h2osoi_liq
     integer :: index_l2e_state_h2osoi_ice

     integer :: index_e2l_state_h2osoi_liq
     integer :: index_e2l_state_h2osoi_ice
     integer :: index_e2l_state_smp
     integer :: index_e2l_state_wtd

     integer :: index_l2e_flux_infil
     integer :: index_l2e_flux_et
     integer :: index_l2e_flux_dew
     integer :: index_l2e_flux_snow_sub
     integer :: index_l2e_flux_snowlyr
     integer :: index_l2e_flux_drainage

     integer :: index_e2l_flux_qrecharge

     integer :: index_l2e_filter_hydrologyc
     integer :: index_l2e_filter_num_hydrologyc

     integer :: index_l2e_column_zi

     ! IDs to indentify the conditions for ATS
     integer :: ats_cond_id_for_infil
     integer :: ats_cond_id_for_et
     integer :: ats_cond_id_for_dew
     integer :: ats_cond_id_for_drainage
     integer :: ats_cond_id_for_snow
     integer :: ats_cond_id_for_sublimation
     integer :: ats_cond_id_for_lateral_flux

     !type(ats_elm_type) :: ats_elm

   contains

     procedure, public :: Populate_L2E_Init_List  => EM_ATS_Populate_L2E_Init_List
     procedure, public :: Populate_E2L_Init_List  => EM_ATS_Populate_E2L_Init_List
     procedure, public :: Populate_L2E_List       => EM_ATS_Populate_L2E_List
     procedure, public :: Populate_E2L_List       => EM_ATS_Populate_E2L_List
     procedure, public :: Init                    => EM_ATS_Init
     procedure, public :: Solve                   => EM_ATS_OneStep
     procedure, public :: Finalize                => EM_ATS_Finalize
  end type em_ats_type

contains

  !------------------------------------------------------------------------
  subroutine EM_ATS_Populate_L2E_Init_List(this, l2e_init_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by ATS from ELM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_ats_type)                  :: this
    class(emi_data_list), intent(inout) :: l2e_init_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index
    !
    !
    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_INITIALIZATION_STAGE

    id                                         = L2E_STATE_WTD
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_wtd              = index

    id                                         = L2E_FLUX_RESTART_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_flux_mflx_snowlyr_col  = index

    id                                         = L2E_COLUMN_ACTIVE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_active             = index

    id                                         = L2E_COLUMN_TYPE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_type               = index

    id                                         = L2E_COLUMN_LANDUNIT_INDEX
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_landunit_index     = index

    id                                         = L2E_COLUMN_ZI
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_zi                 = index

    id                                         = L2E_COLUMN_DZ
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_dz                 = index

    id                                         = L2E_COLUMN_Z
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_z                  = index

    id                                         = L2E_COLUMN_AREA
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_area               = index

    id                                         = L2E_LANDUNIT_TYPE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_type          = index

    id                                         = L2E_LANDUNIT_LAKEPOINT
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_lakepoint     = index

    id                                         = L2E_LANDUNIT_URBANPOINT
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_urbanpoint    = index

    id                                         = L2E_PARAMETER_WATSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_watsatc      = index

    id                                         = L2E_PARAMETER_HKSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_hksatc       = index

    id                                         = L2E_PARAMETER_BSWC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_bswc         = index

    id                                         = L2E_PARAMETER_SUCSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_sucsatc      = index

    id                                         = L2E_PARAMETER_EFFPOROSITYC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_effporosityc = index

    id                                        = L2E_STATE_H2OSOI_LIQ_NLEVGRND
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_h2osoi_liq            = index

    id                                        = L2E_STATE_H2OSOI_ICE_NLEVGRND
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_h2osoi_ice            = index

    deallocate(em_stages)

  end subroutine EM_ATS_Populate_L2E_Init_List

  !------------------------------------------------------------------------
  subroutine EM_ATS_Populate_E2L_Init_List(this, e2l_init_list)
    !
    !
    ! !DESCRIPTION:
    ! Create a list of all variables to be returned by ATS from ELM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_ats_type)                   :: this
    class(emi_data_list) , intent(inout) :: e2l_init_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data)      , pointer       :: data
    integer              , pointer       :: em_stages(:)
    integer                              :: number_em_stages
    integer                              :: id
    integer                              :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_INITIALIZATION_STAGE

    id                                        = E2L_STATE_H2OSOI_LIQ
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_liq      = index

    id                                        = E2L_STATE_H2OSOI_ICE
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_ice      = index

    id                                        = E2L_STATE_SOIL_MATRIC_POTENTIAL
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_smp             = index

    id                                        = E2L_STATE_WTD
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_wtd             = index

    id                                        = E2L_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_flux_mflx_snowlyr_col = index

    deallocate(em_stages)

  end subroutine EM_ATS_Populate_E2L_Init_List

  !------------------------------------------------------------------------
  subroutine EM_ATS_Populate_L2E_List(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by ATS from ELM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_ats_type)                  :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index
    !
    !
    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_ATS_SOIL_HYDRO_STAGE

    id                                   = L2E_STATE_TSOIL_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_tsoil           = index

    id                                   = L2E_STATE_H2OSOI_LIQ_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_liq      = index

    id                                   = L2E_STATE_H2OSOI_ICE_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_ice      = index

    id                                   = L2E_FLUX_INFIL_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_infil            = index

    id                                   = L2E_FLUX_VERTICAL_ET_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_et               = index

    id                                   = L2E_FLUX_DEW_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_dew              = index

    id                                   = L2E_FLUX_SNOW_SUBLIMATION_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_snow_sub         = index

    id                                   = L2E_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_snowlyr          = index

    id                                   = L2E_FLUX_DRAINAGE_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_drainage         = index

    id                                   = L2E_FILTER_HYDROLOGYC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_hydrologyc     = index

    id                                   = L2E_FILTER_NUM_HYDROLOGYC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_num_hydrologyc = index

    id                                   = L2E_COLUMN_ZI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_column_zi             = index

    deallocate(em_stages)

  end subroutine EM_ATS_Populate_L2E_List

  !------------------------------------------------------------------------
  subroutine EM_ATS_Populate_E2L_List(this, e2l_list)
    !
    !
    ! !DESCRIPTION:
    ! Create a list of all variables to be returned by ATS from ELM
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_ats_type)                   :: this
    class(emi_data_list) , intent(inout) :: e2l_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index
    !
    !
    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_ATS_SOIL_HYDRO_STAGE

    id                              = E2L_STATE_H2OSOI_LIQ
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_h2osoi_liq = index

    id                              = E2L_STATE_H2OSOI_ICE
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_h2osoi_ice = index

    id                              = E2L_STATE_SOIL_MATRIC_POTENTIAL
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_smp        = index

    id                              = E2L_STATE_WTD
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_wtd        = index

    deallocate(em_stages)

end subroutine EM_ATS_Populate_E2L_List

  !------------------------------------------------------------------------
  subroutine EM_ATS_Init(this, l2e_init_list, e2l_init_list, iam, bounds_clump)

    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_ats_type)                   :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    integer              , intent(in)    :: iam
    type(bounds_type)    , intent(in)    :: bounds_clump

    ! (TODO)

  end subroutine EM_ATS_Init

  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  subroutine set_material_properties(this, l2e_init_list, bounds_clump)
    !
    ! !DESCRIPTION:
    !
    use elm_instMod               , only : soilstate_vars
    use elm_instMod               , only : soilhydrology_vars
    use elm_varpar                , only : nlevgrnd, nlevsoi
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_ats_type)               :: this
    class(emi_data_list), intent(in) :: l2e_init_list
    type(bounds_type)   , intent(in) :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer    :: elm_watsat(:,:)
    real(r8), pointer    :: elm_hksat(:,:)
    real(r8), pointer    :: elm_bsw(:,:)
    real(r8), pointer    :: elm_sucsat(:,:)
    real(r8), pointer    :: elm_eff_porosity(:,:)
    real(r8), pointer    :: elm_zwt(:)
    real(r8), pointer    :: ats_watsat(:,:)
    real(r8), pointer    :: ats_hksat(:,:)
    real(r8), pointer    :: ats_bsw(:,:)
    real(r8), pointer    :: ats_sucsat(:,:)
    real(r8), pointer    :: ats_eff_porosity(:,:)
    real(r8), pointer    :: ats_residual_sat(:,:)
    integer, pointer     :: ats_filter(:)
    !
    integer              :: c,g,fc,j,l            ! do loop indices
    integer  , pointer                   :: col_active(:)
    integer  , pointer                   :: col_type(:)
    integer  , pointer                   :: col_landunit(:)
    integer  , pointer                   :: lun_type(:)

    integer :: bounds_proc_begc_all, bounds_proc_endc_all
    integer :: bounds_proc_begc, bounds_proc_endc

    !-----------------------------------------------------------------------

    bounds_proc_begc_all = bounds_clump%begc_all
    bounds_proc_endc_all = bounds_clump%endc_all
    bounds_proc_begc     = bounds_clump%begc
    bounds_proc_endc     = bounds_clump%endc

    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_watsatc      , elm_watsat       )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_hksatc       , elm_hksat        )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_bswc         , elm_bsw          )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_sucsatc      , elm_sucsat       )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_effporosityc , elm_eff_porosity )
    call l2e_init_list%GetPointerToReal1D(this%index_l2e_init_state_wtd              , elm_zwt          )

    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_active              , col_active       )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_type                , col_type         )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_landunit_index      , col_landunit     )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_landunit_type           , lun_type         )

    ! Allocate memory
    allocate(ats_filter       (bounds_proc_begc_all:bounds_proc_endc_all           ))
    allocate(ats_watsat       (bounds_proc_begc_all:bounds_proc_endc_all, nlevgrnd ))
    allocate(ats_hksat        (bounds_proc_begc_all:bounds_proc_endc_all, nlevgrnd ))
    allocate(ats_bsw          (bounds_proc_begc_all:bounds_proc_endc_all, nlevgrnd ))
    allocate(ats_sucsat       (bounds_proc_begc_all:bounds_proc_endc_all, nlevgrnd ))
    allocate(ats_eff_porosity (bounds_proc_begc_all:bounds_proc_endc_all, nlevgrnd ))
    allocate(ats_residual_sat (bounds_proc_begc_all:bounds_proc_endc_all, nlevgrnd ))

    ! Initialize
    ats_filter       (:)   = 0
    ats_watsat       (:,:) = 0._r8
    ats_hksat        (:,:) = 0._r8
    ats_bsw          (:,:) = 0._r8
    ats_sucsat       (:,:) = 0._r8
    ats_residual_sat (:,:) = 0._r8

    ! Save data to initialize ATS
    do c = bounds_proc_begc, bounds_proc_endc
       l = col_landunit(c)

       if (col_active(c) == 1) then

          ats_filter    (c) = 1

          do j = 1 ,nlevgrnd
             ats_watsat(c,j)       = elm_watsat(c,j)
             ats_hksat(c,j)        = elm_hksat(c,j)
             ats_bsw(c,j)          = elm_bsw(c,j)
             ats_sucsat(c,j)       = elm_sucsat(c,j)
             ats_eff_porosity(c,j) = elm_eff_porosity(c,j)
          enddo

       endif
    enddo

    ! (TODO) pass data to ATS

    ! Free up memory
    deallocate(ats_filter       )
    deallocate(ats_watsat       )
    deallocate(ats_hksat        )
    deallocate(ats_bsw          )
    deallocate(ats_sucsat       )
    deallocate(ats_eff_porosity )
    deallocate(ats_residual_sat )

  end subroutine set_material_properties

  !------------------------------------------------------------------------
  subroutine set_initial_conditions(this, l2e_init_list, bounds_clump)
    !
    ! !DESCRIPTION:
    !
    use elm_instMod               , only : soilstate_vars
    use elm_instMod               , only : soilhydrology_vars
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_ats_type)              :: this
    class(emi_data_list), intent(in) :: l2e_init_list
    type(bounds_type)   , intent(in) :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    integer              :: c,g,fc,j,l            ! do loop indices
    !-----------------------------------------------------------------------

    !call l2e_init_list%GetPointerToReal1D(this%index_l2e_init_state_wtd         , elm_zwt      )
    !call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_zi            , elm_zi       )

  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine set_boundary_conditions(this, l2e_init_list, bounds_clump)
    !
    ! !DESCRIPTION:
    !
    use elm_instMod               , only : soilstate_vars
    use elm_instMod               , only : soilhydrology_vars
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_ats_type)               :: this
    class(emi_data_list), intent(in) :: l2e_init_list
    type(bounds_type)   , intent(in) :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    integer              :: c,g,fc,j,l            ! do loop indices
    !-----------------------------------------------------------------------


  end subroutine set_boundary_conditions

  !-----------------------------------------------------------------------


    !------------------------------------------------------------------------
   subroutine EM_ATS_OneStep(this, em_stage, dt, nstep, clump_rank, l2e_list, e2l_list, &
        bounds_clump)
    !
    ! !DESCRIPTION:
    ! The ATS driver subroutine
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_ats_type)                   :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump


  end subroutine EM_ATS_OneStep
  !-----------------------------------------------------------------------
  subroutine get_data_for_elm(this, l2e_init_list, e2l_init_list, bounds_clump)
    !
    !DESCRIPTION
    !  Saves
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_ats_type)                   :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    type(bounds_type)    , intent(in)    :: bounds_clump

    ! !LOCAL VARIABLES:
    integer                              :: p,c,fc,j,g  ! do loop indices

    real(r8)  , pointer                  :: l2e_col_zi(:,:)

    real(r8)  , pointer                  :: l2e_soilp(:,:)
    real(r8)  , pointer                  :: l2e_mflx_snowlyr_col(:)

    real(r8)  , pointer                  :: ats_soilp_col_1d(:)
    real(r8)  , pointer                  :: ats_mass_col_1d(:)
    real(r8)  , pointer                  :: ats_smpl_col_1d(:)
    real(r8)  , pointer                  :: e2l_h2osoi_liq(:,:)
    real(r8)  , pointer                  :: e2l_h2osoi_ice(:,:)
    real(r8)  , pointer                  :: e2l_smp_l(:,:)
    real(r8)  , pointer                  :: e2l_zwt(:)
    real(r8)  , pointer                  :: e2l_mflx_snowlyr_col(:)
    real(r8)  , pointer                  :: l2e_h2osoi_liq(:,:)
    real(r8)  , pointer                  :: l2e_h2osoi_ice(:,:)
    integer                              :: jwt
    integer                              :: idx
    integer                              :: soe_auxvar_id
    real(r8)                             :: z_up, z_dn
    integer :: bounds_proc_begc, bounds_proc_endc
    !-----------------------------------------------------------------------

    bounds_proc_begc     = bounds_clump%begc
    bounds_proc_endc     = bounds_clump%endc

    call l2e_init_list%GetPointerToReal1D(this%index_l2e_init_flux_mflx_snowlyr_col , l2e_mflx_snowlyr_col )

    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_state_soilp           , l2e_soilp            )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_zi                , l2e_col_zi           )

    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_h2osoi_liq            , l2e_h2osoi_liq       )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_h2osoi_ice            , l2e_h2osoi_ice       )

    call e2l_init_list%GetPointerToReal1D(this%index_e2l_init_state_wtd             , e2l_zwt              )
    call e2l_init_list%GetPointerToReal1D(this%index_e2l_init_flux_mflx_snowlyr_col , e2l_mflx_snowlyr_col )

    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_liq      , e2l_h2osoi_liq       )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_ice      , e2l_h2osoi_ice       )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_smp             , e2l_smp_l            )

   end subroutine get_data_for_elm

    !------------------------------------------------------------------------
   subroutine EM_ATS_Finalize(this, dt, nstep, clump_rank, l2e_list, e2l_list, &
        bounds_clump)
    !
    ! !DESCRIPTION:
    ! The ATS driver subroutine
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_ats_type)                   :: this
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump


  end subroutine EM_ATS_Finalize


    !------------------------------------------------------------------------

#endif

end module ExternalModelATSMod
