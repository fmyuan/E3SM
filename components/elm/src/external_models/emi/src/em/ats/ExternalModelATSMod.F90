module ExternalModelATSMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !

  ! ELM modules
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use spmdMod                      , only : masterproc, mpicom
  use decompMod                    , only : bounds_type
  use elm_varpar                   , only : nlevgrnd
  ! a few ats coupling options
  use elm_varctl                   , only : use_ats, use_ats_mesh
  use elm_varctl                   , only : ats_hmode, ats_thmode, ats_thcmode
  use elm_varctl                   , only : ats_chkout

  ! a few constants
  use elm_varcon                   , only : denh2o, denice, tfrz, grav
  use histFileMod                  , only : hist_nhtfrq


  ! EMI modules
  use EMI_DataMod                  , only : emi_data_list, emi_data
  use EMI_ColumnType_Constants
  use EMI_Filter_Constants

  use EMI_Atm2LndType_Constants
  use EMI_CanopyStateType_Constants
  use EMI_SoilStateType_Constants
  use EMI_SoilHydrologyType_Constants

  !use EMI_ColumnEnergyFluxType_Constants  ! need redo the list
  use EMI_EnergyFluxType_Constants
  use EMI_ColumnEnergyStateType_Constants
  use EMI_WaterFluxType_Constants
  use EMI_WaterStateType_Constants
  !
  use ExternalModelBaseType        , only : em_base_type
  use ExternalModelConstants

  use ExternalModelATS_createMod

  ! C-F interface
  use ELM_ATS_InterfaceMod


  !
  implicit none
  !
  ! EM data-type and procedures for elm-ats interface
  type, public, extends(em_base_type) :: em_ats_type
     ! ----------------------------------------------------------------------
     ! Indicies required during the initialization
     ! ----------------------------------------------------------------------
     integer :: index_l2e_init_col_active
     integer :: index_l2e_init_col_type
     integer :: index_l2e_init_col_filter
     integer :: index_l2e_init_col_filter_num
     integer :: index_l2e_init_col_zi
     integer :: index_l2e_init_col_dz
     integer :: index_l2e_init_col_z
     integer :: index_l2e_init_col_area

     integer :: index_l2e_init_col_patch_i_beg
     integer :: index_l2e_init_col_patch_i_end
     integer :: index_l2e_init_col_num_patch
     integer :: index_l2e_init_col_pft_type

     integer :: index_l2e_init_state_forc_pbot
     integer :: index_l2e_init_state_h2osoi_liq
     integer :: index_l2e_init_state_h2osoi_ice
     integer :: index_l2e_init_state_soilp
     integer :: index_l2e_init_state_smp
     integer :: index_l2e_init_state_frac_h2osfc
     integer :: index_l2e_init_state_h2osfc
     integer :: index_l2e_init_state_zwt

     integer :: index_l2e_init_state_ts_soil
     integer :: index_l2e_init_state_ts_snow
     integer :: index_l2e_init_state_ts_h2osfc

     integer :: index_l2e_init_parameter_watsatc
     integer :: index_l2e_init_parameter_hksatc
     integer :: index_l2e_init_parameter_bswc
     integer :: index_l2e_init_parameter_sucsatc
     integer :: index_l2e_init_parameter_effporosityc

     ! ----------------------------------------------------------------------
     ! Indicies required during timestepping
     ! ----------------------------------------------------------------------
     integer :: index_l2e_filter
     integer :: index_l2e_filter_num
     integer :: index_l2e_column_zi
     integer :: index_l2e_column_dz

     integer :: index_l2e_state_forc_pbot
     integer :: index_l2e_state_h2osoi_liq
     integer :: index_l2e_state_h2osoi_ice
     integer :: index_l2e_state_soilp
     integer :: index_l2e_state_smp
     integer :: index_l2e_state_zwt

     integer :: index_l2e_parameter_effporosityc
     integer :: index_l2e_rootfrac_col
     integer :: index_l2e_rootfrac_patch

     integer :: index_e2l_state_h2osoi_liq
     integer :: index_e2l_state_h2osoi_ice
     integer :: index_e2l_state_soilp
     integer :: index_e2l_state_smp
     integer :: index_e2l_state_zwt
     integer :: index_e2l_state_h2osfc

     integer :: index_l2e_flux_gross_infil
     integer :: index_l2e_flux_gross_evap
     integer :: index_l2e_flux_pot_tran_pft
     integer :: index_l2e_flux_pot_tran_vr
     integer :: index_l2e_flux_drain_vr
     integer :: index_l2e_flux_surf

     integer :: index_e2l_flux_infil
     integer :: index_e2l_flux_evap
     integer :: index_e2l_flux_tran_veg
     integer :: index_e2l_flux_rootsoi_col
     integer :: index_e2l_flux_rootsoi_patch

     integer :: index_l2e_state_ts_soil
     integer :: index_l2e_state_ts_snow
     integer :: index_l2e_state_ts_h2osfc

     integer :: index_e2l_state_tsoil
     integer :: index_l2e_flux_hs_soil

     ! save col and patch filters
     integer                      :: filter_col_num
     integer   , pointer          :: filter_col(:), filter_pft(:)

     ! elm nstep at start/restart
     integer :: elmnstep0

     type(elm_ats_interface_type) :: ats_interface

   contains

     procedure, public :: Populate_L2E_Init_List  => EM_ATS_Populate_L2E_Init_List
     procedure, public :: Populate_E2L_Init_List  => EM_ATS_Populate_E2L_Init_List
     procedure, public :: Populate_L2E_List       => EM_ATS_Populate_L2E_List
     procedure, public :: Populate_E2L_List       => EM_ATS_Populate_E2L_List
     procedure, public :: Init                    => EM_ATS_Init
     procedure, public :: Solve                   => EM_ATS_OneStep
     procedure, public :: Finalize                => EM_ATS_Finalize
  end type em_ats_type

  !---------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine EM_ATS_Populate_L2E_Init_List(this, l2e_init_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by ATS from ELM during initialization
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
    id                                         = L2E_COLUMN_ACTIVE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_active             = index

    id                                         = L2E_COLUMN_TYPE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_type               = index

    id                                         = L2E_FILTER_SOILC              ! in future, this may be flexible
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_filter             = index

    id                                         = L2E_FILTER_NUM_SOILC          ! in future, this may be flexible
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_filter_num         = index

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

    id                                         = L2E_COLUMN_PATCH_INDEX_BEGIN
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_patch_i_beg        = index

    id                                         = L2E_COLUMN_PATCH_INDEX_END
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_patch_i_end        = index

    id                                         = L2E_COLUMN_NUM_PATCH
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_num_patch          = index

    id                                         = L2E_COLUMN_PFT_TYPE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_pft_type           = index
    !-------------

    id                                        = L2E_STATE_VSFM_PROGNOSTIC_SOILP
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_soilp           = index

    id                                        = L2E_STATE_SOIL_MATRIC_POTENTIAL_NLEVSOI
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_smp             = index

    id                                        = L2E_STATE_H2OSOI_LIQ_NLEVGRND
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_h2osoi_liq      = index

    id                                        = L2E_STATE_H2OSOI_ICE_NLEVGRND
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_h2osoi_ice      = index

    !id                                        = L2E_STATE_FORC_PBOT_DOWNSCALED
    !call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    !this%index_l2e_init_state_forc_pbot       = index

    id                                        = L2E_STATE_TSOIL_NLEVGRND_COL
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_ts_soil         = index

    id                                        = L2E_STATE_TSNOW_COL
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_ts_snow         = index

    id                                        = L2E_STATE_TH2OSFC_COL
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_ts_h2osfc       = index

    id                                        = L2E_STATE_FRAC_H2OSFC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_frac_h2osfc     = index

    id                                        = L2E_STATE_H2OSFC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_h2osfc          = index

    id                                        = L2E_STATE_WTD
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_zwt             = index

    !--------
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

#ifdef ATS_READY
  ! DON'T INITIALIZE ELM by ATS's states
    id                                        = E2L_STATE_H2OSOI_LIQ
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_liq      = index

    id                                        = E2L_STATE_H2OSOI_ICE
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_ice      = index

    id                                        = E2L_STATE_SOIL_MATRIC_POTENTIAL
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_smp             = index

    id                                        = E2L_STATE_VSFM_PROGNOSTIC_SOILP
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_soilp           = index

    id                                        = E2L_STATE_WTD
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_zwt             = index

    if (em_stages(1) == EM_ATS_SOIL_THYDRO_STAGE .or. &
        em_stages(1) == EM_ATS_SOIL_THBGC_STAGE) then
      id                                      = E2L_STATE_TSOIL_NLEVGRND_COL
      call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
      this%index_e2l_init_state_tsoil         = index
    end if
#endif

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

    id                                   = L2E_FILTER_SOILC              ! in future, this may be flexible
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter                = index

    id                                   = L2E_FILTER_NUM_SOILC          ! in future, this may be flexible
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_num            = index

    id                                   = L2E_COLUMN_ZI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_column_zi             = index

    id                                   = L2E_COLUMN_DZ
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_column_dz             = index

    !-------------
    id                                   = L2E_STATE_VSFM_PROGNOSTIC_SOILP
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_soilp           = index

    id                                   = L2E_STATE_SOIL_MATRIC_POTENTIAL_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_smp             = index

    id                                   = L2E_STATE_H2OSOI_LIQ_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_liq      = index

    id                                   = L2E_STATE_H2OSOI_ICE_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_ice      = index

    id                                   = L2E_FLUX_GROSS_INFL_SOIL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_gross_infil      = index

    id                                   = L2E_FLUX_GROSS_EVAP_SOIL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_gross_evap       = index

    id                                   = L2E_FLUX_TRAN_VEG  ! may not need, if 'L2E_FLUX_ROOTSOI' (root-layered) used
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_pot_tran_pft     = index

    id                                   = L2E_FLUX_ROOTSOI   ! may not need, if 'L2E_FLUX_TRAN_VEG' (vegcanopy-integrated) used
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_pot_tran_vr      = index

    id                                   = L2E_FLUX_DRAIN_VR
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_drain_vr         = index

    id                                   = L2E_FLUX_SURF
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_surf             = index

    id                                   = L2E_PARAMETER_EFFPOROSITYC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_parameter_effporosityc = index

    id                                   = L2E_PARAMETER_ROOTFR_COL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_rootfrac_col          = index

    id                                   = L2E_PARAMETER_ROOTFR_PATCH
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_rootfrac_patch        = index

    !-----------------
    if (em_stages(1) == EM_ATS_SOIL_THYDRO_STAGE .or. &
        em_stages(1) == EM_ATS_SOIL_THBGC_STAGE) then

        id                                   = L2E_STATE_TSOIL_NLEVGRND_COL
        call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
        this%index_l2e_state_ts_soil         = index

        id                                   = L2E_STATE_TSNOW_COL
        call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
        this%index_l2e_state_ts_snow         = index

        id                                   = L2E_STATE_TH2OSFC_COL
        call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
        this%index_l2e_state_ts_h2osfc       = index

        id                                   = L2E_FLUX_SOIL_HEAT_FLUX
        call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
        this%index_l2e_flux_hs_soil          = index

    endif

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

    id                              = E2L_STATE_SOIL_MATRIC_POTENTIAL_COL
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_smp        = index

    id                              = E2L_STATE_VSFM_PROGNOSTIC_SOILP
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_soilp      = index

    id                              = E2L_STATE_H2OSFC
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_h2osfc     = index

    id                              = E2L_STATE_WTD
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_zwt        = index

    id                              = E2L_FLUX_INFL_SOIL
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_infil = index
    
    id                              = E2L_FLUX_EVAP_SOIL
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_evap  = index

    id                              = E2L_FLUX_TRAN_VEG
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_tran_veg    = index

    id                              = E2L_FLUX_ROOTSOI
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_rootsoi_col = index

    id                              = E2L_FLUX_ROOTSOI_FRAC
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_rootsoi_patch= index
    !-----------------
    if (em_stages(1) == EM_ATS_SOIL_THYDRO_STAGE .or. &
        em_stages(1) == EM_ATS_SOIL_THBGC_STAGE) then

      id                              = E2L_STATE_TSOIL_NLEVGRND_COL
      call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
      this%index_e2l_state_tsoil      = index
    end if

    deallocate(em_stages)

  end subroutine EM_ATS_Populate_E2L_List

  !------------------------------------------------------------------------
  subroutine EM_ATS_Init(this, l2e_init_list, e2l_init_list, iam, bounds_clump)

    use elm_time_manager, only: get_nstep

    implicit none
    class(em_ats_type)                   :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    integer              , intent(in)    :: iam
    type(bounds_type)    , intent(in)    :: bounds_clump
    integer :: filternum, fc, c, p
    integer, pointer :: filtercol(:), pft_type(:), numpatch(:), pft_beg(:), pft_end(:)

    !-----------------------------------------------------------------------

    this%elmnstep0 = get_nstep()
    !
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_num_patch, numpatch )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_patch_i_beg, pft_beg )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_patch_i_end, pft_end )
    call l2e_init_list%GetIntValue(this%index_l2e_init_col_filter_num, filternum )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_filter, filtercol )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_pft_type, pft_type )
    
    this%filter_col_num = filternum
    
    allocate(this%filter_col(filternum))
    this%filter_col(1:filternum) = filtercol(1:filternum)

    ! 1 pft per column, so use column indexing
    allocate(this%filter_pft(filternum))

    do fc=1, filternum
       c = this%filter_col(fc)
       this%filter_pft(fc) = pft_beg(c) + pft_type(c)
    end do

    if (use_ats) then

      ! create an ATS driver object (via mpicom)
      ! if not use_ats_mesh (or, elm will pass its mesh info to ats)
      if (.not.use_ats_mesh) then
        call EM_ATS_create(this%ats_interface)

        ! pass mesh data to ATS prior to ATS setup (as long as ATS driver object ready)
        call set_mesh(this, l2e_init_list, bounds_clump)
      end if

      ! ats setup
      call this%ats_interface%setup()

      ! if use_ats_mesh, ATS material properties will be passed into ELM as well (TODO)
      if (.not.use_ats_mesh) then
        ! pass material properties to ATS
        call set_material_properties(this, l2e_init_list, bounds_clump)
      end if

      print *, 'EM_ATS_Init @nstep: ', this%elmnstep0

    end if

  end subroutine EM_ATS_Init

  !------------------------------------------------------------------------


  subroutine set_mesh(this, l2e_init_list, bounds_clump)
    implicit none
    class(em_ats_type)                   :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    type(bounds_type)    , intent(in)    :: bounds_clump
    integer                              :: c, fc, j

    real(r8)  , pointer                  :: surf_xi(:)
    real(r8)  , pointer                  :: surf_yi(:)
    real(r8)  , pointer                  :: surf_zi(:,:)
    real(r8)  , pointer                  :: surf_area(:,:)

    real(r8)  , pointer                  :: col_zi(:,:)
    real(r8)  , pointer                  :: col_dz(:,:)
    real(r8)  , pointer                  :: col_z (:,:)

    integer                              :: nz_nodes
    real(r8)  , pointer                  :: col_nodes (:)

    !-----------------------------------------------------------------------

    allocate(surf_xi(0:this%filter_col_num))           ! 0 index starting, so ready for c++ ats data
    allocate(surf_yi(0:1))                             ! assuming columns arranged along X-axis, so gridY nodes is 2
    allocate(surf_zi(0:this%filter_col_num,0:1))       ! surf-grid nodes in 2-D
    allocate(surf_area(0:this%filter_col_num-1,0:0))   ! surf-grid cells in 2-D

    ! hard-wired now as unit distance (m) (TO-FIX)
    surf_xi (0) = 0.0d0
    do c=1, this%filter_col_num
        surf_xi(c) = surf_xi(c-1)+1.0d0
    end do
    surf_yi (0) = 0.0d0; surf_yi(1) = 1.0d0

    surf_zi (:,:) = 0.0d0
    surf_area(:,:)= 1.0d0

    !call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_area           , surf_areas)
    !surf_area(:,0) = surf_areas

    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_zi             , col_zi)   ! layer top/bottom-node (interface) depth (0:nlevgrnd)
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_dz             , col_dz)   ! layer thickness (1:nlevgrnd)
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_z              , col_z )   ! layer centroid depth (1:nlevgrnd)

    nz_nodes = size(col_zi(1,:))
    surf_zi(:,:) = 0.0         ! (TO-FIX) temporarily set to 0

    allocate(col_nodes(0:nz_nodes-1))
    c = 1                                       ! here assuming soil column node coords are same for all ELM columns
    j = 0
    col_nodes(:) = surf_zi(c,0) - col_zi(c,:)   ! elevation (m) for all vertical nodes

    !call this%ats_interface%setmesh(surf_xi, surf_yi, surf_zi, col_nodes)

    deallocate(surf_xi)
    deallocate(surf_yi)
    deallocate(surf_zi)
    deallocate(surf_area)
    deallocate(col_nodes)

  end subroutine set_mesh


  !------------------------------------------------------------------------
  subroutine set_material_properties(this, l2e_init_list, bounds_clump)

    use elm_varpar                , only : nlevgrnd, nlevsoi
    implicit none
    class(em_ats_type)               :: this
    class(emi_data_list), intent(in) :: l2e_init_list
    type(bounds_type)   , intent(in) :: bounds_clump
    real(r8), pointer    :: elm_watsat(:,:)
    real(r8), pointer    :: elm_hksat(:,:)
    real(r8), pointer    :: elm_bsw(:,:)
    real(r8), pointer    :: elm_sucsat(:,:)
    real(r8), pointer    :: ats_watsat(:,:)
    real(r8), pointer    :: ats_hksat(:,:)
    real(r8), pointer    :: ats_bsw(:,:)
    real(r8), pointer    :: ats_sucsat(:,:)
    real(r8), pointer    :: ats_residual_sat(:,:)
    integer              :: fc, c, j
    integer              :: nz

    !-----------------------------------------------------------------------
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_watsatc      , elm_watsat       )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_hksatc       , elm_hksat        )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_bswc         , elm_bsw          )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_sucsatc      , elm_sucsat       )

    nz = size(elm_watsat(1,:))
    allocate(ats_watsat       (this%filter_col_num, nz ))
    allocate(ats_hksat        (this%filter_col_num, nz ))
    allocate(ats_bsw          (this%filter_col_num, nz ))
    allocate(ats_sucsat       (this%filter_col_num, nz ))
    allocate(ats_residual_sat (this%filter_col_num, nz ))

    ats_watsat       (:,:) = 0._r8
    ats_hksat        (:,:) = 1.0e-6_r8 ! this is the value that will be assigned to the 'bedrock' layers in ATS 
    ats_bsw          (:,:) = 0._r8
    ats_sucsat       (:,:) = 0._r8
    ats_residual_sat (:,:) = 0._r8    ! alway zero from ELM

    do fc = 1, this%filter_col_num
      c = this%filter_col(fc)
      do j = 1, nz
         ats_watsat(fc,j)       = elm_watsat(c,j)
         ats_hksat(fc,j)        = ats_hksat(fc,j) + elm_hksat(c,j) * 0.001_r8
         ats_bsw(fc,j)          = elm_bsw(c,j)
         ats_sucsat(fc,j)       = elm_sucsat(c,j)
      end do
    end do

    ! pass data to ATS
    call this%ats_interface%setsoilveg(ats_watsat, ats_hksat, &
         ats_bsw, ats_sucsat, ats_residual_sat)

    deallocate(ats_watsat)
    deallocate(ats_hksat)
    deallocate(ats_bsw)
    deallocate(ats_sucsat)
    deallocate(ats_residual_sat)
  end subroutine set_material_properties

  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
   subroutine EM_ATS_OneStep(this, em_stage, dt, nstep, clump_rank, l2e_list, e2l_list, &
        bounds_clump)

    implicit none
    class(em_ats_type)                   :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump

    ! local variables
    logical :: elm_histout
    logical :: elm_restout
    !-----------------------

    ! initialize all fields in ATS
    if (nstep==this%elmnstep0) then
      call set_initial_conditions(this, nstep, dt, l2e_list, bounds_clump)
    end if

    ! set rates of infiltration and potential evaporation/transpiration
    call set_bc_ss(this, l2e_list, bounds_clump)

    ! set root fraction distribution 
    call set_dyn_properties(this, l2e_list, bounds_clump)

    elm_histout = .false.
    elm_restout = .false.
    if (hist_nhtfrq(1)>0) then
      if (mod(nstep,hist_nhtfrq(1)) == 0) elm_histout = .true.
    end if
    if (ats_chkout) elm_restout = .true.

    ! hardwired for now
    elm_histout = .true.

    ! Advance ATS to time t+dt 
    call this%ats_interface%onestep(dt, elm_histout, elm_restout)

    ! get new water state and fluxes from ATS
    call get_data_for_elm(this, e2l_list, bounds_clump)

  end subroutine EM_ATS_OneStep

  !------------------------------------------------------------------------
  subroutine set_initial_conditions(this, nstep, dt, l2e_list, bounds_clump)

    implicit none
    class(em_ats_type)               :: this
    class(emi_data_list), intent(in) :: l2e_list
    real(r8)            , intent(in) :: dt             ! unit: seconds
    integer             , intent(in) :: nstep
    type(bounds_type)   , intent(in) :: bounds_clump
    integer           :: fc, c, j
    integer           :: nz
    real(r8)          :: starting_time
    real(r8), pointer :: soilp(:,:), soilp_ats(:,:)
    real(r8), pointer :: wtd(:), wtd_ats(:)
    real(r8), pointer :: h2osoi_liq(:,:), h2oliq_ats(:,:)

    call l2e_list%GetPointerToReal2D(this%index_l2e_state_soilp            , soilp )
    nz = size(soilp(1,:))
    allocate(soilp_ats(this%filter_col_num, nz))

    allocate(h2oliq_ats(this%filter_col_num, nz))
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_liq       , h2osoi_liq )

    do fc = 1, this%filter_col_num
      c = this%filter_col(fc)
      do j = 1, nz
        h2oliq_ats(fc,j) = h2osoi_liq(c,j)
      end do
    end do

    starting_time = nstep*dt
    call this%ats_interface%init(starting_time, h2oliq_ats, soilp_ats)
    deallocate(h2oliq_ats)
    deallocate(soilp_ats)
  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine set_dyn_properties(this, l2e_list, bounds_clump)

    use elm_varpar, only : nlevgrnd, nlevsoi
    implicit none
    class(em_ats_type)               :: this
    class(emi_data_list), intent(in) :: l2e_list
    type(bounds_type)   , intent(in) :: bounds_clump
    integer           :: fc, c, j
    integer           :: nc, nz
    real(r8), pointer :: soil_effporo(:,:), soil_effporo_ats(:,:)
    !real(r8), pointer :: pft_rootfrac(:,:), pft_rootfrac_ats(:,:)  ! if 'col_rootfrac' NOT used
    real(r8), pointer :: col_rootfrac(:,:), col_rootfrac_ats(:,:)
    integer,  pointer :: filter_patch(:)
    integer, dimension (5) :: rootfrac_idx ! quick and dirty hardwired pft index for 2D run
    !-----------------------------------------------------------------------
    
    nc = this%filter_col_num
    nz = nlevgrnd

    allocate(soil_effporo_ats(nc, nz))
    call l2e_list%GetPointerToReal2D(this%index_l2e_parameter_effporosityc, soil_effporo)

    allocate(col_rootfrac_ats(nc, nz))
    call l2e_list%GetPointerToReal2D(this%index_l2e_rootfrac_col, col_rootfrac)

    do fc = 1, nc
      c = this%filter_col(fc)
      do j = 1, nz
         col_rootfrac_ats(fc,j) = col_rootfrac(c,j)
      end do
    end do

    call this%ats_interface%setsoilveg_dyn(soil_effporo_ats, col_rootfrac_ats)

    deallocate(col_rootfrac_ats)
    deallocate(soil_effporo_ats)

  end subroutine set_dyn_properties

  !------------------------------------------------------------------------
  subroutine set_bc_ss(this, l2e_list, bounds_clump)

    implicit none
    class(em_ats_type)               :: this
    class(emi_data_list), intent(in) :: l2e_list
    type(bounds_type)   , intent(in) :: bounds_clump
    integer              :: c, fc,  j,  p            ! do loop indices
    integer              :: nc, nz
    real(r8), pointer    :: col_dz(:,:)
    real(r8), pointer    :: soilevap_flux(:),    soilinfl_flux(:),       soilevap_flux_ats(:),    soilinfl_flux_ats(:)
    real(r8), pointer    :: soilbot_flux(:),                             soilbot_flux_ats(:)
    real(r8), pointer    :: tran_veg_flux(:),                            tran_veg_flux_ats(:)      ! summed veg. transpiration for column
    !-----------------------------------------------------------------------
    nc = this%filter_col_num

    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_gross_infil , soilinfl_flux )   ! mmH2O/s, + to soil
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_gross_evap  , soilevap_flux )   ! mmH2O/s, + to atm
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_pot_tran_pft, tran_veg_flux )   ! mmH2O/s, + to atm
    allocate(soilinfl_flux_ats(nc))
    allocate(soilevap_flux_ats(nc))
    allocate(tran_veg_flux_ats(nc))

    do fc = 1, nc
      c = this%filter_col(fc)
      soilinfl_flux_ats(fc) = soilinfl_flux(c)
      soilevap_flux_ats(fc) = -soilevap_flux(c)
      tran_veg_flux_ats(fc) = -tran_veg_flux(c)
    end do

    call this%ats_interface%setss_hydro(soilinfl_flux_ats, soilevap_flux_ats, tran_veg_flux_ats)

    deallocate(soilinfl_flux_ats)
    deallocate(soilevap_flux_ats)
    deallocate(tran_veg_flux_ats)

  end subroutine set_bc_ss

  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------

  subroutine get_data_for_elm(this, e2l_list, bounds_clump)

    implicit none
    class(em_ats_type)                   :: this
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump
    integer              :: c, fc,  j,  p            ! do loop indices
    integer              :: nc, nz
    real(r8)  , pointer  :: e2l_h2osfc(:),         h2osfc_ats(:)
    real(r8)  , pointer  :: e2l_zwt(:),            zwt_ats(:)
    real(r8)  , pointer  :: e2l_h2osoi_liq(:,:),   h2oliq_ats(:,:)
    real(r8)  , pointer  :: e2l_h2osoi_ice(:,:),   h2oice_ats(:,:)
    real(r8)  , pointer  :: e2l_smp_l(:,:), e2l_soilp(:,:), soilpsi_ats(:,:)
    real(r8)  , pointer  :: e2l_soilinfl_flux(:),  soilinfl_flux_ats(:)
    real(r8)  , pointer  :: e2l_evap_flux(:),      evap_flux_ats(:)
    real(r8)  , pointer  :: e2l_root_flux(:,:),    root_flux_ats(:,:)
    real(r8)  , pointer  :: e2l_tran_flux(:),      tran_flux_ats(:)
    real(r8)  , pointer  :: e2l_net_sub_flux(:),   net_sub_flux_ats(:)
    real(r8)  , pointer  :: e2l_net_srf_flux(:),   net_srf_flux_ats(:)
    real(r8)  , pointer  :: e2l_rootsoi_frac(:,:)

    !-----------------------------------------------------------------------
    !
    call e2l_list%GetPointerToReal1D(this%index_e2l_state_h2osfc     , e2l_h2osfc     )
    call e2l_list%GetPointerToReal1D(this%index_e2l_state_zwt        , e2l_zwt        )
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_h2osoi_liq , e2l_h2osoi_liq )
    !call e2l_list%GetPointerToReal2D(this%index_e2l_state_h2osoi_ice , e2l_h2osoi_ice )  ! TODO
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_smp        , e2l_smp_l      )
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_soilp      , e2l_soilp      )

    ! when coupling with ATS, following 2 vars are actual (not gross) returned from ATS
    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_infil , e2l_soilinfl_flux )
    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_evap  , e2l_evap_flux     )

    call e2l_list%GetPointerToReal2D(this%index_e2l_flux_rootsoi_col, e2l_root_flux)  ! if by root-layered transpiration of each column outputs from ATS
    !call e2l_list%GetPointerToReal2D(this%index_e2l_flux_rootsoi_patch, e2l_rootsoi_frac  )  ! if by root-layered transpiration of each PFT outputs from ATS
    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_tran_veg    , e2l_tran_flux     )

    allocate(h2osfc_ats(this%filter_col_num))
    allocate(zwt_ats(this%filter_col_num))
    allocate(net_sub_flux_ats(this%filter_col_num))
    allocate(net_srf_flux_ats(this%filter_col_num))
    allocate(soilinfl_flux_ats(this%filter_col_num))
    allocate(evap_flux_ats(this%filter_col_num))
    allocate(tran_flux_ats(this%filter_col_num))

    nz = size(e2l_h2osoi_liq(1,:))
    allocate(h2oliq_ats(this%filter_col_num, nz))
    allocate(h2oice_ats(this%filter_col_num, nz))
    allocate(soilpsi_ats(this%filter_col_num, nz))
    allocate(root_flux_ats(this%filter_col_num, nz))

    call this%ats_interface%getdata_hydro(h2osfc_ats, zwt_ats, h2oliq_ats, h2oice_ats, soilpsi_ats, &
                                          soilinfl_flux_ats, evap_flux_ats, tran_flux_ats, &
                                          root_flux_ats, net_sub_flux_ats, net_srf_flux_ats)

    ! Conversions occur in ATS
    ! h2osfc_ats[m]      ->  e2l_h2osfc[mm]
    ! zwt_ats[m]         ->  e2l_zwt[m]
    ! h2oliq_ats[-]      ->  e2l_h2osoi_liq[kg/m2]
    do fc = 1, this%filter_col_num
       c = this%filter_col(fc)
       !p = this%filter_pft(fc)
       e2l_h2osfc(c)          = h2osfc_ats(fc)
       e2l_zwt(c)             = zwt_ats(fc)
       e2l_soilinfl_flux(c)   = soilinfl_flux_ats(fc)
       e2l_evap_flux(c)       = evap_flux_ats(fc)
       e2l_tran_flux(c)       = tran_flux_ats(fc)
       do j = 1, nz
          e2l_h2osoi_liq(c,j) = h2oliq_ats(fc,j)
          !e2l_h2osoi_ice(c,j) = h2oice_ats(fc,j)                   ! TODO when thermal coupling is ready
          e2l_smp_l(c,j)      = soilpsi_ats(fc,j)/(-9.80665_r8)    ! Pa ---> -mmH2O
          e2l_soilp(c,j)      = -soilpsi_ats(fc,j)*1.0e-6_r8 + 0.1013250_r8    ! Pa ---> MPa
          e2l_root_flux(c,j)  = root_flux_ats(fc,j)
       end do
    end do

    ! STILL NEED net_sub_flux, net_srf_flux

    deallocate(h2osfc_ats)
    deallocate(zwt_ats)
    deallocate(h2oliq_ats)
    deallocate(h2oice_ats)
    deallocate(soilpsi_ats)
    deallocate(soilinfl_flux_ats)
    deallocate(evap_flux_ats)
    deallocate(tran_flux_ats)
    deallocate(root_flux_ats)
    deallocate(net_sub_flux_ats)
    deallocate(net_srf_flux_ats)

   end subroutine get_data_for_elm

  !------------------------------------------------------------------------
   subroutine EM_ATS_Finalize(this, dt, nstep, clump_rank, l2e_list, e2l_list, &
        bounds_clump)

    implicit none
    class(em_ats_type)                   :: this
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump

    ! (TODO)
    if (associated(this%filter_col)) deallocate(this%filter_col)
    if (associated(this%filter_pft)) deallocate(this%filter_pft)

  end subroutine EM_ATS_Finalize

 !------------------------------------------------------------------------

end module ExternalModelATSMod
