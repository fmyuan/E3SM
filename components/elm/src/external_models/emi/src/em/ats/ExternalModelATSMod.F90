module ExternalModelATSMod

#ifdef USE_ATS_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !

  ! ELM modules
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use spmdMod                      , only : masterproc, mpicom
  use decompMod                    , only : bounds_type
  use elm_varpar                   , only : nlevgrnd
  ! a few ats coupling options
  use elm_varctl                   , only : use_ats, ats_hmode, ats_thmode, ats_thcmode, ats_gmode
  use elm_varcon                   , only : denh2o, denice, tfrz, grav

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

  use ExternalModelATS_readnlMod   , only : ats_inputdir, ats_inputfile

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

     integer :: index_l2e_init_state_forc_pbot
     integer :: index_l2e_init_state_h2osoi_liq
     integer :: index_l2e_init_state_h2osoi_ice
     integer :: index_l2e_init_state_soilp
     integer :: index_l2e_init_state_smp

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

     integer :: index_e2l_state_h2osoi_liq
     integer :: index_e2l_state_h2osoi_ice
     integer :: index_e2l_state_soilp
     integer :: index_e2l_state_smp

     integer :: index_l2e_flux_gross_infil
     integer :: index_l2e_flux_gross_evap
     integer :: index_l2e_flux_tran_vr
     integer :: index_l2e_flux_drain_vr
     integer :: index_l2e_flux_surf

     integer :: index_e2l_flux_qrecharge

     integer :: index_l2e_state_ts_soil
     integer :: index_l2e_state_ts_snow
     integer :: index_l2e_state_ts_h2osfc

     integer :: index_e2l_state_tsoil

     integer :: index_l2e_flux_hs_soil
     integer :: index_l2e_flux_hs_top_snow
     integer :: index_l2e_flux_hs_h2osfc
     integer :: index_l2e_flux_dhsdT
     integer :: index_l2e_flux_sabg_lyr

     ! IDs to indentify the conditions for ATS
     integer :: ats_cond_id_for_infil
     integer :: ats_cond_id_for_et
     integer :: ats_cond_id_for_dew
     integer :: ats_cond_id_for_drainage
     integer :: ats_cond_id_for_snow
     integer :: ats_cond_id_for_sublimation
     integer :: ats_cond_id_for_lateral_flux

     !
     integer                      :: filter_col_num
     integer   , pointer          :: filter_col(:)


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
    id                                        = E2L_STATE_H2OSOI_LIQ
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_liq      = index

    id                                        = E2L_STATE_H2OSOI_ICE
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_ice      = index

    id                                        = E2L_STATE_SOIL_MATRIC_POTENTIAL_NLEVSOI
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_smp             = index

    id                                        = E2L_STATE_WTD
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_wtd             = index

    id                                        = E2L_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_flux_mflx_snowlyr_col = index
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

    id                                   = L2E_FLUX_ROOTSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_tran_vr          = index

    id                                   = L2E_FLUX_GROSS_EVAP_SOIL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_gross_evap       = index

    id                                   = L2E_FLUX_DRAIN_VR
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_drain_vr         = index

    id                                   = L2E_FLUX_SURF
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_surf             = index

    !-----------------
    if (em_stages(1) == EM_ATS_SOIL_THYDRO_STAGE .or. &
        em_stages(1) == EM_ATS_SOIL_THBGC_STAGE) then

        id                                   = L2E_STATE_FORC_PBOT_DOWNSCALED
        call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
        this%index_l2e_state_forc_pbot       = index

        id                                   = L2E_STATE_TSOIL_NLEVGRND_COL
        call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
        this%index_l2e_state_ts_soil         = index

        id                                   = L2E_STATE_TSNOW_COL
        call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
        this%index_l2e_state_ts_snow         = index

        id                                   = L2E_STATE_TH2OSFC_COL
        call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
        this%index_l2e_state_ts_h2osfc       = index

        id                                   = L2E_FLUX_ABSORBED_SOLAR_RADIATION
        call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
        this%index_l2e_flux_sabg_lyr         = index

        id                                   = L2E_FLUX_SOIL_HEAT_FLUX
        call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
        this%index_l2e_flux_hs_soil          = index

        id                                   = L2E_FLUX_SNOW_HEAT_FLUX
        call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
        this%index_l2e_flux_hs_top_snow      = index

        id                                   = L2E_FLUX_H2OSFC_HEAT_FLUX
        call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
        this%index_l2e_flux_hs_h2osfc                  = index

        id                                   = L2E_FLUX_DERIVATIVE_OF_HEAT_FLUX
        call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
        this%index_l2e_flux_dhsdT            = index

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
#ifdef ATS_READY
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
#endif
    deallocate(em_stages)

  end subroutine EM_ATS_Populate_E2L_List

  !------------------------------------------------------------------------
  subroutine EM_ATS_Init(this, l2e_init_list, e2l_init_list, iam, bounds_clump)

    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use timeinfoMod

    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_ats_type)                   :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    integer              , intent(in)    :: iam
    type(bounds_type)    , intent(in)    :: bounds_clump

    ! local
    integer :: filternum
    integer, pointer :: filtercol(:)
    real(r8) :: init_timesecond = 0.0 ! (TODO - shall be set by really initial or restarted time)

    !-----------------------------------------------------------------------

    !
    call l2e_init_list%GetIntValue(this%index_l2e_init_col_filter_num   , filternum )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_filter , filtercol )
    this%filter_col_num = filternum
    allocate(this%filter_col(filternum))
    this%filter_col(1:filternum) = filtercol(1:filternum)

    !
    ! create an ATS driver object
    this%ats_interface = ats_drv()


    if (use_ats) then
      ! pass mesh data to ATS prior to ATS setup (as long as ATS driver object ready)
      call set_mesh(this, l2e_init_list, bounds_clump)

      ! pass material properties to ATS prior to ATS setup
      call set_material_properties(this, l2e_init_list, bounds_clump)

      ! 'mpicom' is communicator group id for land component
      print *, 'EM_ATS_Init: ats inputs - ', trim(ats_inputdir), ' ', trim(ats_inputfile)
      print *, 'communicator id: ', mpicom
      call this%ats_interface%setup(ats_inputdir, ats_inputfile, mpicom)
      !
      ! reset ICs for ATS
      ! call set_initial_conditions(this, l2e_init_list, bounds_clump)  ! soilp not yet available at this stage ?

      ! after setting up mesh, materials
      init_timesecond = nstep_mod * dtime_mod
      call this%ats_interface%init()

    end if

  end subroutine EM_ATS_Init

  !------------------------------------------------------------------------


  subroutine set_mesh(this, l2e_init_list, bounds_clump)
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
    type(bounds_type)    , intent(in)    :: bounds_clump

    ! !LOCAL VARIABLES:
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
    !allocate(surf_xi (bounds_clump%bounds_proc_begc_all:bounds_clump%bounds_proc_endc_all))
    !allocate(surf_yi (bounds_clump%bounds_proc_begc_all:bounds_clump%bounds_proc_endc_all))
    !allocate(surf_zi (bounds_clump%bounds_proc_begc_all:bounds_clump%bounds_proc_endc_all))

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

    call this%ats_interface%setmesh(surf_xi, surf_yi, surf_zi, col_nodes)

    deallocate(surf_xi)
    deallocate(surf_yi)
    deallocate(surf_zi)
    deallocate(surf_area)
    deallocate(col_nodes)

  end subroutine set_mesh

  !------------------------------------------------------------------------
  subroutine set_material_properties(this, l2e_init_list, bounds_clump)
    !
    ! !DESCRIPTION:
    !
    !use elm_instMod               , only : soilstate_vars
    !use elm_instMod               , only : soilhydrology_vars
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
    real(r8), pointer    :: elm_hksat(:,:)          ! [mmH2O/s]
    real(r8), pointer    :: elm_bsw(:,:)
    real(r8), pointer    :: elm_sucsat(:,:)
    real(r8), pointer    :: elm_eff_porosity(:,:)
    real(r8), pointer    :: elm_zwt(:)
    real(r8), pointer    :: ats_watsat(:,:)
    real(r8), pointer    :: ats_hksat(:,:)         !
    real(r8), pointer    :: ats_bsw(:,:)
    real(r8), pointer    :: ats_sucsat(:,:)
    real(r8), pointer    :: ats_eff_porosity(:,:)
    real(r8), pointer    :: ats_residual_sat(:,:)
    !
    integer              :: fc, c, j
    integer              :: nz

    !-----------------------------------------------------------------------

    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_watsatc      , elm_watsat       )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_hksatc       , elm_hksat        )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_bswc         , elm_bsw          )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_sucsatc      , elm_sucsat       )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_parameter_effporosityc , elm_eff_porosity )

    ! Allocate memory (on local rank only), note all dimensions are 0-based
    nz = size(elm_watsat(1,:))
    allocate(ats_watsat       (0:this%filter_col_num-1, 0:nz-1 ))
    allocate(ats_hksat        (0:this%filter_col_num-1, 0:nz-1 ))
    allocate(ats_bsw          (0:this%filter_col_num-1, 0:nz-1 ))
    allocate(ats_sucsat       (0:this%filter_col_num-1, 0:nz-1 ))
    allocate(ats_eff_porosity (0:this%filter_col_num-1, 0:nz-1 ))
    allocate(ats_residual_sat (0:this%filter_col_num-1, 0:nz-1 ))

    ! Initialize
    ats_watsat       (:,:) = 0._r8
    ats_hksat        (:,:) = 0._r8
    ats_bsw          (:,:) = 0._r8
    ats_sucsat       (:,:) = 0._r8
    ats_eff_porosity (:,:) = 0._r8
    ats_residual_sat (:,:) = 0._r8  ! alway zero from ELM

    do fc = 1, this%filter_col_num
      c = this%filter_col(fc)
      do j = 1, nz
         ats_watsat(fc-1,j-1)       = elm_watsat(c,j)
         ats_hksat(fc-1,j-1)        = elm_hksat(c,j)
         ats_bsw(fc-1,j-1)          = elm_bsw(c,j)
         ats_sucsat(fc-1,j-1)       = elm_sucsat(c,j)
         ats_eff_porosity(fc-1,j-1) = elm_eff_porosity(c,j)
      end do
    end do

    ! pass data to ATS
    call this%ats_interface%setmat(ats_watsat, ats_hksat, ats_bsw, ats_sucsat, &
             ats_residual_sat, ats_eff_porosity)

    deallocate(ats_watsat)
    deallocate(ats_hksat)
    deallocate(ats_bsw)
    deallocate(ats_sucsat)
    deallocate(ats_eff_porosity)
    deallocate(ats_residual_sat)

  end subroutine set_material_properties

  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
   subroutine EM_ATS_OneStep(this, em_stage, dt, nstep, clump_rank, l2e_list, e2l_list, &
        bounds_clump)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_ats_type)                   :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt             ! unit: seconds
    integer              , intent(in)    :: nstep
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump

    !-----------------------

    ! reset ICs
    if (nstep==0) then
      call set_initial_conditions(this, l2e_list, bounds_clump)
    end if

    ! BCs and source-sink terms for each ELM-timestep
    call set_bc_ss(this, l2e_list, bounds_clump)

    ! run one ELM-timestep
    call  this%ats_interface%onestep(nstep, dt)

    ! NOT-yet finished
    call get_data_for_elm(this, e2l_list, bounds_clump)

  end subroutine EM_ATS_OneStep

  !------------------------------------------------------------------------
  subroutine set_initial_conditions(this, l2e_list, bounds_clump)
    !
    ! !DESCRIPTION:
    !
    !use elm_instMod               , only : soilstate_vars
    !use elm_instMod               , only : soilhydrology_vars
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_ats_type)               :: this
    class(emi_data_list), intent(in) :: l2e_list
    type(bounds_type)   , intent(in) :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    integer           :: fc, c, j
    integer           :: nz
    real(r8), pointer :: patm(:), patm_ats(:)
    real(r8), pointer :: soilp(:,:), soilp_ats(:,:)
    real(r8), pointer :: wtd(:), wtd_ats(:)
    !-----------------------------------------------------------------------
    !call l2e_list%GetPointerToReal1D(this%index_l2e_state_forc_pbot        , patm  ) ! TODO
    allocate(patm_ats(0:this%filter_col_num-1))
    patm_ats(:) = 101325.0_r8

    call l2e_list%GetPointerToReal2D(this%index_l2e_state_soilp            , soilp )
    nz = size(soilp(1,:))
    allocate(soilp_ats(0:this%filter_col_num-1, 0:nz-1))                               ! index starting from 0 at soil bottom

    !call l2e_list%GetPointerToReal1D(this%index_l2e_state_wtd             , wtd   )   ! TODO
    allocate(wtd_ats(0:this%filter_col_num-1))
    wtd_ats(:) = 100.0_r8 ! TODO

    do fc = 1, this%filter_col_num
      c = this%filter_col(fc)
      do j = 1, nz
        soilp_ats(fc-1,j-1) = soilp(c,j)
      end do
    end do
    call this%ats_interface%setic(patm_ats, soilp_ats, wtd_ats)

    deallocate(patm_ats)
    deallocate(soilp_ats)
    deallocate(wtd_ats)
  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine set_bc_ss(this, l2e_list, bounds_clump)
    !
    ! !DESCRIPTION:
    !
    !use elm_instMod               , only : soilstate_vars
    !use elm_instMod               , only : soilhydrology_vars
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_ats_type)               :: this
    class(emi_data_list), intent(in) :: l2e_list
    type(bounds_type)   , intent(in) :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    integer              :: c, fc,  j            ! do loop indices
    integer              :: nc, nz
    real(r8), pointer    :: col_dz(:,:)
    real(r8), pointer    :: soilevap_flux(:),    soilinfl_flux(:),       soilevap_flux_ats(:),    soilinfl_flux_ats(:)
    real(r8), pointer    :: soilbot_flux(:),                             soilbot_flux_ats(:)
    real(r8), pointer    :: roottran_flux(:,:),                          roottran_flux_ats(:,:)
    real(r8), pointer    :: soildrain_flux(:,:),                         soildrain_flux_ats(:,:)
    !-----------------------------------------------------------------------
    nc = this%filter_col_num

    call l2e_list%GetPointerToReal2D(this%index_l2e_column_dz                   , col_dz        )    ! layer thickness (1:nlevgrnd)

    ! soil water source/sink terms. NOTE: here all variables are the potential rather than actual.
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_gross_infil            , soilinfl_flux )    ! mmH2O/s, + to soil
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_gross_evap             , soilevap_flux )    ! mmH2O/s, + to atm
    allocate(soilinfl_flux_ats(0:nc-1))                            ! index starting from 0, unit: kgH2O/m3/s
    allocate(soilevap_flux_ats(0:nc-1))

    allocate(soilbot_flux_ats(0:nc-1))
    soilbot_flux_ats(:) = 0._r8

    call l2e_list%GetPointerToReal2D(this%index_l2e_flux_tran_vr                , roottran_flux )   ! mmH2O/s, + to atm
    nz = size(roottran_flux(1,:))
    allocate(roottran_flux_ats(0:nc-1,0:nz-1))                    ! index starting from 0, unit: kgH2O/m3/s

    call l2e_list%GetPointerToReal2D(this%index_l2e_flux_drain_vr               , soildrain_flux )
    nz = size(soildrain_flux(1,:))
    allocate(soildrain_flux_ats(0:nc-1,0:nz-1))                   ! index starting from 0, unit: kgH2O/m3/s

    do fc = 1, nc
      c = this%filter_col(fc)
      soilinfl_flux_ats(fc-1) = soilinfl_flux(c)*denh2o/col_dz(c,1)/1000.0        !  mm/s = 0.001m3/m2/s *kg/m3 /m
      soilevap_flux_ats(fc-1) = soilevap_flux(c)*denh2o/col_dz(c,1)/1000.0
      do j = 1, nz
        roottran_flux_ats(fc-1,j-1) = roottran_flux(c,j)*denh2o/col_dz(c,j)/1000.0
        soildrain_flux_ats(fc-1,j-1) = soildrain_flux(c,j)*denh2o/col_dz(c,j)/1000.0
      end do
    end do

    call this%ats_interface%setss(soilinfl_flux_ats, soilevap_flux_ats, soilbot_flux_ats, &
                                      roottran_flux_ats, soildrain_flux_ats)

  end subroutine set_bc_ss

  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------

  subroutine get_data_for_elm(this, e2l_list, bounds_clump)
    !
    !DESCRIPTION
    !  Saves
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_ats_type)                   :: this
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump

    ! !LOCAL VARIABLES:
    integer :: bounds_proc_begc, bounds_proc_endc
    !-----------------------------------------------------------------------

    bounds_proc_begc     = bounds_clump%begc
    bounds_proc_endc     = bounds_clump%endc

    ! currently only checking the ATS data as it is, i.e. NO returning data for ELM
    call this%ats_interface%getdata()

   end subroutine get_data_for_elm

  !------------------------------------------------------------------------
   subroutine EM_ATS_Finalize(this, dt, nstep, clump_rank, l2e_list, e2l_list, &
        bounds_clump)
    !
    ! !DESCRIPTION:
    !
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

    ! (TODO)
    if (associated(this%filter_col)) deallocate(this%filter_col)

  end subroutine EM_ATS_Finalize

 !------------------------------------------------------------------------

#endif

end module ExternalModelATSMod
