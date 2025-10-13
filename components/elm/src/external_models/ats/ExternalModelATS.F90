module ExternalModelATS

  ! Design description
  ! ==================
  !
  ! ATS implements surface and subsurface hydrology, including:
  !  - lateral surface flows,
  !  - lateral and vertical subsurface flows,
  !  - in/ex-filtration, e.g. fluxes between surface/subsurface,
  !  - downregulation of potential evaporation based on free water
  !  - downregulation (based on stomatal closure and water
  !    availability) and distribution to vertical cells of potential
  !    transpiration
  !
  ! To do this, ATS expects ELM to provide the following quantities
  ! (note that all fluxes and state variables are per unit column
  ! area):
  !
  ! Parameters (in ATS units/verbiage)
  ! ----------------------------------
  !  - base porosity [-]
  !    The soil porosity, unaffected by water-based compressibility or ice
  !  - hydraulic conductivity [mm s^-1]
  !  - Clapp & Hornberger b [-]
  !  - Clapp & Hornberger Psi_sat [Pa]
  !    AKA air entry pressure
  !  - residual saturation [-]
  !  - effective porosity [-]
  !    Pore-space excluding ice? Including pore expansion /
  !    compressibility? (See #1).
  !  - PFT [index]
  !    The primary PFT on this grid cell.
  !  - root fraction [-]
  !    Fraction of roots in each vertical cell of PFT.
  !
  ! State Variables
  ! ---------------
  !
  !  - surface water content [mm]
  !    Per unit column area 
  !  - subsurface water content [mm]
  !    Per unit column area
  ! 
  ! Fluxes
  ! ------
  !  - surface water source [mm s^-1]

  !    Positive toward the ground surface, the sum of all gross
  !    influxes of water to the top of the surface water/bare ground
  !    but below the snowpack, including throughfall, canopy
  !    dripping/drainage, irrigation, snowmelt/snow water drainage,
  !    but NOT evaporation or sublimation of ice.  Does this include
  !    condensation/dew as a positive flux (see #2 below)?
  !
  !
  !  - potential evaporation [mm s^-1]
  !
  !    Positive toward the atmosphere, the potential flux of water
  !    from both bare soil and surface water.  Does this include
  !    condensation as a negative flux (see #2 below)?
  !
  !  - potential transpiration [mm s^-1]
  !
  !    Positive toward the atmosphere, the potential flux of water
  !    through the vegatation and used in photosynthesis.
  !
  ! Given these quantities, ATS will pass back ELM the following
  ! quantities:
  !
  ! State Variables
  ! ---------------
  !  - surface water content [mm]
  !    Per unit column area 
  !  - subsurface water content [m]
  !    Per unit column area
  !
  ! Fluxes
  ! ------
  !  - actual evaporation [mm s^-1]
  !    Downregulated evaporation, per unit column area
  !
  !  - actual transpiration [mm s^-1]
  !    Downregulated transpiration, per unit column area, in each
  !    subsurface cell.
  !
  !  - runoff generation [mm s^-1]
  !    Flux of water running off to a downstream cell
  !
  !  - baseflow [mm s^-1]
  !    Column-integrated flux of water leaving the subsurface column
  !    of cells due to lateral or bottom-boundary fluxes.
  !
  ! Diagnostics
  ! -----------
  ! ??
  !
  !
  ! Notes & Open questions
  ! ----------------------
  ! 1. How is ice dealt with?  If porosity does not include ice, then
  !    evaporation should not include sublimation of ice.  But do we
  !    get in and write back ELM's liquid water content or total water
  !    content?
  !
  ! 2. How is condensation/dew dealt with?  Does it go in surface water
  !    source or as negative evaporation?  Probably doesn't matter as
  !    long as we get it somewhere?
  !
  ! 3. Where are area fractions of surface water and bare ground
  !    computed?  ATS or ELM?  And where is surface water converted to
  !    ponded depth?  Need this be consistent or not?
  !
  ! 4. Molar density is assumed to be 55500 mol / m^3.
  !
  !
  !
  ! Map to ELM
  ! ==========
  !
  ! The above quantities are defined at the interface to the two
  ! codes.  Mapping these to ELM requires the following
  ! transformations and renaming conventions:
  !
  ! Parameters
  ! ----------
  ! - base porosity             : soilstate_vars%watsat_col
  ! - hydraulic conductivity    : soilstate_vars%hksat_col
  ! - C&H b                     : soilstate_vars%bsw_col
  ! - C&H Psi_sat               : soilstate_vars%sucsat_col
  ! - residual saturation       : <N/A> (always 0 in C&H)
  ! - effective porosity        : soilstate_vars%eff_porosity_col
  ! - PFT                       : p such that veg_pp%active(p) and maximum
  !                               over p of veg_pp%wtcol(p)
  ! - root fraction             : soilstate_vars%rootfr_patch(PFT)
  !
  !
  ! ELM to ATS State Variables
  ! --------------------------
  ! - surface water content     : col_ws%h2osfc
  ! - subsurface water content  : col_ws%h2osoi_liq
  !
  ! ELM to ATS Fluxes
  ! -----------------
  ! - surface water source      : col_wf%qflx_top_soil
  !                             : Computed:
  !                             : Used: SoilHydrologyMod.F90:431
  ! - potential evaporation     : = (1 - fsno - frac_h2osfc) * qflx_evap +
  !                             :   frac_h2osfc * qflx_ev_h2osfc
  !                             : Computed: SoilHydrologyMod.F90:436 and
  !                             :      SoilHydrologyMod:F90:403, but what
  !                             :      is qflx_evap_grnd vs qflx_ev_soil?
  !                             : Used: as qflx_in_soil and qflx_in_h2osfc
  ! - potential transpiration   : col_wf%qflx_tran_veg
  !
  !
  ! ATS to ELM State Variables
  ! --------------------------
  ! same as ELM to ATS
  !
  !
  ! ATS to ELM Fluxes
  ! -----------------
  !
  ! - actual evaporation        : must be allocated back somehow?  What
  !                             : happens if downregulation changes the
  !                             : energy balance?
  !
  ! - actual transpiration      : given in a vertical-cell basis?  Or just
  !                             : total?  Must ramp back photosynthesis/carbon
  !                             : cycle implications?
  !
  ! Diagnostics
  ! -----------
  ! ??
  !
  ! Progression of Capability
  ! =========================
  ! Step 0: compilation
  ! Step 1: skip passing soil properties, ignore initial conditions,
  !         ignore pressure, ignore downregulation, ignore ELM -> ATS
  !         water state
  ! Step 2: add back in ELM --> ATS water state -- confirm equal?
  !         Likely needs porosity.
  ! Step 3: add back downregulation
  ! Step 4: add back initial conditions, consistent pressure + water
  !         content, at this point water balance checks should pass
  ! Step 5: consistent energy and carbon cycle checks
  ! Step 6: pass soil properties
  
#ifdef USE_ATS_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !

  ! ELM modules
  use shr_kind_mod                 , only : r8 => shr_kind_r8

  use histFileMod                  , only : hist_nhtfrq
  use elm_varpar                   , only : nlevgrnd


  use ExternalModelATS_readnlMod   , only : ats_inputdir, ats_inputfile
  use spmdMod                      , only : mpicom
  use ExternalModelATS_interface   , only : ats_create, ats_setup, ats_initialize, &
                                            ats_advance, ats_delete, &
                                            ats_set_field, ats_get_field

  use iso_c_binding

  ! fortran API for C++ ATS code
  use ExternalModelATS_variables                , only : ats_var_id
  
  !
  implicit none
  !
  ! EM data-type and procedures for elm-ats interface
  
  type, public :: em_ats_type
     ! ----------------------------------------------------------------------
     ! pointer to ATS instance
     ! ----------------------------------------------------------------------
     type(C_PTR) :: ats

   contains
     procedure, public :: Init                    => EM_ATS_Init
     procedure, public :: Advance                 => EM_ATS_Advance
     procedure, public :: Finalize                => EM_ATS_Finalize

  end type em_ats_type

  type(em_ats_type), public :: em_ats

contains
  !------------------------------------------------------------------------
  subroutine EM_ATS_Create(this)
    implicit none
    type(em_ats_type) :: this

    this = EM_ATS_Create2(ats_inputdir, ats_inputfile, mpicom)
  end subroutine EM_ATS_Create

  !------------------------------------------------------------------------
  function EM_ATS_Create2(input_dir, input_file, mpi_comm)
    implicit none
    type(em_ats_type) :: EM_ATS_Create2

    character(len=*), intent(in) :: input_dir
    character(len=*), intent(in) :: input_file
    integer, intent(in) :: mpi_comm

    ! local variables
    character(kind=C_CHAR) :: c_input_file(len_trim(input_dir)+len_trim(input_file)+2)
    integer :: i, n1, n2
    integer :: ierr

    print *, ''
    print *, '============================================================='
    print *,''
    print *, ' -------- ELM-ATS Coupled Mode ------------------------------'
    print *, ''
    print *, 'EM_ATS_Init: ats inputs - ', trim(input_dir), ' ', trim(input_file)
    print *, 'communicator id: ', mpi_comm
    
    ! ----------------------------------------------------------
    ! Converting Fortran-type input filename, incl. dir, to C-type
    n1 = len_trim(input_dir)
    do i = 1, n1
       c_input_file(i) = input_dir(i:i)
    end do
    c_input_file(n1+1:n1+1) = '/'
    n2 = len_trim(input_file)
    do i = 1, n2
       c_input_file(n1+1+i) = input_file(i:i)
    end do
    c_input_file(n1+n2+2) = C_NULL_CHAR

    EM_ATS_Create2%ats = ats_create(mpi_comm, c_input_file)
  end function EM_ATS_Create2
    
  !------------------------------------------------------------------------
  subroutine EM_ATS_Init(this)
    implicit none

    class(em_ats_type)                   :: this

    !! begin setup portion of this call
    ! check PFT assumption -- 1 PFT per column, 1 column per grid cell
    ! TODO: ETC -- Step 0
    ! call EM_ATS_check_heirarchy(this)

    ! check ATS -- ELM mesh consistency
    ! TODO: ETC -- Step 0
    ! call EM_ATS_check_mesh(this)

    ! ATS setup
    call ats_setup(this%ats)

    !! begin init portion of this call
    !
    ! TODO: ETC -- Step 4
    !
    ! initial parameters
    ! call ats_set_field(this%ats, ats_var_id%BASE_POROSITY, ...)
    ! call ats_set_field(this%ats, ats_var_id%HYDRAULIC_CONDUCTIVITY, ...)
    ! call ats_set_field(this%ats, ats_var_id%CLAPP_HORN_B, ...)
    ! call ats_set_field(this%ats, ats_var_id%CLAPP_HORN_PSISAT, ...)
    ! call ats_set_field(this%ats, ats_var_id%RESIDUAL_SATURATION, ...)

    ! initial conditions
    ! call ats_set_scalar(this%ats, ats_var_id%TIME, ...)
    ! call ats_set_field(this%ats, ats_var_id%WATER_CONTENT, ...)
    ! call ats_set_field(this%ats, ats_var_id%SURFACE_WATER_CONTENT, ...)

    ! ats init
    call ats_initialize(this%ats)
    
  end subroutine EM_ATS_Init



  !------------------------------------------------------------------------
  subroutine EM_ATS_Advance(this, dt, nstep, soilstate_vars, col_wf, col_ws)
    use SoilStateType              , only : soilstate_type
    use ColumnDataType             , only : column_water_flux, column_water_state

    implicit none

    class(em_ats_type)                       :: this
    real(r8)                 , intent(in)    :: dt
    integer                  , intent(in)    :: nstep
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(column_water_flux)  , intent(in)    :: col_wf
    type(column_water_state) , intent(in)    :: col_ws

    
    ! local variables
    logical(C_BOOL) :: do_vis
    logical(C_BOOL) :: do_checkpoint

    print *, 'EM_ATS_Advance: advancing for time', dt

    ! pass dynamic parameters to ATS
    call ats_set_field(this%ats, ats_var_id%EFFECTIVE_POROSITY, soilstate_vars%eff_porosity_col)

    ! is this correct?  Is ELM allocated with:
    !      num_patches = 17 * num_columns, or
    !      sum_columns(n_active_patches_in_column)
    ! if the former, a special parser needs to be implemented to
    ! downselect to num_columns, either on ATS or ELM side
    !
    ! Probably this should be done on the ELM side (since this might
    ! allow us to relax some assumptions later) but without memory
    ! allocation, so should be a combination of:
    !   ats_get_field_ptr_w()
    !   fill_with_rooting_info()

    ! TODO: ETC -- Step 4
    ! call ats_set_field(this%ats, ats_var_id%ROOT_FRACTION, soilstate_vars%rootfr_patch)

    ! pass state to ATS
    ! TODO: ETC -- Step 2
    ! call ats_set_field(this%ats, ats_var_id%WATER_CONTENT, col_ws%h2osoi_liq)
    ! call ats_set_field(this%ats, ats_var_id%SURFACE_WATER_CONTENT, col_ws%h2osfc)

    ! pass fluxes to ATS
    call ats_set_field(this%ats, ats_var_id%GROSS_SURFACE_WATER_SOURCE, col_wf%qflx_top_soil)
    call ats_set_field(this%ats, ats_var_id%POTENTIAL_TRANSPIRATION, col_wf%qflx_tran_veg)

    ! This one needs to be computed, and how it is computed is really an ELM thing.
    ! Therefore, we should really use ats_get_field_ptr_w() and fill.
    ! = (1 - fsno - frac_h2osfc) * qflx_evap +
    !                          frac_h2osfc * qflx_ev_h2osfc
    !
    ! note, fsno and qflx_evap are locally computed things themselves,
    ! see SoilHydrologyMod.F90:403.  This awaits input from Peter Thornton.
    ! TODO: ETC -- Step 1
    ! evap = ats_get_field_ptr_w(this%ats, ats_var_id%POTENTIAL_EVAPORATION)
    ! call EM_ATS_compute_evaporation()
    
    ! call advance
    ! -- get whether to vis and checkpoint
    do_vis = .false.
    do_checkpoint = .false.
    if (hist_nhtfrq(1) > 0) then
       if (mod(nstep, hist_nhtfrq(1)) == 0) do_vis = .true.
    end if

    ! -- call the advance method
    call ats_advance(this%ats, dt, do_checkpoint, do_vis)

    ! pass state back to ELM
    call ats_get_field(this%ats, ats_var_id%WATER_CONTENT, col_ws%h2osoi_liq)
    ! TODO: ETC -- Step 4
    !call ats_get_field(this%ats, PRESSURE, ...)
    call ats_get_field(this%ats, ats_var_id%SURFACE_WATER_CONTENT, col_ws%h2osfc)

    ! get actual fluxes back

    ! -- evaporation downregulated by soil water availability... what
    !    happens if this is less than requested?
    ! TODO: ETC -- Step 3
    !call ats_get_field(this%ats, ats_var_id%EVAPORATION, ...)

    ! -- transpiration downregulated by soil water availability... what
    !    happens if this is less than requested?
    ! TODO: ETC -- Step 3
    ! call ats_get_field(this%ats, ats_var_id%COLUMN_TRANSPIRATION, col_wf%qflx_tran_veg)

    ! diagnostics?
    ! ...
  end subroutine EM_ATS_Advance
  
  
  !------------------------------------------------------------------------
  subroutine EM_ATS_Finalize(this)
    implicit none
    class(em_ats_type) :: this

    print *, 'EM_ATS_Finalize: cleaning up'
    call ats_delete(this%ats)
  end subroutine EM_ATS_Finalize


  !========================================================================
  ! private functions
  !

  ! -----------------------------------------------------------------------
  ! Checks that assumptions about the ELM scale heirarchy required for
  ! use of ATS are satsified
  !
  subroutine EM_ATS_check_heirarchy(this)

    implicit none

    class(em_ats_type)                   :: this

    ! checks ELM-only things, e.g. num active PFTs on each water
    ! column is 1, num water columns on each grid cell is 1, no
    ! filters are active, etc
  end subroutine EM_ATS_check_heirarchy
  
  ! -----------------------------------------------------------------------
  ! Checks that the ATS and ELM mesh representations are consistent
  !
  subroutine EM_ATS_check_mesh(this)

    implicit none

    class(em_ats_type)                   :: this

    ! local variables
    integer :: ats_ncols_local
    integer :: ats_ncols_global
    integer :: ats_nlevgrnd

    real(c_double), pointer :: ats_dzs(:)
    real(c_double), pointer :: ats_areas(:)
    real(c_double), pointer :: ats_elevations(:)

    ! compare to mesh info
    ! call ats_get_mesh_info(this%ats, ats_ncols_local, ats_ncols_global, ats_nlevgrnd, ats_dzs)

    ! assertions on shapes
    ! TODO: ETC -- Step 1
    ! assert(ats_ncols_local .eq. ...)
    ! assert(ats_ncols_global .eq. ...)
    ! assert(ats_nlevgrnd .eq. nlevgrnd)
    
    ! calls ats_get_mesh_info, compares nlevgrnd, compares ncolumns,
    ! compares areas, elevations? lat-lon?
    ! TODO: ETC -- Step 1
    ! ats_elevations = ats_get_field_ptr(this%ats, ats_var_id%ELEVATION)
    ! for (...) compare elevation

    ! TODO: ETC -- Step 1
    ! ats_areas = ats_get_field_ptr(this%ats, ats_var_id%AREA)
    ! for (...) compare area
    
  end subroutine EM_ATS_check_mesh


  ! ! -----------------------------------------------------------------------
  ! ! Helper function to hide how we do copies in.  Only fortran knows
  ! ! both real types -- how much has to be done here?
  ! !
  ! subroutine EM_ATS_copy_field_to_ats(this, var_id, field)
  !   implicit none
  !   class(em_ats_type)         :: this
  !   integer(c_int), intent(in) :: var_id
  !   real(r8), intent(in)       :: field(:)

  !   ! local variables
  !   real(c_double)        :: ats_field(:)

  !   ats_field = ats_get_field_ptr_w(this%ats, var_id)
  !   ats_field(:) = field(:)    
  ! end subroutine EM_ATS_copy_field_to_ats

  ! subroutine EM_ATS_copy_field_from_ats(this, var_id, field)
  !   implicit none
  !   class(em_ats_type)         :: this
  !   integer(c_int), intent(in) :: var_id
  !   real(r8), intent(out)      :: field(:)

  !   ! local variables
  !   real(c_double)        :: ats_field(:)

  !   ats_field = ats_get_field_ptr_w(this%ats, var_id)
  !   field(:) = ats_field(:)    
  ! end subroutine EM_ATS_copy_field_to_ats

#endif
  
end module ExternalModelATS
