module ExternalModelATS
  
#ifdef USE_ATS_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !

  ! ELM modules
  use iso_c_binding
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use spmdMod                      , only : mpicom
  use abortutils                   , only : endrun
  use histFileMod                  , only : hist_nhtfrq
  use elm_varpar                   , only : nlevgrnd
  use elm_varctl                   , only : iulog
  
  use ExternalModelATS_interface   , only : ats_create, ats_setup, ats_parse_parameter_list, &
                                            ats_initialize, ats_advance, ats_delete, &
                                            ats_get_field, ats_get_field_ptr_w, &
                                            ats_get_mesh_info, ats_set_scalar

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
     type(c_ptr) :: ats
     integer, pointer :: col_filter(:)
     integer :: ncolumns
     integer :: nlevgrnd
     
   contains
     procedure, public :: Init                          => EM_ATS_Init
     procedure, public :: Advance                       => EM_ATS_Advance
     procedure, public :: Finalize                      => EM_ATS_Finalize
     
  end type em_ats_type

  type(em_ats_type), public :: em_ats

contains
  !------------------------------------------------------------------------
  subroutine EM_ATS_Create(this)
    use ExternalModelATS_readnlMod   , only : ats_inputdir, ats_inputfile

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

    write(iulog,*) ''
    write(iulog,*) '============================================================='
    print *,''
    write(iulog,*) ' -------- ELM-ATS Coupled Mode ------------------------------'
    write(iulog,*) ''
    write(iulog,*) 'EM_ATS_Init: ats inputs - ', trim(input_dir), ' ', trim(input_file)
    write(iulog,*) 'communicator id: ', mpi_comm
    
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
  subroutine EM_ATS_Init(this, time0, filter, nclumps, ncolumns, nlevgrnd)
    use filterMod, only : clumpfilter

    implicit none

    class(em_ats_type)  :: this
    real(r8), intent(in):: time0
    type(clumpfilter), intent(in) :: filter(:)
    integer, intent(in) :: nclumps
    integer, intent(in) :: ncolumns
    integer, intent(in) :: nlevgrnd

    !! begin setup portion of this call
    ! check PFT assumption -- 1 PFT per column, 1 column per grid cell
    ! TODO: ETC -- Step 0
    call EM_ATS_CheckHeirarchy(this, filter, nclumps, ncolumns)

    ! keep the filter?  probably should be 1:1 now and not needed...
    this%ncolumns = filter(1)%num_hydrologyc
    this%col_filter => filter(1)%hydrologyc
    this%nlevgrnd = nlevgrnd

    call ats_parse_parameter_list(this%ats)

    ! check ATS -- ELM mesh consistency
    ! TODO: ETC -- Step 0
    call EM_ATS_CheckMesh(this)

    ! ATS setup
    call ats_setup(this%ats)
    call ats_set_scalar(this%ats, ats_var_id%TIME, time0)

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
  subroutine EM_ATS_Advance(this, dt, nstep, col_pp, soilstate_vars, col_wf, col_ws)
    use ColumnType                 , only : column_physical_properties
    use SoilStateType              , only : soilstate_type
    use ColumnDataType             , only : column_water_flux, column_water_state

    implicit none

    class(em_ats_type)                       :: this
    real(r8)                 , intent(in)    :: dt
    integer                  , intent(in)    :: nstep
    type(column_physical_properties)     , intent(in)    :: col_pp
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(column_water_flux)  , intent(in)    :: col_wf
    type(column_water_state) , intent(in)    :: col_ws

    
    ! local variables
    logical(C_BOOL) :: do_vis
    logical(C_BOOL) :: do_checkpoint

    write(iulog,*) 'EM_ATS_Advance: advancing for time', dt

    ! pass dynamic parameters to ATS
    ! call EM_ATS_SetField_CopySubsurface(this, "effective porosity", ats_var_id%EFFECTIVE_POROSITY, &
    !      soilstate_vars%eff_porosity_col)

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

    ! pass fluxes to ATS -- unit change from [mm H2O / s] to [m H2O / s]
    call EM_ATS_SetField_ScalarMultiply(this, "gross surface water source", ats_var_id%GROSS_SURFACE_WATER_SOURCE, &
         col_wf%qflx_top_soil, 1.e-3_r8)
    call EM_ATS_SetField_ScalarMultiply(this, "potential transpiration", ats_var_id%POTENTIAL_TRANSPIRATION, &
         col_wf%qflx_tran_veg, 1.e-3_r8)

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
    ! call ats_get_field(this%ats, ats_var_id%WATER_CONTENT, col_ws%h2osoi_liq)
    ! TODO: ETC -- Step 4
    !call ats_get_field(this%ats, PRESSURE, ...)
    ! call ats_get_field(this%ats, ats_var_id%SURFACE_WATER_CONTENT, col_ws%h2osfc)

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

    write(iulog,*) 'EM_ATS_Finalize: cleaning up'
    call ats_delete(this%ats)
  end subroutine EM_ATS_Finalize

  !------------------------------------------------------------------------
  ! straight copy
  !
  subroutine EM_ATS_SetField_Copy(this, var_id, field)
    implicit none
    class(em_ats_type)         :: this
    integer(c_int), intent(in) :: var_id
    real(r8), intent(in)       :: field(:)

    ! local variables
    integer :: i, ii
    real(c_double), pointer    :: ats_field(:)

    call EM_ATS_GetFieldPtrW(this, var_id, this%ncolumns, ats_field)
    do i=1,this%ncolumns
       ii = this%col_filter(i)
       ats_field(i) = field(ii)
    end do
  end subroutine EM_ATS_SetField_Copy

  !------------------------------------------------------------------------
  ! straight copy, subsurface
  !
  subroutine EM_ATS_SetField_CopySubsurface(this, var_name, var_id, field)
    implicit none
    class(em_ats_type)         :: this
    character(len=*), intent(in) :: var_name
    integer(c_int), intent(in) :: var_id
    real(r8), intent(in)       :: field(:,:)

    ! local variables
    integer :: i, ii, j
    real(c_double), pointer    :: ats_field(:)

    call EM_ATS_GetFieldPtrW(this, var_id, this%ncolumns, ats_field)
    do i=1,this%ncolumns
       ii = this%col_filter(i)
       do j=1,this%nlevgrnd
          ats_field((i-1) * this%nlevgrnd + j) = field(ii,j)
          write(iulog,*) var_name, " : column ", i, " cell ", j, " = ", field(ii,j)
       end do
    end do
  end subroutine EM_ATS_SetField_CopySubsurface

  
  !------------------------------------------------------------------------
  ! unit change
  !
  subroutine EM_ATS_SetField_ScalarMultiply(this, var_name, var_id, field, factor)
    implicit none
    class(em_ats_type)         :: this
    character(len=*), intent(in) :: var_name
    integer(c_int), intent(in) :: var_id
    real(r8), intent(in)       :: field(:)
    real(r8), intent(in)       :: factor

    ! local variables
    real(c_double), pointer    :: ats_field(:)
    integer :: i, ii

    call EM_ATS_GetFieldPtrW(this, var_id, this%ncolumns, ats_field)
    do i=1,this%ncolumns
       ii = this%col_filter(i)
       ats_field(i) = factor * field(ii)
       write(iulog,*) var_name, " : column ", i, " = ", ats_field(i)
    end do
  end subroutine EM_ATS_SetField_ScalarMultiply


  !========================================================================
  ! private functions
  !

  ! -----------------------------------------------------------------------
  ! Checks that assumptions about the ELM scale heirarchy required for
  ! use of ATS are satsified
  !
  subroutine EM_ATS_CheckHeirarchy(this, filter, nclumps, ncolumns_all)
    use filterMod, only : clumpfilter
    implicit none

    class(em_ats_type)                   :: this
    type(clumpfilter), intent(in) :: filter(:)
    integer, intent(in) :: nclumps
    integer, intent(in) :: ncolumns_all

    integer :: ncolumns
    
    ! only 1 clump or memory is not contiguous
    if (nclumps /= 1) then
       call endrun("ATS only works with 1 clump.")
    end if

    if (ncolumns /= filter(1)%num_hydrologyc) then
       write(iulog,*) "WARNING: ATS does not support non-hydrology columns in spatially explicit mode.  Presuming implicit-mode or URBAN_REGION_ID /= 0."
    end if

    ncolumns = filter(1)%num_hydrologyc

    ! number of soil columns == number of hydrology columns == number
    ! of soil PFTs
    if (ncolumns /= filter(1)%num_soilc) then
       call endrun("ATS expects all soil columns to be hydrology columns.")
    end if

    if (ncolumns /= filter(1)%num_soilp) then
       call endrun("ATS expects one PFT per column.")
    end if

    if (filter(1)%num_lakec /= 0) then
       call endrun("ATS does not support lakes.")
    end if
    if (filter(1)%num_urbanc /= 0) then
       call endrun("ATS does not support urban.")
    end if

    ! checks ELM-only things, e.g. num active PFTs on each water
    ! column is 1, num water columns on each grid cell is 1, no
    ! filters are active, etc
  end subroutine EM_ATS_CheckHeirarchy
  
  ! -----------------------------------------------------------------------
  ! Checks that the ATS and ELM mesh representations are consistent
  !
  subroutine EM_ATS_CheckMesh(this)
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
    call ats_get_mesh_info(this%ats, ats_ncols_local, ats_ncols_global, ats_nlevgrnd, ats_dzs)

    ! assertions on shapes
    if (ats_ncols_local /= this%ncolumns) then
       call endrun("ATS local ncolumns does not match requested ncolumns")
    end if

    if (ats_nlevgrnd /= this%nlevgrnd) then
       print*, "ATS: NLEVGRND = ", ats_nlevgrnd
       print*, "ELM: NLEVGRND = ", this%nlevgrnd
       call endrun("ATS local cells in the vertical does not match ELM nlevgrnd")
    end if

    ! check dzs?  Note this will need more from ELM

    ! check areas of each column? Also needs more from ELM
    
    ! calls ats_get_mesh_info, compares nlevgrnd, compares ncolumns,
    ! compares areas, elevations? lat-lon?
    ! TODO: ETC -- Step 1
    ! ats_elevations = ats_get_field_ptr(this%ats, ats_var_id%ELEVATION)
    ! for (...) compare elevation

    ! TODO: ETC -- Step 1
    ! ats_areas = ats_get_field_ptr(this%ats, ats_var_id%AREA)
    ! for (...) compare area
    
  end subroutine EM_ATS_CheckMesh


  ! -----------------------------------------------------------------------
  ! Assign a field to the c_ptr
  !
  subroutine EM_ATS_GetFieldPtrW(this, var_id, size, field)
    implicit none

    class(em_ats_type)                   :: this
    integer(c_int), intent(in) :: var_id
    integer, intent(in) :: size
    real(c_double), intent(inout), pointer :: field(:)

    ! local variables
    type(c_ptr) :: field_ptr
    
    
    field_ptr = ats_get_field_ptr_w(this%ats, var_id)
    call c_f_pointer(field_ptr, field, [size])
  end subroutine EM_ATS_GetFieldPtrW
  
#endif
  
end module ExternalModelATS
