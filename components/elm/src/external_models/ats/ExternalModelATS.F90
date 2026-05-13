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
  use elm_varctl                   , only : use_ats, use_ats_ic

  use ColumnType                 , only : column_physical_properties
  use GridcellType               , only : gridcell_physical_properties_type
  use PhotosynthesisType         , only : photosyns_type
  use SoilHydrologyType          , only : soilhydrology_type
  use SoilStateType              , only : soilstate_type
  use ColumnDataType             , only : column_water_state, column_water_flux
  use ExternalModelATS_interface

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
     integer, pointer :: pft_filter(:)
     integer :: ncolumns
     integer :: nlevgrnd
     real(r8) :: p_atm
     integer :: verbosity
     
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
    this%verbosity = 1
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
    character(kind=C_CHAR) :: c_log_file(1)
    integer :: i, n1, n2
    integer :: ierr

    write(iulog,*) ''
    write(iulog,*) '============================================================='
    print *,''
    write(iulog,*) ' -------- ELM-ATS Coupled Mode ------------------------------'
    write(iulog,*) ''
    write(iulog,*) 'EM_ATS_Init: ats inputs - ', trim(input_dir), ' ', trim(input_file)
    write(iulog,*) 'communicator id: ', mpi_comm
    
    EM_ATS_Create2%p_atm = 101325._r8

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

    c_log_file(1) = C_NULL_CHAR

    EM_ATS_Create2%ats = ats_create(mpi_comm, c_input_file, c_log_file)
    call ats_parse_parameter_list(EM_ATS_Create2%ats)
  end function EM_ATS_Create2

  !------------------------------------------------------------------------
  subroutine EM_ATS_GetMeshInfo(this, ncolumns)
    implicit none
    class(em_ats_type)  :: this
    integer, intent(out) :: ncolumns
    call ats_get_mesh_info(this%ats, ncolumns)
  end subroutine EM_ATS_GetMeshInfo


  !------------------------------------------------------------------------
  subroutine EM_ATS_Init(this, time0, filter, nclumps, ncolumns, nlevgrnd, &
       grc_pp, col_pp, soilstate_vars, soilhydrology_vars, col_ws)
    use filterMod, only : clumpfilter

    implicit none

    class(em_ats_type)  :: this
    real(r8), intent(in):: time0
    type(clumpfilter), intent(in) :: filter(:)
    integer, intent(in) :: nclumps
    integer, intent(in) :: ncolumns
    integer, intent(in) :: nlevgrnd
    type(gridcell_physical_properties_type) , intent(in) :: grc_pp
    type(column_physical_properties)     , intent(in)    :: col_pp
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(soilstate_type)     , intent(inout)    :: soilstate_vars
    type(column_water_state) , intent(inout)    :: col_ws

    ! locals
    integer :: i, ii

    !! begin setup portion of this call
    ! check PFT assumption -- 1 PFT per column, 1 column per grid cell
    ! TODO: ETC -- Step 0
    call EM_ATS_CheckHierarchy(this, filter, nclumps, ncolumns)

    ! keep the filter?  probably should be 1:1 now and not needed...
    this%ncolumns = filter(1)%num_hydrologyc
    this%col_filter => filter(1)%hydrologyc
    this%pft_filter => filter(1)%soilp
    this%nlevgrnd = nlevgrnd

    ! check ATS -- ELM mesh consistency
    ! TODO: ETC -- Step 0
    call EM_ATS_CheckMesh(this, grc_pp, col_pp)

    ! ATS setup
    call ats_setup(this%ats)
    call ats_set_scalar(this%ats, ats_var_id%TIME, time0)

    ! zero out things that are currently nan, do not get set (because
    ! they are set in blocks of code turned off via use_ats), but do
    ! get included in various water balance calls.
    ! soilhydrology_vars%wa_col(:) = 0._r8
    
    
    !! begin init portion of this call
    !
    ! TODO: ETC -- Step 4
    !
    ! initial parameters
    !
    ! note we initially copy over total porosity, not effective
    ! porosity, because it is uninitialized at the first step.
    call EM_ATS_SetField_CopySubsurface(this, "base porosity", ats_var_id%BASE_POROSITY, &
         soilstate_vars%watsat_col)
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

    ! copy back initial water contents
    call EM_ATS_GetField_WaterContent(this, col_pp, soilstate_vars, col_ws, soilhydrology_vars)

    ! additionally re-save old water
    do i = 1, this%ncolumns
       ii = this%col_filter(i)
       col_ws%h2osoi_liq_old(ii,:) = col_ws%h2osoi_liq(ii,:)
    end do
  end subroutine EM_ATS_Init



  !------------------------------------------------------------------------
  subroutine EM_ATS_Advance(this, dt, nstep, col_pp, soilstate_vars, col_ws, col_wf, soilhydrology_vars, photosyns_vars)
    implicit none

    class(em_ats_type)                       :: this
    real(r8)                 , intent(in)    :: dt
    integer                  , intent(in)    :: nstep
    type(column_physical_properties)  , intent(in)    :: col_pp
    type(soilstate_type)     , intent(inout)    :: soilstate_vars
    type(column_water_state) , intent(inout)    :: col_ws
    type(column_water_flux)  , intent(inout)    :: col_wf
    type(photosyns_type)     , intent(inout) :: photosyns_vars
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
  
    
    ! local variables
    logical(C_BOOL) :: do_vis
    logical(C_BOOL) :: do_checkpoint
    real(r8) :: h2osoi_liq_total
    integer :: i, ii, j

    if (this%verbosity >= 1) then
       write(iulog,*) 'EM_ATS_Advance: advancing for time', dt
    end if

    ! == DEBUG CODE ====
    ! how much has ELM drifted from ATS?  Any changes to h2osoi_liq at
    ! this point will be overwritten, and if they are compensating
    ! changes elsewhere, they will result in mass imbalance.
    do i=1,this%ncolumns
       ii = this%col_filter(i)
       h2osoi_liq_total = 0._r8
       do j=1,this%nlevgrnd
          h2osoi_liq_total = h2osoi_liq_total + col_ws%h2osoi_liq(ii,j)
       end do
       if (this%verbosity >= 1) then
          write(iulog,*) "ELM: h2osoi_liq total : column ", i, " total = ", h2osoi_liq_total
          write(iulog,*) "ELM: h2osfc : column ", i, " total = ", col_ws%h2osfc(ii)
       end if
    end do
    ! == END DEBUG CODE ====
    
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

    ! pass state to ATS -- this is required in hydro-only mode because
    ! freeze-thaw may change water to ice
    call EM_ATS_SetField_WaterContent(this, col_pp, col_ws, soilstate_vars)

    ! pass fluxes to ATS -- unit change from [mm H2O / s] to [m H2O / s]
    call EM_ATS_SetField_ScalarMultiply(this, "gross surface water source", ats_var_id%GROSS_SURFACE_WATER_SOURCE, &
         col_wf%qflx_gross_infl_soil, 1.e-3_r8, 1)
    call EM_ATS_SetField_ScalarMultiply(this, "potential transpiration", ats_var_id%POTENTIAL_TRANSPIRATION, &
         col_wf%qflx_tran_veg, 1.e-3_r8, 1)
    call EM_ATS_SetField_ScalarMultiply(this, "potential evaporation", ats_var_id%POTENTIAL_EVAPORATION, &
         col_wf%qflx_gross_evap_soil, 1.e-3_r8, 1)
    
    ! call advance
    ! -- get whether to vis and checkpoint
    do_vis = .false.
    do_checkpoint = .false.
    if (hist_nhtfrq(1) > 0) then
       if (mod(nstep, hist_nhtfrq(1)) == 0) do_vis = .true.
    end if

    ! -- call the advance method
    write(iulog,*) "CALLING ats_advance"
    call ats_advance(this%ats, dt, do_checkpoint, do_vis)

    ! pass state back to ELM
    ! TODO: ETC -- Step 4
    ! note we use the same porosity as set above
    call EM_ATS_GetField_WaterContent(this, col_pp, soilstate_vars, col_ws, soilhydrology_vars)
    
    ! get actual fluxes back and apply the downregulation where needed
    call EM_ATS_GetField_ActualEvaporation(this, col_wf)
    call EM_ATS_GetField_ActualTranspiration(this, col_wf, photosyns_vars)

    ! get runoff and baseflow back
    call EM_ATS_GetField_ScalarMultiply(this, "baseflow", ats_var_id%BASEFLOW, col_wf%qflx_drain, 1.e3_r8, 1)
    call EM_ATS_GetField_ScalarMultiply(this, "runoff", ats_var_id%RUNOFF, col_wf%qflx_surf, 1.e3_r8, 1)
    
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
  subroutine EM_ATS_GetField_Copy(this, var_name, var_id, field, verbosity)
    implicit none
    class(em_ats_type)         :: this
    character(len=*), intent(in) :: var_name
    integer(c_int), intent(in) :: var_id
    real(r8), intent(inout)       :: field(:)
    integer, intent(in), optional :: verbosity

    ! local variables
    integer :: i, ii
    real(c_double), pointer    :: ats_field(:)
    integer :: lverbosity

    if (present(verbosity)) then
       lverbosity = verbosity
    else
       lverbosity = 2
    end if

    call EM_ATS_GetFieldPtr(this, var_id, this%ncolumns, ats_field)
    do i=1,this%ncolumns
       ii = this%col_filter(i)
       field(i) = ats_field(i)
       if (this%verbosity >= lverbosity) then
          write(iulog,*) "GET: ", var_name, " : column ", i, " = ", field(i)
       end if
    end do
  end subroutine EM_ATS_GetField_Copy

  subroutine EM_ATS_SetField_Copy(this, var_name, var_id, field, verbosity)
    implicit none
    class(em_ats_type)         :: this
    character(len=*), intent(in) :: var_name
    integer(c_int), intent(in) :: var_id
    real(r8), intent(in)       :: field(:)
    integer, intent(in), optional :: verbosity

    ! local variables
    integer :: i, ii
    real(c_double), pointer    :: ats_field(:)
    integer :: lverbosity

    if (present(verbosity)) then
       lverbosity = verbosity
    else
       lverbosity = 2
    end if
    
    call EM_ATS_GetFieldPtrW(this, var_id, this%ncolumns, ats_field)
    do i=1,this%ncolumns
       ii = this%col_filter(i)
       ats_field(i) = field(i)
       if (this%verbosity >= lverbosity) then
          write(iulog,*) "SET: ", var_name, " : column ", i, " = ", field(i)
       end if
    end do
  end subroutine EM_ATS_SetField_Copy

  !------------------------------------------------------------------------
  ! straight copy, subsurface
  !
  subroutine EM_ATS_GetField_CopySubsurface(this, var_name, var_id, field, verbosity)
    implicit none
    class(em_ats_type)         :: this
    character(len=*), intent(in) :: var_name
    integer(c_int), intent(in) :: var_id
    real(r8), intent(inout)       :: field(:,:)
    integer, intent(in), optional :: verbosity

    ! local variables
    integer :: i, ii, j
    real(c_double), pointer    :: ats_field(:)
    integer :: lverbosity

    if (present(verbosity)) then
       lverbosity = verbosity
    else
       lverbosity = 2
    end if

    call EM_ATS_GetFieldPtr(this, var_id, this%ncolumns * this%nlevgrnd, ats_field)
    do i=1,this%ncolumns
       ii = this%col_filter(i)
       do j=1,this%nlevgrnd
          field(i, j) = ats_field((i-1) * this%nlevgrnd + j)
          if (this%verbosity >= lverbosity) then
             write(iulog,*) "GET: ", var_name, " : column ", i, " cell ", j, " = ", field(i,j)
          end if
       end do
    end do
  end subroutine EM_ATS_GetField_CopySubsurface

  subroutine EM_ATS_SetField_CopySubsurface(this, var_name, var_id, field, verbosity)
    implicit none
    class(em_ats_type)         :: this
    character(len=*), intent(in) :: var_name
    integer(c_int), intent(in) :: var_id
    real(r8), intent(in)       :: field(:,:)
    integer, intent(in), optional :: verbosity

    ! local variables
    integer :: i, ii, j
    real(c_double), pointer    :: ats_field(:)
    integer :: lverbosity

    if (present(verbosity)) then
       lverbosity = verbosity
    else
       lverbosity = 2
    end if

    call EM_ATS_GetFieldPtrW(this, var_id, this%ncolumns * this%nlevgrnd, ats_field)
    do i=1,this%ncolumns
       ii = this%col_filter(i)
       do j=1,this%nlevgrnd
          ats_field((i-1) * this%nlevgrnd + j) = field(i,j)
          if (this%verbosity >= lverbosity) then
             write(iulog,*) "SET:", var_name, " : column ", i, " cell ", j, " = ", field(i,j)
          end if
       end do
    end do
  end subroutine EM_ATS_SetField_CopySubsurface

  
  !------------------------------------------------------------------------
  ! unit change
  !
  subroutine EM_ATS_SetField_ScalarMultiply(this, var_name, var_id, field, factor, verbosity)
    implicit none
    class(em_ats_type)         :: this
    character(len=*), intent(in) :: var_name
    integer(c_int), intent(in) :: var_id
    real(r8), intent(in)       :: field(:)
    real(r8), intent(in)       :: factor
    integer, intent(in), optional :: verbosity

    ! local variables
    real(c_double), pointer    :: ats_field(:)
    integer :: i, ii
    integer :: lverbosity

    if (present(verbosity)) then
       lverbosity = verbosity
    else
       lverbosity = 2
    end if

    call EM_ATS_GetFieldPtrW(this, var_id, this%ncolumns, ats_field)
    do i=1,this%ncolumns
       ii = this%col_filter(i)
       ats_field(i) = factor * field(i)
       if (this%verbosity >= lverbosity) then
          write(iulog,*) "SET: ", var_name, " : column ", i, " = ", field(i)
       end if
    end do
  end subroutine EM_ATS_SetField_ScalarMultiply

  subroutine EM_ATS_GetField_ScalarMultiply(this, var_name, var_id, field, factor, verbosity)
    implicit none
    class(em_ats_type)         :: this
    character(len=*), intent(in) :: var_name
    integer(c_int), intent(in) :: var_id
    real(r8), intent(inout)       :: field(:)
    real(r8), intent(in)       :: factor
    integer, intent(in), optional :: verbosity

    ! local variables
    real(c_double), pointer    :: ats_field(:)
    integer :: i, ii
    integer :: lverbosity

    if (present(verbosity)) then
       lverbosity = verbosity
    else
       lverbosity = 2
    end if

    call EM_ATS_GetFieldPtr(this, var_id, this%ncolumns, ats_field)
    do i=1,this%ncolumns
       ii = this%col_filter(i)
       field(i) = factor * ats_field(i)
       if (this%verbosity >= lverbosity) then
          write(iulog,*) "GET: ", var_name, " : column ", i, " = ", field(i)
       end if
    end do
  end subroutine EM_ATS_GetField_ScalarMultiply
  

  !========================================================================
  ! private functions
  !

  ! -----------------------------------------------------------------------
  ! Checks that assumptions about the ELM scale Hierarchy required for
  ! use of ATS are satsified
  !
  subroutine EM_ATS_CheckHierarchy(this, filter, nclumps, ncolumns_all)
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

    ncolumns = filter(1)%num_hydrologyc
    if (ncolumns_all /= ncolumns) then
       write(iulog,*) "ERROR: ATS does not support non-hydrology columns. ncolumns_all = ", ncolumns_all, " num_hydrologyc = ", ncolumns
       call endrun("ATS does not support non-hydrology columns in spatially explicit mode. Check that the land use dataset has no urban, lake, or other non-natural land cover types (URBAN_REGION_ID == 0).")
    end if

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
  end subroutine EM_ATS_CheckHierarchy
  
  ! -----------------------------------------------------------------------
  ! Checks that the ATS and ELM mesh representations are consistent
  !
  subroutine EM_ATS_CheckMesh(this, grc_pp, col_pp)
    implicit none

    class(em_ats_type)                   :: this
    type(column_physical_properties)     , intent(in)    :: col_pp
    type(gridcell_physical_properties_type) , intent(in) :: grc_pp

    ! local variables
    integer :: ats_ncols_local
    integer :: ats_ncols_global
    integer :: ats_nlevgrnd
    integer :: i, j, ii

    real(c_double) :: ats_dzs(this%nlevgrnd)
    real(c_double) :: ats_areas(this%ncolumns)
    real(c_double) :: ats_lats(this%ncolumns)
    real(c_double) :: ats_lons(this%ncolumns)

    ! compare to mesh info
    call ats_get_mesh_info2(this%ats, ats_ncols_local, ats_ncols_global, ats_nlevgrnd, ats_dzs, ats_areas, ats_lats, ats_lons)

    ! assertions on shapes
    if (ats_ncols_local /= this%ncolumns) then
       call endrun("ATS local ncolumns does not match requested ncolumns")
    end if

    if (ats_nlevgrnd /= this%nlevgrnd) then
       print*, "ATS: NLEVGRND = ", ats_nlevgrnd
       print*, "ELM: NLEVGRND = ", this%nlevgrnd
       call endrun("ATS local cells in the vertical does not match ELM nlevgrnd")
    end if

    ! check dzs
    ii = this%col_filter(1)
    do j=1,this%nlevgrnd
       if (abs(ats_dzs(j) - col_pp%dz(ii,j)) > 1.e-10_r8) then
          call endrun("ATS dzs do not match ELM dzs.")
       end if
    end do

    ! check column areas
    do i=1,this%ncolumns
       ii = this%col_filter(i)
       if (abs(ats_areas(i) - grc_pp%area(col_pp%gridcell(ii))) > 1.e-10_r8) then
          write(iulog,*) "WARNING: ATS column areas do not match ELM grid cell areas -- perhaps incorrect ordering."
          ! call endrun("ATS column areas do not match ELM grid cell areas -- perhaps incorrect ordering.")
       end if
    end do

    ! check lat lon
    if (ats_lons(1) >= 0.) then
       do i=1,this%ncolumns
          if (abs(ats_lats(i) - 0.) > 1.e-8 .or. abs(ats_lons(i) - 0.) > 1.e-8) then
             call endrun("WARNING: ATS column lat/lon does not match ELM grid cell lat/lon -- perhaps incorrect ordering.")
          end if
       end do
    end if
    
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
  subroutine EM_ATS_GetFieldPtr(this, var_id, size, field)
    implicit none

    class(em_ats_type)                   :: this
    integer(c_int), intent(in) :: var_id
    integer, intent(in) :: size
    real(c_double), intent(inout), pointer :: field(:)

    ! local variables
    type(c_ptr) :: field_ptr
    
    
    field_ptr = ats_get_field_ptr(this%ats, var_id)
    call c_f_pointer(field_ptr, field, [size])
  end subroutine EM_ATS_GetFieldPtr

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

  ! -----------------------------------------------------------------------
  ! Special purpose -- get the water content fields
  !
  subroutine EM_ATS_GetField_WaterContent(this, col_pp, soilstate_vars, col_ws, soilhydrology_vars)
    use elm_varcon                 , only : denh2o

    implicit none
    class(em_ats_type)                   :: this
    type(column_physical_properties)     , intent(in)    :: col_pp
    type(soilstate_type)     , intent(inout)    :: soilstate_vars
    type(column_water_state) , intent(inout)    :: col_ws
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
  

    ! locals
    integer :: i, ii, j
    real(r8) :: h2osoi_liq_total
    integer :: lverbosity = 2
    
    ! Pressure -- note, we don't set soilp here, just soilpsi.  Is
    ! soilp used?  GDB suggests not...
    call EM_ATS_GetField_CopySubsurface(this, "pressure", ats_var_id%PRESSURE, &
         soilstate_vars%soilpsi_col)
    ! convert from Pa to MPa suction
    do j=1,this%nlevgrnd
       do i=1,this%ncolumns
          ii = this%col_filter(i)
          ! subtract ATS atmospheric pressure, convert to MPa, suction is always negative
          !
          ! note the ELM calculation of soilpsi relies on a sucsat to
          ! bound this away from 0 (see HydrologyNoDrainageMod.F90:509):
          ! ! use the same contants used in the supercool so that psi for frozen soils is consistent
          ! fsattmp = max(vwc/watsat(c,j), 0.001_r8)
          ! psi = sucsat(c,j) * (-9.8e-6_r8) * (fsattmp)**(-bsw(c,j))  ! Mpa
          ! soilpsi(c,j) = min(max(psi,-15.0_r8),0._r8)
          !
          ! Here we set an arbitrary min value of 100Pa, van Genuchten in general has no specific min val
          soilstate_vars%soilpsi_col(ii,j) = max(min(-1.e-4, (soilstate_vars%soilpsi_col(ii,j) - this%p_atm) * 1.e-6_r8), -15.0_r8)
          
       end do
    end do

    ! Water content
    call EM_ATS_GetField_CopySubsurface(this, "water_content", ats_var_id%WATER_CONTENT, &
         col_ws%h2osoi_liq(:,1:this%nlevgrnd))
    ! convert from volumetric water content to kg / m^2
    do i=1,this%ncolumns
       ii = this%col_filter(i)
       h2osoi_liq_total = 0._r8

       do j=1,this%nlevgrnd
          if (this%verbosity >= (lverbosity+1)) then
             write(iulog,*) "GET: volumetric_water_content [-] : column ", i, " cell ", j, " = ", col_ws%h2osoi_liq(ii,j)
          endif

          col_ws%h2osoi_liq(ii,j) = col_ws%h2osoi_liq(ii,j) * denh2o * col_pp%dz(ii,j)
          h2osoi_liq_total = h2osoi_liq_total + col_ws%h2osoi_liq(ii,j)

          if (this%verbosity >= (lverbosity+1)) then
             write(iulog,*) "ELM: dz : column ", i, " cell ", j, " = ", col_pp%dz(ii,j)
             write(iulog,*) "ELM: density = ", denh2o
             write(iulog,*) "GET: h2osoi_liq [kg m^-2] : column ", i, " cell ", j, " = ", col_ws%h2osoi_liq(ii,j)
          end if
       end do
       if (this%verbosity >= lverbosity) then
          write(iulog,*) "GET: h2osoi_liq total : column ", i, " total = ", h2osoi_liq_total
       end if
    end do

    ! depth to water table
    call EM_ATS_GetField_Copy(this, "depth_to_water_table", ats_var_id%DEPTH_TO_WATER_TABLE, soilhydrology_vars%zwt_col)

    ! ponded water
    call EM_ATS_GetField_Copy(this, "ponded depth", ats_var_id%PONDED_DEPTH, col_ws%h2osfc)
    ! convert from m to kg / m^2
    do i=1,this%ncolumns
       ii = this%col_filter(i)
       col_ws%h2osfc(ii) = col_ws%h2osfc(ii) * denh2o
       if (this%verbosity >= lverbosity) then
          write(iulog,*) "GET: h2osfc total : column ", i, " total = ", col_ws%h2osfc(ii)
       end if
    end do
  end subroutine EM_ATS_GetField_WaterContent


  ! -----------------------------------------------------------------------
  ! Special purpose -- set water content to the OLD timestep
  !
  subroutine EM_ATS_SetField_WaterContent(this, col_pp, col_ws, soilstate_vars)
    use elm_varcon                 , only : denh2o
    implicit none

    class(em_ats_type)                                  :: this
    type(column_physical_properties)    , intent(in)    :: col_pp
    type(column_water_state)            , intent(in)    :: col_ws
    type(soilstate_type)                , intent(in)    :: soilstate_vars

    ! local variables
    integer :: i, ii, j
    real(r8) :: h2osoi_liq_tot
    real(c_double), pointer :: ats_surf_wc(:)
    real(c_double), pointer :: ats_wc(:)
    integer :: lverbosity = 2    

    ! set base porosity to be the liquid porosity (e.g. ice = soil)
    call EM_ATS_SetField_CopySubsurface(this, "base porosity", ats_var_id%BASE_POROSITY, &
         soilstate_vars%eff_porosity_col)

    ! set surface water content
    call EM_ATS_GetFieldPtrW(this, ats_var_id%SURFACE_WATER_CONTENT_OLD, this%ncolumns, ats_surf_wc)
    do i=1,this%ncolumns
       ii = this%col_filter(i)
       ! h2osfc in mm --> m
       ats_surf_wc(i) = col_ws%h2osfc(ii) / denh2o
       if (this%verbosity >= lverbosity) then
          write(iulog,*) "SET: h2osfc total : column ", i, " total = ", col_ws%h2osfc(ii)
       end if
    end do

    ! set subsurface water content
    call EM_ATS_GetFieldPtrW(this, ats_var_id%WATER_CONTENT_OLD, this%ncolumns*this%nlevgrnd, ats_wc)
    do i=1,this%ncolumns
       ii = this%col_filter(i)
       h2osoi_liq_tot = 0._r8

       do j=1,this%nlevgrnd
          ! convert to volumetric water content [-]
          ats_wc((i-1) * this%nlevgrnd + j) = col_ws%h2osoi_liq(ii,j) / denh2o / col_pp%dz(ii,j)
          h2osoi_liq_tot = h2osoi_liq_tot + col_ws%h2osoi_liq(ii,j)
       end do

       if (this%verbosity >= lverbosity) then
          write(iulog,*) "SET: h2osoi_liq total : column ", i, " total = ", h2osoi_liq_tot
       end if
    end do

  end subroutine EM_ATS_SetField_WaterContent
  
  ! -----------------------------------------------------------------------
  ! Special purpose -- downregulate transpiration
  !
  subroutine EM_ATS_GetField_ActualTranspiration(this, col_wf, photosyns_vars)
    implicit none

    class(em_ats_type)                                  :: this
    type(column_water_flux)             , intent(inout) :: col_wf
    type(photosyns_type), intent(inout) :: photosyns_vars

    ! locals
    integer :: i,ii, pp
    real(r8) :: downreg, downreg_eps, tot_trans, diff
    real(c_double), pointer :: ats_tot_trans(:)
    integer :: lverbosity = 2
    
    ! get the total transpiration from ATS
    call EM_ATS_GetFieldPtr(this, ats_var_id%COLUMN_TRANSPIRATION, this%ncolumns, ats_tot_trans)

    do i=1,this%ncolumns
       ii = this%col_filter(i)
       pp = this%pft_filter(i)

       if (this%verbosity >= lverbosity) then
          write(iulog,*) "     pot T : column ", i, " total = ", col_wf%qflx_tran_veg(ii)
       end if

       
       downreg = 1._r8
       if (col_wf%qflx_tran_veg(ii) > 0.) then
          tot_trans = ats_tot_trans(i) * 1.e3_r8 ! m/s -> mm/s

          ! note, these are distinct -- ATS transpiration distribution
          ! has a nonlinear solve, which has its own tolerance.  When
          ! potential T is very small, the distribution solution may
          ! have error comparable to the potential T, and then the
          ! actual total_trans may have significant error in it.  This
          ! doesn't affect water balance because the absolute error is
          ! small even if the relative error is large.  So we compute
          ! two downregulation factors -- one to throw an error for
          ! debugging that adds an absolution portion, and one for
          ! downregulating carbon that is bounded below 1.
          downreg = min(tot_trans / col_wf%qflx_tran_veg(ii), 1.0_r8)
          downreg_eps = tot_trans / max(col_wf%qflx_tran_veg(ii), 1.e-10)
          diff = max(col_wf%qflx_tran_veg(ii) - tot_trans, 0._r8)

          if (downreg_eps > 1.01_r8 .OR. downreg < 0._r8) then
             write(iulog,*) "WARNING: ATS transpiration downregulation is out of expected bounds."
             call endrun("ATS transpiration downregulation is out of expected bounds.")
          end if

          ! evap_veg, evap_tot
          col_wf%qflx_evap_veg(ii) = col_wf%qflx_evap_veg(ii) - diff
          col_wf%qflx_evap_tot(ii) = col_wf%qflx_evap_tot(ii) - diff
          
          ! CO2 scales linearly with water? CHECK THIS! --ETC
          if (downreg < 1.0_r8) then
             photosyns_vars%fpsn_patch(pp) = photosyns_vars%fpsn_patch(pp) * downreg
             photosyns_vars%fpsn_wc_patch(pp) = photosyns_vars%fpsn_wc_patch(pp) * downreg
             photosyns_vars%fpsn_wj_patch(pp) = photosyns_vars%fpsn_wj_patch(pp) * downreg
             photosyns_vars%fpsn_wp_patch(pp) = photosyns_vars%fpsn_wp_patch(pp) * downreg
          end if

          ! reduce latent heat, increase sensible
          ! ... TODO: ETC

          ! correct pft_tran_veg
          col_wf%qflx_tran_veg(ii) = tot_trans
       else
          tot_trans = 0.0_r8
       end if

       if (this%verbosity >= lverbosity) then
          write(iulog,*) "     T actual : column ", i, " total = ", tot_trans
          write(iulog,*) "     T downreg : column ", i, " total = ", downreg
          write(iulog,*) "GET: actual transpiration : column ", i, " total = ", col_wf%qflx_tran_veg(ii)
       end if
    end do
  end subroutine EM_ATS_GetField_ActualTranspiration


  ! -----------------------------------------------------------------------
  ! Special purpose -- downregulate evaporation
  !
  subroutine EM_ATS_GetField_ActualEvaporation(this, col_wf)
    implicit none

    class(em_ats_type)                                  :: this
    type(column_water_flux)             , intent(in)    :: col_wf

    ! locals
    integer :: i,ii
    real(r8) :: downreg, evap, pot_evap, diff
    real(c_double), pointer :: ats_pot_evap(:), ats_evap(:)
    integer :: lverbosity = 2
    
    ! get the total evaporation from ATS
    call EM_ATS_GetFieldPtr(this, ats_var_id%EVAPORATION, this%ncolumns, ats_evap)
    call EM_ATS_GetFieldPtr(this, ats_var_id%POTENTIAL_EVAPORATION, this%ncolumns, ats_pot_evap)

    do i=1,this%ncolumns
       ii = this%col_filter(i)

       ! units: m/s --> mm/s
       evap = ats_evap(i) * 1.e3_r8
       pot_evap = ats_pot_evap(i) * 1.e3_r8
       diff = pot_evap - evap
       downreg = 1.0

       if (this%verbosity >= lverbosity) then
          write(iulog,*) "     pot evaporation : column ", i, " = ", pot_evap
          write(iulog,*) "     qflx_evap_tot : column ", i, " = ", col_wf%qflx_evap_tot(i)
       end if

       if (pot_evap > 0. .and. diff > 0.) then
          downreg = evap / pot_evap

          if (downreg > 1.0 .OR. downreg < 0.0) then
             write(iulog,*) "WARNING: ATS evaporation downregulation is out of expected bounds."
             call endrun("ATS evaporation downregulation is out of expected bounds.")
          end if


          ! no distinction between h2osfc and h2osoi evap
          !
          ! take it all from soil evap
          ! note these are all per unit column area
          col_wf%qflx_evap_soi(ii) = col_wf%qflx_evap_soi(ii) - diff
          col_wf%qflx_evap_tot(ii) = col_wf%qflx_evap_tot(ii) - diff
       end if

       if (this%verbosity >= lverbosity) then
          write(iulog,*) "     pot evaporation : column ", i, " = ", pot_evap
          write(iulog,*) "GET: actual evaporation : column ", i, " = ", evap
          write(iulog,*) "     E downreg : column ", i, " total = ", downreg
          write(iulog,*) "     E difference : column ", i, " = ", diff
       end if
    end do
  end subroutine EM_ATS_GetField_ActualEvaporation
  
  
#endif
  
end module ExternalModelATS
