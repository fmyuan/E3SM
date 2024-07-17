module ELM_ATS_InterfaceMod


  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! The ELM Interface to ATS Fortran interface,
  ! which corresponding to ATS's ats_interface (extern "C" interface)
  !
  use iso_c_binding
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  !
  !
  !------------------------------------------------------------------------
  ! c++-fortran interface
  include "ELM_ATS_Interface_c2f.inc"

  ! c++-fortran interface
  !------------------------------------------------------------------------
  !
  type, public :: elm_ats_interface_type
     private
     type(C_PTR) :: ptr                ! pointer to ats driver

     ! ats mesh information
     integer(C_INT), public :: ncols_local
     integer(C_INT), public :: ncols_global
     integer(C_INT), public :: ncells_per_col
     real(r8), pointer, public :: lat(:)
     real(r8), pointer, public :: lon(:)
     real(r8), pointer, public :: elev(:)
     real(r8), pointer, public :: surf_area(:)
     integer(C_INT), pointer, public :: pft(:)
     real(r8), pointer, public :: depth(:,:)

   contains

     final :: ats_delete
     procedure, public :: setup         => ats_setup
     procedure, public :: getmesh       => ats_getmesh
     procedure, public :: setsoilveg    => ats_setsoilveg_parameters
     procedure, public :: init          => ats_init
     procedure, public :: onestep       => ats_advance
     procedure, public :: setsoilveg_dyn=> ats_setsoilveg_properties
     procedure, public :: setss_hydro   => ats_setss_hydro

     procedure, public :: getdata_hydro => ats_getdata_hydro

  end type elm_ats_interface_type
  !------------------------------------------------------------------------
  !
  type(elm_ats_interface_type), public:: ats_drv

  interface ats_drv
     procedure ats_create
  end interface ats_drv

  !------------------------------------------------------------------------
contains

  ! wrap the C++ functions/classes of ATS and data passing
  function ats_create(input_dir, input_file, mpicomm)
    implicit none
    type(elm_ats_interface_type) :: ats_create

    character(len=*), intent(in) :: input_dir
    character(len=*), intent(in) :: input_file
    integer, intent(in) :: mpicomm  ! mpi communicator group id (i.e. MPI_COMM_WORLD for lnd model???)

    ! local variables
    character(kind=C_CHAR) :: c_input_file(len_trim(input_dir)+len_trim(input_file)+2)
    integer :: i, n1, n2
    integer :: ierr

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

    ats_create%ptr = ats_create_c(mpicomm, c_input_file)
  end function ats_create

  !------------------------------------------------------------------------

  subroutine ats_delete(this)
    implicit none
    type(elm_ats_interface_type) :: this

    call ats_delete_c(this%ptr)

    deallocate(this%lat)
    deallocate(this%lon)
    deallocate(this%elev)
    deallocate(this%surf_area)
    deallocate(this%pft)
    deallocate(this%depth)

  end subroutine ats_delete

  !------------------------------------------------------------------------

  subroutine ats_setup(this)
    implicit none
    class(elm_ats_interface_type) :: this
    call ats_setup_c(this%ptr)
  end subroutine ats_setup

   !------------------------------------------------------------------------

  subroutine ats_getmesh(this)
    implicit none
    class(elm_ats_interface_type) :: this

    call ats_get_mesh_info_c(this%ptr, this%ncols_local, this%ncols_global, &
         this%lat, this%lon, this%elev, this%surf_area, &
         this%pft, this%ncells_per_col, this%depth)

    !
  end subroutine ats_getmesh

  !------------------------------------------------------------------------

  subroutine ats_setsoilveg_parameters(this, porosity, hksat, CH_bsw, CH_sucsat, &
       CH_residual_sat)
    implicit none
    class(elm_ats_interface_type) :: this
    real(r8), pointer, intent(in) :: porosity(:,:)
    real(r8), pointer, intent(in) :: hksat (:,:)
    real(r8), pointer, intent(in) :: CH_bsw   (:,:)
    real(r8), pointer, intent(in) :: CH_sucsat(:,:)
    real(r8), pointer, intent(in) :: CH_residual_sat(:,:)

    call ats_set_soil_hydrologic_parameters_c(this%ptr, porosity, hksat, CH_bsw, CH_sucsat, &
         CH_residual_sat)

  end subroutine ats_setsoilveg_parameters

  !----------------------------------------------------------------------------------

  subroutine ats_init(this, starting_time, soil_water_content, soil_pressure)
    implicit none
    class(elm_ats_interface_type) :: this

    real(r8), optional,intent(in) :: starting_time                  ! ELM starting time (in second, 0 by default)
    real(r8), pointer, intent(in) :: soil_water_content(:,:)
    real(r8), pointer, intent(in) :: soil_pressure(:,:)
    call ats_initialize_c(this%ptr, starting_time, soil_water_content, soil_pressure)
  end subroutine ats_init

  !------------------------------------------------------------------------

  subroutine ats_setss_hydro(this, soilinfl_flux, soilevap_flux, pfttran_flux)
    implicit none
    class(elm_ats_interface_type) :: this
    real(r8), pointer, intent(in) :: soilinfl_flux(:)       ! unit: kgH2O/m3/s
    real(r8), pointer, intent(in) :: soilevap_flux(:)
    real(r8), pointer, intent(in) :: pfttran_flux(:)      ! col-level? pft-level (root-fraction summed) transpiration [col, pft]
    call ats_set_potential_sources_c(this%ptr, soilinfl_flux, soilevap_flux, pfttran_flux)
  end subroutine ats_setss_hydro

  !------------------------------------------------------------------------

  subroutine ats_setsoilveg_properties(this, eff_porosity, dyn_rootfrac)
    implicit none
    class(elm_ats_interface_type) :: this
    real(r8), pointer, intent(in) :: eff_porosity(:,:)
    real(r8), pointer, intent(in) :: dyn_rootfrac(:,:)
    
    ! soil effective porosity
    call ats_set_soil_hydrologic_properties_c(this%ptr, eff_porosity)

    ! veg rooting fraction
    call ats_set_veg_properties_c(this%ptr, dyn_rootfrac)

  end subroutine ats_setsoilveg_properties

  !------------------------------------------------------------------------
  
  subroutine ats_advance(this, dt, ats_visout, ats_chkout)
    implicit none
    class(elm_ats_interface_type) :: this
    real(r8), intent(in) :: dt                                      ! one ELM timestep interval (in seconds)
    logical, optional, intent(in) :: ats_visout                     ! instruct ATS output data
    logical, optional, intent(in) :: ats_chkout                     ! instruct ATS output checkpoint
    !
    ! local variables
    logical(KIND=C_BOOL) :: visout = .true.
    logical(KIND=C_BOOL) :: chkout = .false.
    ! ----------------------------------------------------------
    !

    if (present(ats_visout)) visout = ats_visout
    if (present(ats_chkout)) chkout = ats_chkout

    call ats_advance_c(this%ptr, dt, visout, chkout)
  end subroutine ats_advance

  !------------------------------------------------------------------------
  subroutine ats_getdata_hydro(this, h2osfc, zwt, h2oliq, h2oice, soilpsi, &
                               soilinfl_flux, evap_flux, tran_flux, root_flux, netsub_flux, netsrf_runon)
    implicit none
    class(elm_ats_interface_type) :: this

    real(r8), pointer, intent(out) :: h2osfc(:)       ! mm H2O
    real(r8), pointer, intent(out) :: zwt(:)          ! water table depth (m below surface - positive downward)
    real(r8), pointer, intent(out) :: h2oliq(:,:)     ! kgH2O/m2
    real(r8), pointer, intent(out) :: h2oice(:,:)     ! kgH2O/m2
    real(r8), pointer, intent(out) :: soilpsi(:,:)    ! non-negative Pa
    real(r8), pointer, intent(out) :: soilinfl_flux(:)       ! unit: mm/s
    real(r8), pointer, intent(out) :: evap_flux(:)           ! unit: mm/s
    real(r8), pointer, intent(out) :: tran_flux(:)           ! unit: mm/s - transpiration rate at leaves [col]
    real(r8), pointer, intent(out) :: root_flux(:, :)        ! unit: mm/s - transpiration rate at roots [col, nlevgrnd] (TODO)
    real(r8), pointer, intent(out) :: netsub_flux(:)         ! unit: mm/s
    real(r8), pointer, intent(out) :: netsrf_runon(:)        ! unit: mm/s

    ! note: zwt NOT really what is now (TODO)
    call ats_get_waterstate_c(this%ptr, h2osfc, zwt, h2oliq, h2oice, soilpsi)

    ! not yet (TODO)
    call ats_get_water_fluxes_c(this%ptr, soilinfl_flux, evap_flux, tran_flux,  &
         root_flux, netsub_flux, netsrf_runon)

    !
  end subroutine ats_getdata_hydro

  !------------------------------------------------------------------------


end module ELM_ATS_InterfaceMod
