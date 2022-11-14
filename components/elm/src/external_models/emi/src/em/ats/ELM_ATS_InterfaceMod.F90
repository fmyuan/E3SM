module ELM_ATS_InterfaceMod

#ifdef USE_ATS_LIB

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

    ats_create%ptr = ats_create_f(mpicomm, c_input_file)
  end function ats_create

  !------------------------------------------------------------------------

  subroutine ats_delete(this)
    implicit none
    type(elm_ats_interface_type) :: this
    call ats_delete_f(this%ptr)
  end subroutine ats_delete

  !------------------------------------------------------------------------

  subroutine ats_setup(this)
    implicit none
    class(elm_ats_interface_type) :: this

    ! setup
    call ats_setup_f(this%ptr)
  end subroutine ats_setup

   !------------------------------------------------------------------------

  subroutine ats_getmesh(this, ncols_local, ncols_global, ncells_per_col, &
       lat, lon, elev, surf_area, depth)
    implicit none
    class(elm_ats_interface_type) :: this
    integer(C_INT) :: ncols_local
    integer(C_INT) :: ncols_global
    integer(C_INT) :: ncells_per_col

    real(r8), pointer :: lat(:)
    real(r8), pointer :: lon(:)
    real(r8), pointer :: elev(:)
    real(r8), pointer :: surf_area(:)

    real(r8), pointer :: depth(:,:)

    ! ---------------------------------------------------------------------

    call ats_get_mesh_info_f(this%ptr, ncols_local, ncols_global, &
         lat, lon, elev, surf_area, ncells_per_col, depth)

    !
  end subroutine ats_getmesh

  !------------------------------------------------------------------------

  subroutine ats_setsoilveg_parameters(this, porosity, hksat, CH_bsw, CH_sucsat, &
       CH_residual_sat, pft_full_turgor, pft_wilting_point)
    implicit none
    class(elm_ats_interface_type) :: this
    real(r8), pointer, intent(in) :: porosity(:,:)
    real(r8), pointer, intent(in) :: hksat (:,:)
    real(r8), pointer, intent(in) :: CH_bsw   (:,:)
    real(r8), pointer, intent(in) :: CH_sucsat(:,:)
    real(r8), pointer, intent(in) :: CH_residual_sat(:,:)

    real(r8), pointer, intent(in) :: pft_full_turgor(:)
    real(r8), pointer, intent(in) :: pft_wilting_point(:)

    ! soils
    call ats_set_soil_hydrologic_parameters_f(this%ptr, porosity, hksat, CH_bsw, CH_sucsat, &
         CH_residual_sat)

    ! veg
    call ats_set_veg_parameters_f(this%ptr, pft_full_turgor, pft_wilting_point)


  end subroutine ats_setsoilveg_parameters

  !----------------------------------------------------------------------------------

  subroutine ats_init(this, starting_time, patm, soilp)
    implicit none
    class(elm_ats_interface_type) :: this

    real(r8), optional,intent(in) :: starting_time                  ! ELM starting time (in second, 0 by default)
    real(r8), pointer, intent(in) :: patm(:)
    real(r8), pointer, intent(in) :: soilp(:,:)
    call ats_initialize_f(this%ptr, starting_time, patm, soilp)
  end subroutine ats_init

  !------------------------------------------------------------------------

  subroutine ats_setss_hydro(this, soilinfl_flux, soilevap_flux, pfttran_flux)
    implicit none
    class(elm_ats_interface_type) :: this
    real(r8), pointer, intent(in) :: soilinfl_flux(:)       ! unit: kgH2O/m3/s
    real(r8), pointer, intent(in) :: soilevap_flux(:)
    real(r8), pointer, intent(in) :: pfttran_flux(:,:)      ! pft-level (root-fraction summed) transpiration [col, pft]
    call ats_set_potential_sources_f(this%ptr, soilinfl_flux, soilevap_flux, pfttran_flux)
  end subroutine ats_setss_hydro

  !------------------------------------------------------------------------

  subroutine ats_setsoilveg_properties(this, eff_porosity, dyn_rootfrac)
    implicit none
    class(elm_ats_interface_type) :: this
    real(r8), pointer, intent(in) :: eff_porosity(:,:)

    real(r8), pointer, intent(in) :: dyn_rootfrac(:,:,:)

    ! soil effective porosity
    call ats_set_soil_hydrologic_properties_f(this%ptr, eff_porosity)

    ! veg rooting fraction
    call ats_set_veg_properties_f(this%ptr, dyn_rootfrac)


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

    call ats_advance_f(this%ptr, dt, visout, chkout)
  end subroutine ats_advance

  !------------------------------------------------------------------------
  subroutine ats_getdata_hydro(this, zwt, soilp, psi, h2oliq, h2oice, &
                               soilinfl_flux, evap_flux, pfttran_flux, netsub_flux, netsrf_runon)
    implicit none
    class(elm_ats_interface_type) :: this

    real(r8), pointer, intent(out) :: zwt(:)          ! meters
    real(r8), pointer, intent(out) :: soilp(:,:)      ! Pa (atm-pressue additive)
    real(r8), pointer, intent(out) :: psi(:,:)        ! non-negative Pa
    real(r8), pointer, intent(out) :: h2oliq(:,:)     ! kgH2O/m3
    real(r8), pointer, intent(out) :: h2oice(:,:)     ! kgH2O/m3

    real(r8), pointer, intent(out) :: soilinfl_flux(:)       ! unit: mm/s
    real(r8), pointer, intent(out) :: evap_flux(:)           ! unit: mm/s
    real(r8), pointer, intent(out) :: pfttran_flux(:, :, :)  ! pft-level (root-fraction summed) transpiration [col, nlevgrnd, pft]
    real(r8), pointer, intent(out) :: netsub_flux(:)         ! unit: mm/s
    real(r8), pointer, intent(out) :: netsrf_runon(:)        ! unit: mm/s

    ! note: zwt NOT really what is now (TODO)
    call ats_get_waterstate_f(this%ptr, zwt, soilp, psi, h2oliq, h2oice)


    ! not yet (TODO)
    call ats_get_water_fluxes_f(this%ptr, soilinfl_flux, evap_flux, pfttran_flux, &
         netsub_flux, netsrf_runon)

    !
  end subroutine ats_getdata_hydro

  !------------------------------------------------------------------------


#endif

end module ELM_ATS_InterfaceMod
