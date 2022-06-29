module ELM_ATS_InterfaceMod

#ifdef USE_ATS_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! The ELM Interface to ATS Fortran interface,
  ! which corresponding to ATS's ats_elm_interface (extern "C" interface)
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
      procedure, public :: setup        => ats_setup_f
      procedure, public :: init         => ats_init_f
      procedure, public :: onestep      => ats_advance_f

      procedure, public :: setmesh      => ats_setmesh_f
      procedure, public :: setmat       => ats_setmats_f
      procedure, public :: setic        => ats_setIC_f
      procedure, public :: setbc        => ats_setBC_f
      procedure, public :: setss        => ats_setSS_f

      procedure, public :: getmesh      => ats_getmesh_f
      procedure, public :: getdata      => ats_getdata_f

  end type elm_ats_interface_type
  !------------------------------------------------------------------------
  !
  type(elm_ats_interface_type), public:: ats_drv

  interface ats_drv
      procedure ats_create
  end interface

  !------------------------------------------------------------------------
contains

  ! wrap the C++ functions/classes of ATS and data passing

  function ats_create()
      implicit none
      type(elm_ats_interface_type) :: ats_create
      ats_create%ptr = ats_elm_create()
  end function

 !------------------------------------------------------------------------

  subroutine ats_delete(this)
      implicit none
      type(elm_ats_interface_type) :: this
      call ats_elm_delete(this%ptr)
  end subroutine

 !------------------------------------------------------------------------

  subroutine ats_setup_f(this, input_dir, input_file, mpicomm)
      implicit none
      class(elm_ats_interface_type) :: this
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

      ! setup
      call ats_elm_setup(this%ptr, mpicomm, c_input_file) ! mpicomm is from 'spmdMod'

   end subroutine ats_setup_f

 !------------------------------------------------------------------------

   subroutine ats_init_f(this)

      implicit none
      class(elm_ats_interface_type) :: this

        ! initialize
        call ats_elm_initialize(this%ptr)

    end subroutine ats_init_f

 !------------------------------------------------------------------------
    subroutine ats_setmesh_f(this, surfgX, surfgY, surfgZ, soilcol_nodes)

        implicit none
        class(elm_ats_interface_type) :: this
        real(r8), pointer, intent(in) :: surfgX(:), surfgY(:), surfgZ(:,:)
        real(r8), pointer, intent(in) :: soilcol_nodes(:)

        ! local variables
        integer(C_INT) :: lx=2, ly=2, lz=16

        ! ----------------------------------------------------------
        lx = size(surfgX)
        ly = size(surfgY)
        lz = size(soilcol_nodes)
        ! elm surface-grids/soil column dimensions pass to ats
        call ats_elm_setmesh(this%ptr, surfgX, surfgY, surfgZ, soilcol_nodes, lx, ly, lz)

    end subroutine ats_setmesh_f

  !----------------------------------------------------------------------------------

    subroutine ats_setmats_f(this, ats_watsat, ats_hksat, ats_bsw, ats_sucsat, &
             ats_residual_sat, ats_eff_porosity, ats_initzwt)

        implicit none
        class(elm_ats_interface_type) :: this
        real(r8), pointer, intent(in) :: ats_watsat(:,:)
        real(r8), pointer, intent(in) :: ats_hksat (:,:)
        real(r8), pointer, intent(in) :: ats_bsw   (:,:)
        real(r8), pointer, intent(in) :: ats_sucsat(:,:)
        real(r8), pointer, intent(in) :: ats_residual_sat(:,:)
        real(r8), pointer, intent(in) :: ats_eff_porosity(:,:)
        real(r8), pointer, intent(in) :: ats_initzwt(:)
        !
        ! ----------------------------------------------------------
        !
        call ats_elm_setmats(this%ptr, ats_watsat, ats_hksat, ats_bsw, ats_sucsat, &
             ats_residual_sat, ats_eff_porosity, ats_initzwt)
        !
    end subroutine ats_setmats_f

 !------------------------------------------------------------------------

    subroutine ats_setIC_f(this, patm, soilp, wtd, starting_time, ats_visout)
        implicit none
        class(elm_ats_interface_type) :: this
        real(r8), pointer, intent(in) :: patm(:)
        real(r8), pointer, intent(in) :: soilp(:,:)
        real(r8), pointer, intent(in) :: wtd(:)
        real(r8), optional,intent(in) :: starting_time                  ! ELM starting time (in second, 0 by default)
        logical,  optional,intent(in) :: ats_visout                     ! instruct ATS output initial data
        !

        ! Local variables
        real(KIND=C_DOUBLE)  :: stime = 0.0
        logical(KIND=C_BOOL) :: visout = .true.                         ! a note here: fortran logical is an integer, but C/C++ bool is in byte. So this is a must
        ! ----------------------------------------------------------
        !
        if (present(starting_time)) stime = starting_time
        if (present(ats_visout)) visout = ats_visout

        !
        call ats_elm_setIC(this%ptr, stime, patm, soilp, wtd, visout)
        !
    end subroutine ats_setIC_f

 !------------------------------------------------------------------------

    subroutine ats_setBC_f(this, soilbot_flux, soildrain_flux)
        implicit none
        class(elm_ats_interface_type) :: this
        !
        real(r8), pointer, intent(in) :: soilbot_flux(:)
        real(r8), pointer, intent(in) :: soildrain_flux(:,:)

        ! Local variables

        ! ----------------------------------------------------------
        ! TODO
        call ats_elm_setBC(this%ptr)
        !
    end subroutine ats_setBC_f

 !------------------------------------------------------------------------

    subroutine ats_setSS_f(this, soilinfl_flux, soilevap_flux, pfttran_flux, soilbot_flux, &
                           roottran_flux, soildrain_flux)

        implicit none
        class(elm_ats_interface_type) :: this
        real(r8), pointer, intent(in) :: soilinfl_flux(:)       ! unit: kgH2O/m3/s
        real(r8), pointer, intent(in) :: soilevap_flux(:)
        real(r8), pointer, intent(in) :: pfttran_flux(:,:)      ! pft-level (root-fraction summed) transpiration [col, pft]
        real(r8), pointer, intent(in) :: soilbot_flux(:)
        real(r8), pointer, intent(in) :: roottran_flux(:,:)     ! root-level (pft-summed) transpiration [col, soil-layer]
        real(r8), pointer, intent(in) :: soildrain_flux(:,:)
        !

        ! Local variables
        integer(C_INT) :: ncols
        integer(C_INT) :: ncells

        ! ----------------------------------------------------------
        !
        ncols = size(soilevap_flux)
        ncells= size(roottran_flux)

        ! maybe add-up of roottran_flux + soildrain_flux later on

        call ats_elm_setSS(this%ptr, soilinfl_flux, soilevap_flux, &
          pfttran_flux, roottran_flux, soildrain_flux, &
          ncols, ncells)

    end subroutine ats_setSS_f
 !------------------------------------------------------------------------

    subroutine ats_advance_f(this, dt, ats_visout, ats_chkout)
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

        call ats_elm_advance(this%ptr, dt, visout, chkout)

    end subroutine ats_advance_f

 !------------------------------------------------------------------------

    subroutine ats_getmesh_f(this)
        implicit none
        class(elm_ats_interface_type) :: this
        real(r8), pointer :: dz(:)
        real(r8), pointer :: depth(:)
        real(r8), pointer :: elev(:)
        real(r8), pointer :: surf_area_m2(:)
        real(r8), pointer :: lat(:)
        real(r8), pointer :: lon(:)
        integer(C_INT) :: ncols_local
        integer(C_INT) :: ncols_global
        integer(C_INT) :: ncells_per_col

        ! ----------------------------------------------------------

        call ats_elm_getmesh(this%ptr, ncols_local, ncols_global, ncells_per_col, dz, depth, elev, &
          surf_area_m2, lat, lon)
        !

    end subroutine ats_getmesh_f

 !------------------------------------------------------------------------
    subroutine ats_getdata_f(this)
        implicit none
        class(elm_ats_interface_type) :: this

        real(r8), pointer :: surf_pres(:)
        real(r8), pointer :: soil_pres(:,:)
        real(r8), pointer :: satliq(:,:)
        integer(C_INT) :: ncols
        integer(C_INT) :: ncells

        ! ----------------------------------------------------------
        ! (TODO: not finished)
        ncols = 1 !size(surf_pres)
        ncells= 15 !size(soil_pres)
        !call ats_elm_getwaterstate(this%ptr, surf_pres, soil_pres, satliq, ncols, ncells)
        !

    end subroutine ats_getdata_f

 !------------------------------------------------------------------------


#endif

end module ELM_ATS_InterfaceMod
