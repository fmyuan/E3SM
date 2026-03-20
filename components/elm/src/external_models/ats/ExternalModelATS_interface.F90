module ExternalModelATS_interface
  use iso_c_binding
  implicit none

  ! ---------------------------------------------------------------------------
  ! > C function interfaces
  ! ---------------------------------------------------------------------------

  interface
    ! allocate, call constructor, and cast ptr to opaque ELM_ATSDriver_ptr
    function ats_create(f_comm, input_filename, logfile_filename) bind(C, name="ats_create")
      import :: c_ptr, c_char, c_int
      integer(c_int) :: f_comm         ! MPI_Fint is equivalent to int in C
      character(kind=c_char), intent(in) :: input_filename(*)  ! null-terminated string
      character(kind=c_char), intent(in) :: logfile_filename(*)  ! null-terminated string
      type(c_ptr) :: ats_create
    end function ats_create

    ! delete and deallocate
    subroutine ats_delete(ats) bind(C, name="ats_delete")
      import :: c_ptr
      type(c_ptr), value :: ats
    end subroutine ats_delete

    ! -------------------------------------------------------------------------
    ! Control routines
    ! -------------------------------------------------------------------------
    subroutine ats_parse_parameter_list(ats) bind(C, name="ats_parse_parameter_list")
      import :: c_ptr
      type(c_ptr), value :: ats
    end subroutine ats_parse_parameter_list

    subroutine ats_get_mesh_info(ats, ncols_local) bind(C, name="ats_get_mesh_info")
      import :: c_ptr, c_int
      type(c_ptr), value :: ats
      integer(c_int), intent(out) :: ncols_local
    end subroutine ats_get_mesh_info

    subroutine ats_get_mesh_info2(ats, ncols_local, ncols_global, nlevgrnd, dzs, areas, lat, lon) bind(C, name="ats_get_mesh_info2")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: ats
      integer(c_int), intent(out) :: ncols_local
      integer(c_int), intent(out) :: ncols_global
      integer(c_int), intent(out) :: nlevgrnd
      real(c_double), intent(inout) :: dzs(*)   ! array of size nlevgrnd
      real(c_double), intent(inout) :: areas(*) ! array of size ncols_local
      real(c_double), intent(inout) :: lat(*)   ! array of size ncols_local
      real(c_double), intent(inout) :: lon(*)   ! array of size ncols_local
    end subroutine ats_get_mesh_info2

    subroutine ats_setup(ats) bind(C, name="ats_setup")
      import :: c_ptr
      type(c_ptr), value :: ats
    end subroutine ats_setup

    subroutine ats_initialize(ats) bind(C, name="ats_initialize")
      import :: c_ptr
      type(c_ptr), value :: ats
    end subroutine ats_initialize

    subroutine ats_advance(ats, dt, checkpoint, visualize) bind(C, name="ats_advance")
      import :: c_ptr, c_double, c_bool
      type(c_ptr), value :: ats
      real(c_double), value :: dt
      logical(c_bool), value :: checkpoint
      logical(c_bool), value :: visualize
    end subroutine ats_advance

    subroutine ats_advance_test(ats) bind(C, name="ats_advance_test")
      import :: c_ptr
      type(c_ptr), value :: ats
    end subroutine ats_advance_test

    ! -------------------------------------------------------------------------
    ! Memory movement / data passing
    ! -------------------------------------------------------------------------
    subroutine ats_set_scalar(ats, var_id, in) bind(C, name="ats_set_scalar")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: ats
      integer(c_int), value :: var_id
      real(c_double), intent(in) :: in
    end subroutine ats_set_scalar

    subroutine ats_get_field(ats, var_id, out) bind(C, name="ats_get_field")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: ats
      integer(c_int), value :: var_id
      real(c_double), intent(out) :: out(*)
    end subroutine ats_get_field

    function ats_get_field_ptr(ats, var_id) result(ptr) bind(C, name="ats_get_field_ptr")
      import :: c_ptr, c_int
      type(c_ptr), value :: ats
      integer(c_int), value :: var_id
      type(c_ptr) :: ptr
    end function ats_get_field_ptr

    function ats_get_field_ptr_w(ats, var_id) result(ptr) bind(C, name="ats_get_field_ptr_w")
      import :: c_ptr, c_int
      type(c_ptr), value :: ats
      integer(c_int), value :: var_id
      type(c_ptr) :: ptr
    end function ats_get_field_ptr_w

    subroutine ats_set_field(ats, var_id, in) bind(C, name="ats_set_field")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: ats
      integer(c_int), value :: var_id
      real(c_double), intent(in) :: in(*)
    end subroutine ats_set_field

  end interface

end module ExternalModelATS_interface
