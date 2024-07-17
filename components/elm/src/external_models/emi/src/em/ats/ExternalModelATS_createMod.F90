module ExternalModelATS_createMod

#ifdef USE_ATS_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !

  ! ELM modules
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use spmdMod                      , only : masterproc, mpicom
  use decompMod                    , only : bounds_type
  ! a few ats coupling options
  use elm_varctl                   , only : use_ats, use_ats_mesh

  !
  use ExternalModelATS_readnlMod   , only : ats_inputdir, ats_inputfile

  ! C-F interface
  use ELM_ATS_InterfaceMod

  !
  implicit none
  !
  !---------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine EM_ATS_create(elmats_interface)

    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:

    !
    implicit none
    !
    ! !ARGUMENTS:
    type(elm_ats_interface_type), intent(inout) :: elmats_interface

    ! local

    !-----------------------------------------------------------------------
    !
    ! create an ATS driver object
    elmats_interface = ats_create(ats_inputdir, ats_inputfile, mpicom)

    ! 'mpicom' is communicator group id for land component
    print *, ''
    print *, '============================================================='
    print *,''
    print *, ' -------- ELM-ATS Coupled Mode ------------------------------'
    print *, ''
    print *, 'EM_ATS_Init: ats inputs - ', trim(ats_inputdir), ' ', trim(ats_inputfile)
    print *, 'communicator id: ', mpicom
    print *, '============================================================='
    print *, ''
    !


    if (use_ats .and. use_ats_mesh) then
      ! pass mesh from ATS to ELM locally
      call get_mesh_local(elmats_interface)  ! in progress .......

    end if

  end subroutine EM_ATS_create

  !------------------------------------------------------------------------


  subroutine get_mesh_local(elmats_interface)
    !
    !DESCRIPTION
    !  mesh from ATS to ELM, locally
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    type(elm_ats_interface_type), intent(inout) :: elmats_interface

    ! !LOCAL VARIABLES:
    integer                              :: c, fc, j
    integer                              :: col_num, nz

    !-----------------------------------------------------------------------
    col_num = elmats_interface%ncols_local
    nz      = elmats_interface%ncells_per_col
    allocate(elmats_interface%lon(1:col_num))              !
    allocate(elmats_interface%lat(1:col_num))              !
    allocate(elmats_interface%elev(1:col_num))             !
    allocate(elmats_interface%surf_area(1:col_num))        !

    allocate(elmats_interface%depth(1:col_num, 1:nz))      ! col. cell bottom, assuming top-face of 1st cell is the surf (0.0m)
    allocate(elmats_interface%pft(1:col_num))

    ! in progress ...
    call elmats_interface%getmesh()


print *, '----- checking mesh from ATS: '
print *, 'lat:', elmats_interface%lat
print *, 'lon:', elmats_interface%lon
print *, 'elev:', elmats_interface%elev
print *, 'depth:', elmats_interface%depth
print *, 'surf_area:', elmats_interface%surf_area
print *, 'pft index:', elmats_interface%pft

  end subroutine get_mesh_local


 !------------------------------------------------------------------------

#endif

end module ExternalModelATS_createMod
