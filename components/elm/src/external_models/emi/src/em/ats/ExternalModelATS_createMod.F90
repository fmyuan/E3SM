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
  use elm_varctl                   , only : use_ats

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
  subroutine EM_ATS_create(elmats_interface, iam, bounds_clump, col_num)

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
    integer                     , intent(in)    :: iam
    type(bounds_type)           , intent(in)    :: bounds_clump
    integer                     , intent(in)    :: col_num    ! actually elm's grid no.

    ! local

    !-----------------------------------------------------------------------
    !
    ! create an ATS driver object
    elmats_interface = ats_create(ats_inputdir, ats_inputfile, mpicom)


    if (use_ats) then
      ! pass mesh from ATS to ELM locally
      call get_mesh_local(elmats_interface, bounds_clump, col_num)  ! in progress .......

      ! 'mpicom' is communicator group id for land component
      print *, ''
      print *, '============================================================='
      print *,''
      print *, ' -------- ELM-ATS Coupled Mode ------------------------------'
      print *, ''
      print *, 'EM_ATS_Init: ats inputs - ', trim(ats_inputdir), ' ', trim(ats_inputfile)
      print *, 'communicator id: ', mpicom
      !

    end if

  end subroutine EM_ATS_create

  !------------------------------------------------------------------------


  subroutine get_mesh_local(elmats_interface, bounds_clump, col_num)
    !
    !DESCRIPTION
    !  mesh from ATS to ELM, locally
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    type(elm_ats_interface_type), intent(in) :: elmats_interface
    type(bounds_type)           , intent(in) :: bounds_clump
    integer                     , intent(in) :: col_num

    ! !LOCAL VARIABLES:
    integer                              :: c, fc, j

    integer                              :: ncols_local, ncols_global
    integer                              :: ncells_per_col
    real(r8)  , pointer                  :: surf_xi(:)
    real(r8)  , pointer                  :: surf_yi(:)
    real(r8)  , pointer                  :: surf_zi(:)
    real(r8)  , pointer                  :: surf_area(:)
    real(r8)  , pointer                  :: col_zi(:,:)
    real(r8)  , pointer                  :: col_dz(:,:)
    integer(C_INT), pointer              :: pft(:)

    !-----------------------------------------------------------------------

    allocate(surf_xi(1:col_num))           !
    allocate(surf_yi(1:col_num))           !
    allocate(surf_zi(1:col_num))           !
    allocate(surf_area(1:col_num))         !

    allocate(col_zi(1:col_num, 1:15))      ! col. cell bottom, assuming top-face of 1st cell is the surf (0.0m)
    allocate(col_dz(1:col_num, 1:15))      ! col. cell thickness (m)
    allocate(pft(1:col_num))

    ! in progress ...
    call elmats_interface%getmesh(ncols_local, ncols_global, ncells_per_col,  &
       surf_yi, surf_xi, surf_zi, surf_area, pft, col_zi)


print *, '----- checking mesh from ATS: '
print *, 'lat:', surf_yi
print *, 'lon:', surf_xi
print *, 'elev:', surf_zi
print *, 'depth:', col_zi
print *, 'surf_area:', surf_area
print *, 'pft index:', pft


    deallocate(surf_xi)
    deallocate(surf_yi)
    deallocate(surf_zi)
    deallocate(surf_area)
    deallocate(col_zi)
    deallocate(col_dz)
    deallocate(pft)

  end subroutine get_mesh_local


 !------------------------------------------------------------------------

#endif

end module ExternalModelATS_createMod
