module elm_finalizeMod

  !-----------------------------------------------------------------------
  ! Performs land model cleanup
  !
  !
  use SoilLittVertTranspMod, only : cleanupLitterTransportList

  implicit none
  save
  public   ! By default everything is public
  !
  !-----------------------------------------
  ! Instances of component types
  !-----------------------------------------
  !
  public :: final
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine final( )
    !
    ! !DESCRIPTION:
    ! Finalize land surface model
    !
#ifdef HAVE_MOAB
    use MOABGridType, only : elm_moab_finalize
#endif
#ifdef USE_PETSC_LIB
#include <petsc/finclude/petsc.h>
#endif
    ! !USES:
    use elm_varctl             , only : use_vsfm, use_cn, use_ats, use_ats_ic
#ifdef USE_PETSC_LIB
    use petscsys
#endif

#ifdef USE_ATS_LIB
    use ExternalModelATS       , only : em_ats
#endif

    !
    ! !ARGUMENTS
    implicit none
    !

#ifdef HAVE_MOAB
    call elm_moab_finalize()
#endif

#ifdef USE_PETSC_LIB
    PetscErrorCode        :: ierr

    if (use_vsfm) then
       call PetscFinalize(ierr)
    endif
#endif
    if (use_cn) then
       call cleanupLitterTransportList()
    endif

#ifdef USE_ATS_LIB
    if (use_ats .or. use_ats_ic) then
       call em_ats%Finalize()
    endif
#endif
    
  end subroutine final

end module elm_finalizeMod
