module ExternalModelATS_readnlMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! NOTE - this is separated from rest of ExternalModelATS modules so that
  !        avoid cycled calling when initializing ELM.

  ! ELM module use
  use abortutils                   , only : endrun
  use shr_log_mod                  , only : errMsg => shr_log_errMsg
  use spmdMod                      , only : masterproc, mpicom
  use elm_varctl                   , only : iulog
  use elm_varctl                   , only : use_ats
  use elm_varctl                   , only : ats_hmode, ats_thmode, ats_thcmode, ats_gmode
  !
  implicit none
  !
  ! read elm-ats namelist
  character(len=256), public:: ats_inputdir = ''
  character(len=256), public:: ats_inputfile = ''
  public :: elm_ats_readnl

  !---------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
   subroutine elm_ats_readnl( NLFilename )
  !
  ! !DESCRIPTION:
  ! Read namelist for elm-ats interface
  !
  ! !USES:
    use fileutils     , only : getavu, relavu, opnfil
    use elm_nlUtilsMod, only : find_nlgroup_name
    use shr_nl_mod    , only : shr_nl_find_group_name
    use shr_mpi_mod   , only : shr_mpi_bcast

    implicit none

  ! !ARGUMENTS:
    character(len=*), intent(in) :: NLFilename    ! Namelist ats input filename (*.xml)
  ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=32) :: subname = 'elm_ats_readnl'  ! subroutine name
  !
  !-----------------------------------------------------------------------
    namelist / elm_ats_inparm / ats_inputdir, ats_inputfile

    ! ----------------------------------------------------------------------
    ! Read namelist from standard namelist file.
    ! ----------------------------------------------------------------------

    if ( masterproc )then

       unitn = getavu()
       write(iulog,*) 'Read in elm-ats namelist'
       open (unitn, file=trim(NLFilename), status='old', iostat=ierr)
       call shr_nl_find_group_name(unitn, 'elm_ats_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, elm_ats_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg=subname //':: ERROR: reading elm_ats_inparm namelist.'//&
                         errMsg(__FILE__, __LINE__))
          end if
       end if
       close (unitn)
       call relavu(unitn)
       write(iulog, '(/, A)') " elm-ats namelist:"
       write(iulog, '(A, " : ", A,/)') "   ats_inputdir", trim(ats_inputdir)
       write(iulog, '(A, " : ", A,/)') "   ats_inputfile ", trim(ats_inputfile)
    end if

    ! Broadcast namelist variables read in
    call shr_mpi_bcast(ats_inputdir, mpicom)
    call shr_mpi_bcast(ats_inputfile, mpicom)

    ! by default, ats subsurface hydrology is on, if use_ats = .true.
    if (use_ats) ats_hmode = .true.
    call shr_mpi_bcast(ats_hmode, mpicom)
    call shr_mpi_bcast(ats_thmode, mpicom)
    call shr_mpi_bcast(ats_thcmode, mpicom)
    call shr_mpi_bcast(ats_gmode, mpicom)

  end subroutine elm_ats_readnl

end module ExternalModelATS_readnlMod
