module EnhancedWeatheringMod

  !-----------------------------------------------------------------------
  ! !MODULE: EnhancedWeatheringMod
  !
  ! !DESCRIPTION:
  ! Module for rock powder dynamics (application, reaction, leaching)
  ! for coupled carbon-nitrogen(-phosphorus) code.
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use elm_varctl          , only : iulog, year_start_erw, nyear_erw_calibrate, mixing_layer
  use elm_varcon          , only : log_keq_co3, log_keq_hco3, log_keq_sio2am
  use elm_varcon          , only : mass_co3, mass_hco3, mass_co2, mass_h2o, mass_sio2, mass_h
  use elm_varcon          , only : passivation_phi, passivation_tau
  use elm_varcon          , only : zisoi
  use elm_varpar          , only : mixing_depth, nlevgrnd, nlevsoi
  use elm_varpar          , only : nminerals, ncations, nminsecs, nks
  use decompMod           , only : bounds_type
  use ColumnDataType      , only : col_ew, col_ms, col_mf, col_es, col_ws, col_wf
  use ColumnType          , only : col_pp
  use LandunitType        , only : lun_pp
  use TopounitDataType    , only : top_as
  use SoilStateType       , only : soilstate_type
  use ewutils             , only : logmol_to_mass, mol_to_mass, meq_to_mass, mass_to_mol, mass_to_meq, advection_diffusion
  use ewutils             , only : mass_to_logmol, find_net_charge, solve_eq, ph_to_hco3, hco3_to_co3
  use domainMod           , only : ldomain ! debug print
  use shr_sys_mod         , only : shr_sys_flush
  use spmdMod             , only : iam

  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: elm_erw_readnl
  public :: readEnhancedWeatheringParams
  public :: MineralInit
  public :: MineralBackground
  public :: MineralPrimary
  public :: MineralSecondary
  public :: MineralPassivation
  public :: MineralVerticalMovement
  public :: MineralLeaching
  public :: MineralEquilibria

  type, public :: EWParamsType
     character(len=40), pointer  :: minerals_name      (:)      => null()
     real(r8), pointer  :: primary_mass                (:)      => null()   ! molar mass of the primary mineral, g/mol, 1:nminerals (e.g. Mg2SiO4 = 140.6931 g/mol)
     real(r8), pointer  :: log_keq_primary             (:)      => null()   ! log10 of equilibrium constants for primary mineral dissolution
     real(r8), pointer  :: log_k_primary               (:, :)   => null()   ! log10 of primary mineral reaction rate constant at 298.15K (mol m-2 mineral surface area s-1), 1:nminerals x [H+, H2O, OH-]
     real(r8), pointer  :: e_primary                   (:, :)   => null()   ! primary mineral reaction activation energy constant (KJ mol-1), 1:nminerals x [H+, H2O, OH-]
     real(r8), pointer  :: n_primary                   (:, :)   => null()   ! reaction order of H+ and OH- catalyzed weathering, 1:nminerals x [H+, OH-]

     ! reaction stoichiometry: suppose the equation is 
     ! primary mineral + proton + (water) = cations + SiO2 + (water)
     ! coefficient before the mineral is always 1

     real(r8), pointer  :: primary_stoi_proton        (:)      => null()   ! reaction stoichiometry coefficient in front of H+, 1:nminerals
     real(r8), pointer  :: primary_stoi_cations       (:, :)   => null()   ! reaction stoichiometry coefficient in front of cations, 1:nminerals x 1:ncations
     real(r8), pointer  :: primary_stoi_h2o           (:)      => null()   ! reaction stoichiometry coefficient in front of water (posive = produced, negative = consumed), 1:nminerals
     real(r8), pointer  :: primary_stoi_sio2          (:)      => null()   ! reaction stoichiometry coefficient in front of SiO2, 1:nminerals
     real(r8), pointer  :: primary_stoi_hco3          (:)      => null()   ! reaction stoichiometry coefficient in front of HCO3-, 1:nminerals

     character(len=40), pointer  :: minsecs_name      (:)      => null()   ! names of the secondary minerals for the record
     real(r8), pointer  :: minsecs_mass               (:)      => null()   ! molar mass of the secondary mineral, g/mol, 1:nminsecss (e.g. Mg2SiO4 = 140.6931 g/mol)
     real(r8), pointer  :: log_keq_minsecs            (:)      => null()   ! log10 of equilibrium constants for secondary mineral precipitation
     real(r8), pointer  :: ssa_minsecs                (:)      => null()   ! alpha constants for secondary mineral precipitation
     real(r8), pointer  :: k_precip_minsecs           (:, :)   => null()  ! precipitation rate constant for secondary minerals, 1:nminsecs x [H+, H2O, OH-/HCO3-]
     real(r8), pointer  :: e_precip_minsecs           (:, :)   => null()  ! activation energy for secondary mineral precipitation, 1:nminsecs x [H+, H2O, OH-/HCO3-]
     real(r8), pointer  :: ph2o_precip_minsecs        (:)      => null()  ! p power on top of H2O of the precipitation rate law, 1:nminsecs
     real(r8), pointer  :: qh2o_precip_minsecs        (:)      => null()  ! q power on top of H2O of the precipitation rate law, 1:nminsecs
     real(r8), pointer  :: n_precip_minsecs           (:)      => null()  ! n power on top of H+ or HCO3- of the precipitation rate law, 1:nminsecs
     real(r8), pointer  :: k_dissolv_minsecs          (:, :)   => null()  ! dissolution rate constant for secondary minerals, 1:nminsecs x [H+, H2O, OH-/HCO3-]
     real(r8), pointer  :: e_dissolv_minsecs          (:, :)   => null()  ! activation energy for secondary mineral dissolution, 1:nminsecs x [H+, H2O, OH-/HCO3-]
     real(r8), pointer  :: n_dissolv_minsecs          (:, :)   => null()  ! reaction order of H+, H2O, and OH-/HCO3- catalyzed dissolution, 1:nminsecs x [H+, H2O, OH-/HCO3-]

     character(len=40), pointer  :: cations_name       (:)      => null()
     real(r8), pointer  :: cations_mass                (:)      => null()   ! molar masses of the cation species, g/mol
     real(r8), pointer  :: cations_valence             (:)      => null()   ! valence of the cations
     real(r8), pointer  :: cations_diffusivity         (:)      => null()   ! diffusion coefficient of the cations in water, m2/s

     real(r8), pointer  :: bicarbonate_diffusivity              => null()   ! diffusion coefficient of HCO3- in water, m2/s
     real(r8), pointer  :: carbonate_diffusivity                => null()   ! diffusion coefficient of CO3-- in water, m2/s

  end type EWParamsType

  type(EWParamsType), public ::  EWParamsInst
  !$acc declare create(EWParamsInst)

contains

  !-----------------------------------------------------------------------
  !
  ! !IROUTINE: elm_erw_readnl
  !
  ! !INTERFACE:
  subroutine elm_erw_readnl( NLFilename )
  !
  ! !DESCRIPTION:
  ! Read namelist for elm-pflotran interface
  !
  ! !USES:
    use elm_varctl    , only : iulog
    use elm_varctl    , only : year_start_erw, elm_erw_paramfile, use_erw_verbose, builtin_site
    use spmdMod       , only : masterproc, mpicom, MPI_CHARACTER, MPI_INTEGER
    use shr_log_mod   , only : errMsg => shr_log_errMsg
    use fileutils     , only : getavu, relavu, opnfil
    use abortutils    , only : endrun
    use elm_nlUtilsMod, only : find_nlgroup_name
    use shr_nl_mod    , only : shr_nl_find_group_name
    use shr_mpi_mod   , only : shr_mpi_bcast

    implicit none

  ! !ARGUMENTS:
    character(len=*), intent(IN) :: NLFilename ! Namelist filename

  ! !LOCAL VARIABLES:
    character(len=256) :: locfn     ! local file name
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=256):: errline
    character(len=32) :: subname = 'elm_erw_readnl'  ! subroutine name
  !EOP
  !-----------------------------------------------------------------------
    namelist / elm_erw_inparm / year_start_erw, nyear_erw_calibrate, elm_erw_paramfile, use_erw_verbose, builtin_site, mixing_layer

    ! ----------------------------------------------------------------------
    ! Read namelist from standard namelist file.
    ! ----------------------------------------------------------------------

    if ( masterproc)then

       unitn = getavu()
       write(iulog,*) 'Read in elm-erw namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'elm_erw_inparm', status=ierr)
       if ( ierr == 0 ) then
         read(unitn, elm_erw_inparm, iostat=ierr)
         if (ierr /= 0) then
           ! get the error line of namelist
           backspace(unitn)
           read(unitn,fmt='(A)') errline
           print *, 'Invalid line: ', trim(errline), ' in namelist file: ', trim(NLFilename)

           call endrun(msg=subname //':: ERROR: reading elm_erw_inparm namelist.'//&
                       errMsg(__FILE__, __LINE__))
         end if
       end if
       call relavu( unitn )
       write(iulog, '(/, A)') " elm-erw namelist:"
       write(iulog, '(A, " : ", I0,/)') "   elm-erw beginning year ", year_start_erw
       write(iulog, '(A, " : ", I0,/)') "   elm-erw calibration years ", nyear_erw_calibrate
       write(iulog, '(A, " : ", A,/)')  "   elm-erw parameter file ", trim(elm_erw_paramfile)
       write(iulog, '(A, " : ", I0,/)') "   verbose logs ", use_erw_verbose
       write(iulog, '(A, " : ", I0,/)') "   built-in validation site ", builtin_site
       write(iulog, '(A, " : ", I0,/)') "   number of soil layers to mix rock powder", mixing_layer
    end if

    ! Broadcast namelist variables read in
    call mpi_bcast (year_start_erw, 1, MPI_INTEGER, 0, mpicom, ierr)
    call mpi_bcast (nyear_erw_calibrate, 1, MPI_INTEGER, 0, mpicom, ierr)
    call mpi_bcast (elm_erw_paramfile, len(elm_erw_paramfile), MPI_CHARACTER, 0, mpicom, ierr)
    call mpi_bcast (use_erw_verbose, 1, MPI_INTEGER, 0, mpicom, ierr)
    call mpi_bcast (builtin_site, 1, MPI_INTEGER, 0, mpicom, ierr)
    call mpi_bcast (mixing_layer, 1, MPI_INTEGER, 0, mpicom, ierr)
  end subroutine elm_erw_readnl

  !-----------------------------------------------------------------------
  subroutine readEnhancedWeatheringParams ( elm_erw_paramfile )
    !
    ! !DESCRIPTION:
    ! Read in parameters
    !
    ! !USES:
    use fileutils   , only : getfil
    use ncdio_pio   , only : ncd_pio_closefile, ncd_pio_openfile
    use ncdio_pio   , only : file_desc_t, ncd_inqdid, ncd_inqdlen
    use ncdio_pio   , only : ncd_io
    use abortutils  , only : endrun
    use shr_log_mod , only : errMsg => shr_log_errMsg
    use elm_varpar  , only : nminerals, nks, ncations, nminsecs
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: elm_erw_paramfile ! parameter filename

    ! !LOCAL VARIABLES:
    character(len=256) :: locfn   ! local file name
    type(file_desc_t)  :: ncid    ! pio netCDF file id
    integer            :: dimid   ! netCDF dimension id

    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv   ! has variable been read in or not
    character(len=100) :: tString ! temp. var for reading

    character(len=32)  :: subname = 'EWParamsType'

    !-----------------------------------------------------------------------

    !
    call getfil (elm_erw_paramfile, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    call ncd_inqdid(ncid,'nminerals',dimid)
    call ncd_inqdlen(ncid,dimid,nminerals)       ! note this will override value from 'elm_varpar' initials

    call ncd_inqdid(ncid,'ncations',dimid)
    call ncd_inqdlen(ncid,dimid,ncations)        ! note this will override value from 'elm_varpar' initials

    call ncd_inqdid(ncid,'nminsecs',dimid)
    call ncd_inqdlen(ncid,dimid,nminsecs)       ! note this will override value from 'elm_varpar' initials

    allocate(character(40) :: EWParamsInst%minerals_name(1:nminerals))
    allocate(EWParamsInst%primary_mass(1:nminerals))
    allocate(EWParamsInst%log_keq_primary(1:nminerals))
    allocate(EWParamsInst%log_k_primary(1:nminerals, 1:nks))   ! 'nks' NOT read-in as above
    allocate(EWParamsInst%e_primary(1:nminerals, 1:nks))
    allocate(EWParamsInst%n_primary(1:nminerals, 1:nks))

    allocate(EWParamsInst%primary_stoi_proton(1:nminerals))
    allocate(EWParamsInst%primary_stoi_cations(1:nminerals,1:ncations))
    allocate(EWParamsInst%primary_stoi_h2o(1:nminerals))
    allocate(EWParamsInst%primary_stoi_sio2(1:nminerals))
    allocate(EWParamsInst%primary_stoi_hco3(1:nminerals))

    allocate(character(40) :: EWParamsInst%minsecs_name(1:nminsecs))
    allocate(EWParamsInst%minsecs_mass(1:nminsecs))
    allocate(EWParamsInst%log_keq_minsecs(1:nminsecs))
    allocate(EWParamsInst%ssa_minsecs(1:nminsecs))
    allocate(EWParamsInst%k_precip_minsecs(1:nminsecs, 1:nks))  ! 'nks' NOT read-in as above
    allocate(EWParamsInst%e_precip_minsecs(1:nminsecs, 1:nks))  ! 'nks' NOT read-in as above
    allocate(EWParamsInst%ph2o_precip_minsecs(1:nminsecs))
    allocate(EWParamsInst%qh2o_precip_minsecs(1:nminsecs))
    allocate(EWParamsInst%n_precip_minsecs(1:nminsecs))
    allocate(EWParamsInst%k_dissolv_minsecs(1:nminsecs, 1:nks))
    allocate(EWParamsInst%e_dissolv_minsecs(1:nminsecs, 1:nks))
    allocate(EWParamsInst%n_dissolv_minsecs(1:nminsecs, 1:nks))

    allocate(character(40) :: EWParamsInst%cations_name(1:ncations))
    allocate(EWParamsInst%cations_mass(1:ncations))
    allocate(EWParamsInst%cations_valence(1:ncations))
    allocate(EWParamsInst%cations_diffusivity(1:ncations))

    allocate(EWParamsInst%bicarbonate_diffusivity)
    allocate(EWParamsInst%carbonate_diffusivity)

    ! read in parameters
    tString='minerals_name'
    call ncd_io(varname=trim(tString),data=EWParamsInst%minerals_name, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_mass'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_mass, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='log_keq_primary'
    call ncd_io(varname=trim(tString),data=EWParamsInst%log_keq_primary, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='log_k_primary'
    call ncd_io(varname=trim(tString),data=EWParamsInst%log_k_primary, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='e_primary'
    call ncd_io(varname=trim(tString),data=EWParamsInst%e_primary, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='n_primary'
    call ncd_io(varname=trim(tString),data=EWParamsInst%n_primary, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    ! for primary mineral's dissolution reactions (product)
    tString='primary_stoi_proton'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_proton, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_stoi_cations'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_cations, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_stoi_h2o'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_h2o, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_stoi_sio2'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_sio2, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_stoi_hco3'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_hco3, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    ! for secondary mineral precipitions
    tString='minsecs_name'
    call ncd_io(varname=trim(tString),data=EWParamsInst%minsecs_name, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='minsecs_mass'
    call ncd_io(varname=trim(tString),data=EWParamsInst%minsecs_mass, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='log_keq_minsecs'
    call ncd_io(varname=trim(tString),data=EWParamsInst%log_keq_minsecs, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='ssa_minsecs'
    call ncd_io(varname=trim(tString),data=EWParamsInst%ssa_minsecs, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='k_precip_minsecs'
    call ncd_io(varname=trim(tString),data=EWParamsInst%k_precip_minsecs, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='e_precip_minsecs'
    call ncd_io(varname=trim(tString),data=EWParamsInst%e_precip_minsecs, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='ph2o_precip_minsecs'
    call ncd_io(varname=trim(tString),data=EWParamsInst%ph2o_precip_minsecs, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='qh2o_precip_minsecs'
    call ncd_io(varname=trim(tString),data=EWParamsInst%qh2o_precip_minsecs, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='n_precip_minsecs'
    call ncd_io(varname=trim(tString),data=EWParamsInst%n_precip_minsecs, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='k_dissolv_minsecs'
    call ncd_io(varname=trim(tString),data=EWParamsInst%k_dissolv_minsecs, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='e_dissolv_minsecs'
    call ncd_io(varname=trim(tString),data=EWParamsInst%e_dissolv_minsecs, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='n_dissolv_minsecs'
    call ncd_io(varname=trim(tString),data=EWParamsInst%n_dissolv_minsecs, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    ! for cations involving in primary mineral's dissolutions
    tString='cations_name'
    call ncd_io(varname=trim(tString),data=EWParamsInst%cations_name, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='cations_mass'
    call ncd_io(varname=trim(tString),data=EWParamsInst%cations_mass, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='cations_valence'
    call ncd_io(varname=trim(tString),data=EWParamsInst%cations_valence, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='cations_diffusivity'
    call ncd_io(varname=trim(tString),data=EWParamsInst%cations_diffusivity, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='bicarbonate_diffusivity'
    call ncd_io(varname=trim(tString),data=EWParamsInst%bicarbonate_diffusivity, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='carbonate_diffusivity'
    call ncd_io(varname=trim(tString),data=EWParamsInst%carbonate_diffusivity, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    ! close nc file
    call ncd_pio_closefile(ncid)

  end subroutine readEnhancedWeatheringParams


  !-----------------------------------------------------------------------
  subroutine MineralInit(bounds, num_soilc, filter_soilc, soilstate_vars)
    !
    ! !DESCRIPTION: 
    ! Calculate initial cation concentration from background CEC and soil pH
    ! after soil hydrology is already initialized
    ! 
    ! !USES:
    use elm_varctl, only : use_erw_verbose
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: fc,c,j,t,g,l,nlevbed
    integer  :: icat                   ! indices
    real(r8) :: co2_atm                ! CO2 partial pressure in atm
    character(2) :: j_str
    character(4) :: j_lev
    real(r8) :: beta_list(1:ncations)  ! temporary container for cece/cec_tot
    real(r8) :: beta_h,h               ! temporary container for ceca/cec_tot and [H+]
    real(r8) :: keq_list(1:ncations)   ! temporary container for exchange coefficients between H+ and cations

    associate( &
         soil_ph                        => col_ms%soil_ph                 , & ! Input:  [real(r8) (:,:)] calculated soil pH (1:nlevgrnd)
         h2osoi_vol                     => col_ws%h2osoi_vol              , & ! Input:  [real(r8) (:)] volumetric soil water content, ice + water (m3 m-3)
         nlev2bed                       => col_pp%nlevbed                 , & ! Input:  [integer  (:)   ] number of layers to bedrock
         proton_vr                      => col_ms%proton_vr               , & ! Output: [real(r8)(:)   ] calculated soil H+ concentration in soil water each soil layer (1:nlevgrnd) (g m-3 soil [not water])
         bicarbonate_vr                 => col_ms%bicarbonate_vr          , & ! Output:  [real(r8) (:,:)] calculated HCO3- concentration in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)
         carbonate_vr                   => col_ms%carbonate_vr            , & ! Output:  [real(r8) (:,:)] calculated CO3 2- concentration in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)
         cec_cation_vr                  => col_ms%cec_cation_vr           , & ! Output: [real(r8) (:,:,:)] adsorbed cation concentration each soil layer (1:nlevgrnd,1:ncations) (g m-3 soil [not dry soil])
         cec_proton_vr                  => col_ms%cec_proton_vr           , & ! Output: [real(r8) (:,:,:)] adsorbed H+ concentration each soil layer (1:nlevgrnd) (g m-3 soil [not dry soil])
         net_charge_vr                  => col_ms%net_charge_vr           , & ! Output:  [real(r8) (:,:)] net charge of the tracked ions in the soil solution system, constant over time (1:nlevgrnd) (mol kg-1)
         equilibria_conc                => col_ms%equilibria_conc         , & ! Output:  [real(r8) (:,:,:)] soil pore water cation concentration implied by the input soil CEC status and exchange coefficients (mol kg-1)
         cation_vr                      => col_ms%cation_vr               , & ! Output [real(r8) (:,:,:)] cation mass in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd, 1:ncations)
         cect_dyn                       => col_ms%cect_dyn                , & ! Input:  [real(r8) (:,:)] pH-dependent total cation exchange capacity (1:nlevgrnd)
         secondary_mineral_vr           => col_ms%secondary_mineral_vr    & ! Output: [real(r8) (:,:)] secondary mineral content in each soil layer (1:nlevgrnd, 1:nminsecs)
    )

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      t = col_pp%topounit(c)
      g = col_pp%gridcell(c)
      l = col_pp%landunit(c)
      nlevbed = min(nlev2bed(c), nlevsoi)

      ! 2x CO2 run
      !co2_atm = 2 * top_as%pco2bot(t) / 101325
      co2_atm = top_as%pco2bot(t) / 101325

      do j = 1,nlevbed
        ! soil pH has been initialized at cold-start
        proton_vr(c,j) = logmol_to_mass(-soil_ph(c,j), mass_h, h2osoi_vol(c,j))

        ! for the sake of initial charge balance, calculate from soil pH
        bicarbonate_vr(c,j) = mol_to_mass( ph_to_hco3(soil_ph(c,j), co2_atm), mass_hco3, &
                                           h2osoi_vol(c,j) )
        carbonate_vr(c,j) = mol_to_mass( hco3_to_co3(ph_to_hco3(soil_ph(c,j), co2_atm), &
                                           soil_ph(c,j)), mass_co3, h2osoi_vol(c,j) )

        ! initial CEC H+ and cation status
        cec_proton_vr(c,j) = meq_to_mass(soilstate_vars%ceca_col(c,j), 1._r8, mass_h, &
                                         soilstate_vars%bd_col(c,j))
        do icat = 1, ncations
          cec_cation_vr(c,j,icat) = meq_to_mass(soilstate_vars%cece_col(c,j,icat), &
                                                EWParamsInst%cations_valence(icat), &
                                                EWParamsInst%cations_mass(icat), &
                                                soilstate_vars%bd_col(c,j))
        end do

        ! calculate initial soil solution concentration using the equilibrium with CEC
        h = 10**(-soilstate_vars%sph(c,j))
        beta_h = soilstate_vars%ceca_col(c,j) / soilstate_vars%cect_col(c,j)
        do icat = 1,ncations
          beta_list(icat) = soilstate_vars%cece_col(c,j,icat) / soilstate_vars%cect_col(c,j)
          keq_list(icat) = 10**soilstate_vars%log_km_col(c,j,icat)

          equilibria_conc(c,j,icat) = beta_list(icat) / &
            (beta_h * keq_list(icat) / h)**EWParamsInst%cations_valence(icat)
          cation_vr(c,j,icat) = mol_to_mass(equilibria_conc(c,j,icat), &
                                            EWParamsInst%cations_mass(icat), &
                                            h2osoi_vol(c,j))
        end do

        ! calculate the net charge balance during initializaiong
        ! mol/kg
        net_charge_vr(c,j) = find_net_charge(soilstate_vars%sph(c,j), co2_atm, beta_list, &
                                             keq_list, EWParamsInst%cations_valence)
        
        ! reset the total cation exchange capacity
        cect_dyn(c,j) = soilstate_vars%cect_col(c,j)

        ! initial secondary mineral contents
        ! convert from g 100g-1 soil to g m-3 soil
        secondary_mineral_vr(c,j,1) = soilstate_vars%calcite_col(c,j) * 10._r8 * soilstate_vars%bd_col(c,j)
        secondary_mineral_vr(c,j,2) = soilstate_vars%kaolinite_col(c,j) * 10._r8 * soilstate_vars%bd_col(c,j)
        secondary_mineral_vr(c,j,3) = soilstate_vars%gibbsite_col(c,j) * 10._r8 * soilstate_vars%bd_col(c,j)
      end do
    end do

    if (use_erw_verbose > 0) then
      write (100+iam, *) '*************************************************************************'
      write (100+iam, *) '*** Soil Initialization for Enhanced Weathering                       ***'
      do fc = 1,num_soilc
        c = filter_soilc(fc)
        g = col_pp%gridcell(c)

        ! Write diagnostics
        write (100+iam, *) '-----------------------------------------------------------------------'
        write (100+iam, *) 'grid and column: ', ldomain%latc(g), ldomain%lonc(g), g, c
        write (100+iam, *) 'soil pH, proton (g m-3), cation (g m-3):'
        do j = 1,nlevbed
          write (j_str, '(I2)') j
          j_lev = 'j='//j_str
          ! write (100+iam, *) 'j=', j, soil_ph(c,j), h2osoi_vol(c,j), nlevbed, nlevsoi
          write (100+iam, *) j_lev, soil_ph(c,j), proton_vr(c,j), cation_vr(c,j,1:ncations)
          ! write (100+iam, *) 'cation_vr',  j, icat, cation_vr(c,j,icat), soil_ph(c,j), soilstate_vars%log_km_col(c,j,icat), soilstate_vars%ceca_col(c,j), soilstate_vars%cect_col(c,j), EWParamsInst%cations_valence(icat), soilstate_vars%cece_col(c,j,icat), h2osoi_vol(c,j)    
        end do
        write (100+iam, *) 'cec total (meq 100g-1 dry soil), beta H+, beta cation:'
        do j = 1,nlevbed
          write (j_str, '(I2)') j
          j_lev = 'j='//j_str
          write (100+iam, *) j_lev, soilstate_vars%cect_col(c,j), soilstate_vars%ceca_col(c,j) / soilstate_vars%cect_col(c,j), soilstate_vars%cece_col(c,j,1:ncations) / soilstate_vars%cect_col(c,j)
        end do
        write (100+iam, *) 'calcite, kaolinite, gibbsite (g 100g-1 soil):'
        do j = 1,nlevbed
          write (j_str, '(I2)') j
          j_lev = 'j='//j_str
          write (100+iam, *) j_lev, secondary_mineral_vr(c,j,1:nminsecs)
        end do
      end do
      write (100+iam, *) '*************************************************************************'
      call shr_sys_flush(100+iam)
    end if

    end associate
  end subroutine MineralInit


  !-----------------------------------------------------------------------
  subroutine MineralBackground(bounds, num_soilc, filter_soilc, soilstate_vars)
    !
    ! !DESCRIPTION: 
    ! Calculate the background weathering fluxes. 
    ! 
    ! !USES:
    use elm_varcon       , only : spval, secspday
    use elm_varctl       , only : builtin_site
    use elm_time_manager , only : get_step_size, get_curr_date
    use abortutils       , only : endrun
    use timeinfoMod
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilstate_type)     , intent(in)    :: soilstate_vars

    !
    ! !LOCAL VARIABLES:
    integer  :: fc,c,j,g,nlevbed
    integer  :: m, isec, icat          ! indices
    real(r8) :: dt
    real(r8) :: fracday                ! fractional day per time step
    real(r8) :: temp_sum(1:ncations), temp_cec_sum(1:ncations), temp_minsecs_sum(1:nminsecs)
    real(r8) :: temp_flux_sum(1:ncations), temp_flux_cec_sum(1:ncations), temp_flux_minsecs_sum(1:nminsecs)
    integer  :: kyr                    ! current year
    integer  :: kmo                    ! month of year  (1, ..., 12)
    integer  :: kda                    ! day of month   (1, ..., 31)
    integer  :: mcsec                  ! seconds 
    integer  :: current_date
    logical  :: is_hour                ! apply rock only at 10th hour of day

    associate( &
         nlev2bed                       => col_pp%nlevbed                 , & ! Input:  [integer  (:)   ]  number of layers to bedrock
         dz                             => col_pp%dz                      , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)

         rain_ph                        => col_ew%rain_ph                 , & ! Input:  [real(r8) (:)] pH of rain water
         rain_chem                      => col_ew%rain_chem               , & ! Input:  [real(r8) (:,:)] cation concentration in rain water (excluding H+) (g m-3 rain water) (1:ncations)

         annavg_tot_delta               => col_mf%annavg_tot_delta        , & ! Output: [real(r8) (:)] annual average rate of change in total (solute and adsorbed phases) cation concentration before mineral application (g m-3 soil s-1 [not dry soil])
         annavg_cec_delta               => col_mf%annavg_cec_delta        , & ! Output: [real(r8) (:)] annual average rate of change in adsorbed cation concentration before mineral application (g m-3 soil s-1 [not dry soil])
         annavg_minsecs_delta           => col_mf%annavg_minsecs_delta    , & ! Output: [real(r8) (:,:)] annual average rate of change in secondary mineral content before mineral application (1:nminsecs) (g m-3 soil s-1 [not dry soil])

         background_flux_vr             => col_mf%background_flux_vr      , & ! Output: [real(r8) (:)] background flux rate in/out of soil solution (g m-3 s-1)
         background_cec_vr              => col_mf%background_cec_vr       , & ! Output: [real(r8) (:)] background flux rate in/out of adsorbed cations (g m-3 s-1)
         background_minsecs_vr          => col_mf%background_minsecs_vr   , & ! Output: [real(r8) (:,:)] background flux rate in/out of secondary minerals (1:nminsecs) (g m-3 s-1)

         !
         ! Forcing variables for built-in validation sites
         !
         forc_app                       => col_ew%forc_app                 , & ! Input:  [real(r8) (:)] application rate (kg rock m-2 year-1)
         forc_min                       => col_ew%forc_min                 , & ! Input:  [real(r8) (:,:) weight percentage of minerals in rock (1:nminerals) (kg mineral kg-1 rock)
         forc_pho                       => col_ew%forc_pho                 , & ! Input:  [real(r8) (:)] weight percentage of phosphorus content in rock (gP kg-1 rock)
         forc_gra                       => col_ew%forc_gra                 , & ! Input:  [real(r8) (:,:)] grain size (1:nminerals) (um diameter)

         ! initial soil secondary mineral contents
         secondary_mineral_vr          => col_ms%secondary_mineral_vr     & ! Input:  [real(r8) (:,:)] secondary mineral content in each soil layer (1:nlevgrnd, 1:nminsecs) (g m-3 soil [not dry soil])
    )

    dt      = real( get_step_size(), r8 )
    fracday = dt / secspday

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      g = col_pp%gridcell(c)
      nlevbed = min(nlev2bed(c), nlevsoi)

      !------------------------------------------------------------------------------
      ! Background weathering
      ! - goal: prevent forever loss of cations from soil
      ! - solution: replenish both the solute and cation exchange (adsorbed) phase
      !             cations using calibrated long-term average; also replenish
      !             the natural dissolution of initial secondary mineral contents
      !------------------------------------------------------------------------------
      if (year_curr < (year_start_erw + nyear_erw_calibrate)) then
      !! performs less well: let background_flux co-evolve with cumulative loss
      !! if (year_curr < year_start_erw) then

        background_flux_vr(c,1:nlevbed,1:ncations) = 0._r8
        background_cec_vr(c,1:nlevbed,1:ncations) = 0._r8
        background_minsecs_vr(c,1:nlevbed,1:nminsecs) = 0._r8

      else

        temp_cec_sum(1:ncations) = 0._r8 ! true sum of abs(CEC loss) over all soil layers
        temp_flux_cec_sum(1:ncations) = 0._r8 ! sum of max(abs(CEC loss), 0) over all soil layers
        do j = 1,nlevbed
          do icat = 1,ncations
            temp_cec_sum(icat) = temp_cec_sum(icat) - annavg_cec_delta(c,j,icat)*dz(c,j)
            background_cec_vr(c,j,icat) = max(0._r8, - annavg_cec_delta(c,j,icat))
            temp_flux_cec_sum(icat) = temp_flux_cec_sum(icat) + background_cec_vr(c,j,icat)*dz(c,j)
          end do
        end do

        ! re-scale the background weathering rate so that the column total equals
        ! the column total of annavg_cec_delta, which may contain negative numbers
        do icat = 1,ncations
          if (temp_cec_sum(icat) <= 0._r8 .or. temp_flux_cec_sum(icat) <= 0._r8) then
            background_cec_vr(c,1:nlevbed,icat) = 0._r8
          else
            do j = 1,nlevbed
              background_cec_vr(c,j,icat) = background_cec_vr(c,j,icat) * temp_cec_sum(icat) / temp_flux_cec_sum(icat)
            end do
          end if
        end do

        temp_sum(1:ncations) = 0._r8 ! true sum of abs(cation loss) over all soil layers
        temp_flux_sum(1:ncations) = 0._r8 ! sum of max(abs(cation loss), 0) over all soil layers
        do j = 1,nlevbed
          do icat = 1,ncations
            temp_sum(icat) = temp_sum(icat) - annavg_tot_delta(c,j,icat)*dz(c,j) - background_cec_vr(c,j,icat)*dz(c,j)
            background_flux_vr(c,j,icat) = max(0._r8, - annavg_tot_delta(c,j,icat) - background_cec_vr(c,j,icat))
            temp_flux_sum(icat) = temp_flux_sum(icat) + background_flux_vr(c,j,icat)*dz(c,j)
          end do
        end do

        ! re-scale the background weathering rate so that the column total equals
        ! the column total of (annavg_tot_delta - background_cec), which may 
        ! contain negative numbers
        temp_minsecs_sum(1:nminsecs) = 0._r8 ! true sum of abs(minsec loss) over all soil layers
        temp_flux_minsecs_sum(1:nminsecs) = 0._r8 ! sum of max(abs(minsec loss), 0) over all soil layers
        do j = 1,nlevbed
          do isec = 1,nminsecs
            temp_minsecs_sum(isec) = temp_minsecs_sum(isec) - annavg_minsecs_delta(c,j,isec)*dz(c,j)
            background_minsecs_vr(c,j,isec) = max(0._r8, - annavg_minsecs_delta(c,j,isec))
            temp_flux_minsecs_sum(isec) = temp_flux_minsecs_sum(isec) + background_minsecs_vr(c,j,isec)*dz(c,j)
          end do
        end do

        do isec = 1,nminsecs
          if (temp_minsecs_sum(isec) <= 0._r8 .or. temp_flux_minsecs_sum(isec) <= 0._r8) then
            background_minsecs_vr(c,1:nlevbed,isec) = 0._r8
          else
            do j = 1,nlevbed
              background_minsecs_vr(c,j,isec) = background_minsecs_vr(c,j,isec) * temp_minsecs_sum(isec) / temp_flux_minsecs_sum(isec)
            end do
          end if
        end do
      end if

      !------------------------------------------------------------------------------
      ! Define rainfall chemistry
      !------------------------------------------------------------------------------
      !if (builtin_site == 1) then
        ! Hubbard Brook
        ! rain pH data from the monitoring station in Hubbard Brook, 
        ! National Atmospheric Deposition Program
        ! https://nadp.slh.wisc.edu/sites/ntn-NH02/

      !  rain_ph(c) = 4.8_r8
        ! Ca = 0.055 mg/L, Mg = 0.015 mg/L, Na = 0.075 mg/L, K = 0.014 mg/L, Al = 0
      !  rain_chem(c, 1) = 0.055_r8
      !  rain_chem(c, 2) = 0.015_r8
      !  rain_chem(c, 3) = 0.075_r8
      !  rain_chem(c, 4) = 0.014_r8
      !  rain_chem(c, 5) = 0._r8

      !else if (builtin_site == 2) then
        ! U.C. Davis
        ! rain pH data from the monitoring station in Davis,  
        ! National Atmospheric Deposition Program
        ! https://nadp.slh.wisc.edu/sites/ntn-CA88/
      !  rain_ph(c) = 6.2_r8
      !  ! Ca = 0.055 mg/L, Mg = 0.015 mg/L, Na = 0.075 mg/L, K = 0.014 mg/L, Al = 0
      !  rain_chem(c, 1) = 0.06_r8
      !  rain_chem(c, 2) = 0.06_r8
      !  rain_chem(c, 3) = 0.025_r8
      !  rain_chem(c, 4) = 0.04_r8
      !  rain_chem(c, 5) = 0._r8

     !else
        !
        rain_ph(c) = 5.6_r8
        rain_chem(c, :) = 0.0_r8 ! in the new setup, rain_chem should no longer matter
      !end if


      ! ---------------------------------------------------------------
      ! site-specific over-write of forcing
      ! ---------------------------------------------------------------
      call get_curr_date(kyr, kmo, kda, mcsec)
      current_date = kyr*10000 + kmo*100 + kda

      ! if builtin_site > 0, manually overwrite the read-in values
      if (builtin_site == 1) then
        is_hour = (secs_curr-36000.0_r8) > -dt/2.0_r8 .and. (secs_curr-36000.0_r8) <= dt/2.0_r8
        if (current_date .eq. 19991019 .and. is_hour) then
            ! 55 tons / 11.8 ha = 0.466 kg / m2, applied over one day
            ! purity in Table 2: 89.9%
            ! Johnson, C. E., Driscoll, C. T., Blum, J. D., Fahey, T. J., & Battles, J. J. (2014). Soil Chemical Dynamics after Calcium Silicate Addition to a Northern Hardwood Forest. Soil Science Society of America Journal, 78(4), 1458–1468. https://doi.org/10.2136/sssaj2014.03.0114
            forc_app(c) = 0.466_r8 * 0.9
        else
            forc_app(c) = 0._r8
        end if
        forc_min(c, 1:nminerals) = 0._r8
        forc_min(c, 1) = 1._r8

        ! BET surface area is 1.6 m2/g; inversely calculate to make sure 
        ! SSA matches, instead of using the provided grain size
        ! 
        ! Johnson, C. E., Driscoll, C. T., Blum, J. D., Fahey, T. J., & Battles, J. J. (2014). Soil Chemical Dynamics after Calcium Silicate Addition to a Northern Hardwood Forest. Soil Science Society of America Journal, 78(4), 1458–1468. https://doi.org/10.2136/sssaj2014.03.0114
        if ((current_date .lt. 19991019) .or. (current_date .eq. 19991019 .and. is_hour)) then
          forc_gra(c, 1:nminerals) = 20.8_r8 ! 9.6 um
        end if

        forc_pho(c   ) = 0._r8

      else if (builtin_site == 2) then
        is_hour = (secs_curr-36000.0_r8) > -dt/2.0_r8 .and. (secs_curr-36000.0_r8) <= dt/2.0_r8

        if (is_hour .and. ((kyr .eq. 2019) .or. (kyr .eq. 2020)) .and. &
            ((kmo .eq. 9) .or. (kmo .eq. 10) .or. (kmo .eq. 11))) then
          ! 40 t ha-1 = 4 kg / m2, applied over 3 months, convert to per day
          forc_app(c) = 4._r8 / 90._r8
        else
          forc_app(c) = 0._r8
        end if
        forc_min(c,1:nminerals) = 0._r8
        forc_min(c,3) = 0.167_r8
        forc_min(c,4) = 0.167_r8
        forc_min(c,5) = 0.143_r8
        forc_gra(c, 1:nminerals) = 107._r8

        forc_pho(c   ) = 0.0005365_r8

      else if (builtin_site == 3) then
        ! University of Illinois Energy Farm
        is_hour = (secs_curr-36000.0_r8) > -dt/2.0_r8 .and. (secs_curr-36000.0_r8) <= dt/2.0_r8

        if (is_hour .and. (kyr >= 2016) .and. (kyr <= 2019) .and. (kmo == 11)) then
        ! if (is_hour .and. (kyr >= 1850) .and. (kyr <= 1852) .and. (kmo == 11)) then
          ! 50 t ha−1 y−1 = 5 kg / m2, applied over 1 month, convert to per day
          forc_app(c) = 5._r8 / 30._r8
        else
          forc_app(c) = 0._r8
        end if

        !!write (iulog, *) kyr, kmo, kda, c, secs_curr, forc_app(c)

        ! The minerals in the parameter file are different from above
        forc_min(c,1:nminerals) = 0._r8
        forc_min(c,3) = 0.233_r8
        forc_min(c,5) = 0.178_r8
        forc_min(c,6) = 0.026_r8
        forc_min(c,7) = 0.116_r8
        forc_min(c,8) = 0.340_r8
        ! p80 = 267 um => use exponential distribution to convert 
        ! P50 = P80 / 1.218 / 4.395 = 73.994
        forc_gra(c, 1:nminerals) = 73.994_r8

        forc_pho(c   ) = 0._r8

      else

        ! from read-in data of soil amendment application
        !write (iulog, *), nstep_mod, jday_mod, secs_curr, forc_app(c)
        !if (forc_app(c) > 0._r8) then
        !do m = 1, nminerals
        !    write (iulog, *), m, forc_min(c, m), forc_gra(c, m)
        !end do
        !end if

        !! do some checking
        !do m = 1, nminerals
        !  if (forc_gra(c,m) /= forc_gra(c,m)) then
        !    !'nan' cause math issue in following calculation
        !    write (iulog, *), c, m, forc_gra(c, m)
        !    call endrun(msg='Nan of powder grain size')
        !  end if
        !end do

        ! Need P content to affect NEE
        forc_pho(c) = 0._r8  ! (TODO) 0~namendnutr element(s) from soil amend application, by given index of phosphrous

      end if

    end do ! end column loop

    end associate

  end subroutine MineralBackground

  !-----------------------------------------------------------------------
  subroutine MineralPrimary(bounds, num_soilc, filter_soilc, soilstate_vars)
    !
    ! !DESCRIPTION: 
    ! Calculate primary mineral dissolution fluxes.
    ! 
    ! !USES:
    ! rgas = universal gas constant [= 8314.467 J/K/kmole]
    use elm_varcon       , only : spval, rgas, secspday, D_h
    use elm_varctl       , only : use_erw_verbose
    use elm_time_manager , only : get_step_size, get_curr_date
    use abortutils       , only : endrun
    use SharedParamsMod  , only : ParamsShareInst
    use timeinfoMod
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilstate_type)     , intent(in)    :: soilstate_vars

    !
    ! !LOCAL VARIABLES:
    integer  :: fc,c,j,g,nlevbed
    integer  :: m, isec, icat          ! indices
    real(r8) :: dt
    real(r8) :: log_k_dissolve_acid, log_k_dissolve_neutral, log_k_dissolve_base
    real(r8) :: k_tot
    real(r8) :: phi                    ! porosity of the passivation layer
    real(r8) :: tau                    ! tortuosity of the passivation layer
    real(r8) :: D_h_eff                ! effective diffusivity of H+ in the passivation layer (m2/s)
    real(r8) :: Jflux                  ! H+ sink strength due to previous step's dissolution (mol m-2 s-1)
    real(r8) :: dNb_dt                 ! H+ diffusion limited dissolution rate (mol m-3 s-1)

    associate( &
         !
         ! Forcing variables for built-in validation sites
         !
         forc_app                       => col_ew%forc_app                 , & ! Input:  [real(r8) (:)] application rate (kg rock m-2 year-1)
         forc_min                       => col_ew%forc_min                 , & ! Input:  [real(r8) (:,:) weight percentage of minerals in rock (1:nminerals) (kg mineral kg-1 rock)
         forc_pho                       => col_ew%forc_pho                 , & ! Input:  [real(r8) (:)] weight percentage of phosphorus content in rock (gP kg-1 rock)
         forc_gra                       => col_ew%forc_gra                 , & ! Input:  [real(r8) (:,:)] grain size (1:nminerals) (um diameter)

         !
         ! soil pH and ionic states 
         !
         soil_ph                        => col_ms%soil_ph                 , & ! Input: [real(r8) (:,:)] calculated soil pH (1:nlevgrnd)
         bicarbonate_vr                 => col_ms%bicarbonate_vr          , & ! Output: [real(r8) (:,:)] calculated HCO3- concentration in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)
         carbonate_vr                   => col_ms%carbonate_vr            , & ! Output: [real(r8) (:,:)] calculated CO3 2- concentration in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)
         cation_vr                      => col_ms%cation_vr               , & ! Output [real(r8) (:,:,:)] cation mass in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd, 1:ncations)

         !
         ! Primary mineral state
         !
         primary_mineral_vr     => col_ms%primary_mineral_vr             , & ! Output [real(r8) (:,:,:)] primary mineral mass in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd, 1:nminerals)

         silica_vr              => col_ms%silica_vr                      , & ! Output [real(r8) (:,:)] silica mass in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)
         passivation_thickness  => col_ms%passivation_thickness          , & ! Output [real(r8) (:,:)] thickness of the armoring layer on the primary mineral (m) (1:nlevgrnd)
         ssa_dyn                => col_ms%ssa_dyn                        , & ! Output [real(r8) (:,:,:)] specific surface area of the primary minerals (m2 g-1 mineral) (1:nlevgrnd,1:nminerals)

         !
         ! Primary mineral flux
         !
         primary_added_vr               => col_mf%primary_added_vr       , & ! Output [real(r8) (:,:,:)] primary mineral addition through rock powder application (g m-3 s-1) (1:nlevgrnd, 1:nminerals)
         primary_dissolve_vr            => col_mf%primary_dissolve_vr    , & ! Output [real(r8) (:,:,:)] primary mineral loss through dissolution reaction (g m-3 s-1) (1:nlevgrnd, 1:nminerals)

         primary_cation_flux_vr         => col_mf%primary_cation_flux_vr , & ! Output [real(r8) (:,:,:) cations produced due to all the dissolution reactions (g m-3 s-1) (1:nlevgrnd, 1:ncations)
         primary_h2o_flux_vr            => col_mf%primary_h2o_flux_vr    , & ! Output [real(r8) (:,:)] net of water produced and consumed due to all the dissolution reaction (g m-3 s-1) (1:nlevgrnd)
         primary_silica_flux_vr         => col_mf%primary_silica_flux_vr , & ! Output [real(r8) (:,:)] SiO2 produced due to all the dissolution reaction (g m-3 s-1) (1:nlevgrnd)
         primary_residue_flux_vr        => col_mf%primary_residue_flux_vr, & ! Output [real(r8) (:,:)] Non-SiO2 solides produced due to all the dissolution reaction (g m-3 s-1) (1:nlevgrnd,1:nminerals)

         primary_prelease_vr            => col_mf%primary_prelease_vr    , & ! Output [real(r8) (:,:)] release of soluble phosphorus due to all the dissolution reaction (g m-3 s-1) (1:nlevgrnd)

         r_dissolve_vr                  => col_mf%r_dissolve_vr          , & ! Output [real(r8) (:,:)] rate at which the dissolution reaction happens (mol m-3 s-1) (1:nlevgrnd, 1:nminerals)
         log_omega_vr                   => col_mf%log_omega_vr           , & ! Output [real(r8) (:,:)] omega parameter in the dissolution equation (1:nlevgrnd, 1:nminerals)

         !
         ! Other related
         !
         dz                             => col_pp%dz                       , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)
         tsoi                           => col_es%t_soisno                 , & ! Input: [real(r8) (:,:) ] soil temperature [K]
         h2osoi_vol                     => col_ws%h2osoi_vol               , & ! Input:  [real(r8) (:)] volumetric soil water content, ice + water (m3 m-3)
         h2osoi_liqvol                  => col_ws%h2osoi_liqvol            , & ! Input:  [real(r8) (:)] volumetric soil water content, liquid only (m3 m-3)
         nlev2bed                       => col_pp%nlevbed      & ! Input:  [integer  (:)   ]  number of layers to bedrock
    )

    dt      = real( get_step_size(), r8 )

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      g = col_pp%gridcell(c)
      !topo = col_pp%topounit(c)
      nlevbed = min(nlev2bed(c), nlevsoi)

      ! ---------------------------------------------------------------
      ! Apply the primary minerals
      ! ---------------------------------------------------------------
      do j = 1,nlevbed
        do m = 1,nminerals
          ! evenly distributed in the mixing_depth
          ! 20241206: forc_app is the annual applied amount; distribute it within one time step
          ! 20250714: round the application to the nearest discrete ELM vertical layer
          if (j <= mixing_layer) then
            primary_added_vr(c,j,m) = 1000._r8 * forc_app(c) * forc_min(c,m) / sum(dz(c,1:mixing_layer)) / dt
          else
            primary_added_vr(c,j,m) = 0._r8
          end if
        end do
      end do

      ! ---------------------------------------------------------------
      ! Primary mineral dissolution
      ! ---------------------------------------------------------------
      do j = 1,nlevbed
        if (j > mixing_layer .or. h2osoi_liqvol(c,j) < 1e-6) then
          r_dissolve_vr(c,j,1:nminerals) = 0._r8
        else
          ! Primary mineral dissolution
          do m = 1,nminerals
            !-------------------------------------------------------------------
            ! Calculate the dissolution rate based on soil pH
            !-------------------------------------------------------------------
            if (primary_mineral_vr(c,j,m) <= 0._r8) then
              r_dissolve_vr(c,j,m) = 0._r8

            else

              ! log10 of ion activity product divided by equilibrium constant
              log_omega_vr(c,j,m) = - EWParamsInst%log_keq_primary(m) + &
                                    soil_ph(c,j) * EWParamsInst%primary_stoi_proton(m)
              do icat = 1,ncations
                log_omega_vr(c,j,m) = log_omega_vr(c,j,m) + & 
                  EWParamsInst%primary_stoi_cations(m,icat) * &
                  mass_to_logmol(max(1.0e-15, cation_vr(c,j,icat)), &
                                 EWParamsInst%cations_mass(icat), h2osoi_vol(c,j))
              end do
              if (EWParamsInst%primary_stoi_hco3(m) > 0._r8) then
                log_omega_vr(c,j,m) = log_omega_vr(c,j,m) + EWParamsInst%primary_stoi_hco3(m) * &
                  mass_to_logmol(max(1.0e-15, bicarbonate_vr(c,j)), mass_hco3, h2osoi_vol(c,j))
              end if
              if (EWParamsInst%primary_stoi_sio2(m) > 0._r8) then
                log_omega_vr(c,j,m) = log_omega_vr(c,j,m) + EWParamsInst%primary_stoi_sio2(m) * &
                  mass_to_logmol(max(1.0e-15, silica_vr(c,j)), mass_sio2, h2osoi_vol(c,j))
              end if

              !!if (m == 6) then
              !!  write (iulog, *) c, j, m, 'log_omega', log_omega_vr(c,j,m)
              !!end if

              ! check the reaction rate is not negative
              if (log_omega_vr(c,j,m) >= 0._r8) then
                r_dissolve_vr(c,j,m) = 0._r8

                if (use_erw_verbose > 0) then
                  write (100+iam,*) ' WARNING! Omega > 1 meaning dissolution reaction cannot proceed', ldomain%latc(g), ldomain%lonc(g), g, c, j, m, log_omega_vr(c,j,m)

                  !write (100+iam, *) 'log_omega_vr part 1', soil_ph(c,j) * EWParamsInst%primary_stoi_proton(m) - EWParamsInst%log_keq_primary(m)
                  !do icat = 1,ncations
                  !  write (100+iam, *) 'log_omega_vr cation', icat, EWParamsInst%cations_name(icat), EWParamsInst%primary_stoi_cations(m,icat) * &
                  !  mass_to_logmol(cation_vr(c,j,icat), EWParamsInst%cations_mass(icat), h2osoi_vol(c,j))
                  !end do
                end if

              else
                ! log10 of the reaction rate constant by individual weathering agents (log10 mol m-2 s-1)
                k_tot = 0._r8

                ! only add this part if the rate is > -9999 (the NULL value, since base reaction is 
                ! sometimes not reported)
                if (EWParamsInst%log_k_primary(m,1) > -9000) then
                  log_k_dissolve_acid = EWParamsInst%log_k_primary(m,1) + log10(exp(1.0)) * & 
                    (-1.0e6_r8 * EWParamsInst%e_primary(m,1) / rgas * (1/tsoi(c,j) - 1/298.15_r8)) - &
                    EWParamsInst%n_primary(m,1) * soil_ph(c,j) + log10(1 - 10**log_omega_vr(c,j,m))
                  k_tot = k_tot + 10**log_k_dissolve_acid

                  !!if (m == 6) then
                  !!    write (iulog, *) c, j, m, 'log_k_dissolve_acid', log_k_dissolve_acid, &
                  !!          EWParamsInst%log_k_primary(m,1), log10(exp(1.0)) * & 
                  !!      (-1.0e6_r8 * EWParamsInst%e_primary(m,1) / rgas * (1/tsoi(c,j) - 1/298.15_r8)), EWParamsInst%n_primary(m,1) * soil_ph(c,j), log10(1 - 10**log_omega_vr(c,j,m))
                  !!end if

                  !write (100+iam, *) 'log_k_dissolve_acid', ldomain%latc(g), ldomain%lonc(g), g, c, j, m, log_k_dissolve_acid, EWParamsInst%log_k_primary(m,1), EWParamsInst%e_primary(m,1), rgas, tsoi(c,j), EWParamsInst%n_primary(m,1), soil_ph(c,j), log_omega_vr(c,j,m)
                end if

                if (EWParamsInst%log_k_primary(m,2) > -9000) then
                  log_k_dissolve_neutral = EWParamsInst%log_k_primary(m,2) + log10(exp(1.0)) * & 
                    (-1.0e6_r8 * EWParamsInst%e_primary(m,2) / rgas * (1/tsoi(c,j) - 1/298.15_r8)) + &
                    log10(1 - 10**log_omega_vr(c,j,m))
                  k_tot = k_tot + 10**log_k_dissolve_neutral
                  !write (100+iam, *) 'log_k_dissolve_neutral', ldomain%latc(g), ldomain%lonc(g), g, c, j, m, log_k_dissolve_neutral, EWParamsInst%log_k_primary(m,1), EWParamsInst%e_primary(m,1), rgas, tsoi(c,j), log_omega_vr(c,j,m)
                end if

                if (EWParamsInst%log_k_primary(m,3) > -9000) then
                  log_k_dissolve_base = EWParamsInst%log_k_primary(m,3) + log10(exp(1.0)) * & 
                    (-1.0e6_r8 * EWParamsInst%e_primary(m,3) / rgas * (1/tsoi(c,j) - 1/298.15_r8)) - &
                    EWParamsInst%n_primary(m,3) * (14 - soil_ph(c,j)) + log10(1 - 10**log_omega_vr(c,j,m))
                  k_tot = k_tot + 10**log_k_dissolve_base
                  ! write (100+iam, *) 'log_k_dissolve_base', ldomain%latc(g), ldomain%lonc(g), g, c, j, m, log_k_dissolve_base, EWParamsInst%log_k_primary(m,3), EWParamsInst%e_primary(m,3), rgas, tsoi(c,j), EWParamsInst%n_primary(m,3), soil_ph(c,j), log_omega_vr(c,j,m)
                end if

                !!if (m == 6) then
                !!  write (iulog, *) c, j, m, 'k_tot', k_tot
                !!end if

                ! further scale down the reaction rate using the 2/3 of saturation ratio
                ! Bao, C., Li, L., Shi, Y., & Duffy, C. (2017). Understanding watershed hydrogeochemistry: 1. Development of RT-Flux-PIHM. Water Resources Research, 53(3), 2328–2345. https://doi.org/10.1002/2016WR018934
                k_tot = k_tot * (h2osoi_liqvol(c,j) / soilstate_vars%watsat_col(c,j))**0.67_r8

                !!if (m == 6) then
                !!  write (iulog, *) c, j, m, 'k_tot h2o', k_tot * h2osoi_liqvol(c,j)
                !!end if

                ! calculate dissolution rate in mol m-3 s-1
                r_dissolve_vr(c,j,m) = ssa_dyn(c,j,m) * primary_mineral_vr(c,j,m) * k_tot

                !!if (m == 6) then
                !!  write (iulog, *) c, j, m, 'r_dissolve_vr', ssa_dyn(c,j,m), primary_mineral_vr(c,j,m), k_tot, r_dissolve_vr(c,j,m)
                !!end if

                !write (100+iam, *) c, j, m, 'r_dissolve_vr', r_dissolve_vr(c,j,m), k_tot, ssa_dyn(c,j,m), primary_mineral_vr(c,j,m)
              end if
            end if

            !-------------------------------------------------------------------
            ! Take the minimum between reaction rate and diffusion rate, following
            ! shrinking core model, simplified
            ! p362-366 of
            ! Levenspiel, O. (1972). Chapter 12. Fluid-Particle Reactions. In Chemical 
            ! Reaction Engineering (2nd ed., pp. 357–408). USA: John Wiley & Sons, Inc.
            ! 
            ! When the reaction rate is diffusion-rate controlled, then the key
            ! assumption is [H+] = 0 at the passivation layer-primary mineral interface
            !-------------------------------------------------------------------
            if (passivation_thickness(c,j,m) > 0._r8) then
              D_h_eff = D_h * passivation_phi / passivation_tau
              Jflux = D_h_eff * 10**(-soil_ph(c,j)) * 1e3_r8 / passivation_thickness(c,j,m)
              dNb_dt = ssa_dyn(c,j,m) * primary_mineral_vr(c,j,m) * Jflux / EWParamsInst%primary_stoi_proton(m)

              !!write (iulog, *) 'D_h_eff', c, j, m, D_h_eff, D_h, passivation_phi, passivation_tau
              !!write (iulog, *) 'Jflux', c, j, m, Jflux, 10**(-soil_ph(c,j)), forc_gra(c,m), passivation_thickness(c,j,m)
              !!write (iulog, *) 'dNb_dt', c, j, m, dNb_dt, ssa_dyn(c,m), primary_mineral_vr(c,j,m), EWParamsInst%primary_stoi_proton(m)
              !!write (iulog, *) 'r_dissolve_vr', c, j, m, r_dissolve_vr(c,j,m)

              r_dissolve_vr(c,j,m) = min(r_dissolve_vr(c,j,m), dNb_dt)
            end if
          end do
        end if

        ! Limit the dissolution rate to prevent primary mineral from going negative
        do m = 1,nminerals
          r_dissolve_vr(c,j,m) = min(r_dissolve_vr(c,j,m), primary_mineral_vr(c,j,m) / dt)
        end do

        ! Update the mineral and cation fluxes based on the reaction rates
        do m = 1,nminerals
          primary_dissolve_vr(c,j,m) = r_dissolve_vr(c,j,m) * EWParamsInst%primary_mass(m)
        end do

        ! do not consider the proton consumption by the mineral - it's supplied by constant
        ! delivery of CO2 and typically exceeds proton concentration in water
        !primary_proton_flux_vr(c,j) = 0._r8

        do icat = 1,ncations
          primary_cation_flux_vr(c,j,icat) = 0._r8
          do m = 1,nminerals
            primary_cation_flux_vr(c,j,icat) = primary_cation_flux_vr(c,j,icat) + &
              r_dissolve_vr(c,j,m) * EWParamsInst%primary_stoi_cations(m,icat) * EWParamsInst%cations_mass(icat)
          end do
        end do

        primary_h2o_flux_vr(c,j) = 0._r8
        do m = 1,nminerals
          primary_h2o_flux_vr(c,j) = primary_h2o_flux_vr(c,j) + & 
            r_dissolve_vr(c,j,m) * EWParamsInst%primary_stoi_h2o(m) * mass_h2o
        end do

        primary_silica_flux_vr(c,j) = 0._r8
        do m = 1,nminerals
          primary_silica_flux_vr(c,j) = primary_silica_flux_vr(c,j) + &
            r_dissolve_vr(c,j,m) * EWParamsInst%primary_stoi_sio2(m) * mass_sio2
        end do

        do m = 1,nminerals
          primary_residue_flux_vr(c,j,m) = EWParamsInst%primary_mass(m) + &
            EWParamsInst%primary_stoi_proton(m)*mass_h - & 
            EWParamsInst%primary_stoi_h2o(m)*mass_h2o - & 
            EWParamsInst%primary_stoi_sio2(m) * mass_sio2
          do icat = 1,ncations
            primary_residue_flux_vr(c,j,m) = primary_residue_flux_vr(c,j,m) - &
              EWParamsInst%primary_stoi_cations(m,icat) * EWParamsInst%cations_mass(icat)
          end do
          primary_residue_flux_vr(c,j,m) = primary_residue_flux_vr(c,j,m) * r_dissolve_vr(c,j,m)
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!! DOUBLE CHECK P release calculation
        primary_prelease_vr(c,j) = 0._r8
        do m = 1,nminerals
          primary_prelease_vr(c,j) = primary_prelease_vr(c,j) + &
              primary_dissolve_vr(c,j,m) * forc_pho(c)
        end do
      end do
    end do

    end associate
  end subroutine MineralPrimary

  !-----------------------------------------------------------------------
  subroutine MineralSecondary(bounds, num_soilc, filter_soilc, soilstate_vars)
    !
    ! !DESCRIPTION: 
    ! Calculate secondary mineral precipitation and dissolution fluxes. 
    ! 
    ! !USES:
    ! rgas = universal gas constant [= 8314.467 J/K/kmole]
    use elm_varcon       , only : spval, rgas, secspday, D_h
    use elm_varctl       , only : use_erw_verbose
    use elm_time_manager , only : get_step_size, get_curr_date
    use abortutils       , only : endrun
    use SharedParamsMod  , only : ParamsShareInst
    use timeinfoMod
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilstate_type)     , intent(in)    :: soilstate_vars

    !
    ! !LOCAL VARIABLES:
    integer  :: fc,c,j,g,nlevbed
    integer  :: m, isec, icat          ! indices
    real(r8) :: saturation_ratio, log_silica, log_carbonate
    real(r8) :: k_tot
    real(r8) :: dt 

    associate( &
         !
         ! soil pH and ionic states 
         !
         soil_ph                        => col_ms%soil_ph                 , & ! Input: [real(r8) (:,:)] calculated soil pH (1:nlevgrnd)
         bicarbonate_vr                 => col_ms%bicarbonate_vr          , & ! Output: [real(r8) (:,:)] calculated HCO3- concentration in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)
         carbonate_vr                   => col_ms%carbonate_vr            , & ! Output: [real(r8) (:,:)] calculated CO3 2- concentration in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)
         cation_vr                      => col_ms%cation_vr               , & ! Output [real(r8) (:,:,:)] cation mass in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd, 1:ncations)

         ! 
         ! Secondary mineral flux
         ! 
         secondary_mineral_vr           => col_ms%secondary_mineral_vr      , & ! Output [real(r8) (:,:,:)] secondary mineral mass in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd, 1:nminsecs)
         secondary_mineral_flux_vr      => col_mf%secondary_mineral_flux_vr , & ! Output [real(r8) (:,:,:) secondary mineral precipitated (g m-3 s-1) (1:nlevgrnd, 1:nminsecs)
         secondary_cation_flux_vr       => col_mf%secondary_cation_flux_vr  , & ! Output [real(r8) (:,:,:) cations consumed due to precipitation of secondary minerals (g m-3 s-1) (1:nlevgrnd, 1:ncations)
         secondary_silica_flux_vr       => col_mf%secondary_silica_flux_vr  , & ! Output [real(r8) (:,:) sio2 consumed due to precipitation of secondary minerals (g m-3 s-1) (1:nlevgrnd)
         r_precip_vr                    => col_mf%r_precip_vr               , & ! Output [real(r8) (:,:)] rate at which the precipitation of secondary mineral happens (mol m-3 s-1) (1:nlevgrnd, 1:nminsecs)

         !
         ! Other related
         !
         silica_vr                      => col_ms%silica_vr                , & ! Output [real(r8) (:,:)] silica mass in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)
         dz                             => col_pp%dz                       , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)
         tsoi                           => col_es%t_soisno                 , & ! Input: [real(r8) (:,:) ] soil temperature [K]
         h2osoi_vol                     => col_ws%h2osoi_vol               , & ! Input:  [real(r8) (:)] volumetric soil water content, ice + water (m3 m-3)
         h2osoi_liqvol                  => col_ws%h2osoi_liqvol            , & ! Input:  [real(r8) (:)] volumetric soil water content, liquid only (m3 m-3)
         nlev2bed                       => col_pp%nlevbed      & ! Input:  [integer  (:)   ]  number of layers to bedrock
    )

    dt      = real( get_step_size(), r8 )

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      g = col_pp%gridcell(c)
      !topo = col_pp%topounit(c)
      nlevbed = min(nlev2bed(c), nlevsoi)

      ! positive means going into the solution, negative means precipitating out
      secondary_cation_flux_vr(c,:,:) = 0._r8
      secondary_mineral_flux_vr(c,:,:)= 0._r8
      r_precip_vr(c,:,:)              = 0._r8
      secondary_silica_flux_vr(c,:)   = 0._r8

      do j = 1,nlevbed

        if (h2osoi_liqvol(c,j) > 1e-6) then

          ! General precipitation rate law
          ! 
          ! Marty, N. C. M., Claret, F., Lassin, A., Tremosa, J., Blanc, P., 
          ! Madé, B., et al. (2015). A database of dissolution and precipitation 
          ! rates for clay-rocks minerals. Applied Geochemistry, 55, 108–118. 
          ! https://doi.org/10.1016/j.apgeochem.2014.10.012
          !
          ! r [mol m-3 s-1] = k [mol m-2 s-1] * S [m2 m-3] * (\Omega^\theta - 1)^\eta

          do isec = 1,nminsecs
            if (isec == 1) then
              ! ---------------------------------------------------------------
              ! Calcite precipitation (Ca2+ is cation #1)
              ! CaCO3 +1.0000 H+  =  + 1.0000 Ca++ + 1.0000 HCO3- (llnl.dat)
              ! ---------------------------------------------------------------
              icat = 1
              saturation_ratio = &
                mass_to_mol(bicarbonate_vr(c,j), mass_hco3, h2osoi_vol(c,j)) * &
                mass_to_mol(cation_vr(c,j,icat), EWParamsInst%cations_mass(icat), h2osoi_vol(c,j)) * &
                10**(soil_ph(c,j) - EWParamsInst%log_keq_minsecs(icat))
            else if (isec == 2) then
              ! ---------------------------------------------------------------
              ! Kaolinite formation (Al3+ is cation #5)
              ! Al2Si2O5(OH)4 +6.0000 H+  =  + 2.0000 Al+++ + 2.0000 SiO2 + 5.0000 H2O (llnl.dat)
              ! ---------------------------------------------------------------
              icat = 5
              ! check silica concentration - if supersaturated, reduce to saturation point
              log_silica = mass_to_logmol(silica_vr(c,j), mass_sio2, h2osoi_vol(c,j))
              log_silica = min(log_silica, log_keq_sio2am)

              saturation_ratio = 10**(2*log_silica + &
                2 * mass_to_logmol(cation_vr(c,j,icat), EWParamsInst%cations_mass(icat), h2osoi_vol(c,j)) + &
                6 * soil_ph(c,j) - EWParamsInst%log_keq_minsecs(isec))
            else if (isec == 3) then
              ! ---------------------------------------------------------------
              ! Kaolinite formation (Al3+ is cation #5)
              ! Al(OH)3 +3.0000 H+  =  + 1.0000 Al+++ + 3.0000 H2O (llnl.dat)
              ! ---------------------------------------------------------------
              icat = 5
              saturation_ratio = &
                mass_to_mol(cation_vr(c,j,icat), EWParamsInst%cations_mass(icat), h2osoi_vol(c,j)) * &
                10**(3 * soil_ph(c,j) - EWParamsInst%log_keq_minsecs(isec))
            end if

            if (saturation_ratio > 1._r8) then
              if (cation_vr(c,j,icat) > 0._r8) then
                ! run the precipitation reaction

                !! No longer used for consistency with gibbsite
                !! r = \alpha * (\Omega - 1)
                !! \alpha = 9*1e-10 mol dm-3 (solution) s-1
                !! Kirk, G. J. D., Versteegen, A., Ritz, K. & Milodowski, A. E. A simple reactive-transport model of calcite precipitation in soils and other porous media. Geochimica et Cosmochimica Acta 165, 108–122 (2015).

                !! No longer used for consistency with gibbsite
                ! r [mol m-3 s-1] = A_{bulk} [m2 m-3] * k * (\Omega - 1)
                ! Perez-Fodich, A., & Derry, L. A. (2020). A model for germanium-silicon equilibrium fractionation in kaolinite. Geochimica et Cosmochimica Acta, 288, 199–213. https://doi.org/10.1016/j.gca.2020.07.046

                ! Reaction rate constant is
                ! k = k_precip[H2O] * exp[ - E/R * (1/T - 1/298.15)] + 
                !     k_precip[H+] * 10**(-pH) * exp[ - E/R * (1/T - 1/298.15)] + 
                !     k_precip[OH-] * 10**(pH - 14) * exp[ - E/R * (1/T - 1/298.15)]
                ! Marty et al. (2015). A database of dissolution and precipitation rates for clay-rocks minerals. Applied Geochemistry, 55, 108–118. https://doi.org/10.1016/j.apgeochem.2014.10.012
                k_tot = 0._r8
                if (EWParamsInst%k_precip_minsecs(isec,1) > -9000._r8) then
                  k_tot = k_tot + EWParamsInst%k_precip_minsecs(isec,1) * & 
                    10**(-soil_ph(c,j) * EWParamsInst%n_precip_minsecs(1)) * &
                    exp(-1e6_r8 * EWParamsInst%e_precip_minsecs(isec,1) / rgas * (1/tsoi(c,j) - 1/298.15_r8)) * &
                    (saturation_ratio - 1._r8)
                end if

                if (EWParamsInst%k_precip_minsecs(isec,2) > -9000._r8) then
                  k_tot = k_tot + EWParamsInst%k_precip_minsecs(isec,2) * &
                    exp(-1e6_r8 * EWParamsInst%e_precip_minsecs(isec,2) / rgas * (1/tsoi(c,j) - 1/298.15_r8)) * &
                    (saturation_ratio**EWParamsInst%ph2o_precip_minsecs(isec) - 1._r8)**EWParamsInst%qh2o_precip_minsecs(isec)
                end if

                if (EWParamsInst%k_precip_minsecs(isec,3) > -9000._r8) then
                  k_tot = k_tot + EWParamsInst%k_precip_minsecs(isec,3) * & 
                    10**((soil_ph(c,j) - 14) * EWParamsInst%n_precip_minsecs(3)) * &
                    exp(-1e6_r8 * EWParamsInst%e_precip_minsecs(isec,3) / rgas * (1/tsoi(c,j) - 1/298.15_r8)) * &
                    (saturation_ratio - 1._r8)
                end if

                r_precip_vr(c,j,isec) = k_tot * secondary_mineral_vr(c,j,isec) * EWParamsInst%ssa_minsecs(isec)

                ! convert to mol kg-1 water s-1
                r_precip_vr(c,j,isec) = r_precip_vr(c,j,isec) / h2osoi_liqvol(c,j) * 1e-3_r8

                ! limit the precipitation rate by the reactant's concentration
                if (isec == 1) then
                  r_precip_vr(c,j,isec) = min( &
                    mass_to_mol(bicarbonate_vr(c,j), mass_hco3, h2osoi_vol(c,j)) / dt, &
                    mass_to_mol(cation_vr(c,j,icat), EWParamsInst%cations_mass(icat), h2osoi_vol(c,j)) / dt, &
                    r_precip_vr(c,j,isec) )
                else if (isec == 2) then
                  r_precip_vr(c,j,isec) = min( &
                    10**log_silica / 2._r8 / dt, &
                    mass_to_mol(cation_vr(c,j,icat), EWParamsInst%cations_mass(icat), h2osoi_vol(c,j)) / 2._r8 / dt, &
                    r_precip_vr(c,j,isec) )
                else if (isec == 3) then
                  ! limit the precipitation rate by the reactant's concentration
                  r_precip_vr(c,j,isec) = min( &
                    mass_to_mol(cation_vr(c,j,icat), EWParamsInst%cations_mass(icat), h2osoi_vol(c,j)) / dt, &
                    r_precip_vr(c,j,isec) )
                else
                  call endrun('MineralSecondary: isec > 3, this is out of range')
                end if

                ! convert back from mol kg-1 water s-1 to g m-3 s-1
                ! (reaction is in liquid part only)
                r_precip_vr(c,j,isec) = r_precip_vr(c,j,isec) * h2osoi_liqvol(c,j) * 1e3_r8

                ! switch sign to mean negative = out of solution
                r_precip_vr(c,j,isec) = - r_precip_vr(c,j,isec)
              end if

            else
              ! run the dissolution reaction

              ! Reaction rate constant is
              ! k = k_precip[H2O] * exp[ - E/R * (1/T - 1/298.15)] + 
              !     k_precip[H+] * 10**(-pH) * exp[ - E/R * (1/T - 1/298.15)] + 
              !     k_precip[OH-] * 10**(pH - 14) * exp[ - E/R * (1/T - 1/298.15)]
              ! Marty et al. (2015). A database of dissolution and precipitation rates for clay-rocks minerals. Applied Geochemistry, 55, 108–118. https://doi.org/10.1016/j.apgeochem.2014.10.012
              k_tot = 0._r8
              if (EWParamsInst%k_dissolv_minsecs(isec,1) > -9000._r8) then
                k_tot = k_tot + EWParamsInst%k_dissolv_minsecs(isec,1) * & 
                  10**(-soil_ph(c,j) * EWParamsInst%n_dissolv_minsecs(isec,1)) * &
                  exp(-1e6_r8 * EWParamsInst%e_dissolv_minsecs(isec,1) / rgas * (1/tsoi(c,j) - 1/298.15_r8)) * &
                  (1._r8 - saturation_ratio)
              end if

              if (EWParamsInst%k_dissolv_minsecs(isec,2) > -9000._r8) then
                k_tot = k_tot + EWParamsInst%k_dissolv_minsecs(isec,2) * &
                  exp(-1e6_r8 * EWParamsInst%e_dissolv_minsecs(isec,2) / rgas * (1/tsoi(c,j) - 1/298.15_r8)) * &
                  (1._r8 - saturation_ratio)
              end if

              if (EWParamsInst%k_dissolv_minsecs(isec,3) > -9000._r8) then
                k_tot = k_tot + EWParamsInst%k_dissolv_minsecs(isec,3) * & 
                  10**((soil_ph(c,j) - 14) * EWParamsInst%n_dissolv_minsecs(isec,3)) * &
                  exp(-1e6_r8 * EWParamsInst%e_dissolv_minsecs(isec,3) / rgas * (1/tsoi(c,j) - 1/298.15_r8)) * &
                  (1._r8 - saturation_ratio)
              end if

              r_precip_vr(c,j,isec) = k_tot * secondary_mineral_vr(c,j,isec) * EWParamsInst%ssa_minsecs(isec)

              ! limit the precipitation rate by the reactant's concentration
              r_precip_vr(c,j,isec) = min(secondary_mineral_vr(c,j,isec) / dt / EWParamsInst%minsecs_mass(isec), &
                                          r_precip_vr(c,j,isec))
            end if

            ! update the fluxes for operative sec. minerals/cations
            secondary_cation_flux_vr(c,j,icat) = r_precip_vr(c,j,isec) * EWParamsInst%cations_mass(icat)
            secondary_mineral_flux_vr(c,j,isec) = r_precip_vr(c,j,isec) * EWParamsInst%minsecs_mass(isec)
            if (isec == 2) then
              secondary_silica_flux_vr(c,j) = r_precip_vr(c,j,isec) * mass_sio2
            end if

          end do ! isec
        end if ! h2osoi_liqvol > 1e-6
      end do ! soil layer
    end do ! end of the soil column loop

    end associate

  end subroutine MineralSecondary


  !-----------------------------------------------------------------------
  subroutine MineralPassivation(bounds, num_soilc, dt, filter_soilc)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the rate at which secondary minerals
    ! accumulate in the passivation layer
    !
    ! !USES:
    use elm_varcon           , only: vol_sio2
    use ewutils              , only: get_r
    !$acc routine seq
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    real(r8)                 , intent(in)    :: dt              ! time step size (s)
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    !
    ! !LOCAL VARIABLES:
    integer  :: fc,c,j,g,nlevbed,m
    real(r8) :: tot_surf               ! total surface area of primary minerals in the soil column (m2 m-3)
    real(r8) :: mass0, mass            ! mineral mass
    real(r8) :: r0, r                  ! mineral radius

    associate( &
         nlev2bed                       => col_pp%nlevbed                   , & ! Input:  [integer  (:)   ] number of layers to bedrock
         primary_mineral_vr             => col_ms%primary_mineral_vr        , & ! Input:  [real(r8) (:,:,:) ] primary mineral content in the soil column (g m-3) (1:ncol, 1:nlevgrnd, 1:nminerals)
         forc_gra                       => col_ew%forc_gra                  , & ! Input:  [real(r8) (:,:)] grain size (1:nminerals) (um diameter)
         ssa_dyn                        => col_ms%ssa_dyn                   , & ! Input:  [real(r8) (:,:,:)] specific surface area of primary minerals (m2 g-1) (1:ncol, 1:nlevgrnd, 1:nminerals)
         r_dissolve_vr                  => col_mf%r_dissolve_vr             , & ! Input [real(r8) (:,:)] rate at which primary mineral dissolves (mol m-3 s-1) (1:nlevgrnd, 1:nminerals)
         passivation_rate               => col_mf%passivation_rate          , & ! Output: [real(r8) (:)] rate at which the passivation layer accumulates (m s-1) (1:nlevgrnd)
         passivation_thickness          => col_ms%passivation_thickness     & ! Output: [real(r8) (:)] thickness of the passivation layer (m) (1:nlevgrnd)
    )

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      g = col_pp%gridcell(c)
      !topo = col_pp%topounit(c)
      nlevbed = min(nlev2bed(c), nlevsoi)

      do j = 1,nlevbed
        do m = 1,nminerals
          if (primary_mineral_vr(c,j,m) > 0._r8) then
            passivation_rate(c,j,m) = r_dissolve_vr(c,j,m) * EWParamsInst%primary_mass(m) / (2.9e6_r8) / (primary_mineral_vr(c,j,m) * ssa_dyn(c,j,m))

            !!write (iulog, *) 'passivation', c, j, m, passivation_rate(c,j,m), r_dissolve_vr(c,j,m), EWParamsInst%primary_mass(m), primary_mineral_vr(c,j,m), ssa_dyn(c,m)

            passivation_rate(c,j,m) = min(passivation_rate(c,j,m), (forc_gra(c,m)*1e-6_r8 - passivation_thickness(c,j,m)) / dt)
          else
            passivation_rate(c,j,m) = 0._r8
          end if

        end do
      end do
    end do

    end associate
  end subroutine MineralPassivation


  !-----------------------------------------------------------------------
  subroutine MineralVerticalMovement(bounds, num_soilc, filter_soilc, dt, soilstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the boundary conditions of
    ! soil ions (H+, HCO3-, CO3--, cations) caused by vertical movement of soil water
    !
    ! !USES:
    !$acc routine seq
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    real(r8)                 , intent(in)    :: dt              ! radiation time step (seconds)

    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,fc,g,l,t
    integer  :: icat                                   ! indices
    integer  :: nlevbed                                ! number of layers to bedrock
    real(r8) :: frac_thickness                         ! deal with the fractional layer between last layer and max allowed depth
    real(r8) :: tot_water                              ! total column liquid water (kg water/m2)
    real(r8) :: surface_water                          ! liquid water to shallow surface depth (kg water/m2)
    real(r8) :: drain_tot                              ! total drainage flux (mm H2O /s)
    real(r8), parameter :: depth_runoff_Mloss = 0.05   ! (m) depth over which runoff mixes with soil water for ions loss to runoff; same as nitrogen runoff depth
    real(r8) :: co2_atm                                ! CO2 partial pressure in atm
    real(r8) :: rain_carbonate, rain_bicarbonate       ! surface boundary condition (g m-3 H2O)
    real(r8) :: rain_cations(1:ncations)               ! surface boundary condition (g m-3 H2O)
    real(r8) :: sourcesink_zero(1:nlevsoi)             ! src/sink term (g m-3 soil s-1)
    real(r8) :: sourcesink_cations(1:nlevsoi,1:ncations) ! src/sink term (g m-3 soil s-1)
    real(r8) :: adv_water(1:nlevsoi+1)                 ! m H2O / s, negative downward
    real(r8) :: dcation_dt(1:nlevsoi, 1:ncations)      ! cation concentration change rate, g m-3 s-1
    real(r8) :: dhco3_dt(1:nlevsoi)                    ! HCO3- concentration change rate, g m-3 s-1
    real(r8) :: dco3_dt(1:nlevsoi)                     ! CO3-- concentration change rate, g m-3 s-1

    !-----------------------------------------------------------------------

    associate( &
         qin                            => col_wf%qin                     , & ! Input: [real(r8) (:,:) ] flux of water into soil layer [mm h2o/s]
         qout                           => col_wf%qout                    , & ! Input: [real(r8) (:,:) ] flux of water out of soil layer [mm h2o/s]
         h2osoi_vol                     => col_ws%h2osoi_vol              , & ! Input:  [real(r8) (:)] volumetric soil water content, ice + water (m3 m-3)

         nlev2bed                       => col_pp%nlevbed                 , & ! Input:  [integer  (:)   ]  number of layers to bedrock
         dz                             => col_pp%dz                      , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)

         rain_ph                        => col_ew%rain_ph                 , & ! Output: [real(r8) (:)] pH of rain water
         rain_chem                      => col_ew%rain_chem               , & ! Output: [real(r8) (:,:)] cation concentration in rain water (excluding H+) (g m-3 rain water) (1:ncations)

         annavg_qin_col                 => col_wf%annavg_qin_col          , & ! Input:  [real(r8) (:,:) ]  annual average rate of vertical water movement (mm H2O / s)
         mixing_fraction                => col_mf%mixing_fraction         , & ! Input:  [real(r8) (:,:) ] fraction of vertical water flow that participated in mineral reactions (-)

         cation_vr                      => col_ms%cation_vr               , & ! Output [real(r8) (:,:,:)] cation mass in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd, 1:ncations)
         bicarbonate_vr                 => col_ms%bicarbonate_vr          , & ! Output: [real(r8) (:,:)] calculated HCO3- concentration in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)
         carbonate_vr                   => col_ms%carbonate_vr            , & ! Output: [real(r8) (:,:)] calculated CO3 2- concentration in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)

         cation_infl_vr                 => col_mf%cation_infl_vr          , & ! Output: [real(r8) (:,:,:)] cation flux carried from infiltration above (g m-3 soil s-1 [not water]) (1:nlevgrnd, 1:ncations)
         cation_oufl_vr                 => col_mf%cation_oufl_vr          , & ! Output: [real(r8) (:,:,:)] cation flux carried away by infiltration (g m-3 soil s-1 [not water]) (1:nlevgrnd, 1:ncations)
         cation_uptake_vr               => col_mf%cation_uptake_vr        , & ! Output: [real(r8) (:,:,:)] cation flux uptaken by plants (g m-3 soil s-1 [not water]) (1:nlevgrnd, 1:ncations)
         cec_cation_flux_vr             => col_mf%cec_cation_flux_vr      , & ! Output: [real(r8) (:,:,:)] rate at which adsorbed cation is released into water (negative for adsorption into soil) (vertically resolved) (1:nlevgrnd, 1:ncations) (g m-3 s-1)
         cec_cation_flux2_vr            => col_mf%cec_cation_flux2_vr     , & ! Output: [real(r8) (:,:,:)] rate at which adsorbed cation is released into water due to change in total cation exchange capacity (always >= 0) (vertically resolved) (1:nlevgrnd, 1:ncations) (g m-3 s-1)

         bicarbonate_drainage           => col_mf%bicarbonate_drainage    , & ! Output: [real(r8) (:,:)] bottom drainage of HCO3- due to vertical infiltration (positive for increase) (g m-3 soil s-1 [not water]) (1:nlevgrnd)
         carbonate_drainage             => col_mf%carbonate_drainage      , & ! Output: [real(r8) (:,:)] bottom drainage of CO3-- due to vertical infiltration (positive for increase) (g m-3 soil s-1 [not water]) (1:nlevgrnd)

         primary_cation_flux_vr         => col_mf%primary_cation_flux_vr  , & ! Output [real(r8) (:,:,:) cations produced due to all the dissolution reactions (g m-3 s-1) (1:nlevgrnd, 1:ncations)
         background_flux_vr             => col_mf%background_flux_vr      , & ! Output: [real(r8) (:)] background weathering rate (g m-3 s-1)
         secondary_cation_flux_vr       => col_mf%secondary_cation_flux_vr & ! Output [real(r8) (:,:,:) cations consumed due to precipitation of secondary minerals (g m-3 s-1) (1:nlevgrnd, 1:ncations)
    )

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      t = col_pp%topounit(c)
      g = col_pp%gridcell(c)
      nlevbed = min(nlev2bed(c), nlevsoi)

      co2_atm = top_as%pco2bot(t) / 101325

      !------------------------------------------------------------------------------
      ! Collect the water flow boundary conditions (g m-3)
      !------------------------------------------------------------------------------
      ! calculate from rain pH to g/m3 rain
      rain_bicarbonate = mol_to_mass(ph_to_hco3(rain_ph(c), co2_atm), mass_hco3, 1._r8)
      rain_carbonate = mol_to_mass(hco3_to_co3(ph_to_hco3(rain_ph(c), co2_atm), &
                                               rain_ph(c)), mass_co3, 1._r8)
      do icat = 1,ncations
        ! mg/L rain => m/m3 rain
        rain_cations(icat) = rain_chem(c,icat)
      end do

      !------------------------------------------------------------------------------
      ! Calculate the net source sink due to plant uptake, background weathering, 
      ! mineral reactions, and cation exchange (g m-3 s-1)
      !------------------------------------------------------------------------------
      do j = 1,nlevbed
        !------------------------------------------------------------------------------
        ! uptake by vegetation - set to zero, assuming litterfall balances out uptake
        !------------------------------------------------------------------------------
        do icat = 1,ncations
          cation_uptake_vr(c,j,icat) = 0._r8
        end do

        do icat = 1,ncations
          sourcesink_cations(j,icat) = background_flux_vr(c,j,icat) + & 
            primary_cation_flux_vr(c,j,icat) + cec_cation_flux_vr(c,j,icat) + &
            cec_cation_flux2_vr(c,j,icat) - secondary_cation_flux_vr(c,j,icat) - &
            cation_uptake_vr(c,j,icat)
        end do

        ! for HCO3-, CO3--, calculate the pure advection diffusion at equilibrium concentrations
        sourcesink_zero(j) = 0._r8
      end do

      !------------------------------------------------------------------------------
      ! The fraction of incoming water that mixes with soil solution to carry away
      ! cations
      !------------------------------------------------------------------------------
      ! use placeholder => 1 for now
      mixing_fraction(c,1:nlevbed+1) = 1._r8

      !------------------------------------------------------------------------------
      ! Calculate the vertical transport
      !------------------------------------------------------------------------------
      do j = 1,nlevbed
        ! note the flux rate is positive downward
        adv_water(j) = 1.0e-3_r8 * qin(c,j) * mixing_fraction(c,j)
      end do
      adv_water(nlevbed + 1) = 1.0e-3_r8 * qin(c,nlevbed+1) * mixing_fraction(c,nlevbed+1)

      do icat = 1,ncations
        !!write (iulog, *) 'cation_vr', c, icat, cation_vr(c, 1:nlevsoi, icat)
        !!write (iulog, *) 'rain_cations', c, icat, rain_cations(icat)
        !!write (iulog, *) 'adv_water', c, icat, adv_water(1:nlevsoi+1)
        !!write (iulog, *) 'h2osoi_vol', c, icat, h2osoi_vol(c, 1:nlevsoi)
        !!write (iulog, *) 'watsat_col', c, icat, soilstate_vars%watsat_col(c,1:nlevsoi)
        !!write (iulog, *) 'sourcesink_cations', c, icat, sourcesink_cations(1:nlevsoi, icat)
        !!write (iulog, *) 'cations_diffusivity', c, icat, EWParamsInst%cations_diffusivity(icat)
        !!write (iulog, *) 'dz', c, icat, dt, nlevbed, dz(c,1:nlevsoi)

        call advection_diffusion(cation_vr(c, 1:nlevsoi, icat), rain_cations(icat), &
                                 adv_water(1:nlevsoi+1), h2osoi_vol(c, 1:nlevsoi), &
                                 soilstate_vars%watsat_col(c,1:nlevsoi), &
                                 sourcesink_cations(1:nlevsoi, icat), &
                                 EWParamsInst%cations_diffusivity(icat), &
                                 dt, dz(c,1:nlevsoi), nlevbed, dcation_dt(1:nlevsoi, icat))

        !!write (iulog, *) 'dcation_dt', c, j, icat, dcation_dt(1:nlevsoi, icat)

      end do

      call advection_diffusion(bicarbonate_vr(c, 1:nlevsoi), rain_bicarbonate, &
                               adv_water(1:nlevsoi+1), h2osoi_vol(c, 1:nlevsoi), &
                               soilstate_vars%watsat_col(c,1:nlevsoi), &
                               sourcesink_zero(1:nlevsoi), &
                               EWParamsInst%bicarbonate_diffusivity, &
                               dt, dz(c,1:nlevsoi), nlevbed, dhco3_dt(1:nlevsoi))

      call advection_diffusion(carbonate_vr(c, 1:nlevsoi), rain_carbonate, &
                               adv_water(1:nlevsoi+1), h2osoi_vol(c, 1:nlevsoi), &
                               soilstate_vars%watsat_col(c,1:nlevsoi), &
                               sourcesink_zero(1:nlevsoi), &
                               EWParamsInst%carbonate_diffusivity, &
                               dt, dz(c,1:nlevsoi), nlevbed, dco3_dt(1:nlevsoi))

      !------------------------------------------------------------------------------
      ! Update the cation concentrations using the vertical transport
      !------------------------------------------------------------------------------
      do j = 1, nlevbed
        do icat = 1, ncations
          cation_infl_vr(c,j,icat) = dcation_dt(j, icat) ! - sourcesink_cations(j,icat)
          cation_oufl_vr(c,j,icat) = 0._r8
        end do
      end do

      !DEBUG
      !do icat = 1,ncations
      !  write (iulog, *), icat, 'cation_vr', cation_vr(c, 1:nlevsoi, icat)
      !end do
      !do icat = 1,ncations
      !  write (iulog, *), icat, 'sourcesink', sourcesink_cations(1:nlevsoi, icat) * dt
      !end do
      !do icat = 1,ncations
      !  write (iulog, *), icat, 'dcation_dt', dcation_dt(1:nlevsoi, icat) * dt
      !end do

      !------------------------------------------------------------------------------
      ! Calculate the bottom layer drainage of HCO3- and CO3--
      ! 
      ! Rain - Outflow = ΔStorage
      ! => Outflow = Rain - ΔStorage
      !------------------------------------------------------------------------------
      bicarbonate_drainage(c) = rain_bicarbonate * 1.0e-3_r8 * qin(c,1)
      carbonate_drainage(c) = rain_carbonate * 1.0e-3_r8 * qin(c,1)
      do j = 1, nlevbed
        bicarbonate_drainage(c) = bicarbonate_drainage(c) - dhco3_dt(j) * dz(c,j)
        carbonate_drainage(c) = carbonate_drainage(c) - dco3_dt(j) * dz(c,j)
      end do

    end do ! end soil column loop

    end associate
  end subroutine MineralVerticalMovement


  !-----------------------------------------------------------------------
  subroutine MineralLeaching(bounds, num_soilc, filter_soilc, dt)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the boundary conditions of
    ! soil ions (H+, cations) caused by subsurface & surface runoff
    !
    ! !USES:
    !$acc routine seq
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    real(r8)                 , intent(in)    :: dt              ! radiation time step (seconds)

    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,fc,g,l
    integer  :: icat                                   ! indices
    integer  :: nlevbed                                ! number of layers to bedrock
    integer  :: resolve_drain                          ! 1 = vertically resolved drainage; 0 = no
    real(r8) :: frac_thickness                         ! deal with the fractional layer between last layer and max allowed depth
    real(r8) :: tot_water                              ! total column liquid water (kg water/m2)
    real(r8) :: surface_water                          ! liquid water to shallow surface depth (kg water/m2)
    real(r8) :: drain_tot                              ! total drainage flux (mm H2O /s)
    real(r8) :: temp_drain_pct, temp_surf_pct          ! percentage loss per second
    real(r8) :: temp_tot_pct, tot_cation_flush
    real(r8), parameter :: depth_runoff_Mloss = 0.05   ! (m) depth over which runoff mixes with soil water for ions loss to runoff; same as nitrogen runoff depth

    !-----------------------------------------------------------------------

    associate( &
         h2osoi_liq             => col_ws%h2osoi_liq                      , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
         qflx_drain             => col_wf%qflx_drain                      , & ! Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)                    
         qflx_surf              => col_wf%qflx_surf                       , & ! Input:  [real(r8) (:)   ]  surface runoff (mm H2O /s)

         qout_external          => col_wf%qout_external                   , & ! Input:  [real(r8) (:,:) ]  vertically resolved sub-surface runoff (mm H2O /s)                    

         cation_leached_vr      => col_mf%cation_leached_vr               , & ! Output: [real(r8) (:,:,:) ]  rate of cation leaching (g m-3 s-1)
         cation_runoff_vr       => col_mf%cation_runoff_vr                , & ! Output: [real(r8) (:,:,:) ]  rate of cation loss with runoff (g m-3 s-1)

         bicarbonate_vr         => col_ms%bicarbonate_vr                  , & ! Output: [real(r8) (:,:)] calculated HCO3- concentration in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)
         bicarbonate_leached_vr => col_mf%bicarbonate_leached_vr          , & ! Output: [real(r8) (:,:) ] rate of HCO3- leaching (g m-3 s-1) (1:nlevgrnd)

         carbonate_vr           => col_ms%carbonate_vr                    , & ! Output: [real(r8) (:,:)] calculated CO3 2- concentration in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)
         carbonate_leached_vr   => col_mf%carbonate_leached_vr            , & ! Output: [real(r8) (:,:) ] rate of CO3-- leaching (g m-3 s-1) (1:nlevgrnd)

         dz                             => col_pp%dz                      , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)
         nlev2bed                       => col_pp%nlevbed                 , & ! Input:  [integer  (:)   ]  number of layers to bedrock

         rain_ph                        => col_ew%rain_ph                 , & ! Output: [real(r8) (:)] pH of rain water
         rain_chem                      => col_ew%rain_chem               , & ! Output: [real(r8) (:,:)] cation concentration in rain water (excluding H+) (g m-3 rain water) (1:ncations)

         proton_vr                      => col_ms%proton_vr               , & ! Input: calculated soil H+ concentration in soil water each soil layer (1:nlevgrnd) (g m-3 soil [not water])
         cation_vr                      => col_ms%cation_vr               , & ! Output [real(r8) (:,:,:)] cation mass in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd, 1:ncations)

         qflx_rootsoi_col               =>    col_wf%qflx_rootsoi         & ! Input: [real(r8) (:,:) ]  vegetation/soil water exchange (mm H2O/s) (+ = to atm)
    )

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      g = col_pp%gridcell(c)
      nlevbed = min(nlev2bed(c), nlevsoi)

      resolve_drain = 1 ! use vertically resolved subsurface runoff QOUT_EXTERNAL

      !------------------------------------------------------------------------------
      ! Leaching (subsurface runoff) and surface runoff losses
      !------------------------------------------------------------------------------

      ! assume zero mixing between surface runoff and soil pore water
      !! for runoff calculation; calculate total water to a given depth
      !surface_water = 0._r8
      !do j = 1,nlevbed
      !  if ( zisoi(j) <= depth_runoff_Mloss)  then
      !      surface_water = surface_water + h2osoi_liq(c,j)
      !  elseif ( zisoi(j-1) < depth_runoff_Mloss)  then
      !      frac_thickness = (depth_runoff_Mloss - zisoi(j-1)) / dz(c,j)
      !      surface_water = surface_water + h2osoi_liq(c,j) * frac_thickness
      !  end if
      !end do
      !! (qflx_surf / tot_water) is the fraction water lost per second
      !temp_surf_pct = qflx_surf(c) / surface_water
      temp_surf_pct = 0._r8

      if (resolve_drain == 0) then
        ! calculate the total soil water
        tot_water = 0._r8
        do j = 1,nlevbed
          tot_water = tot_water + h2osoi_liq(c,j)
        end do
        ! (drain_tot / tot_water) is the fraction water lost per second
        temp_drain_pct = qflx_drain(c) / tot_water
      end if

      do j = 1,nlevbed
        if (h2osoi_liq(c,j) > 0._r8) then
          ! use the analytical solution if the flushing rate is too large
          ! need to calculate the total leaching and runoff flux, otherwise
          ! they can be individually small, but still in total too large

          ! calculate the loss from surface runoff, assuming a shallow mixing of surface waters into soil and removal based on runoff
          ! &
          ! calculate the leaching flux as a function of the dissolved
          ! concentration (g cation/kg water) and the sub-surface drainage flux

          if (resolve_drain == 1) then
            temp_drain_pct = - qout_external(c,j) / h2osoi_liq(c,j)
          end if

          if ( zisoi(j-1) < depth_runoff_Mloss ) then

            if ( zisoi(j) <= depth_runoff_Mloss )  then
              frac_thickness = 1._r8
            else
              frac_thickness = (depth_runoff_Mloss - zisoi(j-1)) / dz(c,j)
            end if
            temp_tot_pct = temp_surf_pct*frac_thickness + temp_drain_pct

            if ( temp_tot_pct*dt > 0.1_r8 ) then
              do icat = 1,ncations
                tot_cation_flush = cation_vr(c,j,icat)*(1._r8-exp(-temp_tot_pct*dt))/dt
                cation_leached_vr(c,j,icat) = tot_cation_flush * temp_drain_pct / temp_tot_pct
                cation_runoff_vr(c,j,icat) = tot_cation_flush - cation_leached_vr(c,j,icat)
              end do
            else
              do icat = 1,ncations
                cation_leached_vr(c,j,icat) = cation_vr(c,j,icat) * temp_drain_pct
                cation_runoff_vr(c,j,icat) = cation_vr(c,j,icat) * temp_surf_pct * frac_thickness
              end do
            end if

          else

            cation_runoff_vr(c,j,1:ncations) = 0._r8
            if (temp_drain_pct*dt > 0.1_r8) then
              do icat = 1,ncations
                cation_leached_vr(c,j,icat) = cation_vr(c,j,icat)*(1._r8-exp(-temp_drain_pct*dt))/dt
              end do
            else
              do icat = 1,ncations
                cation_leached_vr(c,j,icat) = cation_vr(c,j,icat) * temp_drain_pct
              end do
            end if

          end if

          ! calculate the leaching flux as a function of the HCO3-- and CO3-- concentration and
          ! the sub-surface drainage flux
          if (temp_drain_pct*dt > 0.1_r8) then
            bicarbonate_leached_vr(c,j) = bicarbonate_vr(c,j) * (1._r8 - exp(-temp_drain_pct*dt))/dt
            carbonate_leached_vr(c,j) = carbonate_vr(c,j) * (1._r8 - exp(-temp_drain_pct*dt))/dt
          else
            bicarbonate_leached_vr(c,j) = bicarbonate_vr(c,j) * temp_drain_pct
            carbonate_leached_vr(c,j) = carbonate_vr(c,j) * temp_drain_pct
          end if

        else
          do icat = 1,ncations
            cation_leached_vr(c,j,icat) = 0._r8
            cation_runoff_vr(c,j,icat) = 0._r8
          end do

          bicarbonate_leached_vr(c,j) = 0._r8
          carbonate_leached_vr(c,j) = 0._r8
        end if

      end do ! end soil level loop

    end do ! end soil column loop

    end associate
  end subroutine MineralLeaching

  !-----------------------------------------------------------------------
  subroutine MineralEquilibria(bounds, num_soilc, filter_soilc, soilstate_vars)
    !
    ! !DESCRIPTION: 
    ! Calculate the dynamic pH value from the following set of equations
    ! 
    ! eq1 = sp.Eq(h * hco3 / co2_atm, 10**(-7.8136))
    ! eq2 = sp.Eq(h * co3 / hco3, 10**(-10.3288))
    ! eq3 = sp.Eq(h * oh, 1e-14)
    ! eq4 = sp.Eq(h / beta_h * (beta1 / ca)**(1/valence_Ca2), kex1) # 10**(3.4*(1-beta_h)) *  
    ! eq5 = sp.Eq(h / beta_h * (beta2 / mg)**(1/valence_Mg2), kex2) # 10**(3.4*(1-beta_h)) *  
    ! eq6 = sp.Eq(h / beta_h * (beta3 / na)**(1/valence_Na), kex3) # 10**(3.4*(1-beta_h)) * 
    ! eq7 = sp.Eq(h / beta_h * (beta4 / k)**(1/valence_K), kex4) # 10**(3.4*(1-beta_h)) * 
    ! eq8 = sp.Eq(h / beta_h * (beta5 / al)**(1/valence_Al3), kex5) # 10**(3.4*(1-beta_h)) * 
    ! eq9 = sp.Eq(h - oh - hco3 - 2*co3 + 2*ca + 2*mg + na + k + 3*al, b0)
    ! 
    ! !USES:
    use elm_time_manager , only : get_step_size, get_curr_date
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilstate_type)     , intent(in)    :: soilstate_vars

    !
    ! !LOCAL VARIABLES:
    integer  :: fc,c,t,j,nlevbed
    integer  :: icat                   ! indices
    real(r8) :: co2_atm                ! CO2 partial pressure in atm
    real(r8) :: cece(1:ncations)       ! temporary container (meq 100g-1 soil)
    real(r8) :: beta_list(1:ncations)  ! temporary container for cece/cec_tot
    real(r8) :: beta_h                 ! temporary container for ceca/cec_tot
    real(r8) :: keq_list(1:ncations)   ! temporary container for exchange coefficients between H+ and cations
    real(r8) :: conc(1:ncations)       ! temporary container for cation concentration (mol/kg)
    real(r8) :: prev_pH                ! soil pH before solving the CEC balance
    real(r8) :: oc_frac                ! fraction of organic carbon
    real(r8) :: dt
    real(r8) :: beta_for_release(1:ncations), beta_h_for_release

    associate( &
        net_charge_vr                       => col_ms%net_charge_vr           , & ! Input:  [real(r8) (:,:)] net charge of the tracked ions in the soil solution system, constant over time (1:nlevgrnd) (mol kg-1)
        nlev2bed                            => col_pp%nlevbed                 , & ! Input:  [integer (:)    ]  number of layers to bedrock

        soil_ph                             => col_ms%soil_ph                 , & ! Input:  [real(r8) (:,:)] calculated soil pH (1:nlevgrnd)
        proton_vr                           => col_ms%proton_vr               , & ! Input: [real (r8) (:,:)] calculated soil H+ concentration in soil water each soil layer (1:nlevgrnd) (g m-3 soil [not water])
        cation_vr                           => col_ms%cation_vr               , & ! Input: [real(r8) (:,:,:)] cation concentration in soil water in each soil layer (1:nlevgrnd,1:ncations) (g m-3 soil [not water])
        cec_cation_vr                       => col_ms%cec_cation_vr           , & ! Input: [real(r8) (:,:,:)] adsorbed cation concentration each soil layer (1:nlevgrnd,1:ncations) (g m-3 soil [not dry soil])

        cect_dyn                            => col_ms%cect_dyn                , & ! Input:  [real(r8) (:,:)] pH-dependent total cation exchange capacity (1:nlevgrnd)
        cect_delta                          => col_mf%cect_delta              , & ! Input:  [real(r8) (:,:)] pH-dependent change in cation exchange capacity (1:nlevgrnd)
        cece_delta                          => col_mf%cece_delta              , & ! Input:  [real(r8) (:,:)] pH-dependent change in individual cations occupied exchange capacity (1:nlevgrnd)

        h2osoi_vol                          => col_ws%h2osoi_vol              , & ! Input: [real(r8) (:,:)] soil water volume, liquid + ice (m3 m-3)
        h2osoi_liqvol                       => col_ws%h2osoi_liqvol           , & ! Input: [real(r8) (:,:)] soil water volume, liquid (m3 m-3)

        bicarbonate_vr                      => col_ms%bicarbonate_vr          , & ! Output: [real(r8) (:,:)] calculated HCO3- concentration in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)
        carbonate_vr                        => col_ms%carbonate_vr            , & ! Output: [real(r8) (:,:)] calculated CO3 2- concentration in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)

        cec_cation_flux_vr                  => col_mf%cec_cation_flux_vr      , & ! Output: [real(r8) (:,:,:)] rate at which adsorbed cation is released into water (negative for adsorption into soil) (vertically resolved) (1:nlevgrnd, 1:ncations) (g m-3 s-1)
        cec_cation_flux2_vr                 => col_mf%cec_cation_flux2_vr      & ! Output: [real(r8) (:,:,:)] rate at which adsorbed cation is released into water due to change in total CEC (zero if cect increases, positive if cect decreases) (vertically resolved) (1:nlevgrnd, 1:ncations) (g m-3 s-1)
    )

    dt      = real( get_step_size(), r8 )

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      t = col_pp%topounit(c)
      nlevbed = min(nlev2bed(c), nlevsoi)

      co2_atm = top_as%pco2bot(t) / 101325

      do j = 1,nlevbed

        ! before pH change
        prev_pH = soil_ph(c,j)

        ! use grid search to find the pH
        do icat = 1,ncations
          cece(icat) = mass_to_meq(cec_cation_vr(c,j,icat), EWParamsInst%cations_valence(icat), &
                                   EWParamsInst%cations_mass(icat), soilstate_vars%bd_col(c,j))
          beta_list(icat) = cece(icat) / cect_dyn(c,j)
          ! if too tiny or zero, e.g. due to tiny or zero cec_cation_vr,
          ! the following 'solve_eq' function for soil pH won't work well.
          beta_list(icat) = max(2.0e-4_r8, beta_list(icat))
          keq_list(icat) = 10**soilstate_vars%log_km_col(c,j,icat)
        end do

        if (sum(beta_list)>=(1.0_r8-2.0e-4_r8)) then
          ! if sum(beta_list)>=1.0, it will cause mathematic issue in following 'beta_h' calculation;
          !   and it will also cause 'solve_eq' issue.
          ! the limit of '1.-2e-4' would likely constrain soil pH of up to ~8.6.
          beta_list = beta_list * ((1.0_r8-2.0e-4_r8)/sum(beta_list))
        end if
        beta_h = 1._r8 - sum(beta_list(1:ncations))

        soil_ph(c,j) = solve_eq(net_charge_vr(c,j), co2_atm, beta_list, keq_list, &
                                EWParamsInst%cations_valence)

        ! calculate the implications on HCO3- & CO3 --
        bicarbonate_vr(c,j) = ph_to_hco3(soil_ph(c,j), co2_atm)
        carbonate_vr(c,j) = hco3_to_co3(bicarbonate_vr(c,j), soil_ph(c,j))

        bicarbonate_vr(c,j) = mol_to_mass(bicarbonate_vr(c,j), mass_hco3, h2osoi_vol(c,j))
        carbonate_vr(c,j) = mol_to_mass(carbonate_vr(c,j), mass_co3, h2osoi_vol(c,j))

        ! calculate the change in aqueous concentration due to equilibrium reaction
        do icat = 1,ncations
          conc(icat) = beta_list(icat) / (beta_h * keq_list(icat) / &
                       10**(-soil_ph(c,j)))**EWParamsInst%cations_valence(icat)

          !!write (100+iam, *) 'conc', j,icat, conc(icat), beta_list(icat), beta_h, soil_ph(c,j)

          cec_cation_flux_vr(c,j,icat) = & 
            ( mol_to_mass(conc(icat), EWParamsInst%cations_mass(icat), &
                          h2osoi_liqvol(c,j)) - & 
              cation_vr(c,j,icat) * h2osoi_liqvol(c,j) / h2osoi_vol(c,j) ) / dt

          !!write (100+iam, *) 'cec_cation_flux_vr', j,icat, conc(icat), h2osoi_liqvol(c,j), cation_vr(c,j,icat), h2osoi_vol(c,j)
        end do

        ! calculate the implications on cation exchange capacity
        ! 
        ! Figure 13.12: Y_org,C = -59 + 51*pH, unit: meq 100g-1 org C
        ! Stevenson, F. J. (1982), Chapter 13: Electrochemical and ion-exchange properties, 
        ! in 'Hummus Chemistry: Genesis, Composition, Reactions’, John Wiley & Sons, Inc., p. 328.
        ! The cited original study: Wisconsin soil
        ! https://acsess.onlinelibrary.wiley.com/doi/abs/10.2136/sssaj1964.03615995002800040020x
        !
        ! Also, pH-dependent CEC will not decline when pH <= 4, because charge already becomes 
        ! positive. See p59 of
        ! Kroll, E. S., Okjemstad, J. O., & Baldock, J. A. (2004). Functions of soil organic matter 
        ! and the effect on soil properties (pp. iii, 107 p. : ill.; 30 cm.). Canberra, A.C.T.: 
        ! Cooperative Research Centre for Greenhouse Accounting.
        oc_frac = 0.58_r8 * soilstate_vars%cellorg_col(c,j) / soilstate_vars%bd_col(c,j)
        if (prev_pH > 4._r8) then
          if (soil_pH(c,j) > 4._r8) then
            cect_delta(c,j) = (soil_pH(c,j) - prev_pH) * oc_frac * 51
          else
            cect_delta(c,j) = (4._r8 - prev_pH) * oc_frac * 51
          end if
        else
          if (soil_pH(c,j) > 4._r8) then
            cect_delta(c,j) = (soil_pH(c,j) - 4._r8) * oc_frac * 51
          else
            cect_delta(c,j) = 0._r8
          end if
        end if

        if (cect_delta(c,j) > 0._r8) then
          ! do not release any cation
          cece_delta(c,j,1:ncations) = 0._r8
        else
          ! equally release all the cations (negative)
          ! (but, similar to the numerical catch on beta_list above, do not release
          !  a cation if it is already too tiny; adjust the total release accordingly)
          ! (this calculation is for cece_delta & cect_delta only; no effect on
          !  cec_cation_flux_vr)
          beta_for_release(1:ncations) = cece(1:ncations) / cect_dyn(c,j)
          beta_h_for_release = 1._r8 - sum(beta_for_release)
          if (beta_h_for_release < 0.01_r8) then
            beta_h_for_release = 0._r8
          end if
          do icat = 1,ncations
            if (beta_for_release(icat) < 0.01_r8) then
              beta_for_release(icat) = 0._r8
            end if
            cece_delta(c,j,icat) = cect_delta(c,j) * beta_for_release(icat)
          end do
          cect_delta(c,j) = cect_delta(c,j) * (beta_h_for_release + &
                            sum(beta_for_release(1:ncations)))

          !!write (iam+100, *) 'beta_for_release', c, j, beta_for_release(1:ncations), beta_h_for_release
        end if

        !!write (iam+100, *) 'cec_delta', c, j, cece_delta(c,j,1:ncations), sum(cece_delta(c,j,1:ncations)), cect_delta(c,j)
        !!write (iam+100, *) 'cec_state', c, j, cece(1:ncations), sum(cece(1:ncations)), cect_dyn(c,j)

        do icat = 1,ncations
          ! calculate the change in aqueuous concentration due to change in total cation 
          ! exchange capacity
          ! this flux is guaranteed >= 0
          !  - is = 0, when cect_delta(c,j) > 0
          !  - is > 0 (flow into solution), when cect_delta(c,j) < 0
          cec_cation_flux2_vr(c,j,icat) = &
            meq_to_mass(-cece_delta(c,j,icat)/dt, EWParamsInst%cations_valence(icat), &
                        EWParamsInst%cations_mass(icat), soilstate_vars%bd_col(c,j))
        end do

      end do
    end do

    end associate

  end subroutine MineralEquilibria

end module EnhancedWeatheringMod
