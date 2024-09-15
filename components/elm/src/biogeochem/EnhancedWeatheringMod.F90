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
  use elm_varctl          , only : iulog
  use elm_varcon          , only : log_keq_co3, log_keq_hco3, log_keq_sio2am
  use elm_varcon          , only : mass_co3, mass_hco3, mass_co2, mass_h2o, mass_sio2, mass_h
  use elm_varcon          , only : zisoi
  use elm_varpar          , only : mixing_layer, mixing_depth, nlevgrnd, nlevsoi
  use elm_varpar          , only : nminerals, ncations, nminsecs, nks
  use decompMod           , only : bounds_type
  use ColumnDataType      , only : col_ew, col_ms, col_mf, col_es, col_ws, col_wf
  use ColumnType          , only : col_pp
  use LandunitType        , only : lun_pp
  use TopounitDataType    , only : top_as
  use SoilStateType       , only : soilstate_type
  use ewutils             , only : logmol_to_mass, mol_to_mass, meq_to_mass, mass_to_mol, mass_to_meq, advection_diffusion
  use ewutils             , only : mass_to_logmol, objective_solveq, solve_eq, ph_to_hco3, hco3_to_co3
  use domainMod           , only : ldomain ! debug print
  use shr_sys_mod         , only : shr_sys_flush
  use landunit_varcon , only : istsoil, istcrop

  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: MineralInit
  public :: MineralDynamics
  public :: MineralEquilibria
  public :: MineralLeaching
  public :: elm_erw_readnl
  public :: readEnhancedWeatheringParams

  type, public :: EWParamsType
     character(len=40), pointer  :: minerals_name      (:)      => null()
     real(r8), pointer  :: log_k_primary               (:, :)   => null()   ! log10 of primary mineral reaction rate constant at 298.15K (mol m-2 mineral surface area s-1), 1:nminerals x [H+, H2O, OH-]
     real(r8), pointer  :: e_primary                   (:, :)   => null()   ! primary mineral reaction activation energy constant (KJ mol-1), 1:nminerals x [H+, H2O, OH-]
     real(r8), pointer  :: n_primary                   (:, :)   => null()   ! reaction order of H+ and OH- catalyzed weathering, 1:nminerals x [H+, OH-]
     real(r8), pointer  :: log_keq_primary             (:)      => null()   ! log10 of equilibrium constants for primary mineral dissolution
     real(r8), pointer  :: primary_mass                (:)      => null()   ! molar mass of the primary mineral, g/mol, 1:nminerals (e.g. Mg2SiO4 = 140.6931 g/mol)

     character(len=40), pointer  :: cations_name       (:)      => null()
     real(r8), pointer  :: cations_mass                (:)      => null()   ! molar masses of the cation species, g/mol
     real(r8), pointer  :: cations_valence             (:)      => null()   ! valence of the cations
     real(r8), pointer  :: cations_diffusivity         (:)      => null()   ! diffusion coefficient of the cations in water, m2/s

     character(len=40), pointer  :: minsecs_name       (:)      => null()   ! names of the secondary minerals for the record
     real(r8), pointer  :: minsecs_mass                (:)      => null()   ! molar mass of the secondary mineral, g/mol, 1:nminsecss (e.g. Mg2SiO4 = 140.6931 g/mol)
     real(r8), pointer  :: log_keq_minsecs             (:)      => null()   ! log10 of equilibrium constants for secondary mineral precipitation
     real(r8), pointer  :: alpha_minsecs               (:)      => null()   ! alpha constants for secondary mineral precipitation

     ! reaction stoichiometry: suppose the equation is 
     ! primary mineral + proton + (water) = cations + SiO2 + (water)
     ! coefficient before the mineral is always 1

     real(r8), pointer  :: primary_stoi_proton       (:)      => null()   ! reaction stoichiometry coefficient in front of H+, 1:nminerals
     real(r8), pointer  :: primary_stoi_cations      (:, :)   => null()   ! reaction stoichiometry coefficient in front of cations, 1:nminerals x 1:ncations
     real(r8), pointer  :: primary_stoi_sio2         (:)      => null()   ! reaction stoichiometry coefficient in front of SiO2, 1:nminerals
     real(r8), pointer  :: primary_stoi_h2o          (:)      => null()   ! reaction stoichiometry coefficient in front of water (posive = produced, negative = consumed), 1:nminerals
     real(r8), pointer  :: primary_stoi_hco3         (:)      => null()   ! reaction stoichiometry coefficient in front of HCO3-, 1:nminerals

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
    use elm_varctl    , only : elm_erw_paramfile
    use spmdMod       , only : masterproc, mpicom, MPI_CHARACTER
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
    namelist / elm_erw_inparm / elm_erw_paramfile

    ! ----------------------------------------------------------------------
    ! Read namelist from standard namelist file.
    ! ----------------------------------------------------------------------

    if ( masterproc)then

       unitn = getavu()
       write(iulog,*) 'Read in elm-erw namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'elm_erw_inparm', status=ierr)
       if (ierr == 0) then
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
       write(iulog, '(A, " : ", A,/)') "   elm-erw parameter file ", trim(elm_erw_paramfile)
    end if

    ! Broadcast namelist variables read in
    call mpi_bcast (elm_erw_paramfile, len(elm_erw_paramfile) , MPI_CHARACTER, 0, mpicom, ierr)
    !
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

    allocate(character(40) :: EWParamsInst%minerals_name(1:nminerals))
    allocate(EWParamsInst%log_k_primary(1:nminerals, 1:nks))   ! 'nks' NOT read-in as above
    allocate(EWParamsInst%e_primary(1:nminerals, 1:nks))
    allocate(EWParamsInst%n_primary(1:nminerals, 1:nks))
    allocate(EWParamsInst%log_keq_primary(1:nminerals))
    allocate(EWParamsInst%primary_mass(1:nminerals))

    call ncd_inqdid(ncid,'ncations',dimid)
    call ncd_inqdlen(ncid,dimid,ncations)        ! note this will override value from 'elm_varpar' initials
    allocate(character(40) :: EWParamsInst%cations_name(1:ncations))
    allocate(EWParamsInst%cations_mass(1:ncations))
    allocate(EWParamsInst%cations_valence(1:ncations))
    allocate(EWParamsInst%cations_diffusivity(1:ncations))

    allocate(EWParamsInst%primary_stoi_proton(1:nminerals))
    allocate(EWParamsInst%primary_stoi_h2o(1:nminerals))
    allocate(EWParamsInst%primary_stoi_sio2(1:nminerals))
    allocate(EWParamsInst%primary_stoi_cations(1:nminerals,1:ncations))
    allocate(EWParamsInst%primary_stoi_hco3(1:nminerals))

    call ncd_inqdid(ncid,'nminsecs',dimid)
    call ncd_inqdlen(ncid,dimid,nminsecs)       ! note this will override value from 'elm_varpar' initials
    allocate(character(40) :: EWParamsInst%minsecs_name(1:nminsecs))
    allocate(EWParamsInst%minsecs_mass(1:nminsecs))
    allocate(EWParamsInst%log_keq_minsecs(1:nminsecs))
    allocate(EWParamsInst%alpha_minsecs(1:nminsecs))

    ! read in parameters
    tString='minerals_name'
    call ncd_io(varname=trim(tString),data=EWParamsInst%minerals_name, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_mass'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_mass, flag='read', ncid=ncid, readvar=readv)
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

    tString='log_keq_primary'
    call ncd_io(varname=trim(tString),data=EWParamsInst%log_keq_primary, flag='read', ncid=ncid, readvar=readv)
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

    ! for primary mineral's dissolution reactions (product)
    tString='primary_stoi_proton'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_proton, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_stoi_h2o'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_h2o, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_stoi_cations'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_cations, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_stoi_sio2'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_sio2, flag='read', ncid=ncid, readvar=readv)
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

    tString='alpha_minsecs'
    call ncd_io(varname=trim(tString),data=EWParamsInst%alpha_minsecs, flag='read', ncid=ncid, readvar=readv)
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
    use spmdMod, only : masterproc
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

    associate( &
         soil_ph                        => col_ms%soil_ph                 , & ! Input:  [real(r8) (:,:)] calculated soil pH (1:nlevgrnd)
         h2osoi_vol                     => col_ws%h2osoi_vol              , & ! Input:  [real(r8) (:)] volumetric soil water content, ice + water (m3 m-3)
         nlev2bed                       => col_pp%nlevbed                 , & ! Input:  [integer  (:)   ]  number of layers to bedrock
         proton_vr                      => col_ms%proton_vr               , & ! Output: calculated soil H+ concentration in soil water each soil layer (1:nlevgrnd) (g m-3 soil [not water])
         bicarbonate_vr                 => col_ms%bicarbonate_vr          , & ! Output: [real(r8) (:,:)] calculated HCO3- concentration in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)
         carbonate_vr                   => col_ms%carbonate_vr            , & ! Output: [real(r8) (:,:)] calculated CO3 2- concentration in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)
         cec_cation_vr                  => col_ms%cec_cation_vr           , & ! Output: [real(r8) (:,:,:)] adsorbed cation concentration each soil layer (1:nlevgrnd,1:ncations) (g m-3 soil [not dry soil])
         cec_proton_vr                  => col_ms%cec_proton_vr           , & ! Output: [real(r8) (:,:,:)] adsorbed H+ concentration each soil layer (1:nlevgrnd) (g m-3 soil [not dry soil])
         net_charge_vr                  => col_ms%net_charge_vr           , & ! Output:  [real(r8) (:,:)] net charge of the tracked ions in the soil solution system, constant over time (1:nlevgrnd) (mol kg-1)
         cation_vr                      => col_ms%cation_vr               & ! Output [real(r8) (:,:,:)] cation mass in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd, 1:ncations)
    )

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      t = col_pp%topounit(c)
      g = col_pp%gridcell(c)
      l = col_pp%landunit(c)
      !write (iulog, *) lun_pp%itype(l) == istsoil
      nlevbed = min(nlev2bed(c), nlevsoi)

      co2_atm = top_as%pco2bot(t) / 101325

      do j = 1,nlevbed
        ! soil pH has been initialized at cold-start
        proton_vr(c,j) = logmol_to_mass(-soil_ph(c,j), mass_h, h2osoi_vol(c,j))

        ! for the sake of initial charge balance, calculate from soil pH
        bicarbonate_vr(c,j) = mol_to_mass( ph_to_hco3(soil_ph(c,j), co2_atm), mass_hco3, h2osoi_vol(c,j) )
        carbonate_vr(c,j) = mol_to_mass( hco3_to_co3(ph_to_hco3(soil_ph(c,j), co2_atm), soil_ph(c,j)), mass_co3, h2osoi_vol(c,j) )

        cec_proton_vr(c,j) = meq_to_mass(soilstate_vars%ceca_col(c,j), 1._r8, mass_h, soilstate_vars%bd_col(c,j))
        do icat = 1, ncations
          cec_cation_vr(c,j,icat) = meq_to_mass(soilstate_vars%cece_col(c,j,icat), EWParamsInst%cations_valence(icat), EWParamsInst%cations_mass(icat), soilstate_vars%bd_col(c,j))
        end do

        ! calculate initial soil solution concentration using the equilibrium with CEC
        do icat = 1,ncations
          cation_vr(c,j,icat) = ( 10**(-soil_ph(c,j)-soilstate_vars%log_km_col(c,j,icat)) / &
               (soilstate_vars%ceca_col(c,j)/soilstate_vars%cect_col(c,j)) ) &
              **EWParamsInst%cations_valence(icat)
          cation_vr(c,j,icat) = cation_vr(c,j,icat) * &
            (soilstate_vars%cece_col(c,j,icat)/soilstate_vars%cect_col(c,j))
          cation_vr(c,j,icat) = mol_to_mass(cation_vr(c,j,icat), EWParamsInst%cations_mass(icat), h2osoi_vol(c,j))
        end do

        ! calculate the net charge balance during initializaiong
        ! mol/kg
        net_charge_vr(c,j) = 10**(-soil_ph(c,j)) - 10**(-14_r8+soil_ph(c,j)) - &
            mass_to_mol(bicarbonate_vr(c,j), mass_hco3, h2osoi_vol(c,j)) - & 
            2_r8 * mass_to_mol(carbonate_vr(c,j), mass_co3, h2osoi_vol(c,j))
        do icat = 1, ncations
          net_charge_vr(c,j) = net_charge_vr(c,j) + EWParamsInst%cations_valence(icat) * &
             mass_to_mol(cation_vr(c,j,icat), EWParamsInst%cations_mass(icat), h2osoi_vol(c,j))
        end do
      end do
    end do

    if ( masterproc)then
      write (iulog, *) '***************************************************************************'
      write (iulog, *) '*** Soil Initialization for Enhanced Weathering                       *****'
      do fc = 1,num_soilc
        c = filter_soilc(fc)
        g = col_pp%gridcell(c)

        ! Write diagnostics
        write (iulog, *) '-------------------------------------------------------------------------'
        write (iulog, *) 'grid and column: ', ldomain%latc(g), ldomain%lonc(g), g, c
        write (iulog, *) 'soil pH, proton (g m-3), cation (g m-3):'
        do j = 1,nlevbed
          write (j_str, '(I2)') j
          j_lev = 'j='//j_str
          ! write (iulog, *) 'j=', j, soil_ph(c,j), h2osoi_vol(c,j), nlevbed, nlevsoi
          write (iulog, *) j_lev, soil_ph(c,j), proton_vr(c,j), cation_vr(c,j,1:ncations)
          ! write (iulog, *) 'cation_vr',  j, icat, cation_vr(c,j,icat), soil_ph(c,j), soilstate_vars%log_km_col(c,j,icat), soilstate_vars%ceca_col(c,j), soilstate_vars%cect_col(c,j), EWParamsInst%cations_valence(icat), soilstate_vars%cece_col(c,j,icat), h2osoi_vol(c,j)    
        end do
        write (iulog, *) 'cec total, beta H+, beta cation:'
        do j = 1,nlevbed
          write (j_str, '(I2)') j
          j_lev = 'j='//j_str
          write (iulog, *) j_lev, soilstate_vars%cect_col(c,j), soilstate_vars%ceca_col(c,j) / soilstate_vars%cect_col(c,j), soilstate_vars%cece_col(c,j,1:ncations) / soilstate_vars%cect_col(c,j)
        end do
      end do
      write (iulog, *) '***************************************************************************'
    end if

    end associate
  end subroutine MineralInit

  !-----------------------------------------------------------------------
  subroutine MineralDynamics(bounds, num_soilc, filter_soilc, soilstate_vars)
    !
    ! !DESCRIPTION: 
    ! Calculate the background weathering, primary mineral dissolution, and
    ! secondary mineral precipitation fluxes. 
    ! 
    ! !USES:
    ! rgas = universal gas constant [= 8314.467 J/K/kmole]
    use elm_varcon       , only : spval, rgas, secspday
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
    integer  :: kyr                    ! current year
    integer  :: kmo                    ! month of year  (1, ..., 12)
    integer  :: kda                    ! day of month   (1, ..., 31)
    integer  :: mcsec                  ! seconds 
    integer  :: current_date
    real(r8) :: equilibria_conc(1:nlevsoi)
    real(r8) :: beta_h, beta_cation, keq, h
    real(r8) :: dt
    real(r8) :: log_k_dissolve_acid, log_k_dissolve_neutral, log_k_dissolve_base
    real(r8) :: saturation_ratio, log_silica, log_carbonate
    real(r8) :: k_tot
    real(r8) :: qflx_drain_layer, qflx_surf_layer, dzsum
    real(r8), parameter :: depth_runoff_Mloss = 0.05   ! (m) depth over which runoff mixes with soil water for ions loss to runoff; same as nitrogen runoff depth

    ! TEMPORARY - pick site
    integer :: site_id

    associate( &
         !
         ! Forcing variables
         !
         forc_app                       => col_ew%forc_app                 , & ! Input:  [real(r8) (:)] application rate (kg rock m-2 year-1)
         forc_min                       => col_ew%forc_min                 , & ! Input:  [real(r8) (:,:) weight percentage of minerals in rock (1:nminerals) (kg mineral kg-1 rock)
         forc_pho                       => col_ew%forc_pho                 , & ! Input:  [real(r8) (:)] weight percentage of phosphorus content in rock (gP kg-1 rock)
         forc_gra                       => col_ew%forc_gra                 , & ! Input:  [real(r8) (:,:)] grain size (1:nminerals) (um diameter)
         rain_ph                        => col_ew%rain_ph                  , & ! Input:  [real(r8) (:)] pH of rain water
         rain_chem                      => col_ew%rain_chem                , & ! Input:  [real(r8) (:,:)] cation concentration in rain water (excluding H+) (g m-3 rain water) (1:ncations)

         !
         ! Background weathering flux
         !
         background_weathering_vr       => col_mf%background_weathering_vr , & ! Output: [real(r8) (:)] background weathering rate (g m-3 s-1)

         !
         ! soil pH and ionic states 
         !
         soil_ph                        => col_ms%soil_ph                 , & ! Output: [real(r8) (:,:)] calculated soil pH (1:nlevgrnd)
         bicarbonate_vr                 => col_ms%bicarbonate_vr          , & ! Output: [real(r8) (:,:)] calculated HCO3- concentration in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)
         carbonate_vr                   => col_ms%carbonate_vr            , & ! Output: [real(r8) (:,:)] calculated CO3 2- concentration in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)
         cation_vr                      => col_ms%cation_vr               , & ! Output [real(r8) (:,:,:)] cation mass in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd, 1:ncations)

         !
         ! Primary mineral state
         !
         primary_mineral_vr     => col_ms%primary_mineral_vr             , & ! Output [real(r8) (:,:,:)] primary mineral mass in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd, 1:nminerals)

         silica_vr              => col_ms%silica_vr                      , & ! Output [real(r8) (:,:)] silica mass in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)
         armor_thickness_vr     => col_ms%armor_thickness_vr             , & ! Output [real(r8) (:,:,:)] thickness of the armoring layer on the primary mineral (um) (1:nlevgrnd, 1:nminerals)
         ssa                    => col_ms%ssa                            , & ! Output [real(r8) (:,:)] specific surface area of the primary minerals (m2 g-1 mineral) (1:nminerals)

         !
         ! Primary mineral flux
         !
         primary_added_vr               => col_mf%primary_added_vr       , & ! Output [real(r8) (:,:,:)] primary mineral addition through rock powder application (g m-3 s-1) (1:nlevgrnd, 1:nminerals)
         primary_dissolve_vr            => col_mf%primary_dissolve_vr    , & ! Output [real(r8) (:,:,:)] primary mineral loss through dissolution reaction (g m-3 s-1) (1:nlevgrnd, 1:nminerals)

         primary_proton_flux_vr         => col_mf%primary_proton_flux_vr , & ! Output [real(r8) (:,:)] consumed H+ due to all the dissolution reactions (g m-3 s-1) (1:nlevgrnd)
         primary_cation_flux_vr         => col_mf%primary_cation_flux_vr , & ! Output [real(r8) (:,:,:) cations produced due to all the dissolution reactions (g m-3 s-1) (1:nlevgrnd, 1:ncations)
         primary_h2o_flux_vr            => col_mf%primary_h2o_flux_vr    , & ! Output [real(r8) (:,:)] net of water produced and consumed due to all the dissolution reaction (g m-3 s-1) (1:nlevgrnd)
         primary_silica_flux_vr         => col_mf%primary_silica_flux_vr , & ! Output [real(r8) (:,:)] SiO2 produced due to all the dissolution reaction (g m-3 s-1) (1:nlevgrnd)
         primary_residue_flux_vr        => col_mf%primary_residue_flux_vr, & ! Output [real(r8) (:,:)] Non-SiO2 solides produced due to all the dissolution reaction (g m-3 s-1) (1:nlevgrnd,1:nminerals)

         primary_prelease_vr            => col_mf%primary_prelease_vr    , & ! Output [real(r8) (:,:)] release of soluble phosphorus due to all the dissolution reaction (g m-3 s-1) (1:nlevgrnd)

         r_dissolve_vr                  => col_mf%r_dissolve_vr          , & ! Output [real(r8) (:,:)] rate at which the dissolution reaction happens (mol m-3 s-1) (1:nlevgrnd, 1:nminerals)
         log_omega_vr                   => col_mf%log_omega_vr           , & ! Output [real(r8) (:,:)] omega parameter in the dissolution equation (1:nlevgrnd, 1:nminerals)

         !
         ! Secondary mineral state
         ! 
         secondary_mineral_vr           => col_ms%secondary_mineral_vr      , & ! Output [real(r8) (:,:,:)] secondary mineral mass in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd, 1:nminsecs)

         ! 
         ! Secondary mineral flux
         ! 
         secondary_mineral_flux_vr      => col_mf%secondary_mineral_flux_vr , & ! Output [real(r8) (:,:,:) secondary mineral precipitated (g m-3 s-1) (1:nlevgrnd, 1:nminsecs)
         secondary_cation_flux_vr       => col_mf%secondary_cation_flux_vr  , & ! Output [real(r8) (:,:,:) cations consumed due to precipitation of secondary minerals (g m-3 s-1) (1:nlevgrnd, 1:ncations)
         secondary_silica_flux_vr       => col_mf%secondary_silica_flux_vr  , & ! Output [real(r8) (:,:) sio2 consumed due to precipitation of secondary minerals (g m-3 s-1) (1:nlevgrnd)
         r_precip_vr                    => col_mf%r_precip_vr               , & ! Output [real(r8) (:,:)] rate at which the precipitation of secondary mineral happens (mol m-3 s-1) (1:nlevgrnd, 1:nminsecs)
         !
         ! Other related
         !
         tsoi                          => col_es%t_soisno     , & ! Input: [real(r8) (:,:) ] soil temperature [K]
         qin                           => col_wf%qin          , & ! Input: [real(r8) (:,:) ] flux of water into soil layer [mm h2o/s]
         qout                          => col_wf%qout         , & ! Input: [real(r8) (:,:) ] flux of water out of soil layer [mm h2o/s]
         qflx_drain                    => col_wf%qflx_drain   , & ! Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)                    
         qflx_surf                     => col_wf%qflx_surf    , & ! Input:  [real(r8) (:)   ]  surface runoff (mm H2O /s)
         qflx_rootsoi_col              => col_wf%qflx_rootsoi , & ! Input: [real(r8) (:,:) ]  vegetation/soil water exchange (mm H2O/s) (+ = to atm)
         h2osoi_vol                    => col_ws%h2osoi_vol   , & ! Input:  [real(r8) (:)] volumetric soil water content, ice + water (m3 m-3)
         h2osoi_liqvol                 => col_ws%h2osoi_liqvol, & ! Input:  [real(r8) (:)] volumetric soil water content, liquid only (m3 m-3)
         nlev2bed                       => col_pp%nlevbed     , & ! Input:  [integer  (:)   ]  number of layers to bedrock
         dz                            => col_pp%dz             & ! Input:  [real(r8) (:,:) ]  layer thickness (m)
    )

    dt      = real( get_step_size(), r8 )

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      g = col_pp%gridcell(c)
      !topo = col_pp%topounit(c)
      nlevbed = min(nlev2bed(c), nlevsoi)

      !------------------------------------------------------------------------------
      ! Propagate inputs
      !------------------------------------------------------------------------------
      call get_curr_date(kyr, kmo, kda, mcsec)
      current_date = kyr*10000 + kmo*100 + kda

      site_id = 0  ! 0: read-in, 1: HB site hard-wired, 2: uc-davis site hard-wired

      if (site_id == 1) then
        ! Hubbard Brook
        ! rain pH data from the monitoring station in Hubbard Brook, 
        ! National Atmospheric Deposition Program
        ! https://nadp.slh.wisc.edu/sites/ntn-NH02/
        rain_ph(c) = 4.8_r8
        ! Ca = 0.055 mg/L, Mg = 0.015 mg/L, Na = 0.075 mg/L, K = 0.014 mg/L, Al = 0
        rain_chem(c, 1) = 0.055_r8
        rain_chem(c, 2) = 0.015_r8
        rain_chem(c, 3) = 0.075_r8
        rain_chem(c, 4) = 0.014_r8
        rain_chem(c, 5) = 0._r8

        if (current_date .eq. 19991019) then
            ! 55 tons / 11.8 ha = 0.466 kg / m2, applied over one day
            forc_app(c) = 0.466_r8
        else
            forc_app(c) = 0._r8
        end if
        forc_min(c, 1:nminerals) = 0._r8
        forc_min(c, 1) = 1._r8
        forc_pho(c   ) = 0._r8
        forc_gra(c, 1:nminerals) = 9.6_r8 ! 9.6 um

      else if (site_id == 2) then
        ! U.C. Davis
        ! rain pH data from the monitoring station in Davis,  
        ! National Atmospheric Deposition Program
        ! https://nadp.slh.wisc.edu/sites/ntn-CA88/
        rain_ph(c) = 6.2_r8
        ! Ca = 0.055 mg/L, Mg = 0.015 mg/L, Na = 0.075 mg/L, K = 0.014 mg/L, Al = 0
        rain_chem(c, 1) = 0.06_r8
        rain_chem(c, 2) = 0.06_r8
        rain_chem(c, 3) = 0.025_r8
        rain_chem(c, 4) = 0.04_r8
        rain_chem(c, 5) = 0._r8

        if ((kyr .eq. 2019) .or. (kyr .eq. 2020)) then
          if ((kmo .eq. 9) .or. (kmo .eq. 10) .or. (kmo .eq. 11)) then
           ! 40 t ha-1 = 4 kg / m2, applied over 3 months, convert to per day
            forc_app(c) = 4._r8 / 90._r8
          else
            forc_app(c) = 0._r8
          end if
        else
          forc_app(c) = 0._r8
        end if

        !! overwrite in the control case
        forc_app(c) = 0._r8

        ! manually overwrite the mineral content
        forc_min(c,1:nminerals) = 0._r8
        forc_min(c,3) = 0.334_r8
        forc_min(c,5) = 0.143_r8
        forc_min(c,4) = 0.334_r8
        forc_gra(c, 1:nminerals) = 105._r8

     else
        !
        rain_ph(c) = 5.6_r8
        rain_chem(c, :) = 0.0_r8 ! in the new setup, rain_chem should no longer matter

        ! from read-in data of soil amendment application
        !print *, nstep_mod, jday_mod, secs_curr, forc_app(c)
        !do m = 1, nminerals
        !  print *, m, forc_min(c, m), forc_gra(c, m)
        !end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Need P content to affect NEE
        forc_pho(c) = 0._r8  ! (TODO) 0~namendnutr element(s) from soil amend application, by given index of phosphrous
      end if

      !------------------------------------------------------------------------------
      ! Pure background weathering
      ! - At long-term equilibrium, assume soil solution concentration is balanced
      !   with cation exchange. Also, assume no inter-soil-layer movement of cations.
      ! - There is only loss via (plant uptake + surface & subsurface runoff). 
      ! - Background weathering rate only replenishes that loss. 
      !------------------------------------------------------------------------------
      dzsum = 0._r8
      do j = 1,nlevbed
        dzsum = dzsum + dz(c,j)
      end do

      do icat = 1,ncations
        ! rainwater, mg/L => mol/kg
        ! equilibria_conc(0) = rain_chem(c,icat) * 1e-3 / EWParamsInst%cations_mass(icat)
        do j = 1,nlevbed
          h = 10**(-soilstate_vars%sph(c,j))
          beta_h = soilstate_vars%ceca_col(c,j) / soilstate_vars%cect_col(c,j)
          beta_cation = soilstate_vars%cece_col(c,j,icat) / soilstate_vars%cect_col(c,j)
          keq = 10**soilstate_vars%log_km_col(c,j,icat)
          equilibria_conc(j) = beta_cation/(beta_h*keq/h)**EWParamsInst%cations_valence(icat)
        end do

        ! mol kg-1 * g mol-1 * mm s-1 = g m-2 s-1,divide by dz to get g m-3 s-1
        do j = 1,nlevbed
          qflx_drain_layer = qflx_drain(c) / dzsum
          if (zisoi(j) < depth_runoff_Mloss) then
            qflx_surf_layer = qflx_surf(c) / depth_runoff_Mloss
          else
            qflx_surf_layer = 0._r8
          end if

          background_weathering_vr(c,j,icat) = &
              mol_to_mass(equilibria_conc(j), & !  - equilibria_conc(j-1)
                EWParamsInst%cations_mass(icat), 1e-3_r8 * (qflx_drain_layer + & 
                qflx_surf_layer + qflx_rootsoi_col(c,j)))
        end do
      end do

      !------------------------------------------------------------------------------
      ! Add the weathering of pre-existing calcite in the soil
      !------------------------------------------------------------------------------
      ! soilstate_vars%bd_col(c,j) * soilstate_vars%calcite_col(c,j)/100._r8

      ! ---------------------------------------------------------------
      ! Apply the primary minerals
      ! ---------------------------------------------------------------
      do j = 1,nlevbed
        do m = 1,nminerals
          ! evenly distributed in the top 30 centimeters
          if (j <= mixing_layer) then
            primary_added_vr(c,j,m) = 1000._r8 * forc_app(c) * forc_min(c,m) / mixing_depth / secspday
          else
            primary_added_vr(c,j,m) = 0._r8
          end if
        end do
      end do

      ! ---------------------------------------------------------------
      ! Primary mineral dissolution
      ! ---------------------------------------------------------------
      ! Specific surface area depends on the grain size of the mineral, following
      !    Strefler, J., Amann, T., Bauer, N., Kriegler, E., and Hartmann, J.: Potential and 
      !       costs of carbon dioxide removal by enhanced weathering of rocks, Environ. Res.
      !       Lett., 13, 034010, https://doi.org/10.1088/1748-9326/aaa9c4, 2018.
      ! TODO: more accurate method from geochemistry 
      !   at ~100 um magnitude, the grains are individual minerals
      !   weighted average of the specific surface area of each mineral (m^2 g-1)
      do m = 1,nminerals
        ssa(c,m) = 69.18_r8 * (forc_gra(c,m) ** (-1.24_r8)) ! unit: m^2 g-1
      end do

      do j = 1,nlevbed
        if (j > mixing_layer .or. h2osoi_liqvol(c,j) < 1e-6) then
          do m = 1,nminerals
            r_dissolve_vr(c,j,m) = 0._r8
          end do
        else
          ! Primary mineral dissolution
          do m = 1,nminerals
            if (primary_mineral_vr(c,j,m) == 0._r8) then
              r_dissolve_vr(c,j,m) = 0._r8
            else
              ! log10 of ion activity product divided by equilibrium constant
              log_omega_vr(c,j,m) = soil_ph(c,j) * EWParamsInst%primary_stoi_proton(m) - &
                EWParamsInst%log_keq_primary(m)
              do icat = 1,ncations
                log_omega_vr(c,j,m) = log_omega_vr(c,j,m) + & 
                  EWParamsInst%primary_stoi_cations(m,icat) * &
                  mass_to_logmol(cation_vr(c,j,icat), EWParamsInst%cations_mass(icat), h2osoi_vol(c,j))
              end do

              ! check the reaction rate is not negative
              if (log_omega_vr(c,j,m) >= 0._r8) then
                r_dissolve_vr(c,j,m) = 0._r8

                write (iulog,*) ' WARNING! Omega > 1 meaning dissolution reaction cannot proceed', ldomain%latc(g), ldomain%lonc(g), g, c, j, m, log_omega_vr(c,j,m)
                !write (iulog, *) 'log_omega_vr part 1', soil_ph(c,j) * EWParamsInst%primary_stoi_proton(m) - EWParamsInst%log_keq_primary(m)
                !do icat = 1,ncations
                !  write (iulog, *) 'log_omega_vr cation', icat, EWParamsInst%cations_name(icat), EWParamsInst%primary_stoi_cations(m,icat) * &
                !  mass_to_logmol(cation_vr(c,j,icat), EWParamsInst%cations_mass(icat), h2osoi_vol(c,j))
                !end do

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

                  !write (iulog, *) 'log_k_dissolve_acid', ldomain%latc(g), ldomain%lonc(g), g, c, j, m, log_k_dissolve_acid, EWParamsInst%log_k_primary(m,1), EWParamsInst%e_primary(m,1), rgas, tsoi(c,j), EWParamsInst%n_primary(m,1), soil_ph(c,j), log_omega_vr(c,j,m)
                end if

                if (EWParamsInst%log_k_primary(m,2) > -9000) then
                  log_k_dissolve_neutral = EWParamsInst%log_k_primary(m,2) + log10(exp(1.0)) * & 
                    (-1.0e6_r8 * EWParamsInst%e_primary(m,2) / rgas * (1/tsoi(c,j) - 1/298.15_r8)) + &
                    log10(1 - 10**log_omega_vr(c,j,m))
                  k_tot = k_tot + 10**log_k_dissolve_neutral
                  !write (iulog, *) 'log_k_dissolve_neutral', ldomain%latc(g), ldomain%lonc(g), g, c, j, m, log_k_dissolve_neutral, EWParamsInst%log_k_primary(m,1), EWParamsInst%e_primary(m,1), rgas, tsoi(c,j), log_omega_vr(c,j,m)
                end if

                if (EWParamsInst%log_k_primary(m,3) > -9000) then
                  log_k_dissolve_base = EWParamsInst%log_k_primary(m,3) + log10(exp(1.0)) * & 
                    (-1.0e6_r8 * EWParamsInst%e_primary(m,3) / rgas * (1/tsoi(c,j) - 1/298.15_r8)) - &
                    EWParamsInst%n_primary(m,3) * (14 - soil_ph(c,j)) + log10(1 - 10**log_omega_vr(c,j,m))
                  k_tot = k_tot + 10**log_k_dissolve_base
                  ! write (iulog, *) 'log_k_dissolve_base', ldomain%latc(g), ldomain%lonc(g), g, c, j, m, log_k_dissolve_base, EWParamsInst%log_k_primary(m,3), EWParamsInst%e_primary(m,3), rgas, tsoi(c,j), EWParamsInst%n_primary(m,3), soil_ph(c,j), log_omega_vr(c,j,m)
                end if

                ! further scale down the reaction rate by soil moisture, use liquid only
                ! may try more complex power law 
                ! Bao, C., Li, L., Shi, Y., & Duffy, C. (2017). Understanding watershed hydrogeochemistry: 1. Development of RT-Flux-PIHM. Water Resources Research, 53(3), 2328–2345. https://doi.org/10.1002/2016WR018934
                k_tot = k_tot * h2osoi_liqvol(c,j)

                ! calculate dissolution rate in mol m-3 s-1
                r_dissolve_vr(c,j,m) = ssa(c,m) * primary_mineral_vr(c,j,m) * k_tot

                !write (iulog, *) c, j, m, 'r_dissolve_vr', r_dissolve_vr(c,j,m), k_tot, ssa(c,m), primary_mineral_vr(c,j,m)
              end if
            end if
          end do
        end if

        ! Update the mineral and cation fluxes based on the reaction rates
        do m = 1,nminerals
          primary_dissolve_vr(c,j,m) = r_dissolve_vr(c,j,m) * EWParamsInst%primary_mass(m)
        end do

        ! do not consider the proton consumption by the mineral - it's supplied by constant
        ! delivery of CO2 and typically exceeds proton concentration in water
        primary_proton_flux_vr(c,j) = 0._r8
        !do m = 1,nminerals
        !  primary_proton_flux_vr(c,j) = primary_proton_flux_vr(c,j) + & 
        !    r_dissolve_vr(c,j,m) * EWParamsInst%primary_stoi_proton(m) * mass_h
        !end do

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

        primary_prelease_vr(c,j) = 0._r8
        do m = 1,nminerals
          primary_prelease_vr(c,j) = primary_prelease_vr(c,j) + &
              primary_dissolve_vr(c,j,m) * forc_pho(c)
        end do
      end do

      ! ---------------------------------------------------------------
      ! Secondary mineral precipitation
      ! ---------------------------------------------------------------
      secondary_cation_flux_vr(c,:,:) = 0._r8
      secondary_mineral_flux_vr(c,:,:)= 0._r8
      r_precip_vr(c,:,:)              = 0._r8
      secondary_silica_flux_vr(c,:)   = 0._r8

      do j = 1,nlevbed

        if (h2osoi_liqvol(c,j) > 1e-6) then

         ! Calcite precipitation (Ca2+ is cation #1)
         isec = 1
         icat = 1
         if (cation_vr(c,j,icat)>0._r8) then
          saturation_ratio = &
            mass_to_mol(carbonate_vr(c,j), mass_hco3, h2osoi_vol(c,j)) * &
            mass_to_mol(cation_vr(c,j,icat), EWParamsInst%cations_mass(icat), h2osoi_vol(c,j)) * &
            10**(EWParamsInst%log_keq_minsecs(icat))

          ! Reaction rate is 
          ! r = \alpha * (\Omega - 1)
          ! \alpha = 9*1e-10 mol dm-3 (solution) s-1
          ! 
          ! Kirk, G. J. D., Versteegen, A., Ritz, K. & Milodowski, A. E. A simple reactive-transport model of calcite precipitation in soils and other porous media. Geochimica et Cosmochimica Acta 165, 108–122 (2015).
          r_precip_vr(c,j,isec) = EWParamsInst%alpha_minsecs(isec) * max(saturation_ratio - 1._r8, 0._r8)

          ! limit the precipitation rate by the reactant's concentration
          r_precip_vr(c,j,isec) = min( &
             mass_to_mol(carbonate_vr(c,j), mass_co3, h2osoi_vol(c,j)) / dt, &
             mass_to_mol(cation_vr(c,j,icat), mass_co3, h2osoi_vol(c,j)) / dt, &
             r_precip_vr(c,j,isec) )

          ! convert from per kg solution to per m3 soil
          ! reaction is in liquid part only
          r_precip_vr(c,j,isec) = r_precip_vr(c,j,isec) * h2osoi_liqvol(c,j) * 1e3_r8

          ! update the fluxes for operative sec. minerals/cations
          secondary_cation_flux_vr(c,j,icat) = r_precip_vr(c,j,isec) * EWParamsInst%cations_mass(icat)
          secondary_mineral_flux_vr(c,j,isec) = r_precip_vr(c,j,isec) * EWParamsInst%minsecs_mass(isec)
         endif


         ! Kaolinite formation (Al3+ is cation #5)
         isec = 2
         icat = 5
         if (cation_vr(c,j,icat)>0._r8) then
          ! check silica concentration - if supersaturated, reduce to saturation point
          log_silica = mass_to_logmol(silica_vr(c,j), mass_sio2, h2osoi_vol(c,j))
          log_silica = min(log_silica, log_keq_sio2am)

          saturation_ratio = 2 * log_silica + &
                             2 * mass_to_logmol(cation_vr(c,j,icat), EWParamsInst%cations_mass(icat), h2osoi_vol(c,j)) + &
                             6 * soil_ph(c,j) + EWParamsInst%log_keq_minsecs(isec)

          ! Reaction rate is 
          ! r [mol m-3 s-1] = A_{bulk} [m2 m-3] * k * (\Omega - 1)
          ! Perez-Fodich, A., & Derry, L. A. (2020). A model for germanium-silicon equilibrium fractionation in kaolinite. Geochimica et Cosmochimica Acta, 288, 199–213. https://doi.org/10.1016/j.gca.2020.07.046
          r_precip_vr(c,j,isec) = EWParamsInst%alpha_minsecs(isec) * &
            (soilstate_vars%bd_col(c,j)*1e3*(1-soilstate_vars%cellorg_col(c,j)/ &
             ParamsShareInst%organic_max)*soilstate_vars%kaolinite_col(c,j)/100._r8) * &
            max(10**saturation_ratio - 1._r8, 0._r8)

          ! convert to mol kg-1 water s-1
          r_precip_vr(c,j,isec) = r_precip_vr(c,j,isec) / h2osoi_liqvol(c,j) * 1e-3_r8

          ! limit the precipitation rate by the reactant's concentration
          r_precip_vr(c,j,isec) = min( 10**log_silica / 2 / dt,  &
            mass_to_mol(cation_vr(c,j,icat), EWParamsInst%cations_mass(icat), h2osoi_vol(c,j)) / 2 / dt, &
            r_precip_vr(c,j,isec) )

          ! convert from per kg solution to per m3 soil
          ! reaction is in liquid part only
          r_precip_vr(c,j,isec) = r_precip_vr(c,j,isec) * h2osoi_liqvol(c,j) * 1e3_r8

          ! update the fluxes for operative sec. minerals/cations
          secondary_cation_flux_vr(c,j,icat) = r_precip_vr(c,j,isec) * EWParamsInst%cations_mass(icat)
          secondary_mineral_flux_vr(c,j,isec) = r_precip_vr(c,j,isec) * EWParamsInst%minsecs_mass(isec)

          secondary_silica_flux_vr(c,j) = secondary_silica_flux_vr(c,j) + &
               r_precip_vr(c,j,isec) * mass_sio2
         end if
        end if
      end do ! soil layer
    end do ! end of the soil column loop

    end associate

  end subroutine MineralDynamics


  !-----------------------------------------------------------------------
  subroutine MineralLeaching(bounds, num_soilc, filter_soilc, dt)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the boundary conditions of
    ! soil ions (H+, cations) caused by vertical movement and
    ! subsurface & surface runoff
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
    real(r8) :: frac_thickness                         ! deal with the fractional layer between last layer and max allowed depth
    real(r8) :: tot_water                              ! total column liquid water (kg water/m2)
    real(r8) :: surface_water                          ! liquid water to shallow surface depth (kg water/m2)
    real(r8) :: drain_tot                              ! total drainage flux (mm H2O /s)
    real(r8), parameter :: depth_runoff_Mloss = 0.05   ! (m) depth over which runoff mixes with soil water for ions loss to runoff; same as nitrogen runoff depth
    real(r8) :: rain_proton, rain_cations(1:ncations)  ! surface boundary condition (g m-3 H2O)
    real(r8) :: sourcesink_proton(1:nlevsoi), sourcesink_cations(1:nlevsoi,1:ncations) ! (g m-3 soil s-1)
    real(r8) :: adv_water(1:nlevsoi+1)                 ! m H2O / s, negative downward
    real(r8) :: diffus(1:nlevsoi)                      ! m2/s
    real(r8) :: rho(1:nlevsoi)                         ! "density" factor using soil water content
    real(r8) :: dcation_dt(1:nlevsoi, 1:ncations)      ! cation concentration rate, g m-3 s-1

    !-----------------------------------------------------------------------

    associate( &
         h2osoi_liq             => col_ws%h2osoi_liq                      , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
         qflx_drain             => col_wf%qflx_drain                      , & ! Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)                    
         qflx_surf              => col_wf%qflx_surf                       , & ! Input:  [real(r8) (:)   ]  surface runoff (mm H2O /s)

         cation_leached_vr      => col_mf%cation_leached_vr               , & ! Output: [real(r8) (:,:,:) ]  rate of cation leaching (g m-3 s-1)
         cation_runoff_vr       => col_mf%cation_runoff_vr                , & ! Output: [real(r8) (:,:,:) ]  rate of cation loss with runoff (g m-3 s-1)
         proton_leached_vr      => col_mf%proton_leached_vr               , & ! Output: [real(r8) (:,:,:) ]  rate of H+ leaching (g m-3 s-1)
         proton_runoff_vr       => col_mf%proton_runoff_vr                , & ! Output: [real(r8) (:,:,:) ]  rate of H+ loss with runoff (g m-3 s-1)

         qin                            => col_wf%qin                     , & ! Input: [real(r8) (:,:) ] flux of water into soil layer [mm h2o/s]
         qout                           => col_wf%qout                    , & ! Input: [real(r8) (:,:) ] flux of water out of soil layer [mm h2o/s]
         h2osoi_liqvol                  => col_ws%h2osoi_liqvol           , & ! Input:  [real(r8) (:)] volumetric soil water content, liquid only (m3 m-3)
         h2osoi_vol                     => col_ws%h2osoi_vol              , & ! Input:  [real(r8) (:)] volumetric soil water content, ice + water (m3 m-3)

         dz                             => col_pp%dz                      , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)
         nlev2bed                       => col_pp%nlevbed                 , & ! Input:  [integer  (:)   ]  number of layers to bedrock

         rain_ph                        => col_ew%rain_ph                 , & ! Output: [real(r8) (:)] pH of rain water
         rain_chem                      => col_ew%rain_chem               , & ! Output: [real(r8) (:,:)] cation concentration in rain water (excluding H+) (g m-3 rain water) (1:ncations)

         proton_vr                      => col_ms%proton_vr               , & ! Input: calculated soil H+ concentration in soil water each soil layer (1:nlevgrnd) (g m-3 soil [not water])
         cation_vr                      => col_ms%cation_vr               , & ! Output [real(r8) (:,:,:)] cation mass in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd, 1:ncations)

         proton_infl_vr                 => col_mf%proton_infl_vr          , & ! Output: [real(r8) (:,:)] proton flux carried from infiltration above (g m-3 soil s-1 [not water]) (1:nlevgrnd)
         cation_infl_vr                 => col_mf%cation_infl_vr          , & ! Output: [real(r8) (:,:,:)] cation flux carried from infiltration above (g m-3 soil s-1 [not water]) (1:nlevgrnd, 1:ncations)
         proton_uptake_vr               => col_mf%proton_uptake_vr        , & ! Output: [real(r8) (:,:)] proton flux uptake by plants (g m-3 soil s-1 [not water]) (1:nlevgrnd)

         proton_oufl_vr                 => col_mf%proton_oufl_vr          , & ! Output: [real(r8) (:,:)] proton flux carried away by infiltration (g m-3 soil s-1 [not water]) (1:nlevgrnd)
         cation_oufl_vr                 => col_mf%cation_oufl_vr          , & ! Output: [real(r8) (:,:,:)] cation flux carried away by infiltration (g m-3 soil s-1 [not water]) (1:nlevgrnd, 1:ncations)

         cation_uptake_vr               => col_mf%cation_uptake_vr        , & ! Output: [real(r8) (:,:,:)] cation flux uptaken by plants (g m-3 soil s-1 [not water]) (1:nlevgrnd, 1:ncations)
         cec_cation_flux_vr             => col_mf%cec_cation_flux_vr      , & ! Output: [real(r8) (:,:,:)] rate at which adsorbed cation is released into water (negative for adsorption into soil) (vertically resolved) (1:nlevgrnd, 1:ncations) (g m-3 s-1)

         primary_proton_flux_vr         => col_mf%primary_proton_flux_vr  , & ! Output [real(r8) (:,:)] consumed H+ due to all the dissolution reactions (g m-3 s-1) (1:nlevgrnd)
         primary_cation_flux_vr         => col_mf%primary_cation_flux_vr  , & ! Output [real(r8) (:,:,:) cations produced due to all the dissolution reactions (g m-3 s-1) (1:nlevgrnd, 1:ncations)
         background_weathering_vr       => col_mf%background_weathering_vr, & ! Output: [real(r8) (:)] background weathering rate (g m-3 s-1)
         secondary_cation_flux_vr       => col_mf%secondary_cation_flux_vr, & ! Output [real(r8) (:,:,:) cations consumed due to precipitation of secondary minerals (g m-3 s-1) (1:nlevgrnd, 1:ncations)

         qflx_rootsoi_col               =>    col_wf%qflx_rootsoi         & ! Input: [real(r8) (:,:) ]  vegetation/soil water exchange (mm H2O/s) (+ = to atm)
    )

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      g = col_pp%gridcell(c)
      ! l = col_pp%landunit(c)
      ! write (iulog, *) lun_pp%itype(l) == istsoil
      nlevbed = min(nlev2bed(c), nlevsoi)

      !------------------------------------------------------------------------------
      ! Collect the water flow boundary conditions (g m-3)
      !------------------------------------------------------------------------------
      ! mol/kg rain => g/m3 rain
      rain_proton = 10**(-rain_ph(c)) * mass_h * 1e3
      do icat = 1,ncations
        ! mg/L rain => g/m3 rain
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
          cation_uptake_vr(c,j,icat) = 0._r8 ! mol_to_mass(mass_to_mol(cation_vr(c,j,icat), EWParamsInst%cations_mass(icat), h2osoi_vol(c,j)), EWParamsInst%cations_mass(icat), 1e-3_r8 * qflx_rootsoi_col(c,j) / dz(c,j))
        end do

        do icat = 1,ncations
          sourcesink_cations(j,icat) = background_weathering_vr(c,j,icat) + & 
            primary_cation_flux_vr(c,j,icat) + cec_cation_flux_vr(c,j,icat) - &
            secondary_cation_flux_vr(c,j,icat) - cation_uptake_vr(c,j,icat)
        end do
      end do

      !------------------------------------------------------------------------------
      ! Calculate the vertical transport
      !------------------------------------------------------------------------------
      do icat = 1,ncations
        diffus(1:nlevsoi) = EWParamsInst%cations_diffusivity(icat)
        do j = 1,nlevbed
          ! note the flux rate is negative downward
          adv_water(j) = -1.0e-3_r8 * qin(c,j)
        end do
        adv_water(nlevbed + 1) = qin(c,j+1)

        !write (iulog, *) 'pre-adv', c, icat, cation_vr(c,1:mixing_layer, icat)
        !write (iulog, *) 'adv_water', c, adv_water(1:mixing_layer)
        !write (iulog, *) 'diffus', c, icat, diffus(1:mixing_layer)
        !write (iulog, *) 'sourcesink', c, icat, sourcesink_cations(1:mixing_layer,icat)
        !write (iulog, *) 'rain_cations', c, icat, rain_cations(icat)
        !write (iulog, *) 'rho', c, icat, rho(1:mixing_layer)

        call advection_diffusion( & 
          cation_vr(c,1:nlevsoi, icat), adv_water(1:nlevsoi+1), diffus(1:nlevsoi), &
          sourcesink_cations(1:nlevsoi,icat), rain_cations(icat), nlevbed, dt, &
          h2osoi_liqvol(c,1:nlevsoi), &
          dcation_dt(1:nlevsoi, icat) &
        )
      end do

      !------------------------------------------------------------------------------
      ! Update the cation concentrations using the vertical transport
      !------------------------------------------------------------------------------
      do j = 1, nlevbed
        proton_infl_vr(c,j) = 0._r8
        proton_oufl_vr(c,j) = 0._r8
      end do

      do j = 1, nlevbed
        do icat = 1, ncations
          cation_infl_vr(c,j,icat) = dcation_dt(j, icat) - sourcesink_cations(j,icat)
          cation_oufl_vr(c,j,icat) = 0._r8
        end do
      end do

      !------------------------------------------------------------------------------
      ! Update the cation concentrations using the vertical transport
      !------------------------------------------------------------------------------
      do icat = 1,ncations
        do j = 1,nlevbed
          cation_vr(c, j, icat) = cation_vr(c, j, icat) + dcation_dt(j, icat) * dt
        end do
        !write (iulog, *) 'post-adv', c, icat, cation_vr(c,1:mixing_layer, icat)
      end do

      !------------------------------------------------------------------------------
      ! Leaching (subsurface runoff) and surface runoff losses
      !------------------------------------------------------------------------------
      ! calculate the total soil water
      tot_water = 0._r8
      do j = 1,nlevbed
        tot_water = tot_water + h2osoi_liq(c,j)
      end do

      ! for runoff calculation; calculate total water to a given depth
      surface_water = 0._r8
      do j = 1,nlevbed
        if ( zisoi(j) <= depth_runoff_Mloss)  then
            surface_water = surface_water + h2osoi_liq(c,j)
        elseif ( zisoi(j-1) < depth_runoff_Mloss)  then
            frac_thickness = (depth_runoff_Mloss - zisoi(j-1)) / dz(c,j)
            surface_water = surface_water + h2osoi_liq(c,j) * frac_thickness
        end if
      end do

      do j = 1,nlevbed
        ! calculate the leaching flux as a function of the dissolved
        ! concentration (g cation/kg water) and the sub-surface drainage flux

        if (h2osoi_liq(c,j) > 0._r8) then
          ! (drain_tot / tot_water) is the fraction water lost per second
          proton_leached_vr(c,j) = proton_vr(c,j) * qflx_drain(c) / tot_water
          ! ensure the rate is not larger than the soil pool and positive
          proton_leached_vr(c,j) = max(min(proton_vr(c,j) / dt, proton_leached_vr(c,j)), 0._r8)
        else
          proton_leached_vr(c,j) = 0._r8
        end if

        do icat = 1,ncations
          if (h2osoi_liq(c,j) > 0._r8) then
            ! (drain_tot / tot_water) is the fraction water lost per second
            cation_leached_vr(c,j,icat) = cation_vr(c,j,icat) * qflx_drain(c) / tot_water
            ! ensure the rate is not larger than the soil pool and positive
            cation_leached_vr(c,j,icat) = max(min(cation_vr(c,j,icat) / dt, cation_leached_vr(c,j,icat)), 0._r8)
          else
            cation_leached_vr(c,j,icat) = 0._r8
          end if
        end do

        ! calculate the loss from surface runoff, assuming a shallow mixing of surface waters into soil and removal based on runoff

        if (h2osoi_liq(c,j) > 0._r8) then
          if ( zisoi(j) <= depth_runoff_Mloss )  then
            proton_runoff_vr(c,j) = proton_vr(c,j) * qflx_surf(c) / surface_water
          else if ( zisoi(j-1) < depth_runoff_Mloss )  then
            frac_thickness = (depth_runoff_Mloss - zisoi(j-1)) / dz(c,j)
            proton_runoff_vr(c,j) = proton_vr(c,j) * qflx_surf(c) / surface_water * frac_thickness
          end if
          ! ensure the rate is not larger than the soil pool and positive
          proton_runoff_vr(c,j) = max(min(proton_vr(c,j) / dt, proton_runoff_vr(c,j)), 0._r8)
        else
          proton_runoff_vr(c,j) = 0._r8
        end if

        do icat = 1,ncations
          if (h2osoi_liq(c,j) > 0._r8) then
            if ( zisoi(j) <= depth_runoff_Mloss )  then
              cation_runoff_vr(c,j,icat) = cation_vr(c,j,icat) * qflx_surf(c) / surface_water
            else if ( zisoi(j-1) < depth_runoff_Mloss )  then
              frac_thickness = (depth_runoff_Mloss - zisoi(j-1)) / dz(c,j)
              cation_runoff_vr(c,j,icat) = cation_vr(c,j,icat) * qflx_surf(c) / surface_water * frac_thickness
            end if

            ! ensure the rate is not larger than the soil pool and positive
            cation_runoff_vr(c,j,icat) = max(min(cation_vr(c,j,icat) / dt, cation_runoff_vr(c,j,icat)), 0._r8)
          else
            cation_runoff_vr(c,j,icat) = 0._r8
          end if
        end do

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
    real(r8) :: dt

    associate( &
        net_charge_vr                       => col_ms%net_charge_vr           , & ! Input:  [real(r8) (:,:)] net charge of the tracked ions in the soil solution system, constant over time (1:nlevgrnd) (mol kg-1)
        nlev2bed                            => col_pp%nlevbed                 , & ! Input:  [integer (:)    ]  number of layers to bedrock

        soil_ph                             => col_ms%soil_ph                 , & ! Input:  [real(r8) (:,:)] calculated soil pH (1:nlevgrnd)
        proton_vr                           => col_ms%proton_vr               , & ! Input: [real (r8) (:,:)] calculated soil H+ concentration in soil water each soil layer (1:nlevgrnd) (g m-3 soil [not water])
        cation_vr                           => col_ms%cation_vr               , & ! Input: [real(r8) (:,:,:)] cation concentration in soil water in each soil layer (1:nlevgrnd,1:ncations) (g m-3 soil [not water])
        cec_cation_vr                       => col_ms%cec_cation_vr           , & ! Input: [real(r8) (:,:,:)] adsorbed cation concentration each soil layer (1:nlevgrnd,1:ncations) (g m-3 soil [not dry soil])
        h2osoi_vol                          => col_ws%h2osoi_vol              , & ! Input: [real(r8) (:,:)] soil water volume, liquid + ice (m3 m-3)
        h2osoi_liqvol                       => col_ws%h2osoi_liqvol           , & ! Input: [real(r8) (:,:)] soil water volume, liquid (m3 m-3)

        bicarbonate_vr                      => col_ms%bicarbonate_vr          , & ! Output: [real(r8) (:,:)] calculated HCO3- concentration in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)
        carbonate_vr                        => col_ms%carbonate_vr            , & ! Output: [real(r8) (:,:)] calculated CO3 2- concentration in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd)

        cec_cation_flux_vr                  => col_mf%cec_cation_flux_vr      , & ! Output: [real(r8) (:,:,:)] rate at which adsorbed cation is released into water (negative for adsorption into soil) (vertically resolved) (1:nlevgrnd, 1:ncations) (g m-3 s-1)
        cec_proton_flux_vr                  => col_mf%cec_proton_flux_vr       & ! Output: [real(r8) (:,:)] rate at which adsorbed H+ is released into water (negative for adsorption into soil) (vertically resolved) (g m-3 s-1)
    )

    dt      = real( get_step_size(), r8 )

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      t = col_pp%topounit(c)
      nlevbed = min(nlev2bed(c), nlevsoi)

      co2_atm = top_as%pco2bot(t) / 101325

      do j = 1,nlevbed

        ! use grid search to find the pH
        do icat = 1,ncations
          cece(icat) = mass_to_meq(cec_cation_vr(c,j,icat), EWParamsInst%cations_valence(icat), &
            EWParamsInst%cations_mass(icat), soilstate_vars%bd_col(c,j))
          beta_list(icat) = cece(icat) / soilstate_vars%cect_col(c,j)
          keq_list(icat) = 10**soilstate_vars%log_km_col(c,j,icat)
        end do
        soil_ph(c,j) = solve_eq(net_charge_vr(c,j), co2_atm, beta_list, keq_list, EWParamsInst%cations_valence)
        ! write (iulog, *) 'solve_eq', c, j, beta_list(1), beta_list(2), beta_list(3), beta_list(4), beta_list(5)

        ! calculate the implications on HCO3- & CO3 --
        bicarbonate_vr(c,j) = ph_to_hco3(soil_ph(c,j), co2_atm)
        carbonate_vr(c,j) = hco3_to_co3(bicarbonate_vr(c,j), soil_ph(c,j))

        bicarbonate_vr(c,j) = mol_to_mass(bicarbonate_vr(c,j), mass_hco3, h2osoi_vol(c,j))
        carbonate_vr(c,j) = mol_to_mass(carbonate_vr(c,j), mass_co3, h2osoi_vol(c,j))

        ! calculate the implications on H+
        ! the reaction happens in the liquid water part only
        cec_proton_flux_vr(c,j) = (mol_to_mass(10**(-soil_ph(c,j)), mass_h, h2osoi_liqvol(c,j)) - &
                                   proton_vr(c,j)*h2osoi_liqvol(c,j)/h2osoi_vol(c,j)) / dt

        ! calculate the implications on cations
        beta_h = 1._r8
        do icat = 1,ncations
          beta_h = beta_h - beta_list(icat)
        end do

        do icat = 1,ncations
          conc(icat) = beta_list(icat)/(beta_h*keq_list(icat)/10**(-soil_ph(c,j)))**EWParamsInst%cations_valence(icat)
          ! the reaction happens in the liquid water part only
          cec_cation_flux_vr(c,j,icat) = (mol_to_mass(conc(icat), EWParamsInst%cations_mass(icat), &
            h2osoi_liqvol(c,j)) - cation_vr(c,j,icat)*h2osoi_liqvol(c,j)/h2osoi_vol(c,j)) / dt
        end do

      end do
    end do

    end associate

  end subroutine MineralEquilibria

end module EnhancedWeatheringMod