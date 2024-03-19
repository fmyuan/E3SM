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
  use elm_varctl          , only : iulog, spinup_state, nyears_before_ew
  use elm_varcon          , only : log_keq_co3, log_keq_hco3, log_keq_minsec, log_keq_sio2am, alpha_minsec
  use elm_varcon          , only : mass_co3, mass_hco3, mass_co2, mass_h2o, mass_sio2, mass_h, mass_minsec
  use elm_varpar          , only : cation_mass, cation_valence, mixing_layer, mixing_depth, cation_names
  use elm_varpar          , only : nminerals, ncations, nminsec, nks
  use decompMod           , only : bounds_type
  use ColumnDataType      , only : col_ew, col_ms, col_mf, col_es, col_ws, col_wf
  use ColumnType          , only : col_pp
  use TopounitDataType    , only : top_as
  use SoilStateType       , only : soilstate_type
  use ewutils             , only : logmol_to_mass, mol_to_mass, meq_to_mass, mass_to_mol, mass_to_meq, mass_to_logmol, objective_solveq, solve_eq, ph_to_hco3, hco3_to_co3

  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: MineralInit
  public :: MineralDynamics
  public :: MineralEquilibria
  public :: MineralLeaching
  public :: readEnhancedWeatheringParams

  type, public :: EWParamsType
     real(r8) :: log_k_primary               (1:nminerals, 1:nks)  ! log10 of primary mineral reaction rate constant at 298.15K (mol m-2 mineral surface area s-1), 1:nminerals x [H+, H2O, OH-]
     real(r8) :: e_primary                   (1:nminerals, 1:nks)  ! primary mineral reaction activation energy constant (KJ mol-1), 1:nminerals x [H+, H2O, OH-]
     real(r8) :: n_primary                   (1:nminerals, 1:nks)  ! reaction order of H+ and OH- catalyzed weathering, 1:nminerals x [H+, OH-]
     real(r8) :: log_keq_primary             (1:nminerals)         ! log10 of equilibrium constants for primary mineral dissolution 

     ! reaction stoichiometry: suppose the equation is 
     ! primary mineral + proton + (water) = cations + SiO2 + (water)
     ! coefficient before the mineral is always 1
     real(r8) :: primary_stoi_proton         (1:nminerals)    ! reaction stoichiometry coefficient in front of H+, 1:nminerals
     real(r8) :: primary_stoi_h2o_in         (1:nminerals)    ! reaction stoichiometry coefficient in front of water consumed, 1:nminerals
     real(r8) :: primary_stoi_cations        (1:nminerals, 1:ncations)    ! reaction stoichiometry coefficient in front of cations, 1:nminerals x 1:ncations
     real(r8) :: primary_stoi_silica         (1:nminerals)    ! reaction stoichiometry coefficient in front of SiO2, 1:nminerals
     real(r8) :: primary_stoi_h2o_out        (1:nminerals)    ! reaction stoichiometry coefficient in front of water produced, 1:nminerals

     real(r8) :: primary_mass                (1:nminerals)    ! molar mass of the primary mineral, g/mol, 1:nminerals (e.g. Mg2SiO4 = 140.6931 g/mol)
  end type EWParamsType

  type(EWParamsType), public ::  EWParamsInst
  !$acc declare create(EWParamsInst)

contains

  !-----------------------------------------------------------------------
  subroutine readEnhancedWeatheringParams ( ncid )
    !
    ! !DESCRIPTION:
    ! Read in parameters
    !
    ! !USES:
    use ncdio_pio   , only : file_desc_t,ncd_io
    use abortutils  , only : endrun
    use shr_log_mod , only : errMsg => shr_log_errMsg
    use elm_varpar  , only : nminerals, ncations, nminsec
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'EWParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

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

    tString='primary_stoi_proton'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_proton, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_stoi_h2o_in'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_h2o_in, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_stoi_cations'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_cations, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_stoi_silica'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_silica, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_stoi_h2o_out'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_h2o_out, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_mass'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_mass, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

  end subroutine readEnhancedWeatheringParams


  !-----------------------------------------------------------------------
  subroutine MineralInit(bounds, num_soilc, filter_soilc, soilstate_vars)
    !
    ! !DESCRIPTION: 
    ! Calculate initial cation concentration from background CEC and soil pH
    ! after soil hydrology is already initialized
    ! 
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: fc,c,j,t,a             ! indices
    real(r8) :: co2_atm                ! CO2 partial pressure in atm

    associate( &
         soil_ph                        => col_ms%soil_ph                 , & ! Input:  [real(r8) (:,:)] calculated soil pH (1:nlevgrnd)
         h2osoi_vol                     => col_ws%h2osoi_vol              , & ! Input:  [real(r8) (:)] volumetric soil water content, ice + water (m3 m-3)

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

      co2_atm = top_as%pco2bot(t) / 101325

      do j = 1,mixing_layer
        ! soil pH has been initialized at cold-start
        proton_vr(c,j) = logmol_to_mass(-soil_ph(c,j), mass_h, h2osoi_vol(c,j))

        !write (iulog, *) 'initial proton_vr', c, j, proton_vr(c,j)

        ! for the sake of initial charge balance, calculate from soil pH
        bicarbonate_vr(c,j) = mol_to_mass( ph_to_hco3(soil_ph(c,j), co2_atm), mass_hco3, h2osoi_vol(c,j) )
        carbonate_vr(c,j) = mol_to_mass( hco3_to_co3(ph_to_hco3(soil_ph(c,j), co2_atm), soil_ph(c,j)), mass_co3, h2osoi_vol(c,j) )

        cec_proton_vr(c,j) = meq_to_mass(soilstate_vars%ceca_col(c,j), 1._r8, mass_h, soilstate_vars%bd_col(c,j))
        do a = 1,ncations
          cec_cation_vr(c,j,a) = meq_to_mass(soilstate_vars%cece_col(c,j,a), cation_valence(a), cation_mass(a), soilstate_vars%bd_col(c,j))
        end do

        ! calculate soil solution concentration using the equilibrium with CEC
        do a = 1,ncations
          cation_vr(c,j,a) = mol_to_mass( &
              (10**(-soil_ph(c,j)-soilstate_vars%log_km_col(c,j,a)) / & 
               (soilstate_vars%ceca_col(c,j) / soilstate_vars%cect_col(c,j)) & 
              )**cation_valence(a) * &
              soilstate_vars%cece_col(c,j,a) / soilstate_vars%cect_col(c,j) &
          , cation_mass(a), h2osoi_vol(c,j) )
        end do

        ! write (iulog, *) 'initial cation_vr', c, j, mass_to_mol(cation_vr(c,j,1), cation_mass(1), h2osoi_vol(c,j)), mass_to_mol(cation_vr(c,j,2), cation_mass(2), h2osoi_vol(c,j)), mass_to_mol(cation_vr(c,j,3), cation_mass(3), h2osoi_vol(c,j)), mass_to_mol(cation_vr(c,j,4), cation_mass(4), h2osoi_vol(c,j)), mass_to_mol(cation_vr(c,j,5), cation_mass(5), h2osoi_vol(c,j))

        ! calculate the net charge balance at the first time step
        ! mol/kg
        net_charge_vr(c,j) = 10**(-soil_ph(c,j)) - 10**(-14_r8+soil_ph(c,j)) - &
            mass_to_mol(bicarbonate_vr(c,j), mass_hco3, h2osoi_vol(c,j)) - & 
            2_r8 * mass_to_mol(carbonate_vr(c,j), mass_co3, h2osoi_vol(c,j))
        do a = 1,ncations
          net_charge_vr(c,j) = net_charge_vr(c,j) + cation_valence(a) * mass_to_mol(cation_vr(c,j,a), cation_mass(a), h2osoi_vol(c,j))
        end do

        ! write (iulog, *) 'initial net charge', c, j, net_charge_vr(c,j), 10**(-soil_ph(c,j)), - 10**(-14_r8+soil_ph(c,j)), - mass_to_mol(bicarbonate_vr(c,j), mass_hco3, h2osoi_vol(c,j)), - 2_r8 * mass_to_mol(carbonate_vr(c,j), mass_co3, h2osoi_vol(c,j)), cation_valence(1) * mass_to_mol(cation_vr(c,j,1), cation_mass(1), h2osoi_vol(c,j)), cation_valence(2) * mass_to_mol(cation_vr(c,j,2), cation_mass(2), h2osoi_vol(c,j)), cation_valence(3) * mass_to_mol(cation_vr(c,j,3), cation_mass(3), h2osoi_vol(c,j)), cation_valence(4) * mass_to_mol(cation_vr(c,j,4), cation_mass(4), h2osoi_vol(c,j)), cation_valence(5) * mass_to_mol(cation_vr(c,j,5), cation_mass(5), h2osoi_vol(c,j))

      end do
    end do

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
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilstate_type)     , intent(in)    :: soilstate_vars

    !
    ! !LOCAL VARIABLES:
    integer  :: fc,c,j,a,m             ! indices
    integer  :: kyr                    ! current year
    integer  :: kmo                    ! month of year  (1, ..., 12)
    integer  :: kda                    ! day of month   (1, ..., 31)
    integer  :: mcsec                  ! seconds 
    integer  :: current_date
    real(r8) :: equilibria_conc(0:mixing_layer)
    real(r8) :: beta_h, beta_cation, keq, h
    real(r8) :: dt
    real(r8) :: log_k_dissolve_acid, log_k_dissolve_neutral, log_k_dissolve_base
    real(r8) :: saturation_ratio, log_silica, log_carbonate
    real(r8) :: k_tot
    character(6):: n_str, m_str, c_str
    real(r8) :: frac_kaolinite

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
         secondary_mineral_vr           => col_ms%secondary_mineral_vr      , & ! Output [real(r8) (:,:,:)] secondary mineral mass in each layer of the soil (g m-3 soil [not water]) (1:nlevgrnd, 1:nminsec)

         ! 
         ! Secondary mineral flux
         ! 
         secondary_mineral_flux_vr      => col_mf%secondary_mineral_flux_vr , & ! Output [real(r8) (:,:,:) secondary mineral precipitated (g m-3 s-1) (1:nlevgrnd, 1:nminsec)
         secondary_cation_flux_vr       => col_mf%secondary_cation_flux_vr  , & ! Output [real(r8) (:,:,:) cations consumed due to precipitation of secondary minerals (g m-3 s-1) (1:nlevgrnd, 1:ncations)
         secondary_silica_flux_vr       => col_mf%secondary_silica_flux_vr  , & ! Output [real(r8) (:,:) sio2 consumed due to precipitation of secondary minerals (g m-3 s-1) (1:nlevgrnd)
         r_precip_vr                    => col_mf%r_precip_vr               , & ! Output [real(r8) (:,:)] rate at which the precipitation of secondary mineral happens (mol m-3 s-1) (1:nlevgrnd, 1:nminsec)

         !
         ! Other related
         !
         tsoi                          => col_es%t_soisno     , &
         qin                           => col_wf%qin          , & ! Input: [real(r8) (:,:) ] flux of water into soil layer [mm h2o/s]
         qout                          => col_wf%qout         , & ! Input: [real(r8) (:,:) ] flux of water out of soil layer [mm h2o/s]
         qflx_rootsoi_col              => col_wf%qflx_rootsoi , & ! Input: [real(r8) (:,:) ]  vegetation/soil water exchange (mm H2O/s) (+ = to atm)
         h2osoi_vol                    => col_ws%h2osoi_vol   , & ! Input:  [real(r8) (:)] volumetric soil water content, ice + water (m3 m-3)
         h2osoi_liqvol                 => col_ws%h2osoi_liqvol, & ! Input:  [real(r8) (:)] volumetric soil water content, liquid only (m3 m-3)
         dz                            => col_pp%dz             & ! Input:  [real(r8) (:,:) ]  layer thickness (m)
    )

    dt      = real( get_step_size(), r8 )

    do fc = 1,num_soilc
      c = filter_soilc(fc)

      !------------------------------------------------------------------------------
      ! Propagate inputs
      !------------------------------------------------------------------------------
      call get_curr_date(kyr, kmo, kda, mcsec)
      current_date = kyr*10000 + kmo*100 + kda

      site_id = 2

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
        forc_min(c, 1) = 1._r8
        forc_min(c, 2) = 0._r8
        forc_min(c, 3) = 0._r8
        forc_min(c, 4) = 0._r8
        forc_min(c, 5) = 0._r8
        forc_pho(c   ) = 0._r8
        forc_gra(c, 1:5) = 9.6_r8 ! 9.6 um

        ! weight fraction of kaolinite in soil, g g-1 soil
        frac_kaolinite = 0.004_r8
      else
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
        forc_min(c, 1) = 0._r8
        forc_min(c, 2) = 0._r8
        forc_min(c, 3) = 0.334_r8 ! albite
        forc_min(c, 4) = 0.334_r8 ! anorthite
        forc_min(c, 5) = 0.143_r8 ! epidote
        forc_pho(c   ) = 0._r8
        forc_gra(c, 1:5) = 105._r8 ! 102-107 um

        ! weight fraction of kaolinite in soil, g g-1 soil
        frac_kaolinite = 0.2_r8
      end if

      !------------------------------------------------------------------------------
      ! At long-term equilibrium, water influx at top should be equal to 
      ! (outflux at bottom + plant uptake + surface & subsurface runoff)
      ! However, plant uptake is zero in the cation balance. 
      ! Therefore, set background weathering rate equal to 
      ! (Initial cation concentration - Above layer's concentration) * 
      ! (Influx from above - Plant uptake)
      !------------------------------------------------------------------------------
      do a = 1, ncations
        ! rainwater, mg/L => mol/kg
        equilibria_conc(0) = rain_chem(c,a) * 1e-3 / cation_mass(a)
        do j = 1,mixing_layer
          h = 10**(-soilstate_vars%sph(c,j))
          beta_h = soilstate_vars%ceca_col(c,j) / soilstate_vars%cect_col(c,j)
          beta_cation = soilstate_vars%cece_col(c,j,a) / soilstate_vars%cect_col(c,j)
          keq = 10**soilstate_vars%log_km_col(c,j,a)
          equilibria_conc(j) = beta_cation/(beta_h*keq/h)**cation_valence(a)
        end do

        ! mol kg-1 * g mol-1 * mm s-1 = g m-2 s-1, divide by dz to get g m-3 s-1
        do j = 1,mixing_layer
          background_weathering_vr(c,j,a) = mol_to_mass(equilibria_conc(j) - equilibria_conc(j-1), cation_mass(a), 1e-3_r8 * (qin(c,j) - qflx_rootsoi_col(c,j)) / dz(c,j))

          ! write (iulog, *) c, j, a, 'background_weathering', equilibria_conc(j), equilibria_conc(j) - equilibria_conc(j-1), qin(c,j), qin(c,j) - qflx_rootsoi_col(c,j)
        end do
      end do

      ! ---------------------------------------------------------------
      ! Apply the primary minerals
      ! ---------------------------------------------------------------
      do j = 1,mixing_layer
        do m = 1,nminerals
          ! evenly distributed in the top 30 centimeters
          primary_added_vr(c,j,m) = 1000._r8 * forc_app(c) * forc_min(c,m) / mixing_depth / secspday
        end do
      end do

      ! ---------------------------------------------------------------
      ! Dissolution reaction
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

      do j = 1,mixing_layer
        if (h2osoi_liqvol(c,j) < 1e-6) then
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
              do a = 1,ncations
                log_omega_vr(c,j,m) = log_omega_vr(c,j,m) + & 
                  EWParamsInst%primary_stoi_cations(m,a) * & 
                  mass_to_logmol(cation_vr(c,j,a), cation_mass(a), h2osoi_vol(c,j))
              end do

              !write (iulog, *) 'log_omega_vr part 1', c, j, soil_ph(c,j)*EWParamsInst%primary_stoi_proton(m) - EWParamsInst%log_keq_primary(m), soil_ph(c,j), EWParamsInst%primary_stoi_proton(m), EWParamsInst%log_keq_primary(m)
              !do a = 1,ncations
              !  write (iulog, *) 'log_omega_vr cation', c, j, a, cation_names(a), EWParamsInst%primary_stoi_cations(m,a) * & 
              !  mass_to_logmol(cation_vr(c,j,a), cation_mass(a), h2osoi_vol(c,j)), EWParamsInst%primary_stoi_cations(m,a), mass_to_logmol(cation_vr(c,j,a), cation_mass(a), h2osoi_vol(c,j)), cation_vr(c,j,a)
              !end do

              ! check the reaction rate is not negative
              if (log_omega_vr(c,j,m) >= 0._r8) then
                r_dissolve_vr(c,j,m) = 0._r8

                write (m_str, '(I6)') m
                write (c_str, '(I6)') c
                write (n_str, '(I6)') j

                write (iulog,*) ' WARNING! Omega > 1 meaning dissolution reaction cannot proceed for mineral '//m_str//' column '//c_str//' layer '//n_str//':'

                !write (iulog, *) 'h2o_vol = ', h2osoi_vol(c,j)
                !do a = 1,ncations
                !  if (EWParamsInst%primary_stoi_cations(m,a) > 0._r8) then
                !    write (iulog, *) 'cation ', a, '=', cation_vr(c,j,a)
                !  end if
                !end do

                write (iulog, *) 'log_omega_vr part 1', soil_ph(c,j) * EWParamsInst%primary_stoi_proton(m) - EWParamsInst%log_keq_primary(m)
                do a = 1,ncations
                  write (iulog, *) 'log_omega_vr cation', a, cation_names(a), EWParamsInst%primary_stoi_cations(m,a) * & 
                  mass_to_logmol(cation_vr(c,j,a), cation_mass(a), h2osoi_vol(c,j))
                end do

              else
                ! log10 of the reaction rate constant by individual weathering agents (log10 mol m-2 s-1)
                k_tot = 0._r8

                ! only add this part if the rate is > -9999 (the NULL value, since base reaction is 
                ! sometimes not reported)
                if (EWParamsInst%log_k_primary(m,1) > -9999) then
                  log_k_dissolve_acid = EWParamsInst%log_k_primary(m,1) + log10(exp(1.0)) * & 
                    (-1.0e6_r8 * EWParamsInst%e_primary(m,1) / rgas * (1/tsoi(c,j) - 1/298.15_r8)) - &
                    EWParamsInst%n_primary(m,1) * soil_ph(c,j) + log10(1 - 10**log_omega_vr(c,j,m))
                  k_tot = k_tot + 10**log_k_dissolve_acid

                  write (iulog, *) c, j, m, 'log_k_dissolve_acid', log_k_dissolve_acid, EWParamsInst%log_k_primary(m,1), EWParamsInst%e_primary(m,1), rgas, tsoi(c,j), EWParamsInst%n_primary(m,1), soil_ph(c,j), log_omega_vr(c,j,m)
                end if

                if (EWParamsInst%log_k_primary(m,2) > -9999) then
                  log_k_dissolve_neutral = EWParamsInst%log_k_primary(m,2) + log10(exp(1.0)) * & 
                    (-1.0e6_r8 * EWParamsInst%e_primary(m,2) / rgas * (1/tsoi(c,j) - 1/298.15_r8)) + &
                    log10(1 - 10**log_omega_vr(c,j,m))
                  k_tot = k_tot + 10**log_k_dissolve_neutral
                  write (iulog, *) c, j, m, 'log_k_dissolve_neutral', log_k_dissolve_neutral, EWParamsInst%log_k_primary(m,1), EWParamsInst%e_primary(m,1), rgas, tsoi(c,j), log_omega_vr(c,j,m)
                end if

                if (EWParamsInst%log_k_primary(m,3) > -9999) then
                  log_k_dissolve_base = EWParamsInst%log_k_primary(m,3) + log10(exp(1.0)) * & 
                    (-1.0e6_r8 * EWParamsInst%e_primary(m,3) / rgas * (1/tsoi(c,j) - 1/298.15_r8)) - &
                    EWParamsInst%n_primary(m,3) * (14 - soil_ph(c,j)) + log10(1 - 10**log_omega_vr(c,j,m))
                  k_tot = k_tot + 10**log_k_dissolve_base
                  write (iulog, *) c, j, m, 'log_k_dissolve_base', log_k_dissolve_base, EWParamsInst%log_k_primary(m,3), EWParamsInst%e_primary(m,3), rgas, tsoi(c,j), EWParamsInst%n_primary(m,3), soil_ph(c,j), log_omega_vr(c,j,m)
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

        do a = 1,ncations
          primary_cation_flux_vr(c,j,a) = 0._r8
          do m = 1,nminerals
            primary_cation_flux_vr(c,j,a) = primary_cation_flux_vr(c,j,a) + &
              r_dissolve_vr(c,j,m) * EWParamsInst%primary_stoi_cations(m,a) * cation_mass(a)
          end do
        end do

        primary_h2o_flux_vr(c,j) = 0._r8
        do m = 1,nminerals
          primary_h2o_flux_vr(c,j) = primary_h2o_flux_vr(c,j) + & 
            r_dissolve_vr(c,j,m) * (EWParamsInst%primary_stoi_h2o_out(m) - EWParamsInst%primary_stoi_h2o_in(m)) * mass_h2o
        end do

        primary_silica_flux_vr(c,j) = 0._r8
        do m = 1,nminerals
          primary_silica_flux_vr(c,j) = primary_silica_flux_vr(c,j) + &
            r_dissolve_vr(c,j,m) * EWParamsInst%primary_stoi_silica(m) * mass_sio2
        end do

        do m = 1,nminerals
          primary_residue_flux_vr(c,j,m) = EWParamsInst%primary_mass(m) + &
            EWParamsInst%primary_stoi_proton(m)*mass_h - & 
            (EWParamsInst%primary_stoi_h2o_out(m)-EWParamsInst%primary_stoi_h2o_in(m))*mass_h2o - & 
            EWParamsInst%primary_stoi_silica(m) * mass_sio2
          do a = 1,ncations
            primary_residue_flux_vr(c,j,m) = primary_residue_flux_vr(c,j,m) - &
              EWParamsInst%primary_stoi_cations(m,a) * cation_mass(a)
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
      do j = 1,mixing_layer
        do m = 1,nminsec
          r_precip_vr(c,j,m) = 0._r8
        end do

        if (h2osoi_liqvol(c,j) > 1e-6) then
          ! Calcite formation (Ca2+ is cation #1)
          saturation_ratio = &
            mass_to_mol(carbonate_vr(c,j), mass_hco3, h2osoi_vol(c,j)) * &
            mass_to_mol(cation_vr(c,j,1), cation_mass(1), h2osoi_vol(c,j)) * &
            10**(-log_keq_minsec(1))

          ! Reaction rate is 
          ! r = \alpha * (\Omega - 1)
          ! \alpha = 9*1e-10 mol dm-3 (solution) s-1
          ! 
          ! Kirk, G. J. D., Versteegen, A., Ritz, K. & Milodowski, A. E. A simple reactive-transport model of calcite precipitation in soils and other porous media. Geochimica et Cosmochimica Acta 165, 108–122 (2015).
          r_precip_vr(c,j,1) = alpha_minsec(1) * max(saturation_ratio - 1._r8, 0._r8)

          ! limit the precipitation rate by the reactant's concentration
          r_precip_vr(c,j,1) = min( & 
            mass_to_mol(carbonate_vr(c,j), mass_co3, h2osoi_vol(c,j)) / dt, &
            mass_to_mol(cation_vr(c,j,1), mass_co3, h2osoi_vol(c,j)) / dt, &
            r_precip_vr(c,j,1) &
          )

          ! convert from per kg solution to per m3 soil
          ! reaction is in liquid part only
          r_precip_vr(c,j,1) = r_precip_vr(c,j,1) * h2osoi_liqvol(c,j) * 1e3_r8

          ! Kaolinite formation (Al3+ is cation #5)
          ! check silica concentration - if supersaturated, reduce to saturation point
          log_silica = mass_to_logmol(silica_vr(c,j), mass_sio2, h2osoi_vol(c,j))
          log_silica = min(log_silica, log_keq_sio2am)

          saturation_ratio = 2 * log_silica + 2 * mass_to_logmol(cation_vr(c,j,5), cation_mass(5), h2osoi_vol(c,j)) + 6 * soil_ph(c,j) - log_keq_minsec(2)

          ! Reaction rate is 
          ! r [mol m-3 s-1] = A_{bulk} [m2 m-3] * k * (\Omega - 1)
          ! Perez-Fodich, A., & Derry, L. A. (2020). A model for germanium-silicon equilibrium fractionation in kaolinite. Geochimica et Cosmochimica Acta, 288, 199–213. https://doi.org/10.1016/j.gca.2020.07.046
          r_precip_vr(c,j,2) = alpha_minsec(2) * (soilstate_vars%bd_col(c,j)*1e3*(1-soilstate_vars%cellorg_col(c,j)/ParamsShareInst%organic_max)*frac_kaolinite) * max(10**saturation_ratio - 1._r8, 0._r8)

          ! convert to mol kg-1 water s-1
          r_precip_vr(c,j,2) = r_precip_vr(c,j,2) / h2osoi_liqvol(c,j) * 1e-3_r8

          !write (iulog, *) 'r_precip_vr', c, j, r_precip_vr(c,j,1), r_precip_vr(c,j,2), h2osoi_liqvol(c,j), alpha_minsec(2), soilstate_vars%bd_col(c,j), (1._r8-soilstate_vars%cellorg_col(c,j)/ParamsShareInst%organic_max), frac_kaolinite, saturation_ratio

          ! limit the precipitation rate by the reactant's concentration
          r_precip_vr(c,j,2) = min( 10**log_silica / 2 / dt,  &
            mass_to_mol(cation_vr(c,j,5), cation_mass(5), h2osoi_vol(c,j)) / 2 / dt, &
            r_precip_vr(c,j,2) )

          ! convert from per kg solution to per m3 soil
          ! reaction is in liquid part only
          r_precip_vr(c,j,2) = r_precip_vr(c,j,2) * h2osoi_liqvol(c,j) * 1e3_r8

          ! update the fluxes
          secondary_cation_flux_vr(c,j,1) = r_precip_vr(c,j,1) * cation_mass(1)
          secondary_cation_flux_vr(c,j,5) = r_precip_vr(c,j,2) * cation_mass(5)
          do a = 2,4
            secondary_cation_flux_vr(c,j,a) = 0._r8
          end do
          secondary_mineral_flux_vr(c,j,1) = r_precip_vr(c,j,1) * mass_minsec(1)
          secondary_mineral_flux_vr(c,j,2) = r_precip_vr(c,j,2) * mass_minsec(2)
          secondary_silica_flux_vr(c,j) = r_precip_vr(c,j,2) * mass_sio2
        end if
      end do
    end do ! end of the soil column loop

    end associate

  end subroutine MineralDynamics


  !-----------------------------------------------------------------------
  subroutine MineralLeaching(bounds, num_soilc, filter_soilc, dt)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the boundary conditions of
    ! soil ions (H+, cations) caused by rain water infiltration and
    ! subsurface & surface runoff
    !
    ! !USES:
    !$acc routine seq
    use elm_varcon       , only : zisoi
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    real(r8)                 , intent(in)    :: dt              ! radiation time step (seconds)

    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,fc,m,a                             ! indices
    integer  :: nlevbed				                         ! number of layers to bedrock
    real(r8) :: frac_thickness                         ! deal with the fractional layer between last layer and max allowed depth
    real(r8) :: tot_water(bounds%begc:bounds%endc)     ! total column liquid water (kg water/m2)
    real(r8) :: surface_water(bounds%begc:bounds%endc) ! liquid water to shallow surface depth (kg water/m2)
    real(r8) :: drain_tot(bounds%begc:bounds%endc)     ! total drainage flux (mm H2O /s)
    real(r8), parameter :: depth_runoff_Mloss = 0.05   ! (m) depth over which runoff mixes with soil water for ions loss to runoff; same as nitrogen runoff depth
    !-----------------------------------------------------------------------

    associate( &
         nlev2bed               => col_pp%nlevbed                         , & ! Input:  [integer (:)    ]  number of layers to bedrock
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

         qflx_rootsoi_col               =>    col_wf%qflx_rootsoi         & ! Input: [real(r8) (:,:) ]  vegetation/soil water exchange (mm H2O/s) (+ = to atm)
    )

    !------------------------------------------------------------------------------
    ! Calculate the amount of proton and cation in/out flux from infiltration
    !------------------------------------------------------------------------------
    do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! mol L-1 water * g mol-1 * mm s-1 = g m-2 s-1, divide by dz to get g m-3 s-1
      j = 1
      if (qin(c,j) > 0._r8) then
        proton_infl_vr(c,j) = 10**(-rain_ph(c)) * mass_h * qin(c,j) / dz(c,j)
        do a = 1,ncations
          ! 1e-3 g L-1 water * mm s-1 = 1e-3 g m-2 s-1, divide by dz to get g m-3 s-1
          cation_infl_vr(c,j,a) = rain_chem(c,a) * 1e-3_r8 * qin(c,j) / dz(c,j)
        end do
      end if

      do j = 2,mixing_layer
        ! mol kg-1 water * g mol-1 * mm s-1, divide by dz to get g m-3 s-1
        if (qin(c,j) > 0._r8) then
          ! flow from above layer
          proton_infl_vr(c,j) = mol_to_mass(mass_to_mol(proton_vr(c,j-1), mass_h, h2osoi_vol(c,j-1)), mass_h, 1e-3_r8 * qin(c,j) / dz(c,j))
          do a = 1,ncations
            cation_infl_vr(c,j,a) = mol_to_mass(mass_to_mol(cation_vr(c,j-1,a), cation_mass(a), &
              h2osoi_vol(c,j-1)), cation_mass(a), 1e-3_r8 * qin(c,j) / dz(c,j))
          end do
        else
          ! flow into above layer
          proton_infl_vr(c,j) = mol_to_mass(mass_to_mol(proton_vr(c,j), mass_h, h2osoi_vol(c,j)), mass_h, 1e-3_r8 * qin(c,j) / dz(c,j))
          do a = 1,ncations
            cation_infl_vr(c,j,a) = mol_to_mass(mass_to_mol(cation_vr(c,j,a), cation_mass(a), &
              h2osoi_vol(c,j)), cation_mass(a), 1e-3_r8 * qin(c,j) / dz(c,j))
          end do
        end if
      end do

      ! change mixing_layer to nlevsoi and let groundwater boundary condition to zero
      do j = 1,mixing_layer
        ! mol kg-1 water * g mol-1 * mm s-1, divide by dz to get g m-3 s-1
        if (qout(c,j) > 0._r8) then
          ! flow into below layer
          proton_oufl_vr(c,j) = mol_to_mass(mass_to_mol(proton_vr(c,j), mass_h, h2osoi_vol(c,j)), mass_h, 1e-3_r8 * qout(c,j) / dz(c,j))
          do a = 1,ncations
            cation_oufl_vr(c,j,a) = mol_to_mass(mass_to_mol(cation_vr(c,j,a), cation_mass(a), &
              h2osoi_vol(c,j)), cation_mass(a), 1e-3_r8 * qout(c,j) / dz(c,j))
          end do
        else
          ! flow from below layer
          if (j < mixing_layer) then
            proton_oufl_vr(c,j) = mol_to_mass(mass_to_mol(proton_vr(c,j+1), mass_h, h2osoi_vol(c,j+1)), mass_h, 1e-3_r8 * qout(c,j) / dz(c,j))
            do a = 1,ncations
              cation_oufl_vr(c,j,a) = mol_to_mass(mass_to_mol(cation_vr(c,j+1,a), cation_mass(a), h2osoi_vol(c,j+1)), cation_mass(a), 1e-3_r8 * qout(c,j) / dz(c,j))
            end do
          else
            proton_oufl_vr(c,j) = mol_to_mass(mass_to_mol(proton_vr(c,j), mass_h, h2osoi_vol(c,j+1)), mass_h, 1e-3_r8 * qout(c,j) / dz(c,j))
            do a = 1,ncations
              cation_oufl_vr(c,j,a) = mol_to_mass(mass_to_mol(cation_vr(c,j,a), cation_mass(a), &
                h2osoi_vol(c,j+1)), cation_mass(a), 1e-3_r8 * qout(c,j) / dz(c,j))
            end do
          end if
        end if

        ! uptake by vegetation
        ! set to zero, assuming litterfall balances out uptake
        proton_uptake_vr(c,j) = 0._r8 ! mol_to_mass(mass_to_mol(proton_vr(c,j), mass_h, h2osoi_vol(c,j)), mass_h, 1e-3_r8 * qflx_rootsoi_col(c,j) / dz(c,j))
        do a = 1,ncations
          cation_uptake_vr(c,j,a) = 0._r8 ! mol_to_mass(mass_to_mol(cation_vr(c,j,a), cation_mass(a), h2osoi_vol(c,j)), cation_mass(a), 1e-3_r8 * qflx_rootsoi_col(c,j) / dz(c,j))
        end do
      end do
    end do

    !do fc = 1,num_soilc
    !  c = filter_soilc(fc)
    !  do j = 1,mixing_layer
        !write (iulog, *) c, j, 'Reaction water flux', col_ws%h2osoi_vol(c,j) * dz(c,j), qin(c,j)*dt, - qout(c,j)*dt, -qflx_rootsoi_col(c,j)*dt
        !write (iulog, *) c, j, 'Reaction proton flux', proton_vr(c,j) * dz(c,j), proton_infl_vr(c,j) * dz(c,j) * dt, -proton_oufl_vr(c,j) * dz(c,j) * dt, -proton_uptake_vr(c,j)*dt
    !  end do
    !end do

    !------------------------------------------------------------------------------
    ! Leaching (subsurface runoff) and surface runoff losses
    !------------------------------------------------------------------------------
    ! calculate the total soil water
    tot_water(bounds%begc:bounds%endc) = 0._r8
    do fc = 1,num_soilc
        c = filter_soilc(fc)
        nlevbed = nlev2bed(c)
        do j = 1,nlevbed
          tot_water(c) = tot_water(c) + h2osoi_liq(c,j)
        end do
    end do

    ! for runoff calculation; calculate total water to a given depth
    surface_water(bounds%begc:bounds%endc) = 0._r8
    do fc = 1,num_soilc
        c = filter_soilc(fc)
        do j = 1,mixing_layer
          if ( zisoi(j) <= depth_runoff_Mloss)  then
              surface_water(c) = surface_water(c) + h2osoi_liq(c,j)
          elseif ( zisoi(j-1) < depth_runoff_Mloss)  then
              frac_thickness = (depth_runoff_Mloss - zisoi(j-1)) / dz(c,j)
              surface_water(c) = surface_water(c) + h2osoi_liq(c,j) * frac_thickness
          end if
        end do
    end do

    ! Loop through columns
    do fc = 1,num_soilc
      c = filter_soilc(fc)

      do j = 1,mixing_layer
        ! calculate the leaching flux as a function of the dissolved
        ! concentration (g cation/kg water) and the sub-surface drainage flux

        if (h2osoi_liq(c,j) > 0._r8) then
          ! (drain_tot / tot_water) is the fraction water lost per second
          proton_leached_vr(c,j) = proton_vr(c,j) * qflx_drain(c) / tot_water(c)
          ! ensure the rate is not larger than the soil pool and positive
          proton_leached_vr(c,j) = max(min(proton_vr(c,j) / dt, proton_leached_vr(c,j)), 0._r8)
        else
          proton_leached_vr(c,j) = 0._r8
        end if

        do a = 1,ncations
          if (h2osoi_liq(c,j) > 0._r8) then
            ! (drain_tot / tot_water) is the fraction water lost per second
            cation_leached_vr(c,j,a) = cation_vr(c,j,a) * qflx_drain(c) / tot_water(c)
            ! ensure the rate is not larger than the soil pool and positive
            cation_leached_vr(c,j,a) = max(min(cation_vr(c,j,a) / dt, cation_leached_vr(c,j,a)), 0._r8)
          else
            cation_leached_vr(c,j,a) = 0._r8
          end if
        end do

        ! calculate the loss from surface runoff, assuming a shallow mixing of surface waters into soil and removal based on runoff

        if (h2osoi_liq(c,j) > 0._r8) then
          if ( zisoi(j) <= depth_runoff_Mloss )  then
            proton_runoff_vr(c,j) = proton_vr(c,j) * qflx_surf(c) / surface_water(c)
          else if ( zisoi(j-1) < depth_runoff_Mloss )  then
            frac_thickness = (depth_runoff_Mloss - zisoi(j-1)) / dz(c,j)
            proton_runoff_vr(c,j) = proton_vr(c,j) * qflx_surf(c) / surface_water(c) * frac_thickness
          end if
          ! ensure the rate is not larger than the soil pool and positive
          proton_runoff_vr(c,j) = max(min(proton_vr(c,j) / dt, proton_runoff_vr(c,j)), 0._r8)
        else
          proton_runoff_vr(c,j) = 0._r8
        end if

        do a = 1,ncations
          if (h2osoi_liq(c,j) > 0._r8) then
            if ( zisoi(j) <= depth_runoff_Mloss )  then
              cation_runoff_vr(c,j,a) = cation_vr(c,j,a) * qflx_surf(c) / surface_water(c)
            else if ( zisoi(j-1) < depth_runoff_Mloss )  then
              frac_thickness = (depth_runoff_Mloss - zisoi(j-1)) / dz(c,j)
              cation_runoff_vr(c,j,a) = cation_vr(c,j,a) * qflx_surf(c) / surface_water(c) * frac_thickness
            end if

            ! ensure the rate is not larger than the soil pool and positive
            cation_runoff_vr(c,j,a) = max(min(cation_vr(c,j,a) / dt, cation_runoff_vr(c,j,a)), 0._r8)
          else
            cation_runoff_vr(c,j,a) = 0._r8
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
    integer  :: fc,c,t,j,a             ! indices
    real(r8) :: co2_atm                ! CO2 partial pressure in atm
    real(r8) :: cece(1:ncations)       ! temporary container (meq 100g-1 soil)
    real(r8) :: beta_list(1:ncations)  ! temporary container for cece/cec_tot
    real(r8) :: beta_h                 ! temporary container for ceca/cec_tot
    real(r8) :: ph                     ! temporary container for soil pH at equilibrium
    real(r8) :: keq_list(1:ncations)   ! temporary container for exchange coefficients between H+ and cations
    real(r8) :: conc(1:ncations)       ! temporary container for cation concentration (mol/kg)
    real(r8) :: dt

    associate( &
        net_charge_vr                       => col_ms%net_charge_vr           , & ! Input:  [real(r8) (:,:)] net charge of the tracked ions in the soil solution system, constant over time (1:nlevgrnd) (mol kg-1)

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

      co2_atm = top_as%pco2bot(t) / 101325

      do j = 1,mixing_layer

        ! use grid search to find the pH
        do a = 1,ncations
          cece(a) = mass_to_meq(cec_cation_vr(c,j,a), cation_valence(a), cation_mass(a), soilstate_vars%bd_col(c,j))
          beta_list(a) = cece(a) / soilstate_vars%cect_col(c,j)
          keq_list(a) = 10**soilstate_vars%log_km_col(c,j,a)
        end do
        ph = solve_eq(net_charge_vr(c,j), co2_atm, beta_list, keq_list)

        ! calculate the implications on HCO3- & CO3 --
        bicarbonate_vr(c,j) = ph_to_hco3(ph, co2_atm)
        carbonate_vr(c,j) = hco3_to_co3(bicarbonate_vr(c,j), ph)

        bicarbonate_vr(c,j) = mol_to_mass(bicarbonate_vr(c,j), mass_hco3, h2osoi_vol(c,j))
        carbonate_vr(c,j) = mol_to_mass(carbonate_vr(c,j), mass_co3, h2osoi_vol(c,j))

        ! calculate the implications on H+
        ! the reaction happens in the liquid water part only
        cec_proton_flux_vr(c,j) = (mol_to_mass(10**(-ph), mass_h, h2osoi_liqvol(c,j)) - &
                                   proton_vr(c,j)*h2osoi_liqvol(c,j)/h2osoi_vol(c,j)) / dt

        ! calculate the implications on cations
        beta_h = 1._r8
        do a = 1,ncations
          beta_h = beta_h - beta_list(a)
        end do

        do a = 1,ncations
          conc(a) = beta_list(a)/(beta_h*keq_list(a)/10**(-ph))**cation_valence(a)
          ! the reaction happens in the liquid water part only
          cec_cation_flux_vr(c,j,a) = (mol_to_mass(conc(a), cation_mass(a), h2osoi_liqvol(c,j)) - &
                                       cation_vr(c,j,a)*h2osoi_liqvol(c,j)/h2osoi_vol(c,j)) / dt
        end do

        !write (iulog, *) 'cece', c, j, cece(1), cece(2), cece(3), cece(4), cece(5), soilstate_vars%cect_col(c,j)
        !write (iulog, *) 'beta_list', c, j, beta_list(1), beta_list(2), beta_list(3), beta_list(4), beta_list(5), soilstate_vars%log_km_col(c,j,1), soilstate_vars%log_km_col(c,j,2), soilstate_vars%log_km_col(c,j,3), soilstate_vars%log_km_col(c,j,4), soilstate_vars%log_km_col(c,j,5)
        !write (iulog, *) 'new_ph', c,j, ph, beta_h, co2_atm, net_charge_vr(c,j)
        !write (iulog, *) 'new cation concentration', c,j, conc(1), conc(2), conc(3), conc(4), conc(5)

        !do a = 1,ncations
        !  write (iulog, *) 'cec', c, j, a, beta_list(a), beta_h, keq_list(a), ph, cation_valence(a), conc(a)
        !end do

      end do
    end do

    end associate

  end subroutine MineralEquilibria

end module EnhancedWeatheringMod