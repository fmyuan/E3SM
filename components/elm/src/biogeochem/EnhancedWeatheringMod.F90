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
  use elm_varpar          , only : cation_mass, cation_valence
  use elm_varcon          , only : log_keq_co3, log_keq_hco3, log_keq_caco3
  use elm_varcon          , only : mass_caco3, mass_co3, mass_hco3, mass_co2, mass_h2o, mass_sio2, mass_h
  use elm_varpar          , only : nminerals, ncations, nminsec, nlevgrnd, mixing_depth
  use decompMod           , only : bounds_type
  use ColumnDataType      , only : col_ew, col_ms, col_mf, col_es, col_ws, col_wf
  use ColumnType          , only : col_pp
  use TopounitDataType    , only : top_as

  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: MineralReaction
  public :: MineralLeaching
  public :: readEnhancedWeatheringParams

  type, public :: EWParamsType
     real(r8) :: log_k_primary               (1:nminerals, 1:3)    ! log10 of primary mineral reaction rate constant at 298.15K (mol m-2 mineral surface area s-1), 1:nminerals x [H+, H2O, OH-]
     real(r8) :: e_primary                   (1:nminerals, 1:3)    ! primary mineral reaction activation energy constant (KJ mol-1), 1:nminerals x [H+, H2O, OH-]
     real(r8) :: n_primary                   (1:nminerals, 1:2)    ! reaction order of H+ and OH- catalyzed weathering, 1:nminerals x [H+, OH-]
     real(r8) :: log_keq_primary             (1:nminerals)         ! log10 of equilibrium constants for primary mineral dissolution 

     ! reaction stoichiometry: suppose the equation is 
     ! primary mineral + proton + (water) = cations + SiO2 + (water)
     ! coefficient before the mineral is always 1
     real(r8) :: primary_stoi_proton         (1:nminerals)    ! reaction stoichiometry coefficient in front of H+, 1:nminerals
     real(r8) :: primary_stoi_h2o_in         (1:nminerals)    ! reaction stoichiometry coefficient in front of water consumed, 1:nminerals
     real(r8) :: primary_stoi_cations        (1:nminerals, 1:ncations)    ! reaction stoichiometry coefficient in front of cations, 1:nminerals x 1:ncations
     real(r8) :: primary_stoi_silicate       (1:nminerals)    ! reaction stoichiometry coefficient in front of SiO2, 1:nminerals
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

    tString='primary_stoi_silicate'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_silicate, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_stoi_h2o_out'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_h2o_out, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

  end subroutine readEnhancedWeatheringParams


  !-----------------------------------------------------------------------
  subroutine MineralReaction(bounds, num_soilc, filter_soilc)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, calculate the dissolution of primary minerals
    ! and precipitation of secondary minerals. 
    !
    ! ! USES
    ! rgas = universal gas constant [= 8314.467 J/K/kmole]
    use elm_varcon      , only : spval, rgas, secspday, zisoi, dzsoi
    use elm_time_manager, only : get_step_size, get_curr_date
    use landunit_varcon , only : istice, istwet, istsoil, istdlak, istcrop, istice_mec
    !
    ! !ARGUMENTS:
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    !
    ! !LOCAL VARIABLES:
    integer :: c, fc, topo, lu, m, n, a
    integer :: maxlayer               ! the deepest decomposition layer
    integer :: kyr                    ! current year
    integer :: kmo                    ! month of year  (1, ..., 12)
    integer :: kda                    ! day of month   (1, ..., 31)
    integer :: mcsec                  ! seconds 
    integer :: current_date
    real(r8):: dt
    real(r8):: co2_ppmv_val
    real(r8):: log_bicarbonate, log_carbonate
    real(r8):: log_k_dissolve_acid, log_k_dissolve_neutral, log_k_dissolve_base, log_omega, k_tot
    !-----------------------------------------------------------------------
    associate( &
         !
         ! Forcing variables
         !
         forc_app               => col_ew%forc_app, & ! Input:  [real(r8) (:)] application rate (kg rock m-2 year-1)
         forc_min               => col_ew%forc_min, & ! Input:  [real(r8) (:, 1:nminerals) weight percentage of minerals in rock (1:nminerals) (kg mineral kg-1 rock)
         forc_pho               => col_ew%forc_pho, & ! Input:  [real(r8) (:)] weight percentage of phosphorus content in rock (gP kg-1 rock)
         forc_gra               => col_ew%forc_gra, & ! Input:  [real(r8) (:)] grain size (um diameter)
         forc_sph               => col_ew%forc_sph, & ! Input:  [real(r8) (:)] soil pH
         !
         ! soil pH and ionic states 
         !
         soil_ph                => col_ms%soil_ph              , & ! Output: [real(r8) (:,:)] calculated soil pH, at present equal to the forcing (1:nlevgrnd)
         bicarbonate_vr         => col_ms%bicarbonate_vr       , & ! Output: [real(r8) (:,:)] calculated HCO3- concentration in each layer of the soil (g m-3) (1:nlevgrnd)
         carbonate_vr           => col_ms%carbonate_vr         , & ! Output: [real(r8) (:,:)] calculated CO3-- concentration in each layer of the soil (g m-3) (1:nlevgrnd)
         cation_vr              => col_ms%cation_vr            , & ! Output [real(r8) (:,:,:)] cation mass in each layer of the soil (g m-3) (1:nlevgrnd, 1:ncations)
         cation                 => col_ms%cation               , & ! Output [real(r8) (:,:)] vertically integrated cation mass (g m-2) (1:ncations)
         !
         ! Primary mineral state
         !
         primary_mineral_vr     => col_ms%primary_mineral_vr   , & ! Output [real(r8) (:,:,:)] primary mineral mass in each layer of the soil (g m-3) (1:nlevgrnd, 1:nminerals)

         silicate_vr            => col_ms%silicate_vr            , & ! Output [real(r8) (:,:)] silica mass in each layer of the soil (g m-3) (1:nlevgrnd)
         armor_thickness_vr     => col_ms%armor_thickness_vr   , & ! Output [real(r8) (:,:,:)] thickness of the armoring layer on the primary mineral (um) (1:nlevgrnd, 1:nminerals)
         ssa                    => col_ms%ssa                  , & ! Output [real(r8) (:,:)] specific surface area of the primary minerals (m2 kg-1 mineral) (1:nminerals)

         primary_mineral        => col_ms%primary_mineral      , & ! Output [real(r8) (:,:)] vertically integrated primary mineral mass (g m-2) (1:nminerals)
         silicate               => col_ms%silicate             , & ! Output [real(r8) (:)] vertically integrated SiO2 mass (g m-2)

         !
         ! Primary mineral flux
         !
         primary_added_vr               => col_mf%primary_added_vr       , & ! Output [real(r8) (:,:,:)] primary mineral addition through rock powder application (g m-3 s-1) (1:nlevgrnd, 1:nminerals)
         primary_dissolve_vr            => col_mf%primary_dissolve_vr    , & ! Output [real(r8) (:,:,:)] primary mineral loss through dissolution reaction (g m-3 s-1) (1:nlevgrnd, 1:nminerals)

         primary_proton_flux_vr         => col_mf%primary_proton_flux_vr , & ! Output [real(r8) (:,:)] consumed H+ due to all the dissolution reactions (g m-3 s-1) (1:nlevgrnd)
         primary_cation_flux_vr         => col_mf%primary_cation_flux_vr , & ! Output [real(r8) (:,:,:) cations produced due to all the dissolution reactions (g m-3 s-1) (1:nlevgrnd, 1:ncations)
         primary_h2o_flux_vr            => col_mf%primary_h2o_flux_vr    , & ! Output [real(r8) (:,:)] net of water produced and consumed due to all the dissolution reaction (g m-3 s-1) (1:nlevgrnd)
         primary_silicate_flux_vr       => col_mf%primary_silicate_flux_vr, & ! Output [real(r8) (:,:)] SiO2 produced due to all the dissolution reaction (g m-3 s-1) (1:nlevgrnd)

         !primary_prelease_vr            => col_mf%mineral_prelease_vr    , & ! Output [real(r8) (:,:)] release of soluble phosphorus due to all the dissolution reaction (g m-3 s-1) (1:nlevgrnd)

         r_dissolve_vr                  => col_mf%r_dissolve_vr          , & ! Output [real(r8) (:,:)] rate at which the dissolution reaction happens (1e-7 mol m-3 s-1) (1:nlevgrnd, 1:nminerals)
         primary_co2_flux_vr            => col_mf%primary_co2_flux_vr    , &
         secondary_co2_flux_vr          => col_mf%secondary_co2_flux_vr  , &

         primary_added                  => col_mf%primary_added          , & ! Output [real(r8) (:,:)] vertially integrated primary mineral addition (g m-2 s-1) (1:nminerals)
         primary_dissolve               => col_mf%primary_dissolve       , & ! Output [real(r8) (:,:)] vertically integrated primary mineral loss (g m-2 s-1) (1:nminerals)
         primary_proton_flux            => col_mf%primary_proton_flux    , & ! Output [real(r8) (:)] vertically integrated consumed H+ due to all the dissolution reactions (1e-7 mol m-2 s-1)
         primary_cation_flux            => col_mf%primary_cation_flux    , & ! Output [real(r8) (:,:) vertically integrated cations produced due to all the dissolution reactions (g m-2 s-1) (1:ncations)
         primary_h2o_flux               => col_mf%primary_h2o_flux       , & ! Output [real(r8) (:)] vertically integrated net water produced and consumed due to all the dissolution reaction (g m-2 s-1)
         primary_silicate_flux          => col_mf%primary_silicate_flux  , & ! Output [real(r8) (:)] vertically integrated SiO2 production (g m-2 s-1)
         !primary_prelease_vr            => col_mf%primary_prelease       , & ! Output [real(r8) (:,:)] vertically integrated release of soluble phosphorus (g m-2 s-1)

         !
         ! Secondary mineral state
         ! 
         secondary_mineral_vr           => col_ms%secondary_mineral_vr      , & ! Output [real(r8) (:,:,:)] secondary mineral mass in each layer of the soil (g m-3) (1:nlevgrnd, 1:nminsec)
         secondary_mineral              => col_ms%secondary_mineral         , & ! Output [real(r8) (:,:)] vertically integrated secondary mineral mass (g m-2) (1:nminsec)

         ! 
         ! Secondary mineral flux
         ! 
         secondary_cation_flux_vr       => col_mf%secondary_cation_flux_vr  , & ! Output [real(r8) (:,:,:) cations consumed due to precipitation of secondary minerals (g m-3 s-1) (1:nlevgrnd, 1:ncations)
         r_precip_vr                    => col_mf%r_precip_vr               , & ! Output [real(r8) (:,:)] rate at which the precipitation of secondary mineral happens (1e-7 mol m-3 s-1) (1:nlevgrnd, 1:nminsec)
         secondary_cation_flux          => col_mf%secondary_cation_flux     , & ! Output [real(r8) (:,:)] vertically integrated cations consumed due to formation of secondary minerals (g m-2 s-1)

         !
         ! Other related
         !
         tsoi                           => col_es%t_soisno                 , &
         h2osoi_liqvol                  => col_ws%h2osoi_liqvol            , &
         dz                             => col_pp%dz                         & ! Input:  [real(r8) (:,:) ]  layer thickness (m)
    )

    dt      = real( get_step_size(), r8 )

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      topo = col_pp%topounit(c)

      ! ---------------------------------------------------------------
      ! TODO: hard code the forcing values here
      ! Hubbard Brook
      call get_curr_date(kyr, kmo, kda, mcsec)
      current_date = kyr*10000 + kmo*100 + kda

      if (current_date .eq. 19991019) then
          ! 55 tons / 11.8 ha = 0.466 kg / m2, applied over one day
          forc_app(c) = 0.466_r8
      else
          forc_app(c) = 0._r8
      end if
      forc_min(c, 1) = 1._r8
      forc_min(c, 2) = 0._r8
      forc_pho(c   ) = 0._r8
      forc_gra(c   ) = 9.6_r8 ! 9.6 um
      if (current_date .lt. 19991019) then
          forc_sph(c) = 3.5
      else
          forc_sph(c) = 4.5
      end if
      ! ---------------------------------------------------------------

      ! Specific surface area depends on the grain size of the mineral, following
      !    Strefler, J., Amann, T., Bauer, N., Kriegler, E., and Hartmann, J.: Potential and 
      !       costs of carbon dioxide removal by enhanced weathering of rocks, Environ. Res.
      !       Lett., 13, 034010, https://doi.org/10.1088/1748-9326/aaa9c4, 2018.
      ! TODO: more accurate method from geochemistry 
      !   at ~100 um magnitude, the grains are individual minerals
      !   weighted average of the specific surface area of each mineral (m^2 g-1)
      ssa(c) = 69.18_r8 * (forc_gra(c) ** (-1.24_r8)) ! unit: m^2 g-1

      ! find the maximum layer of reaction
      maxlayer = 1
      do n = 2,nlevgrnd
        if (zisoi(n-1) < mixing_depth) then
          maxlayer = n
        end if
      end do

      do n = 1,maxlayer

        ! TODO: will calculate soil pH dynamically instead of using prescribed values
        soil_ph(c, n) = forc_sph(c)

        ! TODO: calculate CO2 ppmv from soil heterotrohpic respiration and diffusion rates
        ! Step 1. Get the CO2 concentration
        co2_ppmv_val = top_as%pco2bot(topo) / top_as%pbot(topo) * 1.0e6_r8

        ! Step 2-3. Calculate HCO3- and CO3 2- concentrations
        log_bicarbonate = log_keq_hco3 + 4.7854_r8 + log10(co2_ppmv_val) + soil_ph(c,n)
        log_carbonate = log_keq_co3 - 0.007240 + soil_ph(c,n) + log_bicarbonate
        bicarbonate_vr(c,n) = 10**log_bicarbonate
        carbonate_vr(c,n) = 10**log_carbonate

        ! Step 4. Mineral addition and dissolution reaction
        do m = 1,nminerals
          ! evenly distributed in the top 30 centimeters
          primary_added_vr(c,n,m) = 1000._r8 * forc_app(c) * forc_min(c,m) / mixing_depth / secspday

          ! log10 of ion activity product divided by equilibrium constant
          log_omega = soil_ph(c,n) * EWParamsInst%primary_stoi_proton(m) - &
                      EWParamsInst%log_keq_primary(m)
          do a = 1,ncations
            log_omega = log_omega + EWParamsInst%primary_stoi_cations(m,a) * &
              (log10(cation_vr(c,n,a)) - log10(cation_mass(a)) - 3._r8)
          end do

          ! log10 of the reaction rate constant by individual weathering agents (log10 mol m-2 s-1)
          log_k_dissolve_acid = EWParamsInst%log_k_primary(m,1) + log10(exp(1.0)) * & 
            (-1.0e6_r8 * EWParamsInst%e_primary(m,1) / rgas * (1/tsoi(c,n) - 1/298.15_r8)) - &
            EWParamsInst%n_primary(m,1) * soil_ph(c,n) + log10(1 - 10**log_omega)

          log_k_dissolve_neutral = EWParamsInst%log_k_primary(m,2) + log10(exp(1.0)) * & 
            (-1.0e6_r8 * EWParamsInst%e_primary(m,2) / rgas * (1/tsoi(c,n) - 1/298.15_r8)) + &
            log10(1 - 10**log_omega)

          ! base reaction is sometimes not reported ...
          if (EWParamsInst%log_k_primary(m,3) >= 1e-6) then
            log_k_dissolve_base = EWParamsInst%log_k_primary(m,3) + log10(exp(1.0)) * & 
              (-1.0e6_r8 * EWParamsInst%e_primary(m,3) / rgas * (1/tsoi(c,n) - 1/298.15_r8)) - &
              EWParamsInst%n_primary(m,2) * (14 - soil_ph(c,n)) + log10(1 - 10**log_omega)
          else
            log_k_dissolve_base = 0._r8
          end if

          ! sum up the reaction rates and calculate dissolution rate in mol m-3 s-1
          ! the following is not right??? - fmy
          if (log_k_dissolve_base == 0._r8) then
            k_tot = 10**log_k_dissolve_acid + 10**log_k_dissolve_neutral
          else
            k_tot = 10**log_k_dissolve_acid + 10**log_k_dissolve_neutral + 10**log_k_dissolve_base
          end if
          r_dissolve_vr(c,n,m) = ssa(c) * primary_mineral_vr(c,n,m) * k_tot
        end do

        ! convert the reaction rate to mineral and cation fluxes
        primary_dissolve_vr(c,n,m) = r_dissolve_vr(c,n,m) * EWParamsInst%primary_mass(m)
        primary_proton_flux_vr(c,n) = 0._r8
        do m = 1,nminerals
          primary_proton_flux_vr(c,n) = primary_proton_flux_vr(c,n) + & 
            r_dissolve_vr(c,n,m) * EWParamsInst%primary_stoi_proton(m) * mass_h
        end do
        do a = 1,ncations
          primary_cation_flux_vr(c,n,a) = 0._r8
          do m = 1,nminerals
            primary_cation_flux_vr(c,n,a) = primary_cation_flux_vr(c,n,a) + &
              r_dissolve_vr(c,n,m) * EWParamsInst%primary_stoi_cations(m,a) * cation_mass(a)
          end do
        end do
        primary_h2o_flux_vr(c,n) = 0._r8
        do m = 1,nminerals
          primary_h2o_flux_vr(c,n) = primary_h2o_flux_vr(c,n) + & 
            r_dissolve_vr(c,n,m) * (EWParamsInst%primary_stoi_h2o_out(m) - EWParamsInst%primary_stoi_h2o_in(m)) * mass_h2o
        end do
        primary_silicate_flux_vr(c,n) = 0._r8
        do m = 1,nminerals
          primary_silicate_flux_vr(c,n) = primary_silicate_flux_vr(c,n) + &
            r_dissolve_vr(c,n,m) * EWParamsInst%primary_stoi_silicate(m) * mass_sio2
        end do

        !primary_prelease_vr(c,n) = primary_dissolve_vr(c,n,m) * forc_pho(c)


        ! Step 5. Secondary mineral precipitation
        ! 5.1. Calcite

      end do

    end do

    end associate

  end subroutine MineralReaction


  !-----------------------------------------------------------------------
  subroutine MineralLeaching(bounds, num_soilc, filter_soilc, dt)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the cation and bicarbonate
    ! leaching rate as a function of concentration and total 
    ! soil water outflow.
    !
    ! !USES:
    !$acc routine seq
    use elm_varpar       , only : nlevdecomp, nlevsoi, nlevgrnd
    use elm_varcon       , only : secspday, zisoi
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    real(r8)                 , intent(in)    :: dt              ! radiation time step (seconds)

    !
    ! !LOCAL VARIABLES:
    integer  :: j,c,fc,m,a                             ! indices
    real(r8) :: conc_cation_temp(1:ncations)           ! temporary variable cation concentration in kg/kg soil water to follow the nitrogen routine
    real(r8) :: conc_bicarbonate_temp                  ! temporary variable cation concentration in kg/kg soil water to follow the nitrogen routine
    real(r8) :: frac_thickness                         ! deal with the fractional layer between last layer and max allowed depth
    real(r8) :: tot_water(bounds%begc:bounds%endc)     ! total column liquid water (kg water/m2)
    real(r8) :: surface_water(bounds%begc:bounds%endc) ! liquid water to shallow surface depth (kg water/m2)
    real(r8) :: drain_tot(bounds%begc:bounds%endc)     ! total drainage flux (mm H2O /s)
    real(r8), parameter :: depth_runoff_Mloss = 0.05   ! (m) depth over which runoff mixes with soil water for ions loss to runoff; same as nitrogen runoff depth
    !-----------------------------------------------------------------------

    associate( &
         h2osoi_liq             => col_ws%h2osoi_liq                      , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
         qflx_drain             => col_wf%qflx_drain                      , & ! Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)                    
         qflx_surf              => col_wf%qflx_surf                       , & ! Input:  [real(r8) (:)   ]  surface runoff (mm H2O /s)
         dz                     => col_pp%dz                              , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)
         !
         cation_vr              => col_ms%cation_vr                       , & ! Input:  [real(r8) (:,:,:) ] cation in each soil layer (kg m-2)
         bicarbonate_vr         => col_ms%bicarbonate_vr                  , & ! Input:  [real(r8) (:,:,:) ] bicarbonate in each soil layer (kg m-2)
         ! TODO: change the unit to m^-3
         cation_leached_vr      => col_mf%cation_leached_vr               , & ! Output: [real(r8) (:,:,:) ]  rate of cation leaching (kg m-2 s-1) ! compare to N (g m-3 s-1)
         bicarbonate_leached_vr => col_mf%bicarbonate_leached_vr          , & ! Output: [real(r8) (:,:) ]  rate of cation leaching (kg m-2 s-1) ! compare to N (g m-3 s-1)
         cation_runoff_vr       => col_mf%cation_runoff_vr                , & ! Output: [real(r8) (:,:,:) ]  rate of cation loss with runoff (kg m-2 s-1) ! compare to N (g m-3 s-1)
         bicarbonate_runoff_vr  => col_mf%bicarbonate_runoff_vr           , & ! Output: [real(r8) (:,:) ]  rate of bicarbonate loss with runoff (kg m-2 s-1) ! compare to N (g m-3 s-1)
         ! 
         primary_co2_flux_vr    => col_mf%primary_co2_flux_vr             , &
         secondary_co2_flux_vr  => col_mf%secondary_co2_flux_vr           , &
         r_sequestration        => col_mf%r_sequestration                   &
    )

    ! for runoff calculation; calculate total water to a given depth
    surface_water(bounds%begc:bounds%endc) = 0._r8
    do fc = 1,num_soilc
        c = filter_soilc(fc)
        do j = 1,nlevgrnd
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
        drain_tot(c) = qflx_drain(c)
    end do

    do j = 1,nlevgrnd
        ! Loop through columns
        do fc = 1,num_soilc
          c = filter_soilc(fc)

          ! calculate the dissolved cation and bicarbonate concentration (kg/kg water)
          do a = 1,ncations
            conc_cation_temp(a) = 0._r8
            if (h2osoi_liq(c,j) > 0._r8) then
              conc_cation_temp(a) = cation_vr(c,j,a) / h2osoi_liq(c,j)
            end if
          end do
          if (h2osoi_liq(c,j) > 0._r8) then
            conc_bicarbonate_temp = bicarbonate_vr(c,j) / h2osoi_liq(c,j)
          else
            conc_bicarbonate_temp = 0._r8
          end if

          ! calculate the leaching flux as a function of the dissolved
          ! concentration and the sub-surface drainage flux
          do a = 1,ncations
            cation_leached_vr(c,j,a) = conc_cation_temp(a) * drain_tot(c) * h2osoi_liq(c,j) / tot_water(c)
          end do
          bicarbonate_leached_vr(c,j) = conc_bicarbonate_temp * drain_tot(c) * h2osoi_liq(c,j) / tot_water(c)

          !
          ! ensure that leaching rate isn't larger than soil pool size
          do a = 1,ncations
            cation_leached_vr(c,j,a) = min(cation_vr(c,j,a) / dt, cation_leached_vr(c,j,a))
          end do
          bicarbonate_leached_vr(c,j) = min(bicarbonate_vr(c,j) / dt, bicarbonate_leached_vr(c,j))

          !
          ! limit the leaching flux to a positive value
          do a = 1,ncations
            cation_leached_vr(c,j,a) = max(cation_leached_vr(c,j,a), 0._r8)
          end do
          bicarbonate_leached_vr(c,j) = max(bicarbonate_leached_vr(c,j), 0._r8)

          !
          !
          ! calculate the loss from surface runoff, assuming a shallow mixing of surface waters into soil and removal based on runoff
          if ( zisoi(j) <= depth_runoff_Mloss )  then
            do a = 1,ncations
              cation_runoff_vr(c,j,a) = conc_cation_temp(a) * qflx_surf(c) * h2osoi_liq(c,j) / surface_water(c)
            end do
            bicarbonate_runoff_vr(c,j) = conc_bicarbonate_temp * qflx_surf(c) * h2osoi_liq(c,j) / surface_water(c)

          else if ( zisoi(j-1) < depth_runoff_Mloss )  then
            do a = 1,ncations
              cation_runoff_vr(c,j,a) = conc_cation_temp(a) * qflx_surf(c) * h2osoi_liq(c,j) / surface_water(c) * frac_thickness
            end do
            bicarbonate_runoff_vr(c,j) = conc_bicarbonate_temp * qflx_surf(c) * h2osoi_liq(c,j) / surface_water(c) * frac_thickness

          else
            do a = 1,ncations
              cation_runoff_vr(c,j,a) = 0._r8
            end do
            bicarbonate_runoff_vr(c,j) = 0._r8
          end if

          !
          ! ensure that runoff rate isn't larger than soil cations pool
          do a = 1,ncations
            cation_runoff_vr(c,j,a) = min(cation_vr(c,j,a) / dt - cation_leached_vr(c,j,a), cation_runoff_vr(c,j,a))
          end do
          bicarbonate_runoff_vr(c,j) = min(bicarbonate_vr(c,j) / dt - bicarbonate_leached_vr(c,j), bicarbonate_runoff_vr(c,j))

          !
          ! limit the flux to a positive value
          do a = 1,ncations
            cation_runoff_vr(c,j,a) = max(cation_runoff_vr(c,j,a), 0._r8)
          end do
          bicarbonate_runoff_vr(c,j) = max(bicarbonate_runoff_vr(c,j), 0._r8)

       end do

    end do

    ! calculate the total CO2 sequestration rate
    do fc = 1,num_soilc
        c = filter_soilc(fc)
        r_sequestration(c) = 0._r8
        do j = 1,nlevgrnd
          r_sequestration(c) = r_sequestration(c) + primary_co2_flux_vr(c,j) - secondary_co2_flux_vr(c,j)
        end do
        ! convert from kg CO2 to gC
        r_sequestration(c) = r_sequestration(c) * 1000._r8 * 12._r8 / 44.0_r8

      ! TODO: use cation valence to further reduce by ocean efficiency
    end do

    end associate
  end subroutine MineralLeaching

end module EnhancedWeatheringMod
