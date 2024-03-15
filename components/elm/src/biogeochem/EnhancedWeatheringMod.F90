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
  use shr_const_mod       , only : SHR_CONST_RGAS ! universal gas constant [= 8314.467 J/K/kmole]
  use elm_varpar          , only : nminerals, ncations, nminsec, nlevgrnd, mixing_depth
  use decompMod           , only : bounds_type
  use ColumnDataType      , only : col_ew, col_ms, col_mf, col_es, col_ws, col_wf
  use ColumnType          , only : col_pp

  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: MineralReaction
  public :: MineralLeaching
  public :: readEnhancedWeatheringParams

  !
  ! !PRIVATE DATA:
  real(r8), parameter :: co2_weight = 44.0_r8      ! molecular weight of CO2, 44.01 g/mol
  real(r8), parameter :: h2o_weight = 18.0_r8      ! molecular weight of H2O, 18.01 g/mol
  real(r8), parameter :: bicarbonate_weight = 61.0_r8   ! molecular weight of HCO3-, 61.02 g/mol
  real(r8), parameter :: silicate_weight = 96.0_r8 ! molecular weight of H4SiO4, 60. (SiO2) + 2*18 (2*H2O)

  type, public :: EWParamsType
     real(r8) :: k_primary                   (1:nminerals)    ! primary mineral reaction constant at 298.15K (mol m-2 mineral surface area s-1), 1:nminerals
     real(r8) :: e_primary                   (1:nminerals)    ! primary mineral reaction activation energy constant at 298.15K (J mol-1), 1:nminerals

     ! reaction stoichiometry: suppose the equation is 
     ! mineral + co2 + water -> cation(s) + bicarbonate + silicate (H4SiO4)
     ! coefficient before the mineral is always 1
     real(r8) :: primary_stoi_co2            (1:nminerals)    ! reaction stoichiometry coefficient in front of CO2, 1:nminerals
     real(r8) :: primary_stoi_h2o            (1:nminerals)    ! reaction stoichiometry coefficient in front of water, 1:nminerals
     real(r8) :: primary_stoi_cations        (1:nminerals, 1:ncations)    ! reaction stoichiometry coefficient in front of cations, 1:nminerals x 1:ncations
     real(r8) :: primary_stoi_bicarbonate    (1:nminerals)    ! reaction stoichiometry coefficient in front of bicarbonate, 1:nminerals
     real(r8) :: primary_stoi_silicate       (1:nminerals)    ! reaction stoichiometry coefficient in front of silicate, 1:nminerals

     real(r8) :: primary_weight              (1:nminerals)    ! molecular weight of the primary mineral, g/mol, 1:nminerals (e.g. Mg2SiO4 = 140.6931 g/mol)
     real(r8) :: cation_weight               (1:ncations)    ! molecular weight of the cations in the minerals, g/mol, 1:ncations
     real(r8) :: cation_valence              (1:ncations)    ! valence of the cations (1+, 2+) in the minerals, g/mol, 1:ncations
     real(r8) :: secondary_weight            (1:nminsec)    ! molecular weight of the secondary mineral, g/mol, 1:nminsec (e.g. CaCO3 = 140.6931 g/mol)

     real(r8) :: k_secondary                 (1:nminsec)    ! secondary mineral precipitation reaction constant at 298.15K (the goal is to get reaction rate in the unit of mol reaction m-3 solution s-1), 1:nminsec
     real(r8) :: e_secondary                 (1:nminsec)    ! secondary mineral precipitation reaction constant at 298.15K (J mol-1), 1:nminsec

     ! reaction stoichiometry: suppose the equation is 
     ! pathway1: cation + bicarbonate -> carbonate + co2 + water
     ! coefficient before the cation is always 1
     ! pathway2 - clay formation: Al3+ + silicate + cations -> clay
     real(r8) :: secondary_stoi_cations     (1:nminsec,1:ncations)    ! reaction stoiciometry in front of the cations, 1:nminsec, 1:ncations
     real(r8) :: secondary_stoi_bicarbonate (1:nminsec)    ! reaction stoichiometry coefficient in front of bicarbonate, 1:nminsec
     real(r8) :: secondary_stoi_minsec      (1:nminsec)    ! reaction stoichiometry coefficient in front of secondary mineral, 1:nminsec
     real(r8) :: secondary_stoi_co2         (1:nminsec)    ! reaction stoichiometry coefficient in front of CO2, 1:nminsec
     real(r8) :: secondary_stoi_h2o         (1:nminsec)    ! reaction stoichiometry coefficient in front of water, 1:nminsec

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

    tString='k_primary'
    call ncd_io(varname=trim(tString),data=EWParamsInst%k_primary, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='e_primary'
    call ncd_io(varname=trim(tString),data=EWParamsInst%e_primary, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_stoi_co2'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_co2, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_stoi_h2o'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_h2o, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_stoi_cations'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_cations, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_stoi_bicarbonate'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_bicarbonate, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_stoi_silicate'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_stoi_silicate, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='primary_weight'
    call ncd_io(varname=trim(tString),data=EWParamsInst%primary_weight, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='cation_weight'
    call ncd_io(varname=trim(tString),data=EWParamsInst%cation_weight, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='cation_valence'
    call ncd_io(varname=trim(tString),data=EWParamsInst%cation_valence, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='secondary_weight'
    call ncd_io(varname=trim(tString),data=EWParamsInst%secondary_weight, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='k_secondary'
    call ncd_io(varname=trim(tString),data=EWParamsInst%k_secondary, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='e_secondary'
    call ncd_io(varname=trim(tString),data=EWParamsInst%e_secondary, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='secondary_stoi_cations'
    call ncd_io(varname=trim(tString),data=EWParamsInst%secondary_stoi_cations, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='secondary_stoi_bicarbonate'
    call ncd_io(varname=trim(tString),data=EWParamsInst%secondary_stoi_bicarbonate, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='secondary_stoi_minsec'
    call ncd_io(varname=trim(tString),data=EWParamsInst%secondary_stoi_minsec, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='secondary_stoi_co2'
    call ncd_io(varname=trim(tString),data=EWParamsInst%secondary_stoi_co2, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))

    tString='secondary_stoi_h2o'
    call ncd_io(varname=trim(tString),data=EWParamsInst%secondary_stoi_h2o, flag='read', ncid=ncid, readvar=readv)
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
    !$acc routine seq
    use elm_varcon      , only : spval, rgas, secspday, zisoi, dzsoi
    use elm_time_manager, only : get_step_size
    use landunit_varcon , only : istice, istwet, istsoil, istdlak, istcrop, istice_mec
    !
    ! !ARGUMENTS:
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    !
    ! !LOCAL VARIABLES:
    integer :: c, fc, l, lu, m, a
    integer :: maxlayer ! the deepest decomposition layer
    real(r8):: dt
    real(r8):: conc_cation_temp(1:ncations) ! temporary variable cation concentration in mol m-3 soil water
    real(r8):: conc_bicarbonate_temp ! temporary variable cation concentration in mol m-3 soil water    real(r8), pointer :: co2_primary_temp(1:nminerals) => null()
    real(r8) :: co2_primary_temp(1:nminerals)
    real(r8) :: h2o_primary_temp(1:nminerals)
    real(r8) :: cation_primary_temp(1:nminerals, 1:ncations)
    real(r8) :: bicarbonate_primary_temp(1:nminerals)
    real(r8) :: silicate_primary_temp(1:nminerals)
    real(r8) :: scale_r_precip_bycation(1:ncations) ! scaling down factors for precipitation
    real(r8) :: scale_r_precip_byminsec(1:nminsec) ! the above, minimum over all the used cations
    real(r8) :: cation_secondary_temp(1:nminsec, 1:ncations)
    real(r8) :: cation_secondary_tot_temp(1:ncations)
    real(r8) :: bicarbonate_secondary_temp(1:nminsec)
    real(r8) :: h2o_secondary_temp(1:nminsec)
    real(r8) :: co2_secondary_temp(1:nminsec)

    !-----------------------------------------------------------------------
    associate( &
         !
         ! Enhanced weathering parameter values
         !
         k_primary => EWParamsInst%k_primary, &
         e_primary => EWParamsInst%e_primary, &
         primary_stoi_co2 => EWParamsInst%primary_stoi_co2, &
         primary_stoi_h2o => EWParamsInst%primary_stoi_h2o, &
         primary_stoi_cations => EWParamsInst%primary_stoi_cations, &
         primary_stoi_bicarbonate => EWParamsInst%primary_stoi_bicarbonate, &
         primary_stoi_silicate => EWParamsInst%primary_stoi_silicate, &
         primary_weight => EWParamsInst%primary_weight, &
         cation_weight => EWParamsInst%cation_weight, &
         cation_valence => EWParamsInst%cation_valence, &
         secondary_weight => EWParamsInst%secondary_weight, &
         k_secondary => EWParamsInst%k_secondary, &
         e_secondary => EWParamsInst%e_secondary, &
         secondary_stoi_cations => EWParamsInst%secondary_stoi_cations, &
         secondary_stoi_bicarbonate => EWParamsInst%secondary_stoi_bicarbonate, &
         secondary_stoi_minsec => EWParamsInst%secondary_stoi_minsec, &
         secondary_stoi_co2 => EWParamsInst%secondary_stoi_co2, &
         secondary_stoi_h2o => EWParamsInst%secondary_stoi_h2o, &

         !
         ! Forcing variables
         !
         forc_app               => col_ew%forc_app, & ! Input:  [real(r8) (:)] application rate (kg rock m-2 year-1)
         forc_min               => col_ew%forc_min, & ! Input:  [real(r8) (:, 1:nminerals) weight percentage of minerals in rock (1:nminerals) (kg mineral kg-1 rock)
         forc_pho               => col_ew%forc_pho, & ! Input:  [real(r8) (:)] weight percentage of phosphorus content in rock (gP kg-1 rock)
         forc_gra               => col_ew%forc_gra, & ! Input:  [real(r8) (:)] grain size (um diameter)
         forc_sph               => col_ew%forc_sph, & ! Input:  [real(r8) (:)] soil pH

         !
         ! Mineral state
         !
         primary_mineral_vr     => col_ms%primary_mineral_vr   , &
         cation_vr              => col_ms%cation_vr            , &
         bicarbonate_vr         => col_ms%bicarbonate_vr       , &
         armor_thickness_vr     => col_ms%armor_thickness_vr   , &
         ssa                    => col_ms%ssa                  , &

         !
         ! Mineral flux
         !
         primary_added_vr               => col_mf%primary_added_vr     , &
         primary_dissolve_vr            => col_mf%primary_dissolve_vr  , &
         primary_co2_flux_vr            => col_mf%primary_co2_flux_vr  , &
         primary_h2o_flux_vr            => col_mf%primary_h2o_flux_vr  , &

         primary_cation_flux_vr         => col_mf%primary_cation_flux_vr     , &
         primary_bicarbonate_flux_vr    => col_mf%primary_bicarbonate_flux_vr, &
         primary_silicate_flux_vr       => col_mf%primary_silicate_flux_vr   , &

         r_dissolve_vr                  => col_mf%r_dissolve_vr, &

         secondary_cation_flux_vr       => col_mf%secondary_cation_flux_vr     , &
         secondary_bicarbonate_flux_vr  => col_mf%secondary_bicarbonate_flux_vr, &
         secondary_precip_vr            => col_mf%secondary_precip_vr          , &
         secondary_co2_flux_vr          => col_mf%secondary_co2_flux_vr        , &
         secondary_h2o_flux_vr          => col_mf%secondary_h2o_flux_vr        , &

         r_precip_vr                    => col_mf%r_precip_vr                  , &

         cation_leached_vr              => col_mf%cation_leached_vr           , &
         bicarbonate_leached_vr         => col_mf%bicarbonate_leached_vr      , &

         mineral_prelease_vr            => col_mf%mineral_prelease_vr          , &
         r_sequestration                => col_mf%r_sequestration              , &

         !
         ! Other related
         !
         tsoi                           => col_es%t_soisno, &
         h2osoi_liqvol                  => col_ws%h2osoi_liqvol, &
         dz                             => col_pp%dz & ! Input:  [real(r8) (:,:) ]  layer thickness (m)
    )

    dt      = real( get_step_size(), r8 )

    do fc = 1,num_soilc
      c = filter_soilc(fc)

      !! TODO: hard code the forcing values here
      !forc_app(c) = 0.004_r8 / 365._r8 / 86400._r8 ! 40 ton/ha = 4 kg m-2, convert from yr-1 to s-1
      !forc_min(c,1) = 1._r8
      !forc_pho(c) = 0._r8 ! no P content yet
      !forc_gra(c) = 100._r8 ! 100 microns in diameter
      !forc_sph(c) = 7._r8 ! soil pH, not used yet

      write(iulog, *) forc_app(c), forc_min(c, 1), forc_pho(c), forc_gra(c), forc_sph(c)

      ! Specific surface area depends on the grain size of the mineral
      ! Strefler, J., Amann, T., Bauer, N., Kriegler, E., and Hartmann, J.: Potential and costs of carbon dioxide removal by enhanced weathering of rocks, Environ. Res. Lett., 13, 034010, https://doi.org/10.1088/1748-9326/aaa9c4, 2018.
      ! More accurate method from geochemistry: 
      !   at ~100 um magnitude, the grains are individual minerals
      !   weighted average of the specific surface area of each mineral (m^2 g-1)
      ssa(c) = 69.18_r8 * (forc_gra(c) ** (-1.24_r8)) * 1000._r8 ! unit: m^2 kg-1 mineral

      ! find the maximum layer of reaction
      maxlayer = 1
      do l = 2,nlevgrnd
        if (zisoi(l-1) < mixing_depth) then
          maxlayer = l
        end if
      end do

      do l = 1,maxlayer
        !--------------------------------------------------------------------------------------
        ! Dissolution reaction
        !--------------------------------------------------------------------------------------
        do m = 1,nminerals
          ! TODO: the reaction rates equation may use fitted temperature response function 
          !       instead of activation energy

          ! evenly distributed in the top 30 centimeters
          if (l < maxlayer) then
            primary_added_vr(c,l,m) = forc_app(c) * forc_min(c,m) * dz(c,l) / mixing_depth
          else
            primary_added_vr(c,l,m) = forc_app(c) * forc_min(c,m) * & 
                (mixing_depth - zisoi(l-1)) / mixing_depth
          end if

          ! reaction rate in the unit of mol reaction kg-1 mineral s-1
          ! rgas need to be divided by 1000 because the unit in ELM is J/kmol/k
          ! not using soil pH & water for now
          r_dissolve_vr(c,l,m) = ssa(c) * k_primary(m) * exp( - e_primary(m) / (rgas/1000.) * (1/tsoi(c, l) - 1/298.15_r8) )

          ! convert to individual fluxes
          ! kg m-2 area s-1 = stoichiometry * mol reaction kg-1 mineral s-1 * kg mineral m-2 area s-1 * g mol-1 / 1000. for unit conversion
          primary_dissolve_vr(c,l,m) = r_dissolve_vr(c,l,m) * primary_mineral_vr(c,l,m) * primary_weight(m) / 1000._r8

          co2_primary_temp(m) = primary_stoi_co2(m) * r_dissolve_vr(c,l,m) * primary_mineral_vr(c,l,m) * co2_weight / 1000._r8
          h2o_primary_temp(m) = primary_stoi_h2o(m) * r_dissolve_vr(c,l,m) * primary_mineral_vr(c,l,m) * h2o_weight / 1000._r8

          do a = 1,ncations
            cation_primary_temp(m,a) = primary_stoi_cations(m,a) * r_dissolve_vr(c,l,m) * primary_mineral_vr(c,l,m) * cation_weight(a) / 1000._r8
          end do

          bicarbonate_primary_temp(m) = primary_stoi_bicarbonate(m) * r_dissolve_vr(c,l,m) * primary_mineral_vr(c,l,m) * bicarbonate_weight / 1000._r8

          silicate_primary_temp(m) = primary_stoi_silicate(m) * r_dissolve_vr(c,l,m) * primary_mineral_vr(c,l,m) * silicate_weight / 1000._r8
        end do

        primary_co2_flux_vr(c,l) = 0._r8
        primary_h2o_flux_vr(c,l) = 0._r8
        primary_cation_flux_vr(c,l,:) = 0._r8
        primary_bicarbonate_flux_vr(c,l) = 0._r8
        primary_silicate_flux_vr(c,l) = 0._r8
        do m = 1,nminerals
          primary_co2_flux_vr(c,l) = primary_co2_flux_vr(c,l) + co2_primary_temp(m)
          primary_h2o_flux_vr(c,l) = primary_h2o_flux_vr(c,l) + h2o_primary_temp(m)
          do a = 1,ncations
            primary_cation_flux_vr(c,l,a) = primary_cation_flux_vr(c,l,a) + cation_primary_temp(m,a)
          end do
          primary_bicarbonate_flux_vr(c,l) = primary_bicarbonate_flux_vr(c,l) + bicarbonate_primary_temp(m)
          primary_silicate_flux_vr(c,l) = primary_silicate_flux_vr(c,l) + silicate_primary_temp(m)
        end do

        !--------------------------------------------------------------------------------------
        ! Secondary mineral reaction          !--------------------------------------------------------------------------------------
        ! convert cation concentration in kg m-2 to concentration in soil solution (mol m-3)
        ! need a catch for zero soil water
        do a = 1,ncations
          if (h2osoi_liqvol(c,l) > 0._r8) then
            conc_cation_temp(a) = cation_vr(c,l,a) * 1000._r8 / cation_weight(a) / h2osoi_liqvol(c,l) / dz(c,l)
          else
            conc_cation_temp(a) = 0._r8
          end if
        end do
        if (h2osoi_liqvol(c,l) > 0._r8) then
          conc_bicarbonate_temp = bicarbonate_vr(c,l) * 1000._r8 / bicarbonate_weight / h2osoi_liqvol(c,l) / dz(c,l)
        else
          conc_bicarbonate_temp = 0._r8
        end if

        ! First pass: get the potential rate
        do m = 1,nminsec

          ! rgas need to be divided by 1000 because the unit in ELM is J/kmol/k
          ! TODO: not using soil pH & water & armoring layer thickness for now
          ! suppose the reaction rate is proportional to surface area instead of concentration in the solution
          r_precip_vr(c,l,m) = ssa(c) * k_secondary(m) * exp( - e_secondary(m) / (rgas/1000.) * (1/tsoi(c, l) - 1/298.15_r8) ) ! * conc_bicarbonate_temp

          !! suppose the reaction rate is in the unit of mol reaction m-3 soil solution s-1
          !do a = 1,ncations
          !  if (secondary_stoi_cations(m,a) > 0._r8) then
          !    r_precip_vr(c,l,m) = r_precip_vr(c,l,m) * conc_cation_temp(a)
          !  end if
          !end do

          ! convert to individual fluxes
          ! kg m-2 area s-1 = stoichiometry * mol reaction m-3 soil solution s-1 * soil water content (m3 m-3) * soil layer thickness (m) * g mol-1 / 1000. for unit conversion
          do a = 1,ncations
            cation_secondary_temp(m,a) = secondary_stoi_cations(m,a) * r_precip_vr(c,l,m) * h2osoi_liqvol(c,l) * dz(c,l) * cation_weight(a)
          end do
          bicarbonate_secondary_temp(m) = secondary_stoi_bicarbonate(m) * r_precip_vr(c,l,m) * h2osoi_liqvol(c,l) * dz(c,l) * bicarbonate_weight
        end do

        ! Calculate the total amount of cations required
        do a = 1,ncations
          cation_secondary_tot_temp(a) = 0._r8
          do m = 1,nminsec
            cation_secondary_tot_temp(a) = cation_secondary_tot_temp(a) + cation_secondary_temp(m,a)
          end do

          ! scale down if the total needed is larger than available
          if ((cation_secondary_tot_temp(a) * dt) > cation_vr(c,l,a)) then
            scale_r_precip_bycation(a) = max(0._r8, cation_vr(c,l,a) / cation_secondary_tot_temp(a) / dt)
          else
            scale_r_precip_bycation(a) = 1._r8 ! no cation rqeuired
          end if
        end do

        ! Compare to the available cations
        do m = 1,nminsec
          scale_r_precip_byminsec(m) = 1._r8
          do a = 1,ncations
            if (secondary_stoi_cations(m,a) > 0._r8) then
              scale_r_precip_byminsec(m) = min(scale_r_precip_byminsec(m), &
                                               scale_r_precip_bycation(a))
            end if
          end do
        end do

        ! Scale down
        do m = 1,nminsec
          r_precip_vr(c,l,m) = r_precip_vr(c,l,m) * scale_r_precip_byminsec(m)
          do a = 1,ncations
            cation_secondary_temp(m,a) = secondary_stoi_cations(m,a) * r_precip_vr(c,l,m) * h2osoi_liqvol(c,l) * dz(c,l) * cation_weight(a)
          end do
          bicarbonate_secondary_temp(m) = secondary_stoi_bicarbonate(m) * r_precip_vr(c,l,m) * h2osoi_liqvol(c,l) * dz(c,l) * bicarbonate_weight
        end do

        ! Integrate over minerals
        do a = 1,ncations
          secondary_cation_flux_vr(c,l,a) = 0._r8
          do m = 1,nminsec
            secondary_cation_flux_vr(c,l,a) = cation_secondary_temp(m,a)
          end do
        end do
        secondary_bicarbonate_flux_vr(c,l) = 0._r8
        secondary_co2_flux_vr(c,l) = 0._r8
        secondary_h2o_flux_vr(c,l) = 0._r8
        do m = 1,nminsec
          secondary_bicarbonate_flux_vr(c,l) = secondary_bicarbonate_flux_vr(c,l) + bicarbonate_secondary_temp(m)
          secondary_co2_flux_vr(c,l) = secondary_co2_flux_vr(c,l) + co2_secondary_temp(m)
          secondary_h2o_flux_vr(c,l) = secondary_h2o_flux_vr(c,l) + h2o_secondary_temp(m)
        end do
      end do

      ! calculate the total CO2 sequestration rate
      r_sequestration(c) = 0._r8
      do l = 1,nlevgrnd
        r_sequestration(c) = r_sequestration(c) + primary_co2_flux_vr(c,l) - secondary_co2_flux_vr(c,l)
      end do
      ! convert from kg CO2 to gC
      r_sequestration(c) = r_sequestration(c) * 1000._r8 * 12._r8 / co2_weight

      ! TODO: use cation valence to further reduce by ocean efficiency

      ! calculate the P-release rate (gP m-2 s-1)
      do l = 1,nlevgrnd
        mineral_prelease_vr(c,l) = 0._r8
        do m = 1, nminerals
          mineral_prelease_vr(c,l) = primary_dissolve_vr(c,l,m) * forc_pho(c)
        end do
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

         cation_vr              => col_ms%cation_vr                       , & ! Input:  [real(r8) (:,:,:) ] cation in each soil layer (kg m-2)
         bicarbonate_vr         => col_ms%bicarbonate_vr                  , & ! Input:  [real(r8) (:,:,:) ] bicarbonate in each soil layer (kg m-2)

         ! TODO: change the unit to m^-3
         cation_leached_vr      => col_mf%cation_leached_vr               , & ! Output: [real(r8) (:,:,:) ]  rate of cation leaching (kg m-2 s-1) ! compare to N (g m-3 s-1)
         bicarbonate_leached_vr => col_mf%bicarbonate_leached_vr          , & ! Output: [real(r8) (:,:) ]  rate of cation leaching (kg m-2 s-1) ! compare to N (g m-3 s-1)

         cation_runoff_vr       => col_mf%cation_runoff_vr                , & ! Output: [real(r8) (:,:,:) ]  rate of cation loss with runoff (kg m-2 s-1) ! compare to N (g m-3 s-1)
         bicarbonate_runoff_vr  => col_mf%bicarbonate_runoff_vr             & ! Output: [real(r8) (:,:) ]  rate of bicarbonate loss with runoff (kg m-2 s-1) ! compare to N (g m-3 s-1)
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

    end associate
  end subroutine MineralLeaching

end module EnhancedWeatheringMod