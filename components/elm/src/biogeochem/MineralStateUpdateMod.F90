module MineralStateUpdateMod

  !-----------------------------------------------------------------------
  ! Module for enhanced weathering state variable update.
  !
  ! !USES:
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use decompMod               , only : bounds_type
  use spmdMod                 , only : iam
  use elm_varpar              , only : nminerals, ncations, nminsecs, nlevgrnd, nlevsoi
  use elm_varcon              , only : zisoi, dzsoi, mass_h, mass_hco3, mass_co3, secspday
  use elm_varctl              , only : iulog, use_erw_verbose
  use shr_sys_mod             , only : shr_sys_flush
  use spmdMod                 , only : masterproc
  use abortutils              , only : endrun
  use shr_log_mod             , only : errMsg => shr_log_errMsg
  use ewutils                 , only : mass_to_mol, mass_to_meq, mol_to_mass, meq_to_mass
  use ColumnDataType          , only : col_ws, col_ms, col_mf, col_pp, col_ew
  use ColumnDataType          , only : column_mineral_state, column_mineral_flux, column_water_flux, column_erw_forcing
  use SoilStateType           , only : soilstate_type
  use CNStateType             , only : cnstate_type
  use EnhancedWeatheringMod   , only : EWParamsInst
  use domainMod               , only : ldomain ! debug print
  use elm_time_manager        , only : get_curr_time_string
  use timeinfoMod
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: MineralFluxLimit
  public :: MineralStateUpdate1
  public :: MineralStateUpdate2
  public :: MineralStateUpdate3
  public :: MineralSelfCalibrate
  public :: MineralStateDiags
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine MineralStateUpdate1(num_soilc, filter_soilc, col_ms, col_mf, col_ew, dt, soilstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the mineral state variables that are not
    ! affected by vertical or horizontal soil water movement, update pH-dependent CEC.
    !
    !$acc routine seq
    !
    ! !USES:
    use ewutils          , only : u_pdf, get_ssa
    !
    ! !ARGUMENTS:
    integer                      , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                      , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(column_erw_forcing)     , intent(inout) :: col_ew
    type(column_mineral_state)   , intent(inout) :: col_ms
    type(column_mineral_flux)    , intent(inout) :: col_mf
    real(r8)                     , intent(in)    :: dt              ! radiation time step (seconds)
    type(soilstate_type)         , intent(in)    :: soilstate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,k,icat,m,g       ! indices
    integer  :: fp,fc                  ! lake filter indices
    integer  :: nlevbed
    real(r8) :: temp_gra               ! temporary variable for grain size
    real(r8) :: theta                  ! parameter of the gamma distribution
    real(r8) :: num, prob              ! helper variables to integrate over distribution
    real(r8) :: u_min, u_max           ! bounds of integration
    real(r8) :: u, du                  ! discretization step (size)
    integer, parameter :: n_int = 1000 ! number of discretized intervals
    integer  :: ii                     ! iterator for discretization
    !-----------------------------------------------------------------------

    ! Update mineral state
    do fc = 1,num_soilc
      c = filter_soilc(fc)
      nlevbed = min(col_pp%nlevbed(c), nlevsoi)

      do j = 1,nlevbed
        ! soil H+ concentration
        ! only determined by CEC equilibrium
        col_ms%proton_vr(c,j) = mol_to_mass(10**(-col_ms%soil_ph(c,j)), mass_h, col_ws%h2osoi_vol(c,j))

        ! soil cation concentration - not updated here
        ! must be preserved before calling the vertical solute movement solver

        ! pH-dependent CEC
        col_ms%cect_dyn(c,j) = col_mf%cect_delta(c,j) + col_ms%cect_dyn(c,j)

        ! CEC cations - only depends on flux limit
        do icat = 1,ncations
          col_ms%cec_cation_vr(c,j,icat) = col_ms%cec_cation_vr(c,j,icat) + &
                  (- col_mf%cec_cation_flux_vr(c,j,icat) - col_mf%cec_cation_flux2_vr(c,j,icat) &
                   + col_mf%background_cec_vr(c,j,icat))*dt
        end do

        ! CEC H+ with pH dependent changes
        ! the Equilibria subroutine cannot distinguish the effects of CO2 and cation exchange
        ! instead, use charge balance on the mineral surface to get the change in adsorped H+
        ! note the newly exposed surface due to dynamic CEC
        col_ms%cec_proton_vr(c,j) = col_ms%cect_dyn(c,j)
        do icat = 1,ncations
          col_ms%cec_proton_vr(c,j) = col_ms%cec_proton_vr(c,j) - &
            mass_to_meq(col_ms%cec_cation_vr(c,j,icat), EWParamsInst%cations_valence(icat), &
                        EWParamsInst%cations_mass(icat), soilstate_vars%bd_col(c,j))
        end do
        col_ms%cec_proton_vr(c,j) = meq_to_mass(col_ms%cec_proton_vr(c,j), 1._r8, mass_h, soilstate_vars%bd_col(c,j))

        ! primary mineral
        do m = 1,nminerals
          col_ms%primary_mineral_vr(c,j,m) = col_ms%primary_mineral_vr(c,j,m) + &
            col_mf%primary_added_vr(c,j,m)*dt - col_mf%primary_dissolve_vr(c,j,m)*dt

          ! catch tiny negative values
          col_ms%primary_mineral_vr(c,j,m) = max(col_ms%primary_mineral_vr(c,j,m), 0._r8)

          ! non-SiO2 minerals
          col_ms%primary_residue_vr(c,j,m) = col_ms%primary_residue_vr(c,j,m) + &
            col_mf%primary_residue_flux_vr(c,j,m)*dt
        end do

        ! silicate
        col_ms%silica_vr(c,j) = col_ms%silica_vr(c,j) + col_mf%primary_silica_flux_vr(c,j) * dt + col_mf%secondary_silica_flux_vr(c,j) * dt

        ! secondary mineral
        do m = 1,nminsecs
            col_ms%secondary_mineral_vr(c,j,m) = col_ms%secondary_mineral_vr(c,j,m) + & 
              col_mf%secondary_mineral_flux_vr(c,j,m)*dt + col_mf%background_minsecs_vr(c,j,m)*dt
        end do

        ! ignore the effect on soil water content for now

        ! passivation layer thickness
        do m = 1,nminerals
          col_ms%passivation_thickness(c,j,m) = col_ms%passivation_thickness(c,j,m) + &
            col_mf%passivation_rate(c,j,m) * dt
        end do

        ! specific surface area
        ! Specific surface area depends on the grain size of the mineral
        !    Strefler, J., Amann, T., Bauer, N., Kriegler, E., and Hartmann, J.: Potential and 
        !       costs of carbon dioxide removal by enhanced weathering of rocks, Environ. Res.
        !       Lett., 13, 034010, https://doi.org/10.1088/1748-9326/aaa9c4, 2018.
        do m = 1,nminerals
          temp_gra = col_ew%forc_gra(c,m) - col_ms%passivation_thickness(c,j,m)*1e6_r8 ! um

          if (builtin_site == 1) then
            ! Hubbard Brook applied in the form of pellets
            col_ms%ssa_dyn(c,j,m) = get_ssa(temp_gra) ! unit: m^2 g-1
          else
            ! integrate over grain size distribution to calculate the surface area
            ! range from 0.01 * grain size to 100 * grain size
            theta = 4.395_r8 * temp_gra
            u_min = sqrt(0.01_r8 * temp_gra)
            u_max = sqrt(100._r8 * temp_gra)
            du = (u_max - u_min) / n_int

            num = get_ssa(u_min**2) * u_pdf(u_min, theta)
            prob = u_pdf(u_min, theta)

            do ii = 1, n_int-1
              u = u_min + du * ii
              if (mod(ii,2) == 1) then
                num = num + 4_r8 * get_ssa(u**2) * u_pdf(u, theta)
                prob = prob + 4_r8 * u_pdf(u, theta)
              else
                num = num + 2_r8 * get_ssa(u**2) * u_pdf(u, theta)
                prob = prob + 2_r8 * u_pdf(u, theta)
              end if
            end do

            num = num + get_ssa(u_max**2) * u_pdf(u_max, theta)
            prob = prob + u_pdf(u_max, theta)

            num = num * du / 3._r8
            prob = prob * du / 3._r8

            col_ms%ssa_dyn(c,j,m) = num / prob

            !!write (iulog, *) u_min, u_max, du, num, prob, ssa(c,m)

            !! Add a roughness factor: lambda = (10^10 * r [m])**0.33
            !!    Kanzaki et al. (2022) Soil Cycles of Elements simulator for Predicting TERrestrial
            !!        regulation of greenhouse gases: SCEPTER v0.9. 
            !!        https://doi.org/10.5194/gmd-15-4959-2022         Eq. 39
            !!    Beerling, D. J., Kantzas, E. P., Lomas, M. R., Wade, P., Eufrasio, R. M., Renforth, 
            !!        P., et al. (2020). Potential for large-scale CO2 removal via enhanced rock 
            !!        weathering with croplands. Nature, 583(7815), 242–248. 
            !!        https://doi.org/10.1038/s41586-020-2448-9        SI Eq. 8 
            !!ssa_dyn(c,m) = (1e4_r8 * temp_gra) ** 0.33_r8 * ssa_dyn(c,j,m)
          end if
        end do
      end do
    end do
  end subroutine MineralStateUpdate1


  !-----------------------------------------------------------------------
  subroutine MineralStateUpdate2(num_soilc, filter_soilc, col_ms, col_mf, dt)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the mineral state variables
    ! related to vertical water movement
    !$acc routine seq
    ! !ARGUMENTS:
    integer                      , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                      , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(column_mineral_state)   , intent(inout) :: col_ms
    type(column_mineral_flux)    , intent(inout) :: col_mf
    real(r8)                     , intent(in)    :: dt              ! radiation time step (seconds)

    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,k,icat,m,g ! indices
    integer  :: fp,fc         ! lake filter indices
    integer  :: nlevbed
    !-----------------------------------------------------------------------

    ! Update mineral state
    do fc = 1,num_soilc
      c = filter_soilc(fc)
      nlevbed = min(col_pp%nlevbed(c), nlevsoi)
      do icat = 1,ncations
        do j = 1,nlevbed
          ! note the source sink terms are called in the advection_diffusion solver
          col_ms%cation_vr(c,j,icat) = col_ms%cation_vr(c,j,icat) + & 
            ( col_mf%background_flux_vr(c,j,icat) + & 
              col_mf%primary_cation_flux_vr(c,j,icat) + & 
              col_mf%cec_cation_flux_vr(c,j,icat) + &
              col_mf%cec_cation_flux2_vr(c,j,icat) + &
              col_mf%secondary_cation_flux_vr(c,j,icat) - & 
              col_mf%cation_uptake_vr(c,j,icat) ) * dt + & 
            ( col_mf%cation_infl_vr(c,j,icat) - col_mf%cation_oufl_vr(c,j,icat) ) * dt
        end do
      end do
    end do
  end subroutine MineralStateUpdate2

  !-----------------------------------------------------------------------
  subroutine MineralStateUpdate3(num_soilc, filter_soilc, col_ms, col_mf, dt)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the mineral state variables
    ! related to horizontal soil water movement.
    ! Also update carbon sequestration rate.
    !$acc routine seq
    ! !ARGUMENTS:
    integer                      , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                      , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(column_mineral_state)   , intent(inout) :: col_ms
    type(column_mineral_flux)    , intent(inout) :: col_mf
    real(r8)                     , intent(in)    :: dt              ! radiation time step (seconds)

    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,icat,g    ! indices
    integer  :: fc            ! lake filter indices
    integer  :: nlevbed
    integer  :: rmethod       ! 1 - use cation, 2 - use HCO3-- and CO3-- flux
    !-----------------------------------------------------------------------

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      nlevbed = min(col_pp%nlevbed(c), nlevsoi)

      ! Update mineral state
      do j = 1,nlevbed
        do icat = 1,ncations
          col_ms%cation_vr(c,j,icat) = col_ms%cation_vr(c,j,icat) - col_mf%cation_leached_vr(c,j,icat)*dt - col_mf%cation_runoff_vr(c,j,icat)*dt
        end do
      end do

      ! Calculate the total CO2 sequestration rate in mol m-2 s-1
      rmethod = 2
      if (rmethod == 1) then

        col_mf%r_sequestration(c) = 0._r8
        do j = 1,nlevbed
          ! precipitated by calcite: 1 mol CO2 per mol Ca2+
          col_mf%r_sequestration(c) = col_mf%r_sequestration(c) + & 
            col_mf%secondary_cation_flux_vr(c,j,1) / EWParamsInst%cations_mass(1) * col_pp%dz(c,j)
          ! transported to ocean: 2x for 2+ cations, 1x for 1+ cations, multiply by
          ! ocean efficiency (0.86)
          ! - col_mf%background_flux_vr(c,j,icat)
          do icat = 1,ncations-1
            col_mf%r_sequestration(c) = col_mf%r_sequestration(c) + & 
                ( col_mf%cation_leached_vr(c,j,icat) + col_mf%cation_runoff_vr(c,j,icat) - &
                  col_mf%cation_infl_vr(c,j,icat) ) * 0.86_r8 * col_pp%dz(c,j) / &
                EWParamsInst%cations_mass(icat) * EWParamsInst%cations_valence(icat)
          end do
        end do

      else

        ! calculate the total CO2 sequestration rate in mol m-2 s-1 as the 
        ! bottom drainage of HCO3- + 2*CO3--
        col_mf%r_sequestration(c) = (col_mf%bicarbonate_drainage(c) / mass_hco3 + &
          col_mf%carbonate_drainage(c) / mass_co3 * 2._r8) * 0.86_r8
        do j = 1,nlevbed
          ! add the subsurface drainage
          col_mf%r_sequestration(c) = (col_mf%r_sequestration(c) + &
            col_mf%bicarbonate_leached_vr(c,j) * col_pp%dz(c,j) / mass_hco3 + &
            col_mf%carbonate_leached_vr(c,j) * col_pp%dz(c,j) / mass_co3 * 2._r8) * 0.86_r8

          ! add the precipitated by calcite: 1 mol CO2 per mol Ca2+
          col_mf%r_sequestration(c) = col_mf%r_sequestration(c) + & 
            col_mf%secondary_cation_flux_vr(c,j,1) / EWParamsInst%cations_mass(1) * col_pp%dz(c,j)
        end do

      end if

      ! convert from mol m-2 s-1 to gC m-2 s-1
      col_mf%r_sequestration(c) = col_mf%r_sequestration(c) * 12._r8
    end do
  end subroutine MineralStateUpdate3

  !-----------------------------------------------------------------------
  subroutine MineralSelfCalibrate(num_soilc, filter_soilc, col_ms, col_mf, col_wf, dt)
    !
    ! !DESCRIPTION:
    ! During the self-calibration stage, calculate the amount of background weathering flux
    ! needed to replenish the cation lost from the system
    ! !ARGUMENTS:
    integer                      , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                      , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(column_mineral_state)   , intent(inout) :: col_ms
    type(column_mineral_flux)    , intent(inout) :: col_mf
    type(column_water_flux)      , intent(in)    :: col_wf
    real(r8)                     , intent(in)    :: dt              ! radiation time step (seconds)
    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,icat,g,isec    ! indices
    integer  :: fc                 ! lake filter indices
    integer  :: nlevbed
    real(r8) :: step_delta
    real(r8) :: fracday

    fracday = dt / secspday

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      nlevbed = min(col_pp%nlevbed(c), nlevsoi)

      ! Update the annual average flow rate accumulator
      do j = 1,nlevbed+1
        col_wf%tempavg_qin_col(c,j) = col_wf%tempavg_qin_col(c,j) + abs(col_wf%qin(c,j)) * fracday / dayspyr_mod
      end do

      do j = 1,nlevbed
        do icat = 1,ncations
          ! Update the annual total column cation loss rate accumulator
          ! (when background flux does not exist)
          step_delta = col_mf%cec_cation_flux_vr(c,j,icat) + &
                       col_mf%cec_cation_flux2_vr(c,j,icat) + &
                       col_mf%secondary_cation_flux_vr(c,j,icat) - &
                       col_mf%cation_uptake_vr(c,j,icat) + &
                       col_mf%cation_infl_vr(c,j,icat) - & 
                       col_mf%cation_oufl_vr(c,j,icat) - &
                       col_mf%cation_leached_vr(c,j,icat) - &
                       col_mf%cation_runoff_vr(c,j,icat)

          col_mf%tempavg_tot_delta(c,j,icat) = col_mf%tempavg_tot_delta(c,j,icat) + &
            step_delta * fracday / dayspyr_mod

          ! Also calibrate a replenishment term to the cation exchange phase
          ! because it seems to be lost pretty severly
          ! note: cec_cation_flux_vr is defined negative for adsorption into soil
          !       this term is defined positive for adsorption into soil
          col_mf%tempavg_cec_delta(c,j,icat) = col_mf%tempavg_cec_delta(c,j,icat) - &
            (col_mf%cec_cation_flux_vr(c,j,icat) + col_mf%cec_cation_flux2_vr(c,j,icat) &
            ) * fracday / dayspyr_mod
        end do

        do isec = 1,nminsecs
          ! Also update the annual total column secondary mineral loss rate accumulator
          ! (when background flux does not exist)
          col_mf%tempavg_minsecs_delta(c,j,isec) = col_mf%tempavg_minsecs_delta(c,j,isec) + &
            col_mf%secondary_mineral_flux_vr(c,j,isec) * fracday / dayspyr_mod
        end do
      end do
    end do

  end subroutine MineralSelfCalibrate


  !-----------------------------------------------------------------------
  subroutine MineralStateDiags(num_soilc, filter_soilc, col_ms, col_mf, dt, soilstate_vars)
    !
    ! !DESCRIPTION:
    ! Write out diagnostics of mineral state depending on verbosity
    !$acc routine seq
    ! !ARGUMENTS:
    integer                      , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                      , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(column_mineral_state)   , intent(inout) :: col_ms
    type(column_mineral_flux)    , intent(inout) :: col_mf
    real(r8)                     , intent(in)    :: dt              ! radiation time step (seconds)
    type(soilstate_type)         , intent(in)    :: soilstate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,icat,g    ! indices
    integer  :: fc            ! lake filter indices
    integer  :: nlevbed
    character(len=256) :: dateTimeString

    !-----------------------------------------------------------------------
    call get_curr_time_string(dateTimeString)

    !-----------------------------------------------------------------------
    ! Write mass balance diagnostics only if verbose level is set to high

    !if ((ldomain%latc(g) > 44.7) .and. (ldomain%latc(g) < 44.8) .and. &
    !    (ldomain%lonc(g) > -69.8) .and. (ldomain%lonc(g) < -69.7)) then

    if (use_erw_verbose == 2) then
      ! print soil solution proton
      do fc = 1,num_soilc
        c = filter_soilc(fc)
        g = col_pp%gridcell(c)
        nlevbed = min(col_pp%nlevbed(c), nlevsoi)
        write (100+iam, *) 'Post-reaction H+: ', ldomain%latc(g), ldomain%lonc(g), trim(dateTimeString)
        do j = 1,nlevbed
          write (100+iam, *) c, j, col_ms%soil_ph(c,j), col_ms%proton_vr(c,j), & 
             mass_to_mol(col_ms%proton_vr(c,j), mass_h, col_ws%h2osoi_vol(c,j))
        end do
      end do

      ! print soil solution cations
      do fc = 1,num_soilc
        c = filter_soilc(fc)
        g = col_pp%gridcell(c)
        nlevbed = min(col_pp%nlevbed(c), nlevsoi)
        write (100+iam, *) 'Post-reaction cation: ', ldomain%latc(g), ldomain%lonc(g), trim(dateTimeString)
        ! note: sourcesink_cations term in EnhancedWeatheringMod.F90
        !       should be approximately matched to cation_infl_vr
        ! leaching & runoff are the additionals
        do j = 1,nlevbed
          do icat = 1,ncations
            write (100+iam, *) c, j, icat, col_ms%cation_vr(c,j,icat), &
              mass_to_mol(col_ms%cation_vr(c,j,icat), EWParamsInst%cations_mass(icat), &
                          col_ws%h2osoi_vol(c,j)), &
              col_mf%background_flux_vr(c,j,icat)*dt, &
              col_mf%primary_cation_flux_vr(c,j,icat)*dt, &
              col_mf%cec_cation_flux_vr(c,j,icat)*dt, & 
              col_mf%cec_cation_flux2_vr(c,j,icat)*dt, & 
              col_mf%secondary_cation_flux_vr(c,j,icat)*dt, &
              - col_mf%cation_uptake_vr(c,j,icat)*dt, &
              col_mf%cation_infl_vr(c,j,icat)*dt, &
              - col_mf%cation_leached_vr(c,j,icat)*dt, &
              - col_mf%cation_runoff_vr(c,j,icat)*dt
            end do
          end do
        end do

      ! print CEC protons
      do fc = 1,num_soilc
        c = filter_soilc(fc)
        g = col_pp%gridcell(c)
        nlevbed = min(col_pp%nlevbed(c), nlevsoi)
        write (100+iam, *) 'Post-reaction cec H+: ', ldomain%latc(g), ldomain%lonc(g), trim(dateTimeString)
        do j = 1,nlevbed
          ! note: the change in CEC H+ is equal to the sum of other cations' influx
          write (100+iam, *) c, j, col_ms%cec_proton_vr(c,j), &
            mass_to_meq(col_ms%cec_proton_vr(c,j), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
            col_mf%cect_delta(c,j), &
            mass_to_meq(col_mf%cec_cation_flux_vr(c,j,1)*dt/EWParamsInst%cations_mass(1) &
              *mass_h*EWParamsInst%cations_valence(1), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
            mass_to_meq(col_mf%cec_cation_flux_vr(c,j,2)*dt/EWParamsInst%cations_mass(2) &
              *mass_h*EWParamsInst%cations_valence(2), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
            mass_to_meq(col_mf%cec_cation_flux_vr(c,j,3)*dt/EWParamsInst%cations_mass(3) &
              *mass_h*EWParamsInst%cations_valence(3), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
            mass_to_meq(col_mf%cec_cation_flux_vr(c,j,4)*dt/EWParamsInst%cations_mass(4) &
              *mass_h*EWParamsInst%cations_valence(4), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
            mass_to_meq(col_mf%cec_cation_flux_vr(c,j,5)*dt/EWParamsInst%cations_mass(5) &
              *mass_h*EWParamsInst%cations_valence(5), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
            mass_to_meq(col_mf%cec_cation_flux2_vr(c,j,1)*dt/EWParamsInst%cations_mass(1) &
              *mass_h*EWParamsInst%cations_valence(1), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
            mass_to_meq(col_mf%cec_cation_flux2_vr(c,j,2)*dt/EWParamsInst%cations_mass(2) &
              *mass_h*EWParamsInst%cations_valence(2), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
            mass_to_meq(col_mf%cec_cation_flux2_vr(c,j,3)*dt/EWParamsInst%cations_mass(3) &
              *mass_h*EWParamsInst%cations_valence(3), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
            mass_to_meq(col_mf%cec_cation_flux2_vr(c,j,4)*dt/EWParamsInst%cations_mass(4) &
              *mass_h*EWParamsInst%cations_valence(4), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
            mass_to_meq(col_mf%cec_cation_flux2_vr(c,j,5)*dt/EWParamsInst%cations_mass(5) &
              *mass_h*EWParamsInst%cations_valence(5), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
            mass_to_meq(-col_mf%background_cec_vr(c,j,1)*dt/EWParamsInst%cations_mass(1) &
              *mass_h*EWParamsInst%cations_valence(1), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
            mass_to_meq(-col_mf%background_cec_vr(c,j,2)*dt/EWParamsInst%cations_mass(2) &
              *mass_h*EWParamsInst%cations_valence(2), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
            mass_to_meq(-col_mf%background_cec_vr(c,j,3)*dt/EWParamsInst%cations_mass(3) &
              *mass_h*EWParamsInst%cations_valence(3), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
            mass_to_meq(-col_mf%background_cec_vr(c,j,4)*dt/EWParamsInst%cations_mass(4) &
              *mass_h*EWParamsInst%cations_valence(4), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
            mass_to_meq(-col_mf%background_cec_vr(c,j,5)*dt/EWParamsInst%cations_mass(5) &
              *mass_h*EWParamsInst%cations_valence(5), 1._r8, mass_h, soilstate_vars%bd_col(c,j))
        end do
      end do

      ! print CEC cations
      do fc = 1,num_soilc
        c = filter_soilc(fc)
        g = col_pp%gridcell(c)
        nlevbed = min(col_pp%nlevbed(c), nlevsoi)
        write (100+iam, *) 'Post-reaction cec cation: ', ldomain%latc(g), ldomain%lonc(g), trim(dateTimeString)
        do j = 1,nlevbed
          do icat = 1,ncations
            write (100+iam, *) c, j, icat, col_ms%cec_cation_vr(c,j,icat), &
              mass_to_meq(col_ms%cec_cation_vr(c,j,icat), EWParamsInst%cations_valence(icat), &
                          EWParamsInst%cations_mass(icat), soilstate_vars%bd_col(c,j)), &
              -col_mf%cec_cation_flux_vr(c,j,icat)*dt, -col_mf%cec_cation_flux2_vr(c,j,icat)*dt, col_mf%background_cec_vr(c,j,icat)*dt
          end do
        end do
      end do
    end if
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Negative cation concentration check
    do fc = 1,num_soilc
      c = filter_soilc(fc)
      g = col_pp%gridcell(c)
      nlevbed = min(col_pp%nlevbed(c), nlevsoi)
      do j = 1,nlevbed
        do icat = 1,ncations
          if (col_ms%cation_vr(c,j,icat) < 0 .and. col_ms%cation_vr(c,j,icat) > -1e-12_r8) then
            ! Numerical accuracy problems can cause runoff to leach away all the cations
            ! Reset cation_vr to zero in this case. Balance Check will ignore everything 
            ! smaller than 1e-12
            col_ms%cation_vr(c,j,icat) = 0._r8
          end if

          if (col_ms%cation_vr(c,j,icat) < 0 .or. col_ms%cec_cation_vr(c,j,icat) < 0) then
            write (100+iam, *) 'Negative cation_vr/cec_cation_vr diagnostics:', ldomain%latc(g), ldomain%lonc(g), c, j, icat, col_ms%cation_vr(c,j,icat), trim(dateTimeString)

            !------------------------------------------------------------------------------
            ! Print out the mass balance diagnostics like above
            write (100+iam, *) 'Post-reaction H+: ', ldomain%latc(g), ldomain%lonc(g), trim(dateTimeString)
            write (100+iam, *) c, j, col_ms%soil_ph(c,j), col_ms%proton_vr(c,j), & 
              mass_to_mol(col_ms%proton_vr(c,j), mass_h, col_ws%h2osoi_vol(c,j))

            write (100+iam, *) 'Post-reaction cation: ', ldomain%latc(g), ldomain%lonc(g), trim(dateTimeString)
            write (100+iam, *) c, j, icat, col_ms%cation_vr(c,j,icat), &
              mass_to_mol(col_ms%cation_vr(c,j,icat), EWParamsInst%cations_mass(icat), &
                          col_ws%h2osoi_vol(c,j)), &
              col_mf%background_flux_vr(c,j,icat)*dt, &
              col_mf%primary_cation_flux_vr(c,j,icat)*dt, &
              col_mf%cec_cation_flux_vr(c,j,icat)*dt, & 
              col_mf%cec_cation_flux2_vr(c,j,icat)*dt, & 
              col_mf%secondary_cation_flux_vr(c,j,icat)*dt, &
              - col_mf%cation_uptake_vr(c,j,icat)*dt, &
              col_mf%cation_infl_vr(c,j,icat)*dt, &
              - col_mf%cation_leached_vr(c,j,icat)*dt, &
              - col_mf%cation_runoff_vr(c,j,icat)*dt

            write (100+iam, *) 'Post-reaction cec H+: ', ldomain%latc(g), ldomain%lonc(g), trim(dateTimeString)
            write (100+iam, *) c, j, col_ms%cec_proton_vr(c,j), & 
              mass_to_meq(col_ms%cec_proton_vr(c,j), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
              col_mf%cect_delta(c,j), &
              mass_to_meq(col_mf%cec_cation_flux_vr(c,j,1)*dt/EWParamsInst%cations_mass(1) & 
                *mass_h*EWParamsInst%cations_valence(1), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
              mass_to_meq(col_mf%cec_cation_flux_vr(c,j,2)*dt/EWParamsInst%cations_mass(2) & 
                *mass_h*EWParamsInst%cations_valence(2), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
              mass_to_meq(col_mf%cec_cation_flux_vr(c,j,3)*dt/EWParamsInst%cations_mass(3) & 
                *mass_h*EWParamsInst%cations_valence(3), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
              mass_to_meq(col_mf%cec_cation_flux_vr(c,j,4)*dt/EWParamsInst%cations_mass(4) & 
                *mass_h*EWParamsInst%cations_valence(4), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
              mass_to_meq(col_mf%cec_cation_flux_vr(c,j,5)*dt/EWParamsInst%cations_mass(5) & 
              *mass_h*EWParamsInst%cations_valence(5), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
              mass_to_meq(col_mf%cec_cation_flux2_vr(c,j,1)*dt/EWParamsInst%cations_mass(1) & 
                *mass_h*EWParamsInst%cations_valence(1), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
              mass_to_meq(col_mf%cec_cation_flux2_vr(c,j,2)*dt/EWParamsInst%cations_mass(2) & 
                *mass_h*EWParamsInst%cations_valence(2), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
              mass_to_meq(col_mf%cec_cation_flux2_vr(c,j,3)*dt/EWParamsInst%cations_mass(3) & 
                *mass_h*EWParamsInst%cations_valence(3), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
              mass_to_meq(col_mf%cec_cation_flux2_vr(c,j,4)*dt/EWParamsInst%cations_mass(4) & 
                *mass_h*EWParamsInst%cations_valence(4), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
              mass_to_meq(col_mf%cec_cation_flux2_vr(c,j,5)*dt/EWParamsInst%cations_mass(5) & 
                *mass_h*EWParamsInst%cations_valence(5), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
              mass_to_meq(col_mf%background_cec_vr(c,j,1)*dt/EWParamsInst%cations_mass(1) & 
                *mass_h*EWParamsInst%cations_valence(1), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
              mass_to_meq(col_mf%background_cec_vr(c,j,2)*dt/EWParamsInst%cations_mass(2) & 
                *mass_h*EWParamsInst%cations_valence(2), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
              mass_to_meq(col_mf%background_cec_vr(c,j,3)*dt/EWParamsInst%cations_mass(3) & 
                *mass_h*EWParamsInst%cations_valence(3), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
              mass_to_meq(col_mf%background_cec_vr(c,j,4)*dt/EWParamsInst%cations_mass(4) & 
                *mass_h*EWParamsInst%cations_valence(4), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
              mass_to_meq(col_mf%background_cec_vr(c,j,5)*dt/EWParamsInst%cations_mass(5) & 
                *mass_h*EWParamsInst%cations_valence(5), 1._r8, mass_h, soilstate_vars%bd_col(c,j))

            write (100+iam, *) 'Post-reaction cec cation: ', ldomain%latc(g), ldomain%lonc(g), trim(dateTimeString)
            write (100+iam, *) c, j, icat, col_ms%cec_cation_vr(c,j,icat), &
              mass_to_meq(col_ms%cec_cation_vr(c,j,icat), EWParamsInst%cations_valence(icat), &
                          EWParamsInst%cations_mass(icat), soilstate_vars%bd_col(c,j)), &
              -col_mf%cec_cation_flux_vr(c,j,icat)*dt, -col_mf%cec_cation_flux2_vr(c,j,icat)*dt, col_mf%background_cec_vr(c,j,icat)*dt
            !------------------------------------------------------------------------------

            call endrun(msg='cation_vr/cec_cation_vr < 0')
          end if
        end do
      end do
    end do

    !end if

    !------------------------------------------------------------------------------
  end subroutine MineralStateDiags

  !-----------------------------------------------------------------------
  subroutine MineralFluxLimit(num_soilc, filter_soilc, col_ms, col_mf, dt, soilstate_vars)
    ! 
    ! !DESCRIPTION:
    ! Scale down reaction flux rates if they cause negative cation balance or
    ! exceed total cation exchange capacity.
    ! 
    ! Order of priority: 
    !  1. Primary mineral must be > 0 after dissolution (done in EnhancedWeatheringMod.F90)
    !     (fully determines r_dissolve and related terms @(2))
    !  2. Secondary mineral must be > 0 after dissolution, and solution cation must be > 0
    !     after secondary mineral precipitation (done in EnhancedWeatheringMod.F90)
    !     (fully determines r_precip and related terms @(1))
    !  3. pH-induced change in total CEC must not cause negative CEC
    !     (individual cation release already restricted in EnhancedWeatheringMod.F90)
    !     (partly determines cect_delta @(3), fully determines cec_cation_flux2_vr @(3)&(4))
    !  4. When cec_cation_flux_vr > 0, it must not cause negative CEC
    !     (fully determines cec_cation_flux_vr @(4) when cec_cation_flux_vr > 0)
    !  5. When cec_cation_flux_vr < 0, it must not cause negative solution cation
    !     (fully determines cec_cation_flux_vr @(4) when cec_cation_flux_vr < 0)
    !  6. The total cations in CEC must not exceed the total CEC minus minimum CEC[H+]
    !     (fully determines cect_delta @(3))
    ! 
    !     ▲ ┌──────┐  (4) ┌─────────┐
    ! (3) │ │ CEC  │ ◀──▶ │Solution │
    !     ▼ └──────┘      └─────────┘
    !                      ▲       ▲
    !                  (1) │       │ (2)
    !                      ▼       │
    !              ┌─────────┐   ┌──────────┐
    !              │ Minsecs │   │ Primary  │
    !              └─────────┘   └──────────┘
    !
    !$acc routine seq
    !
    ! !USES:
    use SharedParamsMod   , only : ParamsShareInst
    !
    ! !ARGUMENTS:
    integer                      , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                      , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(column_mineral_state)   , intent(inout) :: col_ms
    type(column_mineral_flux)    , intent(inout) :: col_mf
    real(r8)                     , intent(in)    :: dt              ! radiation time step (seconds)
    type(soilstate_type)         , intent(in)    :: soilstate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,icat,m,g ! indices
    integer  :: fc        ! lake filter indices
    integer  :: nlevbed
    character(len=20) :: icat_str, c_str, j_str
    real(r8) :: residual_quantity, residual_cect
    real(r8) :: temp_avail_cece(1:num_soilc, 1:nlevsoi, 1:ncations)
    real(r8) :: temp_delta1_cece(1:num_soilc, 1:nlevsoi, 1:ncations)
    real(r8) :: temp_delta2_cece(1:num_soilc, 1:nlevsoi, 1:ncations)
    real(r8) :: temp_avail_ceca(1:num_soilc, 1:nlevsoi)
    real(r8) :: temp_delta1_cation(1:num_soilc, 1:nlevsoi, 1:ncations)
    real(r8) :: temp_delta2_cation(1:num_soilc, 1:nlevsoi, 1:ncations)
    real(r8) :: temp_avail_cation(1:num_soilc, 1:nlevsoi)
    integer  :: err_fc, err_lev, err_icat, err_col
    character(len=256) :: dateTimeString
    logical  :: print_flux_limit(1:num_soilc, 1:nlevsoi, 1:ncations)
    character(len=32) :: subname = 'elm_erw_mineral_flux_limit'  ! subroutine name
    !-----------------------------------------------------------------------

    call get_curr_time_string(dateTimeString)

    ! ensure a tiny bit of cation is left due to numerical accuracy reasons
    residual_cect = 6e-20_r8 ! minimum amount of total CEC to be left in the soil
    residual_quantity = 1e-20_r8 ! minimum amount of individual CEC or pore water cations to be left in the soil

    print_flux_limit(1:num_soilc, 1:nlevsoi, 1:ncations) = .false.

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      nlevbed = min(col_pp%nlevbed(c), nlevsoi)
      do j = 1,nlevbed

        ! ---------------------------------------------------------------------------
        ! Limit the change in total cation exchange capacity due to the total amount
        ! of sites
        ! ---------------------------------------------------------------------------
        if (col_ms%cect_dyn(c,j) <= residual_cect) then
          ! If the dynamic CEC is too small, do not allow any changes to the CEC
          col_mf%cec_delta_limit(c,j) = 0._r8

          col_mf%cect_delta(c,j) = 0._r8
          col_mf%cece_delta(c,j,1:ncations) = 0._r8
          col_mf%cec_cation_flux2_vr(c,j,1:ncations) = 0._r8

          print_flux_limit(fc,j) = .true.
        else if (col_mf%cect_delta(c,j) + col_ms%cect_dyn(c,j) < residual_cect) then
          col_mf%cec_delta_limit(c,j) = - (col_ms%cect_dyn(c,j) - residual_cect) / col_mf%cect_delta(c,j)

          col_mf%cect_delta(c,j) = col_mf%cect_delta(c,j) * col_mf%cec_delta_limit(c,j)
          do icat = 1,ncations
            col_mf%cece_delta(c,j,icat) = col_mf%cece_delta(c,j,icat) * col_mf%cec_delta_limit(c,j)
            col_mf%cec_cation_flux2_vr(c,j,icat) = col_mf%cec_cation_flux2_vr(c,j,icat) * col_mf%cec_delta_limit(c,j)
          end do

          ! ==========================================
          ! Assert total CEC is above minimum allowed
          ! ==========================================
          if (col_mf%cect_delta(c,j) + col_ms%cect_dyn(c,j) < residual_cect) then
            write(c_str, '(I0)') c
            write(j_str, '(I0)') j
            call endrun(msg='Total CEC < minimum allowed (6e-20 g 100g-1) after pH induced change at column '//trim(c_str)//' level '//trim(j_str))
          end if

          print_flux_limit(fc,j) = .true.

        else
          col_mf%cec_delta_limit(c,j) = 1._r8
        end if

        do icat = 1,ncations

          ! ---------------------------------------------------------------------------
          ! If flow direction is from CEC to solution, limit this flux due to the 
          ! availability of individual cations in the CEC phase
          ! ---------------------------------------------------------------------------
          if (col_mf%cec_cation_flux_vr(c,j,icat) > 0._r8) then

            ! CEC-adsorbed cation content after pH-induced change in CEC and background replenishment
            temp_avail_cece(fc,j,icat) = mass_to_meq(col_ms%cec_cation_vr(c,j,icat) - &
              col_mf%cec_cation_flux2_vr(c,j,icat) * dt + col_mf%background_cec_vr(c,j,icat) * dt, &
              EWParamsInst%cations_valence(icat), EWParamsInst%cations_mass(icat), soilstate_vars%bd_col(c,j))

            ! flux from CEC to solution in meq 100g-1 soil unit
            temp_delta1_cation(fc,j,icat) = mass_to_meq(col_mf%cec_cation_flux_vr(c,j,icat) * dt, &
              EWParamsInst%cations_valence(icat), EWParamsInst%cations_mass(icat), soilstate_vars%bd_col(c,j))

            if (temp_avail_cece(fc,j,icat) <= residual_quantity) then
              ! if the available cation sites are too small, do not allow flow to solution
              col_mf%cec_limit_vr(c,j,icat) = 0._r8
              col_mf%cec_cation_flux_vr(c,j,icat) = 0._r8

              print_flux_limit(fc,j) = .true.
            else
              ! If the CEC-adsorbed cation content becomes too small after flow to solution, 
              ! scale down the flux from CEC to solution
              if (temp_avail_cece(fc,j,icat) - temp_delta1_cation(fc,j,icat) <= residual_quantity) then
                col_mf%cec_limit_vr(c,j,icat) = (temp_avail_cece(fc,j,icat) - residual_quantity) / temp_delta1_cation(fc,j,icat)
                col_mf%cec_cation_flux_vr(c,j,icat) = col_mf%cec_cation_flux_vr(c,j,icat) * col_mf%cec_limit_vr(c,j,icat)

                print_flux_limit(fc,j) = .true.
              else
                col_mf%cec_limit_vr(c,j,icat) = 1._r8
              end if

              ! ===========================================================
              ! Assert the individual cation's CEC is above minimum allowed
              ! ===========================================================
              temp_avail_cece(fc,j,icat) = mass_to_meq(col_ms%cec_cation_vr(c,j,icat) - &
                col_mf%cec_cation_flux2_vr(c,j,icat) * dt + col_mf%background_cec_vr(c,j,icat) * dt - &
                col_mf%cec_cation_flux_vr(c,j,icat) * dt, EWParamsInst%cations_valence(icat), &
                EWParamsInst%cations_mass(icat), soilstate_vars%bd_col(c,j))
              if (temp_avail_cece(fc,j,icat) < residual_quantity) then
                write(c_str, '(I0)') c
                write(j_str, '(I0)') j
                write(icat_str, '(I0)') icat
                call endrun(msg='CEC '//trim(icat_str)//' < minimum allowed (1e-20 g 100g-1) after flow to solution at column '//trim(c_str)//' level '//trim(j_str))
              end if
            end if

          ! ---------------------------------------------------------------------------
          ! If flow direction is from solution to CEC, limit this flux due to the 
          ! availability of individual cations in the solution (after secondary mineral
          ! precipitation or replenishment and background flux replenishment)
          ! ---------------------------------------------------------------------------
          else

            ! Assert the primary mineral cation flux to solution is always >= 0
            if (col_mf%primary_cation_flux_vr(c,j,icat) < 0._r8) then
              write(c_str, '(I0)') c
              write(j_str, '(I0)') j
              write(icat_str, '(I0)') icat
              call endrun(msg='Erroneous negative primary dissolution flux for cation '//trim(icat_str)//' at column '//trim(c_str)//' level '//trim(j_str))
            end if

            ! Assert the background flux to solution is always >= 0
            if (col_mf%background_flux_vr(c,j,icat) < 0._r8) then
              write(c_str, '(I0)') c
              write(j_str, '(I0)') j
              write(icat_str, '(I0)') icat
              call endrun(msg='Erroneous negative background flux for cation '//trim(icat_str)//' at column '//trim(c_str)//' level '//trim(j_str))
            end if

            ! Assert the flow from CEC to solution due to change in total CEC
            ! is always >= 0, because it only happens when the total CEC declines due to pH change
            if (col_mf%cec_cation_flux2_vr(c,j,icat) < 0._r8) then
              write(c_str, '(I0)') c
              write(j_str, '(I0)') j
              write(icat_str, '(I0)') icat
              call endrun(msg='Erroneous pH-induced increase in CEC-phase cation '//trim(icat_str)//' at column '//trim(c_str)//' level '//trim(j_str))
            end if

            ! Assert the restriction on r_precip worked
            if (col_ms%cation_vr(c,j,icat) + col_mf%secondary_cation_flux_vr(c,j,icat)*dt < 0._r8) then
              write(c_str, '(I0)') c
              write(j_str, '(I0)') j
              write(icat_str, '(I0)') icat
              call endrun(msg='Too large secondary mineral precipitation for cation '//trim(icat_str)//' at column '//trim(c_str)//' level '//trim(j_str))
            end if

            ! delta1 should always be positive when the above assertions are true
            temp_delta1_cation(fc,j,icat) = col_mf%primary_cation_flux_vr(c,j,icat)*dt + &
              col_mf%background_flux_vr(c,j,icat)*dt + col_mf%cec_cation_flux2_vr(c,j,icat)*dt + &
              col_mf%secondary_cation_flux_vr(c,j,icat)*dt

            ! delta2 is the absolute value out of the solution
            temp_delta2_cation(fc,j,icat) = (col_mf%cation_uptake_vr(c,j,icat) - col_mf%cec_cation_flux_vr(c,j,icat))*dt

            if (temp_delta1_cation(fc,j,icat) <= residual_quantity) then
              col_mf%flux_limit_vr(c,j,icat) = 0._r8
              col_mf%cation_uptake_vr(c,j,icat) = 0._r8
              col_mf%cec_cation_flux_vr(c,j,icat) = 0._r8

              print_flux_limit(fc,j) = .true.

            else if (temp_delta1_cation(fc,j,icat) - temp_delta2_cation(fc,j,icat) < residual_quantity) then
              col_mf%flux_limit_vr(c,j,icat) = (temp_delta1_cation(fc,j,icat) - residual_quantity) / temp_delta2_cation(fc,j,icat)
              col_mf%cation_uptake_vr(c,j,icat) = col_mf%cation_uptake_vr(c,j,icat) * col_mf%flux_limit_vr(c,j,icat)
              col_mf%cec_cation_flux_vr(c,j,icat) = col_mf%cec_cation_flux_vr(c,j,icat) * col_mf%flux_limit_vr(c,j,icat)

              !===========================================================
              ! Assert the solution cation is above minimum allowed
              !===========================================================
              temp_avail_cation(fc,j,icat) = col_ms%cation_vr(c,j,icat) + col_mf%background_flux_vr(c,j,icat)*dt + &
                col_mf%primary_cation_flux_vr(c,j,icat)*dt + col_mf%cec_cation_flux_vr(c,j,icat)*dt + &
                col_mf%cec_cation_flux2_vr(c,j,icat)*dt + col_mf%secondary_cation_flux_vr*dt - col_mf%cation_uptake_vr(c,j,icat)*dt
              if (temp_avail_cation(fc,j,icat) < residual_quantity) then
                write(c_str, '(I0)') c
                write(j_str, '(I0)') j
                write(icat_str, '(I0)') icat
                call endrun(msg='Solution cation '//trim(icat_str)//' < minimum allowed (1e-20 g m-3 soil) after all reaction fluxes at column '//trim(c_str)//' level '//trim(j_str))
              end if

              print_flux_limit(fc,j) = .true.

            else
              col_mf%flux_limit_vr(c,j,icat) = 1._r8

            end if
          end if

        end do

        ! ---------------------------------------------------------------------------
        ! If there is not enough H+ sites (the difference between total CEC after 
        ! pH-induced change, and the cation-occupied sites after applying various changes),
        ! then increase the pH-induced change of total CEC, i.e. adjust cect_delta,
        ! so that CEC[H+] > residual_quantity. 
        ! Note the change in CEC[X+] for the cations do not need to be adjusted;
        ! DO NOT adjust background_cec_vr to insure control-run v.s. treatment-run. consistency.
        ! ---------------------------------------------------------------------------
        do icat = 1,ncations

          ! Assert the background CEC flux is always >= 0
          if (col_mf%background_cec_vr(c,j,icat) < 0._r8) then
            write(c_str, '(I0)') c
            write(j_str, '(I0)') j
            write(icat_str, '(I0)') icat
            call endrun(msg='Erroneous negative CEC background flux for cation '//trim(icat_str)//' at column '//trim(c_str)//' level '//trim(j_str))
          end if

          ! cation-occupied sites after various cec change (meq 100g-1 soil unit)
          temp_avail_cece(fc,j,icat) = mass_to_meq( &
            col_ms%cec_cation_vr(c,j,icat) + (- col_mf%cec_cation_flux_vr(c,j,icat) &
            - col_mf%cec_cation_flux2_vr(c,j,icat) + col_mf%background_cec_vr(c,j,icat))*dt, &
            EWParamsInst%cations_valence(icat), EWParamsInst%cations_mass(icat), &
            soilstate_vars%bd_col(c,j))
        end do

        ! calculate the expected available H+ sites, after total CEC changes
        temp_avail_ceca(fc,j) = col_ms%cect_dyn(c,j) + col_mf%cect_delta(c,j) - sum(temp_avail_cec(fc,j,1:ncations))
        if (temp_avail_ceca(fc,j) < residual_quantity) then
          col_mf%cect_delta_add(c,j) = residual_quantity - (col_ms%cect_dyn(c,j) + &
            col_mf%cect_delta(c,j) - sum(temp_avail_cec(fc,j,1:ncations)))
          col_mf%cect_delta(c,j) = col_mf%cect_delta(c,j) + col_mf%cect_delta_add(c,j)

          print_flux_limit = .true.
        else
          col_mf%cect_delta_add(c,j) = 0._r8
        end if
        ! re-calculate the available H+ sites, assert it is above minimum allowed
        temp_avail_ceca(fc,j) = col_ms%cect_dyn(c,j) + col_mf%cect_delta(c,j) - sum(temp_avail_cec(fc,j,1:ncations))
        if (temp_avail_ceca(fc,j,icat) < residual_quantity) then
          write(c_str, '(I0)') c
          write(j_str, '(I0)') j
          write(icat_str, '(I0)') icat
          call endrun(msg='CEC H+ '//trim(icat_str)//' < minimum allowed (1e-20 meq 100g-1) after all reaction fluxes at column '//trim(c_str)//' level '//trim(j_str))
        end if

      end do
    end do

    ! -------------------------------------------------------------------------------------------
    ! Print out flux limit factor on the fly if verbose mode
    if (use_erw_verbose > 0) then
      do fc = 1,num_soilc
        c = filter_soilc(fc)
        g = col_pp%gridcell(c)
        do j = 1,nlevbed
          if (print_flux_limit(fc,j)) then

            write (100+iam, *) 'Flux limit diagnostics: ', ldomain%latc(g), ldomain%lonc(g), j, trim(dateTimeString)
            call shr_sys_flush(100+iam)

            if (col_mf%cec_delta_limit(c,j,icat) < 1._r8) then
              write (100+iam, *) '   negative total CEC ', c, j, icat, col_mf%cec_delta_limit(c,j,icat), col_mf%cect_delta(c,j), col_mf%cect_dyn(c,j), col_mf%cece_delta(c,j,1:ncations)
            end if

            do icat = 1,ncations
              if (col_mf%cec_limit_vr(c,j,icat) < 1._r8) then
                write (100+iam, *) '   negative CEC cation ', c, j, icat, col_mf%cec_limit_vr(c,j,icat), col_ms%cec_cation_vr(c,j,icat), col_mf%background_cec_vr(c,j,icat)*dt, - col_mf%cec_cation_flux_vr(c,j,icat)*dt, - col_mf%cec_cation_flux2_vr(c,j,icat)*dt
              end if
            end do

            do icat = 1,ncations
              if (col_mf%flux_limit_vr(c,j,icat) < 1._r8) then
                write (100+iam, *) '   negative solution cation ', c, j, icat, col_mf%flux_limit_vr(c,j,icat), col_ms%cation_vr(c,j,icat), col_mf%primary_cation_flux_vr(c,j,icat)*dt, col_mf%background_flux_vr(c,j,icat)*dt, - col_mf%cation_uptake_vr(c,j,icat)*dt, col_mf%secondary_cation_flux_vr(c,j,icat)*dt, col_mf%cec_cation_flux_vr(c,j,icat)*dt, col_mf%cec_cation_flux2_vr(c,j,icat)*dt
              end if
            end do

            if (col_mf%cect_delta_add(c,j) > 0._r8) then
              write (100+iam, *) '   negative CEC H+ ', c,j, col_mf%cect_delta_add(c,j), temp_avail_ceca(fc,j), col_mf%cect_delta(c,j), col_mf%cect_dyn(c,j), temp_avail_cec(fc,j,1:ncations)
            end if

          end if
        end do
      end do
    end if

  end subroutine MineralFluxLimit

end module MineralStateUpdateMod