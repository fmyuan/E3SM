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
  use ColumnDataType          , only : col_ws
  use ColumnDataType          , only : col_ms, col_mf, col_pp
  use ColumnDataType          , only : column_mineral_state, column_mineral_flux, column_water_flux
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
  subroutine MineralStateUpdate1(num_soilc, filter_soilc, col_ms, col_mf, dt, soilstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the mineral state variables that are not
    ! affected by vertical or horizontal soil water movement, update pH-dependent CEC.
    !
    !$acc routine seq
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
    integer  :: c,p,j,k,icat,m,g ! indices
    integer  :: fp,fc         ! lake filter indices
    integer  :: nlevbed
    !-----------------------------------------------------------------------

    ! Update mineral state
    do fc = 1,num_soilc
      c = filter_soilc(fc)
      nlevbed = min(col_pp%nlevbed(c), nlevsoi)

      do j = 1,nlevbed
        ! soil H+ concentration
        ! only determined by CEC equilibrium
        col_ms%proton_vr(c,j) = mol_to_mass(10**(-col_ms%soil_ph(c,j)), mass_h, & 
                                            col_ws%h2osoi_vol(c,j))

        ! soil cation concentration - not updated here
        ! must be preserved before calling the vertical solute movement solver

        ! pH-dependent CEC
        col_ms%cect_dyn(c,j) = col_mf%cect_delta(c,j) + col_ms%cect_dyn(c,j)

        ! CEC cations - only depends on flux limit

        do icat = 1,ncations

          !!write (iam+100, *) '++++++++++++++++++++++++++'
          !!write (iam+100, *) c, j, icat, col_ms%cec_cation_vr(c,j,icat), col_mf%cec_cation_flux_vr(c,j,icat)*dt, col_mf%cec_cation_flux2_vr(c,j,icat)*dt, col_mf%background_cec_vr(c,j,icat)*dt
          !!write (iam+100, *) '++++++++++++++++++++++++++'

          col_ms%cec_cation_vr(c,j,icat) = col_ms%cec_cation_vr(c,j,icat) - &
                  (col_mf%cec_cation_flux_vr(c,j,icat) + col_mf%cec_cation_flux2_vr(c,j,icat) &
                   - col_mf%background_cec_vr(c,j,icat))*dt

          !!write (iam+100, *) '++++++++++++++++++++++++++'
          !!write (iam+100, *) c, j, icat, col_ms%cec_cation_vr(c,j,icat),  col_mf%cec_cation_flux_vr(c,j,icat)*dt, col_mf%cec_cation_flux2_vr(c,j,icat)*dt, col_mf%background_cec_vr(c,j,icat)*dt
          !!write (iam+100, *) '++++++++++++++++++++++++++'

        end do

        ! CEC H+ with pH dependent changes
        ! the Equilibria subroutine cannot distinguish the effects of CO2 and cation exchange
        ! instead, use charge balance on the mineral surface to get the change in adsorped H+
        ! note the newly exposed surface due to dynamic CEC
        do icat = 1,ncations
          col_ms%cec_proton_vr(c,j) = col_ms%cec_proton_vr(c,j) + &
            (col_mf%cec_cation_flux_vr(c,j,icat) + col_mf%cec_cation_flux2_vr(c,j,icat) &
             - col_mf%background_cec_vr(c,j,icat))*dt & 
            / EWParamsInst%cations_mass(icat) * mass_h * EWParamsInst%cations_valence(icat)
        end do
        col_ms%cec_proton_vr(c,j) = col_ms%cec_proton_vr(c,j) + &
          meq_to_mass(col_mf%cect_delta(c,j), 1._r8, mass_h, soilstate_vars%bd_col(c,j))

        !!DEBUG
        !!write (iulog, *) j, 'cect_dyn', col_ms%cect_dyn(c,j), col_mf%cect_delta(c,j)
        !!do icat = 1,ncations
        !!  write (iulog, *) j, 'cec_cation_vr', &
        !!    mass_to_meq(col_ms%cec_cation_vr(c,j,icat), &
        !!    EWParamsInst%cations_valence(icat), EWParamsInst%cations_mass(icat), &
        !!    soilstate_vars%bd_col(c,j)), mass_to_meq(col_mf%cec_cation_flux2_vr(c,j,icat), &
        !!    EWParamsInst%cations_valence(icat), EWParamsInst%cations_mass(icat), &
        !!    soilstate_vars%bd_col(c,j)), icat
        !!end do
        !!write (iulog, *) j, 'cec_proton', mass_to_meq(col_ms%cec_proton_vr(c,j), &
        !!    1._r8, mass_h, soilstate_vars%bd_col(c,j))
        !!write (iulog, *) j, 'col_mf%cec_delta_limit(c,j)', col_mf%cec_delta_limit(c,j)

        ! primary mineral
        do m = 1,nminerals
          col_ms%primary_mineral_vr(c,j,m) = col_ms%primary_mineral_vr(c,j,m) + &
            col_mf%primary_added_vr(c,j,m)*dt - col_mf%primary_dissolve_vr(c,j,m)*dt

          !!if (m == 6) then
          !!  write (iulog, *) c, j, m, col_ms%primary_mineral_vr(c,j,m), col_mf%primary_added_vr(c,j,m), col_mf%primary_dissolve_vr(c,j,m), dt
          !!end if

          ! non-SiO2 minerals
          col_ms%primary_residue_vr(c,j,m) = col_ms%primary_residue_vr(c,j,m) + &
            col_mf%primary_residue_flux_vr(c,j,m)*dt
        end do

        ! silicate
        col_ms%silica_vr(c,j) = col_ms%silica_vr(c,j) + col_mf%primary_silica_flux_vr(c,j) * dt - col_mf%secondary_silica_flux_vr(c,j) * dt

        ! secondary mineral
        do m = 1,nminsecs
            col_ms%secondary_mineral_vr(c,j,m) = col_ms%secondary_mineral_vr(c,j,m) + col_mf%secondary_mineral_flux_vr(c,j,m) * dt
        end do

        ! ignore the effect on soil water content for now

        ! passivation layer thickness
        col_ms%passivation_thickness(c,j) = col_ms%passivation_thickness(c,j) + &
          col_mf%passivation_rate(c,j) * dt
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
              col_mf%cec_cation_flux2_vr(c,j,icat) - &
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
        col_mf%r_sequestration(c) = col_mf%bicarbonate_drainage(c) / mass_hco3 + &
          col_mf%carbonate_drainage(c) / mass_co3 * 2._r8
        do j = 1,nlevbed
          ! add the subsurface drainage
          col_mf%r_sequestration(c) = col_mf%r_sequestration(c) + &
            col_mf%bicarbonate_leached_vr(c,j) * col_pp%dz(c,j) / mass_hco3 + &
            col_mf%carbonate_leached_vr(c,j) * col_pp%dz(c,j) / mass_co3 * 2._r8

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
    integer  :: c,j,icat,g    ! indices
    integer  :: fc            ! lake filter indices
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
          step_delta = - col_mf%secondary_cation_flux_vr(c,j,icat) - &
                       col_mf%cation_uptake_vr(c,j,icat) + &
                       col_mf%cation_infl_vr(c,j,icat) - & 
                       col_mf%cation_oufl_vr(c,j,icat) - &
                       col_mf%cation_leached_vr(c,j,icat) - &
                       col_mf%cation_runoff_vr(c,j,icat)

          col_mf%tempavg_tot_delta(c,j,icat) = col_mf%tempavg_tot_delta(c,j,icat) + &
            step_delta * fracday / dayspyr_mod

          ! also calibrate a replenishment term to the cation exchange phase
          ! because it seems to be lost pretty severly
          ! note: cec_cation_flux_vr is defined negative for adsorption into soil
          !       this term is defined positive for adsorption into soil
          col_mf%tempavg_cec_delta(c,j,icat) = col_mf%tempavg_cec_delta(c,j,icat) - &
            (col_mf%cec_cation_flux_vr(c,j,icat) + col_mf%cec_cation_flux2_vr(c,j,icat) &
            ) * fracday / dayspyr_mod
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
              - col_mf%secondary_cation_flux_vr(c,j,icat)*dt, &
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
              - col_mf%secondary_cation_flux_vr(c,j,icat)*dt, &
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
    ! exceed total cation exchange capacity
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
    real(r8) :: residual_factor
    real(r8) :: temp_avail_cece(1:num_soilc, 1:nlevsoi, 1:ncations)
    real(r8) :: temp_delta1_cece(1:num_soilc, 1:nlevsoi, 1:ncations)
    real(r8) :: temp_delta2_cece(1:num_soilc, 1:nlevsoi, 1:ncations)
    real(r8) :: temp_delta_ceca(1:num_soilc, 1:nlevsoi)
    real(r8) :: temp_avail_ceca(1:num_soilc, 1:nlevsoi)
    real(r8) :: temp_delta1_cation(1:num_soilc, 1:nlevsoi, 1:ncations)
    real(r8) :: temp_delta2_cation(1:num_soilc, 1:nlevsoi, 1:ncations)
    real(r8) :: min_flux_limit(1:num_soilc, 1:nlevsoi)
    real(r8) :: small_fix_cec
    logical  :: err_found
    integer  :: err_fc, err_lev, err_icat, err_col
    character(len=256) :: dateTimeString
    logical  :: print_flux_limit
    character(len=32) :: subname = 'elm_erw_mineral_flux_limit'  ! subroutine name
    !-----------------------------------------------------------------------

    call get_curr_time_string(dateTimeString)

    ! ensure a tiny bit of cation is left due to numerical accuracy reasons
    residual_factor = 0.99_r8

    min_flux_limit(1:num_soilc, 1:nlevsoi) = 1._r8
    err_found = .false.

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      nlevbed = min(col_pp%nlevbed(c), nlevsoi)
      do j = 1,nlevbed

        ! Limit the change in total cation exchange capacity due to the total amount
        ! of sites
        if (col_mf%cect_delta(c,j) + col_ms%cect_dyn(c,j) < 0._r8) then
          col_mf%cec_delta_limit(c,j) = - col_ms%cect_dyn(c,j) / &
                                          col_mf%cect_delta(c,j) * residual_factor
          col_mf%cect_delta(c,j) = col_mf%cect_delta(c,j) * col_mf%cec_delta_limit(c,j)
          do icat = 1,ncations
            col_mf%cece_delta(c,j,icat) = col_mf%cece_delta(c,j,icat) * col_mf%cec_delta_limit(c,j)
            col_mf%cec_cation_flux2_vr(c,j,icat) = col_mf%cec_cation_flux2_vr(c,j,icat) * &
                                                   col_mf%cec_delta_limit(c,j)
          end do
        else
          col_mf%cec_delta_limit(c,j) = 1._r8
        end if

        !!write (iam+100, *) '---------------------------------------'
        !!write (iam+100, *) 'z', col_mf%cec_delta_limit(c,j) 
        !!write (iam+100, *) '---------------------------------------'

        ! Limit the cec cation flux due to the availability of individual cations in the
        ! cation exchange phase
        do icat = 1,ncations

          ! cation-occupied sites after total cec change
          temp_avail_cece(fc,j,icat) = col_ms%cec_cation_vr(c,j,icat) - &
            col_mf%cec_cation_flux2_vr(c,j,icat) * dt

          !!write (iam+100, *) '-------------------------------------------'
          !!write (iam+100, *) 'a0', c,j,icat, col_ms%cec_cation_vr(c,j,icat),-col_mf%cec_cation_flux2_vr(c,j,icat)*dt, col_mf%cece_delta(c,j,icat)
          !!write (iam+100, *) '-------------------------------------------'
          !!write (iam+100, *) 'a', c,j,icat, temp_avail_cece(fc,j,icat), col_mf%background_cec_vr(c,j,icat)*dt, -col_mf%cec_cation_flux_vr(c,j,icat)*dt
          !!write (iam+100, *) '-------------------------------------------'

          ! background flux will always be > 0
          temp_delta1_cece(fc,j,icat) = col_mf%background_cec_vr(c,j,icat)*dt
          ! cec_cation_flux_vr > 0 := flow from CEC to solution
          temp_delta2_cece(fc,j,icat) = - col_mf%cec_cation_flux_vr(c,j,icat)*dt

          if ((temp_avail_cece(fc,j,icat) + temp_delta1_cece(fc,j,icat) + &
               temp_delta2_cece(fc,j,icat)) < 0._r8) then

            if ((temp_avail_cece(fc,j,icat) + temp_delta1_cece(fc,j,icat)) < 1e-50_r8) then
              ! Do not allow any more flow from CEC to solution to prevent numerical
              ! errors in the calculated cec_limit factor. 
              ! Instead, add a small background flux to increase the cec.
              col_mf%cec_limit_vr(c,j,icat) = 0._r8

              col_mf%background_cec_vr(c,j,icat) = col_mf%background_cec_vr(c,j,icat) + &
                  meq_to_mass(2.0e-4 * col_ms%cect_dyn(c,j), EWParamsInst%cations_valence(icat), &
                              EWParamsInst%cations_mass(icat), soilstate_vars%bd_col(c,j))
            else
              ! Otherwise, calculate the flux limit normally.
              col_mf%cec_limit_vr(c,j,icat) = - (temp_delta1_cece(fc,j,icat) + & 
                                                 temp_avail_cece(fc,j,icat)) / &
                  temp_delta2_cece(fc,j,icat) * residual_factor
            end if
            ! Apply the flux limite factor on cec_cation_flux_vr
            col_mf%cec_cation_flux_vr(c,j,icat) = col_mf%cec_cation_flux_vr(c,j,icat) * & 
                                                  col_mf%cec_limit_vr(c,j,icat)
            min_flux_limit(fc,j) = min(min_flux_limit(fc,j), col_mf%cec_limit_vr(c,j,icat))

          else
            col_mf%cec_limit_vr(c,j,icat) = 1._r8
          end if

          !!write (iam+100, *) '-------------------------------------------'
          !!write (iam+100, *) 'a', c,j,icat, col_mf%background_cec_vr(c,j,icat)*dt, -col_mf%cec_cation_flux_vr(c,j,icat)*dt, -col_mf%cec_cation_flux2_vr(c,j,icat)*dt, col_mf%cece_delta(c,j,icat)
          !!write (iam+100, *) '-------------------------------------------'

        end do


        ! Limit the cec cation flux due to not having enough remaining (H+) sites

        ! - calculate the expected available H+ sites, after total CEC changes
        temp_avail_ceca(fc,j) = mass_to_meq(col_ms%cec_proton_vr(c,j), 1._r8, mass_h, &
                                            soilstate_vars%bd_col(c,j))
        if (col_mf%cect_delta(c,j) > 0._r8) then
          temp_avail_ceca(fc,j) = temp_avail_ceca(fc,j) + col_mf%cect_delta(c,j)
        else
          temp_avail_ceca(fc,j) = temp_avail_ceca(fc,j) + col_mf%cect_delta(c,j) - &
                            sum(col_mf%cece_delta(c,j,1:ncations))
        end if

        ! - calculate the expected H+ sites to be displaced by cations
        !   note col_mf%cec_cation_flux2_vr(c,j,icat) does not displace anything, because
        !   they come from the lost cation-occupied sites, not H+ occupied sites.
        temp_delta_ceca(fc,j) = 0._r8
        do icat = 1,ncations
          temp_delta_ceca(fc,j) = temp_delta_ceca(fc,j) + mass_to_meq( &
            - col_mf%background_cec_vr(c,j,icat) + col_mf%cec_cation_flux_vr(c,j,icat), &
            EWParamsInst%cations_valence(icat), &
            EWParamsInst%cations_mass(icat), soilstate_vars%bd_col(c,j)) * dt
        end do

        ! - compare
        if (temp_avail_ceca(fc,j) <= 0._r8) then

          col_mf%proton_limit_vr(c,j) = 0._r8
          do icat = 1,ncations
            col_mf%background_cec_vr(c,j,icat) = col_mf%background_cec_vr(c,j,icat) * &
                                                 col_mf%proton_limit_vr(c,j)
            col_mf%cec_cation_flux_vr(c,j,icat) = col_mf%cec_cation_flux_vr(c,j,icat) * &
              col_mf%proton_limit_vr(c,j)
          end do

        else if ((temp_avail_ceca(fc,j) + temp_delta_ceca(fc,j)) < 0._r8) then
          col_mf%proton_limit_vr(c,j) = - temp_avail_ceca(fc,j) / temp_delta_ceca(fc,j) &
                                        * residual_factor
          do icat = 1,ncations
            col_mf%background_cec_vr(c,j,icat) = col_mf%background_cec_vr(c,j,icat) * &
                                                 col_mf%proton_limit_vr(c,j)
            col_mf%cec_cation_flux_vr(c,j,icat) = col_mf%cec_cation_flux_vr(c,j,icat) * &
              col_mf%proton_limit_vr(c,j)

            !!write (iam+100, *) '-------------------------------------------'
            !!write (iam+100, *) 'b', c,j,icat, temp_avail_ceca(fc,j), col_mf%cect_delta(c,j), sum(col_mf%cece_delta(c,j,1:ncations)), col_mf%background_cec_vr(c,j,icat)*dt, col_mf%cec_cation_flux_vr(c,j,icat)*dt, col_mf%cec_cation_flux2_vr(c,j,icat)*dt
            !!write (iam+100, *) '-------------------------------------------'

          end do
          min_flux_limit(fc,j) = min(min_flux_limit(fc,j), col_mf%proton_limit_vr(c,j))
        else
          col_mf%proton_limit_vr(c,j) = 1._r8
        end if


        ! Limit due to soil solution cation concentration, after the previous two limits
        ! have been applieds
        do icat = 1,ncations
          ! delta1 is always positive
          temp_delta1_cation(fc,j,icat) = col_mf%primary_cation_flux_vr(c,j,icat)*dt + &
                                          col_mf%background_flux_vr(c,j,icat)*dt + &
                                          col_mf%cec_cation_flux2_vr(c,j,icat)*dt
          ! delta2 may be positive or negative
          temp_delta2_cation(fc,j,icat) = - col_mf%cation_uptake_vr(c,j,icat)*dt - &
                                          col_mf%secondary_cation_flux_vr(c,j,icat)*dt + & 
                                          col_mf%cec_cation_flux_vr(c,j,icat)*dt

          if ((col_ms%cation_vr(c,j,icat) + temp_delta1_cation(fc,j,icat) + & 
              temp_delta2_cation(fc,j,icat)) < 0._r8) then
            ! ensure a tiny bit of cation is left due to numerical accuracy reasons
            col_mf%flux_limit_vr(c,j,icat) = - (temp_delta1_cation(fc,j,icat) + & 
              col_ms%cation_vr(c,j,icat)) / temp_delta2_cation(fc,j,icat) * residual_factor
            col_mf%cation_uptake_vr(c,j,icat) = col_mf%cation_uptake_vr(c,j,icat) * & 
                  col_mf%flux_limit_vr(c,j,icat)
            col_mf%secondary_cation_flux_vr(c,j,icat) = &
                  col_mf%secondary_cation_flux_vr(c,j,icat) * col_mf%flux_limit_vr(c,j,icat)
            col_mf%cec_cation_flux_vr(c,j,icat) = col_mf%cec_cation_flux_vr(c,j,icat) * & 
                  col_mf%flux_limit_vr(c,j,icat)

            min_flux_limit(fc,j) = min(min_flux_limit(fc,j), col_mf%flux_limit_vr(c,j,icat))

          else
            col_mf%flux_limit_vr(c,j,icat) = 1._r8
          end if

          !!write (iam+100, *) '-------------------------------------------'
          !!write (iam+100, *) 'c1', c,j,icat, col_mf%primary_cation_flux_vr(c,j,icat)*dt, col_mf%background_flux_vr(c,j,icat)*dt, col_mf%cec_cation_flux2_vr(c,j,icat)*dt
          !!write (iam+100, *) 'c2', c,j,icat,-col_mf%cation_uptake_vr(c,j,icat)*dt, col_mf%secondary_cation_flux_vr(c,j,icat)*dt, col_mf%cec_cation_flux_vr(c,j,icat)*dt
          !!write (iam+100, *) 'c3', c,j,icat, col_mf%flux_limit_vr(c,j,icat), temp_delta1_cation(fc,j,icat), temp_delta2_cation(fc,j,icat)
          !!write (iam+100, *) '-------------------------------------------'

        end do

      end do
    end do

    ! -------------------------------------------------------------------------------------------
    ! Print out flux limit factor on the fly if verbose mode
    if (use_erw_verbose > 0) then
      do fc = 1,num_soilc
        c = filter_soilc(fc)
        g = col_pp%gridcell(c)
        do j = 1,nlevbed
          if (min_flux_limit(fc,j) < 1._r8) then
            write (100+iam, *) 'Flux limit diagnostics: ', ldomain%latc(g), ldomain%lonc(g), j, trim(dateTimeString)
            call shr_sys_flush(100+iam)

            do icat = 1,ncations
              if (col_mf%cec_limit_vr(c,j,icat) < 1._r8) then
                write (100+iam, *) '   negative CEC cation ', c, j, icat, col_mf%cec_limit_vr(c,j,icat), col_ms%cec_cation_vr(c,j,icat), col_mf%background_cec_vr(c,j,icat)*dt, - col_mf%cec_cation_flux_vr(c,j,icat)*dt, - col_mf%cec_cation_flux2_vr(c,j,icat)*dt
              end if
            end do
            if (col_mf%proton_limit_vr(c,j) < 1._r8) then
              write (100+iam, *) '   negative CEC H+ ', c,j, col_mf%proton_limit_vr(c,j), col_ms%cec_proton_vr(c,j), temp_avail_ceca(fc,j), temp_delta_ceca(fc,j)
            end if
            do icat = 1,ncations
              if (col_mf%flux_limit_vr(c,j,icat) < 1._r8) then
                write (100+iam, *) '   negative solution cation ', c, j, icat, col_mf%flux_limit_vr(c,j,icat), col_ms%cation_vr(c,j,icat), col_mf%primary_cation_flux_vr(c,j,icat)*dt, col_mf%background_flux_vr(c,j,icat)*dt, col_mf%cation_uptake_vr(c,j,icat)*dt, - col_mf%secondary_cation_flux_vr(c,j,icat)*dt, col_mf%cec_cation_flux_vr(c,j,icat)*dt, col_mf%cec_cation_flux2_vr(c,j,icat)*dt
              end if
            end do
          end if
        end do
      end do
    end if

  end subroutine MineralFluxLimit

end module MineralStateUpdateMod