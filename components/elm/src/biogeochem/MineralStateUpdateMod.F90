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
  use elm_varctl              , only : use_erw_verbose
  use shr_sys_mod             , only : shr_sys_flush
  use spmdMod                 , only : masterproc
  use abortutils              , only : endrun
  use shr_log_mod             , only : errMsg => shr_log_errMsg
  use ewutils                 , only : mass_to_mol, mass_to_meq, mol_to_mass
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
  subroutine MineralStateUpdate1(num_soilc, filter_soilc, col_ms, col_mf, dt)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the dissolution/precipitation dynamics and 
    ! cation exchange affected mineral state variables
    !
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

      do j = 1,nlevbed
        ! soil H+ concentration
        ! only determined by CEC equilibrium
        col_ms%proton_vr(c,j) = mol_to_mass(10**(-col_ms%soil_ph(c,j)), mass_h, & 
                                            col_ws%h2osoi_vol(c,j))
        
        ! soil cation concentration - not updated here
        ! must be preserved before calling the vertical solute movement solver

        ! CEC cations - only depends on flux limit
        do icat = 1,ncations
          col_ms%cec_cation_vr(c,j,icat) = col_ms%cec_cation_vr(c,j,icat) - &
                  (col_mf%cec_cation_flux_vr(c,j,icat) - col_mf%background_cec_vr(c,j,icat))*dt
        end do

        ! CEC H+
        ! the Equilibria subroutine cannot distinguish the effects of CO2 and cation exchange
        ! instead, use charge balance on the mineral surface to get the change in adsorped H+
        do icat = 1,ncations
          col_ms%cec_proton_vr(c,j) = col_ms%cec_proton_vr(c,j) + &
            (col_mf%cec_cation_flux_vr(c,j,icat) - col_mf%background_cec_vr(c,j,icat))*dt & 
            / EWParamsInst%cations_mass(icat) * mass_h * EWParamsInst%cations_valence(icat)
        end do

        ! primary mineral
        do m = 1,nminerals
          col_ms%primary_mineral_vr(c,j,m) = col_ms%primary_mineral_vr(c,j,m) + &
            col_mf%primary_added_vr(c,j,m)*dt - col_mf%primary_dissolve_vr(c,j,m)*dt
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

        ! TODO: ignore the effect on soil water for now
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
          col_ms%cation_vr(c, j, icat) = col_ms%cation_vr(c, j, icat) + & 
            ( col_mf%background_flux_vr(c,j,icat) + & 
              col_mf%primary_cation_flux_vr(c,j,icat) + & 
              col_mf%cec_cation_flux_vr(c,j,icat) - &
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
    ! On the radiation time step, update the prognostic mineral state variables
    ! related to leaching and carbon sequestration rate
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

    ! Update mineral state
    do fc = 1,num_soilc
      c = filter_soilc(fc)
      nlevbed = min(col_pp%nlevbed(c), nlevsoi)

      do j = 1,nlevbed
        do icat = 1,ncations
          col_ms%cation_vr(c,j,icat) = col_ms%cation_vr(c,j,icat) - col_mf%cation_leached_vr(c,j,icat)*dt - col_mf%cation_runoff_vr(c,j,icat)*dt
        end do
      end do

      ! Calculate the total CO2 sequestration rate in mol m-2 s-1
      rmethod = 1
      if (rmethod == 1) then

        col_mf%r_sequestration(c) = 0._r8
        do j = 1,nlevbed
          ! precipitated by calcite: 1 mol CO2 per mol Ca2+
          col_mf%r_sequestration(c) = col_mf%r_sequestration(c) + & 
            col_mf%secondary_cation_flux_vr(c,j,1) / EWParamsInst%cations_mass(1) * col_pp%dz(c,j)
          ! transported to ocean: 2x for 2+ cations, 1x for 1+ cations, multiply by
          ! ocean efficiency (0.86)
          ! - col_mf%background_flux_vr(c,j,icat)
          do icat = 1,ncations
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
            col_mf%cec_cation_flux_vr(c,j,icat) * fracday / dayspyr_mod
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

    !------------------------------------------------------------------------------
    ! Write mass balance diagnostics only if verbose level is set to high
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
        do icat = 1,ncations
          do j = 1,nlevbed
            write (100+iam, *) c, j, icat, col_ms%cation_vr(c,j,icat), &
              mass_to_mol(col_ms%cation_vr(c,j,icat), EWParamsInst%cations_mass(icat), &
                          col_ws%h2osoi_vol(c,j)), &
              col_mf%background_flux_vr(c,j,icat)*dt, &
              col_mf%primary_cation_flux_vr(c,j,icat)*dt, &
              col_mf%cec_cation_flux_vr(c,j,icat)*dt, & 
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
            mass_to_meq(col_mf%cec_cation_flux_vr(c,j,1)*dt/EWParamsInst%cations_mass(1) & 
              *mass_h*EWParamsInst%cations_valence(1), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
            mass_to_meq(col_mf%cec_cation_flux_vr(c,j,2)*dt/EWParamsInst%cations_mass(2) & 
              *mass_h*EWParamsInst%cations_valence(2), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
            mass_to_meq(col_mf%cec_cation_flux_vr(c,j,3)*dt/EWParamsInst%cations_mass(3) & 
              *mass_h*EWParamsInst%cations_valence(3), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
            mass_to_meq(col_mf%cec_cation_flux_vr(c,j,4)*dt/EWParamsInst%cations_mass(4) & 
              *mass_h*EWParamsInst%cations_valence(4), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
            mass_to_meq(col_mf%cec_cation_flux_vr(c,j,5)*dt/EWParamsInst%cations_mass(5) & 
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
              -col_mf%cec_cation_flux_vr(c,j,icat)*dt, col_mf%background_cec_vr(c,j,icat)*dt
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
          if (col_ms%cation_vr(c,j,icat) < 0) then
            if (col_ms%cation_vr(c,j,icat) > -1e-12_r8) then
              ! Numerical accuracy problems can cause runoff to leach away all the cations
              ! Reset to zero in this case. Balance Check will ignore everything smaller than 1e-12
              col_ms%cation_vr(c,j,icat) = 0._r8
            else
              ! write (100+iam, *) 'cation_vr diagnostics:', ldomain%latc(g), ldomain%lonc(g), c, j, icat, col_ms%cation_vr(c,j,icat), trim(dateTimeString)

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
                - col_mf%secondary_cation_flux_vr(c,j,icat)*dt, &
                - col_mf%cation_uptake_vr(c,j,icat)*dt, &
                col_mf%cation_infl_vr(c,j,icat)*dt, &
                - col_mf%cation_leached_vr(c,j,icat)*dt, &
                - col_mf%cation_runoff_vr(c,j,icat)*dt

              write (100+iam, *) 'Post-reaction cec H+: ', ldomain%latc(g), ldomain%lonc(g), trim(dateTimeString)
              write (100+iam, *) c, j, col_ms%cec_proton_vr(c,j), & 
                mass_to_meq(col_ms%cec_proton_vr(c,j), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
                mass_to_meq(col_mf%cec_cation_flux_vr(c,j,1)*dt/EWParamsInst%cations_mass(1) & 
                  *mass_h*EWParamsInst%cations_valence(1), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
                mass_to_meq(col_mf%cec_cation_flux_vr(c,j,2)*dt/EWParamsInst%cations_mass(2) & 
                  *mass_h*EWParamsInst%cations_valence(2), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
                mass_to_meq(col_mf%cec_cation_flux_vr(c,j,3)*dt/EWParamsInst%cations_mass(3) & 
                  *mass_h*EWParamsInst%cations_valence(3), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
                mass_to_meq(col_mf%cec_cation_flux_vr(c,j,4)*dt/EWParamsInst%cations_mass(4) & 
                  *mass_h*EWParamsInst%cations_valence(4), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), &
                mass_to_meq(col_mf%cec_cation_flux_vr(c,j,5)*dt/EWParamsInst%cations_mass(5) & 
                *mass_h*EWParamsInst%cations_valence(5), 1._r8, mass_h, soilstate_vars%bd_col(c,j))

              write (100+iam, *) 'Post-reaction cec cation: ', ldomain%latc(g), ldomain%lonc(g), trim(dateTimeString)
              write (100+iam, *) c, j, icat, col_ms%cec_cation_vr(c,j,icat), &
                mass_to_meq(col_ms%cec_cation_vr(c,j,icat), EWParamsInst%cations_valence(icat), &
                            EWParamsInst%cations_mass(icat), soilstate_vars%bd_col(c,j)), &
                -col_mf%cec_cation_flux_vr(c,j,icat)*dt
              !------------------------------------------------------------------------------

              call endrun(msg='cation_vr < 0')
            end if
          end if
        end do
      end do
    end do
    !------------------------------------------------------------------------------
  end subroutine MineralStateDiags

  !-----------------------------------------------------------------------
  subroutine MineralFluxLimit(num_soilc, filter_soilc, col_ms, col_mf, dt)
    !
    ! !DESCRIPTION:
    ! Scale down reaction flux rates if they cause negative cation balance or
    ! exceed total cation exchange capacity
    !
    !$acc routine seq
    ! !ARGUMENTS:
    integer                      , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                      , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(column_mineral_state)   , intent(inout) :: col_ms
    type(column_mineral_flux)    , intent(inout) :: col_mf
    real(r8)                     , intent(in)    :: dt              ! radiation time step (seconds)
    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,icat,m,g ! indices
    integer  :: fc        ! lake filter indices
    integer  :: nlevbed
    real(r8) :: residual_factor
    real(r8) :: temp_delta1_cece(1:num_soilc, 1:nlevsoi, 1:ncations)
    real(r8) :: temp_delta2_cece(1:num_soilc, 1:nlevsoi, 1:ncations)
    real(r8) :: temp_delta_ceca(1:num_soilc, 1:nlevsoi)
    real(r8) :: temp_delta1_cation(1:num_soilc, 1:nlevsoi, 1:ncations)
    real(r8) :: temp_delta2_cation(1:num_soilc, 1:nlevsoi, 1:ncations)
    real(r8) :: min_flux_limit(1:num_soilc, 1:nlevsoi)
    logical  :: err_found
    integer  :: err_fc, err_lev, err_icat, err_col
    character(len=256) :: dateTimeString
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

        ! Limit due to cation exchange capacity of individual non-proton cations
        do icat = 1,ncations
          ! cec_cation_flux_vr > 0 := flow from CEC to solution
          temp_delta1_cece(fc,j,icat) = col_mf%background_cec_vr(c,j,icat)*dt
          temp_delta2_cece(fc,j,icat) = - col_mf%cec_cation_flux_vr(c,j,icat)*dt

          if ((col_ms%cec_cation_vr(c,j,icat) + temp_delta1_cece(fc,j,icat) + &
               temp_delta2_cece(fc,j,icat)) < 0._r8) then
            col_mf%cec_limit_vr(c,j,icat) = - (temp_delta1_cece(fc,j,icat) + & 
                                               col_ms%cec_cation_vr(c,j,icat)) / &
                temp_delta2_cece(fc,j,icat) * residual_factor
            col_mf%cec_cation_flux_vr(c,j,icat) = col_mf%cec_cation_flux_vr(c,j,icat) * & 
                                  col_mf%cec_limit_vr(c,j,icat)
            min_flux_limit(fc,j) = min(min_flux_limit(fc,j), col_mf%cec_limit_vr(c,j,icat))
          else
            col_mf%cec_limit_vr(c,j,icat) = 1._r8
          end if
        end do

        ! Limit due to cation exchange capacity of H+
        temp_delta_ceca(fc,j) = 0._r8
        do icat = 1,ncations
          temp_delta_ceca(fc,j) = temp_delta_ceca(fc,j) + &
           (col_mf%cec_cation_flux_vr(c,j,icat) - col_mf%background_cec_vr(c,j,icat)) * dt & 
           / EWParamsInst%cations_mass(icat) * mass_h * EWParamsInst%cations_valence(icat)
        end do
        if ((col_ms%cec_proton_vr(c,j) + temp_delta_ceca(fc,j)) < 0._r8) then        
          col_mf%proton_limit_vr(c,j) = - col_ms%cec_proton_vr(c,j) / temp_delta_ceca(fc,j) & 
                              * residual_factor
          do icat = 1,ncations
            col_mf%cec_cation_flux_vr(c,j,icat) = col_mf%cec_cation_flux_vr(c,j,icat) * &
                                      col_mf%proton_limit_vr(c,j)
            col_mf%background_cec_vr(c,j,icat) = col_mf%background_cec_vr(c,j,icat) * &
                                      col_mf%proton_limit_vr(c,j)
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
                                          col_mf%background_flux_vr(c,j,icat)*dt
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
                write (100+iam, *) '   negative CEC cation ', icat, col_mf%cec_limit_vr(c,j,icat), col_ms%cec_cation_vr(c,j,icat), - col_mf%cec_cation_flux_vr(c,j,icat)*dt, col_mf%background_cec_vr(c,j,icat)*dt
              end if
            end do
            if (col_mf%proton_limit_vr(c,j) < 1._r8) then
              write (100+iam, *) '   negative CEC H+ ', col_mf%proton_limit_vr(c,j), col_ms%cec_proton_vr(c,j), temp_delta_ceca(fc,j)
            end if
            do icat = 1,ncations
              if (col_mf%flux_limit_vr(c,j,icat) < 1._r8) then
                write (100+iam, *) '   negative solution cation ', icat, col_mf%flux_limit_vr(c,j,icat), col_ms%cation_vr(c,j,icat), col_mf%primary_cation_flux_vr(c,j,icat)*dt, col_mf%background_flux_vr(c,j,icat)*dt, col_mf%cation_uptake_vr(c,j,icat)*dt, -col_mf%secondary_cation_flux_vr(c,j,icat)*dt, col_mf%cec_cation_flux_vr(c,j,icat)*dt
              end if
            end do
          end if
          do icat = 1,ncations
            if (col_mf%flux_limit_vr(c,j,icat) < 1._r8) then
              write (100+iam, *) '   negative solution cation ', icat, col_mf%flux_limit_vr(c,j,icat), col_ms%cation_vr(c,j,icat), temp_delta1_cation(fc,j,icat), temp_delta2_cation(fc,j,icat)
            end if
          end do
        end do
      end do

    end if

  end subroutine MineralFluxLimit

end module MineralStateUpdateMod
