module MineralStateUpdateMod

  !-----------------------------------------------------------------------
  ! Module for enhanced weathering state variable update.
  !
  ! !USES:
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use decompMod               , only : bounds_type
  use elm_varpar              , only : nminerals, ncations, nminsecs, nlevgrnd, mixing_layer
  use elm_varcon              , only : zisoi, dzsoi, mass_h
  use elm_varctl              , only : iulog
  use abortutils              , only : endrun
  use shr_log_mod             , only : errMsg => shr_log_errMsg
  use ewutils                 , only : mass_to_mol, mass_to_meq, mol_to_mass
  use ColumnDataType          , only : col_ws
  use ColumnDataType          , only : col_ms, col_mf, col_pp
  use ColumnDataType          , only : column_mineral_state, column_mineral_flux
  use SoilStateType           , only : soilstate_type
  use EnhancedWeatheringMod   , only : EWParamsInst
  use domainMod               , only : ldomain ! debug print
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: MineralStateUpdate
  public :: MineralFluxLimit
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine MineralStateUpdate(num_soilc, filter_soilc, col_ms, col_mf, dt, soilstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic mineral state variables
    !
    !$acc routine seq
    ! !ARGUMENTS:
    integer                      , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                      , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(column_mineral_state)   , intent(inout) :: col_ms
    type(column_mineral_flux)    , intent(inout) :: col_mf
    real(r8)                     , intent(in)    :: dt              ! radiation time step (seconds)

    ! DEBUG
    type(soilstate_type)         , intent(in)    :: soilstate_vars

    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,k,a,m,g ! indices
    integer  :: fp,fc         ! lake filter indices
    real(r8) :: flux_limit
    real(r8) :: temp_solution_balance
    !-----------------------------------------------------------------------

    ! Update mineral state
    do fc = 1,num_soilc
      c = filter_soilc(fc)

      col_mf%r_sequestration(c) = 0._r8

      do j = 1,mixing_layer

        ! -----------------------------------------------------------------------------------
        ! Balance update
        ! -----------------------------------------------------------------------------------

        ! cations
        do a = 1,ncations
          col_ms%cation_vr(c,j,a) = col_ms%cation_vr(c,j,a) + &
            col_mf%background_weathering_vr(c,j,a)*dt + col_mf%primary_cation_flux_vr(c,j,a)*dt - &
            col_mf%secondary_cation_flux_vr(c,j,a)*dt + col_mf%cec_cation_flux_vr(c,j,a)*dt + &
            col_mf%cation_infl_vr(c,j,a)*dt - col_mf%cation_oufl_vr(c,j,a)*dt - &
            col_mf%cation_uptake_vr(c,j,a)*dt - col_mf%cation_leached_vr(c,j,a)*dt - col_mf%cation_runoff_vr(c,j,a)*dt

          col_ms%cec_cation_vr(c,j,a) = col_ms%cec_cation_vr(c,j,a) - col_mf%cec_cation_flux_vr(c,j,a)*dt
        end do

        ! soil H+ concentration (g m-3 soil)
        !! col_ms%proton_vr(c,j) = col_ms%proton_vr(c,j) - col_mf%primary_proton_flux_vr(c,j)*dt + col_mf%cec_proton_flux_vr(c,j)*dt + col_mf%proton_infl_vr(c,j)*dt - col_mf%proton_oufl_vr(c,j)*dt - col_mf%proton_uptake_vr(c,j)*dt - col_mf%proton_leached_vr(c,j)*dt - col_mf%proton_runoff_vr(c,j)*dt
        !! col_ms%soil_ph(c,j) = - log10(mass_to_mol(col_ms%proton_vr(c,j), 1._r8, col_ws%h2osoi_vol(c,j)))
        col_ms%proton_vr(c,j) = mol_to_mass(10**(-col_ms%soil_ph(c,j)), mass_h, col_ws%h2osoi_liq(c,j))

        !write (iulog, *) 'soil_ph', c, j, col_ms%soil_ph(c,j), col_ms%proton_vr(c,j), col_ws%h2osoi_vol(c,j)

        ! primary mineral
        do m = 1,nminerals
          col_ms%primary_mineral_vr(c,j,m) = col_ms%primary_mineral_vr(c,j,m) + &
            col_mf%primary_added_vr(c,j,m)*dt - col_mf%primary_dissolve_vr(c,j,m)*dt
          ! non-SiO2 minerals
          col_ms%primary_residue_vr(c,j,m) = col_ms%primary_residue_vr(c,j,m) + &
            col_mf%primary_residue_flux_vr(c,j,m)*dt
        end do

        ! note "cec_proton_flux_vr" integrates the impacts from CO2 reactions
        ! instead, use charge balance on the mineral surface to get the change in adsorped H+
        do a = 1,ncations
          col_ms%cec_proton_vr(c,j) = col_ms%cec_proton_vr(c,j) + &
            col_mf%cec_cation_flux_vr(c,j,a)*dt/EWParamsInst%cations_mass(a)*mass_h*EWParamsInst%cations_valence(a)
        end do

        ! silicate
        col_ms%silica_vr(c,j) = col_ms%silica_vr(c,j) + col_mf%primary_silica_flux_vr(c,j) * dt - col_mf%secondary_silica_flux_vr(c,j) * dt

        ! secondary mineral
        do m = 1,nminsecs
            col_ms%secondary_mineral_vr(c,j,m) = col_ms%secondary_mineral_vr(c,j,m) + col_mf%secondary_mineral_flux_vr(c,j,m) * dt
        end do

        ! TODO: ignore the effect on soil water for now

        ! Calculate the total CO2 sequestration rate in mol m-2 s-1
        ! precipitated by calcite: 1 mol CO2 per mol Ca2+
        col_mf%r_sequestration(c) = col_mf%r_sequestration(c) + col_mf%secondary_cation_flux_vr(c,j,1) / EWParamsInst%cations_mass(1)
        ! transported to ocean: 2x for 2+ cations, 1x for 1+ cations
        do a = 1,ncations
          col_mf%r_sequestration(c) = col_mf%r_sequestration(c) + & 
              (col_mf%cation_leached_vr(c,j,a) + col_mf%cation_leached_vr(c,j,a)) / &
              EWParamsInst%cations_mass(a) * EWParamsInst%cations_valence(a)
        end do
      end do

      ! convert from mol m-2 s-1 to gC m-2 s-1
      col_mf%r_sequestration(c) = col_mf%r_sequestration(c) * 12._r8
      ! multiply by oceanic efficiency
      col_mf%r_sequestration(c) = col_mf%r_sequestration(c) * 0.86

    end do

    !write (iulog, *) 'Post-reaction H+'
    !do j = 1,mixing_layer
    !  write (iulog, *) c, j, col_ms%soil_ph(c,j), col_ms%proton_vr(c,j), mass_to_mol(col_ms%proton_vr(c,j), mass_h, col_ws%h2osoi_vol(c,j)), - col_mf%primary_proton_flux_vr(c,j)*dt, col_mf%cec_proton_flux_vr(c,j)*dt, col_mf%proton_infl_vr(c,j)*dt, - col_mf%proton_oufl_vr(c,j)*dt, -col_mf%proton_uptake_vr(c,j)*dt, -col_mf%proton_leached_vr(c,j)*dt, -col_mf%proton_runoff_vr(c,j)*dt
    !end do

    !write (iulog, *) 'Post-reaction cation'
    !do j = 1,mixing_layer
    !  do a = 1, ncations
    !    write (iulog, *) c, j, a, col_ms%cation_vr(c,j,a), mass_to_mol(col_ms%cation_vr(c,j,a), EWParamsInst%cations_mass(a), col_ws%h2osoi_vol(c,j)), col_mf%background_weathering_vr(c,j,a)*dt, col_mf%primary_cation_flux_vr(c,j,a)*dt, - col_mf%secondary_cation_flux_vr(c,j,a)*dt, col_mf%cec_cation_flux_vr(c,j,a)*dt, col_mf%cation_infl_vr(c,j,a)*dt, -col_mf%cation_oufl_vr(c,j,a)*dt, - col_mf%cation_uptake_vr(c,j,a)*dt, - col_mf%cation_leached_vr(c,j,a)*dt, - col_mf%cation_runoff_vr(c,j,a)*dt
    !  end do
    !end do

    !write (iulog, *) 'Post-reaction cec H+'
    !do j = 1,mixing_layer
    !  write (iulog, *) c, j, col_ms%cec_proton_vr(c,j), mass_to_meq(col_ms%cec_proton_vr(c,j), 1._r8, mass_h, soilstate_vars%bd_col(c,j)), mass_to_meq(col_mf%cec_cation_flux_vr(c,j,1)*dt/EWParamsInst%cations_mass(1)*mass_h*EWParamsInst%cations_valence(1) + col_mf%cec_cation_flux_vr(c,j,2)*dt/EWParamsInst%cations_mass(2)*mass_h*EWParamsInst%cations_valence(2) + col_mf%cec_cation_flux_vr(c,j,3)*dt/EWParamsInst%cations_mass(3)*mass_h*EWParamsInst%cations_valence(3) + col_mf%cec_cation_flux_vr(c,j,4)*dt/EWParamsInst%cations_mass(4)*mass_h*EWParamsInst%cations_valence(4) + col_mf%cec_cation_flux_vr(c,j,5)*dt/EWParamsInst%cations_mass(5)*mass_h*EWParamsInst%cations_valence(5), 1._r8, mass_h, soilstate_vars%bd_col(c,j))
    !end do

    !write (iulog, *) 'Post-reaction cec cation'
    !do j = 1,mixing_layer
    !  do a = 1, ncations
    !    write (iulog, *) c, j, a, col_ms%cec_cation_vr(c,j,a), mass_to_meq(col_ms%cec_cation_vr(c,j,a), EWParamsInst%cations_valence(a), EWParamsInst%cations_mass(a), soilstate_vars%bd_col(c,j)), -col_mf%cec_cation_flux_vr(c,j,a)*dt
    !  end do
    !end do

    !do j = 1,mixing_layer
    !  do a = 1,ncations
    !    if (col_ms%cation_vr(c,j,a) < 0) then
    !      write (iulog, *) c, j, a, col_mf%cec_cation_flux_vr(c,j,a)*dt
    !      call endrun(msg='cation_vr < 0')
    !    end if
    !  end do
    !end do

  end subroutine MineralStateUpdate


  !-----------------------------------------------------------------------
  subroutine MineralFluxLimit(num_soilc, filter_soilc, col_ms, col_mf, dt)
    !
    ! !DESCRIPTION:
    ! Scale down reaction flux rates if they cause negative cation balance
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
    integer  :: c,j,a,m,g ! indices
    integer  :: fc        ! lake filter indices
    real(r8) :: flux_limit
    real(r8) :: temp_in, temp_out
    character(len=32) :: subname = 'elm_erw_mineral_flux_limit'  ! subroutine name
    !-----------------------------------------------------------------------

    do fc = 1,num_soilc
      c = filter_soilc(fc)
      g = col_pp%gridcell(c)

      do j = 1,mixing_layer

        do a = 1,ncations
          ! Limit due to cation exchange capacity
          temp_in = 0._r8
          temp_out = - col_mf%cec_cation_flux_vr(c,j,a)*dt
          if ((col_ms%cec_cation_vr(c,j,a) + temp_in + temp_out) < 0._r8) then
            ! ensure a tiny bit of cation is left due to numerical accuracy reasons
            flux_limit = - (temp_in + col_ms%cec_cation_vr(c,j,a)) / temp_out * 0.99_r8

            col_mf%cec_cation_flux_vr(c,j,a) = col_mf%cec_cation_flux_vr(c,j,a) * flux_limit

            write (iulog, *) 'Flux limit due to negative CEC cation; factor = ', c, j, a, flux_limit, col_mf%cec_cation_flux_vr(c,j,a)
          end if
        end do

        do a = 1,ncations
          ! Limit due to soil solution cation concentration
          temp_in = col_mf%background_weathering_vr(c,j,a)*dt + &
                    col_mf%primary_cation_flux_vr(c,j,a)*dt + &
                    col_mf%cation_infl_vr(c,j,a)*dt - & 
                    col_mf%cation_oufl_vr(c,j,a)*dt - &
                    col_mf%cation_uptake_vr(c,j,a)*dt - & 
                    col_mf%cation_leached_vr(c,j,a)*dt - &
                    col_mf%cation_runoff_vr(c,j,a)*dt
          temp_out = - col_mf%secondary_cation_flux_vr(c,j,a)*dt + & 
                       col_mf%cec_cation_flux_vr(c,j,a)*dt

          if ((col_ms%cation_vr(c,j,a) + temp_in + temp_out) < 0._r8) then
            if (col_ms%cation_vr(c,j,a) + temp_in < 0._r8) then
              write (iulog, *) ldomain%latc(g), ldomain%lonc(g), g, c, j, a, 'Problematic flushing rate: ', 'dt=', dt, 'initial cation=', col_ms%cation_vr(c,j,a), 'terms=', col_mf%background_weathering_vr(c,j,a)*dt, col_mf%primary_cation_flux_vr(c,j,a)*dt, col_mf%cation_infl_vr(c,j,a)*dt, col_mf%cation_oufl_vr(c,j,a)*dt, col_mf%cation_uptake_vr(c,j,a)*dt, col_mf%cation_leached_vr(c,j,a)*dt, col_mf%cation_runoff_vr(c,j,a)*dt
              call endrun(msg=subname //':: ERROR: Negative cation balance'//errMsg(__FILE__, __LINE__))
            end if

            ! ensure a tiny bit of cation is left due to numerical accuracy reasons
            flux_limit = - (temp_in + col_ms%cation_vr(c,j,a)) / temp_out * 0.99_r8

            col_mf%secondary_cation_flux_vr(c,j,a) = col_mf%secondary_cation_flux_vr(c,j,a) * flux_limit
            col_mf%cec_cation_flux_vr(c,j,a) = col_mf%cec_cation_flux_vr(c,j,a) * flux_limit
            !col_mf%cation_infl_vr(c,j,a) = col_mf%cation_infl_vr(c,j,a) * flux_limit
            !col_mf%cation_oufl_vr(c,j,a) = col_mf%cation_oufl_vr(c,j,a) * flux_limit
            !col_mf%cation_uptake_vr(c,j,a) = col_mf%cation_uptake_vr(c,j,a) * flux_limit
            !col_mf%cation_leached_vr(c,j,a) = col_mf%cation_leached_vr(c,j,a) * flux_limit
            !col_mf%cation_runoff_vr(c,j,a) = col_mf%cation_runoff_vr(c,j,a) * flux_limit

            write (iulog, *) 'Flux limit due to negative cation concentration; factor = ', ldomain%latc(g), ldomain%lonc(g), g, c, j, a, flux_limit, col_mf%secondary_cation_flux_vr(c,j,a), col_mf%cec_cation_flux_vr(c,j,a)
          end if
        end do

        ! Limit due to acid exchange capacity
        temp_in = 0._r8
        temp_out = 0._r8
        do a = 1,ncations
          temp_out = temp_out + col_mf%cec_cation_flux_vr(c,j,a)*dt/EWParamsInst%cations_mass(a)*mass_h*EWParamsInst%cations_valence(a)
        end do
        if ((col_ms%cec_proton_vr(c,j) + temp_in + temp_out) < 0._r8) then        
            ! ensure a tiny bit of H+ is left due to numerical accuracy reasons
            flux_limit = - (temp_in + col_ms%cec_proton_vr(c,j)) / temp_out * 0.99_r8

            do a = 1,ncations
              col_mf%cec_cation_flux_vr(c,j,a) = col_mf%cec_cation_flux_vr(c,j,a) * flux_limit
            end do

            write (iulog, *) 'Flux limit due to negative CEC H+; factor = ', c, j, flux_limit, col_mf%cec_cation_flux_vr(c,j,1), col_mf%cec_cation_flux_vr(c,j,2), col_mf%cec_cation_flux_vr(c,j,3), col_mf%cec_cation_flux_vr(c,j,4), col_mf%cec_cation_flux_vr(c,j,5)
        end if
      end do
    end do
  end subroutine MineralFluxLimit

end module MineralStateUpdateMod
