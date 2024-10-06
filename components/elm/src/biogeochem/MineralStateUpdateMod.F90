module MineralStateUpdateMod

  !-----------------------------------------------------------------------
  ! Module for enhanced weathering state variable update.
  !
  ! !USES:
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use decompMod               , only : bounds_type
  use spmdMod                 , only : iam
  use elm_varpar              , only : nminerals, ncations, nminsecs, nlevgrnd, nlevsoi
  use elm_varcon              , only : zisoi, dzsoi, mass_h
  use elm_varctl              , only : use_erw_verbose
  use shr_sys_mod             , only : shr_sys_flush
  use spmdMod                 , only : masterproc
  use abortutils              , only : endrun
  use shr_log_mod             , only : errMsg => shr_log_errMsg
  use ewutils                 , only : mass_to_mol, mass_to_meq, mol_to_mass
  use ColumnDataType          , only : col_ws
  use ColumnDataType          , only : col_ms, col_mf, col_pp
  use ColumnDataType          , only : column_mineral_state, column_mineral_flux
  use SoilStateType           , only : soilstate_type
  use EnhancedWeatheringMod   , only : EWParamsInst
  use domainMod               , only : ldomain ! debug print
  use elm_time_manager        , only : get_curr_time_string
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
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine MineralStateUpdate1(num_soilc, filter_soilc, col_ms, col_mf, dt)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the prognostic mineral state variables that do not
    ! depend on vertical water movement or leaching
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
        ! -----------------------------------------------------------------------------------
        ! Balance update
        ! -----------------------------------------------------------------------------------
        ! CEC cations
        do icat = 1,ncations
          col_ms%cec_cation_vr(c,j,icat) = col_ms%cec_cation_vr(c,j,icat) - &
                  col_mf%cec_cation_flux_vr(c,j,icat)*dt
        end do

        ! CEC H+
        ! note "cec_proton_flux_vr" integrates the impacts from CO2 reactions
        ! instead, use charge balance on the mineral surface to get the change in adsorped H+
        do icat = 1,ncations
          col_ms%cec_proton_vr(c,j) = col_ms%cec_proton_vr(c,j) + &
            col_mf%cec_cation_flux_vr(c,j,icat)*dt/EWParamsInst%cations_mass(icat)*mass_h*EWParamsInst%cations_valence(icat)
        end do

        ! soil H+ concentration (g m-3 soil)
        ! only determined by CEC equilibrium
        !! col_ms%proton_vr(c,j) = col_ms%proton_vr(c,j) - col_mf%primary_proton_flux_vr(c,j)*dt + col_mf%cec_proton_flux_vr(c,j)*dt + col_mf%proton_infl_vr(c,j)*dt - col_mf%proton_oufl_vr(c,j)*dt - col_mf%proton_uptake_vr(c,j)*dt - col_mf%proton_leached_vr(c,j)*dt - col_mf%proton_runoff_vr(c,j)*dt
        !! col_ms%soil_ph(c,j) = - log10(mass_to_mol(col_ms%proton_vr(c,j), 1._r8, col_ws%h2osoi_vol(c,j)))
        col_ms%proton_vr(c,j) = mol_to_mass(10**(-col_ms%soil_ph(c,j)), mass_h, col_ws%h2osoi_vol(c,j))

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
    ! On the radiation time step, update the prognostic mineral state variables
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
          col_ms%cation_vr(c, j, icat) = col_ms%cation_vr(c, j, icat) + ( &
            col_mf%background_weathering_vr(c,j,icat) + & 
            col_mf%primary_cation_flux_vr(c,j,icat) + col_mf%cec_cation_flux_vr(c,j,icat) - &
            col_mf%secondary_cation_flux_vr(c,j,icat) - col_mf%cation_uptake_vr(c,j,icat) + & 
            col_mf%cation_infl_vr(c,j,icat) - col_mf%cation_oufl_vr(c,j,icat) ) * dt
        end do
        !write (100+iam, *) 'post-adv', c, icat, cation_vr(c,1:mixing_layer, icat)
      end do
    end do
  end subroutine MineralStateUpdate2

  !-----------------------------------------------------------------------
  subroutine MineralStateUpdate3(num_soilc, filter_soilc, col_ms, col_mf, dt, soilstate_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the prognostic mineral state variables
    ! related to leaching
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
    integer  :: c,p,j,k,icat,m,g ! indices
    integer  :: fp,fc         ! lake filter indices
    integer  :: nlevbed
    logical  :: print_proton, print_cations, print_cec_proton, print_cec_cations ! flags to turn on debug printing
    character(len=256) :: dateTimeString
    !-----------------------------------------------------------------------

    call get_curr_time_string(dateTimeString)

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
      col_mf%r_sequestration(c) = 0._r8
      do j = 1,nlevbed
        ! precipitated by calcite: 1 mol CO2 per mol Ca2+
        col_mf%r_sequestration(c) = col_mf%r_sequestration(c) + & 
          col_mf%secondary_cation_flux_vr(c,j,1) / EWParamsInst%cations_mass(1) * col_pp%dz(c,j)
        ! transported to ocean: 2x for 2+ cations, 1x for 1+ cations, multiply by
        ! ocean efficiency (0.86)
        ! - col_mf%background_weathering_vr(c,j,icat)
        do icat = 1,ncations
          col_mf%r_sequestration(c) = col_mf%r_sequestration(c) + & 
              ( col_mf%cation_leached_vr(c,j,icat) + col_mf%cation_runoff_vr(c,j,icat) ) & 
              * 0.86_r8 * col_pp%dz(c,j) / &
              EWParamsInst%cations_mass(icat) * EWParamsInst%cations_valence(icat)
        end do
      end do
      ! convert from mol m-2 s-1 to gC m-2 s-1
      col_mf%r_sequestration(c) = col_mf%r_sequestration(c) * 12._r8
    end do

    !------------------------------------------------------------------------------
    ! Write mass balance diagnostics only if verbose level is set to high
    if (use_erw_verbose == 2) then
      ! print soil solution proton
      write (100+iam, *) 'Post-reaction H+: ', ldomain%latc(g), ldomain%lonc(g), trim(dateTimeString)
      do fc = 1,num_soilc
        c = filter_soilc(fc)
        nlevbed = min(col_pp%nlevbed(c), nlevsoi)
        do j = 1,nlevbed
          write (100+iam, *) c, j, col_ms%soil_ph(c,j), col_ms%proton_vr(c,j), & 
             mass_to_mol(col_ms%proton_vr(c,j), mass_h, col_ws%h2osoi_vol(c,j))
        end do
      end do
      !!  - col_mf%primary_proton_flux_vr(c,j)*dt, & 
      !!  col_mf%cec_proton_flux_vr(c,j)*dt, &
      !!  col_mf%proton_infl_vr(c,j)*dt, &
      !!  - col_mf%proton_oufl_vr(c,j)*dt, & 
      !!  -col_mf%proton_uptake_vr(c,j)*dt, & 
      !!  -col_mf%proton_leached_vr(c,j)*dt, &
      !!  -col_mf%proton_runoff_vr(c,j)*dt

      ! print soil solution cations
      write (100+iam, *) 'Post-reaction cation: ', ldomain%latc(g), ldomain%lonc(g), trim(dateTimeString)
      do fc = 1,num_soilc
        c = filter_soilc(fc)
        nlevbed = min(col_pp%nlevbed(c), nlevsoi)
        ! note: sourcesink_cations term in EnhancedWeatheringMod.F90
        !       should be approximately matched to cation_infl_vr
        ! leaching & runoff are the additionals
        do icat = 1,ncations
          do j = 1,nlevbed
            write (100+iam, *) c, j, icat, col_ms%cation_vr(c,j,icat), &
              mass_to_mol(col_ms%cation_vr(c,j,icat), EWParamsInst%cations_mass(icat), &
                          col_ws%h2osoi_vol(c,j)), &
              (col_mf%background_weathering_vr(c,j,icat) + col_mf%primary_cation_flux_vr(c,j,icat) & 
              + col_mf%cec_cation_flux_vr(c,j,icat) - col_mf%secondary_cation_flux_vr(c,j,icat) &
              - col_mf%cation_uptake_vr(c,j,icat))*dt, &
              col_mf%cation_infl_vr(c,j,icat)*dt, &
              - col_mf%cation_leached_vr(c,j,icat)*dt, &
              - col_mf%cation_runoff_vr(c,j,icat)*dt
            end do
          end do
        end do

      ! print CEC protons
      write (100+iam, *) 'Post-reaction cec H+: ', ldomain%latc(g), ldomain%lonc(g), trim(dateTimeString)
      do fc = 1,num_soilc
        c = filter_soilc(fc)
        nlevbed = min(col_pp%nlevbed(c), nlevsoi)
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
      write (100+iam, *) 'Post-reaction cec cation: ', ldomain%latc(g), ldomain%lonc(g), trim(dateTimeString)
      do fc = 1,num_soilc
        c = filter_soilc(fc)
        nlevbed = min(col_pp%nlevbed(c), nlevsoi)
        do j = 1,nlevbed
          do icat = 1,ncations
            write (100+iam, *) c, j, icat, col_ms%cec_cation_vr(c,j,icat), &
              mass_to_meq(col_ms%cec_cation_vr(c,j,icat), EWParamsInst%cations_valence(icat), &
                          EWParamsInst%cations_mass(icat), soilstate_vars%bd_col(c,j)), &
              -col_mf%cec_cation_flux_vr(c,j,icat)*dt
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
                (col_mf%background_weathering_vr(c,j,icat) + col_mf%primary_cation_flux_vr(c,j,icat) & 
                + col_mf%cec_cation_flux_vr(c,j,icat) - col_mf%secondary_cation_flux_vr(c,j,icat) &
                - col_mf%cation_uptake_vr(c,j,icat))*dt, &
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

  end subroutine MineralStateUpdate3

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
    integer  :: c,j,icat,m,g ! indices
    integer  :: fc        ! lake filter indices
    integer  :: nlevbed
    real(r8) :: residual_factor
    real(r8) :: temp_delta_cece(1:num_soilc, 1:nlevsoi, 1:ncations)
    real(r8) :: temp_delta_ceca(1:num_soilc, 1:nlevsoi)
    real(r8) :: temp_delta1_cation(1:num_soilc, 1:nlevsoi, 1:ncations)
    real(r8) :: temp_delta2_cation(1:num_soilc, 1:nlevsoi, 1:ncations)
    real(r8) :: min_flux_limit(1:num_soilc, 1:nlevsoi)
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

        ! Limit due to cation exchange capacity of individual non-proton cations
        do icat = 1,ncations
          ! cec_cation_flux_vr > 0 := flow from CEC to solution
          temp_delta_cece(fc,j,icat) = - col_mf%cec_cation_flux_vr(c,j,icat)*dt
          if ((col_ms%cec_cation_vr(c,j,icat) + temp_delta_cece(fc,j,icat)) < 0._r8) then
            col_mf%cec_limit_vr(c,j,icat) = - col_ms%cec_cation_vr(c,j,icat) / &
                                  temp_delta_cece(fc,j,icat) * residual_factor
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
          temp_delta_ceca(fc,j) = temp_delta_ceca(fc,j) + col_mf%cec_cation_flux_vr(c,j,icat) & 
            * dt / EWParamsInst%cations_mass(icat) * mass_h * EWParamsInst%cations_valence(icat)
        end do
        if ((col_ms%cec_proton_vr(c,j) + temp_delta_ceca(fc,j)) < 0._r8) then        
          col_mf%proton_limit_vr(c,j) = - col_ms%cec_proton_vr(c,j) / temp_delta_ceca(fc,j) & 
                              * residual_factor
          do icat = 1,ncations
            col_mf%cec_cation_flux_vr(c,j,icat) = col_mf%cec_cation_flux_vr(c,j,icat) * &
                                      col_mf%proton_limit_vr(c,j)
          end do
          min_flux_limit(fc,j) = min(min_flux_limit(fc,j), col_mf%proton_limit_vr(c,j))
        else
          col_mf%proton_limit_vr(c,j) = 1._r8
        end if

        ! Limit due to soil solution cation concentration, after the previous two limits
        ! have been applieds
        do icat = 1,ncations
          temp_delta1_cation(fc,j,icat) = col_mf%background_weathering_vr(c,j,icat)*dt + &
                    col_mf%primary_cation_flux_vr(c,j,icat)*dt + &
                    col_mf%cation_uptake_vr(c,j,icat)*dt
          temp_delta2_cation(fc,j,icat) = - col_mf%secondary_cation_flux_vr(c,j,icat)*dt + & 
                    col_mf%cec_cation_flux_vr(c,j,icat)*dt

          if (col_ms%cation_vr(c,j,icat) + temp_delta1_cation(fc,j,icat) < 0._r8) then
            err_found = .true.
            err_fc = fc
            err_col = filter_soilc(err_fc)
            err_lev = j
            err_icat = icat
          else if ((col_ms%cation_vr(c,j,icat) + temp_delta1_cation(fc,j,icat) + & 
                    temp_delta2_cation(fc,j,icat)) < 0._r8) then
            ! ensure a tiny bit of cation is left due to numerical accuracy reasons
            col_mf%flux_limit_vr(c,j,icat) = - (temp_delta1_cation(fc,j,icat) + & 
              col_ms%cation_vr(c,j,icat)) / temp_delta2_cation(fc,j,icat) * residual_factor

            col_mf%secondary_cation_flux_vr(c,j,icat) = col_mf%secondary_cation_flux_vr(c,j,icat)*& 
                  col_mf%flux_limit_vr(c,j,icat)
            col_mf%cec_cation_flux_vr(c,j,icat) = col_mf%cec_cation_flux_vr(c,j,icat) * & 
                  col_mf%flux_limit_vr(c,j,icat)

            min_flux_limit(fc,j) = min(min_flux_limit(fc,j), col_mf%flux_limit_vr(c,j,icat))
          else
            col_mf%flux_limit_vr(c,j,icat) = 1._r8
          end if
        end do
      end do
    end do

    if (err_found) then
      g = col_pp%gridcell(err_col)
      write (100+iam, *) 'Flushing rate diagnostics: ', ldomain%latc(g), ldomain%lonc(g), err_col, err_lev, err_icat, trim(dateTimeString)
      write (100+iam, *) ' initial cation=', col_ms%cation_vr(err_col, err_lev, err_icat)
      write (100+iam, *) ' delta1=', temp_delta1_cation(err_fc, err_lev, err_icat)
      write (100+iam, *) ' terms/dt=', &
        col_mf%background_weathering_vr(err_col, err_lev, err_icat), &
        col_mf%primary_cation_flux_vr(err_col, err_lev, err_icat), &
        col_mf%cation_uptake_vr(err_col, err_lev, err_icat)
      call endrun(msg=subname //':: ERROR: Problematic flushing rate'//errMsg(__FILE__, __LINE__))
    end if

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
                write (100+iam, *) '   negative CEC cation ', icat, col_mf%cec_limit_vr(c,j,icat), col_ms%cec_cation_vr(c,j,icat), temp_delta_cece(fc,j,icat)
              end if
            end do
            if (col_mf%proton_limit_vr(c,j) < 1._r8) then
              write (100+iam, *) '   negative CEC H+ ', col_mf%proton_limit_vr(c,j), col_ms%cec_proton_vr(c,j), temp_delta_ceca(fc,j)
            end if
            do icat = 1,ncations
              if (col_mf%flux_limit_vr(c,j,icat) < 1._r8) then
                write (100+iam, *) '   negative solution cation ', icat, col_mf%flux_limit_vr(c,j,icat), col_ms%cation_vr(c,j,icat), temp_delta1_cation(fc,j,icat), temp_delta2_cation(fc,j,icat)
              end if
            end do
          end if
        end do
      end do
    end if

  end subroutine MineralFluxLimit

end module MineralStateUpdateMod