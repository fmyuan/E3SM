module MineralStateUpdateMod

  !-----------------------------------------------------------------------
  ! Module for enhanced weathering state variable update.
  !
  ! !USES:
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use decompMod               , only : bounds_type
  use elm_varpar              , only : nminerals, ncations, nminsec, nlevgrnd, mixing_depth
  use elm_varcon              , only : zisoi, dzsoi
  !use elm_varctl              , only : nu_com
  !use elm_varctl              , only : use_pflotran, pf_cmode, use_fates

  !use GridcellDataType        , only : grc_cs, c13_grc_cs, c14_grc_cs
  !use GridcellDataType        , only : grc_cf, c13_grc_cf, c14_grc_cf
  use ColumnDataType          , only : col_ms, col_mf, col_pp
  use ColumnDataType          , only : column_mineral_state, column_mineral_flux
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: MineralStateUpdate
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine MineralStateUpdate(num_soilc, filter_soilc, col_ms, col_mf, dt)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic carbon state
    ! variables (except for gap-phase mortality and fire fluxes)
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
    integer  :: c,p,j,k,l,a,m ! indices
    integer  :: fp,fc         ! lake filter indices
    integer  :: maxlayer      ! deepest layer for mineral reactions
    !-----------------------------------------------------------------------

    ! Update mineral state
    do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! find the maximum layer of reaction
      maxlayer = 1
      do l = 2,nlevgrnd
        if (zisoi(l-1) < mixing_depth) then
          maxlayer = l
        end if
      end do

      do l = 1,maxlayer
          ! primary mineral
          do m = 1,nminerals
              col_ms%primary_mineral_vr(c,l,m) = col_ms%primary_mineral_vr(c,l,m) + &
                (col_mf%primary_added_vr(c,l,m) - col_mf%primary_dissolve_vr(c,l,m)) * dt
          end do

          ! cations
          do a = 1,ncations
              col_ms%cation_vr(c,l,a) = col_ms%cation_vr(c,l,a) + &
                (col_mf%primary_cation_flux_vr(c,l,a) - &
                 col_mf%secondary_cation_flux_vr(c,l,a) - &
                 col_mf%cation_leached_vr(c,l,a) - &
                 col_mf%cation_runoff_vr(c,l,a)) * dt
          end do

          ! bicarbonate
          col_ms%bicarbonate_vr(c,l) = col_ms%bicarbonate_vr(c,l) + &
            (col_mf%primary_bicarbonate_flux_vr(c,l) - &
             col_mf%secondary_bicarbonate_flux_vr(c,l) - &
             col_mf%bicarbonate_leached_vr(c,l) - &
             col_mf%bicarbonate_runoff_vr(c,l)) * dt

          ! silicate
          col_ms%silicate_vr(c,l) = col_ms%silicate_vr(c,l) + &
            col_mf%primary_silicate_flux_vr(c,l) * dt

          ! secondary mineral
          do m = 1,nminsec
              col_ms%secondary_mineral_vr(c,l,m) = col_ms%secondary_mineral_vr(c,l,m) + col_mf%secondary_precip_vr(c,l,m) * dt
          end do

          ! TODO: ignore the effect on soil water for now
      end do
    end do

  end subroutine MineralStateUpdate

end module MineralStateUpdateMod