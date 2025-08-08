module EMI_WaterFluxType_ExchangeMod
  !
  use shr_kind_mod                          , only : r8 => shr_kind_r8
  use shr_log_mod                           , only : errMsg => shr_log_errMsg
  use abortutils                            , only : endrun
  use elm_varctl                            , only : iulog
  use EMI_DataMod         , only : emi_data_list, emi_data
  use EMI_DataDimensionMod , only : emi_data_dimension_list_type
  use ColumnDataType                        , only : column_water_flux
  use VegetationDataType                    , only : vegetation_water_flux
  use EMI_Atm2LndType_Constants
  use EMI_CanopyStateType_Constants
  use EMI_ChemStateType_Constants
  use EMI_CNCarbonStateType_Constants
  use EMI_CNNitrogenStateType_Constants
  use EMI_CNNitrogenFluxType_Constants
  use EMI_CNCarbonFluxType_Constants
  use EMI_ColumnEnergyStateType_Constants
  use EMI_ColumnWaterStateType_Constants
  use EMI_ColumnWaterFluxType_Constants
  use EMI_EnergyFluxType_Constants
  use EMI_SoilHydrologyType_Constants
  use EMI_SoilStateType_Constants
  use EMI_TemperatureType_Constants
  use EMI_WaterFluxType_Constants
  use EMI_WaterStateType_Constants
  use EMI_Filter_Constants
  use EMI_ColumnType_Constants
  use EMI_Landunit_Constants
  !
  implicit none
  !
  !
  public :: EMI_Pack_WaterFluxType_at_Column_Level_for_EM
  public :: EMI_Unpack_WaterFluxType_at_Column_Level_from_EM

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_Pack_WaterFluxType_at_Column_Level_for_EM(data_list, em_stage, &
        num_filter, filter, col_wf, veg_wf)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM waterflux_vars for EM
    !
    ! !USES:
    use elm_varpar             , only : nlevgrnd
    use elm_varpar             , only : nlevsoi
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list)   , intent(in) :: data_list
    integer                , intent(in) :: em_stage
    integer                , intent(in) :: num_filter
    integer                , intent(in) :: filter(:)
    type(column_water_flux)     , intent(in) :: col_wf
    type(vegetation_water_flux) , intent(in) :: veg_wf
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fc,c,j,k
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         mflx_infl_col            => col_wf%mflx_infl            , &
         mflx_et_col              => col_wf%mflx_et              , &
         mflx_dew_col             => col_wf%mflx_dew             , &
         mflx_sub_snow_col        => col_wf%mflx_sub_snow        , &
         mflx_snowlyr_disp_col    => col_wf%mflx_snowlyr_disp    , &
         mflx_snowlyr_col         => col_wf%mflx_snowlyr         , &
         mflx_drain_col           => col_wf%mflx_drain           , &
         qflx_infl_col            => col_wf%qflx_infl            , &
         qflx_totdrain_col        => col_wf%qflx_totdrain        , &
         qflx_gross_evap_soil_col => col_wf%qflx_gross_evap_soil , &
         qflx_gross_infl_soil_col => col_wf%qflx_gross_infl_soil , &
         qflx_surf_col            => col_wf%qflx_surf            , &
         qflx_dew_grnd_col        => col_wf%qflx_dew_grnd        , &
         qflx_dew_snow_col        => col_wf%qflx_dew_snow        , &
         qflx_h2osfc2topsoi_col   => col_wf%qflx_h2osfc2topsoi   , &
         qflx_sub_snow_col        => col_wf%qflx_sub_snow        , &
         qflx_snow2topsoi_col     => col_wf%qflx_snow2topsoi     , &
         qflx_rootsoi_col         => col_wf%qflx_rootsoi         , &
         qflx_adv_col             => col_wf%qflx_adv             , &
         qflx_drain_vr_col        => col_wf%qflx_drain_vr        , &
         qflx_tran_veg_col        => col_wf%qflx_tran_veg        , &
         qflx_rootsoi_frac_patch  => veg_wf%qflx_rootsoi_frac      &
         )

    count = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit
       count = count + 1

       need_to_pack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_pack = .true.
             exit
          endif
       enddo

       if (need_to_pack) then

          select case (cur_data%id)

          case (L2E_FLUX_INFIL_MASS_FLUX)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = mflx_infl_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_VERTICAL_ET_MASS_FLUX)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = mflx_et_col(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_DEW_MASS_FLUX)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = mflx_dew_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_SNOW_SUBLIMATION_MASS_FLUX)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = mflx_sub_snow_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = mflx_snowlyr_disp_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_RESTART_SNOW_LYR_DISAPPERANCE_MASS_FLUX)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = mflx_snowlyr_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_DRAINAGE_MASS_FLUX)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = mflx_drain_col(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_INFL)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = qflx_infl_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_TOTDRAIN)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = qflx_totdrain_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_GROSS_EVAP_SOIL)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = qflx_gross_evap_soil_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_GROSS_INFL_SOIL)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = qflx_gross_infl_soil_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_SURF)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = qflx_surf_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_DEW_GRND)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = qflx_dew_grnd_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_DEW_SNOW)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = qflx_dew_snow_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_SUB_SNOW_VOL)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = qflx_h2osfc2topsoi_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_SUB_SNOW)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = qflx_sub_snow_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_H2OSFC2TOPSOI)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = qflx_h2osfc2topsoi_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_SNOW2TOPSOI)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = qflx_snow2topsoi_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_ROOTSOI)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = qflx_rootsoi_col(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_ADV)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 0, nlevgrnd
                   cur_data%data_real_2d(c,j) = qflx_adv_col(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_DRAIN_VR)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = qflx_drain_vr_col(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_TRAN_VEG)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = qflx_tran_veg_col(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_ROOTSOI_FRAC)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = qflx_rootsoi_frac_patch(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Pack_WaterFluxType_at_Column_Level_for_EM

!-----------------------------------------------------------------------
  subroutine EMI_Unpack_WaterFluxType_at_Column_Level_from_EM(data_list, em_stage, &
        num_filter, filter, col_wf)
    !
    ! !DESCRIPTION:
    ! Unpack data for ALM waterflux_vars from EM
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list)   , intent(in) :: data_list
    integer                , intent(in) :: em_stage
    integer                , intent(in) :: num_filter
    integer                , intent(in) :: filter(:)
    type(column_water_flux), intent(in) :: col_wf
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fc,c,j,k
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         mflx_snowlyr_col => col_wf%mflx_snowlyr   &
         )

    count = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit
       count = count + 1

       need_to_pack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_pack = .true.
             exit
          endif
       enddo

       if (need_to_pack) then

          select case (cur_data%id)

          case (E2L_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX)
             do fc = 1, num_filter
                c = filter(fc)
                mflx_snowlyr_col(c) = cur_data%data_real_1d(c)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Unpack_WaterFluxType_at_Column_Level_from_EM


end module EMI_WaterFluxType_ExchangeMod
