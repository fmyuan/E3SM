module EMI_EnergyFluxType_ExchangeMod
  !
  use shr_kind_mod                          , only : r8 => shr_kind_r8
  use shr_log_mod                           , only : errMsg => shr_log_errMsg
  use abortutils                            , only : endrun
  use elm_varctl                            , only : iulog
  use EMI_DataMod                           , only : emi_data_list, emi_data
  use EMI_DataDimensionMod                  , only : emi_data_dimension_list_type
  use ColumnDataType                        , only : col_ef
  use EMI_Atm2LndType_Constants
  use EMI_CanopyStateType_Constants
  use EMI_ChemStateType_Constants
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
  public :: EMI_Pack_EnergyFluxType_at_Column_Level_for_EM

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_Pack_EnergyFluxType_at_Column_Level_for_EM(data_list, em_stage, &
        num_filter, filter)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM energyflux_vars for EM
    !
    ! !USES:
    use elm_varpar             , only : nlevsoi, nlevgrnd, nlevsno
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list)   , intent(in) :: data_list
    integer                , intent(in) :: em_stage
    integer                , intent(in) :: num_filter
    integer                , intent(in) :: filter(:)
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fc,c,j
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         eflx_sabg_lyr    => col_ef%eflx_sabg_lyr    , &
         eflx_hs_soil     => col_ef%eflx_hs_soil     , &
         eflx_hs_top_snow => col_ef%eflx_hs_top_snow , &
         eflx_hs_h2osfc   => col_ef%eflx_hs_h2osfc   , &
         eflx_dhsdT       => col_ef%eflx_dhsdT         &
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

          case (L2E_FLUX_ABSORBED_SOLAR_RADIATION)
             do fc = 1, num_filter
                c = filter(fc)
                do j = -nlevsno+1, 1
                   cur_data%data_real_2d(c,j) = eflx_sabg_lyr(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_SOIL_HEAT_FLUX)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = eflx_hs_soil(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_SNOW_HEAT_FLUX)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = eflx_hs_top_snow(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_H2OSFC_HEAT_FLUX)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = eflx_hs_h2osfc(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_DERIVATIVE_OF_HEAT_FLUX)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = eflx_dhsdT(c)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Pack_EnergyFluxType_at_Column_Level_for_EM


end module EMI_EnergyFluxType_ExchangeMod
