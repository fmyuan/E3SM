module EMI_GrndWaterStateType_ExchangeMod
  !
  use shr_kind_mod                          , only : r8 => shr_kind_r8
  use shr_log_mod                           , only : errMsg => shr_log_errMsg
  use abortutils                            , only : endrun
  use elm_varctl                            , only : iulog
  use EMI_DataMod                           , only : emi_data_list, emi_data
  use EMI_DataDimensionMod                  , only : emi_data_dimension_list_type
  use ColumnDataType                        , only : column_water_state
  use EMI_GrndWaterStateType_Constants
  !
  implicit none
  !
  !
  public :: EMI_Pack_GrndWaterStateType_at_Column_Level_for_EM

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_Pack_GrndWaterStateType_at_Column_Level_for_EM(data_list, em_stage, &
        num_filter, filter, col_ws)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM waterstate_vars for EM
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
    type(column_water_state) , intent(in) :: col_ws
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fc,c,j
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         frac_h2osfc   => col_ws%frac_h2osfc   , &
         finundated    => col_ws%finundated    , &
         h2osfc        => col_ws%h2osfc        , &
         h2osno        => col_ws%h2osno        , &
         h2osoi_liq    => col_ws%h2osoi_liq    , &
         h2osoi_ice    => col_ws%h2osoi_ice    , &
         h2osoi_vol    => col_ws%h2osoi_vol    , &
         frac_sno_eff  => col_ws%frac_sno_eff    &
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

          case (L2E_STATE_FRAC_H2OSFC)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = frac_h2osfc(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_FRAC_INUNDATED)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = finundated(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSOI_LIQ_NLEVSNOW)
             do fc = 1, num_filter
                c = filter(fc)
                do j = -nlevsno + 1, 0
                   cur_data%data_real_2d(c,j) = h2osoi_liq(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSOI_ICE_NLEVSNOW)
             do fc = 1, num_filter
                c = filter(fc)
                do j = -nlevsno + 1, 0
                   cur_data%data_real_2d(c,j) = h2osoi_ice(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSNOW)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = h2osno(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSFC)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = h2osfc(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_FRAC_SNOW_EFFECTIVE)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = frac_sno_eff(c)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Pack_GrndWaterStateType_at_Column_Level_for_EM

!-----------------------------------------------------------------------


end module EMI_GrndWaterStateType_ExchangeMod
