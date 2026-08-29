module EMI_ColumnWaterStateType_ExchangeMod
  !
  use shr_kind_mod                          , only : r8 => shr_kind_r8
  use shr_log_mod                           , only : errMsg => shr_log_errMsg
  use abortutils                            , only : endrun
  use elm_varctl                            , only : iulog
  use EMI_DataMod                           , only : emi_data_list, emi_data
  use EMI_DataDimensionMod                  , only : emi_data_dimension_list_type
  use ColumnDataType       , only : column_water_state
  use EMI_ColumnWaterStateType_Constants
  !
  implicit none
  !
  !
  public :: EMI_Pack_ColumnWaterStateType_at_Column_Level_for_EM
  public :: EMI_Unpack_ColumnWaterStateType_at_Column_Level_from_EM

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_Pack_ColumnWaterStateType_at_Column_Level_for_EM(data_list, em_stage, &
        num_filter, filter, col_ws)
    !
    ! !DESCRIPTION:
    ! Pack data from ELM col_ws for EM
    !
    ! !USES:
    use elm_varpar             , only : nlevgrnd
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list)     , intent(in) :: data_list
    integer                  , intent(in) :: em_stage
    integer                  , intent(in) :: num_filter
    integer                  , intent(in) :: filter(:)
    type(column_water_state) , intent(in) :: col_ws
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fc,c,j,k
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         h2osoi_liqvol => col_ws%h2osoi_liqvol , &
         h2osoi_icevol => col_ws%h2osoi_icevol , &
         h2osfc        => col_ws%h2osfc        , &
         salinity      => col_ws%salinity      , &
         nitrate_tide  => col_ws%nitrate_tide  , &
         h2osfc_tide   => col_ws%h2osfc_tide     &
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

          case (L2E_STATE_SOIL_LIQ_VOL_COL)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = h2osoi_liqvol(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_SOIL_ICE_VOL_COL)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = h2osoi_icevol(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSFC_COL)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = h2osfc(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_SALINITY_COL)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = salinity(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_TIDE_NITRATE_COL)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = nitrate_tide(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSFC_TIDE_COL)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = h2osfc_tide(c)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Pack_ColumnWaterStateType_at_Column_Level_for_EM

!-----------------------------------------------------------------------
  subroutine EMI_Unpack_ColumnWaterStateType_at_Column_Level_from_EM(data_list, em_stage, &
        num_filter, filter, col_ws)
    !
    ! !DESCRIPTION:
    ! Unpack data from ELM col_ws for EM
    !
    ! !USES:
    use elm_varpar             , only : nlevgrnd
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list)     , intent(in) :: data_list
    integer                  , intent(in) :: em_stage
    integer                  , intent(in) :: num_filter
    integer                  , intent(in) :: filter(:)
    type(column_water_state) , intent(in) :: col_ws
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fc,c,j,k
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(&
         h2osoi_liq => col_ws%h2osoi_liq , &
         h2osoi_ice => col_ws%h2osoi_ice   &
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

          case (E2L_STATE_H2OSOI_LIQ)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   h2osoi_liq(c,j) = cur_data%data_real_2d(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (E2L_STATE_H2OSOI_ICE)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   h2osoi_ice(c,j) = cur_data%data_real_2d(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Unpack_ColumnWaterStateType_at_Column_Level_from_EM


end module EMI_ColumnWaterStateType_ExchangeMod
