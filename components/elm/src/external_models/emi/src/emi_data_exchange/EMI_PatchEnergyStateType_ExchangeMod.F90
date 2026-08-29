module EMI_PatchEnergyStateType_ExchangeMod
  !
  use shr_kind_mod                          , only : r8 => shr_kind_r8
  use shr_log_mod                           , only : errMsg => shr_log_errMsg
  use abortutils                            , only : endrun
  use elm_varctl                            , only : iulog
  use EMI_DataMod                           , only : emi_data_list, emi_data
  use EMI_DataDimensionMod                  , only : emi_data_dimension_list_type
  use VegetationDataType                    , only : veg_es
  use EMI_GrndEnergyStateType_Constants
  !
  implicit none
  !
  !
  public :: EMI_Pack_EnergyStateType_at_Patch_Level_for_EM

contains
  

!-----------------------------------------------------------------------
  subroutine EMI_Pack_EnergyStateType_at_Patch_Level_for_EM(data_list, em_stage, &
        num_filter, filter)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM temperature_vars for EM
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
    integer                             :: fp,p,j
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         t_veg => veg_es%t_veg   &
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

          case (L2E_STATE_TVEG)
             do fp = 1, num_filter
                p = filter(fp)
                cur_data%data_real_1d(p) = t_veg(p)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Pack_EnergyStateType_at_Patch_Level_for_EM

!-----------------------------------------------------------------------

end module EMI_PatchEnergyStateType_ExchangeMod
