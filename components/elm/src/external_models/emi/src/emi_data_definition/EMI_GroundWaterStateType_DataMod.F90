module EMI_GroundWaterStateType_DataMod
  !
  use EMI_GroundWaterStateType_Constants
  !
  implicit none
  !
  public :: EMI_GroundWaterStateType_DataInfoByID

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_GroundWaterStateType_DataInfoByID(data_id, id_val, name_val, long_name_val,&
        units_val, is_int_type, is_real_type, ndim, &
        dim1_beg_name, dim1_end_name, dim2_beg_name, dim2_end_name, &
        dim3_beg_name, dim3_end_name, dim4_beg_name, dim4_end_name, &
        data_found)
    !
    ! !DESCRIPTION:
    ! Defines information of data exchanged between ELM and EM
    !
    ! !USES: 
    use EMI_DataDimensionMod
    implicit none
    !
    ! !ARGUMENTS:
    integer            , intent(in)  :: data_id
    integer            , intent(out) :: id_val
    character (len=32) , intent(out) :: name_val
    character (len=128), intent(out) :: long_name_val
    character (len=32) , intent(out) :: units_val
    logical            , intent(out) :: is_int_type
    logical            , intent(out) :: is_real_type
    integer            , intent(out) :: ndim
    character (len=32) , intent(out) :: dim1_beg_name
    character (len=32) , intent(out) :: dim1_end_name
    character (len=32) , intent(out) :: dim2_beg_name
    character (len=32) , intent(out) :: dim2_end_name
    character (len=32) , intent(out) :: dim3_beg_name
    character (len=32) , intent(out) :: dim3_end_name
    character (len=32) , intent(out) :: dim4_beg_name
    character (len=32) , intent(out) :: dim4_end_name
    logical            , intent(inout) :: data_found

    is_int_type    = .false.
    is_real_type   = .false.
    dim1_beg_name  = ''
    dim2_beg_name  = ''
    dim3_beg_name  = ''
    dim4_beg_name  = ''
    dim1_end_name  = ''
    dim2_end_name  = ''
    dim3_end_name  = ''
    dim4_end_name  = ''

    select case(data_id)

    case(L2E_STATE_FRAC_H2OSFC)
       id_val         =  L2E_STATE_FRAC_H2OSFC
       name_val       =  'Fraction of standing water'
       long_name_val  =  'Fraction of standing water: ELM to EM'
       units_val      =  '[-]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_STATE_FRAC_INUNDATED)
       id_val         =  L2E_STATE_FRAC_INUNDATED
       name_val       =  'Fraction inundated'
       long_name_val  =  'Fraction inundated: ELM to EM'
       units_val      =  '[-]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_STATE_H2OSOI_LIQ_NLEVSNOW)
       id_val         =  L2E_STATE_H2OSOI_LIQ_NLEVSNOW
       name_val       =  'Soil liquid water'
       long_name_val  =  'Soil liquid water: ELM to EM'
       units_val      =  '[kg/m2]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_nlevsno_plus_one
       dim2_end_name  =  dimname_zero
       data_found   =  .true.

    case(L2E_STATE_H2OSOI_ICE_NLEVSNOW)
       id_val         =  L2E_STATE_H2OSOI_ICE_NLEVSNOW
       name_val       =  'Soil ice water'
       long_name_val  =  'Soil ice water: ELM to EM'
       units_val      =  '[kg/m2]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_nlevsno_plus_one
       dim2_end_name  =  dimname_zero
       data_found   =  .true.

    case(L2E_STATE_H2OSNOW)
       id_val         =  L2E_STATE_H2OSNOW
       name_val       =  'Snow water'
       long_name_val  =  'Snow water: ELM to EM'
       units_val      =  '[mm]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_STATE_H2OSFC)
       id_val         =  L2E_STATE_H2OSFC
       name_val       =  'Standing surface water'
       long_name_val  =  'Standing surface water: ELM to EM'
       units_val      =  '[mm]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_STATE_FRAC_SNOW_EFFECTIVE)
       id_val         =  L2E_STATE_FRAC_SNOW_EFFECTIVE
       name_val       =  'Frac. of grnd covered with snow'
       long_name_val  =  'Frac. of grnd covered with snow: ELM to EM'
       units_val      =  '[-]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    end select
    
  end subroutine EMI_GroundWaterStateType_DataInfoByID
    
end module EMI_GroundWaterStateType_DataMod
