module EMI_ColumnWaterStateType_Constants
  !
  implicit none
  private

  !
  integer, parameter, public :: L2E_STATE_SOIL_LIQ_VOL_COL  = 0401
  integer, parameter, public :: L2E_STATE_SOIL_ICE_VOL_COL  = 0402
  integer, parameter, public :: L2E_STATE_H2OSFC_COL        = 0403
  integer, parameter, public :: L2E_STATE_SALINITY_COL      = 0404
  integer, parameter, public :: L2E_STATE_TIDE_NITRATE_COL  = 0405
  integer, parameter, public :: L2E_STATE_H2OSFC_TIDE_COL   = 0406

  integer, parameter, public :: L2E_STATE_H2OSOI_LIQ_NLEVGRND = 0411
  integer, parameter, public :: L2E_STATE_H2OSOI_ICE_NLEVGRND = 0412
  integer, parameter, public :: L2E_STATE_AIR_VOL_NLEVGRND    = 0413
  integer, parameter, public :: L2E_STATE_RHO_VAP_NLEVGRND    = 0414
  integer, parameter, public :: L2E_STATE_RHVAP_SOI_NLEVGRND  = 0415
  integer, parameter, public :: L2E_STATE_SOIL_MATRIC_POTENTIAL_NLEGRND  = 0416
  integer, parameter, public :: L2E_STATE_H2OSOI_LIQ_NLEVSNOW            = 0417
  integer, parameter, public :: L2E_STATE_H2OSOI_ICE_NLEVSNOW            = 0418

  integer, parameter, public :: E2L_STATE_H2OSOI_LIQ        = 0421
  integer, parameter, public :: E2L_STATE_H2OSOI_ICE        = 0422

end module EMI_ColumnWaterStateType_Constants
