module ExternalModelATS_variables
  use, intrinsic :: iso_c_binding
  implicit none

  type :: ats_var_id_type
     integer(c_int) :: NO_VARIABLE = 0
     integer(c_int) :: BASE_POROSITY = 1
     integer(c_int) :: HYDRAULIC_CONDUCTIVITY = 2
     integer(c_int) :: CLAPP_HORNBERGER_B = 3
     integer(c_int) :: CLAPP_HORNBERGER_PSI_SAT = 4
     integer(c_int) :: RESIDUAL_SATURATION = 5
     integer(c_int) :: EFFECTIVE_POROSITY = 6
     integer(c_int) :: ROOT_FRACTION = 7
     integer(c_int) :: SURFACE_WATER_CONTENT = 8
     integer(c_int) :: WATER_CONTENT = 9
     integer(c_int) :: PRESSURE = 10
     integer(c_int) :: GROSS_SURFACE_WATER_SOURCE = 11
     integer(c_int) :: POTENTIAL_EVAPORATION = 12
     integer(c_int) :: POTENTIAL_TRANSPIRATION = 13
     integer(c_int) :: EVAPORATION = 14
     integer(c_int) :: COLUMN_TRANSPIRATION = 15
     integer(c_int) :: RUNOFF = 16
     integer(c_int) :: BASEFLOW = 17
     integer(c_int) :: ELEVATION = 18
     integer(c_int) :: TIME = 19
  end type ats_var_id_type

  type(ats_var_id_type) :: ats_var_id = ats_var_id_type()
  
end module ExternalModelATS_variables

