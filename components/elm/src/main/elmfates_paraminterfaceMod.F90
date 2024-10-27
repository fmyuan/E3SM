module CLMFatesParamInterfaceMod
  ! NOTE(bja, 2017-01) this code can not go into the main clm-fates
  ! interface module because of circular dependancies with pftvarcon.

  use FatesGlobals, only : fates_log
  use shr_kind_mod, only : r8 => shr_kind_r8
  
  implicit none

  ! NOTE(bja, 2017-01) these methods can NOT be part of the clmi-fates
  ! nterface type because they are called before the instance is
  ! initialized.
  public :: FatesReadParameters
  public :: FatesReadPFTs
  private :: ParametersFromNetCDF
  private :: SetParameterDimensions
  private :: GetUsedDimensionSizes

  logical :: DEBUG  = .false.

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  
contains
  
 !-----------------------------------------------------------------------
 subroutine FatesReadParameters()

   use elm_varctl, only : use_fates, paramfile, fates_paramfile
   use spmdMod, only : masterproc

   use FatesParametersInterface, only : fates_parameters_type

   use EDParamsMod, only : FatesRegisterParams, FatesReceiveParams
   use SFParamsMod, only : SpitFireRegisterParams, SpitFireReceiveParams
   use PRTInitParamsFATESMod, only : PRTRegisterParams, PRTReceiveParams
   use FatesSynchronizedParamsMod, only : FatesSynchronizedParamsInst

   implicit none

   character(len=32)  :: subname = 'FatesReadParameters'
   class(fates_parameters_type), allocatable :: fates_params
   logical :: is_host_file

   if (use_fates) then
      if (masterproc) then
         write(fates_log(), *) 'clmfates_parameterinterfaceMod.F90::'//trim(subname)//' :: CLM reading ED/FATES '//' parameters '
      end if

      allocate(fates_params)
      call fates_params%Init()
      call FatesRegisterParams(fates_params)
      call SpitFireRegisterParams(fates_params)
      call PRTRegisterParams(fates_params)
      call FatesSynchronizedParamsInst%RegisterParams(fates_params)

      is_host_file = .false.
      call ParametersFromNetCDF(fates_paramfile, is_host_file, fates_params)

      is_host_file = .true.
      call ParametersFromNetCDF(paramfile, is_host_file, fates_params)

      call FatesReceiveParams(fates_params)
      call SpitFireReceiveParams(fates_params)
      call PRTReceiveParams(fates_params)
      call FatesSynchronizedParamsInst%ReceiveParams(fates_params)

      call fates_params%Destroy()
      deallocate(fates_params)
   end if

 end subroutine FatesReadParameters

 !-----------------------------------------------------------------------
 subroutine FatesReadPFTs()

   use elm_varctl, only : use_fates, paramfile, fates_paramfile
   use spmdMod, only : masterproc

   use FatesParametersInterface, only : fates_parameters_type
   use EDPftvarcon , only : EDPftvarcon_inst

   use fileutils  , only : getfil
   use ncdio_pio  , only : file_desc_t, ncd_pio_closefile, ncd_pio_openfile

   implicit none

   character(len=32)  :: subname = 'FatesReadPFTs'
   class(fates_parameters_type), allocatable :: fates_params
   logical :: is_host_file

   character(len=256) :: locfn ! local file name
   type(file_desc_t)  :: ncid  ! pio netCDF file id

   if (use_fates) then
      if (masterproc) then
         write(fates_log(), *) 'clmfates_interfaceMod.F90::'//trim(subname)//' :: CLM reading ED/FATES '//' PFTs '
      end if

      allocate(fates_params)
      call fates_params%Init()
      call EDPftvarcon_inst%Init()

      call EDPftvarcon_inst%Register(fates_params)

      is_host_file = .false.
      
      call ParametersFromNetCDF(fates_paramfile, is_host_file, fates_params)

      is_host_file = .true.
      write(fates_log(), *) 'JD  elmfates_paraminterfaceMod.F90 L111'
      write(fates_log(), *) 'host model par file: ', paramfile
      call ParametersFromNetCDF(paramfile, is_host_file, fates_params)

      write(fates_log(), *) 'JD L115 call EDPftvarcon_inst%Receive'
      call EDPftvarcon_inst%Receive(fates_params)
      write(fates_log(), *) 'JD  done EDPftvarcon_inst%Receive'
      
      call fates_params%Destroy()
      deallocate(fates_params)
   end if

 end subroutine FatesReadPFTs

 !-----------------------------------------------------------------------
 subroutine SetParameterDimensions(ncid, is_host_file, fates_params)
   ! Get the list of dimensions used by the fates parameters,
   ! retreive them from the parameter file, then give the information
   ! back to fates.
   use FatesParametersInterface, only : fates_parameters_type, param_string_length, max_dimensions, max_used_dimensions
   use ncdio_pio  , only : file_desc_t

   implicit none

   type(file_desc_t), intent(inout) :: ncid
   logical, intent(in) :: is_host_file
   class(fates_parameters_type), intent(inout) :: fates_params

   integer :: num_used_dimensions
   character(len=param_string_length) :: used_dimension_names(max_used_dimensions)
   integer :: used_dimension_sizes(max_used_dimensions)

   call fates_params%GetUsedDimensions(is_host_file, num_used_dimensions, used_dimension_names)

   call GetUsedDimensionSizes(ncid, num_used_dimensions, used_dimension_names, used_dimension_sizes)

   call fates_params%SetDimensionSizes(is_host_file, num_used_dimensions, used_dimension_names, used_dimension_sizes)

 end subroutine SetParameterDimensions

 !-----------------------------------------------------------------------
 subroutine GetUsedDimensionSizes(ncid, num_used_dimensions, dimension_names, dimension_sizes)

   use ncdio_pio  , only : ncd_inqdid, ncd_inqdlen
   use FatesParametersInterface, only : param_string_length
   use ncdio_pio, only : file_desc_t


   implicit none

   type(file_desc_t), intent(inout) :: ncid
   integer, intent(in) :: num_used_dimensions
   character(len=param_string_length), intent(in) :: dimension_names(:)
   integer, intent(out) :: dimension_sizes(:)

   integer :: d, max_dim_size, num_dims
   integer :: dim_len, dim_id

   dimension_sizes(:) = 0
   max_dim_size = 0

   do d = 1, num_used_dimensions
      call ncd_inqdid(ncid, dimension_names(d), dim_id)
      call ncd_inqdlen(ncid, dim_id, dim_len)
      dimension_sizes(d) = dim_len
      !write(*, *) '--> ', trim(dimension_names(d)), ' setting size ', dimension_sizes(d)
   end do

 end subroutine GetUsedDimensionSizes

 !-----------------------------------------------------------------------
 subroutine ParametersFromNetCDF(filename, is_host_file, fates_params)

   use abortutils, only : endrun
   use fileutils  , only : getfil
   use ncdio_pio  , only : file_desc_t, ncd_pio_closefile, ncd_pio_openfile
   use paramUtilMod, only : readNcdio

   use FatesParametersInterface, only : fates_parameters_type
   use FatesParametersInterface, only : param_string_length, max_dimensions, max_used_dimensions
   use FatesParametersInterface, only : dimension_shape_scalar, dimension_shape_1d, dimension_shape_2d

   implicit none

   character(len=*), intent(in) :: filename
   logical, intent(in) :: is_host_file
   class(fates_parameters_type), intent(inout) :: fates_params

   character(len=40)  :: subname = 'clmfates_interface::ReadParameters'
   character(len=256) :: locfn ! local file name
   type(file_desc_t)  :: ncid  ! pio netCDF file id
   integer            :: dimid ! netCDF dimension id
   integer :: i, num_params, dimension_shape
   integer :: max_dim_size
   real(r8), allocatable :: data(:, :)
   character(len=param_string_length) :: name
   integer :: dimension_sizes(max_dimensions)
   character(len=param_string_length) :: dimension_names(max_dimensions)
   integer :: size_dim_1, size_dim_2
   logical :: is_host_param

   call getfil (filename, locfn, 0)
   write(fates_log(),*) 'JD elmfates_paraminterfaceMod.F90  L 211'
   write(fates_log(),*) 'file name:', filename, ', locfn: ', locfn
   
   call ncd_pio_openfile (ncid, trim(locfn), 0)
   
   write(fates_log(),*) 'call SetParameterDimensions', 'is_host_file: ', is_host_file
   call SetParameterDimensions(ncid, is_host_file, fates_params)
   write(fates_log(),*) 'done call SetParameterDimensions'
    
   max_dim_size = fates_params%GetMaxDimensionSize()
   allocate(data(max_dim_size, max_dim_size))

   num_params = fates_params%num_params()
   write(fates_log(),*) 'num_params: ', num_params   
   do i = 1, num_params
      
      call fates_params%GetMetaData(i, name, dimension_shape, dimension_sizes, dimension_names, is_host_param)
      write(fates_log(),*) 'from call fates_params%GetMetaData'
      write(fates_log(),*) 'i: ', i
      write(fates_log(),*) 'name: ', name
      write(fates_log(),*) 'dimension_shape: ', dimension_shape
      write(fates_log(),*) 'dimension_sizes: ', dimension_sizes      
      write(fates_log(),*) 'dimension_names: ', dimension_names
      write(fates_log(),*) 'is_host_param: ', is_host_param
      write(fates_log(),*) '=================='               
                
      if (is_host_file .eqv. is_host_param) then
         select case(dimension_shape)
         case(dimension_shape_scalar)
            size_dim_1 = 1
            size_dim_2 = 1
         case(dimension_shape_1d)
            size_dim_1 = dimension_sizes(1)
            size_dim_2 = 1
         case(dimension_shape_2d)
            size_dim_1 = dimension_sizes(1)
            size_dim_2 = dimension_sizes(2)
         case default
            write(fates_log(),*) 'dimension shape:',dimension_shape
            call endrun(msg='unsupported number of dimensions reading parameters.')
         end select
         if(DEBUG) then
            write(fates_log(), *) 'clmfates_interfaceMod.F90:: reading '//trim(name)
         end if
         call readNcdio(ncid, name, dimension_shape, dimension_names, subname, data(1:size_dim_1, 1:size_dim_2))
         call fates_params%SetData(i, data(1:size_dim_1, 1:size_dim_2))
      end if
   end do
   write(fates_log(),*) 'JD elmfates_paraminterfaceMod.F90  L 258 (end)'
   deallocate(data)
   call ncd_pio_closefile(ncid)
 end subroutine ParametersFromNetCDF
 !-----------------------------------------------------------------------

end module CLMFatesParamInterfaceMod
