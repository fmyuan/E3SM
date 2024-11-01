module dynSoilAmendmentsMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle reading of the dataset that soil amendment, e.g. basalt rock powder, applications
  ! assuming that:
  !    (1) application on certain days of a year (<max. application no.), onto each soil patches (pft+cft), of each grid cell.
  !    (2) application rate in g/m2/yr, onto each soil patches (pft+cft), of each grid cell.
  !    (3) applied soil amendment is of same quality, grain size and percentage of mixed species, for individual grid cell.
  !        (i.e. NOT patch-level varied like appl. doy and rate)
  !    (4) applied soil amendment may or may not contain nutrient elements (e.g. N, P, K, ... for plant growth/development).
  !        and, if any, its weight percentage is part of species.
  
  ! !USES:
#include "shr_assert.h"
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use decompMod             , only : bounds_type, BOUNDS_LEVEL_PROC
  use dynFileMod            , only : dyn_file_type
  use dynVarTimeUninterpMod , only : dyn_var_time_uninterp_type
  use elm_varctl            , only : iulog
  use elm_varcon            , only : grlnd, namec
  use abortutils            , only : endrun
  use spmdMod               , only : masterproc, mpicom
  use LandunitType          , only : lun_pp
  use ColumnType            , only : col_pp
  use VegetationType        , only : veg_pp
  use topounit_varcon       , only : max_topounits
  
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  public :: dynsoilamendments_init     ! initialize information read from landuse.timeseries dataset
  public :: dynsoilamendments_appl     ! get soil amendment application data for the current time step, if needed
  
  ! ! PRIVATE TYPES
  type(dyn_file_type), target      :: dynsoilamendments_file ! information for the file containing transient soil amendment app data
  type(dyn_var_time_uninterp_type) :: soilamend_doy       ! soil amendment application day of year
  type(dyn_var_time_uninterp_type) :: soilamend_rate      ! soil amendment application rate [g/m2/yr] (yr, pft, lat, lon)
  type(dyn_var_time_uninterp_type) :: soilamend_grainsize ! soil amendment grain size [um] (yr, lat, lon)
  type(dyn_var_time_uninterp_type) :: soilamend_specpct   ! soil amendment chemical/mineral species weight percentage [-] (yr, species, lat, lon)
  type(dyn_var_time_uninterp_type) :: soilamend_nutrpct   ! soil amendment nutrient element weight percentage [um] (yr, n-nutrients, lat, lon)

  integer :: nspecies          ! species no. of amendment composites
  integer :: nnutrients        ! nutrient elmement no. of amendment composites (optional)

  ! Names of variables on file
  character(len=*), parameter :: soilamend_doy_varname  = 'SOIL_AMENDMENTS_DOY'
  character(len=*), parameter :: soilamend_rate_varname = 'SOIL_AMENDMENTS_RATE'
  character(len=*), parameter :: soilamend_grainsize_varname = 'SOIL_AMENDMENTS_GRAINSIZE'
  character(len=*), parameter :: soilamend_specpct_varname   = 'SOIL_AMENDMENTS_PCT'
  character(len=*), parameter :: soilamend_nutrpct_varname   = 'SOIL_AMENDMENTS_NUTRIENT_PCT'

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !---------------------------------------------------------------------------

contains
  !-----------------------------------------------------------------------
  subroutine dynsoilamendments_init(bounds, dynsoilamendments_filename, namendspec, namendnutr)
    
    ! !DESCRIPTION:
    ! Initialize dataset containing transient crop info (position it to the right time
    ! samples that bound the initial model date)
    
    ! !USES:
    use elm_varpar     , only : cft_size, natpft_size
    use ncdio_pio      , only : check_dim
    use dynTimeInfoMod , only : YEAR_POSITION_START_OF_TIMESTEP
    
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds                     ! proc-level bounds
    character(len=*)  , intent(in) :: dynsoilamendments_filename ! name of file containing soil amendment app information
    integer           , intent(in) :: namendspec                 ! max. number of species in a soil amendment
    integer,optional  , intent(in) :: namendnutr                 ! max. number of nutrients, e.g. N,P,K, ... in a soil amendment (optional)
    
    ! !LOCAL VARIABLES:
    integer :: num_points          ! number of spatial points
    integer :: mpft                ! natpft+cft
    integer :: amendappl_shape(3)  ! shape of the amendment application data
    integer :: amendgrain_shape(2)
    integer :: amendspec_shape(3)  ! shape of the amendment species data
    integer :: amendnutr_shape(3)  ! shape of the amendment nutrients data
    character(len=*), parameter :: subname = 'dynsoilamendments_init'
    !-----------------------------------------------------------------------
    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    if (masterproc) then
       write(iulog,*) 'Attempting to read dynamic soil amendment application data .....'
    end if

    !
    mpft = natpft_size + cft_size

    ! Get the year from the START of the timestep
    dynsoilamendments_file = dyn_file_type(dynsoilamendments_filename, YEAR_POSITION_START_OF_TIMESTEP)

    call check_dim(dynsoilamendments_file, 'mpft', mpft)
    call check_dim(dynsoilamendments_file, 'namendspec', namendspec)
    nspecies = namendspec

    if (present(namendnutr)) then
        call check_dim(dynsoilamendments_file, 'namendnutr', namendnutr)
        nnutrients = namendnutr
    end if

    ! read data SOIL_AMENDMENTS_DOY, SOIL_AMENDMENTS_RATE, SOIL_AMENDMENTS_GRAINSIZE, and,
    ! SOIL_AMENDMENTS_PCT corresponding to correct year

    num_points = (bounds%endg - bounds%begg + 1)

    amendappl_shape = [num_points, max_topounits, mpft]
    soilamend_doy = dyn_var_time_uninterp_type( &
        dyn_file = dynsoilamendments_file, varname=soilamend_doy_varname, &
        dim1name=grlnd, conversion_factor=1._r8, &
        do_check_sums_equal_1=.false., data_shape=amendappl_shape, &
        allow_nodata=.true.)
    soilamend_rate = dyn_var_time_uninterp_type( &
        dyn_file = dynsoilamendments_file, varname=soilamend_rate_varname, &
        dim1name=grlnd, conversion_factor=1._r8, &
        do_check_sums_equal_1=.false., data_shape=amendappl_shape, &
        allow_nodata=.true.)

    amendgrain_shape = [num_points, max_topounits]
    soilamend_grainsize = dyn_var_time_uninterp_type( &
        dyn_file = dynsoilamendments_file, varname=soilamend_grainsize_varname, &
        dim1name=grlnd, conversion_factor=1._r8, &
        do_check_sums_equal_1=.false., data_shape=amendgrain_shape, &
        allow_nodata=.true.)

    amendspec_shape = [num_points, max_topounits, namendspec]
    ! not necessary added up to 100%, but must less than 100% (will do check on-fly)
    soilamend_specpct = dyn_var_time_uninterp_type( &
        dyn_file = dynsoilamendments_file, varname=soilamend_specpct_varname, &
        dim1name=grlnd, conversion_factor=100._r8, &
        do_check_sums_equal_1=.false., data_shape=amendspec_shape, &
        allow_nodata=.true.)

    if (present(namendnutr)) then
      amendnutr_shape = [num_points, max_topounits, namendnutr]
      ! not necessary counted separatedly or added up to 100%, but must less than 100% (will do check on-fly)
      soilamend_nutrpct = dyn_var_time_uninterp_type( &
        dyn_file = dynsoilamendments_file, varname=soilamend_nutrpct_varname, &
        dim1name=grlnd, conversion_factor=100._r8, &
        do_check_sums_equal_1=.false., data_shape=amendnutr_shape, &
        allow_nodata=.true.)
    end if

  end subroutine dynsoilamendments_init

  !-----------------------------------------------------------------------
  subroutine dynsoilamendments_appl(bounds, appl_rate, appl_grainsize, appl_specpct, appl_nutrpct)

    ! !DESCRIPTION:
    ! Get soil amendment application for model time exactly (no temporal interp), when needed.

    ! !USES:
    use landunit_varcon   , only : istcrop, istsoil
    use elm_varpar        , only : cft_size, natpft_size
    use elm_varpar        , only : cft_lb, cft_ub, natpft_lb, natpft_ub
    use subgridWeightsMod , only : set_landunit_weight
    use subgridWeightsMod , only : get_landunit_weight
    use GridcellType      , only : grc_pp
    use timeinfoMod

    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    real(r8), pointer, intent(out):: appl_rate(:)                 ! col-level, summed from patches
    real(r8), pointer, intent(out):: appl_grainsize(:,:)          ! col-level, summed from patches
    real(r8), pointer, intent(out):: appl_specpct(:,:)            ! col-level, summed from patches
    real(r8), pointer, optional,intent(out):: appl_nutrpct(:,:)   ! col-level, summed from patches

    ! !LOCAL VARIABLES:
    integer               :: m,p,c,l,g,t,t2,ti,topi      ! indices
    real(r8), allocatable :: appl_doy_cur(:,:,:)         ! current timestep appl doy
    real(r8), allocatable :: appl_rate_cur(:,:,:)        ! current timestep appl rate
    real(r8), allocatable :: appl_grainsize_cur(:,:)     ! current timestep appl grain size
    real(r8), allocatable :: appl_specpct_cur(:,:,:)     ! current timestep appl specpct
    real(r8), allocatable :: appl_nutrpct_cur(:,:,:)     ! current timestep appl nutrpct
    character(len=*), parameter :: subname = 'dynsoilamendments_appl'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    call dynsoilamendments_file%time_info%set_current_year_get_year()

    allocate(appl_doy_cur(bounds%begg:bounds%endg, max_topounits, natpft_lb:cft_ub))
    allocate(appl_rate_cur(bounds%begg:bounds%endg, max_topounits, natpft_lb:cft_ub))
    allocate(appl_grainsize_cur(bounds%begg:bounds%endg, max_topounits))
    allocate(appl_specpct_cur(bounds%begg:bounds%endg, max_topounits, nspecies))
    if (present(appl_nutrpct)) &
      allocate(appl_specpct_cur(bounds%begg:bounds%endg, max_topounits, nspecies))

    ! note here data obtained is instant (not averaged over whole year)
    ! Keep in mind that what needed is an instant value for each variable
    call soilamend_doy%get_current_data(appl_doy_cur)            ! patch-level
    call soilamend_rate%get_current_data(appl_rate_cur)          ! patch-level
    call soilamend_grainsize%get_current_data(appl_grainsize_cur)! grid-level
    call soilamend_specpct%get_current_data(appl_specpct_cur)    ! grid-level
    if (present(appl_nutrpct)) &
      call soilamend_nutrpct%get_current_data(appl_nutrpct_cur)    ! grid-level

    !! output vars (already allocated memory, but may not zeroed)
    !! comment out because zeroed in col_ew_init()
    !appl_rate(bounds%begc:bounds%endc)        = 0._r8
    !appl_grainsize(bounds%begc:bounds%endc,:) = 1._r8 ! non-zero to avoid math issue. it's ok as long as rate is 0
    !appl_specpct(bounds%begc:bounds%endc, :)  = 0._r8
    !!if (present(appl_nutrpct)) &
    !  appl_nutrpct(bounds%begc:bounds%endc, :)  = 0._r8
    do p = bounds%begp, bounds%endp
       g = veg_pp%gridcell(p)
       l = veg_pp%landunit(p)
       c = veg_pp%column(p)
       t = veg_pp%topounit(p)
       topi = grc_pp%topi(g)
       ti = t - topi + 1

       if (lun_pp%itype(l) == istcrop .or. lun_pp%itype(l) == istsoil) then
          m = veg_pp%itype(p)   ! 0-based

          ! assuming application occurs at mid-day of the day of year
          ! this wouldn't make large impact, but if it is then needs more thoughts on when to add into column.
          if (appl_doy_cur(g,ti,m) == jday_mod) then
             ! weighted sum of all active pft.
             ! So if really want to patch-level, must only allow one patch per column.
            appl_rate(c) = appl_rate(c) + appl_rate_cur(g,ti,m) * veg_pp%wtcol(p)
          else
            appl_rate(c) = 0._r8
          end if

          ! assuming grain size is same for its mineral components at application
          ! but which may be changeable during dissolution and upon mineral type
          appl_grainsize(c,:) = appl_grainsize_cur(g,ti)

          appl_specpct(c,:) = appl_specpct_cur(g,ti,:)

          if (present(appl_nutrpct)) &
            appl_nutrpct(c,:) = appl_nutrpct_cur(g,ti,:)

       end if
    end do

    deallocate(appl_doy_cur)
    deallocate(appl_rate_cur)
    deallocate(appl_grainsize_cur)
    deallocate(appl_specpct_cur)
    if (present(appl_nutrpct)) &
      deallocate(appl_nutrpct_cur)

  end subroutine dynsoilamendments_appl

end module dynSoilAmendmentsMod
