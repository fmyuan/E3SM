module C14DecayMod

  !-----------------------------------------------------------------------
  ! Module for 14-carbon flux variable update, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use elm_varpar             , only : ndecomp_cascade_transitions, nlevdecomp, ndecomp_pools
  use elm_varcon             , only : secspday
  use elm_varctl             , only : spinup_state
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CNStateType            , only : cnstate_type
  use elm_varctl             , only : nu_com
  use ColumnDataType         , only : c14_col_cs
  use VegetationDataType     , only : c14_veg_cs

  use timeInfoMod
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: C14Decay
  public:: C14BombSpike
  public:: C14_init_BombSpike

  ! !PUBLIC TYPES:
  logical, public :: use_c14_bombspike = .false.         ! do we use time-varying atmospheric C14?
  !$acc declare copyin(use_c14_bombspike)

  character(len=256), public :: atm_c14_filename = ' '   ! file name of C14 input data

  ! !PRIVATE TYPES:
  real(r8), allocatable, private :: atm_c14file_time(:)
  real(r8), allocatable, private :: atm_delta_c14(:)
  !$acc declare create(atm_c14file_time(:))
  !$acc declare create(atm_delta_c14(:)   )
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine C14Decay( num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnstate_vars )
    !
    ! !DESCRIPTION:
    ! On the radiation time step, calculate the radioactive decay of C14
    !
      !$acc routine seq
    ! !ARGUMENTS:
    integer                , intent(in)    :: num_soilc       ! number of soil columns filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnstate_type)     , intent(in)    :: cnstate_vars
    real(r8)  :: dt , days_per_year          ! radiation time step (seconds)
    integer   :: yr, mon, day, tod, offset

    !
    ! !LOCAL VARIABLES:
    integer  :: fp,j,l,p,fc,c,i
    real(r8) :: half_life
    real(r8) :: decay_const
    real(r8) :: spinup_term   ! spinup accelerated decomposition factor, used to accelerate transport as well
    !-----------------------------------------------------------------------

    associate(                                                                    &
         spinup_factor      =>    decomp_cascade_con%spinup_factor              , & ! Input:   [real(r8) (:)     ]  factor for AD spinup associated with each pool
         decomp_cpools_vr   =>    c14_col_cs%decomp_cpools_vr     , & ! Output:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
         seedc              =>    c14_col_cs%seedc                , & ! Output:  [real(r8) (:)     ]
         cpool              =>    c14_veg_cs%cpool              , & ! Output:  [real(r8) (:)     ]  (gC/m2) temporary photosynthate C pool
         xsmrpool           =>    c14_veg_cs%xsmrpool           , & ! Output:  [real(r8) (:)     ]  (gC/m2) execss maint resp C pool
         deadcrootc         =>    c14_veg_cs%deadcrootc         , & ! Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C
         deadcrootc_storage =>    c14_veg_cs%deadcrootc_storage , & ! Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C storage
         deadcrootc_xfer    =>    c14_veg_cs%deadcrootc_xfer    , & ! Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C transfer
         deadstemc          =>    c14_veg_cs%deadstemc          , & ! Output:  [real(r8) (:)     ]  (gC/m2) dead stem C
         deadstemc_storage  =>    c14_veg_cs%deadstemc_storage  , & ! Output:  [real(r8) (:)     ]  (gC/m2) dead stem C storage
         deadstemc_xfer     =>    c14_veg_cs%deadstemc_xfer     , & ! Output:  [real(r8) (:)     ]  (gC/m2) dead stem C transfer
         frootc             =>    c14_veg_cs%frootc             , & ! Output:  [real(r8) (:)     ]  (gC/m2) fine root C
         frootc_storage     =>    c14_veg_cs%frootc_storage     , & ! Output:  [real(r8) (:)     ]  (gC/m2) fine root C storage
         frootc_xfer        =>    c14_veg_cs%frootc_xfer        , & ! Output:  [real(r8) (:)     ]  (gC/m2) fine root C transfer
         gresp_storage      =>    c14_veg_cs%gresp_storage      , & ! Output:  [real(r8) (:)     ]  (gC/m2) growth respiration storage
         gresp_xfer         =>    c14_veg_cs%gresp_xfer         , & ! Output:  [real(r8) (:)     ]  (gC/m2) growth respiration transfer
         leafc              =>    c14_veg_cs%leafc              , & ! Output:  [real(r8) (:)     ]  (gC/m2) leaf C
         leafc_storage      =>    c14_veg_cs%leafc_storage      , & ! Output:  [real(r8) (:)     ]  (gC/m2) leaf C storage
         leafc_xfer         =>    c14_veg_cs%leafc_xfer         , & ! Output:  [real(r8) (:)     ]  (gC/m2) leaf C transfer
         livecrootc         =>    c14_veg_cs%livecrootc         , & ! Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C
         livecrootc_storage =>    c14_veg_cs%livecrootc_storage , & ! Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C storage
         livecrootc_xfer    =>    c14_veg_cs%livecrootc_xfer    , & ! Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C transfer
         livestemc          =>    c14_veg_cs%livestemc          , & ! Output:  [real(r8) (:)     ]  (gC/m2) live stem C
         livestemc_storage  =>    c14_veg_cs%livestemc_storage  , & ! Output:  [real(r8) (:)     ]  (gC/m2) live stem C storage
         livestemc_xfer     =>    c14_veg_cs%livestemc_xfer     , & ! Output:  [real(r8) (:)     ]  (gC/m2) live stem C transfer
         pft_ctrunc         =>    c14_veg_cs%ctrunc               & ! Output:  [real(r8) (:)     ]  (gC/m2) pft-level sink for C truncation
         )

      ! set time steps
      yr = year_curr; mon = mon_curr; day = day_curr; tod = secs_curr
      dt = dtime_mod
      days_per_year = dayspyr_mod


      half_life = 5568._r8 * secspday * days_per_year  !! libby half-life value, for comparison against ages calculated with this value
      ! half_life = 5730._r8 * secspday * days_per_year  !! recent half-life value
      decay_const = - log(0.5_r8) / half_life

      ! column loop
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         seedc(c) = seedc(c) *  (1._r8 - decay_const * dt)
      end do ! end of columns loop


         do l = 1, ndecomp_pools
            if ( spinup_state .eq. 1) then
               ! speed up radioactive decay by the same factor as decomposition so tat SOM ages prematurely in all respects
               spinup_term = spinup_factor(l)
            else
               spinup_term = 1.
            endif
            do j = 1, nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  if (spinup_term > 1._r8 .and. yr .ge. 40) then
                      decomp_cpools_vr(c,j,l) = decomp_cpools_vr(c,j,l) * (1._r8 - decay_const * &
                          (spinup_term / cnstate_vars%scalaravg_col(c,j)) * dt)
                  else
                      decomp_cpools_vr(c,j,l) = decomp_cpools_vr(c,j,l) * (1._r8 - decay_const * spinup_term * dt)
                  end if
               end do
            end do
         end do ! end of columns loop


      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         cpool(p)              = cpool(p)               * (1._r8 - decay_const * dt)
         xsmrpool(p)           = xsmrpool(p)            * (1._r8 - decay_const * dt)
         deadcrootc(p)         = deadcrootc(p)          * (1._r8 - decay_const * dt)
         deadcrootc_storage(p) = deadcrootc_storage(p)  * (1._r8 - decay_const * dt)
         deadcrootc_xfer(p)    = deadcrootc_xfer(p)     * (1._r8 - decay_const * dt)
         deadstemc(p)          = deadstemc(p)           * (1._r8 - decay_const * dt)
         deadstemc_storage(p)  = deadstemc_storage(p)   * (1._r8 - decay_const * dt)
         deadstemc_xfer(p)     = deadstemc_xfer(p)      * (1._r8 - decay_const * dt)
         frootc(p)             = frootc(p)              * (1._r8 - decay_const * dt)
         frootc_storage(p)     = frootc_storage(p)      * (1._r8 - decay_const * dt)
         frootc_xfer(p)        = frootc_xfer(p)         * (1._r8 - decay_const * dt)
         gresp_storage(p)      = gresp_storage(p)       * (1._r8 - decay_const * dt)
         gresp_xfer(p)         = gresp_xfer(p)          * (1._r8 - decay_const * dt)
         leafc(p)              = leafc(p)               * (1._r8 - decay_const * dt)
         leafc_storage(p)      = leafc_storage(p)       * (1._r8 - decay_const * dt)
         leafc_xfer(p)         = leafc_xfer(p)          * (1._r8 - decay_const * dt)
         livecrootc(p)         = livecrootc(p)          * (1._r8 - decay_const * dt)
         livecrootc_storage(p) = livecrootc_storage(p)  * (1._r8 - decay_const * dt)
         livecrootc_xfer(p)    = livecrootc_xfer(p)     * (1._r8 - decay_const * dt)
         livestemc(p)          = livestemc(p)           * (1._r8 - decay_const * dt)
         livestemc_storage(p)  = livestemc_storage(p)   * (1._r8 - decay_const * dt)
         livestemc_xfer(p)     = livestemc_xfer(p)      * (1._r8 - decay_const * dt)
         pft_ctrunc(p)         = pft_ctrunc(p)          * (1._r8 - decay_const * dt)
      end do

    end associate

  end subroutine C14Decay

  !-----------------------------------------------------------------------
   subroutine C14BombSpike( num_soilp, filter_soilp, &
        cnstate_vars )
     !
     ! !DESCRIPTION:
     ! for transient pulse simulation, impose a simplified bomb spike
     !
     ! !USES:
      !$acc routine seq
     use elm_varcon       , only : c14ratio, secspday
     !
     ! !ARGUMENTS:
     integer, intent(in) :: num_soilp       ! number of soil patches in filter
     integer, intent(in) :: filter_soilp(:) ! filter for soil patches
     type(cnstate_type), intent(inout) :: cnstate_vars
     integer  :: yr, mon, day, tod
     real(r8)  :: days_per_year ! days per year

     !
     ! !LOCAL VARIABLES:

     real(r8) :: dateyear
     real(r8) :: delc14o2_atm
     integer  :: fp, p, nt
     integer  :: ind_below
     integer  :: ntim_atm_ts
     real(r8) :: twt_1, twt_2  ! weighting fractions for interpolating
     !-----------------------------------------------------------------------

     if ( use_c14_bombspike ) then

        ! get current date
        yr = year_curr; mon = mon_curr; day = day_curr; tod = secs_curr
        days_per_year = dayspyr_mod


        dateyear = real(yr) + real(mon)/12._r8 + real(day)/days_per_year + real(tod)/(secspday*days_per_year)

        ! find points in atm timeseries to interpolate between
        ntim_atm_ts = size(atm_c14file_time)
        ind_below = 0
        do nt = 1, ntim_atm_ts
           if (dateyear >= atm_c14file_time(nt) ) then
              ind_below = ind_below+1
           endif
        end do

        ! interpolate between nearest two points in atm c14 timeseries
        if (ind_below .eq. 0 ) then
           delc14o2_atm = atm_delta_c14(1)
        elseif (ind_below .eq. ntim_atm_ts ) then
           delc14o2_atm = atm_delta_c14(ntim_atm_ts)
        else
           twt_2 = min(1._r8, max(0._r8,(dateyear-atm_c14file_time(ind_below)) &
                / (atm_c14file_time(ind_below+1)-atm_c14file_time(ind_below))))
           twt_1 = 1._r8 - twt_2
           delc14o2_atm = atm_delta_c14(ind_below) * twt_1 +  atm_delta_c14(ind_below+1) * twt_2
        endif

        ! change delta units to ratio, put on pft loop
        do fp = 1,num_soilp
           p = filter_soilp(fp)
           cnstate_vars%rc14_atm_patch(p) = (delc14o2_atm * 1.e-3_r8 + 1._r8) * c14ratio
        end do

     else

        ! for constant 14c concentration
        do fp = 1,num_soilp
           p = filter_soilp(fp)
           cnstate_vars%rc14_atm_patch(p) = c14ratio
        end do

     endif

   end subroutine C14BombSpike

  !-----------------------------------------------------------------------
  subroutine C14_init_BombSpike()
    !
    ! !DESCRIPTION:
    ! read netcdf file containing a timeseries of atmospheric delta C14 values; save in module-level array
    !
    ! !USES:
    use ncdio_pio
    use fileutils   , only : getfil
    use abortutils  , only : endrun
    use elm_varctl  , only : iulog
    use spmdMod     , only : masterproc
    use shr_log_mod , only : errMsg => shr_log_errMsg
    !
    ! !LOCAL VARIABLES:
    character(len=256) :: locfn           ! local file name
    type(file_desc_t)  :: ncid            ! netcdf id
    integer :: dimid,varid                ! input netCDF id's
    integer :: ntim                       ! number of input data time samples
    integer :: t
    !-----------------------------------------------------------------------

    if ( use_c14_bombspike ) then

       if ( masterproc ) then
          write(iulog, *) 'C14_init_BombSpike: preparing to open file:'
          write(iulog, *) trim(locfn)
       endif

       call getfil(atm_c14_filename, locfn, 0)

       call ncd_pio_openfile (ncid, trim(locfn), 0)

       call ncd_inqdlen(ncid,dimid,ntim,'time')

       !! allocate arrays based on size of netcdf timeseries
       allocate(atm_c14file_time(ntim))
       allocate(atm_delta_c14(ntim))

       call ncd_io(ncid=ncid, varname='time', flag='read', data=atm_c14file_time)

       call ncd_io(ncid=ncid, varname='atm_delta_c14', flag='read', data=atm_delta_c14)

       call ncd_pio_closefile(ncid)

       ! check to make sure that time dimension is well behaved
       do t = 2, ntim
          if ( atm_c14file_time(t) - atm_c14file_time(t-1) <= 0._r8 ) then
             write(iulog, *) 'C14_init_BombSpike: error.  time axis must be monotonically increasing'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          endif
       end do
    endif

  end subroutine C14_init_BombSpike

end module C14DecayMod
