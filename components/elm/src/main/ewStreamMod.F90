module ewStreamMod

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION: 
  ! Contains methods for reading from dynamic gridded files
  !  - application rate (kg rock m-2 year-1)
  !  - weight percentage of minerals in rock (kg mineral kg-1 rock)
  !  - weight percentage of phosphorus content in rock (gP kg-1 rock)
  !  - grain size (um diameter)
  !  - soil pH
  ! Also includes functions for dynamic file handling and interpolation.

  ! Below will be from parameter files
  ! - reaction equation stoichiometry for each of the minerals

  use shr_kind_mod, only: r8 => shr_kind_r8, CL => shr_kind_cl, CX => shr_kind_CXX
  use shr_strdata_mod
  use shr_stream_mod
  use shr_string_mod
  use shr_sys_mod
  use shr_mct_mod
  use mct_mod
  use spmdMod     , only: mpicom, masterproc, comp_id, iam
  use elm_varctl  , only: iulog
  use shr_sys_mod , only : shr_sys_flush
  use controlMod  , only: NLFilename
  use elm_varpar  , only: nminerals
  use abortutils  , only: endrun
  use fileutils   , only: getavu, relavu
  use decompMod   , only: bounds_type, ldecomp, gsmap_lnd_gdc2glo 
  use domainMod   , only: ldomain
  use ColumnDataType, only: col_ew

  ! !PUBLIC TYPES:
  implicit none
  private
  save

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: ew_init   ! initialize enhanced weathering
  public :: ew_interp ! interpolate enhanced weathering data

  ! ! PRIVATE TYPES
  type(shr_strdata_type)  :: sdat_app       ! input data stream
  type(shr_strdata_type)  :: sdat_min       ! input data stream
  type(shr_strdata_type)  :: sdat_pho       ! input data stream
  type(shr_strdata_type)  :: sdat_gra       ! input data stream
  type(shr_strdata_type)  :: sdat_sph       ! input data stream
  type(shr_strdata_type)  :: sdat_pctcrop   ! input data stream

  ! local variables
  integer :: stream_year_first_ew       ! first year in stream to use
  integer :: stream_year_last_ew        ! last year in stream to use
  integer :: model_year_align_ew        ! align stream_year_first_ew with
  integer :: doy_application_ew         ! date of application of rock powder in each year (1-365)
  !----------------------------------------------------------------------- 

contains
  !-----------------------------------------------------------------------
  subroutine ew_init(bounds)
   !
   ! Initialize enhanced weathering data streams
   !
   ! Uses:
   use elm_varctl       , only : inst_name
   use elm_time_manager , only : get_calendar
   use ncdio_pio        , only : pio_subsystem
   use shr_pio_mod      , only : shr_pio_getiotype
   use shr_nl_mod       , only : shr_nl_find_group_name
   use ndepStreamMod    , only : elm_domain_mct
   use shr_log_mod      , only : errMsg => shr_log_errMsg
   use shr_string_mod   , only : shr_string_listCreateField
   !
   ! arguments
   implicit none
   type(bounds_type), intent(in) :: bounds
   !

   integer            :: nu_nml    ! unit for namelist file
   integer            :: nml_error ! namelist i/o error flag
   type(mct_ggrid)    :: dom_elm   ! domain information 
   character(len=CL)  :: stream_fldFilename_ew
   character(len=CL)  :: ew_mapalgo = 'bilinear'
   character(*), parameter :: shr_strdata_unset = 'NOT_SET'
   character(*), parameter :: subName = "('ew_dyn_init')"
   character(*), parameter :: F00 = "('(ew_dyn_init) ',4a)"
   character(CX)      :: fldList            ! field string for mineral dimensions
   !-----------------------------------------------------------------------

   namelist /ew_streams/         &
        stream_year_first_ew,    &
        stream_year_last_ew,     &
        model_year_align_ew,     &
        ew_mapalgo,              &
        doy_application_ew,      &
        stream_fldFilename_ew

   ! Default values for namelist
   stream_year_first_ew     = 1          ! first year in stream to use
   stream_year_last_ew      = 1          ! last  year in stream to use
   model_year_align_ew      = 1          ! align stream_year_first_ew with this model year
   doy_application_ew       = 1          ! date of application of rock powder in each year (1-365)
   stream_fldFileName_ew    = ' '

   ! Read ew_streams namelist
   if (masterproc) then
      nu_nml = getavu()
      open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call shr_nl_find_group_name(nu_nml, 'ew_streams', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=ew_streams, iostat=nml_error)
         if (nml_error /= 0) then
            call endrun(msg=' ERROR reading ew_streams namelist'//errMsg(__FILE__, __LINE__))
         end if
      end if
      close(nu_nml)
      call relavu( nu_nml )
   endif

   call shr_mpi_bcast(stream_year_first_ew, mpicom)
   call shr_mpi_bcast(stream_year_last_ew, mpicom)
   call shr_mpi_bcast(model_year_align_ew, mpicom)
   call shr_mpi_bcast(doy_application_ew, mpicom)
   call shr_mpi_bcast(stream_fldFileName_ew, mpicom)

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'EW application rate stream settings:'
      write(iulog,*) '  stream_year_first_ew  = ',stream_year_first_ew
      write(iulog,*) '  stream_year_last_ew   = ',stream_year_last_ew
      write(iulog,*) '  model_year_align_ew   = ',model_year_align_ew
      write(iulog,*) '  doy_application_ew    = ',doy_application_ew
      write(iulog,*) '  stream_fldFileName_ew = ',trim(stream_fldFileName_ew)
      write(iulog,*) ' '
   endif

   call elm_domain_mct (bounds, dom_elm)

   ! enhanced weathering application rate
   call shr_strdata_create(sdat_app,name="ew_appdyn",      &
        pio_subsystem=pio_subsystem,               &
        pio_iotype=shr_pio_getiotype(inst_name),   &
        mpicom=mpicom, compid=comp_id,             &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_elm,    &
        nxg=ldomain%ni, nyg=ldomain%nj,            &
        yearFirst=stream_year_first_ew,            &
        yearLast=stream_year_last_ew,              &
        yearAlign=model_year_align_ew,             &
        offset=0,                                  &
        domFilePath='',                            &
        domFileName=trim(stream_fldFileName_ew),   &
        domTvarName='time',                        &
        domXvarName='lon' ,                        &
        domYvarName='lat' ,                        &  
        domAreaName='area',                        &
        domMaskName='mask',                        &
        filePath='',                               &
        filename=(/trim(stream_fldFileName_ew)/),  &
        fldListFile='ew_app',                      & ! this is the field name in the NetCDF file
        fldListModel='ew_app',                     &
        fillalgo='none',                           &
        mapalgo=ew_mapalgo,                        &
        tintalgo='nearest',                        &
        calendar=get_calendar(),                   &
        taxmode='extend'                           )

   if (masterproc) then
      call shr_strdata_print(sdat_app,'enhanced weathering application rate data')
   endif

   ! enhanced weathering weight percentage of mineral in rock (kg mineral kg-1 rock)
   fldList = shr_string_listCreateField( nminerals,  'ew_min')

   call shr_strdata_create(sdat_min,name="ew_mindyn",      &
        pio_subsystem=pio_subsystem,               &
        pio_iotype=shr_pio_getiotype(inst_name),   &
        mpicom=mpicom, compid=comp_id,             &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_elm,    &
        nxg=ldomain%ni, nyg=ldomain%nj,            &
        yearFirst=stream_year_first_ew,            &
        yearLast=stream_year_last_ew,              &
        yearAlign=model_year_align_ew,             &
        offset=0,                                  &
        domFilePath='',                            &
        domFileName=trim(stream_fldFileName_ew),   &
        domTvarName='time',                        &
        domXvarName='lon' ,                        &
        domYvarName='lat' ,                        &  
        domAreaName='area',                        &
        domMaskName='mask',                        &
        filePath='',                               &
        filename=(/trim(stream_fldFileName_ew)/),  &
        fldListFile=fldList,                       & ! this is the field name in the NetCDF file
        fldListModel=fldList,                      &
        fillalgo='none',                           &
        mapalgo=ew_mapalgo,                        &
        tintalgo='nearest',                        &
        calendar=get_calendar(),                   &
        taxmode='extend'                           )

   if (masterproc) then
      call shr_strdata_print(sdat_min,'enhanced weathering weight percentage of minerals in rock (kg mineral kg-1 rock)')
   endif

   ! enhanced weathering weight percentage of phosphorus content in rock (gP kg-1 rock)
   call shr_strdata_create(sdat_pho,name="ew_phodyn",      &
        pio_subsystem=pio_subsystem,               &
        pio_iotype=shr_pio_getiotype(inst_name),   &
        mpicom=mpicom, compid=comp_id,             &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_elm,    &
        nxg=ldomain%ni, nyg=ldomain%nj,            &
        yearFirst=stream_year_first_ew,            &
        yearLast=stream_year_last_ew,              &
        yearAlign=model_year_align_ew,             &
        offset=0,                                  &
        domFilePath='',                            &
        domFileName=trim(stream_fldFileName_ew),   &
        domTvarName='time',                        &
        domXvarName='lon' ,                        &
        domYvarName='lat' ,                        &  
        domAreaName='area',                        &
        domMaskName='mask',                        &
        filePath='',                               &
        filename=(/trim(stream_fldFileName_ew)/),  &
        fldListFile='ew_pho',                      & ! this is the field name in the NetCDF file
        fldListModel='ew_pho',                     &
        fillalgo='none',                           &
        mapalgo=ew_mapalgo,                        &
        tintalgo='nearest',                        &
        calendar=get_calendar(),                   &
        taxmode='extend'                           )

   if (masterproc) then
      call shr_strdata_print(sdat_min,'enhanced weathering weight percentage of phosphorus in rock (gP kg-1 rock)')
   endif

   ! enhanced weathering grain size (um diameter)
   call shr_strdata_create(sdat_gra,name="ew_gradyn",      &
        pio_subsystem=pio_subsystem,               &
        pio_iotype=shr_pio_getiotype(inst_name),   &
        mpicom=mpicom, compid=comp_id,             &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_elm,    &
        nxg=ldomain%ni, nyg=ldomain%nj,            &
        yearFirst=stream_year_first_ew,            &
        yearLast=stream_year_last_ew,              &
        yearAlign=model_year_align_ew,             &
        offset=0,                                  &
        domFilePath='',                            &
        domFileName=trim(stream_fldFileName_ew),   &
        domTvarName='time',                        &
        domXvarName='lon' ,                        &
        domYvarName='lat' ,                        &  
        domAreaName='area',                        &
        domMaskName='mask',                        &
        filePath='',                               &
        filename=(/trim(stream_fldFileName_ew)/),  &
        fldListFile='ew_gra',                      & ! this is the field name in the NetCDF file
        fldListModel='ew_gra',                     &
        fillalgo='none',                           &
        mapalgo=ew_mapalgo,                        &
        tintalgo='nearest',                        &
        calendar=get_calendar(),                   &
        taxmode='extend'                           )

   if (masterproc) then
      call shr_strdata_print(sdat_gra,'enhanced weathering grain size (um diameter)')
   endif

   ! enhanced weathering soil pH
   call shr_strdata_create(sdat_sph,name="ew_sphdyn",      &
        pio_subsystem=pio_subsystem,               &
        pio_iotype=shr_pio_getiotype(inst_name),   &
        mpicom=mpicom, compid=comp_id,             &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_elm,    &
        nxg=ldomain%ni, nyg=ldomain%nj,            &
        yearFirst=stream_year_first_ew,            &
        yearLast=stream_year_last_ew,              &
        yearAlign=model_year_align_ew,             &
        offset=0,                                  &
        domFilePath='',                            &
        domFileName=trim(stream_fldFileName_ew),   &
        domTvarName='time',                        &
        domXvarName='lon' ,                        &
        domYvarName='lat' ,                        &  
        domAreaName='area',                        &
        domMaskName='mask',                        &
        filePath='',                               &
        filename=(/trim(stream_fldFileName_ew)/),  &
        fldListFile='ew_sph',                      & ! this is the field name in the NetCDF file
        fldListModel='ew_sph',                     &
        fillalgo='none',                           &
        mapalgo=ew_mapalgo,                        &
        tintalgo='nearest',                        & ! temporal interpolation algorithm = nearest
        calendar=get_calendar(),                   &
        taxmode='extend'                           ) ! use the last value beyond time range instead of cycling

   if (masterproc) then
      call shr_strdata_print(sdat_sph,'enhanced weathering soil pH')
   endif

   ! percentage of the rock powder applied on croplands (= 100% - percent applieid on natveg)
   call shr_strdata_create(sdat_pctcrop,name="ew_pctcropdyn",      &
        pio_subsystem=pio_subsystem,               &
        pio_iotype=shr_pio_getiotype(inst_name),   &
        mpicom=mpicom, compid=comp_id,             &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_elm,    &
        nxg=ldomain%ni, nyg=ldomain%nj,            &
        yearFirst=stream_year_first_ew,            &
        yearLast=stream_year_last_ew,              &
        yearAlign=model_year_align_ew,             &
        offset=0,                                  &
        domFilePath='',                            &
        domFileName=trim(stream_fldFileName_ew),   &
        domTvarName='time',                        &
        domXvarName='lon' ,                        &
        domYvarName='lat' ,                        &  
        domAreaName='area',                        &
        domMaskName='mask',                        &
        filePath='',                               &
        filename=(/trim(stream_fldFileName_ew)/),  &
        fldListFile='ew_pct_crop',                 &
        fldListModel='ew_pct_crop',                &
        fillalgo='none',                           &
        mapalgo=ew_mapalgo,                        &
        tintalgo='nearest',                        & ! temporal interpolation algorithm = nearest
        calendar=get_calendar(),                   &
        taxmode='extend'                           ) ! use the last value beyond time range instead of cycling

   if (masterproc) then
      call shr_strdata_print(sdat_pctcrop,'enhanced weathering percent applied on croplands')
   endif

 end subroutine ew_init

 !================================================================
 subroutine ew_interp(bounds)
   !-----------------------------------------------------------------------
   ! Application rate should be once a year, on doy_application_app
   !-----------------------------------------------------------------------
   use elm_time_manager, only : get_curr_date, get_days_per_year
   use elm_varcon      , only : secspday
   use landunit_varcon , only : istsoil, istcrop
   use ColumnType      , only : col_pp
   use LandunitType    , only : lun_pp
   !
   ! Arguments
   type(bounds_type)           , intent(in)    :: bounds
   !
   ! Local variables
   integer :: g, gc, l, ic, ig, m, c
   integer :: year    ! year (0, ...) for nstep+1
   integer :: mon     ! month (1, ..., 12) for nstep+1
   integer :: day     ! day of month (1, ..., 31) for nstep+1
   integer :: sec     ! seconds into current date for nstep+1
   integer :: mcdate  ! Current model date (yyyymmdd)
   integer :: dayspyr ! days per year
   real(r8):: pctcol  ! percent rock powder going to crop column or natural vegetated column (%)
   character(len=CL) :: mineral_index
   !-----------------------------------------------------------------------

   call get_curr_date(year, mon, day, sec)
   mcdate = year*10000 + mon*100 + day

   call shr_strdata_advance(sdat_app, mcdate, sec, mpicom, 'ew_appdyn')
   call shr_strdata_advance(sdat_min, mcdate, sec, mpicom, 'ew_mindyn')
   call shr_strdata_advance(sdat_pho, mcdate, sec, mpicom, 'ew_phodyn')
   call shr_strdata_advance(sdat_gra, mcdate, sec, mpicom, 'ew_gradyn')
   call shr_strdata_advance(sdat_sph, mcdate, sec, mpicom, 'ew_sphdyn')
   call shr_strdata_advance(sdat_pctcrop, mcdate, sec, mpicom, 'ew_pctcropdyn')

   do c = bounds%begc, bounds%endc

      ! limit to vegetated columns
      l = col_pp%landunit(c)
      if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then

         ! Determine the vector index corresponding to soil column
         gc = col_pp%gridcell(c)
         ig = 0
         do g = bounds%begg,bounds%endg
            ig = ig+1
            if (g == gc) exit
         end do

         ! only fill in application rates one day per year
         if ( mon*100 + day == doy_application_ew ) then

            ! Determine the percent rock powder given to this column (crop or natveg)
            pctcol = sdat_pctcrop%avs(1)%rAttr(1,ig)
            if (lun_pp%itype(l) == istsoil) then
               pctcol = 100 - pctcol
            end if

            ! kg m-2 yr-1 => kg m-2 s-1
            col_ew%forc_app(c) = pctcol / 100 * sdat_app%avs(1)%rAttr(1,ig) / secspday
         end if

         ! the other values are needed every day

         ! % weight of up to five minerals in rock powder
         do m = 1, nminerals
            ! get the stream index
            write (mineral_index, "(i6)") m
            mineral_index = 'ew_min_'//trim(adjustl(mineral_index))
            ic = mct_aVect_indexRA(sdat_min%avs(1), trim(mineral_index))

            col_ew%forc_min(c,m) = sdat_min%avs(1)%rAttr(ic,ig)
         end do

         col_ew%forc_pho(c) = sdat_pho%avs(1)%rAttr(1,ig)
         col_ew%forc_gra(c) = sdat_gra%avs(1)%rAttr(1,ig)
         col_ew%forc_sph(c) = sdat_sph%avs(1)%rAttr(1,ig)

      else
         col_ew%forc_app(c) = 0._r8
         do m = 1, nminerals
            col_ew%forc_min(c, m) = 0._r8
         end do
         col_ew%forc_pho(c) = 0._r8
         col_ew%forc_gra(c) = 0._r8
         col_ew%forc_sph(c) = 0._r8
      end if

   end do

  end subroutine ew_interp

end module ewStreamMod