module lnd_import_export

  use shr_kind_mod , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use abortutils   , only: endrun
  use decompmod    , only: bounds_type
  use lnd2atmType  , only: lnd2atm_type
  use lnd2glcMod   , only: lnd2glc_type
  use atm2lndType  , only: atm2lnd_type
  use glc2lndMod   , only: glc2lnd_type
  use lnd2iacMod   , only: lnd2iac_type
  use iac2lndMod   , only: iac2lnd_type
  use elm_varctl   , only: iac_present, iulog
  use elm_varpar   , only: numpft, numharvest
  use ocn2lndType  , only: ocn2lnd_type
  use GridcellType , only: grc_pp          ! for access to gridcell topology
  use TopounitDataType , only: top_as, top_af  ! atmospheric state and flux variables
  use elm_cpl_indices
  use mct_mod
  use seq_flds_mod    , only : rof_sed
  !
  implicit none
  !===============================================================================

  ! Constants to compute vapor pressure over water/ice (shared by several
  ! routines below), and Kelvin<->Celsius helper -- module-scope statement
  ! functions can't be used outside the routine that declares them in
  ! Fortran, so these are ordinary functions instead.

contains

  !===============================================================================
  ! Kelvin -> bounded Celsius; and saturation vapor pressure over water/ice (Pa).
  ! Pulled out of lnd_import's local statement functions so both
  ! lnd_import_cpl_bypass and lnd_import_coupler (and the shared downstream
  ! helper) can use them without duplicating the coefficient tables.
  !===============================================================================
  pure function tdc(t)
    use shr_const_mod, only: SHR_CONST_TKFRZ
    real(r8), intent(in) :: t
    real(r8) :: tdc
    tdc = min(50._r8, max(-50._r8, (t - SHR_CONST_TKFRZ)))
  end function tdc

  pure function esatw(t)
    real(r8), intent(in) :: t
    real(r8) :: esatw
    real(r8), parameter :: a0=6.107799961_r8    , a1=4.436518521e-01_r8, &
         a2=1.428945805e-02_r8, a3=2.650648471e-04_r8, &
         a4=3.031240396e-06_r8, a5=2.034080948e-08_r8, &
         a6=6.136820929e-11_r8
    esatw = 100._r8*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
  end function esatw

  pure function esati(t)
    real(r8), intent(in) :: t
    real(r8) :: esati
    real(r8), parameter :: b0=6.109177956_r8    , b1=5.034698970e-01_r8, &
         b2=1.886013408e-02_r8, b3=4.176223716e-04_r8, &
         b4=5.824720280e-06_r8, b5=4.838803174e-08_r8, &
         b6=1.838826904e-10_r8
    esati = 100._r8*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))
  end function esati

  !===============================================================================
  ! lnd_import: dispatcher. Keeps the historical public signature so both
  ! lnd_comp_mct.F90 call sites (HAVE_MOAB vs. MCT -- an axis unrelated to
  ! CPL_BYPASS vs. coupler) are untouched. Under CPL_BYPASS, forcing is read
  ! from files and x2l is unused for meteorology; otherwise forcing comes
  ! from the coupler's x2l array directly.
  !===============================================================================
  subroutine lnd_import( bounds, x2l, atm2lnd_vars, glc2lnd_vars, ocn2lnd_vars, lnd2atm_vars, iac2lnd_vars)
    !
    ! !ARGUMENTS:
    type(bounds_type)  , intent(in)    :: bounds   ! bounds
    real(r8)           , intent(in)    :: x2l(:,:) ! driver import state to land model
    type(atm2lnd_type) , intent(inout) :: atm2lnd_vars      ! elm internal input data type
    type(glc2lnd_type) , intent(inout) :: glc2lnd_vars      ! elm internal input data type
    type(ocn2lnd_type) , intent(inout) :: ocn2lnd_vars
    type(lnd2atm_type) , intent(in)    :: lnd2atm_vars
    type(iac2lnd_type) , intent(inout) :: iac2lnd_vars ! elm iac to land
    !---------------------------------------------------------------------------

#ifdef CPL_BYPASS
    call lnd_import_cpl_bypass( bounds, x2l, atm2lnd_vars, glc2lnd_vars, ocn2lnd_vars, lnd2atm_vars, iac2lnd_vars)
#else
    call lnd_import_coupler( bounds, x2l, atm2lnd_vars, glc2lnd_vars, ocn2lnd_vars, lnd2atm_vars, iac2lnd_vars)
#endif

  end subroutine lnd_import

#ifdef CPL_BYPASS
  !===============================================================================
  ! lnd_import_cpl_bypass: meteorological forcing read directly from NetCDF
  ! files on disk, bypassing the CIME coupler. x2l is still accepted (and used
  ! for the small number of non-meteorological fields -- flood/runoff/ocean/
  ! ice/glc/iac -- that always arrive via the coupler even in CPL_BYPASS runs)
  ! but is not used for meteorology.
  !
  ! Design (see docs/coupler_bypass.md and this file's git history for the
  ! rationale): every variable in a CPL_BYPASS dataset shares one DTIME time
  ! axis (validated once at setup, not assumed). Time indices/interpolation
  ! weights are computed once per timestep, not per gridcell or per variable.
  ! The in-memory page buffer (atm2lnd_vars%atm_input) is sized by the number
  ! of *distinct* source columns actually referenced by this rank's
  ! gridcells, not one slot per gridcell -- the spatial remap from unique
  ! columns to per-gridcell forc_*_grc values happens fresh every timestep,
  ! not at page-read time.
  !===============================================================================
  subroutine lnd_import_cpl_bypass( bounds, x2l, atm2lnd_vars, glc2lnd_vars, ocn2lnd_vars, lnd2atm_vars, iac2lnd_vars)
    !
    ! !USES:
    use elm_varctl       , only: co2_type, co2_ppmv, iulog, use_c13, create_glacier_mec_landunit, &
                                 metdata_type, metdata_bypass, metdata_biases, co2_file, aero_file, use_atm_downscaling_to_topunit
    use elm_varctl       , only: const_climate_hist, add_temperature, add_co2, use_cn, use_fates
    use elm_varctl       , only: startdate_add_temperature, startdate_add_co2
    use elm_varctl       , only: metdata_read_buffer_ntimes, metdata_fsds_interpolant
    use elm_varcon       , only: rair, o2_molar_const, c13ratio
    use elm_time_manager , only: get_nstep, get_step_size, get_curr_calday, get_curr_date
    use controlMod       , only: NLFilename
    use shr_const_mod    , only: SHR_CONST_TKFRZ, SHR_CONST_STEBOL
    use domainMod        , only: ldomain
    use fileutils        , only: getavu, relavu
    use spmdmod          , only: masterproc, mpicom, iam, npes, MPI_REAL8, MPI_INTEGER, MPI_STATUS_SIZE
    use elm_nlUtilsMod   , only : find_nlgroup_name
    use FrictionVelocityMod, only: implicit_stress, atm_gustiness
    use lnd_disagg_forc
    use lnd_downscale_atm_forcing
    use netcdf
    !
    ! !ARGUMENTS:
    type(bounds_type)  , intent(in)    :: bounds
    real(r8)           , intent(in)    :: x2l(:,:)
    type(atm2lnd_type) , intent(inout) :: atm2lnd_vars
    type(glc2lnd_type) , intent(inout) :: glc2lnd_vars
    type(ocn2lnd_type) , intent(inout) :: ocn2lnd_vars
    type(lnd2atm_type) , intent(in)    :: lnd2atm_vars
    type(iac2lnd_type) , intent(inout) :: iac2lnd_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: g,topo,i,m,thism,nstep,ier
    integer status(MPI_STATUS_SIZE)
    real(r8) :: forc_rainc, forc_rainl, forc_snowc, forc_snowl
    real(r8) :: e, ea, qsat, forc_t, forc_q, forc_pbot
    real(r8) :: co2_ppmv_diag, co2_ppmv_prog, co2_ppmv_val
    integer  :: co2_type_idx
    real(r8) :: vp
    integer  :: thisng, np, num, nu_nml, nml_error
    integer, allocatable :: lnfmind_x_all(:), lnfmind_y_all(:)
    integer, allocatable :: lnfmind_x_local(:), lnfmind_y_local(:)
    real(r8), allocatable :: lnfm_sendbuf(:,:), lnfm_recvbuf(:,:)
    integer, allocatable :: recvcounts(:), displs(:)
    real(r8) :: swndf, swndr, swvdf, swvdr, ratio_rvrf, frac
    real(r8) :: swndr2, swndr3, swvdr2, swvdr3
    real(r8) :: thiscosz, avgcosz, szenith
    integer  :: swrad_period_len, swrad_period_start, thishr, thismin
    real(r8) :: latixy(500000), longxy(500000)
    integer ::  ierr, varid, dimid, yr, mon, day, tod, nindex(2), caldaym(13)
    integer ::  ncid, mask_ncid, thisncid, ng, tm
    integer ::  aindex(2), starti(3), counti(3)
    integer ::  grid_map(500000), zone_map(500000)
    integer ::  met_nvars, nyears_spinup, nyears_trans, starti_site, endi_site
    real(r8) :: smap05_lat(360), smap05_lon(720)
    real(r8) :: smapt62_lat(94), smapt62_lon(192)
    real(r8) :: smap2_lat(96), smap2_lon(144)
    real(r8) :: thisdist, mindist, thislon
    real(r8) :: tbot, tempndep(1,1,158), thiscalday, wt1(14), wt2(14), thisdoy
    real(r8) :: site_metdata(14,12)
    real(r8) :: var_month_mean(12)
    integer  :: var_month_count(12)
    integer :: xtoget, ytoget, thisx, thisy, calday_start
    integer :: sdate_addt, sy_addt, sm_addt, sd_addt
    integer :: sdate_addco2, sy_addco2, sm_addco2, sd_addco2
    character(len=200) metsource_str, thisline
    character(len=*), parameter :: sub = 'lnd_import_cpl_bypass'
    integer :: av, v, n, nummetdims, g3, gtoget, ztoget, line, mystart, tod_start, thistimelen
    integer, allocatable :: gtoget_all(:), ztoget_all(:)
    integer :: met_ngrid
    character(len=20) aerovars(14), metvars(14)
    character(len=3) zst
    integer :: stream_year_first_lightng, stream_year_last_lightng, model_year_align_lightng
    integer :: stream_year_first_popdens, stream_year_last_popdens, model_year_align_popdens
    integer :: stream_year_first_ndep,    stream_year_last_ndep,    model_year_align_ndep
    character(len=CL)  :: metdata_fname
    character(len=CL)  :: lightngmapalgo = 'bilinear'
    character(len=CL)  :: popdensmapalgo = 'bilinear'
    character(len=CL)  :: ndepmapalgo    = 'bilinear'
    character(len=CL)  :: stream_fldFileName_lightng
    character(len=CL)  :: stream_fldFileName_popdens
    character(len=CL)  :: stream_fldFileName_ndep
    logical :: use_sitedata, has_zonefile, use_daymet, use_livneh
    integer  :: dtime
    integer  :: col, t1, t2
    integer*2, allocatable :: site_range_buf(:)
    data caldaym / 1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 /
    !---------------------------------------------------------------------------

    namelist /light_streams/         &
        stream_year_first_lightng,  &
        stream_year_last_lightng,   &
        model_year_align_lightng,   &
        lightngmapalgo,             &
        stream_fldFileName_lightng

    namelist /popd_streams/          &
        stream_year_first_popdens,  &
        stream_year_last_popdens,   &
        model_year_align_popdens,   &
        popdensmapalgo,             &
        stream_fldFileName_popdens

    namelist /ndepdyn_nml/        &
        stream_year_first_ndep,  &
        stream_year_last_ndep,   &
        model_year_align_ndep,   &
        ndepmapalgo,             &
        stream_fldFileName_ndep

    stream_fldFileName_lightng = ' '
    stream_fldFileName_popdens = ' '

    co2_type_idx = 0
    if (co2_type == 'prognostic') then
       co2_type_idx = 1
    else if (co2_type == 'diagnostic') then
       co2_type_idx = 2
    end if
    if (co2_type == 'prognostic' .and. index_x2l_Sa_co2prog == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2prog for co2_type equal to prognostic' )
    else if (co2_type == 'diagnostic' .and. index_x2l_Sa_co2diag == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2diag for co2_type equal to diagnostic' )
    end if

    thisng = bounds%endg - bounds%begg + 1

    ! Lightning: read file and bcast lat/lon coords before the g loop so all
    ! ranks participate in collective calls at the same point.
    if (use_cn .or. use_fates) then
      if (atm2lnd_vars%loaded_bypassdata .eq. 0) then
        if (masterproc) then
          nu_nml = getavu()
          open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
          call find_nlgroup_name(nu_nml, 'light_streams', status=nml_error)
          if (nml_error == 0) then
            read(nu_nml, nml=light_streams, iostat=nml_error)
            if (nml_error /= 0) call endrun(msg='ERROR reading light_streams namelist')
          end if
          close(nu_nml)
          call relavu( nu_nml )
          allocate(atm2lnd_vars%lnfm_all(192,94,2920))
          ierr = nf90_open(trim(stream_fldFileName_lightng), NF90_NOWRITE, ncid)
          ierr = nf90_inq_varid(ncid, 'lat', varid)
          ierr = nf90_get_var(ncid, varid, smapt62_lat)
          ierr = nf90_inq_varid(ncid, 'lon', varid)
          ierr = nf90_get_var(ncid, varid, smapt62_lon)
          ierr = nf90_inq_varid(ncid, 'lnfm', varid)
          ierr = nf90_get_var(ncid, varid, atm2lnd_vars%lnfm_all)
          ierr = nf90_close(ncid)
        end if
        call mpi_bcast(smapt62_lon, 192, MPI_REAL8, 0, mpicom, ier)
        call mpi_bcast(smapt62_lat, 94,  MPI_REAL8, 0, mpicom, ier)
        allocate(lnfmind_x_local(thisng))
        allocate(lnfmind_y_local(thisng))
      end if
    end if

    call get_curr_date( yr, mon, day, tod )
    thiscalday = get_curr_calday()
    nstep = get_nstep()
    dtime = get_step_size()

    !this subroutine is called with nstep = 0 and nstep = 2,3,..., but not nstep=1.
    !Similarly, current date is called at begin time and 2x dtime, but not at 1x dtime
    !Conceptually, the rest of this code acts as if the tod is at the END of the step
    if (nstep.eq.0) then
       nstep = 1
       if ((tod + dtime).eq.86400) then
          day = day + 1
          tod = 0
       else
          tod = tod + dtime
       end if
    endif

    met_nvars=7
    if (metdata_type(1:3) == 'cpl') met_nvars=14

    if (atm2lnd_vars%loaded_bypassdata == 0) then
      call cpl_bypass_setup_metsource(atm2lnd_vars, metdata_type, metsource_str, use_livneh, use_daymet)
      call cpl_bypass_setup_varnames(atm2lnd_vars, metvars, met_nvars)

      inquire(file=trim(metdata_biases), exist=use_sitedata)

      !Read zone_mappings.txt once for all cells
      inquire(file=trim(metdata_bypass) // '/zone_mappings.txt', exist=has_zonefile)
      if (has_zonefile) then
        open(unit=13, file=trim(metdata_bypass) // '/zone_mappings.txt')
      else if (atm2lnd_vars%metsource .ne. 2) then
        call endrun( sub//' ERROR: Zone mapping file does not exist for cpl_bypass' )
      end if

      if (atm2lnd_vars%metsource .ne. 2) then
        ng = 0
        do v=1,500000
          read(13,*, end=10), longxy(v), latixy(v), zone_map(v), grid_map(v)
          ng = ng + 1
        end do
10      continue
        close(unit=13)
      end if

      !Pre-compute gtoget and ztoget for every local gridcell
      allocate(gtoget_all(bounds%begg:bounds%endg))
      allocate(ztoget_all(bounds%begg:bounds%endg))
      do g = bounds%begg, bounds%endg
        if (atm2lnd_vars%metsource .ne. 2) then
          mindist = 99999
          do g3 = 1,ng
            thisdist = 100*((latixy(g3) - ldomain%latc(g))**2 + &
                            (longxy(g3) - ldomain%lonc(g))**2)**0.5
            if (thisdist .lt. mindist) then
              mindist = thisdist
              ztoget_all(g) = zone_map(g3)
              gtoget_all(g) = grid_map(g3)
            end if
          end do
        else
          gtoget_all(g) = 1
          ztoget_all(g) = 0
        end if
      end do

      ! Build the unique-column map: atm2lnd_vars%column_index(g) indexes
      ! directly into atm_input's compact column dimension, replacing the
      ! old gtoget_saved/ztoget_saved pair. Each distinct (zone, file-row)
      ! pair referenced by any local gridcell gets exactly one column slot.
      call cpl_bypass_build_column_index(bounds, atm2lnd_vars, gtoget_all, ztoget_all)

      ! Needs a real zone (atm2lnd_vars%col_ztoget(1), just resolved above by
      ! cpl_bypass_build_column_index) to open the right per-zone file and read
      ! start/end year directly from it where the dataset supports that --
      ! never a hard-coded zone, since datasets like this rank's GSWP3
      ! Oakharbor extraction may not even have a zone 1 file on disk.
      call cpl_bypass_setup_year_defaults(atm2lnd_vars, use_livneh, use_daymet, metdata_bypass, &
                                       metsource_str, atm2lnd_vars%col_ztoget(1), &
                                       nyears_spinup, nyears_trans)

      atm2lnd_vars%page_ntimes = metdata_read_buffer_ntimes

      ! One-time dataset setup: opens each variable's file, reads DTIME's
      ! CF units/spacing/length, validates all variables share one time
      ! axis, and does the initial page load. yr/mon/day/tod/dtime are
      ! passed in (already carrying the nstep==0 bootstrap correction
      ! applied above) rather than re-read via get_curr_date here, so the
      ! very first call's initial tindex/page load agrees with every
      ! subsequent per-timestep call's cpl_bypass_update_tindex.
      call cpl_bypass_setup_dataset(bounds, atm2lnd_vars, met_nvars, metvars, &
                                  metsource_str, use_livneh, use_daymet, &
                                  yr, mon, day, tod, dtime)

      !Compute site bias corrections for each local cell
      mystart = atm2lnd_vars%startyear_met

      do g = bounds%begg, bounds%endg
        atm2lnd_vars%var_offset(:,g,:) = 0._r8
        atm2lnd_vars%var_mult(:,g,:)   = 1._r8

        if (use_sitedata) then
          open(unit=9, file=trim(metdata_biases),status='old')
          read(9,*) thisline
          site_metdata(:,:)=-999._r8
          do while ((site_metdata(1,1) .lt. ldomain%lonc(g) - 0.01 .or. &
                  site_metdata(1,1) .gt. ldomain%lonc(g) + 0.01) .and. &
                    (site_metdata(2,1) .lt. ldomain%latc(g) - 0.01 .or. &
                     site_metdata(2,1) .gt. ldomain%latc(g) + 0.01))
            read(9,*) site_metdata(1:7,1)
            if (site_metdata(1,1) .lt. 0) site_metdata(1,1) = site_metdata(1,1)+360._r8
          end do
          do line=2,12
            read(9,*) site_metdata(1:7,line)
          end do
          close(unit=9)

          col = atm2lnd_vars%column_index(g)
          do v=1,met_nvars
            if (v == 1) then
              starti_site = max((nint(site_metdata(4,1))-atm2lnd_vars%startyear_met) * &
                                   365*nint(24./atm2lnd_vars%timeres)+1,1)
              endi_site   = (min(atm2lnd_vars%endyear_met_trans,nint(site_metdata(5,1))) - &
                                   atm2lnd_vars%startyear_met+1)*(365*nint(24./atm2lnd_vars%timeres))
            end if
            ! This bias-correction scan needs a multi-year range [starti_site,
            ! endi_site] that is normally far larger than the rolling page
            ! buffer -- read it directly from file rather than through
            ! atm2lnd_vars%atm_input (which only holds the current page).
            call cpl_bypass_read_range(atm2lnd_vars, metvars(v), metsource_str, use_livneh, use_daymet, &
                                     atm2lnd_vars%col_ztoget(col), atm2lnd_vars%col_gtoget(col), &
                                     starti_site, endi_site, site_range_buf)
            var_month_mean(:)  = 0._r8
            var_month_count(:) = 0
            do i=starti_site, endi_site
              thisdoy = mod(i,365*nint(24./atm2lnd_vars%timeres))/(nint(24./atm2lnd_vars%timeres))+1
              do m=1,12
                if (thisdoy .ge. caldaym(m) .and. thisdoy .lt. caldaym(m+1)) thism = m
              enddo
              var_month_mean(thism) = var_month_mean(thism) + (site_range_buf(i-starti_site+1)* &
                                        atm2lnd_vars%scale_factors(v) + atm2lnd_vars%add_offsets(v))
              var_month_count(thism) = var_month_count(thism)+1
            end do
            deallocate(site_range_buf)
            do m = 1,12
              var_month_mean(m) = var_month_mean(m)/var_month_count(m)
              if (v .eq. 1) atm2lnd_vars%var_offset(v,g,m) = (site_metdata(6,m)+SHR_CONST_TKFRZ) - var_month_mean(m)
              if (v .eq. 5 .and. var_month_mean(m) .gt. 0) &
                    atm2lnd_vars%var_mult(v,g,m) = (site_metdata(7,m))/(caldaym(m+1)-caldaym(m))/24._r8/ &
                                                    3600._r8 / var_month_mean(m)
            end do
          end do
        end if

        if (.not. const_climate_hist) then
          if (yr .lt. mystart) then
            call endrun(msg=' ERROR: lnd_import_cpl_bypass: RUN_STARTDATE is before the start of' // &
                            ' the met forcing record. If running a spinup, set const_climate_hist=.true.' // &
                            ' in user_nl_elm to cycle the met record.')
          end if
          if (yr .gt. atm2lnd_vars%endyear_met_trans) then
            call endrun(msg=' ERROR: lnd_import_cpl_bypass: RUN_STARTDATE is after the end of' // &
                            ' the met forcing record -- check metdata_bypass and RUN_STARTDATE')
          end if
        end if
      end do

      deallocate(gtoget_all)
      deallocate(ztoget_all)

    end if   !end loaded_bypassdata == 0

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)

       ! Determine flooding input, sign convention is positive downward and
       ! hierarchy is atm/glc/lnd/rof/ice/ocn.  so water sent from rof to land is negative,
       ! change the sign to indicate addition of water to system.

       atm2lnd_vars%forc_flood_grc(g)   = 0._r8
       atm2lnd_vars%volr_grc(g)   = 0._r8
       atm2lnd_vars%volrmch_grc(g)= x2l(index_x2l_Flrr_volrmch,i) * (ldomain%area(g) * 1.e6_r8)
       atm2lnd_vars%supply_grc(g) = x2l(index_x2l_Flrr_supply,i)
       atm2lnd_vars%deficit_grc(g) = x2l(index_x2l_Flrr_deficit,i)

       if (index_x2l_Sr_h2orof /= 0) then
         atm2lnd_vars%h2orof_grc(g)      = x2l(index_x2l_Sr_h2orof,i)
         atm2lnd_vars%frac_h2orof_grc(g) = x2l(index_x2l_Sr_frac_h2orof,i)
       endif

       if (index_x2l_So_ssh /= 0) then
         ocn2lnd_vars%ssh_grc(g) = x2l(index_x2l_So_ssh,i)
       endif

       if (index_x2l_So_frac_h2oocn /= 0) then
         ocn2lnd_vars%frac_h2oocn_grc(g) = x2l(index_x2l_So_frac_h2oocn,i)
       endif
    end do

    !=========================================================================
    ! Shared time axis: compute the current bracketing record indices and
    ! interpolation weights ONCE, before any gridcell loop -- every variable
    ! in the dataset shares this same axis (validated at setup).
    !=========================================================================
    call cpl_bypass_update_tindex(atm2lnd_vars, yr, mon, day, tod, dtime, const_climate_hist)

    ! Lazily reload the page if the current bracketing indices fall outside it.
    if (atm2lnd_vars%page_ntimes .gt. 0) then
      call cpl_bypass_ensure_page(bounds, atm2lnd_vars, met_nvars, metvars, metsource_str, &
                                use_livneh, use_daymet, const_climate_hist)
    end if

    t1 = cpl_bypass_buffer_offset(atm2lnd_vars, atm2lnd_vars%tindex(1))
    t2 = cpl_bypass_buffer_offset(atm2lnd_vars, atm2lnd_vars%tindex(2))

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       col = atm2lnd_vars%column_index(g)

        ! wt1/wt2 start from the shared interpolation weight every gridcell
        ! (all variables share one time axis), then individual slots (wt2(4)
        ! for shortwave zenith reshaping, wt1(1)/wt2(1) for fire/aerosol/CO2,
        ! wt2(8:10) for metsource==5) get overwritten below within this same
        ! iteration -- reset fresh each gridcell so no iteration inherits a
        ! stale override from the previous one.
        wt1 = atm2lnd_vars%interp_wt1
        wt2 = atm2lnd_vars%interp_wt2

        !Air temperature
        atm2lnd_vars%forc_t_not_downscaled_grc(g)  = min(((atm2lnd_vars%atm_input(1,col,t1)*atm2lnd_vars%scale_factors(1)+ &
                                                      atm2lnd_vars%add_offsets(1))*wt1(1) + (atm2lnd_vars%atm_input(1,col,t2)* &
                                                      atm2lnd_vars%scale_factors(1)+atm2lnd_vars%add_offsets(1))*wt2(1)) * &
                                                      atm2lnd_vars%var_mult(1,g,mon) + atm2lnd_vars%var_offset(1,g,mon), 323._r8)
        atm2lnd_vars%forc_th_not_downscaled_grc(g) = min(((atm2lnd_vars%atm_input(1,col,t1)*atm2lnd_vars%scale_factors(1)+ &
                                                      atm2lnd_vars%add_offsets(1))*wt1(1) + (atm2lnd_vars%atm_input(1,col,t2)* &
                                                      atm2lnd_vars%scale_factors(1)+atm2lnd_vars%add_offsets(1))*wt2(1)) * &
                                                      atm2lnd_vars%var_mult(1,g,mon) + atm2lnd_vars%var_offset(1,g,mon), 323._r8)

        tbot = atm2lnd_vars%forc_t_not_downscaled_grc(g)

        !Air pressure
        atm2lnd_vars%forc_pbot_not_downscaled_grc(g) = max(((atm2lnd_vars%atm_input(2,col,t1)*atm2lnd_vars%scale_factors(2)+ &
                                                        atm2lnd_vars%add_offsets(2))*wt1(2) + (atm2lnd_vars%atm_input(2,col,t2) &
                                                        *atm2lnd_vars%scale_factors(2)+atm2lnd_vars%add_offsets(2))*wt2(2)) * &
                                                        atm2lnd_vars%var_mult(2,g,mon) + atm2lnd_vars%var_offset(2,g,mon), 4e4_r8)
        !Specific humidity
        atm2lnd_vars%forc_q_not_downscaled_grc(g) = max(((atm2lnd_vars%atm_input(3,col,t1)*atm2lnd_vars%scale_factors(3)+ &
                                                     atm2lnd_vars%add_offsets(3))*wt1(3) + (atm2lnd_vars%atm_input(3,col,t2) &
                                                     *atm2lnd_vars%scale_factors(3)+atm2lnd_vars%add_offsets(3))*wt2(3)) * &
                                                     atm2lnd_vars%var_mult(3,g,mon) + atm2lnd_vars%var_offset(3,g,mon), 1e-9_r8)

        if (atm2lnd_vars%metsource == 2) then  !convert RH to qbot
          if (tbot > SHR_CONST_TKFRZ) then
            e = esatw(tdc(tbot))
          else
            e = esati(tdc(tbot))
          end if
          qsat           = 0.622_r8*e / (atm2lnd_vars%forc_pbot_not_downscaled_grc(g) - 0.378_r8*e)
          atm2lnd_vars%forc_q_not_downscaled_grc(g) = qsat * atm2lnd_vars%forc_q_not_downscaled_grc(g) / 100.0_r8
        end if

        !use longwave from file if provided
        atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) = ((atm2lnd_vars%atm_input(7,col,t1)*atm2lnd_vars%scale_factors(7)+ &
                                                        atm2lnd_vars%add_offsets(7))*wt1(7) + (atm2lnd_vars%atm_input(7,col,t2) &
                                                        *atm2lnd_vars%scale_factors(7)+atm2lnd_vars%add_offsets(7))*wt2(7)) * &
                                                        atm2lnd_vars%var_mult(7,g,mon) + atm2lnd_vars%var_offset(7,g,mon)
        if (atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) .le. 50 .or. atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) .ge. 600) then
        !Longwave radiation (calculated from air temperature, humidity)
            e =  atm2lnd_vars%forc_pbot_not_downscaled_grc(g) * atm2lnd_vars%forc_q_not_downscaled_grc(g) / &
                 (0.622_R8 + 0.378_R8 * atm2lnd_vars%forc_q_not_downscaled_grc(g) )
            ea = 0.70_R8 + 5.95e-05_R8 * 0.01_R8 * e * exp(1500.0_R8/tbot)
            atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) = ea * SHR_CONST_STEBOL * tbot**4
        end if

        !Shortwave radiation (cosine zenith angle interpolation, or plain
        !linear interpolation between forcing records -- see
        !metdata_fsds_interpolant).
        if (trim(metdata_fsds_interpolant) == 'linear') then
          wt2(4) = 1.0_r8
        else
          thishr = (tod-get_step_size()/2)/3600
          if (thishr < 0) thishr=thishr+24
          thismin = mod((tod-get_step_size()/2)/60, 60)
          thiscosz = max(cos(szenith(ldomain%lonc(g),ldomain%latc(g),0,int(thiscalday),thishr,thismin,0)* &
                          3.14159265358979/180.0d0), 0.001d0)
          avgcosz = 0d0
          if (atm2lnd_vars%npf - 1._r8 .gt. 1e-3) then
            swrad_period_len   = get_step_size()*nint(atm2lnd_vars%npf)
            swrad_period_start = ((tod-get_step_size()/2)/swrad_period_len) * swrad_period_len
            !set to last period if first model timestep of the day
            if (tod-get_step_size()/2 < 0) swrad_period_start = ((86400-get_step_size()/2)/swrad_period_len) * swrad_period_len

            do tm=1,nint(atm2lnd_vars%npf)
              !Get the average cosine zenith angle over the time resolution of the input data
              thishr  = (swrad_period_start+(tm-1)*get_step_size()+get_step_size()/2)/3600
              if (thishr > 23) thishr=thishr-24
              thismin = mod((swrad_period_start+(tm-1)*get_step_size()+get_step_size()/2)/60, 60)
              avgcosz  = avgcosz + max(cos(szenith(ldomain%lonc(g),ldomain%latc(g),0,int(thiscalday),thishr, thismin, 0) &
                         *3.14159265358979/180.0d0), 0.001d0)/atm2lnd_vars%npf
            end do
          else
            avgcosz = thiscosz
          end if
          if (thiscosz > 0.001d0) then
            wt2(4) = min(thiscosz/avgcosz, 10.0_r8)
          else
            wt2(4) = 0d0
          end if
        end if

        if (atm2lnd_vars%metsource == 5) then
            wt2(4)=1.0   !cosz interp not working
            wt2(8:10)=1.0
            swndf = max(((atm2lnd_vars%atm_input(4,col,t2)*atm2lnd_vars%scale_factors(4)+ &
                                     atm2lnd_vars%add_offsets(4))*wt2(4)), 0.0_r8)
            swndr = max(((atm2lnd_vars%atm_input(8,col,t2)*atm2lnd_vars%scale_factors(8)+ &
                                     atm2lnd_vars%add_offsets(8))*wt2(8)), 0.0_r8)
            swvdf = max(((atm2lnd_vars%atm_input(9,col,t2)*atm2lnd_vars%scale_factors(9)+ &
                                     atm2lnd_vars%add_offsets(9))*wt2(9)), 0.0_r8)
            swvdr = max(((atm2lnd_vars%atm_input(10,col,t2)*atm2lnd_vars%scale_factors(10)+ &
                                     atm2lnd_vars%add_offsets(10))*wt2(10)), 0.0_r8)
            atm2lnd_vars%forc_solad_grc(g,2) = swndr
            atm2lnd_vars%forc_solad_grc(g,1) = swvdr
            atm2lnd_vars%forc_solai_grc(g,2) = swndf
            atm2lnd_vars%forc_solai_grc(g,1) = swvdf
        else
            ! compute shared SW input once to prevent compiler FP optimization producing NaN/Inf
            frac = max((atm2lnd_vars%atm_input(4,col,t2)*atm2lnd_vars%scale_factors(4) &
                 + atm2lnd_vars%add_offsets(4)) * wt2(4) * 0.50_R8, 0.0_r8)

            swndr = frac
            swndf = frac
            swvdr = frac
            swvdf = frac

            swndr2 = swndr * swndr
            swndr3 = swndr2 * swndr
            ratio_rvrf =   min(0.99_R8, max(0.01_R8, 0.29548_R8 + 0.00504_R8*swndr - 1.4957e-05_R8 * swndr2 + 1.4881e-08_R8 * swndr3))
            atm2lnd_vars%forc_solad_grc(g,2) = ratio_rvrf*swndr
            atm2lnd_vars%forc_solai_grc(g,2) = (1._R8 - ratio_rvrf)*swndf

            swvdr2 = swvdr * swvdr
            swvdr3 = swvdr2 * swvdr
            ratio_rvrf =   min(0.99_R8,max(0.17639_R8 + 0.00380_R8*swvdr - 9.0039e-06_R8*swvdr2 + 8.1351e-09_R8*swvdr3,0.01_R8))
            atm2lnd_vars%forc_solad_grc(g,1) = ratio_rvrf*swvdr
            atm2lnd_vars%forc_solai_grc(g,1) = (1._R8 - ratio_rvrf)*swvdf
        end if
        !Rain and snow
        if (atm2lnd_vars%metsource == 5) then
          forc_rainc = max((((atm2lnd_vars%atm_input(5,col,t2)*atm2lnd_vars%scale_factors(5)+ &
                                        atm2lnd_vars%add_offsets(5)))*atm2lnd_vars%var_mult(5,g,mon) + &
                                        atm2lnd_vars%var_offset(5,g,mon)), 0.0_r8)
          forc_rainl = max((((atm2lnd_vars%atm_input(11,col,t2)*atm2lnd_vars%scale_factors(11)+ &
                                        atm2lnd_vars%add_offsets(11)))*atm2lnd_vars%var_mult(11,g,mon) + &
                                        atm2lnd_vars%var_offset(11,g,mon)), 0.0_r8)
          forc_snowc = max((((atm2lnd_vars%atm_input(12,col,t2)*atm2lnd_vars%scale_factors(12)+ &
                                        atm2lnd_vars%add_offsets(12)))*atm2lnd_vars%var_mult(12,g,mon) + &
                                        atm2lnd_vars%var_offset(12,g,mon)), 0.0_r8)
          forc_snowl = max((((atm2lnd_vars%atm_input(13,col,t2)*atm2lnd_vars%scale_factors(13)+ &
                                        atm2lnd_vars%add_offsets(13)))*atm2lnd_vars%var_mult(13,g,mon) + &
                                          atm2lnd_vars%var_offset(13,g,mon)), 0.0_r8)
        else
          frac = (atm2lnd_vars%forc_t_not_downscaled_grc(g) - SHR_CONST_TKFRZ)*0.5_R8       ! ramp near freezing
          frac = min(1.0_R8,max(0.0_R8,frac))           ! bound in [0,1]
          !Don't interpolate rainfall data
          forc_rainc = 0.1_R8 * frac * max((((atm2lnd_vars%atm_input(5,col,t2)*atm2lnd_vars%scale_factors(5)+ &
                                        atm2lnd_vars%add_offsets(5)))*atm2lnd_vars%var_mult(5,g,mon) + &
                                        atm2lnd_vars%var_offset(5,g,mon)), 0.0_r8)
          forc_rainl = 0.9_R8 * frac * max((((atm2lnd_vars%atm_input(5,col,t2)*atm2lnd_vars%scale_factors(5)+ &
                                         atm2lnd_vars%add_offsets(5)))*atm2lnd_vars%var_mult(5,g,mon) + &
                                         atm2lnd_vars%var_offset(5,g,mon)), 0.0_r8)
          forc_snowc = 0.1_R8 * (1.0_R8 - frac) * max((((atm2lnd_vars%atm_input(5,col,t2)*atm2lnd_vars%scale_factors(5)+ &
                  atm2lnd_vars%add_offsets(5)))*atm2lnd_vars%var_mult(5,g,mon) + atm2lnd_vars%var_offset(5,g,mon)), 0.0_r8)
          forc_snowl = 0.9_R8 * (1.0_R8 - frac) * max((((atm2lnd_vars%atm_input(5,col,t2)*atm2lnd_vars%scale_factors(5)+ &
                  atm2lnd_vars%add_offsets(5))) * atm2lnd_vars%var_mult(5,g,mon) + atm2lnd_vars%var_offset(5,g,mon)), 0.0_r8)
        end if
        !Wind
        atm2lnd_vars%forc_u_grc(g) = (atm2lnd_vars%atm_input(6,col,t1)*atm2lnd_vars%scale_factors(6)+ &
                                     atm2lnd_vars%add_offsets(6))*wt1(6) + (atm2lnd_vars%atm_input(6,col,t2)* &
                                     atm2lnd_vars%scale_factors(6)+atm2lnd_vars%add_offsets(6))*wt2(6)
        if (atm2lnd_vars%metsource == 5) then
          atm2lnd_vars%forc_v_grc(g) = (atm2lnd_vars%atm_input(14,col,t1)*atm2lnd_vars%scale_factors(14)+ &
                                     atm2lnd_vars%add_offsets(14))*wt1(14) + (atm2lnd_vars%atm_input(14,col,t2)* &
                                     atm2lnd_vars%scale_factors(14)+atm2lnd_vars%add_offsets(14))*wt2(14)
        else
            atm2lnd_vars%forc_v_grc(g) = 0.0_R8
        end if
        atm2lnd_vars%forc_hgt_grc(g) = 30.0_R8   ! zgcmxy  Atm state, default=30m

  !------------------------------------Fire data -------------------------------------------------------

        nindex(1) = yr-1848
        nindex(2) = nindex(1)+1
        if (yr .lt. 1850 .or. const_climate_hist) nindex(1:2) = 2
        if (yr .ge. 2010 .and. .not. const_climate_hist) nindex(1:2) = 161

        model_filter: if (use_cn .or. use_fates) then
          if (atm2lnd_vars%loaded_bypassdata == 0 .or. (mon .eq. 1 .and. day .eq. 1 .and. tod .eq. 0)) then
            if (masterproc .and. i .eq. 1) then
              ! Read pop_dens streams namelist to get filename
              nu_nml = getavu()
              open(nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
              call find_nlgroup_name(nu_nml, 'popd_streams', status=nml_error)
              if (nml_error == 0) then
                  read(nu_nml, nml=popd_streams,iostat=nml_error)
                  if (nml_error /= 0) then
                      call endrun(msg='ERROR reading popdens namelist')
                  end if
              end if
              close(nu_nml)
              call relavu( nu_nml )

              ierr = nf90_open(trim(stream_fldFileName_popdens), NF90_NOWRITE, ncid)
              ierr = nf90_inq_varid(ncid, 'lat', varid)
              ierr = nf90_get_var(ncid, varid, smap05_lat)
              ierr = nf90_inq_varid(ncid, 'lon', varid)
              ierr = nf90_get_var(ncid, varid, smap05_lon)
              ierr = nf90_inq_varid(ncid, 'hdm', varid)
              starti(1:2) = 1
              starti(3)   = nindex(1)
              counti(1) = 720
              counti(2) = 360
              counti(3) = 1
              ierr = nf90_get_var(ncid, varid, atm2lnd_vars%hdm1, starti, counti)
              starti(3) = nindex(2)
              if (nindex(1) .ne. nindex(2)) then
                  ierr = nf90_get_var(ncid, varid, atm2lnd_vars%hdm2, starti, counti)
              else
                  atm2lnd_vars%hdm2 = atm2lnd_vars%hdm1
              end if
              ierr = nf90_close(ncid)
            end if

            if (i .eq. 1) then
              call mpi_bcast (atm2lnd_vars%hdm1, 360*720, MPI_REAL8, 0, mpicom, ier)
              call mpi_bcast (atm2lnd_vars%hdm2, 360*720, MPI_REAL8, 0, mpicom, ier)
              call mpi_bcast (smap05_lon, 720, MPI_REAL8, 0, mpicom, ier)
              call mpi_bcast (smap05_lat, 360, MPI_REAL8, 0, mpicom, ier)
            end if
          end if

          !figure out which point to get
          if (atm2lnd_vars%loaded_bypassdata == 0) then
            mindist=99999
            do thisx = 1,720
              do thisy = 1,360
                  if (ldomain%lonc(g) .lt. 0) then
                      if (smap05_lon(thisx) >= 180) smap05_lon(thisx) = smap05_lon(thisx)-360._r8
                  else if (ldomain%lonc(g) .ge. 180) then
                      if (smap05_lon(thisx) < 0) smap05_lon(thisx) = smap05_lon(thisx) + 360._r8
                  end if
                  thisdist = 100*((smap05_lat(thisy) - ldomain%latc(g))**2 + &
                          (smap05_lon(thisx) - ldomain%lonc(g))**2)**0.5
                  if (thisdist .lt. mindist) then
                      mindist = thisdist
                      atm2lnd_vars%hdmind(g,1) = thisx
                      atm2lnd_vars%hdmind(g,2) = thisy
                  end if
              end do
            end do
          end if
          !get weights for interpolation
          wt1(1) = 1._r8 - (thiscalday -1._r8)/365._r8
          wt2(1) = 1._r8 - wt1(1)
          atm2lnd_vars%forc_hdm(g) = atm2lnd_vars%hdm1(atm2lnd_vars%hdmind(g,1),atm2lnd_vars%hdmind(g,2),1)*wt1(1) + &
                                     atm2lnd_vars%hdm2(atm2lnd_vars%hdmind(g,1),atm2lnd_vars%hdmind(g,2),1)*wt2(1)

          if (atm2lnd_vars%loaded_bypassdata .eq. 0) then
            mindist=99999
            do thisx = 1,192
              do thisy = 1,94
                if (ldomain%lonc(g) .lt. 0) then
                  if (smapt62_lon(thisx) >= 180) smapt62_lon(thisx) = smapt62_lon(thisx)-360._r8
                else if (ldomain%lonc(g) .ge. 180) then
                  if (smapt62_lon(thisx) < 0) smapt62_lon(thisx) = smapt62_lon(thisx) + 360._r8
                end if
                thisdist = 100*((smapt62_lat(thisy) - ldomain%latc(g))**2 + &
                            (smapt62_lon(thisx) - ldomain%lonc(g))**2)**0.5
                if (thisdist .lt. mindist) then
                  mindist = thisdist
                  lnfmind_x_local(i) = thisx
                  lnfmind_y_local(i) = thisy
                end if
              end do
            end do
          end if

   !------------------------------------Nitrogen deposition----------------------------------------------

          !DMR note - ndep will NOT be correct if more than 1850 years of model
          !spinup (model year > 1850)
          nindex(1) = min(max(yr-1848,2), 168)
          nindex(2) = min(nindex(1)+1, 168)

          if (atm2lnd_vars%loaded_bypassdata .eq. 0 .or. (mon .eq. 1 .and. day .eq. 1 .and. tod .eq. 0)) then
            if (masterproc .and. i .eq. 1) then
              nu_nml = getavu()
              open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
              call find_nlgroup_name(nu_nml, 'ndepdyn_nml', status=nml_error)
              if (nml_error == 0) then
                read(nu_nml, nml=ndepdyn_nml,iostat=nml_error)
                if (nml_error /= 0) then
                  call endrun(msg='ERROR reading ndep namelist')
                end if
              end if
              close(nu_nml)
              call relavu( nu_nml )

              ierr = nf90_open(trim(stream_fldFileName_ndep), nf90_nowrite, ncid)
              ierr = nf90_inq_varid(ncid, 'lat', varid)
              ierr = nf90_get_var(ncid, varid, smap2_lat)
              ierr = nf90_inq_varid(ncid, 'lon', varid)
              ierr = nf90_get_var(ncid, varid, smap2_lon)
              ierr = nf90_inq_varid(ncid, 'NDEP_year', varid)
              starti(1:2) = 1
              starti(3)   = nindex(1)
              counti(1)   = 144
              counti(2)   = 96
              counti(3)   = 1
              ierr = nf90_get_var(ncid, varid, atm2lnd_vars%ndep1, starti, counti)
              if (nindex(1) .ne. nindex(2)) then
                starti(3) = nindex(2)
                ierr = nf90_get_var(ncid, varid, atm2lnd_vars%ndep2, starti, counti)
              else
                atm2lnd_vars%ndep2 = atm2lnd_vars%ndep1
              end if
              ierr = nf90_close(ncid)
             end if
             if (i .eq. 1) then
               call mpi_bcast (atm2lnd_vars%ndep1, 144*96, MPI_REAL8, 0, mpicom, ier)
               call mpi_bcast (atm2lnd_vars%ndep2, 144*96, MPI_REAL8, 0, mpicom, ier)
               call mpi_bcast (smap2_lon, 144, MPI_REAL8, 0, mpicom, ier)
               call mpi_bcast (smap2_lat, 96, MPI_REAL8, 0, mpicom, ier)
             end if
          end if

          if (atm2lnd_vars%loaded_bypassdata .eq. 0) then
            mindist=99999
            do thisx = 1,144
              do thisy = 1,96
                if (ldomain%lonc(g) .lt. 0) then
                  if (smap2_lon(thisx) >= 180) smap2_lon(thisx) = smap2_lon(thisx)-360._r8
                else if (ldomain%lonc(g) .ge. 180) then
                  if (smap2_lon(thisx) < 0) smap2_lon(thisx) = smap2_lon(thisx) + 360._r8
                end if
                thislon = smap2_lon(thisx)
                thisdist = 100*((smap2_lat(thisy) - ldomain%latc(g))**2 + &
                              (thislon - ldomain%lonc(g))**2)**0.5
                if (thisdist .lt. mindist) then
                  mindist = thisdist
                  atm2lnd_vars%ndepind(g,1) = thisx
                  atm2lnd_vars%ndepind(g,2) = thisy
                end if
              end do
            end do
          end if

          !get weights for interpolation
          wt1(1) = 1._r8 - (thiscalday -1._r8)/365._r8
          wt2(1) = 1._r8 - wt1(1)

          atm2lnd_vars%forc_ndep_grc(g)    = (atm2lnd_vars%ndep1(atm2lnd_vars%ndepind(g,1),atm2lnd_vars%ndepind(g,2),1)*wt1(1) + &
                                              atm2lnd_vars%ndep2(atm2lnd_vars%ndepind(g,1),atm2lnd_vars%ndepind(g,2),1)*wt2(1)) / (365._r8 * 86400._r8)
        end if model_filter

   !------------------------------------Aerosol forcing--------------------------------------------------
        if (atm2lnd_vars%loaded_bypassdata .eq. 0 .or. (mon .eq. 1 .and. day .eq. 1 .and. tod .eq. 0)) then
          if (masterproc .and. i .eq. 1) then
            aerovars(1) = 'BCDEPWET'
            aerovars(2) = 'BCPHODRY'
            aerovars(3) = 'BCPHIDRY'
            aerovars(4) = 'OCDEPWET'
            aerovars(5) = 'OCPHODRY'
            aerovars(6) = 'OCPHIDRY'
            aerovars(7) = 'DSTX01DD'
            aerovars(8) = 'DSTX02DD'
            aerovars(9) = 'DSTX03DD'
            aerovars(10) = 'DSTX04DD'
            aerovars(11) = 'DSTX01WD'
            aerovars(12) = 'DSTX02WD'
            aerovars(13) = 'DSTX03WD'
            aerovars(14) = 'DSTX04WD'
            ierr = nf90_open(trim(aero_file), nf90_nowrite, ncid)
            ierr = nf90_inq_varid(ncid, 'lat', varid)
            ierr = nf90_get_var(ncid, varid, smap2_lat)
            ierr = nf90_inq_varid(ncid, 'lon', varid)
            ierr = nf90_get_var(ncid, varid, smap2_lon)
            starti(1:2) = 1
            starti(3)   = max((min(yr,2100)-1849)*12+1, 13)-1
            counti(1)   = 144
            counti(2)   = 96
            counti(3)   = 14
            do av=1,14
              ierr = nf90_inq_varid(ncid, trim(aerovars(av)), varid)
              ierr = nf90_get_var(ncid, varid, atm2lnd_vars%aerodata(av,:,:,:), starti, counti)
            end do
            ierr = nf90_close(ncid)
          end if
          if (i .eq. 1) then
             call mpi_bcast (atm2lnd_vars%aerodata, 14*144*96*14, MPI_REAL8, 0, mpicom, ier)
          end if
        end if

        !Use ndep grid indices since they're on the same grid
        if (atm2lnd_vars%loaded_bypassdata .eq. 0 .and. (.not. (use_fates .or. use_cn) )   ) then
            mindist=99999
            do thisx = 1,144
              do thisy = 1,96
                if (ldomain%lonc(g) .lt. 0) then
                  if (smap2_lon(thisx) >= 180) smap2_lon(thisx) = smap2_lon(thisx)-360._r8
                else if (ldomain%lonc(g) .ge. 180) then
                  if (smap2_lon(thisx) < 0) smap2_lon(thisx) = smap2_lon(thisx) + 360._r8
                end if
                thislon = smap2_lon(thisx)
                thisdist = 100*((smap2_lat(thisy) - ldomain%latc(g))**2 + &
                              (thislon - ldomain%lonc(g))**2)**0.5
                if (thisdist .lt. mindist) then
                  mindist = thisdist
                  atm2lnd_vars%ndepind(g,1) = thisx
                  atm2lnd_vars%ndepind(g,2) = thisy
                end if
              end do
            end do
        end if

        !get weights for interpolation (note this method doesn't get the month boundaries quite right..)
        aindex(1) = mon+1
        if (thiscalday .le. (caldaym(mon+1)+caldaym(mon))/2._r8) then
           wt1(1) = 0.5_r8 + (thiscalday-caldaym(mon))/(caldaym(mon+1)-caldaym(mon))
           aindex(2) = aindex(1)-1
        else
           wt1(1) = 1.0_r8 - (thiscalday-(caldaym(mon+1)+caldaym(mon))/2._r8)/   &
                          (caldaym(mon+1)-caldaym(mon))
           aindex(2) = aindex(1)+1
        end if
        wt2(1) = 1._r8 - wt1(1)

        do av = 1,14
          atm2lnd_vars%forc_aer_grc(g,av)  =  atm2lnd_vars%aerodata(av,atm2lnd_vars%ndepind(g,1), &
            atm2lnd_vars%ndepind(g,2),aindex(1))*wt1(1)+atm2lnd_vars%aerodata(av,atm2lnd_vars%ndepind(g,1), &
            atm2lnd_vars%ndepind(g,2),aindex(2))*wt2(1)
        end do

       !Parse startdate for adding temperature
       if (startdate_add_temperature .ne. '') then
         call get_curr_date( yr, mon, day, tod )
         read(startdate_add_temperature,*) sdate_addt
         sy_addt     = sdate_addt/10000
         sm_addt     = (sdate_addt-sy_addt*10000)/100
         sd_addt     = sdate_addt-sy_addt*10000-sm_addt*100
       end if
       if (startdate_add_temperature .ne. '') then
         if ((yr == sy_addt .and. mon == sm_addt .and. day >= sd_addt) .or. &
             (yr == sy_addt .and. mon > sm_addt) .or. (yr > sy_addt)) then
           atm2lnd_vars%forc_t_not_downscaled_grc(g) = atm2lnd_vars%forc_t_not_downscaled_grc(g) + add_temperature
           atm2lnd_vars%forc_th_not_downscaled_grc(g) = atm2lnd_vars%forc_th_not_downscaled_grc(g) + add_temperature
         end if
       end if

       ! Adding topographic downscaling capability within CPL_BYPASS code block
       if (use_atm_downscaling_to_topunit) then
        atm2lnd_vars%forc_uovern = x2l(index_x2l_Sa_uovern,i)
        atm2lnd_vars%forc_rain_not_downscaled_grc = forc_rainc + forc_rainl
        atm2lnd_vars%forc_snow_not_downscaled_grc = forc_snowc + forc_snowl

          if(atm_gustiness) then
             call endrun("Error: atm_gustiness not yet supported with multiple topounits (in CPL_BYPASS)")
          end if
         do topo = grc_pp%topi(g) , grc_pp%topf(g)
            top_as%ugust(topo) = 0._r8
         end do

         call downscale_atm_forcing_to_topounit_cpl_bypass(g, atm2lnd_vars, lnd2atm_vars)
       else
         do topo = grc_pp%topi(g), grc_pp%topf(g)
          top_as%tbot(topo)    = atm2lnd_vars%forc_t_not_downscaled_grc(g)      ! forc_txy  Atm state K
          top_as%thbot(topo)   = atm2lnd_vars%forc_th_not_downscaled_grc(g)     ! forc_thxy Atm state K
          top_as%pbot(topo)    = atm2lnd_vars%forc_pbot_not_downscaled_grc(g)   ! ptcmxy    Atm state Pa
          top_as%qbot(topo)    = atm2lnd_vars%forc_q_not_downscaled_grc(g)      ! forc_qxy  Atm state kg/kg
          top_as%ubot(topo)    = atm2lnd_vars%forc_u_grc(g)                     ! forc_uxy  Atm state m/s
          top_as%vbot(topo)    = atm2lnd_vars%forc_v_grc(g)                     ! forc_vxy  Atm state m/s
          top_as%zbot(topo)    = atm2lnd_vars%forc_hgt_grc(g)                   ! zgcmxy    Atm state m
          top_as%windbot(topo) = sqrt(top_as%ubot(topo)**2 + top_as%vbot(topo)**2)
          ! Relative humidity (percent)
          if (top_as%tbot(topo) > SHR_CONST_TKFRZ) then
             e = esatw(tdc(top_as%tbot(topo)))
          else
             e = esati(tdc(top_as%tbot(topo)))
          end if
          qsat           = 0.622_r8*e / (top_as%pbot(topo) - 0.378_r8*e)
          top_as%rhbot(topo) = 100.0_r8*(top_as%qbot(topo) / qsat)
          ! partial pressure of oxygen (Pa)
          top_as%po2bot(topo) = o2_molar_const * top_as%pbot(topo)
          ! air density (kg/m**3) - uses a temporary calculation of water vapor pressure (Pa)
          vp = top_as%qbot(topo) * top_as%pbot(topo)  / (0.622_r8 + 0.378_r8 * top_as%qbot(topo))
          top_as%rhobot(topo) = (top_as%pbot(topo) - 0.378_r8 * vp) / (rair * top_as%tbot(topo))
          top_af%rain(topo)    = forc_rainc + forc_rainl            ! sum of convective and large-scale rain
          top_af%snow(topo)    = forc_snowc + forc_snowl            ! sum of convective and large-scale snow
          top_af%solad(topo,2) = atm2lnd_vars%forc_solad_grc(g,2)   ! forc_sollxy  Atm flux  W/m^2
          top_af%solad(topo,1) = atm2lnd_vars%forc_solad_grc(g,1)   ! forc_solsxy  Atm flux  W/m^2
          top_af%solai(topo,2) = atm2lnd_vars%forc_solai_grc(g,2)   ! forc_solldxy Atm flux  W/m^2
          top_af%solai(topo,1) = atm2lnd_vars%forc_solai_grc(g,1)   ! forc_solsdxy Atm flux  W/m^2
          top_af%lwrad(topo)   = atm2lnd_vars%forc_lwrad_not_downscaled_grc(g)     ! flwdsxy Atm flux  W/m^2
          ! derived flux forcings
          top_af%solar(topo) = top_af%solad(topo,2) + top_af%solad(topo,1) + &
                              top_af%solai(topo,2) + top_af%solai(topo,1)
         end do
       end if

       !set the topounit-level atmospheric variables that are not handled in downscaling code
       do topo = grc_pp%topi(g), grc_pp%topf(g)
        if (implicit_stress) then
            top_as%wsresp(topo)  = 0._r8                                       !           Atm state m/s/Pa
            top_as%tau_est(topo) = 0._r8                                       !           Atm state Pa
         end if
         top_as%ugust(topo) = 0._r8                                            !           Atm state m/s
         if (atm_gustiness) then
            top_as%windbot(topo) = sqrt(top_as%windbot(topo)**2 + top_as%ugust(topo)**2)
         end if
       end do

!---------------------------------- CO2 -------------------------------------------------------------------

     if (co2_type_idx /= 0) then
        !atmospheric CO2 (to be used for transient simulations only)
        if (atm2lnd_vars%loaded_bypassdata .eq. 0) then
          ierr = nf90_open(trim(co2_file), nf90_nowrite, ncid)
          if (ierr .ne. 0) call endrun(msg=' ERROR: Failed to open cpl_bypass input CO2 file' )
          ierr = nf90_inq_dimid(ncid, 'time', dimid)
          ierr = nf90_Inquire_Dimension(ncid, dimid, len = thistimelen)
          ierr = nf90_inq_varid(ncid, 'CO2', varid)
          ierr = nf90_get_var(ncid, varid, atm2lnd_vars%co2_input(:,:,1:thistimelen))
          ierr = nf90_inq_varid(ncid, 'C13O2', varid)
          ierr = nf90_get_var(ncid, varid, atm2lnd_vars%c13o2_input(:,:,1:thistimelen))
          ierr = nf90_close(ncid)
        end if

        !get weights/indices for interpolation (assume values represent annual averages)
        nindex(1) = min(max(yr,1850),2100)-1764
        if (thiscalday .le. 182.5) then
          nindex(2) = nindex(1)-1
        else
          nindex(2) = nindex(1)+1
        end if
        wt1(1) = 1._r8 - abs((182.5 - (thiscalday -1._r8))/365._r8)
        wt2(1) = 1._r8 - wt1(1)

        co2_ppmv_val = atm2lnd_vars%co2_input(1,1,nindex(1))*wt1(1) + atm2lnd_vars%co2_input(1,1,nindex(2))*wt2(1)

        if (startdate_add_co2 .ne. '') then
          read(startdate_add_co2,*) sdate_addco2
          sy_addco2     = sdate_addco2/10000
          sm_addco2     = (sdate_addco2-sy_addco2*10000)/100
          sd_addco2     = sdate_addco2-sy_addco2*10000-sm_addt*100
          if ((yr == sy_addco2 .and. mon == sm_addco2 .and. day >= sd_addco2) .or. &
              (yr == sy_addco2 .and. mon > sm_addco2) .or. (yr > sy_addco2)) then
            co2_ppmv_val=co2_ppmv_val + add_co2
          end if
        end if

        if (use_c13) then
          atm2lnd_vars%forc_pc13o2_grc(g) = (atm2lnd_vars%c13o2_input(1,1,nindex(1))*wt1(1) + &
             atm2lnd_vars%c13o2_input(1,1,nindex(2))*wt2(1)) * 1.e-6_r8 * atm2lnd_vars%forc_pbot_not_downscaled_grc(g)
        end if

        ! bypass mode doesn't receive _prog/_diag from atm, but have to reset them here
        co2_ppmv_prog = co2_ppmv_val
        co2_ppmv_diag = co2_ppmv_val
     else if (co2_type_idx == 0) then

        ! CO2 constant, value from namelist
        co2_ppmv_val = co2_ppmv
        if (use_c13) then
          atm2lnd_vars%forc_pc13o2_grc(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 &
                                           * atm2lnd_vars%forc_pbot_not_downscaled_grc(g)
        end if

     else

        call endrun( sub//' ERROR: Invalid co2_type_idx, must be 0 or not (constant or diagnostic) for CPL_BYPASS' )

     end if

     call finalize_forcing_common(bounds, g, i, x2l, atm2lnd_vars, glc2lnd_vars, iac2lnd_vars, co2_ppmv_val, &
                                 co2_ppmv_prog, co2_ppmv_diag, co2_type_idx, use_c13, forc_rainc, forc_rainl, &
                                 forc_snowc, forc_snowl)

    end do

    ! Lightning: gather lnfmind from all ranks to rank 0, do lookups, scatter lnfm back.
    if (use_cn .or. use_fates) then
      if (atm2lnd_vars%loaded_bypassdata .eq. 0) then
        allocate(recvcounts(0:npes-1))
        allocate(displs(0:npes-1))
        call mpi_allgather(thisng, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, mpicom, ier)
        displs(0) = 0
        do np = 1, npes-1
          displs(np) = displs(np-1) + recvcounts(np-1)
        end do
        if (masterproc) then
          allocate(lnfmind_x_all(sum(recvcounts)))
          allocate(lnfmind_y_all(sum(recvcounts)))
        else
          allocate(lnfmind_x_all(1))  ! dummy on non-master
          allocate(lnfmind_y_all(1))
        end if
        call mpi_gatherv(lnfmind_x_local, thisng, MPI_INTEGER, &
             lnfmind_x_all, recvcounts, displs, MPI_INTEGER, 0, mpicom, ier)
        call mpi_gatherv(lnfmind_y_local, thisng, MPI_INTEGER, &
             lnfmind_y_all, recvcounts, displs, MPI_INTEGER, 0, mpicom, ier)

        ! rank 0 builds send buffer: lnfm_sendbuf(2920, total_ng)
        if (masterproc) then
          allocate(lnfm_sendbuf(2920, sum(recvcounts)))
          do np = 0, npes-1
            do num = 1, recvcounts(np)
              lnfm_sendbuf(:, displs(np)+num) = &
                   atm2lnd_vars%lnfm_all(lnfmind_x_all(displs(np)+num), &
                                         lnfmind_y_all(displs(np)+num), :)
            end do
          end do
          deallocate(atm2lnd_vars%lnfm_all)
        else
          allocate(lnfm_sendbuf(2920, 1))  ! dummy on non-master
        end if

        ! scatter lnfm back into a (2920, thisng) temp buffer then transpose into
        ! lnfm(begg:endg, 2920) since the two layouts are incompatible for direct scatterv
        allocate(lnfm_recvbuf(2920, thisng))
        call mpi_scatterv(lnfm_sendbuf, recvcounts*2920, displs*2920, MPI_REAL8, &
             lnfm_recvbuf, thisng*2920, MPI_REAL8, 0, mpicom, ier)
        do num = 1, thisng
          atm2lnd_vars%lnfm(bounds%begg+num-1, :) = lnfm_recvbuf(:, num)
        end do

        deallocate(lnfm_recvbuf)
        deallocate(lnfmind_x_local)
        deallocate(lnfmind_y_local)
        deallocate(lnfmind_x_all)
        deallocate(lnfmind_y_all)
        deallocate(lnfm_sendbuf)
        deallocate(recvcounts)
        deallocate(displs)

        ! assign forc_lnfm for each local gridcell
        do g = bounds%begg, bounds%endg
          atm2lnd_vars%forc_lnfm(g) = atm2lnd_vars%lnfm(g, ((int(thiscalday)-1)*8+tod/(3600*3))+1)
        end do
      end if
    end if

    atm2lnd_vars%loaded_bypassdata = 1

  end subroutine lnd_import_cpl_bypass

  !===============================================================================
  ! cpl_bypass_setup_metsource: parse metdata_type into the integer metsource code
  ! and the livneh/daymet regional-override flags. Filename-construction
  ! detail only -- never used to infer time metadata (see cpl_bypass_setup_dataset).
  !===============================================================================
  subroutine cpl_bypass_setup_metsource(atm2lnd_vars, metdata_type, metsource_str, use_livneh, use_daymet)
    type(atm2lnd_type), intent(inout) :: atm2lnd_vars
    character(len=*),   intent(in)    :: metdata_type
    character(len=200), intent(out)   :: metsource_str
    logical,             intent(out)  :: use_livneh, use_daymet
    character(len=*), parameter :: sub = 'cpl_bypass_setup_metsource'

    if (index(metdata_type, 'qian') .gt. 0) then
      atm2lnd_vars%metsource = 0
      metsource_str = 'qian'
    else if (index(metdata_type,'cru') .gt. 0) then
      atm2lnd_vars%metsource = 1
      metsource_str = 'cruncep'
    else if (index(metdata_type,'site') .gt. 0) then
      atm2lnd_vars%metsource = 2
      metsource_str = 'site'
    else if (index(metdata_type,'princeton') .gt. 0) then
      atm2lnd_vars%metsource = 3
      metsource_str = 'princeton'
    else if (index(metdata_type,'gswp3') .gt. 0) then
      atm2lnd_vars%metsource = 4
      metsource_str = 'gswp3'
    else if (index(metdata_type,'cpl') .gt. 0) then
      atm2lnd_vars%metsource = 5
      metsource_str = 'cpl'
    else if (index(metdata_type,'era5') .gt. 0) then
      atm2lnd_vars%metsource = 6
      metsource_str = 'era5'
    else if (index(metdata_type,'atssubdaily') .gt. 0) then
      atm2lnd_vars%metsource = 8
      metsource_str = 'atssubdaily'
    else
      call endrun( sub//' ERROR: Invalid met data source for cpl_bypass' )
    end if

    use_livneh = .false.
    use_daymet = .false.
    if (index(metdata_type, 'livneh') .gt. 0) then
      use_livneh = .true.
    else if (index(metdata_type, 'daymet') .gt. 0) then
      use_daymet = .true.
    end if
  end subroutine cpl_bypass_setup_metsource

  !===============================================================================
  ! cpl_bypass_setup_varnames: variable name list per metsource. Every metsource
  ! shares the same 7-variable layout (TBOT, PSRF, QBOT/RH, FSDS, PRECTmms,
  ! WIND, FLDS) except metsource==5 ('cpl'), which provides a genuinely
  ! different 14-variable set (radiative-band decomposition, split
  ! convective/large-scale precip phase, vector wind U/V).
  !===============================================================================
  subroutine cpl_bypass_setup_varnames(atm2lnd_vars, metvars, met_nvars)
    type(atm2lnd_type), intent(in)  :: atm2lnd_vars
    character(len=20),  intent(out) :: metvars(14)
    integer,             intent(in) :: met_nvars

    metvars(1) = 'TBOT'
    metvars(2) = 'PSRF'
    metvars(3) = 'QBOT'
    if (atm2lnd_vars%metsource .eq. 2) metvars(3) = 'RH'
    metvars(4) = 'FSDS'
    metvars(5) = 'PRECTmms'
    metvars(6) = 'WIND'
    metvars(7) = 'FLDS'
    if (atm2lnd_vars%metsource .eq. 5) then
      metvars(6) = 'U'
      metvars(8) = 'SWNDR'
      metvars(9) = 'SWVDF'
      metvars(10) = 'SWVDR'
      metvars(11) = 'RAINL'
      metvars(12) = 'SNOWC'
      metvars(13) = 'SNOWL'
      metvars(14) = 'V'
    end if
  end subroutine cpl_bypass_setup_varnames

  !===============================================================================
  ! cpl_bypass_setup_year_defaults: per-source calendar-year metadata.
  ! Where a dataset provides no start_year/end_year variables (GSWP3;
  ! CRUNCEP/Princeton/ERA5 still fall back to a hard-coded table below --
  ! no such data ships in this repo's inputdata to fix against), the real
  ! start year is instead read from the file's own DTIME:units CF string
  ! (e.g. "Days since 1901-01-01 00:00"), and the end year derived from
  ! DTIME's last value -- never hard-coded, since a dataset's actual
  ! coverage (e.g. a site-specific GSWP3 extraction) need not match the
  ! nominal "1901-2014" full record implied by the filename convention.
  ! site/atssubdaily read an explicit start_year/end_year pair directly
  ! instead, since those formats provide one.
  ! This is otherwise filename-adjacent bookkeeping only -- cpl_bypass_setup_dataset
  ! still always reads the true record length/resolution from DTIME itself.
  !===============================================================================
  subroutine cpl_bypass_setup_year_defaults(atm2lnd_vars, use_livneh, use_daymet, metdata_bypass, &
                                         metsource_str, ztoget, nyears_spinup, nyears_trans)
    use netcdf
    use shr_string_mod, only: shr_string_parseCFtunit, shr_string_toLower
    type(atm2lnd_type), intent(inout) :: atm2lnd_vars
    logical,             intent(in)   :: use_livneh, use_daymet
    character(len=*),    intent(in)   :: metdata_bypass
    character(len=*),    intent(in)   :: metsource_str
    integer,             intent(in)   :: ztoget
    integer,             intent(out)  :: nyears_spinup, nyears_trans
    integer :: ierr, ncid, varid
    character(len=256) :: metdata_fname, dtime_units
    character(len=80) :: dtime_unit_name
    integer :: dtime_bdate, this_startyear
    real(r8) :: dtime_bsec, dtime_first, dtime_last
    integer :: dtime_dimid, dtime_len

    atm2lnd_vars%startyear_met       = 1901
    atm2lnd_vars%endyear_met_spinup  = 1920
    if (atm2lnd_vars%metsource == 0) then
      atm2lnd_vars%startyear_met       = 1948
      atm2lnd_vars%endyear_met_spinup  = 1972
      atm2lnd_vars%endyear_met_trans   = 2004
    else if (atm2lnd_vars%metsource == 1) then
      atm2lnd_vars%endyear_met_trans  = 2016
    else if (atm2lnd_vars%metsource == 2) then
      ierr = nf90_open(trim(metdata_bypass) // '/all_hourly.nc', nf90_nowrite, ncid)
      ierr = nf90_inq_varid(ncid, 'start_year', varid)
      ierr = nf90_get_var(ncid, varid, atm2lnd_vars%startyear_met)
      ierr = nf90_inq_varid(ncid, 'end_year', varid)
      ierr = nf90_get_var(ncid, varid, atm2lnd_vars%endyear_met_spinup)
      ierr = nf90_close(ncid)
      atm2lnd_vars%endyear_met_trans = atm2lnd_vars%endyear_met_spinup
    else if (atm2lnd_vars%metsource == 3) then
      atm2lnd_vars%endyear_met_trans = 2012
    else if (atm2lnd_vars%metsource == 4) then
      ! Read the real start/end year from GSWP3's own DTIME axis instead of
      ! assuming the nominal "1901-2014" full record -- a site extraction
      ! (e.g. this rank's Oakharbor single-column pull) can cover a
      ! different, possibly shorter or offset, span. shr_string_parseCFtunit
      ! (CIME's share/ CF-time-units parser, already used elsewhere in ELM,
      ! e.g. SatellitePhenologyMod) handles the full CF units grammar
      ! (optional time-of-day/timezone, arbitrary whitespace) instead of
      ! hand-parsing the string here.
      metdata_fname = cpl_bypass_filename(atm2lnd_vars, 'TBOT', metsource_str, use_livneh, use_daymet, ztoget)
      ierr = nf90_open(trim(metdata_bypass) // '/' // trim(metdata_fname), nf90_nowrite, ncid)
      if (ierr .ne. 0) call endrun(msg=' ERROR: cpl_bypass_setup_year_defaults: '// &
                                    'failed to open '//trim(metdata_fname)//' to read DTIME:units')
      ierr = nf90_inq_varid(ncid, 'DTIME', varid)
      ierr = nf90_get_att(ncid, varid, 'units', dtime_units)
      ierr = nf90_inq_dimid(ncid, 'DTIME', dtime_dimid)
      ierr = nf90_inquire_dimension(ncid, dtime_dimid, len=dtime_len)
      ierr = nf90_get_var(ncid, varid, dtime_first, (/1/))
      ierr = nf90_get_var(ncid, varid, dtime_last, (/dtime_len/))
      ierr = nf90_close(ncid)

      ! shr_string_parseCFtunit matches lowercase 'days '/' since ' only --
      ! this dataset's DTIME:units is capitalized ("Days since ..."), so
      ! normalize case first rather than duplicating its parsing logic.
      call shr_string_parseCFtunit(shr_string_toLower(dtime_units), dtime_unit_name, dtime_bdate, dtime_bsec)
      this_startyear = dtime_bdate / 10000

      ! Add the whole years spanned by DTIME's range to get the end year
      ! (DTIME is days-since-start; integer division by 365 is exact enough
      ! for this record-selection purpose -- precise dates/leap-days are
      ! handled elsewhere via DTIME's values directly, not through
      ! startyear_met/endyear_met_trans).
      atm2lnd_vars%startyear_met      = this_startyear
      atm2lnd_vars%endyear_met_trans  = this_startyear + int((dtime_last - dtime_first) / 365._r8)
      ! Preserve the pre-existing behavior of always spinup-cycling through
      ! the record's first 20 years (module default endyear_met_spinup=1920
      ! previously did this only because startyear_met also defaulted to
      ! 1901) -- now relative to the real start year, clamped to the actual
      ! record length in case a dataset covers less than 20 years.
      atm2lnd_vars%endyear_met_spinup = min(this_startyear + 19, atm2lnd_vars%endyear_met_trans)
    else if (atm2lnd_vars%metsource == 5) then
      atm2lnd_vars%startyear_met      = 566
      atm2lnd_vars%endyear_met_spinup = 590
      atm2lnd_vars%endyear_met_trans  = 590
    else if (atm2lnd_vars%metsource == 6) then
      atm2lnd_vars%startyear_met      = 1950
      atm2lnd_vars%endyear_met_spinup = 1970
      atm2lnd_vars%endyear_met_trans  = 2025
    else if (atm2lnd_vars%metsource == 8) then
      ierr = nf90_open(trim(metdata_bypass) // '/ATS-subdaily_TBOT_z01.nc', nf90_nowrite, ncid)
      ierr = nf90_inq_varid(ncid, 'start_year', varid)
      ierr = nf90_get_var(ncid, varid, atm2lnd_vars%startyear_met)
      ierr = nf90_inq_varid(ncid, 'end_year', varid)
      ierr = nf90_get_var(ncid, varid, atm2lnd_vars%endyear_met_spinup)
      ierr = nf90_close(ncid)
      atm2lnd_vars%endyear_met_trans = atm2lnd_vars%endyear_met_spinup
    end if

    if (use_livneh) then
      atm2lnd_vars%startyear_met      = 1950
      atm2lnd_vars%endyear_met_spinup = 1969
    else if (use_daymet) then
      atm2lnd_vars%startyear_met      = 1980
      atm2lnd_vars%endyear_met_spinup = atm2lnd_vars%endyear_met_trans
    end if

    nyears_spinup = atm2lnd_vars%endyear_met_spinup - atm2lnd_vars%startyear_met + 1
    nyears_trans  = atm2lnd_vars%endyear_met_trans  - atm2lnd_vars%startyear_met + 1
  end subroutine cpl_bypass_setup_year_defaults

  !===============================================================================
  ! cpl_bypass_filename: filename-construction only, per metsource. Used solely
  ! to locate a file on disk -- never a source of time metadata.
  !===============================================================================
  function cpl_bypass_filename(atm2lnd_vars, metvar, metsource_str, use_livneh, use_daymet, ztoget) result(fname)
    type(atm2lnd_type), intent(in) :: atm2lnd_vars
    character(len=*),   intent(in) :: metvar, metsource_str
    logical,             intent(in) :: use_livneh, use_daymet
    integer,             intent(in) :: ztoget
    character(len=256) :: fname
    character(len=3)   :: zst

    write(zst, '(I3)') 100+ztoget
    if (atm2lnd_vars%metsource == 0) then
      fname = trim(metsource_str) // '_' // trim(metvar) // '_z' // zst(2:3) // '.nc'
    else if (atm2lnd_vars%metsource == 1) then
      fname = 'CRUNCEP.v5_' // trim(metvar) // '_1901-2013_z' // zst(2:3) // '.nc'
      if (use_livneh .and. ztoget .ge. 16 .and. ztoget .le. 20) then
        fname = 'CRUNCEP5_Livneh_' // trim(metvar) // '_1950-2013_z' // zst(2:3) // '.nc'
      else if (use_daymet .and. ztoget .ge. 16 .and. ztoget .le. 20) then
        fname = 'CRUNCEP5_Daymet3_' // trim(metvar) // '_1980-2013_z' // zst(2:3) // '.nc'
      end if
    else if (atm2lnd_vars%metsource == 2) then
      fname = 'all_hourly.nc'
    else if (atm2lnd_vars%metsource == 3) then
      fname = 'Princeton_' // trim(metvar) // '_1901-2012_z' // zst(2:3) // '.nc'
      if (use_livneh .and. ztoget .ge. 16 .and. ztoget .le. 20) then
        fname = 'Princeton_Livneh_' // trim(metvar) // '_1950-2012_z' // zst(2:3) // '.nc'
      else if (use_daymet .and. ztoget .ge. 16 .and. ztoget .le. 20) then
        fname = 'Princeton_Daymet3_' // trim(metvar) // '_1980-2013_z' // zst(2:3) // '.nc'
      end if
    else if (atm2lnd_vars%metsource == 4) then
      fname = 'GSWP3_' // trim(metvar) // '_1901-2014_z' // zst(2:3) // '.nc'
      if (use_livneh .and. ztoget .ge. 16 .and. ztoget .le. 20) then
        fname = 'GSWP3_Livneh_' // trim(metvar) // '_1950-2010_z' // zst(2:3) // '.nc'
      else if (use_daymet .and. ztoget .ge. 16 .and. ztoget .le. 20) then
        fname = 'GSWP3_Daymet3_' // trim(metvar) // '_1980-2010_z' // zst(2:3) // '.nc'
      end if
    else if (atm2lnd_vars%metsource == 5) then
      fname = 'CBGC1850S.ne30_' // trim(metvar) // '_0566-0590_z' // zst(2:3) // '.nc'
    else if (atm2lnd_vars%metsource == 6) then
      fname = 'ERA5_' // trim(metvar) // '_1950-2025_z' // zst(2:3) // '.nc'
    else if (atm2lnd_vars%metsource == 8) then
      fname = 'ATS-subdaily_' // trim(metvar) // '_z' // zst(2:3) // '.nc'
    end if
  end function cpl_bypass_filename

  !===============================================================================
  ! cpl_bypass_build_column_index: assign each local gridcell a flat index into atm_input's
  ! compact column dimension. Distinct local gridcells that share the same
  ! (zone, file-row) pair get the same column index -- computed once at setup,
  ! static for the run (per the design discussion: nothing about paging or
  ! spinup wrapping requires recomputing this later).
  !===============================================================================
  subroutine cpl_bypass_build_column_index(bounds, atm2lnd_vars, gtoget_all, ztoget_all)
    type(bounds_type),  intent(in)    :: bounds
    type(atm2lnd_type), intent(inout) :: atm2lnd_vars
    integer,             intent(in)   :: gtoget_all(bounds%begg:bounds%endg)
    integer,             intent(in)   :: ztoget_all(bounds%begg:bounds%endg)

    integer :: g, g2, ncols

    ncols = 0
    do g = bounds%begg, bounds%endg
      atm2lnd_vars%column_index(g) = 0
      do g2 = bounds%begg, g-1
        if (gtoget_all(g2) == gtoget_all(g) .and. ztoget_all(g2) == ztoget_all(g)) then
          atm2lnd_vars%column_index(g) = atm2lnd_vars%column_index(g2)
          exit
        end if
      end do
      if (atm2lnd_vars%column_index(g) == 0) then
        ncols = ncols + 1
        atm2lnd_vars%column_index(g) = ncols
        ! Record this new column's (zone, file-row) once -- cpl_bypass_load_page
        ! reads col_gtoget/col_ztoget on every (re)load rather than re-deriving
        ! uniqueness from begg:endg gridcells each time.
        atm2lnd_vars%col_gtoget(ncols) = gtoget_all(g)
        atm2lnd_vars%col_ztoget(ncols) = ztoget_all(g)
      end if
    end do
    atm2lnd_vars%ncols_needed = ncols
  end subroutine cpl_bypass_build_column_index

  !===============================================================================
  ! cpl_bypass_setup_dataset: one-time dataset setup. Opens each variable's file
  ! (for each unique zone), reads DTIME's CF units/spacing/length ONCE from
  ! the first variable/zone opened, and validates every other variable's file
  ! agrees (length and spacing) -- endrun with a clear message on mismatch.
  ! Then performs the initial page load.
  !===============================================================================
  subroutine cpl_bypass_setup_dataset(bounds, atm2lnd_vars, met_nvars, metvars, &
                                    metsource_str, use_livneh, use_daymet, &
                                    yr, mon, day, tod, dtime)
    use netcdf
    use elm_varctl,      only: metdata_bypass, const_climate_hist
    type(bounds_type),  intent(in)    :: bounds
    type(atm2lnd_type), intent(inout) :: atm2lnd_vars
    integer,             intent(in)   :: met_nvars
    character(len=20),   intent(in)   :: metvars(14)
    character(len=200),  intent(in)   :: metsource_str
    logical,             intent(in)   :: use_livneh, use_daymet
    ! yr/mon/day/tod/dtime: the caller's already-bootstrap-corrected clock
    ! (see the nstep==0 adjustment in lnd_import_cpl_bypass) -- read fresh
    ! via get_curr_date here instead would disagree with that correction on
    ! the very first call and could compute a negative eval_hr.
    integer,             intent(in)   :: yr, mon, day, tod, dtime

    character(len=*), parameter :: sub = 'cpl_bypass_setup_dataset'
    integer :: v, col, ierr, varid, dimid, ncid
    real(r8) :: timetemp(2), this_timeres
    integer :: this_timelen
    character(len=256) :: metdata_fname
    integer :: nyears_spinup

    nyears_spinup = atm2lnd_vars%endyear_met_spinup - atm2lnd_vars%startyear_met + 1

    ! Open each variable's file once per unique column (cpl_bypass_build_column_index has
    ! already reduced this rank's gridcells down to their distinct source
    ! columns, in col_ztoget/col_gtoget) to establish and validate the shared
    ! DTIME time axis. Reading metadata per-column rather than per-gridcell
    ! avoids reopening the same file redundantly for every gridcell that maps
    ! to it.
    do v = 1, met_nvars
      do col = 1, atm2lnd_vars%ncols_needed
        metdata_fname = cpl_bypass_filename(atm2lnd_vars, metvars(v), metsource_str, use_livneh, use_daymet, &
                                           atm2lnd_vars%col_ztoget(col))
        ierr = nf90_open(trim(metdata_bypass) // '/' // trim(metdata_fname), NF90_NOWRITE, ncid)
        if (ierr .ne. 0) call endrun(msg=' ERROR: '//sub//': Failed to open cpl_bypass input meteorology file: '//trim(metdata_fname))

        ierr = nf90_inq_dimid(ncid, 'DTIME', dimid)
        ierr = nf90_Inquire_Dimension(ncid, dimid, len=this_timelen)
        ierr = nf90_inq_varid(ncid, 'DTIME', varid)
        ierr = nf90_get_var(ncid, varid, timetemp, (/1/), (/2/))
        this_timeres = (timetemp(2)-timetemp(1))*24._r8

        if (v == 1 .and. col == 1) then
          !Reference: establish the shared axis from the first file opened.
          atm2lnd_vars%timelen = this_timelen
          atm2lnd_vars%timeres = this_timeres
          atm2lnd_vars%npf     = 86400d0*(timetemp(2)-timetemp(1))/real(dtime, r8)
          atm2lnd_vars%timelen_spinup = nyears_spinup*(365*nint(24./atm2lnd_vars%timeres))
        else
          !Validate every other variable/column file agrees with the reference axis.
          if (this_timelen /= atm2lnd_vars%timelen .or. &
              abs(this_timeres - atm2lnd_vars%timeres) .gt. 1.e-6_r8*max(this_timeres,atm2lnd_vars%timeres)) then
            call endrun(msg=' ERROR: '//sub//': DTIME axis mismatch between meteorology files -- '// &
                            trim(metdata_fname)//' does not share the reference time axis established by variable 1.'// &
                            ' Every variable in a CPL_BYPASS dataset must share one time axis.')
          end if
        end if

        if (col == 1) then
          atm2lnd_vars%scale_factors(v) = 1.0_r8
          atm2lnd_vars%add_offsets(v)   = 0.0_r8
          ierr = nf90_inq_varid(ncid, trim(metvars(v)), varid)
          ierr = nf90_get_att(ncid, varid, 'scale_factor', atm2lnd_vars%scale_factors(v))
          ierr = nf90_get_att(ncid, varid, 'add_offset',   atm2lnd_vars%add_offsets(v))
        end if

        ierr = nf90_close(ncid)
      end do
    end do

    !Resolve page depth: page_ntimes==-1 means "whole record, no paging."
    if (atm2lnd_vars%page_ntimes == -1) then
      if (.not. const_climate_hist .and. (yr .ge. 1850)) then
        atm2lnd_vars%page_ntimes = atm2lnd_vars%timelen
      else
        atm2lnd_vars%page_ntimes = atm2lnd_vars%timelen_spinup
      end if
    else if (atm2lnd_vars%page_ntimes .lt. 2) then
      ! A page must hold at least tindex(1) and tindex(2) (the one-slice
      ! overlap the paging design requires) -- a depth of 1 would leave no
      ! room for the bracketing "next" record.
      call endrun(msg=' ERROR: '//sub//': metdata_read_buffer_ntimes must be -1 (whole record)'// &
                      ' or >= 2 (a page must hold both bracketing time records).')
    end if

    !Initial time index/weights and page load.
    call cpl_bypass_update_tindex(atm2lnd_vars, yr, mon, day, tod, dtime, const_climate_hist)
    call cpl_bypass_load_page(bounds, atm2lnd_vars, met_nvars, metvars, metsource_str, use_livneh, use_daymet, &
                            atm2lnd_vars%tindex(1))

  end subroutine cpl_bypass_setup_dataset

  !===============================================================================
  ! cpl_bypass_update_tindex: compute the current bracketing absolute record
  ! indices tindex(1:2) and interpolation weights ONCE per timestep, shared
  ! across all variables (they share one time axis) and all gridcells (the
  ! shared clock doesn't depend on gridcell).
  !
  ! Definition (same for every variable, independent of how that variable's
  ! interpolant uses it downstream): the query time query_hr always refers to
  ! the END of the current ELM step, i.e. [elm_step_start, query_hr]. The
  ! evaluation time eval_hr = query_hr - elm_dt/2 (the ELM step's midpoint)
  ! is used for bracketing -- this applies uniformly to both interpolants,
  ! not just the linear one. The data's own dt is always >= elm_dt, so
  ! eval_hr is always bracketed by the data interval it falls in: tindex(1)
  ! ("previous") is the last data slice whose time is < eval_hr, and
  ! tindex(2) ("next") is the first data slice whose time is >= eval_hr --
  ! an exact match is treated as "next" (right-continuous convention). This
  ! single bracket supports both a linear interpolant (wt1/wt2 blend between
  ! the two) and a step/right-continuous-piecewise-constant interpolant
  ! (which just uses tindex(2)'s value outright, wt2=1) -- record selection
  ! never differs by field, only how the two records (and their weights)
  ! are combined downstream.
  !===============================================================================
  subroutine cpl_bypass_update_tindex(atm2lnd_vars, yr, mon, day, tod, dtime, const_climate_hist)
    type(atm2lnd_type), intent(inout) :: atm2lnd_vars
    integer,             intent(in)   :: yr, mon, day, tod, dtime
    logical,             intent(in)   :: const_climate_hist

    character(len=*), parameter :: sub = 'cpl_bypass_update_tindex'
    integer, parameter :: caldaym(13) = &
         (/ 1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 /)
    real(r8) :: query_hr, eval_hr, mystart, record_hr

    mystart = real(atm2lnd_vars%startyear_met, r8)
    query_hr = (real(yr,r8) - mystart) * 8760._r8 &
             + real(caldaym(mon) + day - 2, r8) * 24._r8 &
             + real(tod, r8) / 3600._r8
    eval_hr = query_hr - real(dtime, r8) / 3600._r8 / 2._r8

    ! Record k (1-based) is stored at time (k-1)*timeres. "next" is the first
    ! record whose time >= eval_hr (an exact match counts as "next", per the
    ! right-continuous convention); "previous" is the record immediately
    ! before it, so eval_hr is always in (previous_time, next_time].
    atm2lnd_vars%tindex(2) = ceiling(eval_hr / atm2lnd_vars%timeres) + 1
    atm2lnd_vars%tindex(1) = atm2lnd_vars%tindex(2) - 1

    ! Linear-interpolant weight: how far eval_hr sits between tindex(1)'s and
    ! tindex(2)'s own times, computed from the PRE-WRAP tindex(1) (directly
    ! consistent with eval_hr's own unwrapped timeline) -- not the wrapped
    ! absolute record index below, which for a const_climate_hist run whose
    ! synthetic clock year is far from startyear_met (e.g. spinup year 1
    ! against a forcing file stamped startyear_met=2009) sits nowhere near
    ! eval_hr, producing a huge bogus weight instead of one in [0,1]. Step-
    ! function fields ignore this and use wt2=1 (tindex(2)'s value outright)
    ! instead -- see lnd_import_cpl_bypass.
    record_hr = real(atm2lnd_vars%tindex(1) - 1, r8) * atm2lnd_vars%timeres
    atm2lnd_vars%interp_wt2 = (eval_hr - record_hr) / atm2lnd_vars%timeres
    atm2lnd_vars%interp_wt1 = 1._r8 - atm2lnd_vars%interp_wt2

    if (const_climate_hist) then
      atm2lnd_vars%tindex(1) = modulo(atm2lnd_vars%tindex(1) - 1, atm2lnd_vars%timelen_spinup) + 1
      atm2lnd_vars%tindex(2) = modulo(atm2lnd_vars%tindex(2) - 1, atm2lnd_vars%timelen_spinup) + 1
    else
      if (atm2lnd_vars%tindex(1) .lt. 1 .or. atm2lnd_vars%tindex(1) .gt. atm2lnd_vars%timelen) then
        call endrun(msg=' ERROR: '//sub//': run has advanced past the end of the met forcing record.'// &
                        ' If running a spinup, set const_climate_hist=.true. in user_nl_elm to cycle the met record.')
      end if
      if (atm2lnd_vars%tindex(2) .gt. atm2lnd_vars%timelen) then
        atm2lnd_vars%tindex(2) = atm2lnd_vars%tindex(1)
      end if
    end if

  end subroutine cpl_bypass_update_tindex

  !===============================================================================
  ! cpl_bypass_buffer_offset: convert an absolute record index (a tindex(1) or tindex(2)
  ! value) into its offset within the current page's atm_input buffer.
  ! Handles the spinup-wrap page (page_start > page_end), where the buffer
  ! holds [page_start .. record_len] in the first page_first_part offsets
  ! followed by [1 .. page_end] in the remaining offsets -- see
  ! cpl_bypass_load_page for how the page is populated.
  !===============================================================================
  pure function cpl_bypass_buffer_offset(atm2lnd_vars, abs_index) result(off)
    type(atm2lnd_type), intent(in) :: atm2lnd_vars
    integer,             intent(in) :: abs_index
    integer :: off

    if (atm2lnd_vars%page_start .le. atm2lnd_vars%page_end) then
      ! Ordinary, non-wrapped page.
      off = abs_index - atm2lnd_vars%page_start + 1
    else
      ! Wrapped page: indices >= page_start belong to the first (tail) part;
      ! indices <= page_end belong to the second (head) part.
      if (abs_index .ge. atm2lnd_vars%page_start) then
        off = abs_index - atm2lnd_vars%page_start + 1
      else
        off = atm2lnd_vars%page_first_part + abs_index
      end if
    end if
  end function cpl_bypass_buffer_offset

  !===============================================================================
  ! cpl_bypass_index_in_page: true if the given absolute record index is currently
  ! covered by the page, correctly handling the spinup-wrap case where
  ! page_start > page_end (the page covers [page_start..record_len] and
  ! [1..page_end], not a contiguous [page_start..page_end] range).
  !===============================================================================
  pure function cpl_bypass_index_in_page(atm2lnd_vars, abs_index) result(covered)
    type(atm2lnd_type), intent(in) :: atm2lnd_vars
    integer,             intent(in) :: abs_index
    logical :: covered

    if (atm2lnd_vars%page_start .le. atm2lnd_vars%page_end) then
      covered = (abs_index .ge. atm2lnd_vars%page_start) .and. (abs_index .le. atm2lnd_vars%page_end)
    else
      covered = (abs_index .ge. atm2lnd_vars%page_start) .or. (abs_index .le. atm2lnd_vars%page_end)
    end if
  end function cpl_bypass_index_in_page

  !===============================================================================
  ! cpl_bypass_ensure_page: lazily reload the page if the current bracketing
  ! indices fall outside [page_start, page_end]. Pages overlap by one time
  ! slice (page_end of the old page == page_start of the new page's first
  ! usable interpolation partner) so interpolation across a page boundary
  ! never needs a second read within one timestep.
  !===============================================================================
  subroutine cpl_bypass_ensure_page(bounds, atm2lnd_vars, met_nvars, metvars, metsource_str, &
                                  use_livneh, use_daymet, const_climate_hist)
    type(bounds_type),  intent(in)    :: bounds
    type(atm2lnd_type), intent(inout) :: atm2lnd_vars
    integer,             intent(in)   :: met_nvars
    character(len=20),   intent(in)   :: metvars(14)
    character(len=200),  intent(in)   :: metsource_str
    logical,             intent(in)   :: use_livneh, use_daymet
    logical,             intent(in)   :: const_climate_hist

    logical :: need_reload

    need_reload = (.not. cpl_bypass_index_in_page(atm2lnd_vars, atm2lnd_vars%tindex(1))) .or. &
                  (.not. cpl_bypass_index_in_page(atm2lnd_vars, atm2lnd_vars%tindex(2)))

    if (.not. need_reload) return

    call cpl_bypass_load_page(bounds, atm2lnd_vars, met_nvars, metvars, metsource_str, use_livneh, use_daymet, &
                            atm2lnd_vars%tindex(1))

  end subroutine cpl_bypass_ensure_page

  !===============================================================================
  ! cpl_bypass_load_page: read a page of atm_input starting at absolute record
  ! index start_tidx, sized atm2lnd_vars%page_ntimes (records), reading only
  ! the ncols_needed unique (zone, file-row) columns actually referenced by
  ! this rank's gridcells (atm2lnd_vars%col_gtoget/col_ztoget, populated once
  ! by cpl_bypass_build_column_index) -- never one read per gridcell.
  !
  ! For const_climate_hist spinup, if the requested page would run past the
  ! end of the spinup window, it wraps modulo timelen_spinup: the page reads
  ! the tail of the record and then the head, in one load, so interpolation
  ! across the wrap never needs a second reload right at the boundary.
  ! Pages otherwise overlap the next page's first record by one slice.
  !===============================================================================
  subroutine cpl_bypass_load_page(bounds, atm2lnd_vars, met_nvars, metvars, metsource_str, &
                                use_livneh, use_daymet, start_tidx)
    use netcdf
    use elm_varctl, only: metdata_bypass, const_climate_hist
    type(bounds_type),  intent(in)    :: bounds
    type(atm2lnd_type), intent(inout) :: atm2lnd_vars
    integer,             intent(in)   :: met_nvars
    character(len=20),   intent(in)   :: metvars(14)
    character(len=200),  intent(in)   :: metsource_str
    logical,             intent(in)   :: use_livneh, use_daymet
    integer,             intent(in)   :: start_tidx

    character(len=*), parameter :: sub = 'cpl_bypass_load_page'
    integer :: v, col, gtoget, ztoget, ierr, varid, ncid
    integer :: starti(2), counti(2)
    integer :: nbuf, first_part, second_part
    integer :: record_len
    integer*2, allocatable :: met_buf(:,:)
    character(len=256) :: metdata_fname

    record_len = merge(atm2lnd_vars%timelen_spinup, atm2lnd_vars%timelen, const_climate_hist)
    nbuf = min(atm2lnd_vars%page_ntimes, record_len)

    ! Spinup wrap: if the requested page runs past the record, it wraps --
    ! the page holds [start_tidx .. record_len] (buffer offsets 1..first_part)
    ! followed by [1 .. second_part] (buffer offsets first_part+1..nbuf), so
    ! page_end < page_start in absolute-record terms in that case. See
    ! cpl_bypass_buffer_offset() for how callers convert an absolute tindex to a buffer
    ! offset using page_start/page_first_part/record length. Only const_climate_hist
    ! actually wraps, though -- a plain transient run has nowhere to wrap TO (the
    ! record isn't cyclic), so a "whole record" page (page_ntimes==-1, nbuf==timelen)
    ! starting anywhere but record 1 must clamp nbuf down to what's actually left in
    ! the file instead of requesting a full-length page starting mid-file, which would
    ! run nbuf-first_part records past EOF -- an out-of-bounds nf90_get_var whose
    ! (unchecked) failure previously left atm_input holding uninitialized/garbage
    ! values for every later timestep of a transient run not starting at record 1.
    first_part  = min(nbuf, record_len - start_tidx + 1)
    second_part = nbuf - first_part
    if (.not. const_climate_hist) then
      nbuf = first_part
      second_part = 0
    end if
    atm2lnd_vars%page_first_part = first_part

    if (.not. associated(atm2lnd_vars%atm_input)) then
      allocate(atm2lnd_vars%atm_input(met_nvars, atm2lnd_vars%ncols_needed, 1:nbuf))
    end if

    atm2lnd_vars%page_start = start_tidx
    atm2lnd_vars%page_end = start_tidx + nbuf - 1
    if (const_climate_hist .and. second_part .gt. 0) then
      atm2lnd_vars%page_end = second_part
    end if

    allocate(met_buf(nbuf, 1))

    do v = 1, met_nvars
      do col = 1, atm2lnd_vars%ncols_needed
        gtoget = atm2lnd_vars%col_gtoget(col)
        ztoget = atm2lnd_vars%col_ztoget(col)

        metdata_fname = cpl_bypass_filename(atm2lnd_vars, metvars(v), metsource_str, use_livneh, use_daymet, ztoget)
        ierr = nf90_open(trim(metdata_bypass) // '/' // trim(metdata_fname), NF90_NOWRITE, ncid)
        if (ierr .ne. 0) call endrun(msg=' ERROR: '//sub//': Failed to open cpl_bypass input meteorology file: '//trim(metdata_fname))
        ierr = nf90_inq_varid(ncid, trim(metvars(v)), varid)

        if (const_climate_hist .and. second_part .gt. 0) then
          starti(1) = start_tidx; starti(2) = gtoget
          counti(1) = first_part; counti(2) = 1
          ierr = nf90_get_var(ncid, varid, met_buf(1:first_part,:), starti(1:2), counti(1:2))
          starti(1) = 1
          counti(1) = second_part
          ierr = nf90_get_var(ncid, varid, met_buf(first_part+1:nbuf,:), starti(1:2), counti(1:2))
        else
          starti(1) = start_tidx; starti(2) = gtoget
          counti(1) = nbuf; counti(2) = 1
          ierr = nf90_get_var(ncid, varid, met_buf, starti(1:2), counti(1:2))
        end if
        ierr = nf90_close(ncid)

        atm2lnd_vars%atm_input(v, col, :) = met_buf(:, 1)
      end do
    end do

    deallocate(met_buf)

  end subroutine cpl_bypass_load_page

  !===============================================================================
  ! cpl_bypass_read_range: read one variable/column's raw (still int16-packed)
  ! values for an arbitrary absolute record range [starti, endi] directly
  ! from its file. Used by the one-time site-bias-correction setup, which
  ! needs to scan a multi-year range that is normally far larger than the
  ! rolling page buffer -- so it reads independently of atm2lnd_vars%atm_input
  ! rather than assuming the page happens to cover the requested range.
  !===============================================================================
  subroutine cpl_bypass_read_range(atm2lnd_vars, metvar, metsource_str, use_livneh, use_daymet, &
                                 ztoget, gtoget, starti_rec, endi_rec, buf)
    use netcdf
    use elm_varctl, only: metdata_bypass
    type(atm2lnd_type), intent(in)  :: atm2lnd_vars
    character(len=*),   intent(in)  :: metvar, metsource_str
    logical,             intent(in) :: use_livneh, use_daymet
    integer,             intent(in) :: ztoget, gtoget, starti_rec, endi_rec
    integer*2, allocatable, intent(out) :: buf(:)

    character(len=*), parameter :: sub = 'cpl_bypass_read_range'
    character(len=256) :: metdata_fname
    integer :: ierr, varid, ncid, starti(2), counti(2)

    allocate(buf(endi_rec - starti_rec + 1))
    metdata_fname = cpl_bypass_filename(atm2lnd_vars, metvar, metsource_str, use_livneh, use_daymet, ztoget)
    ierr = nf90_open(trim(metdata_bypass) // '/' // trim(metdata_fname), NF90_NOWRITE, ncid)
    if (ierr .ne. 0) call endrun(msg=' ERROR: '//sub//': Failed to open cpl_bypass input meteorology file: '//trim(metdata_fname))
    ierr = nf90_inq_varid(ncid, trim(metvar), varid)
    starti(1) = starti_rec; starti(2) = gtoget
    counti(1) = endi_rec - starti_rec + 1; counti(2) = 1
    ierr = nf90_get_var(ncid, varid, buf, starti(1:2), counti(1:2))
    ierr = nf90_close(ncid)

  end subroutine cpl_bypass_read_range

  !===============================================================================
  ! finalize_forcing_common: downstream logic shared by both lnd_import_cpl_bypass
  ! and lnd_import_coupler once forc_*_grc (and the CO2 ppmv values) have been
  ! populated for gridcell g -- CO2 partial pressure conversion, CH4, derived
  ! humidity/density/wind-speed/solar quantities, and glc/iac coupling.
  !===============================================================================
  subroutine finalize_forcing_common(bounds, g, i, x2l, atm2lnd_vars, glc2lnd_vars, iac2lnd_vars, co2_ppmv_val, &
                                    co2_ppmv_prog, co2_ppmv_diag, co2_type_idx, use_c13, forc_rainc, forc_rainl, &
                                    forc_snowc, forc_snowl)
    use elm_varctl, only: create_glacier_mec_landunit
    use elm_varcon, only: rair, o2_molar_const, c13ratio
    use shr_const_mod, only: SHR_CONST_TKFRZ
    type(bounds_type),  intent(in)    :: bounds
    integer,             intent(in)   :: g, i
    real(r8),            intent(in)   :: x2l(:,:)
    type(atm2lnd_type), intent(inout) :: atm2lnd_vars
    type(glc2lnd_type), intent(inout) :: glc2lnd_vars
    type(iac2lnd_type), intent(inout) :: iac2lnd_vars
    real(r8),            intent(inout):: co2_ppmv_val
    real(r8),            intent(in)   :: co2_ppmv_prog, co2_ppmv_diag
    integer,             intent(in)   :: co2_type_idx
    logical,             intent(in)   :: use_c13
    real(r8),            intent(in)   :: forc_rainc, forc_rainl, forc_snowc, forc_snowl

    integer :: topo, num
    real(r8) :: forc_t, forc_q, forc_pbot, e, qsat

    do topo = grc_pp%topi(g), grc_pp%topf(g)
      top_as%pco2bot(topo) = co2_ppmv_val * 1.e-6_r8 * top_as%pbot(topo)
      if (use_c13) then
        top_as%pc13o2bot(topo) = top_as%pco2bot(topo) * c13ratio
      end if
    end do
    atm2lnd_vars%forc_pco2_grc(g) = co2_ppmv_val * 1.e-6_r8 * atm2lnd_vars%forc_pbot_not_downscaled_grc(g)

    if (index_x2l_Sa_methane /= 0) then
      do topo = grc_pp%topi(g), grc_pp%topf(g)
        top_as%pch4bot(topo) = x2l(index_x2l_Sa_methane,i)
      end do
      atm2lnd_vars%forc_pch4_grc(g) = x2l(index_x2l_Sa_methane,i)
    endif

    forc_t    = atm2lnd_vars%forc_t_not_downscaled_grc(g)
    forc_q    = atm2lnd_vars%forc_q_not_downscaled_grc(g)
    forc_pbot = atm2lnd_vars%forc_pbot_not_downscaled_grc(g)

    atm2lnd_vars%forc_hgt_u_grc(g) = atm2lnd_vars%forc_hgt_grc(g)
    atm2lnd_vars%forc_hgt_t_grc(g) = atm2lnd_vars%forc_hgt_grc(g)
    atm2lnd_vars%forc_hgt_q_grc(g) = atm2lnd_vars%forc_hgt_grc(g)
    atm2lnd_vars%forc_vp_grc(g)    = forc_q * forc_pbot / (0.622_r8 + 0.378_r8 * forc_q)
    atm2lnd_vars%forc_rho_not_downscaled_grc(g) = &
         (forc_pbot - 0.378_r8 * atm2lnd_vars%forc_vp_grc(g)) / (rair * forc_t)
    atm2lnd_vars%forc_po2_grc(g)   = o2_molar_const * forc_pbot
    atm2lnd_vars%forc_wind_grc(g)  = sqrt(atm2lnd_vars%forc_u_grc(g)**2 + atm2lnd_vars%forc_v_grc(g)**2)
    atm2lnd_vars%forc_solar_grc(g) = atm2lnd_vars%forc_solad_grc(g,1) + atm2lnd_vars%forc_solai_grc(g,1) + &
                                     atm2lnd_vars%forc_solad_grc(g,2) + atm2lnd_vars%forc_solai_grc(g,2)

    atm2lnd_vars%forc_rain_not_downscaled_grc(g) = forc_rainc + forc_rainl
    atm2lnd_vars%forc_snow_not_downscaled_grc(g) = forc_snowc + forc_snowl

    if (forc_t > SHR_CONST_TKFRZ) then
      e = esatw(tdc(forc_t))
    else
      e = esati(tdc(forc_t))
    end if
    qsat = 0.622_r8*e / (forc_pbot - 0.378_r8*e)
    atm2lnd_vars%forc_rh_grc(g) = 100.0_r8*(forc_q / qsat)

    if (co2_type_idx == 1) then
      co2_ppmv_val = co2_ppmv_prog
    else if (co2_type_idx == 2) then
      co2_ppmv_val = co2_ppmv_diag
      if (use_c13) then
        atm2lnd_vars%forc_pc13o2_grc(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * forc_pbot
      end if
    else
      if (use_c13) then
        atm2lnd_vars%forc_pc13o2_grc(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * forc_pbot
      end if
    end if
    atm2lnd_vars%forc_pco2_grc(g) = co2_ppmv_val * 1.e-6_r8 * forc_pbot

    if (create_glacier_mec_landunit) then
      do num = 0,glc_nec
        glc2lnd_vars%frac_grc(g,num)  = x2l(index_x2l_Sg_frac(num),i)
        glc2lnd_vars%topo_grc(g,num)  = x2l(index_x2l_Sg_topo(num),i)
        glc2lnd_vars%hflx_grc(g,num)  = x2l(index_x2l_Flgg_hflx(num),i)
      end do
      glc2lnd_vars%icemask_grc(g)  = x2l(index_x2l_Sg_icemask,i)
      glc2lnd_vars%icemask_coupled_fluxes_grc(g)  = x2l(index_x2l_Sg_icemask_coupled_fluxes,i)
    end if

    if (iac_present) then
      do num = 0,numpft
        iac2lnd_vars%frac_pft(g,num) = x2l(index_x2l_Sz_pct_pft(num),i)
        iac2lnd_vars%frac_pft_prev(g,num) = x2l(index_x2l_Sz_pct_pft_prev(num),i)
        if (num < numharvest) then
          iac2lnd_vars%harvest_frac(g,num) = x2l(index_x2l_Sz_harvest_frac(num),i)
        end if
      end do
    endif

  end subroutine finalize_forcing_common

#else
  !===============================================================================
  ! lnd_import_coupler: forcing taken directly from the CIME coupler's x2l
  ! array. This is today's non-CPL_BYPASS path, unchanged in substance.
  !===============================================================================
  subroutine lnd_import_coupler( bounds, x2l, atm2lnd_vars, glc2lnd_vars, ocn2lnd_vars, lnd2atm_vars, iac2lnd_vars)
    use elm_varctl       , only: co2_type, co2_ppmv, use_c13, create_glacier_mec_landunit, use_atm_downscaling_to_topunit
    use elm_varcon       , only: rair, o2_molar_const, c13ratio
    use shr_const_mod    , only: SHR_CONST_TKFRZ
    use domainMod        , only: ldomain
    use FrictionVelocityMod, only: implicit_stress, atm_gustiness
    use lnd_downscale_atm_forcing
    !
    type(bounds_type)  , intent(in)    :: bounds
    real(r8)           , intent(in)    :: x2l(:,:)
    type(atm2lnd_type) , intent(inout) :: atm2lnd_vars
    type(glc2lnd_type) , intent(inout) :: glc2lnd_vars
    type(ocn2lnd_type) , intent(inout) :: ocn2lnd_vars
    type(lnd2atm_type) , intent(in)    :: lnd2atm_vars
    type(iac2lnd_type) , intent(inout) :: iac2lnd_vars
    !
    integer  :: g, topo, i, co2_type_idx
    real(r8) :: forc_rainc, forc_rainl, forc_snowc, forc_snowl
    real(r8) :: e, qsat, vp
    real(r8) :: co2_ppmv_diag, co2_ppmv_prog, co2_ppmv_val
    character(len=*), parameter :: sub = 'lnd_import_coupler'

    co2_type_idx = 0
    if (co2_type == 'prognostic') then
       co2_type_idx = 1
    else if (co2_type == 'diagnostic') then
       co2_type_idx = 2
    end if
    if (co2_type == 'prognostic' .and. index_x2l_Sa_co2prog == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2prog for co2_type equal to prognostic' )
    else if (co2_type == 'diagnostic' .and. index_x2l_Sa_co2diag == 0) then
       call endrun( sub//' ERROR: must have nonzero index_x2l_Sa_co2diag for co2_type equal to diagnostic' )
    end if

    ! Note that the precipitation fluxes received from the coupler are in
    ! units of kg/s/m^2. To convert to mm/sec, divide by 1000 kg/m^3 and
    ! multiply by 1000 mm/m -- an overall factor of unity. Units below are
    ! therefore given in mm/s.

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)

       atm2lnd_vars%forc_flood_grc(g)   = -x2l(index_x2l_Flrr_flood,i)
       atm2lnd_vars%volr_grc(g)   = x2l(index_x2l_Flrr_volr,i) * (ldomain%area(g) * 1.e6_r8)
       atm2lnd_vars%volrmch_grc(g)= x2l(index_x2l_Flrr_volrmch,i) * (ldomain%area(g) * 1.e6_r8)
       atm2lnd_vars%supply_grc(g) = x2l(index_x2l_Flrr_supply,i)
       atm2lnd_vars%deficit_grc(g) = x2l(index_x2l_Flrr_deficit,i)

       if (index_x2l_Sr_h2orof /= 0) then
         atm2lnd_vars%h2orof_grc(g)      = x2l(index_x2l_Sr_h2orof,i)
         atm2lnd_vars%frac_h2orof_grc(g) = x2l(index_x2l_Sr_frac_h2orof,i)
       endif

       if (index_x2l_So_ssh /= 0) then
         ocn2lnd_vars%ssh_grc(g) = x2l(index_x2l_So_ssh,i)
       endif

       if (index_x2l_So_frac_h2oocn /= 0) then
         ocn2lnd_vars%frac_h2oocn_grc(g) = x2l(index_x2l_So_frac_h2oocn,i)
       endif

       atm2lnd_vars%forc_hgt_grc(g)     = x2l(index_x2l_Sa_z,i)
       atm2lnd_vars%forc_u_grc(g)       = x2l(index_x2l_Sa_u,i)
       atm2lnd_vars%forc_v_grc(g)       = x2l(index_x2l_Sa_v,i)
       atm2lnd_vars%forc_solad_grc(g,2) = x2l(index_x2l_Faxa_swndr,i)
       atm2lnd_vars%forc_solad_grc(g,1) = x2l(index_x2l_Faxa_swvdr,i)
       atm2lnd_vars%forc_solai_grc(g,2) = x2l(index_x2l_Faxa_swndf,i)
       atm2lnd_vars%forc_solai_grc(g,1) = x2l(index_x2l_Faxa_swvdf,i)

       atm2lnd_vars%forc_th_not_downscaled_grc(g)    = x2l(index_x2l_Sa_ptem,i)
       atm2lnd_vars%forc_q_not_downscaled_grc(g)     = x2l(index_x2l_Sa_shum,i)
       atm2lnd_vars%forc_pbot_not_downscaled_grc(g)  = x2l(index_x2l_Sa_pbot,i)
       atm2lnd_vars%forc_t_not_downscaled_grc(g)     = x2l(index_x2l_Sa_tbot,i)
       atm2lnd_vars%forc_lwrad_not_downscaled_grc(g) = x2l(index_x2l_Faxa_lwdn,i)

       forc_rainc = x2l(index_x2l_Faxa_rainc,i)
       forc_rainl = x2l(index_x2l_Faxa_rainl,i)
       forc_snowc = x2l(index_x2l_Faxa_snowc,i)
       forc_snowl = x2l(index_x2l_Faxa_snowl,i)

       atm2lnd_vars%forc_aer_grc(g,1)  =  x2l(index_x2l_Faxa_bcphidry,i)
       atm2lnd_vars%forc_aer_grc(g,2)  =  x2l(index_x2l_Faxa_bcphodry,i)
       atm2lnd_vars%forc_aer_grc(g,3)  =  x2l(index_x2l_Faxa_bcphiwet,i)
       atm2lnd_vars%forc_aer_grc(g,4)  =  x2l(index_x2l_Faxa_ocphidry,i)
       atm2lnd_vars%forc_aer_grc(g,5)  =  x2l(index_x2l_Faxa_ocphodry,i)
       atm2lnd_vars%forc_aer_grc(g,6)  =  x2l(index_x2l_Faxa_ocphiwet,i)
       atm2lnd_vars%forc_aer_grc(g,7)  =  x2l(index_x2l_Faxa_dstwet1,i)
       atm2lnd_vars%forc_aer_grc(g,8)  =  x2l(index_x2l_Faxa_dstdry1,i)
       atm2lnd_vars%forc_aer_grc(g,9)  =  x2l(index_x2l_Faxa_dstwet2,i)
       atm2lnd_vars%forc_aer_grc(g,10) =  x2l(index_x2l_Faxa_dstdry2,i)
       atm2lnd_vars%forc_aer_grc(g,11) =  x2l(index_x2l_Faxa_dstwet3,i)
       atm2lnd_vars%forc_aer_grc(g,12) =  x2l(index_x2l_Faxa_dstdry3,i)
       atm2lnd_vars%forc_aer_grc(g,13) =  x2l(index_x2l_Faxa_dstwet4,i)
       atm2lnd_vars%forc_aer_grc(g,14) =  x2l(index_x2l_Faxa_dstdry4,i)

       !set the topounit-level atmospheric state and flux forcings
       if (use_atm_downscaling_to_topunit) then
          if(atm_gustiness) then
             call endrun("Error: atm_gustiness not yet supported with multiple topounits")
          end if
         do topo = grc_pp%topi(g) , grc_pp%topf(g)
            top_as%ugust(topo) = 0._r8
         end do

         call downscale_atm_forcing_to_topounit(g, i, x2l, lnd2atm_vars)
       else
         do topo = grc_pp%topi(g), grc_pp%topf(g)
           top_as%tbot(topo)    = x2l(index_x2l_Sa_tbot,i)
           top_as%thbot(topo)   = x2l(index_x2l_Sa_ptem,i)
           top_as%pbot(topo)    = x2l(index_x2l_Sa_pbot,i)
           top_as%qbot(topo)    = x2l(index_x2l_Sa_shum,i)
           top_as%ubot(topo)    = x2l(index_x2l_Sa_u,i)
           top_as%vbot(topo)    = x2l(index_x2l_Sa_v,i)
           top_as%zbot(topo)    = x2l(index_x2l_Sa_z,i)
         if (implicit_stress) then
            top_as%wsresp(topo)  = x2l(index_x2l_Sa_wsresp,i)
            top_as%tau_est(topo) = x2l(index_x2l_Sa_tau_est,i)
         end if
         if (atm_gustiness) then
            top_as%ugust(topo)  = x2l(index_x2l_Sa_ugust,i)
         else
            top_as%ugust(topo) = 0._r8
         end if
           top_as%windbot(topo) = sqrt(top_as%ubot(topo)**2 + top_as%vbot(topo)**2)
           if (atm_gustiness) then
              top_as%windbot(topo) = sqrt(top_as%windbot(topo)**2 + top_as%ugust(topo)**2)
           end if
           if (top_as%tbot(topo) > SHR_CONST_TKFRZ) then
            e = esatw(tdc(top_as%tbot(topo)))
           else
            e = esati(tdc(top_as%tbot(topo)))
           end if
           qsat = 0.622_r8*e / (top_as%pbot(topo) - 0.378_r8*e)
           top_as%rhbot(topo) = 100.0_r8*(top_as%qbot(topo) / qsat)
           top_as%po2bot(topo) = o2_molar_const * top_as%pbot(topo)
           vp = top_as%qbot(topo) * top_as%pbot(topo)  / (0.622_r8 + 0.378_r8 * top_as%qbot(topo))
           top_as%rhobot(topo) = (top_as%pbot(topo) - 0.378_r8 * vp) / (rair * top_as%tbot(topo))

           top_af%rain(topo)    = forc_rainc + forc_rainl
           top_af%snow(topo)    = forc_snowc + forc_snowl
           top_af%solad(topo,2) = x2l(index_x2l_Faxa_swndr,i)
           top_af%solad(topo,1) = x2l(index_x2l_Faxa_swvdr,i)
           top_af%solai(topo,2) = x2l(index_x2l_Faxa_swndf,i)
           top_af%solai(topo,1) = x2l(index_x2l_Faxa_swvdf,i)
           top_af%lwrad(topo)   = x2l(index_x2l_Faxa_lwdn,i)
           top_af%solar(topo) = top_af%solad(topo,2) + top_af%solad(topo,1) + &
                              top_af%solai(topo,2) + top_af%solai(topo,1)
         end do
       end if

       ! Determine optional receive fields
       ! CO2 (and C13O2) concentration: constant, prognostic, or diagnostic
       if (co2_type_idx == 0) then
         co2_ppmv_val = co2_ppmv
       else if (co2_type_idx == 1) then
         co2_ppmv_val = x2l(index_x2l_Sa_co2prog,i)
       else if (co2_type_idx == 2) then
         co2_ppmv_val = x2l(index_x2l_Sa_co2diag,i)
       else
         call endrun( sub//' ERROR: Invalid co2_type_idx, must be 0, 1, or 2 (constant, prognostic, or diagnostic)' )
       end if

      if (index_x2l_Sa_co2prog /= 0) then
         co2_ppmv_prog = x2l(index_x2l_Sa_co2prog,i)
      else
         co2_ppmv_prog = co2_ppmv
      end if

      if (index_x2l_Sa_co2diag /= 0) then
         co2_ppmv_diag = x2l(index_x2l_Sa_co2diag,i)
      else
         co2_ppmv_diag = co2_ppmv
      end if

      call finalize_forcing_common(bounds, g, i, x2l, atm2lnd_vars, glc2lnd_vars, iac2lnd_vars, co2_ppmv_val, &
                                  co2_ppmv_prog, co2_ppmv_diag, co2_type_idx, use_c13, forc_rainc, forc_rainl, &
                                  forc_snowc, forc_snowl)

       if (index_x2l_Sa_methane /= 0) then
          atm2lnd_vars%forc_pch4_grc(g) = x2l(index_x2l_Sa_methane,i)
       endif

    end do

  end subroutine lnd_import_coupler

  !===============================================================================
  ! finalize_forcing_common: downstream logic shared by both lnd_import_cpl_bypass
  ! and lnd_import_coupler once forc_*_grc (and the CO2 ppmv values) have been
  ! populated for gridcell g -- CO2 partial pressure conversion, CH4, derived
  ! humidity/density/wind-speed/solar quantities, and glc/iac coupling.
  !===============================================================================
  subroutine finalize_forcing_common(bounds, g, i, x2l, atm2lnd_vars, glc2lnd_vars, iac2lnd_vars, co2_ppmv_val, &
                                    co2_ppmv_prog, co2_ppmv_diag, co2_type_idx, use_c13, forc_rainc, forc_rainl, &
                                    forc_snowc, forc_snowl)
    use elm_varctl, only: create_glacier_mec_landunit
    use elm_varcon, only: rair, o2_molar_const, c13ratio
    use shr_const_mod, only: SHR_CONST_TKFRZ
    type(bounds_type),  intent(in)    :: bounds
    integer,             intent(in)   :: g, i
    real(r8),            intent(in)   :: x2l(:,:)
    type(atm2lnd_type), intent(inout) :: atm2lnd_vars
    type(glc2lnd_type), intent(inout) :: glc2lnd_vars
    type(iac2lnd_type), intent(inout) :: iac2lnd_vars
    real(r8),            intent(inout):: co2_ppmv_val
    real(r8),            intent(in)   :: co2_ppmv_prog, co2_ppmv_diag
    integer,             intent(in)   :: co2_type_idx
    logical,             intent(in)   :: use_c13
    real(r8),            intent(in)   :: forc_rainc, forc_rainl, forc_snowc, forc_snowl

    integer :: topo, num
    real(r8) :: forc_t, forc_q, forc_pbot, e, qsat

    do topo = grc_pp%topi(g), grc_pp%topf(g)
      top_as%pco2bot(topo) = co2_ppmv_val * 1.e-6_r8 * top_as%pbot(topo)
      if (use_c13) then
        top_as%pc13o2bot(topo) = top_as%pco2bot(topo) * c13ratio
      end if
    end do
    atm2lnd_vars%forc_pco2_grc(g) = co2_ppmv_val * 1.e-6_r8 * atm2lnd_vars%forc_pbot_not_downscaled_grc(g)

    forc_t    = atm2lnd_vars%forc_t_not_downscaled_grc(g)
    forc_q    = atm2lnd_vars%forc_q_not_downscaled_grc(g)
    forc_pbot = atm2lnd_vars%forc_pbot_not_downscaled_grc(g)

    atm2lnd_vars%forc_hgt_u_grc(g) = atm2lnd_vars%forc_hgt_grc(g)
    atm2lnd_vars%forc_hgt_t_grc(g) = atm2lnd_vars%forc_hgt_grc(g)
    atm2lnd_vars%forc_hgt_q_grc(g) = atm2lnd_vars%forc_hgt_grc(g)
    atm2lnd_vars%forc_vp_grc(g)    = forc_q * forc_pbot / (0.622_r8 + 0.378_r8 * forc_q)
    atm2lnd_vars%forc_rho_not_downscaled_grc(g) = &
         (forc_pbot - 0.378_r8 * atm2lnd_vars%forc_vp_grc(g)) / (rair * forc_t)
    atm2lnd_vars%forc_po2_grc(g)   = o2_molar_const * forc_pbot
    atm2lnd_vars%forc_wind_grc(g)  = sqrt(atm2lnd_vars%forc_u_grc(g)**2 + atm2lnd_vars%forc_v_grc(g)**2)
    atm2lnd_vars%forc_solar_grc(g) = atm2lnd_vars%forc_solad_grc(g,1) + atm2lnd_vars%forc_solai_grc(g,1) + &
                                     atm2lnd_vars%forc_solad_grc(g,2) + atm2lnd_vars%forc_solai_grc(g,2)

    atm2lnd_vars%forc_rain_not_downscaled_grc(g) = forc_rainc + forc_rainl
    atm2lnd_vars%forc_snow_not_downscaled_grc(g) = forc_snowc + forc_snowl

    if (forc_t > SHR_CONST_TKFRZ) then
      e = esatw(tdc(forc_t))
    else
      e = esati(tdc(forc_t))
    end if
    qsat = 0.622_r8*e / (forc_pbot - 0.378_r8*e)
    atm2lnd_vars%forc_rh_grc(g) = 100.0_r8*(forc_q / qsat)

    if (co2_type_idx == 1) then
      co2_ppmv_val = co2_ppmv_prog
    else if (co2_type_idx == 2) then
      co2_ppmv_val = co2_ppmv_diag
      if (use_c13) then
        atm2lnd_vars%forc_pc13o2_grc(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * forc_pbot
      end if
    else
      if (use_c13) then
        atm2lnd_vars%forc_pc13o2_grc(g) = co2_ppmv_val * c13ratio * 1.e-6_r8 * forc_pbot
      end if
    end if
    atm2lnd_vars%forc_pco2_grc(g) = co2_ppmv_val * 1.e-6_r8 * forc_pbot

    if (create_glacier_mec_landunit) then
      do num = 0,glc_nec
        glc2lnd_vars%frac_grc(g,num)  = x2l(index_x2l_Sg_frac(num),i)
        glc2lnd_vars%topo_grc(g,num)  = x2l(index_x2l_Sg_topo(num),i)
        glc2lnd_vars%hflx_grc(g,num)  = x2l(index_x2l_Flgg_hflx(num),i)
      end do
      glc2lnd_vars%icemask_grc(g)  = x2l(index_x2l_Sg_icemask,i)
      glc2lnd_vars%icemask_coupled_fluxes_grc(g)  = x2l(index_x2l_Sg_icemask_coupled_fluxes,i)
    end if

    if (iac_present) then
      do num = 0,numpft
        iac2lnd_vars%frac_pft(g,num) = x2l(index_x2l_Sz_pct_pft(num),i)
        iac2lnd_vars%frac_pft_prev(g,num) = x2l(index_x2l_Sz_pct_pft_prev(num),i)
        if (num < numharvest) then
          iac2lnd_vars%harvest_frac(g,num) = x2l(index_x2l_Sz_harvest_frac(num),i)
        end if
      end do
    endif

  end subroutine finalize_forcing_common

#endif

  !===============================================================================

  subroutine lnd_export( bounds, lnd2atm_vars, lnd2glc_vars, lnd2iac_vars, l2x)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the data to be sent from the elm model to the coupler
    !
    ! !USES:
    use shr_kind_mod       , only : r8 => shr_kind_r8
    use elm_varctl         , only : iulog, create_glacier_mec_landunit, iac_present, use_finetop_rad
    use elm_time_manager   , only : get_nstep, get_step_size
    use domainMod          , only : ldomain
    use seq_drydep_mod     , only : n_drydep
    use shr_megan_mod      , only : shr_megan_mechcomps_n
    use FireMod            , only : forc_hdm
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in)    :: bounds  ! bounds
    type(lnd2atm_type), intent(inout) :: lnd2atm_vars ! elm land to atmosphere exchange data type
    type(lnd2glc_type), intent(inout) :: lnd2glc_vars ! elm land to atmosphere exchange data type
    type(lnd2iac_type), intent(inout) :: lnd2iac_vars ! elm lnd to gcam exchange vars
    real(r8)          , intent(out)   :: l2x(:,:)! land to coupler export state on land grid
    !
    ! !LOCAL VARIABLES:
    integer  :: g,i,p ! indices
    integer  :: ier   ! error status
    integer  :: nstep ! time step index
    integer  :: dtime ! time step
    integer  :: num   ! counter
    character(len=*), parameter :: sub = 'lnd_export_mct'
    !---------------------------------------------------------------------------

    dtime = get_step_size()

    ! cesm sign convention is that fluxes are positive downward

    l2x(:,:) = 0.0_r8

    do g = bounds%begg,bounds%endg
       i = 1 + (g-bounds%begg)
       l2x(index_l2x_Sl_t,i)        =  lnd2atm_vars%t_rad_grc(g)
       l2x(index_l2x_Sl_snowh,i)    =  lnd2atm_vars%h2osno_grc(g)

       if (use_finetop_rad) then
          l2x(index_l2x_Sl_avsdr,i)    =  lnd2atm_vars%apparent_albd_grc(g,1)
          l2x(index_l2x_Sl_anidr,i)    =  lnd2atm_vars%apparent_albd_grc(g,2)
          l2x(index_l2x_Sl_avsdf,i)    =  lnd2atm_vars%apparent_albi_grc(g,1)
          l2x(index_l2x_Sl_anidf,i)    =  lnd2atm_vars%apparent_albi_grc(g,2)
       else
          l2x(index_l2x_Sl_avsdr,i)    =  lnd2atm_vars%albd_grc(g,1)
          l2x(index_l2x_Sl_anidr,i)    =  lnd2atm_vars%albd_grc(g,2)
          l2x(index_l2x_Sl_avsdf,i)    =  lnd2atm_vars%albi_grc(g,1)
          l2x(index_l2x_Sl_anidf,i)    =  lnd2atm_vars%albi_grc(g,2)
       end if

       l2x(index_l2x_Sl_tref,i)     =  lnd2atm_vars%t_ref2m_grc(g)
       l2x(index_l2x_Sl_qref,i)     =  lnd2atm_vars%q_ref2m_grc(g)
       l2x(index_l2x_Sl_u10,i)      =  lnd2atm_vars%u_ref10m_grc(g)
       l2x(index_l2x_Sl_u10withgusts,i)=lnd2atm_vars%u_ref10m_with_gusts_grc(g)
       l2x(index_l2x_Fall_taux,i)   = -lnd2atm_vars%taux_grc(g)
       l2x(index_l2x_Fall_tauy,i)   = -lnd2atm_vars%tauy_grc(g)
       l2x(index_l2x_Fall_lat,i)    = -lnd2atm_vars%eflx_lh_tot_grc(g)
       l2x(index_l2x_Fall_sen,i)    = -lnd2atm_vars%eflx_sh_tot_grc(g)
       l2x(index_l2x_Fall_lwup,i)   = -lnd2atm_vars%eflx_lwrad_out_grc(g)
       l2x(index_l2x_Fall_evap,i)   = -lnd2atm_vars%qflx_evap_tot_grc(g)
       l2x(index_l2x_Fall_swnet,i)  =  lnd2atm_vars%fsa_grc(g)
       if (index_l2x_Fall_fco2_lnd /= 0) then
          l2x(index_l2x_Fall_fco2_lnd,i) = -lnd2atm_vars%nee_grc(g)
       end if

       ! Additional fields for DUST, PROGSSLT, dry-deposition and VOC
       ! These are now standard fields, but the check on the index makes sure the driver handles them
       if (index_l2x_Sl_ram1      /= 0 )  l2x(index_l2x_Sl_ram1,i)     =  lnd2atm_vars%ram1_grc(g)
       if (index_l2x_Sl_fv        /= 0 )  l2x(index_l2x_Sl_fv,i)       =  lnd2atm_vars%fv_grc(g)
       if (index_l2x_Sl_soilw     /= 0 )  l2x(index_l2x_Sl_soilw,i)    =  lnd2atm_vars%h2osoi_vol_grc(g,1)
       if (index_l2x_Fall_flxdst1 /= 0 )  l2x(index_l2x_Fall_flxdst1,i)= -lnd2atm_vars%flxdst_grc(g,1)
       if (index_l2x_Fall_flxdst2 /= 0 )  l2x(index_l2x_Fall_flxdst2,i)= -lnd2atm_vars%flxdst_grc(g,2)
       if (index_l2x_Fall_flxdst3 /= 0 )  l2x(index_l2x_Fall_flxdst3,i)= -lnd2atm_vars%flxdst_grc(g,3)
       if (index_l2x_Fall_flxdst4 /= 0 )  l2x(index_l2x_Fall_flxdst4,i)= -lnd2atm_vars%flxdst_grc(g,4)


       ! for dry dep velocities
       if (index_l2x_Sl_ddvel     /= 0 )  then
          l2x(index_l2x_Sl_ddvel:index_l2x_Sl_ddvel+n_drydep-1,i) = &
               lnd2atm_vars%ddvel_grc(g,:n_drydep)
       end if

       ! for MEGAN VOC emis fluxes
       if (index_l2x_Fall_flxvoc  /= 0 ) then
          l2x(index_l2x_Fall_flxvoc:index_l2x_Fall_flxvoc+shr_megan_mechcomps_n-1,i) = &
               -lnd2atm_vars%flxvoc_grc(g,:shr_megan_mechcomps_n)
       end if

       if (index_l2x_Fall_methane /= 0) then
          l2x(index_l2x_Fall_methane,i) = -lnd2atm_vars%flux_ch4_grc(g)
       endif

       ! for FAN NH3 fluxes
       if (index_l2x_Fall_flxnh3 /= 0) then
          l2x(index_l2x_Fall_flxnh3,i) = -lnd2atm_vars%flux_nh3_grc(g)
       end if

       ! sign convention is positive downward with
       ! hierarchy of atm/glc/lnd/rof/ice/ocn.  so water sent from land to rof is positive

       l2x(index_l2x_Flrl_rofi,i) = lnd2atm_vars%qflx_rofice_grc(g)
       l2x(index_l2x_Flrl_rofsur,i) = lnd2atm_vars%qflx_rofliq_qsur_grc(g) &
                                    + lnd2atm_vars%qflx_rofliq_qsurp_grc(g)   !  surface ponding
       l2x(index_l2x_Flrl_rofsub,i) = lnd2atm_vars%qflx_rofliq_qsub_grc(g) &
                                    + lnd2atm_vars%qflx_rofliq_qsubp_grc(g)   !  perched drainiage
       l2x(index_l2x_Flrl_rofgwl,i) = lnd2atm_vars%qflx_rofliq_qgwl_grc(g)

       l2x(index_l2x_Flrl_demand,i) =  lnd2atm_vars%qflx_irr_demand_grc(g)   ! needs to be filled in
       if (l2x(index_l2x_Flrl_demand,i) > 0.0_r8) then
           write(iulog,*)'lnd2atm_vars%qflx_irr_demand_grc is',lnd2atm_vars%qflx_irr_demand_grc(g)
           write(iulog,*)'l2x(index_l2x_Flrl_demand,i) is',l2x(index_l2x_Flrl_demand,i)
           call endrun( sub//' ERROR: demand must be <= 0.')
       endif
       l2x(index_l2x_Flrl_Tqsur,i)  = lnd2atm_vars%Tqsur_grc(g)
       l2x(index_l2x_Flrl_Tqsub,i)  = lnd2atm_vars%Tqsub_grc(g)
       l2x(index_l2x_coszen_str,i)  = lnd2atm_vars%coszen_str(g)
	   if (rof_sed) then
           l2x(index_l2x_Flrl_rofmud,i) = lnd2atm_vars%qflx_rofmud_grc(g)
	   end if
       l2x(index_l2x_Flrl_wslake,i) = lnd2atm_vars%wslake_grc(g)/dtime

       if (index_l2x_Flrl_inundinf /= 0) then
          l2x(index_l2x_Flrl_inundinf,i) = lnd2atm_vars%qflx_h2orof_drain_grc(g)
       endif

       ! glc coupling

       if (create_glacier_mec_landunit) then
          do num = 0,glc_nec
             l2x(index_l2x_Sl_tsrf(num),i)   = lnd2glc_vars%tsrf_grc(g,num)
             l2x(index_l2x_Sl_topo(num),i)   = lnd2glc_vars%topo_grc(g,num)
             l2x(index_l2x_Flgl_qice(num),i) = lnd2glc_vars%qice_grc(g,num)
          end do
       end if

       if (iac_present) then
          do p = 0,numpft
             l2x(index_l2x_Sl_hr(p),i)        = lnd2iac_vars%hr(g,p)
             l2x(index_l2x_Sl_npp(p),i)       = lnd2iac_vars%npp(g,p)
             l2x(index_l2x_Sl_pftwgt(p),i)    = lnd2iac_vars%pftwgt(g,p)
          end do
          ! Scalar per-gridcell fields
          ! forc_hdm is time-interpolated from stream file (FireMod in non-CPL_BYPASS,
          ! or mirrored from atm2lnd_vars in CPL_BYPASS during lnd_import);
          ! copy here into lnd2iac_vars and pack into the coupler vector.
          lnd2iac_vars%forc_hdm(g)          = forc_hdm(g)
          l2x(index_l2x_Sl_forc_hdm, i)     = lnd2iac_vars%forc_hdm(g)
          l2x(index_l2x_Sl_hdd,i)           = lnd2iac_vars%hdd(g)
          l2x(index_l2x_Sl_cdd,i)           = lnd2iac_vars%cdd(g)
       end if
    end do

  end subroutine lnd_export

end module lnd_import_export


double precision function szenith(xcoor, ycoor, ltm, jday, hr, min, offset)
  !Function to calcualte solar zenith angle
  !Used in coupler bypass mode to compute inerpolation for incoming solar

  use shr_kind_mod , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  implicit none
  !inputs
  real(r8) xcoor, ycoor, offset_min
  integer jday, hr, min, ltm, offset
  !working variables
  real(r8) d2r, r2d, lsn, latrad, decrad, decdeg, ha
  real(r8) hangle, harad, saltrad, saltdeg, sazirad, sazideg
  real(r8) szendeg,szenrad

  real pi
  parameter(pi = 3.14159265358979)
  offset_min = offset/60d0   !note assumes 1hr or smaller timestep
  min = min - offset_min

  !adjust time for offsets
  if (min < 0) then
    hr = hr - 1
    min = min+60
  end if
  if (min >= 60) then
    hr = hr+1
    min = min-60
  end if
  if (hr < 0) then
    hr = hr+24
    jday = jday-1
  end if
  if (hr >= 24) then
    hr = hr-24
    jday = jday+1
  end if

  if (jday < 1) jday = 1
  if (xcoor > 180d0) xcoor = xcoor-360d0

  d2r     = pi/180d0
  r2d     = 1/d2r
  lsn     = 12.0d0+((ltm-xcoor)/15.0d0)
  latrad  = ycoor*d2r
  decrad  = 23.45*d2r*sin(d2r*360d0*(284d0+jday)/365d0)
  decdeg  = decrad*r2d
  ha      = hr+min/60.0d0
  hangle  = (lsn-ha)*60.0d0
  harad   = hangle*0.0043633d0

  saltrad = asin((sin(latrad)*sin(decrad))+(cos(latrad)*cos(decrad) &
       *cos(harad)))
  saltdeg = saltrad * r2d
  sazirad = asin(cos(decrad)*sin(harad)/cos(saltrad))
  sazideg = sazirad * r2d

  if (saltdeg .lt. 0.0d0 .or. saltrad .gt. 180.0d0) then  ! sun is below horizon
     saltdeg = 0.0d0
     saltrad = 0.0d0
     szendeg = 90.0d0
     szenrad = 90.0d0*d2r
  else
     szendeg = 90d0-saltdeg
     szenrad = szendeg*d2r
  end if
  szenith = szendeg

end function szenith
