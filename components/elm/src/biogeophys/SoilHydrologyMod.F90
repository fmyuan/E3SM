module SoilHydrologyMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate soil hydrology
  !
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
   use elm_varctl        , only : iulog, use_vichydro, use_humhol, obs_zwt_forcing
  use elm_varctl        , only : use_lnd_rof_two_way, lnd_rof_coupling_nstep
  use elm_varctl        , only : use_modified_infil
  use elm_varcon        , only : e_ice, denh2o, denice, rpi, spval
  use EnergyFluxType    , only : energyflux_type
  use SoilHydrologyType , only : soilhydrology_type
  use SoilStateType     , only : soilstate_type
  use WaterfluxType     , only : waterflux_type
  use TopounitType      , only : top_pp
  use TopounitDataType  , only : top_ws
  use GridcellType      , only : grc_pp
  use LandunitType      , only : lun_pp
  use ColumnType        , only : col_pp
  use ColumnDataType    , only : col_es, col_ws, col_wf
  use VegetationType    , only : veg_pp
  use VegetationDataType, only : veg_wf
  use abortutils      , only : endrun

  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SurfaceRunoff        ! Calculate surface runoff
  public :: Infiltration         ! Calculate infiltration into surface soil layer
  public :: WaterTable           ! Calculate water table before imposing drainage
  public :: Drainage             ! Calculate subsurface drainage
  public :: DrainageVSFM         ! Calculate subsurface drainage for VSFM
  public :: ELMVICMap
  private :: rediagnose_zwt_from_liq_profile
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine rediagnose_zwt_from_liq_profile(nlevbed, use_var_soil_thick, &
       h2osoi_liq_col, eff_porosity_col, dzmm_col, zi_col, zwt_col, jwt_col)
    !
    ! !DESCRIPTION:
    ! Infer a shallower water table from a continuous saturated liquid profile.
    !
    use elm_varcon, only : watmin
    !
    ! !ARGUMENTS:
    integer , intent(in)    :: nlevbed
    logical , intent(in)    :: use_var_soil_thick
    real(r8), intent(in)    :: h2osoi_liq_col(:)
    real(r8), intent(in)    :: eff_porosity_col(:)
    real(r8), intent(in)    :: dzmm_col(:)
    real(r8), intent(in)    :: zi_col(0:)
    real(r8), intent(inout) :: zwt_col
    integer , intent(inout) :: jwt_col
    !
    ! !LOCAL VARIABLES:
    integer  :: j
    integer  :: k_sat
    real(r8) :: layer_capacity
    real(r8) :: sat_liq
    real(r8) :: sat_frac
    real(r8) :: zwt_layer_top
    real(r8) :: zwt_layer_bot
    real(r8) :: zwt_profile
    real(r8), parameter :: zwt_interp_min = 0.90_r8
    real(r8), parameter :: zwt_sat_threshold = 0.99_r8
    !-----------------------------------------------------------------------

    k_sat = nlevbed + 1
    do j = nlevbed, 1, -1
       layer_capacity = max(watmin, eff_porosity_col(j)*dzmm_col(j))
       sat_liq = h2osoi_liq_col(j) / layer_capacity
       if (sat_liq >= zwt_sat_threshold) then
          k_sat = j
       else
          exit
       end if
    end do

    if (k_sat == 1) then
       zwt_profile = 0._r8
    else
       j = min(k_sat-1, nlevbed)
       if (j >= 1) then
          layer_capacity = max(watmin, eff_porosity_col(j)*dzmm_col(j))
          sat_liq = h2osoi_liq_col(j) / layer_capacity
          if (sat_liq >= zwt_interp_min) then
             if (j == 1) then
                zwt_layer_top = 0._r8
             else
                zwt_layer_top = zi_col(j-1)
             end if
             zwt_layer_bot = zi_col(j)
             sat_frac = (sat_liq - zwt_interp_min) / (zwt_sat_threshold - zwt_interp_min)
             sat_frac = max(0._r8, min(1._r8, sat_frac))
             zwt_profile = zwt_layer_bot - sat_frac*(zwt_layer_bot - zwt_layer_top)
          else if (k_sat <= nlevbed) then
             zwt_profile = zi_col(k_sat-1)
          else
             zwt_profile = zwt_col
          end if
       else
          zwt_profile = zwt_col
       end if
    end if

    if (zwt_profile < zwt_col) then
       zwt_col = max(0._r8, zwt_profile)
    end if

    jwt_col = nlevbed
    do j = 1, nlevbed
       if (use_var_soil_thick) then
          if (zwt_col <= zi_col(j) .and. zwt_col < zi_col(nlevbed)) then
             jwt_col = j-1
             exit
          end if
       else
          if (zwt_col <= zi_col(j)) then
             jwt_col = j-1
             exit
          end if
       end if
    end do

  end subroutine rediagnose_zwt_from_liq_profile

  subroutine SurfaceRunoff (bounds, num_hydrologyc, filter_hydrologyc, &
       num_urbanc, filter_urbanc, soilhydrology_vars, soilstate_vars, dtime)
    !
    ! !DESCRIPTION:
    ! Calculate surface runoff
    !
    ! !USES:
      !$acc routine seq
    use elm_varcon      , only : denice, denh2o, wimp, pondmx_urban, frac_from_uphill
    use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall
    use column_varcon   , only : icol_road_imperv, icol_road_perv
    use elm_varpar      , only : nlevsoi, nlevgrnd, maxpatch_pft
    use elm_varpar      , only : nlayer, nlayert
    use elm_varctl      , only : use_var_soil_thick, use_IM2_hillslope_hydrology
    use pftvarcon       , only : humhol_ht
    use landunit_varcon  , only : istsoil
    use SoilWaterMovementMod, only : zengdecker_2009_with_var_soil_thick
    use timeinfoMod      , only : nstep_mod
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
    integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    real(r8), intent(in)  :: dtime
    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,fc,g,l,t,i,topo_index                             !indices
    integer  :: tpair                                                 ! paired bog hummock/hollow topounit
    integer  :: nlevbed                                    !# levels to bedrock
    real(r8) :: xs(bounds%begc:bounds%endc)                !excess soil water above urban ponding limit
    real(r8) :: vol_ice(bounds%begc:bounds%endc,1:nlevgrnd) !partial volume of ice lens in layer
    real(r8) :: fff(bounds%begc:bounds%endc)               !decay factor (m-1)
    real(r8) :: s1                                         !variable to calculate qinmax
    real(r8) :: su                                         !variable to calculate qinmax
    real(r8) :: v                                          !variable to calculate qinmax
    real(r8) :: qinmax                                     !maximum infiltration capacity (mm/s)
    real(r8) :: A(bounds%begc:bounds%endc)                 !fraction of the saturated area
    real(r8) :: ex(bounds%begc:bounds%endc)                !temporary variable (exponent)
    real(r8) :: frac_from_uphill_eff                       ! fraction of uphill storage applied this timestep
    real(r8) :: from_uphill_used(bounds%begt:bounds%endt)  ! topounit from_uphill storage consumed this timestep (mm)
    real(r8) :: humhol_ht_eff                              ! hummock/hollow height used by this topounit (m)
    real(r8) :: qflx_top_soil_at_surf(bounds%begc:bounds%endc) ! qflx_top_soil used to form qflx_surf (mm/s)
    real(r8) :: qflx_top_soil_after_snow(bounds%begc:bounds%endc) ! qflx_top_soil after snow/surface water additions (mm/s)
    real(r8) :: top_moist(bounds%begc:bounds%endc)         !temporary, soil moisture in top VIC layers
    real(r8) :: top_max_moist(bounds%begc:bounds%endc)     !temporary, maximum soil moisture in top VIC layers
    real(r8) :: top_ice(bounds%begc:bounds%endc)           !temporary, ice len in top VIC layers
    logical  :: print_surface_runoff_diag
    logical, parameter :: debug_surface_runoff_diag = .false.
    character(len=32) :: subname = 'SurfaceRunoff'         !subroutine name
    !-----------------------------------------------------------------------

    associate(                                                        &
         snl              =>    col_pp%snl                             , & ! Input:  [integer  (:)   ]  minus number of snow layers
         dz               =>    col_pp%dz                              , & ! Input:  [real(r8) (:,:) ]  layer depth (m)
         nlev2bed         =>    col_pp%nlevbed                         , & ! Input:  [integer  (:)   ]  number of layers to bedrock

         sucsat           =>    soilstate_vars%sucsat_col           , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
         watsat           =>    soilstate_vars%watsat_col           , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         wtfact           =>    soilstate_vars%wtfact_col           , & ! Input:  [real(r8) (:)   ]  maximum saturated fraction for a gridcell
         hksat            =>    soilstate_vars%hksat_col            , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
         bsw              =>    soilstate_vars%bsw_col              , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"

         h2osoi_ice       =>    col_ws%h2osoi_ice      , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)                                
         h2osoi_liq       =>    col_ws%h2osoi_liq      , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)                            
         h2osfc           =>    col_ws%h2osfc          , & !Output: [real(r8) (:)   ]  surface water (mm) 

         qflx_rain_grnd   =>    col_wf%qflx_rain_grnd   , & ! Input:  [real(r8) (:)   ]  rain on ground after interception (mm/s)
         qflx_snow_grnd   =>    col_wf%qflx_snow_grnd   , & ! Input:  [real(r8) (:)   ]  snow on ground after interception (mm/s)
         qflx_snow_h2osfc =>    col_wf%qflx_snow_h2osfc , & ! Input:  [real(r8) (:)   ]  snow falling on surface water (mm/s)
         qflx_floodc      =>    col_wf%qflx_floodc      , & ! Input:  [real(r8) (:)   ]  column flux of flood water from RTM
         qflx_evap_grnd   =>    col_wf%qflx_evap_grnd   , & ! Input:  [real(r8) (:)   ]  ground surface evaporation rate (mm H2O/s) [+]
         qflx_top_soil    =>    col_wf%qflx_top_soil    , & ! Output: [real(r8) (:)   ]  net water input into soil from top (mm/s)
         qflx_from_uphill =>    col_wf%qflx_from_uphill , & ! Output: [real(r8) (:)   ]  water received from uphill topounit (mm/s)
         qflx_surf        =>    col_wf%qflx_surf        , & ! Output: [real(r8) (:)   ]  surface runoff (mm H2O /s)
         qflx_irrig       =>    col_wf%qflx_irrig       , & ! Input:  [real(r8) (:)   ]  irrigation flux (mm H2O /s)
         irrig_rate       =>    veg_wf%irrig_rate       , & ! Input:  [real(r8) (:)   ]  current irrigation rate (applied if !n_irrig_steps_left > 0) [mm/s]

         zwt              =>    soilhydrology_vars%zwt_col          , & ! Input:  [real(r8) (:)   ]  water table depth (m)
         max_moist        =>    soilhydrology_vars%max_moist_col    , & ! Input:  [real(r8) (:,:) ]  maximum soil moisture (ice + liq, mm)
         frost_table      =>    soilhydrology_vars%frost_table_col  , & ! Input:  [real(r8) (:)   ]  frost table depth (m)
         zwt_perched      =>    soilhydrology_vars%zwt_perched_col  , & ! Input:  [real(r8) (:)   ]  perched water table depth (m)
         b_infil          =>    soilhydrology_vars%b_infil_col      , & ! Input:  [real(r8) (:)   ]  VIC b infiltration parameter
         moist            =>    soilhydrology_vars%moist_col        , & ! Input:  [real(r8) (:,:) ]  soil moisture in each VIC layers (liq, mm)
         hkdepth          =>    soilhydrology_vars%hkdepth_col      , & ! Input:  [real(r8) (:)   ]  decay factor (m)
         fover            =>    soilhydrology_vars%fover            , & ! Input:  [real(r8) (:)   ]  decay factor for saturation fraction (m)
         origflag         =>    soilhydrology_vars%origflag         , & ! Input:  logical
         fcov             =>    soilhydrology_vars%fcov_col         , & ! Output: [real(r8) (:)   ]  fractional impermeable area
         fsat             =>    soilhydrology_vars%fsat_col         , & ! Output: [real(r8) (:)   ]  fractional area with water table at surface
         fracice          =>    soilhydrology_vars%fracice_col      , & ! Output: [real(r8) (:,:) ]  fractional impermeability (-)
         icefrac          =>    soilhydrology_vars%icefrac_col      , & ! Output: [real(r8) (:,:) ]
         ice              =>    soilhydrology_vars%ice_col          , & ! Output: [real(r8) (:,:) ]  ice len in each VIC layers(ice, mm)
         max_infil        =>    soilhydrology_vars%max_infil_col    , & ! Output: [real(r8) (:)   ]  maximum infiltration capacity in VIC (mm)
         i_0              =>    soilhydrology_vars%i_0_col            & ! Output: [real(r8) (:)   ]  column average soil moisture in top VIC layers (mm)
         )

      ! Get time step

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         nlevbed = nlev2bed(c)
         qflx_top_soil_at_surf(c) = spval
         qflx_top_soil_after_snow(c) = spval
         do j = 1,nlevbed

            ! Porosity of soil, partial volume of ice and liquid, fraction of ice in each layer,
            ! fractional impermeability

            vol_ice(c,j) = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
            if (origflag == 1) then
               icefrac(c,j) = min(1._r8,h2osoi_ice(c,j)/(h2osoi_ice(c,j)+h2osoi_liq(c,j)))
            else
               icefrac(c,j) = min(1._r8,vol_ice(c,j)/watsat(c,j))
            endif

            fracice(c,j) = max(0._r8,exp(-3._r8*(1._r8-icefrac(c,j)))- exp(-3._r8))/(1.0_r8-exp(-3._r8))
         end do
      end do

      ! Saturated fraction

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         g = col_pp%gridcell(c)
         t = col_pp%topounit(c)
         topo_index = top_pp%topo_grc_ind(t)
         humhol_ht_eff = humhol_ht
         if (use_humhol .and. top_pp%is_bog(t) .and. top_pp%peat_depth(t) > 0._r8) then
            do tpair = grc_pp%topi(g), grc_pp%topf(g)
               if (tpair /= t .and. top_pp%active(tpair) .and. top_pp%is_bog(tpair) .and. &
                    top_pp%peat_depth(tpair) > 0._r8) then
                  humhol_ht_eff = abs(top_pp%elevation(tpair) - top_pp%elevation(t))
                  exit
               end if
            end do
            if (humhol_ht_eff <= 0._r8) humhol_ht_eff = humhol_ht
         end if
         fff(c) = fover(g)
         if (zengdecker_2009_with_var_soil_thick) then
            nlevbed = nlev2bed(c)
            fff(c) = 0.5_r8 * col_pp%zi(c,nlevsoi) / min(col_pp%zi(c,nlevbed), col_pp%zi(c,nlevsoi))
         end if
         if (use_vichydro) then
            top_moist(c) = 0._r8
            top_ice(c) = 0._r8
            top_max_moist(c) = 0._r8
            do j = 1, nlayer - 1
               top_ice(c) = top_ice(c) + ice(c,j)
               top_moist(c) =  top_moist(c) + moist(c,j) + ice(c,j)
               top_max_moist(c) = top_max_moist(c) + max_moist(c,j)
            end do
            if(top_moist(c)> top_max_moist(c)) top_moist(c)= top_max_moist(c)
            top_ice(c)     = max(0._r8,top_ice(c))
            max_infil(c)   = (1._r8+b_infil(c)) * top_max_moist(c)
            ex(c)          = b_infil(c) / (1._r8 + b_infil(c))
            A(c)           = 1._r8 - (1._r8 - top_moist(c) / top_max_moist(c))**ex(c)
            i_0(c)         = max_infil(c) * (1._r8 - (1._r8 - A(c))**(1._r8/b_infil(c)))
            fsat(c)        = A(c)  !for output
         else
            fsat(c) = wtfact(c) * exp(-0.5_r8*fff(c)*zwt(c))
         end if
         if (use_humhol .and. top_pp%peat_depth(t) > 0._r8) then
            fsat(c) = exp(-3.0_r8/humhol_ht_eff*zwt(c))
         end if

#if (defined MARSH)
         fsat(c) = 1.0 * exp(-3.0_r8/humhol_ht*(zwt(c)))   !at 30cm, hummock saturated at 5% changed to 0.1 TAO
!          if (topo_index .eq. 2) fsat(c) = min(1.0 * exp(-3.0_r8/humhol_ht*(zwt(c)-h2osfc(c)/1000.+humhol_ht)), 1._r8) !TAO 0.3 t0 0.1, 0.15 to 0.35 !bsulman: what does 0.15 represent?
#endif
         ! use perched water table to determine fsat (if present)
         if ( frost_table(c) > zwt(c)) then
            if (use_vichydro) then
               fsat(c) =  A(c)
            else
               fsat(c) = wtfact(c) * exp(-0.5_r8*fff(c)*zwt(c))
            end if
            if (use_humhol .and. top_pp%peat_depth(t) > 0._r8) then
               fsat(c) = exp(-3.0_r8/humhol_ht_eff*zwt(c))
            end if

#if (defined MARSH)
            fsat(c) = 1.0 * exp(-3.0_r8/humhol_ht*(zwt(c)))   !at 30cm, hummock saturated at 5%
!             if (topo_index .eq. 2) fsat(c) = min(1.0 * exp(-3.0_r8/humhol_ht*(zwt(c)-h2osfc(c)/1000.+humhol_ht)), 1._r8) !TAO 0.3 t 0.1, 0.15 to 0.35
#endif

         else
            if ( frost_table(c) > zwt_perched(c)) then
               fsat(c) = wtfact(c) * exp(-0.5_r8*fff(c)*zwt_perched(c))!*( frost_table(c) - zwt_perched(c))/4.0
            endif
            if (use_humhol .and. top_pp%peat_depth(t) > 0._r8) then
               fsat(c) = exp(-3.0_r8/humhol_ht_eff*zwt(c))
            end if

#if (defined MARSH)
            fsat(c) = 1.0 * exp(-3.0_r8/humhol_ht*(zwt(c))) !at 30cm, hummock saturated at 5%
!             if (topo_index .eq. 2) fsat(c) = min(1.0 * exp(-3.0_r8/humhol_ht*(zwt(c)-h2osfc(c)/1000.+humhol_ht)), 1._r8) !TAO 0.3 t 1.5, 0.15 to 0.35
#endif 
         endif
         if (origflag == 1) then
            if (use_vichydro) then
#ifndef _OPENACC
               call endrun(msg="VICHYDRO is not available for origflag=1"//errmsg(__FILE__, __LINE__))
#endif
            else
               fcov(c) = (1._r8 - fracice(c,1)) * fsat(c) + fracice(c,1)
            end if
         else
            fcov(c) = fsat(c)
         endif
      end do

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         l = col_pp%landunit(c)
         g = col_pp%gridcell(c)
         t = col_pp%topounit(c)
         topo_index = top_pp%topo_grc_ind(t)
         qflx_top_soil_at_surf(c) = qflx_top_soil(c)
         ! no qflx_surf in polygonal ground
         if (lun_pp%ispolygon(l)) then
            qflx_surf(c) = 0._r8

         ! assume qinmax large relative to qflx_top_soil in control
         else if (origflag == 1) then
           if (use_humhol) then
              if (lun_pp%itype(l) == istsoil) then
                ! Fen/lagg and hollow topounits should not generate fast
                ! saturation-excess runoff; keep that pathway on hummocks/uplands.
                if (grc_pp%ntopounits(g) > 1 .and. (topo_index == 1 .or. topo_index == 2)) then
                   qflx_surf(c) = 0._r8
                else
                   qflx_surf(c) = fcov(c) * qflx_top_soil(c) !TAO
                end if
              else
                qflx_surf(c) = 0._r8
              endif
! #elif (defined MARSH)
!             if (topo_index .eq. 1) then  !XS - only compute sfc runoff from hummock, send to hollow
!              qflx_surf(c) = fcov(c) * qflx_top_soil(c) !TAO
!             else
!              qflx_surf(c) = 0._r8   !turn off surface runoff for hollow
!             endif
           else
              qflx_surf(c) = fcov(c) * qflx_top_soil(c)
           endif

         else
            ! only send fast runoff directly to streams
            qflx_surf(c) = fsat(c) * qflx_top_soil(c) !TAO
         endif




      end do

      ! Determine water in excess of ponding limit for urban roof and impervious road.
      ! Excess goes to surface runoff. No surface runoff for sunwall and shadewall.

      do fc = 1, num_urbanc
         c = filter_urbanc(fc)
         if (col_pp%itype(c) == icol_roof .or. col_pp%itype(c) == icol_road_imperv) then

            ! If there are snow layers then all qflx_top_soil goes to surface runoff
            if (snl(c) < 0) then
               qflx_surf(c) = max(0._r8,qflx_top_soil(c))
            else
               xs(c) = max(0._r8, &
                    h2osoi_liq(c,1)/dtime + qflx_top_soil(c) - qflx_evap_grnd(c) - &
                    pondmx_urban/dtime)
               if (xs(c) > 0.) then
                  h2osoi_liq(c,1) = pondmx_urban
               else
                  h2osoi_liq(c,1) = max(0._r8,h2osoi_liq(c,1)+ &
                       (qflx_top_soil(c)-qflx_evap_grnd(c))*dtime)
               end if
               qflx_surf(c) = xs(c)
            end if
         else if (col_pp%itype(c) == icol_sunwall .or. col_pp%itype(c) == icol_shadewall) then
            qflx_surf(c) = 0._r8
         end if
         ! send flood water flux to runoff for all urban columns
         qflx_surf(c) = qflx_surf(c)  + qflx_floodc(c)

      end do

      ! remove stormflow and snow on h2osfc from qflx_top_soil
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         t = col_pp%topounit(c)
         ! add flood water flux to qflx_top_soil
         qflx_top_soil(c) = qflx_top_soil(c) + qflx_snow_h2osfc(c) + qflx_floodc(c)
         qflx_top_soil_after_snow(c) = qflx_top_soil(c)
      end do

      ! when using the subgrid hillslope lateral flow mechanism (IM2 from NGEE Arctic):
      ! 1. calculate fraction of topounit flow that goes to each column
      ! 2. Calculate the weighted flux from uphill (qflx_from_uphill) and add to qflx_top_soil 
      if (use_IM2_hillslope_hydrology .or. use_humhol) then
         ! calculate the sum of column weights on each topounit for columns in the hydrologyc filter
         ! This will be the istsoil, istcrop, and icol_road_perv subset of urban columns
         ! First zero the topounit sum of weights. This zeros some multiple times, but no harm done.
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            t = col_pp%topounit(c)
            top_pp%uphill_wt(t) = 0._r8
            from_uphill_used(t) = 0._r8
            if (top_ws%from_uphill(t) < 1.e-20_r8) then
               top_ws%from_uphill(t) = 0._r8
            endif
         end do
         ! Next sum the relevant column weights
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            t = col_pp%topounit(c)
            top_pp%uphill_wt(t) = top_pp%uphill_wt(t) + col_pp%wttopounit(c)
         end do

         ! flow from uphill topounit goes to top of soil for soil, crop, and pervious road columns
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            t = col_pp%topounit(c)
            if (top_pp%uphill_wt(t) > 0._r8) then
               frac_from_uphill_eff = frac_from_uphill
               if (top_pp%peat_depth(t) > 0._r8) frac_from_uphill_eff = 1.0_r8
               qflx_from_uphill(c) = (col_pp%wttopounit(c)/top_pp%uphill_wt(t)) * &
                    (frac_from_uphill_eff * top_ws%from_uphill(t)) / dtime
               qflx_top_soil(c) = qflx_top_soil(c) + qflx_from_uphill(c)
               from_uphill_used(t) = from_uphill_used(t) + qflx_from_uphill(c) * dtime
            else
               qflx_from_uphill(c) = 0._r8
            endif
         end do

         ! Consume the topounit reservoir after all receiver columns have used
         ! the same pre-timestep state.
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            t = col_pp%topounit(c)
            if (from_uphill_used(t) > 0._r8) then
               top_ws%from_uphill(t) = max(0._r8, top_ws%from_uphill(t) - from_uphill_used(t))
               if (top_ws%from_uphill(t) < 1.e-20_r8) then
                  top_ws%from_uphill(t) = 0._r8
               endif
               from_uphill_used(t) = 0._r8
            endif
         end do
      endif

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         t = col_pp%topounit(c)
         print_surface_runoff_diag = debug_surface_runoff_diag .and. &
              (top_pp%is_bog(t) .and. top_pp%topo_grc_ind(t) == 3 &
              .and. abs(top_pp%wtgcell(t) - 2.1752511754418422e-2_r8) < 1.e-10_r8 &
              .and. abs(top_pp%peat_depth(t) - 1.8641646375629948_r8) < 1.e-8_r8 &
              .and. nstep_mod >= 7670 .and. nstep_mod <= 7680)
         if (print_surface_runoff_diag) then
            write(iulog,*)'surface runoff diagnostic'
            write(iulog,*)'nstep                      = ',nstep_mod
            write(iulog,*)'column index               = ',c
            write(iulog,*)'gridcell index             = ',col_pp%gridcell(c)
            write(iulog,*)'topounit index             = ',t
            write(iulog,*)'topounit topo_grc_ind      = ',top_pp%topo_grc_ind(t)
            write(iulog,*)'topounit weight gridcell   = ',top_pp%wtgcell(t)
            write(iulog,*)'topounit peat_depth        = ',top_pp%peat_depth(t)
            write(iulog,*)'fsat                       = ',fsat(c)
            write(iulog,*)'fcov                       = ',fcov(c)
            write(iulog,*)'zwt                        = ',zwt(c)
            write(iulog,*)'zwt_perched                = ',zwt_perched(c)
            write(iulog,*)'h2osfc                     = ',h2osfc(c)
            write(iulog,*)'qflx_rain_grnd             = ',qflx_rain_grnd(c)
            write(iulog,*)'qflx_snow_grnd             = ',qflx_snow_grnd(c)
            write(iulog,*)'qflx_snow_h2osfc           = ',qflx_snow_h2osfc(c)
            write(iulog,*)'qflx_floodc                = ',qflx_floodc(c)
            write(iulog,*)'qflx_evap_grnd             = ',qflx_evap_grnd(c)
            write(iulog,*)'qflx_from_uphill           = ',qflx_from_uphill(c)
            write(iulog,*)'qflx_top_soil at surf      = ',qflx_top_soil_at_surf(c)
            write(iulog,*)'qflx_top_soil after snow   = ',qflx_top_soil_after_snow(c)
            write(iulog,*)'qflx_top_soil final        = ',qflx_top_soil(c)
            write(iulog,*)'qflx_surf                  = ',qflx_surf(c)
            write(iulog,*)'qflx_surf dt               = ',qflx_surf(c)*dtime
         endif
      end do

    end associate

   end subroutine SurfaceRunoff

   !-----------------------------------------------------------------------
   subroutine Infiltration(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc, &
        atm2lnd_vars, lnd2atm_vars, energyflux_vars, soilhydrology_vars, soilstate_vars, dtime)
     !
     ! !DESCRIPTION:
     ! Calculate infiltration into surface soil layer (minus the evaporation)
     !
     ! !USES:
      !$acc routine seq
     use shr_const_mod    , only : shr_const_pi
     use elm_varpar       , only : nlayer, nlayert
     use elm_varpar       , only : nlevsoi, nlevgrnd
     use elm_varcon       , only : denh2o, denice, roverg, wimp, mu, tfrz
     use elm_varcon       , only : pondmx, watmin
     use column_varcon    , only : icol_roof, icol_road_imperv, icol_sunwall, icol_shadewall, icol_road_perv
     use landunit_varcon  , only : istsoil, istcrop, ilowcenpoly, iflatcenpoly, ihighcenpoly
     use elm_time_manager , only : get_step_size, get_nstep
     use atm2lndType      , only : atm2lnd_type ! land river two way coupling
     use lnd2atmType      , only : lnd2atm_type
     use subgridAveMod    , only : c2g
     use abortutils       , only : endrun
     use pftvarcon        , only : humhol_ht, qflx_h2osfc_surfrate
#if defined MARSH
     use pftvarcon        , only : humhol_dist
#endif
     use elm_time_manager , only : get_step_size, get_nstep, get_curr_date, get_curr_time
#if defined MARSH
      use pftvarcon       , only : num_tide_comps, tide_baseline,tide_coeff_period, tide_coeff_phase, tide_coeff_amp,sfcflow_ratescale
      use elm_varctl      , only : tide_file
#endif
     use elm_varcon       , only : secspday
     use timeinfoMod      , only : nstep_mod
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds
     integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
     integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
     type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars         ! land river two way coupling
     type(lnd2atm_type)       , intent(in)    :: lnd2atm_vars 
     type(energyflux_type)    , intent(in)    :: energyflux_vars
     type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
     type(soilstate_type)     , intent(inout) :: soilstate_vars
     real(r8), intent(in)  :: dtime
     !
     ! !LOCAL VARIABLES:
     integer  :: c,j,l,fc,g,t,topo_index                                 ! indices
     integer  :: nlevbed                                    !# levels to bedrock
     real(r8) :: s1,su,v                                    ! variable to calculate qinmax
     real(r8) :: qinmax                                     ! maximum infiltration capacity (mm/s)
     real(r8) :: vol_ice(bounds%begc:bounds%endc,1:nlevgrnd) ! partial volume of ice lens in layer
     real(r8) :: alpha_evap(bounds%begc:bounds%endc)        ! fraction of total evap from h2osfc
     real(r8) :: qflx_evap(bounds%begc:bounds%endc)         ! local evaporation array
     real(r8) :: qflx_h2osfc_drain(bounds%begc:bounds%endc) ! bottom drainage from h2osfc
     real(r8) :: qflx_in_h2osfc(bounds%begc:bounds%endc)    ! surface input to h2osfc
     real(r8) :: qflx_in_soil(bounds%begc:bounds%endc)      ! surface input to soil
     real(r8) :: qflx_infl_excess(bounds%begc:bounds%endc)  ! infiltration excess runoff -> h2osfc
     real(r8) :: frac_infclust                              ! fraction of submerged area that is connected
     real(r8) :: fsno                                       ! copy of frac_sno
     real(r8) :: k_wet                                      ! linear reservoir coefficient for h2osfc
     real(r8) :: fac                                        ! soil wetness of surface layer
     real(r8) :: psit                                       ! negative potential of soil
     real(r8) :: hr                                         ! relative humidity
     real(r8) :: wx                                         ! partial volume of ice and water of surface layer
     real(r8) :: z_avg
     real(r8) :: rho_avg
     real(r8) :: fmelt
     real(r8) :: f_sno
     real(r8) :: imped
     real(r8) :: d
     !real(r8) :: h2osoi_vol
     real(r8) :: basis                                      ! temporary, variable soil moisture holding capacity
     real(r8) :: vdep                                       ! temporary, ice wedge polygon volumetric depression depth (m)
     real(r8) :: phi_eff                                    ! temporary, polygonal ground effective subsidence (maxes out at 0.4) (m)
     real(r8) :: swc                                        ! temporary, polygonal surface water content in m
     ! in top VIC layers for runoff calculation
     real(r8) :: rsurf_vic                                  ! temp VIC surface runoff
     real(r8) :: top_moist(bounds%begc:bounds%endc)         ! temporary, soil moisture in top VIC layers
     real(r8) :: top_max_moist(bounds%begc:bounds%endc)     ! temporary, maximum soil moisture in top VIC layers
     real(r8) :: top_ice(bounds%begc:bounds%endc)           ! temporary, ice len in top VIC layers
     real(r8) :: top_icefrac                                ! temporary, ice fraction in top VIC layers
     real(r8) :: h2osoi_left_vol1                           ! temporary, available volume in the first soil layer
     real(r8) :: pc                                         ! temporary, threhold for surface water storage to outflow
     ! variables for HUM_HOL and MARSH
     real(r8) :: dzmm(bounds%begc:bounds%endc,1:nlevsoi)   ! layer thickness (mm)
     real(r8) :: ka_ho                                     ! hydraulic conductivity terms at saturation for hummock (mmH2O/s)
     real(r8) :: ka_hu                                     ! hydraulic conductivity terms at saturation for hollow(mmH2O/s)
     real(r8) :: zwt_ho, zwt_hu                            ! water table depth for hollows and hummocks respectively (m)
     real(r8) :: s_node
     integer  :: jwt(bounds%begc:bounds%endc)
     integer  :: yr, mon, day, tod               !
     integer  :: days, seconds               !
     integer  :: ii
     integer  :: topi, topf, t_ref                        ! topounit bounds / lower-neighbor reference within a gridcell
   !   real(r8) :: h2osfc_tide
     real(r8) :: h2osfc_before
     real(r8) :: dryness_factor
     real(r8) :: ka_col(bounds%begc:bounds%endc)         ! column-scale lateral conductivity used for HUM_HOL exchange
     real(r8) :: head_depth_lat(bounds%begc:bounds%endc) ! column total-head depth below the local surface used for lateral exchange
     real(r8) :: ka_top(bounds%begt:bounds%endt)         ! natural-vegetation topounit conductivity
     real(r8) :: head_depth_top(bounds%begt:bounds%endt) ! natural-vegetation topounit total-head depth for lateral exchange
     real(r8) :: sy_top(bounds%begt:bounds%endt)         ! effective specific yield for lateral head-response limiting
     real(r8) :: qflx_lat_aqu_top(bounds%begt:bounds%endt) ! topounit lateral aquifer flux before mapping back to columns
     real(r8) :: elev_offset                             ! elevation difference between adjacent / reference topounits (m)
     real(r8) :: head_diff_pair                          ! water-surface head difference for a lateral exchange pair (m)
     real(r8) :: flux_pair                               ! lateral flux exchanged between a pair of topounits (mm/s)
     real(r8) :: pair_area_scale_ref                     ! pair-flux area scaling applied to the reference topounit
     real(r8) :: pair_area_scale_t                       ! pair-flux area scaling applied to the current topounit
     real(r8) :: lateral_cfl_denom                       ! pair head-response denominator for explicit flux limiting
     real(r8) :: max_flux_pair                           ! maximum pair flux allowed by the lateral CFL limiter (mm/s)
     real(r8) :: infil_ice_imped(1:3)                    ! frozen-soil impedance factors for top infiltration layers
     real(r8) :: h2osfc_runoff_depth                     ! surface water available for HUM/HOL runoff (mm)
     real(r8) :: regional_lateral_scale                  ! active connection fraction for regional aquifer exchange
     real(r8) :: receiver_increment                      ! local lateral aquifer inflow into the receiver topounit (mm/s)
     real(r8) :: receiver_scale                          ! area scaling from pair flux to receiver-column flux
     real(r8) :: obs_ref_zwt                             ! observed water table referenced to the first topounit in the gridcell
     real(r8), parameter :: min_full_lateral_receiver_frac = 0.10_r8 ! receiver area needed for full aquifer connection
     real(r8), parameter :: lateral_head_relax_frac = 0.25_r8 ! maximum pair-head relaxation per explicit timestep
     real(r8), parameter :: h2osfc_surf_exp_depth = 50._r8 ! depth scale for exponential HUM/HOL surface-water export (mm)
     real(r8), parameter :: humhol_frozen_infil_min_imped = 0.05_r8 ! minimum ice impedance for HUM/HOL infiltration
     real(r8), parameter :: humhol_winter_drain_min_imped = 0.02_r8 ! minimum ice impedance for HUM/HOL ZWT drainage
     integer  :: t_recv                                  ! receiver topounit index for lateral exchange
     integer  :: natveg_col_top(bounds%begt:bounds%endt) ! natural vegetation hydrology column associated with each topounit
     logical  :: bog_nonbog_pair                         ! true when a lateral exchange pair crosses bog/non-bog peat
     logical  :: ice_block_top(bounds%begt:bounds%endt)  ! true when ice suppresses HUM_HOL lateral exchange on a topounit
     logical  :: print_infiltration_diag
     logical  :: spruce_three_topounit_groundwater        ! true for the SPRUCE fen/hollow/hummock column set
     logical, parameter :: debug_infiltration_diag = .false.
   ! obs_zwt_forcing is set via namelist (elm_inparm) in elm_varctl
     !-----------------------------------------------------------------------

     associate(                                                    &
          snl                  =>    col_pp%snl                  , & ! Input:  [integer  (:)   ]  minus number of snow layers
          dz                   =>    col_pp%dz                   , & ! Input:  [real(r8) (:,:) ]  layer depth (m)
          z                    =>    col_pp%z                    , & ! Input:  [real(r8) (:,:) ]  layer depth (m)
          nlev2bed             =>    col_pp%nlevbed              , & ! Input:  [integer  (:)   ]  number of layers to bedrock
          cgridcell            =>    col_pp%gridcell             , & ! Input:  [integer  (:)   ]  column's gridcell    
          wtgcell              =>    col_pp%wtgcell              , & ! Input:  [real(r8) (:)   ]  weight (relative to gridcell)
          meangradz            =>    col_pp%meangradz            , & ! Input:  [real(r8) (:)   ]  mean topographic gradient at the column level (unitless)

          t_soisno             =>    col_es%t_soisno             , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)

          frac_h2osfc          =>    col_ws%frac_h2osfc          , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1)
          frac_h2osfc_act      =>    col_ws%frac_h2osfc_act      , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by surface water (0 to 1) without adjustment from snow fraction
          frac_sno             =>    col_ws%frac_sno_eff         , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
          h2osoi_ice           =>    col_ws%h2osoi_ice           , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)
          h2osoi_liq           =>    col_ws%h2osoi_liq           , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
          h2osno               =>    col_ws%h2osno               , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)
          snow_depth           =>    col_ws%snow_depth           , & ! Input:  [real(r8) (:)   ]  snow height (m)
          h2osfc               =>    col_ws%h2osfc               , & ! Output: [real(r8) (:)   ]  surface water (mm)
          h2orof               =>    col_ws%h2orof               , & ! Output:  [real(r8) (:)   ]  floodplain inudntion volume (mm)
          frac_h2orof          =>    col_ws%frac_h2orof          , & ! Output:  [real(r8) (:)   ]  floodplain inudntion fraction (-)
          iwp_microrel         =>    col_ws%iwp_microrel         , & ! Input:  [real(r8) (:)   ]  ice wedge polygon microtopographic relief (m)
          iwp_exclvol          =>    col_ws%iwp_exclvol          , & ! Input:  [real(r8) (:)   ]  ice wedge polygon excluded volume (m)
          iwp_ddep             =>    col_ws%iwp_ddep             , & ! Input:  [real(r8) (:)   ]  ice wedge polygon depression depth (m)
          iwp_subsidence       =>    col_ws%iwp_subsidence       , & ! Input:  [real(r8) (:)   ]  ice wedge polygon ground subsidence (m)
          salinity             =>    col_ws%salinity            , & ! Input:  [real(r8) (:,:)   ] salinity concentration (ppt)
          h2osfc_tide          =>    col_ws%h2osfc_tide         , & ! Output: [real(r8) (:) ] Tide height above surface (mm)
         
          qflx_ev_soil         =>    col_wf%qflx_ev_soil         , & ! Input:  [real(r8) (:)   ]  evaporation flux from soil (W/m**2) [+ to atm]
          qflx_evap_soi        =>    col_wf%qflx_evap_soi        , & ! Input:  [real(r8) (:)   ]  ground surface evaporation rate (mm H2O/s) [+]
          qflx_evap_grnd       =>    col_wf%qflx_evap_grnd       , & ! Input:  [real(r8) (:)   ]  ground surface evaporation rate (mm H2O/s) [+]
          qflx_top_soil        =>    col_wf%qflx_top_soil        , & ! Input:  [real(r8) (:)   ]  net water input into soil from top (mm/s)
          qflx_ev_h2osfc       =>    col_wf%qflx_ev_h2osfc       , & ! Input:  [real(r8) (:)   ]  evaporation flux from h2osfc (W/m**2) [+ to atm]
          qflx_surf            =>    col_wf%qflx_surf            , & ! Output: [real(r8) (:)   ]  surface runoff (mm H2O /s)
          qflx_h2osfc_surf     =>    col_wf%qflx_h2osfc_surf     , & ! Output: [real(r8) (:)   ]  surface water runoff (mm/s)
          qflx_infl            =>    col_wf%qflx_infl            , & ! Output: [real(r8) (:)   ] infiltration (mm H2O /s)
          qflx_gross_infl_soil =>    col_wf%qflx_gross_infl_soil , & ! Output: [real(r8) (:)] gross infiltration (mm H2O/s)
          qflx_gross_evap_soil =>    col_wf%qflx_gross_evap_soil , & ! Output: [real(r8) (:)] gross evaporation (mm H2O/s)
          qflx_h2orof_drain    =>    col_wf%qflx_h2orof_drain    , & ! Output: [real(r8) (:)] drainange from floodplain inundation volume (mm H2O/s) 
          qflx_lat_aqu         =>    col_wf%qflx_lat_aqu         , & ! Output: [real(r8) (:,:) ] total lateral flow 
          qflx_tide            =>    col_wf%qflx_tide            , & ! Output: [real(r8) (:,:) ]
          qflx_lat_aqu_layer   =>    col_wf%qflx_lat_aqu_layer   , & ! Output: [real(r8) (:,:) ] lateral flow for each layer
          qflx_surf_input      =>    col_wf%qflx_surf_input      , & ! Output: [real(r8) (:,:) ] Input to hollowInput to hollow from hummock surface runoff

          smpmin               =>    soilstate_vars%smpmin_col               , & ! Input:  [real(r8) (:)   ]  restriction for min of soil potential (mm)
          sucsat               =>    soilstate_vars%sucsat_col               , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
          watsat               =>    soilstate_vars%watsat_col               , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
          bsw                  =>    soilstate_vars%bsw_col                  , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
          hksat                =>    soilstate_vars%hksat_col                , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
          eff_porosity         =>    soilstate_vars%eff_porosity_col         , & ! Output: [real(r8) (:,:) ]  effective porosity = porosity - vol_ice
          watfc_col            =>    soilstate_vars%watfc_col                , & ! Input:  [real(r8) (:)   ]  field capacity (volumetric)
          h2osfc_thresh        =>    soilhydrology_vars%h2osfc_thresh_col    , & ! Input:  [real(r8) (:)   ]  level at which h2osfc "percolates"
          zwt                  =>    soilhydrology_vars%zwt_col              , & ! Input:  [real(r8) (:)   ]  water table depth (m)
          zwt_perched          =>    soilhydrology_vars%zwt_perched_col      , & ! Input:  [real(r8) (:)   ]  perched water table depth (m)
          fcov                 =>    soilhydrology_vars%fcov_col             , & ! Input:  [real(r8) (:)   ]  fractional area with water table at surface
          b_infil              =>    soilhydrology_vars%b_infil_col          , & ! Input:  [real(r8) (:)   ]  VIC b infiltration parameter
          frost_table          =>    soilhydrology_vars%frost_table_col      , & ! Input:  [real(r8) (:)   ]  frost table depth (m)
          fsat                 =>    soilhydrology_vars%fsat_col             , & ! Input:  [real(r8) (:)   ]  fractional area with water table at surface
          moist                =>    soilhydrology_vars%moist_col            , & ! Input:  [real(r8) (:,:) ]  soil moisture in each VIC layers (liq, mm)
          max_moist            =>    soilhydrology_vars%max_moist_col        , & ! Input:  [real(r8) (:,:) ]  maximum soil moisture (ice + liq, mm)
          max_infil            =>    soilhydrology_vars%max_infil_col        , & ! Input:  [real(r8) (:)   ]  maximum infiltration capacity in VIC (mm)
          ice                  =>    soilhydrology_vars%ice_col              , & ! Input:  [real(r8) (:,:) ]  ice len in each VIC layers(ice, mm)
          i_0                  =>    soilhydrology_vars%i_0_col              , & ! Input:  [real(r8) (:)   ]  column average soil moisture in top VIC layers (mm)
          h2osfcflag           =>    soilhydrology_vars%h2osfcflag           , & ! Input:  logical
          pc_grid              =>    soilhydrology_vars%pc                   , & ! Input:  [real(r8) (:)   ]  threshold for outflow from surface water storage
          icefrac              =>    soilhydrology_vars%icefrac_col          , & ! Output: [real(r8) (:,:) ]  fraction of ice
          h2osoi_vol           =>    col_ws%h2osoi_vol                       , & ! Input: [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
          zi                   =>    col_pp%zi                                 & ! Input: [real(r8) (:,:) ]  interface level below a "z" level (m)
              )


       ! Infiltration into surface soil layer (minus the evaporation)
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          nlevbed = nlev2bed(c)
          do j = 1,nlevbed
             ! Porosity of soil, partial volume of ice and liquid
             vol_ice(c,j) = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
             eff_porosity(c,j) = max(0.01_r8, watsat(c,j)-vol_ice(c,j))
             icefrac(c,j) = min(1._r8, vol_ice(c,j)/watsat(c,j))
          end do
       end do

       do j = 1,nlevbed
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             dzmm(c,j) = dz(c,j)*1.e3_r8
          end do
       end do
       ka_hu = 0._r8
       ka_ho = 0._r8

       !compute jwt index
       ! The layer index of the first unsaturated layer, i.e., the layer right
       ! above
       ! the water table

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          jwt(c) = nlevbed
          ! allow jwt to equal zero when zwt is in top layer
          do j = 1,nlevbed
             if(zwt(c) <= zi(c,j)) then
                jwt(c) = j-1
                exit
             end if
          enddo
       enddo

       qflx_surf_input(:) = 0._r8


       do fc = 1, num_hydrologyc
          c  = filter_hydrologyc(fc)
          g  = cgridcell(c)
          t = col_pp%topounit(c)
          topo_index = top_pp%topo_grc_ind(t)
          pc = pc_grid(g)
          
          ! partition moisture fluxes between soil and h2osfc
          if (lun_pp%itype(col_pp%landunit(c)) == istsoil .or. lun_pp%itype(col_pp%landunit(c))==istcrop) then

             ! explicitly use frac_sno=0 if snl=0
             if (snl(c) >= 0) then
                fsno=0._r8
                ! if no snow layers, sublimation is removed from h2osoi_ice in drainage
                qflx_evap(c)=qflx_evap_grnd(c)
             else
                fsno=frac_sno(c)
                qflx_evap(c)=qflx_ev_soil(c)
             endif

             !0. partition grid-level floodplain inundation volume and fraction to each column
             if (use_lnd_rof_two_way) then
                if (mod(get_nstep(),lnd_rof_coupling_nstep) == 1 .or. get_nstep() <= 1 .or. lnd_rof_coupling_nstep == 1) then
                   h2orof(c)      = atm2lnd_vars%h2orof_grc(g) * wtgcell(c)
                   frac_h2orof(c) = atm2lnd_vars%frac_h2orof_grc(g) * wtgcell(c)
                endif
                ! TODO: add inundfrac from ocean 
                if ( frac_h2orof(c) > 1.0_r8 - fsno - frac_h2osfc(c) ) then
                  frac_h2orof(c) = 1.0_r8 - fsno - frac_h2osfc(c)
                endif
             endif

             !1. partition surface inputs between soil and h2osfc
             qflx_in_soil(c) = (1._r8 - frac_h2osfc(c)) * (qflx_top_soil(c)  - qflx_surf(c))
             qflx_in_h2osfc(c) = frac_h2osfc(c) * (qflx_top_soil(c)  - qflx_surf(c))
             qflx_gross_infl_soil(c) = (1._r8 - frac_h2osfc(c)) * (qflx_top_soil(c)  - qflx_surf(c))
             !2. remove evaporation (snow treated in SnowHydrology)
             qflx_in_soil(c) = qflx_in_soil(c) - (1.0_r8 - fsno - frac_h2osfc(c))*qflx_evap(c)
             qflx_in_h2osfc(c) = qflx_in_h2osfc(c) - frac_h2osfc(c) * qflx_ev_h2osfc(c)

             if (qflx_evap(c)>0._r8) then
                qflx_gross_evap_soil(c) = (1.0_r8 - fsno - frac_h2osfc(c))*qflx_evap(c)
             else
                qflx_gross_evap_soil(c) = 0._r8
                qflx_gross_infl_soil(c) = qflx_gross_infl_soil(c)-(1.0_r8 - fsno - frac_h2osfc(c))*qflx_evap(c)
             endif

             !3. determine maximum infiltration rate
             if (use_vichydro) then
                top_moist(c)= 0._r8
                top_ice(c)=0._r8
                top_max_moist(c)= 0._r8
                do j = 1, nlayer - 1
                   top_ice(c) = top_ice(c) + ice(c,j)
                   top_moist(c) =  top_moist(c) + moist(c,j) + ice(c,j)
                   top_max_moist(c) = top_max_moist(c) + max_moist(c,j)
                end do
                top_icefrac = min(1._r8,top_ice(c)/top_max_moist(c))
                if(qflx_in_soil(c) <= 0._r8) then
                   rsurf_vic = 0._r8
                else if(max_infil(c) <= 0._r8) then
                   rsurf_vic = qflx_in_soil(c)
                else if((i_0(c) + qflx_in_soil(c)*dtime) > max_infil(c)) then             !(Eq.(3a) Wood et al. 1992)
                   rsurf_vic = (qflx_in_soil(c)*dtime - top_max_moist(c) + top_moist(c))/dtime
                else                                                                      !(Eq.(3b) Wood et al. 1992)
                   basis = 1._r8 - (i_0(c) + qflx_in_soil(c)*dtime)/max_infil(c)
                   rsurf_vic = (qflx_in_soil(c)*dtime - top_max_moist(c) + top_moist(c)    &
                        + top_max_moist(c) * basis**(1._r8 + b_infil(c)))/dtime
                end if
                rsurf_vic = min(qflx_in_soil(c), rsurf_vic)
                qinmax = (1._r8 - fsat(c)) * 10._r8**(-e_ice*top_icefrac)*(qflx_in_soil(c) - rsurf_vic)
             else
                infil_ice_imped(1:3) = 10._r8**(-e_ice*(icefrac(c,1:3)))
                if (use_humhol) then
                   infil_ice_imped(1:3) = max(infil_ice_imped(1:3), humhol_frozen_infil_min_imped)
                endif
                if ( use_modified_infil ) then
                  qinmax=minval(infil_ice_imped(1:3)*hksat(c,1:3))
                else
                  dryness_factor = 1.0_r8 !max(min((h2osoi_vol(c,1)/watfc_col(c,1))**4, 1.0_r8), 0.01_r8)
                  qinmax=(1._r8 - fsat(c)) * dryness_factor * minval(infil_ice_imped(1:3)*hksat(c,1:3))
                end if
             end if
             
             if ( use_modified_infil ) then
                ! Assume frac_h2osfc occurs on fsat
                if ( frac_h2osfc(c) >= fsat(c) ) then
                   qflx_infl_excess(c) = max(0._r8,qflx_in_soil(c) -  (1.0_r8 - frac_h2osfc(c))*qinmax)
                else
                   qflx_infl_excess(c) = max(0._r8,qflx_in_soil(c) -  (1.0_r8 - fsat(c))*qinmax)
                end if
             else
                qflx_infl_excess(c) = max(0._r8,qflx_in_soil(c) -  (1.0_r8 - frac_h2osfc(c))*qinmax)
             end if
             
             if (use_lnd_rof_two_way) then
                qflx_infl_excess(c) = max(0._r8,qflx_in_soil(c) -  (1.0_r8 - frac_h2osfc(c) - frac_h2orof(c))*qinmax)
             else
                qflx_infl_excess(c) = max(0._r8,qflx_in_soil(c) -  (1.0_r8 - frac_h2osfc(c))*qinmax)
             endif

             !4. soil infiltration and h2osfc "run-on"
             qflx_infl(c) = qflx_in_soil(c) - qflx_infl_excess(c)
             qflx_in_h2osfc(c) =  qflx_in_h2osfc(c) + qflx_infl_excess(c)
             qflx_gross_infl_soil(c) = qflx_gross_infl_soil(c)- qflx_infl_excess(c)
             !5. surface runoff from h2osfc
             if (h2osfcflag==1) then
                ! calculate runoff from h2osfc  -------------------------------------
                if (use_modified_infil) then
                  if (frac_h2osfc_act(c) <= pc .and. frac_h2osfc(c) <= pc) then 
                     frac_infclust=0.0_r8
                  else
                      if (frac_h2osfc(c) <= pc) then
                        frac_infclust=(frac_h2osfc_act(c)-pc)**mu
                      else
                        frac_infclust=(frac_h2osfc(c)-pc)**mu
                      endif
                  endif
                else
                  if (frac_h2osfc(c) <= pc) then
                    frac_infclust=0.0_r8
                  else
                    frac_infclust=(frac_h2osfc(c)-pc)**mu
                  endif
                endif
             endif

             if (lun_pp%ispolygon(col_pp%landunit(c))) then
                vdep = (2_r8*iwp_exclvol(c) - iwp_microrel(c)) * (iwp_ddep(c)/iwp_microrel(c))**3_r8 &
                       + (2_r8*iwp_microrel(c) - 3_r8*iwp_exclvol(c)) * (iwp_ddep(c)/iwp_microrel(c))**2_r8
                phi_eff = min(iwp_subsidence(c), 0.4_r8)  !fix this variable when available to pull from alt calculations
                swc = h2osfc(c)/1000_r8 ! convert to m
                
                if (swc >= vdep) then
                   if (lun_pp%polygontype(col_pp%landunit(c)) == ilowcenpoly) then
                      k_wet = (2890_r8*phi_eff**4_r8 - 1171.1_r8*phi_eff**3_r8 + 144.94_r8*phi_eff**2_r8 + 1.682_r8*phi_eff + 2.028_r8) &
                              * (710.3_r8*meangradz(c)**2_r8 - 28.736_r8*meangradz(c) + 12.74_r8)
                   else
                      k_wet = 24.925_r8 * (710.3_r8*meangradz(c)**2_r8 - 28.736_r8*meangradz(c) + 12.74_r8)
                   endif
                   qflx_h2osfc_surf(c) = k_wet * (swc - vdep) / 86400_r8 ! coefficients estimated for mm/day; convert from mm/day -> mm/s
                   qflx_h2osfc_surf(c) = min(qflx_h2osfc_surf(c), (swc - vdep)*1000_r8/dtime)
                else
                   qflx_h2osfc_surf(c) = 0._r8
                endif
                
             else
               ! limit runoff to value of storage above S(pc)
               if (use_humhol) then
                  ! HUM_HOL surface water can leave the lowest topounit and
                  ! peat topounits, while groundwater routing remains separate.
                  if (topo_index == 1 .or. top_pp%peat_depth(t) > 0._r8) then
                     h2osfc_runoff_depth = max(0._r8, h2osfc(c))
                     if (h2osfc_runoff_depth .gt. 0._r8) then
                        qflx_h2osfc_surf(c) = qflx_h2osfc_surfrate*h2osfc_surf_exp_depth**2.0_r8 * &
                             (1._r8 - exp(-h2osfc_runoff_depth/h2osfc_surf_exp_depth))
                        qflx_h2osfc_surf(c) = min(qflx_h2osfc_surf(c), h2osfc_runoff_depth / dtime)
                     else
                        qflx_h2osfc_surf(c) = 0._r8
                     endif
                  else
                     qflx_h2osfc_surf(c) = 0._r8
                  endif
#if defined MARSH
               else
                  qflx_tide(c) = 0._r8
                  if (h2osfc(c) .gt. 0._r8) then
                     qflx_h2osfc_surf(c) = min(qflx_h2osfc_surfrate*h2osfc(c)**2.0_r8,h2osfc(c) / dtime)
                  else
                     qflx_h2osfc_surf(c) = 0._r8
                  endif
#else
               else
                  if(h2osfc(c) >= h2osfc_thresh(c) .and. h2osfcflag/=0) then
                     ! spatially variable k_wet
                     k_wet=1.0_r8 * sin((rpi/180._r8) * col_pp%topo_slope(c))
                     qflx_h2osfc_surf(c) = k_wet * frac_infclust * (h2osfc(c) - h2osfc_thresh(c))
                     qflx_h2osfc_surf(c)=min(qflx_h2osfc_surf(c),(h2osfc(c) - h2osfc_thresh(c))/dtime)
                  else
                     qflx_h2osfc_surf(c)= 0._r8
                  endif
#endif
               end if
             endif

             ! cutoff lower limit
             if ( qflx_h2osfc_surf(c) < 1.0e-8) qflx_h2osfc_surf(c) = 0._r8

             ! use this for non-h2osfc code
             if(h2osfcflag==0) then
                qflx_h2osfc_surf(c)= 0._r8
                ! shift infiltration excess from h2osfc input to surface runoff
                qflx_in_h2osfc(c) =  qflx_in_h2osfc(c) - qflx_infl_excess(c)
                qflx_surf(c)= qflx_surf(c) + qflx_infl_excess(c)
                qflx_infl_excess(c) = 0._r8
             endif

             qflx_in_h2osfc(c) =  qflx_in_h2osfc(c) - qflx_h2osfc_surf(c)

             !6. update h2osfc prior to calculating bottom drainage from h2osfc
             h2osfc_before = h2osfc(c)
             h2osfc(c) = h2osfc(c) + qflx_in_h2osfc(c) * dtime

             !--  if all water evaporates, there will be no bottom drainage
             if (h2osfc(c) < 0.0) then
                qflx_infl(c) = qflx_infl(c) + h2osfc(c)/dtime
                qflx_gross_evap_soil(c) = qflx_gross_evap_soil(c) - h2osfc(c)/dtime
                h2osfc(c) = 0.0
                qflx_h2osfc_drain(c)= 0._r8
             else
                if ( use_modified_infil ) then
                   ! Assume frac_h2osfc occurs on top of fsat
                   if (frac_h2osfc(c) <= fsat(c)) then
                     qflx_h2osfc_drain(c)=0._r8
                   else
                     qflx_h2osfc_drain(c)=min((frac_h2osfc(c)-fsat(c))*qinmax,h2osfc(c)/dtime)
                   endif
                else
                   ! Original scheme
                   qflx_h2osfc_drain(c)=min(frac_h2osfc(c)*qinmax,h2osfc(c)/dtime)
                end if
             endif

             if(h2osfcflag==0) then
                qflx_h2osfc_drain(c)= max(0._r8,h2osfc(c)/dtime) !ensure no h2osfc
             endif

             if (use_humhol) then
                ! Compute column-local states first, then aggregate to topounits.
                ka_col(c) = 0._r8
                ! Use the water-surface elevation as lateral head. H2OSFC
                ! remains a separate store, but ponded water contributes to
                ! groundwater exchange head for bog and non-bog topounits.
                head_depth_lat(c) = zwt(c) - h2osfc(c)/1000._r8
                if (jwt(c) .lt. nlevbed) then
                   do j=nlevbed,jwt(c)+1,-1
                     s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
                     s_node = min(1.0_r8, s_node)
                     s1 = 0.5_r8*(1.0+s_node)
                     s1 = min(1._r8, s1)
                     ka_col(c) = ka_col(c) + (hksat(c,j)*s1**(2._r8*bsw(c,j)+3._r8))* &
                                  dzmm(c,j)/sum(dzmm(c,jwt(c)+1:nlevbed))
                   end do
                else
                     s_node = max(h2osoi_vol(c,jwt(c))/watsat(c,jwt(c)), 0.01_r8)
                     s_node = min(1.0_r8, s_node)
                     s1 = 0.5_r8*(1.0+s_node)
                     s1 = min(1._r8, s1)
                     ka_col(c) = ka_col(c) + (hksat(c,jwt(c))*s1**(2._r8*bsw(c,jwt(c))+3._r8))
                end if

             end if

#ifdef MARSH
               call get_curr_time(days, seconds)
               h2osfc_tide(c) = 0.0_r8
               if(tide_file .ne. ' ') then
#ifdef CPL_BYPASS
                  ! If external forcing tide file is specified then use that via coupler bypass
                  ! Indexing assumes that tide forcing is a time series of hourly values
                  ! write(iulog,*) h2osfc_tide
                  h2osfc_tide(c) = (atm2lnd_vars%tide_height(1,1+mod(int((days*secspday+seconds)/3600),atm2lnd_vars%tide_forcing_len)))*1000 !convert m to mm
                  col_ws%salinity(c) = atm2lnd_vars%tide_salinity(1,1+mod(int((days*secspday+seconds)/3600),atm2lnd_vars%tide_forcing_len))
                  salinity(c) = col_ws%salinity(c)
                  ! write(iulog,*) h2osfc_tide
#endif
               else
                  do ii=1,num_tide_comps
                     h2osfc_tide(c) =    h2osfc_tide(c)    +  tide_coeff_amp(ii) * sin(2.0_r8*SHR_CONST_PI*(1/tide_coeff_period(ii)*(days*secspday+seconds) + tide_coeff_phase(ii)))
                  enddo
               endif

                h2osfc_tide(c) = h2osfc_tide(c) + tide_baseline - humhol_ht*1000._r8

               !compute lateral subsurface flux
               if (jwt(c) .lt. nlevbed) then
                  do j=nlevbed,jwt(c)+1,-1
                  s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
                  s_node = min(1.0_r8, s_node)
                  s1 = 0.5_r8*(1.0+s_node)
                  s1 = min(1._r8, s1)
                  ! hu (hummock) is marsh and hollow is tidal channel
                  ka_hu = ka_hu+(hksat(c,j)*s1**(2._r8*bsw(c,j)+3._r8))* &
                           dzmm(c,j)/sum(dzmm(c,jwt(c)+1:nlevbed))
                  end do
               else
                  s_node = max(h2osoi_vol(c,jwt(c))/watsat(c,jwt(c)), 0.01_r8)
                  s_node = min(1.0_r8, s_node)
                  s1 = 0.5_r8*(1.0+s_node)
                  s1 = min(1._r8, s1)
                  ka_hu = ka_hu+(hksat(c,jwt(c))*s1**(2._r8*bsw(c,jwt(c))+3._r8))
               end if

               zwt_hu = zwt(c)
               zwt_hu = zwt_hu - h2osfc(c)/1000._r8

               zwt_ho = -(h2osfc_tide(c)+humhol_ht)/1000._r8 ! Tidal channel assumed to be saturated at bottom
               ka_hu = max(ka_hu, 1e-5_r8)

               !DMR 12/4/2015
               if (icefrac(c,min(jwt(1)+1,nlevbed)) .ge. 0.90_r8) then
                  !turn off lateral transport if any ice is present at or below,
                  !changed from 0.01 to 0.90 TAO 6/4/2021
                  !water table
                  qflx_lat_aqu(c) = 0._r8
               else
                  qflx_lat_aqu(c) =  2._r8*ka_hu * (h2osfc_tide(c)/1000._r8 - (h2osfc(c)/1000._r8 - zwt(c))) / humhol_dist
               endif
             
                ! If flooded water surface of one column is higher than the other, add faster flow since aquifer transfer (ka parameters) is slow
                if(h2osfc_tide(c)>0 .and. h2osfc_tide(c)>h2osfc(c)) then
                  ! qflx_lat_aqu(2) = qflx_lat_aqu(2) - min((h2osfc(2)-(h2osfc(1)+humhol_ht*1000.0))*sfcflow_ratescale,h2osfc(2)*0.5/dtime)
                  qflx_lat_aqu(c) = qflx_lat_aqu(c) + min((h2osfc_tide(c)-h2osfc(c))*sfcflow_ratescale,h2osfc_tide(c)*0.5/dtime)
                elseif(h2osfc(c)>0 .and. h2osfc(c) > h2osfc_tide(c)) then
                  ! qflx_lat_aqu(2) = qflx_lat_aqu(2) + min((h2osfc(1)-(h2osfc(2)-humhol_ht*1000.0))*sfcflow_ratescale,h2osfc(1)*0.5/dtime)
                  qflx_lat_aqu(c) = qflx_lat_aqu(c) - min((h2osfc(c)-h2osfc_tide(c))*sfcflow_ratescale,h2osfc(c)*0.5/dtime)
                endif
               !  write(iulog,*), 'qflx_lat_aqu(c) after', qflx_lat_aqu(c)               
               !  write(iulog,*), 'h2osfc(c) after', h2osfc(c) 
#endif


             !7. remove drainage from h2osfc and add to qflx_infl
             h2osfc(c) = h2osfc(c) - qflx_h2osfc_drain(c) * dtime
             qflx_infl(c) = qflx_infl(c) + qflx_h2osfc_drain(c)
             
             !8. add drainage from river inundation to qflx_infl (land river two way coupling)
             if (use_lnd_rof_two_way) then

               ! estimate the available volume [mm H2O] in the first soil layer for floodplain infiltration
               h2osoi_left_vol1 = max(0._r8,(pondmx+watsat(c,1)*dz(c,1)*1.e3_r8-h2osoi_ice(c,1)-watmin)) - &
                                  max(0._r8,h2osoi_liq(c,1)-watmin)
               if (h2osoi_left_vol1 < 0._r8) then
                   h2osoi_left_vol1 = 0._r8
               endif

               if (frac_h2orof(c) > 0._r8) then
                  h2osoi_left_vol1 = frac_h2orof(c) * h2osoi_left_vol1
                  qflx_h2orof_drain(c)=min(frac_h2orof(c)*qinmax,h2orof(c)/dtime)
                  ! no drainage from rof inundation if the 1st layer soil is saturated
                  qflx_h2orof_drain(c)=min(qflx_h2orof_drain(c), h2osoi_left_vol1/dtime)
               else
                  qflx_h2orof_drain(c)=0._r8
               endif

               ! remove drainage from inundation volume
               h2orof(c) = h2orof(c) - qflx_h2orof_drain(c) * dtime

               qflx_infl(c) = qflx_infl(c) + qflx_h2orof_drain(c) 
               qflx_gross_infl_soil(c) = qflx_gross_infl_soil(c) + qflx_h2osfc_drain(c) + qflx_h2orof_drain(c) 

             else
               qflx_gross_infl_soil(c) = qflx_gross_infl_soil(c) + qflx_h2osfc_drain(c)
             endif

             print_infiltration_diag = debug_infiltration_diag .and. &
                  (top_pp%is_bog(t) .and. topo_index == 3 &
                  .and. abs(top_pp%wtgcell(t) - 2.1752511754418422e-2_r8) < 1.e-10_r8 &
                  .and. abs(top_pp%peat_depth(t) - 1.8641646375629948_r8) < 1.e-8_r8 &
                  .and. nstep_mod >= 7670 .and. nstep_mod <= 7680)
             if (print_infiltration_diag) then
                write(iulog,*)'infiltration diagnostic'
                write(iulog,*)'nstep                      = ',nstep_mod
                write(iulog,*)'column index               = ',c
                write(iulog,*)'gridcell index             = ',col_pp%gridcell(c)
                write(iulog,*)'topounit index             = ',t
                write(iulog,*)'topounit topo_grc_ind      = ',top_pp%topo_grc_ind(t)
                write(iulog,*)'topounit weight gridcell   = ',top_pp%wtgcell(t)
                write(iulog,*)'topounit peat_depth        = ',top_pp%peat_depth(t)
                write(iulog,*)'fsat                       = ',fsat(c)
                write(iulog,*)'fcov                       = ',fcov(c)
                write(iulog,*)'frac_h2osfc                = ',frac_h2osfc(c)
                write(iulog,*)'frac_h2osfc_act            = ',frac_h2osfc_act(c)
                write(iulog,*)'frac_sno_eff               = ',frac_sno(c)
                write(iulog,*)'zwt                        = ',zwt(c)
                write(iulog,*)'zwt_perched                = ',zwt_perched(c)
                write(iulog,*)'h2osfc before update       = ',h2osfc_before
                write(iulog,*)'h2osfc after update        = ',h2osfc(c)
                write(iulog,*)'qflx_top_soil              = ',qflx_top_soil(c)
                write(iulog,*)'qflx_surf                  = ',qflx_surf(c)
                write(iulog,*)'qflx_in_soil               = ',qflx_in_soil(c)
                write(iulog,*)'qflx_in_h2osfc             = ',qflx_in_h2osfc(c)
                write(iulog,*)'qflx_infl                  = ',qflx_infl(c)
                write(iulog,*)'qflx_infl_excess           = ',qflx_infl_excess(c)
                write(iulog,*)'qflx_h2osfc_surf           = ',qflx_h2osfc_surf(c)
                write(iulog,*)'qflx_h2osfc_drain          = ',qflx_h2osfc_drain(c)
                write(iulog,*)'qflx_gross_infl_soil       = ',qflx_gross_infl_soil(c)
                write(iulog,*)'qflx_gross_evap_soil       = ',qflx_gross_evap_soil(c)
                write(iulog,*)'qflx_evap_grnd             = ',qflx_evap_grnd(c)
                write(iulog,*)'qflx_ev_h2osfc             = ',qflx_ev_h2osfc(c)
                write(iulog,*)'qinmax                     = ',qinmax
             endif

          else
             ! non-vegetated landunits (i.e. urban) use original CLM4 code
             if (snl(c) >= 0) then
                ! when no snow present, sublimation is removed in Drainage
                qflx_infl(c) = qflx_top_soil(c) - qflx_surf(c) - qflx_evap_grnd(c)
                qflx_gross_infl_soil(c) = qflx_top_soil(c) - qflx_surf(c)
                qflx_gross_evap_soil(c) = qflx_evap_grnd(c)
             else
                qflx_infl(c) = qflx_top_soil(c) - qflx_surf(c) &
                     - (1.0_r8 - frac_sno(c)) * qflx_ev_soil(c)
                qflx_gross_infl_soil(c) = qflx_top_soil(c) - qflx_surf(c)
                qflx_gross_evap_soil(c) = (1.0_r8 - frac_sno(c)) * qflx_ev_soil(c)
             end if
             qflx_h2osfc_surf(c) = 0._r8
          endif

       enddo

       if (use_humhol) then
          ka_top(:) = 0._r8
          head_depth_top(:) = 0._r8
          sy_top(:) = 0._r8
          qflx_lat_aqu_top(:) = 0._r8
          natveg_col_top(:) = 0
          ice_block_top(:) = .false.
          qflx_lat_aqu(:) = 0._r8

          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             t = col_pp%topounit(c)
             nlevbed = nlev2bed(c)
             l = col_pp%landunit(c)

             if (lun_pp%itype(l) /= istsoil) cycle

             if (natveg_col_top(t) /= 0 .and. natveg_col_top(t) /= c) then
                call endrun(msg='HUM_HOL expects one natural vegetation hydrology column per topounit'//errMsg(__FILE__, __LINE__))
             end if

             natveg_col_top(t) = c
             ka_top(t) = ka_col(c)
             head_depth_top(t) = head_depth_lat(c)
             if (jwt(c) == nlevbed) then
                sy_top(t) = watsat(c,nlevbed) * &
                     (1._r8 - (1._r8 + 1.e3_r8*zwt(c)/sucsat(c,nlevbed))**(-1._r8/bsw(c,nlevbed)))
             else
                j = min(jwt(c)+1,nlevbed)
                sy_top(t) = eff_porosity(c,j) * &
                     (1._r8 - (1._r8 + 1.e3_r8*zwt(c)/sucsat(c,j))**(-1._r8/bsw(c,j)))
             endif
             sy_top(t) = max(sy_top(t), 0.02_r8)
             if (icefrac(c,min(jwt(c)+1,nlevbed)) .ge. 0.01_r8) ice_block_top(t) = .true.
          end do

          do g = bounds%begg, bounds%endg
             if (grc_pp%ntopounits(g) <= 1) cycle

             topi = grc_pp%topi(g)
             topf = grc_pp%topf(g)
             spruce_three_topounit_groundwater = .false.
             if (grc_pp%ntopounits(g) == 3) then
                spruce_three_topounit_groundwater = &
                     top_pp%active(topi) .and. top_pp%active(topi+1) .and. top_pp%active(topf) .and. &
                     top_pp%topo_grc_ind(topi) == 1 .and. top_pp%topo_grc_ind(topi+1) == 2 .and. &
                     top_pp%topo_grc_ind(topf) == 3 .and. &
                     .not. top_pp%is_bog(topi) .and. top_pp%peat_depth(topi) > 0._r8 .and. &
                     top_pp%is_bog(topi+1) .and. top_pp%peat_depth(topi+1) > 0._r8 .and. &
                     top_pp%is_bog(topf) .and. top_pp%peat_depth(topf) > 0._r8
             endif

             if (obs_zwt_forcing) then
                obs_ref_zwt = atm2lnd_vars%forc_zwt_not_downscaled_grc(g)
                do t = topi, topf
                   t_ref = top_pp%regional_target_ti(t)
                   if (t_ref < topi .or. t_ref > topf .or. t_ref == t) cycle
                   if (natveg_col_top(t_ref) == 0 .or. natveg_col_top(t) == 0) cycle
                   ! Regional runs keep bog groundwater isolated from fen/upland
                   ! groundwater. In the standalone SPRUCE setup, topounit 1 is
                   ! the under-boardwalk fen in the enclosure, so include it in
                   ! ordinary ZWT exchange with the hollow and hummock.
                   bog_nonbog_pair = (top_pp%is_bog(t_ref) .and. top_pp%peat_depth(t_ref) > 0._r8) .neqv. &
                        (top_pp%is_bog(t) .and. top_pp%peat_depth(t) > 0._r8)
                   if (bog_nonbog_pair .and. .not. spruce_three_topounit_groundwater) cycle
                   ! Temporarily allow lateral exchange through partially frozen topounits.
                   ! if (ice_block_top(t_ref) .or. ice_block_top(t)) cycle

                   ka_hu = max(ka_top(t_ref), 1.e-5_r8)
                   ka_ho = max(ka_top(t), 1.e-5_r8)
                   elev_offset = top_pp%elevation(t_ref) - top_pp%elevation(topi)
                   qflx_lat_aqu_top(t_ref) = qflx_lat_aqu_top(t_ref) + &
                        ka_hu * (head_depth_top(t_ref) - (obs_ref_zwt + elev_offset))
                   elev_offset = top_pp%elevation(t) - top_pp%elevation(topi)
                   qflx_lat_aqu_top(t) = qflx_lat_aqu_top(t) + &
                        ka_ho * (head_depth_top(t) - (obs_ref_zwt + elev_offset))
                end do
             else
                do t = topi, topf
                   t_ref = top_pp%regional_target_ti(t)
                   if (t_ref < topi .or. t_ref > topf .or. t_ref == t) cycle
                   if (natveg_col_top(t_ref) == 0 .or. natveg_col_top(t) == 0) cycle
                   ! Regional runs keep bog groundwater isolated from fen/upland
                   ! groundwater. In the standalone SPRUCE setup, topounit 1 is
                   ! the under-boardwalk fen in the enclosure, so include it in
                   ! ordinary ZWT exchange with the hollow and hummock.
                   bog_nonbog_pair = (top_pp%is_bog(t_ref) .and. top_pp%peat_depth(t_ref) > 0._r8) .neqv. &
                        (top_pp%is_bog(t) .and. top_pp%peat_depth(t) > 0._r8)
                   if (bog_nonbog_pair .and. .not. spruce_three_topounit_groundwater) cycle
                   ! Temporarily allow lateral exchange through partially frozen topounits.
                   ! if (ice_block_top(t_ref) .or. ice_block_top(t)) cycle

                   ka_hu = max(ka_top(t_ref), 1.e-5_r8)
                   ka_ho = max(ka_top(t), 1.e-5_r8)
                   elev_offset = top_pp%elevation(t) - top_pp%elevation(t_ref)
                   head_diff_pair = elev_offset - (head_depth_top(t) - head_depth_top(t_ref))
                   flux_pair = 2._r8/(1._r8/ka_hu + 1._r8/ka_ho) * &
                        head_diff_pair / top_pp%lateral_dist(t)

                   if (flux_pair > 0._r8) then
                      t_recv = t_ref
                      receiver_scale = sqrt(top_pp%wtgcell(t)/top_pp%wtgcell(t_ref))
                      receiver_increment = flux_pair * receiver_scale
                   else if (flux_pair < 0._r8) then
                      t_recv = t
                      receiver_scale = sqrt(top_pp%wtgcell(t_ref)/top_pp%wtgcell(t))
                      receiver_increment = -flux_pair * receiver_scale
                   else
                      t_recv = -1
                      receiver_scale = 0._r8
                      receiver_increment = 0._r8
                   endif

                   ! Small receiver topounits should only exchange directly
                   ! with a proportional fraction of the source topounit. Apply
                   ! this to the pair flux after identifying the receiver and
                   ! before capacity checks, while preserving the conservative
                   ! source/receiver area mapping below.
                   regional_lateral_scale = 1._r8
                   if (t_recv /= -1 .and. top_pp%wtgcell(t_recv) < min_full_lateral_receiver_frac) then
                      regional_lateral_scale = max(0._r8, min(1._r8, &
                           top_pp%wtgcell(t_recv) / min_full_lateral_receiver_frac))
                      flux_pair = flux_pair * regional_lateral_scale
                      receiver_increment = receiver_increment * regional_lateral_scale
                   endif

                   pair_area_scale_ref = sqrt(top_pp%wtgcell(t)/top_pp%wtgcell(t_ref))
                   pair_area_scale_t = sqrt(top_pp%wtgcell(t_ref)/top_pp%wtgcell(t))
                   if (abs(head_diff_pair) > 0._r8) then
                      lateral_cfl_denom = pair_area_scale_t/max(sy_top(t), 0.02_r8) + &
                           pair_area_scale_ref/max(sy_top(t_ref), 0.02_r8)
                      if (lateral_cfl_denom > 0._r8) then
                         max_flux_pair = lateral_head_relax_frac * abs(head_diff_pair) * &
                              1000._r8 / (dtime * lateral_cfl_denom)
                         flux_pair = sign(min(abs(flux_pair), max_flux_pair), flux_pair)
                         if (t_recv /= -1) receiver_increment = abs(flux_pair) * receiver_scale
                      endif
                   endif

                   qflx_lat_aqu_top(t_ref) = qflx_lat_aqu_top(t_ref) + flux_pair * pair_area_scale_ref
                   qflx_lat_aqu_top(t) = qflx_lat_aqu_top(t) - flux_pair * pair_area_scale_t
                end do
             end if
          end do

          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             t = col_pp%topounit(c)
             if (lun_pp%itype(col_pp%landunit(c)) == istsoil .and. natveg_col_top(t) == c) then
                qflx_lat_aqu(c) = qflx_lat_aqu_top(t)
             else
                qflx_lat_aqu(c) = 0._r8
             endif
          end do
       end if

       ! No infiltration for impervious urban surfaces

       do fc = 1, num_urbanc
          c = filter_urbanc(fc)
          if (col_pp%itype(c) /= icol_road_perv) then
             qflx_infl(c) = 0._r8
             if (use_lnd_rof_two_way) then
                qflx_h2orof_drain(c) = 0._r8
             endif
          end if
       end do

    end associate

   end subroutine Infiltration

   !-----------------------------------------------------------------------
   subroutine WaterTable(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc, &
        soilhydrology_vars, soilstate_vars, dtime)
     !
     ! !DESCRIPTION:
     ! Calculate watertable, considering aquifer recharge but no drainage.
     !
     ! !USES:
      !$acc routine seq
     use elm_varcon       , only : pondmx, tfrz, watmin,denice,denh2o
     use elm_varpar       , only : nlevsoi, nlevgrnd
     use elm_time_manager , only : get_step_size, get_curr_date, get_curr_time
     use column_varcon    , only : icol_roof, icol_road_imperv
     use elm_varctl       , only : use_vsfm, use_var_soil_thick
     use domainMod        , only : ldomain
     use SoilWaterMovementMod, only : zengdecker_2009_with_var_soil_thick
     use timeinfoMod      , only : nstep_mod
     use landunit_varcon  , only : istsoil
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds
     integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
     integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
     type(soilstate_type)     , intent(in)    :: soilstate_vars
     real(r8), intent(in)  :: dtime
     !
     ! !LOCAL VARIABLES:
     integer  :: c,j,fc,i,l,g,t,topo_index                            ! indices
     integer  :: topi, topf, t_fen, t_hollow, c_fen, c_hollow         ! SPRUCE surface-spill indices
     integer  :: nlevbed                                 ! # layers to bedrock
     real(r8) :: xs(bounds%begc:bounds%endc)             ! water needed to bring soil moisture to watmin (mm)
     real(r8) :: dzmm(bounds%begc:bounds%endc,1:nlevgrnd) ! layer thickness (mm)
     integer  :: jwt(bounds%begc:bounds%endc)            ! index of the soil layer right above the water table (-)
     real(r8) :: rsub_bot(bounds%begc:bounds%endc)       ! subsurface runoff - bottom drainage (mm/s)
     real(r8) :: rsub_top(bounds%begc:bounds%endc)       ! subsurface runoff - topographic control (mm/s)
     real(r8) :: fff(bounds%begc:bounds%endc)            ! decay factor (m-1)
     real(r8) :: xsi(bounds%begc:bounds%endc)            ! excess soil water above saturation at layer i (mm)
     real(r8) :: rous                                    ! aquifer yield (-)
     real(r8) :: wh                                      ! smpfz(jwt)-z(jwt) (mm)
     real(r8) :: ws                                      ! summation of pore space of layers below water table (mm)
     real(r8) :: s_node                                  ! soil wetness (-)
     real(r8) :: dzsum                                   ! summation of dzmm of layers below water table (mm)
     real(r8) :: icefracsum                              ! summation of icefrac*dzmm of layers below water table (-)
     real(r8) :: fracice_rsub(bounds%begc:bounds%endc)   ! fractional impermeability of soil layers (-)
     real(r8) :: ka                                      ! hydraulic conductivity of the aquifer (mm/s)
     real(r8) :: dza                                     ! fff*(zwt-z(jwt)) (-)
     real(r8) :: available_h2osoi_liq                    ! available soil liquid water in a layer
     real(r8) :: imped
     real(r8) :: rsub_top_tot
     real(r8) :: rsub_top_layer
     real(r8) :: qcharge_tot
     real(r8) :: qcharge_layer
     real(r8) :: theta_unsat
     real(r8) :: f_unsat
     real(r8) :: s_y
     integer  :: k,k_frz,k_perch
     real(r8) :: sat_lev
     real(r8) :: s1
     real(r8) :: s2
     real(r8) :: m
     real(r8) :: b
     real(r8) :: q_perch
     real(r8) :: q_perch_max
     real(r8) :: dflag=0._r8
     real(r8) :: qcharge_temp
     real(r8) :: qflx_lat_aqu_tot
     real(r8) :: lat_aqu_store_beg
     real(r8) :: lat_aqu_store_end
     real(r8) :: qflx_lat_aqu_requested
     real(r8) :: qflx_lat_aqu_actual
     real(r8) :: layer_capacity
     real(r8) :: aquifer_excess
     real(r8) :: aquifer_bottom_add
     real(r8) :: surface_head_diff                                    ! ponded surface-water head difference (m)
     real(r8) :: surface_area_scale                                   ! source-to-receiver topounit area ratio (-)
     real(r8) :: surface_spill_source                                 ! surface water removed from source topounit (mm)
     real(r8) :: surface_spill_receiver                               ! surface water added to receiver topounit (mm)
     real(r8), parameter :: surface_head_relax_frac = 0.25_r8         ! max explicit relaxation of ponded head difference
     real(r8), parameter :: surface_spill_tol = 1.e-12_r8             ! negligible ponded surface transfer (mm)
     integer  :: natveg_col_top(bounds%begt:bounds%endt)              ! natural vegetation hydrology column by topounit
     integer  :: days, seconds  
     logical  :: print_lat_aqu_diag
     logical  :: spruce_three_topounit_surface                        ! true for SPRUCE fen/hollow/hummock setup
     logical, parameter :: debug_lat_aqu_diag = .false.
     !-----------------------------------------------------------------------

     associate(                                                            &
          snl                =>    col_pp%snl                               , & ! Input:  [integer  (:)   ]  number of snow layers
          dz                 =>    col_pp%dz                                , & ! Input:  [real(r8) (:,:) ]  layer depth (m)
          z                  =>    col_pp%z                                 , & ! Input:  [real(r8) (:,:) ]  layer depth (m)
          zi                 =>    col_pp%zi                                , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)
         nlev2bed            =>    col_pp%nlevbed                           , & ! Input:  [integer  (:)   ]  number of layers to bedrock

          t_soisno           =>    col_es%t_soisno         , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)

          h2osfc             =>    col_ws%h2osfc            , & ! Input:  [real(r8) (:)   ]  surface water (mm)
          h2osoi_liq         =>    col_ws%h2osoi_liq        , & ! Output: [real(r8) (:,:) ]  liquid water (kg/m2)
          h2osoi_ice         =>    col_ws%h2osoi_ice        , & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)
          h2osoi_vol         =>    col_ws%h2osoi_vol        , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
          frac_h2osfc        =>    col_ws%frac_h2osfc       , & ! Input:  [real(r8) (:)   ]

          qflx_dew_grnd      =>    col_wf%qflx_dew_grnd      , & ! Input:  [real(r8) (:)   ]  ground surface dew formation (mm H2O /s) [+]
          qflx_dew_snow      =>    col_wf%qflx_dew_snow      , & ! Input:  [real(r8) (:)   ]  surface dew added to snow pack (mm H2O /s) [+]

          bsw                =>    soilstate_vars%bsw_col                , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
          hksat              =>    soilstate_vars%hksat_col              , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
          sucsat             =>    soilstate_vars%sucsat_col             , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
          watsat             =>    soilstate_vars%watsat_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
          eff_porosity       =>    soilstate_vars%eff_porosity_col       , & ! Input:  [real(r8) (:,:) ]  effective porosity = porosity - vol_ice

          qflx_irrig         =>    col_wf%qflx_irrig         , & ! Input:  [real(r8) (:)   ]  irrigation flux (mm H2O /s)
          qflx_grnd_irrig_col=>    col_wf%qflx_grnd_irrig    , & ! Output: [real(r8) (:)   ]  col real groundwater irrigation flux (mm H2O /s)                                                                                                                                                               
                        
          icefrac            =>    soilhydrology_vars%icefrac_col      , &  !Output: [real(r8) (:,:) ]      
          qflx_surf_input    =>    col_wf%qflx_surf_input    , & ! Output: [real(r8) (:,:) ] surface runoff input to hollow (mmH2O/s)
          qflx_lat_aqu       =>    col_wf%qflx_lat_aqu       , & ! Output: [real(r8) (:,:) ] total lateral flow
          qflx_lat_aqu_layer =>    col_wf%qflx_lat_aqu_layer , & ! Output: [real(r8) (:,:) ] lateral flow for each layer
#ifdef MARSH
          qflx_tide          =>    col_wf%qflx_tide          , & ! Output: [real(r8) (:,:) ]
#endif

          zwt                =>    soilhydrology_vars%zwt_col            , & ! Output: [real(r8) (:)   ]  water table depth (m)
          zwt_perched        =>    soilhydrology_vars%zwt_perched_col    , & ! Output: [real(r8) (:)   ]  perched water table depth (m)
          zwt_perched_k      =>    soilhydrology_vars%zwt_perched_k_col  , & ! Output: [real(r8) (:)   ]  perched water table crossing layer
          zwt_perched_kbarrier =>  soilhydrology_vars%zwt_perched_kbarrier_col, & ! Output: [real(r8) (:)   ]  perched water table barrier layer
          zwt_perched_sat1   =>    soilhydrology_vars%zwt_perched_sat1_col, & ! Output: [real(r8) (:)   ]  saturation ratio at crossing layer
          zwt_perched_sat2   =>    soilhydrology_vars%zwt_perched_sat2_col, & ! Output: [real(r8) (:)   ]  saturation ratio below crossing layer
          frost_table        =>    soilhydrology_vars%frost_table_col    , & ! Output: [real(r8) (:)   ]  frost table depth (m)
          wa                 =>    soilhydrology_vars%wa_col             , & ! Output: [real(r8) (:)   ]  water in the unconfined aquifer (mm)
          qcharge            =>    soilhydrology_vars%qcharge_col        , & ! Input:  [real(r8) (:)   ]  aquifer recharge rate (mm/s)
          origflag           =>    soilhydrology_vars%origflag           , & ! Input:  logical
          h2osfcflag         =>    soilhydrology_vars%h2osfcflag         , & ! Input:  logical

          qflx_sub_snow      =>    col_wf%qflx_sub_snow      , & ! Output: [real(r8) (:)   ]  sublimation rate from snow pack (mm H2O /s) [+]
          qflx_drain         =>    col_wf%qflx_drain         , & ! Output: [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)
          qflx_drain_perched =>    col_wf%qflx_drain_perched , & ! Output: [real(r8) (:)   ]  perched wt sub-surface runoff (mm H2O /s)
          qflx_qrgwl         =>    col_wf%qflx_qrgwl         , & ! Output: [real(r8) (:)   ]  excess groundwater export (mm H2O /s)
          qflx_rsub_sat      =>    col_wf%qflx_rsub_sat        & ! Output: [real(r8) (:)   ]  soil saturation excess [mm h2o/s]
          )


       ! Convert layer thicknesses from m to mm

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          nlevbed = nlev2bed(c)
          do j = 1,nlevbed
             dzmm(c,j) = dz(c,j)*1.e3_r8
          end do
       end do

       if (.not.use_vsfm) then
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             qflx_drain(c)    = 0._r8
             qflx_rsub_sat(c) = 0._r8
             qflx_drain_perched(c)  = 0._r8
          end do
       endif

       ! The layer index of the first unsaturated layer, i.e., the layer right above
       ! the water table

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          nlevbed = nlev2bed(c)
          jwt(c) = nlevbed
          ! allow jwt to equal zero when zwt is in top layer
          do j = 1,nlevbed
             if(zwt(c) <= zi(c,j)) then
                if (zengdecker_2009_with_var_soil_thick .and. zwt(c) == zi(c,nlevbed)) then
                   exit
                else
                   jwt(c) = j-1
                   exit
                end if
             end if
          enddo
       end do

       !============================== QCHARGE =========================================
       ! Water table changes due to qcharge
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          nlevbed = nlev2bed(c)

          !scs: use analytical expression for aquifer specific yield
          rous = watsat(c,nlevbed) &
               * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevbed))**(-1./bsw(c,nlevbed)))
          rous=max(rous,0.02_r8)

          !--  water table is below the soil column  --------------------------------------
          g = col_pp%gridcell(c)
          l = col_pp%landunit(c)
          t = col_pp%topounit(c)
          qcharge_temp = qcharge(c)

          wa(c)  = wa(c) - qflx_grnd_irrig_col(c) * dtime
          zwt(c) = zwt(c) + (qflx_grnd_irrig_col(c) * dtime)/1000._r8/rous

          if(jwt(c) == nlevbed) then
            if (.not. (zengdecker_2009_with_var_soil_thick)) then
              wa(c)  = wa(c) + qcharge(c)  * dtime
              zwt(c) = zwt(c) - (qcharge(c)  * dtime)/1000._r8/rous
            end if
            !-- recompute jwt for following calculations  ---------------------------------
            ! allow jwt to equal zero when zwt is in top layer
            jwt(c) = nlevbed
            do j = 1,nlevbed
               if(zwt(c) <= zi(c,j)) then
                  if (zengdecker_2009_with_var_soil_thick .and. zwt(c) == zi(c,nlevbed)) then
                     exit
                  else
                     jwt(c) = j-1
                     exit
                  end if
               end if
            enddo
          else
             !-- water table within soil layers 1-9  -------------------------------------
             ! try to raise water table to account for qcharge
             qcharge_tot = qcharge(c) * dtime
             if(qcharge_tot > 0.) then !rising water table
                do j = jwt(c)+1, 1,-1
                   !scs: use analytical expression for specific yield
                   s_y = watsat(c,j) &
                        * ( 1. -  (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                   s_y=max(s_y,0.02_r8)

                   qcharge_layer=min(qcharge_tot,(s_y*(zwt(c) - zi(c,j-1))*1.e3))
                   qcharge_layer=max(qcharge_layer,0._r8)

                   if(s_y > 0._r8) zwt(c) = zwt(c) - qcharge_layer/s_y/1000._r8

                   qcharge_tot = qcharge_tot - qcharge_layer
                   if (qcharge_tot <= 0.) exit
                enddo
             else ! deepening water table (negative qcharge)
                do j = jwt(c)+1, nlevbed
                   !scs: use analytical expression for specific yield
                   s_y = watsat(c,j) &
                        * ( 1. -  (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                   s_y=max(s_y,0.02_r8)

                   qcharge_layer=max(qcharge_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                   qcharge_layer=min(qcharge_layer,0._r8)
                   qcharge_tot = qcharge_tot - qcharge_layer
                   if (qcharge_tot >= 0.) then
                      zwt(c) = zwt(c) - qcharge_layer/s_y/1000._r8
                      exit
                   else
                      zwt(c) = zi(c,j)
                   endif

                enddo
                if (qcharge_tot > 0.) zwt(c) = zwt(c) - qcharge_tot/1000._r8/rous
             endif

             !-- recompute jwt for following calculations  ---------------------------------
             ! allow jwt to equal zero when zwt is in top layer
             jwt(c) = nlevbed
             do j = 1,nlevbed
                if(zwt(c) <= zi(c,j)) then
                   if (zengdecker_2009_with_var_soil_thick .and. zwt(c) == zi(c,nlevbed)) then
                      exit
                   else
                      jwt(c) = j-1
                      exit
                   end if
                end if
             enddo
          endif
          qcharge(c) = qcharge_temp
       enddo

       if (use_humhol) then
!Compute lateral fluxes between Hummock and Hollow (code by XS)
! Water table changes due to lateral flux between hommock and hollow in aquifer 
          do fc = 1, num_hydrologyc
              c = filter_hydrologyc(fc)
              nlevbed = nlev2bed(c)
              t = col_pp%topounit(c)
              topo_index = top_pp%topo_grc_ind(t)

              qflx_lat_aqu_layer(c,:) = 0.0

       ! use analytical expression for aquifer specific yield
                rous = watsat(c,nlevbed) &
                     * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevbed))**(-1./bsw(c,nlevbed)))
                rous=max(rous,0.02_r8)

	                qflx_rsub_sat(c) = 0._r8
	                lat_aqu_store_beg = wa(c) + h2osfc(c)
	                do j = 1, nlevgrnd
	                   lat_aqu_store_beg = lat_aqu_store_beg + h2osoi_liq(c,j)
	                enddo
	                qflx_lat_aqu_requested = qflx_lat_aqu(c)
	          !-  water table is below the soil column -----------------------------------
	                qflx_lat_aqu_tot = qflx_lat_aqu(c) * dtime
                if(jwt(c) == nlevbed) then
                 if (qflx_lat_aqu_tot <= 0._r8) then
                   qflx_lat_aqu_tot = max(qflx_lat_aqu_tot, -max(0._r8, wa(c)))
                 endif
                 wa(c)  = wa(c) + qflx_lat_aqu_tot
                 zwt(c) = zwt(c) - qflx_lat_aqu_tot/1000._r8/rous
                 aquifer_excess = max(0._r8, wa(c) - 5000._r8)
                 if (aquifer_excess > 0._r8) then
                    do j = nlevbed, 1, -1
                       layer_capacity = max(0._r8, &
                            eff_porosity(c,j)*dz(c,j)*denh2o - h2osoi_liq(c,j))
                       aquifer_bottom_add = min(aquifer_excess, layer_capacity)
                       if (aquifer_bottom_add > 0._r8) then
                          h2osoi_liq(c,j) = h2osoi_liq(c,j) + aquifer_bottom_add
                          aquifer_excess = aquifer_excess - aquifer_bottom_add
                       endif
                       if (aquifer_excess <= 0._r8) exit
                    enddo
                    if (aquifer_excess > 0._r8) then
                       if (h2osfcflag == 1) then
                          h2osfc(c) = h2osfc(c) + aquifer_excess
                       else
                          qflx_qrgwl(c) = qflx_qrgwl(c) + aquifer_excess/dtime
                       endif
                    endif
                    call rediagnose_zwt_from_liq_profile(nlevbed, zengdecker_2009_with_var_soil_thick, &
                         h2osoi_liq(c,1:), eff_porosity(c,1:), dzmm(c,:), zi(c,0:), zwt(c), jwt(c))
                 endif
                 wa(c)  = min(wa(c), 5000._r8)
             else
      !-- water table within soil layers 1-9  -------------------------------------
      ! try to raise water table to account for lateral flux qflx_lat_aqu
                qflx_lat_aqu_tot = qflx_lat_aqu(c)  * dtime
                if(qflx_lat_aqu_tot > 0.) then !rising water table
                   do j = jwt(c)+1, 1,-1
      ! use analytical expression for specific yield
                      s_y = eff_porosity(c,j) &
                          * ( 1. -  (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                      s_y=max(s_y,0.02_r8)

                      qflx_lat_aqu_layer(c,j)=min(qflx_lat_aqu_tot,(s_y*(zwt(c) - zi(c,j-1))*1.e3))
                      qflx_lat_aqu_layer(c,j)=max(qflx_lat_aqu_layer(c,j),0._r8)

                      h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) &
                        + h2osoi_ice(c,j)/(dz(c,j)*denice)
                      ! Limit lateral inflow by the liquid storage available in the unfrozen pore space.
                      qflx_lat_aqu_layer(c,j) = min(qflx_lat_aqu_layer(c,j), &
                           max(0._r8, eff_porosity(c,j)*dz(c,j)*denh2o - h2osoi_liq(c,j)))
                      if (qflx_lat_aqu_layer(c,j) > 0._r8) then
                        h2osoi_liq(c,j) = h2osoi_liq(c,j) + qflx_lat_aqu_layer(c,j)
                        qflx_lat_aqu_tot = qflx_lat_aqu_tot - qflx_lat_aqu_layer(c,j)
                      end if

                      if(s_y > 0._r8) zwt(c) = zwt(c) - qflx_lat_aqu_layer(c,j)/s_y/1000._r8
                      if (qflx_lat_aqu_tot <= 0.) then
                        exit
                      else if (j .eq. 1) then
                         h2osfc(c) = h2osfc(c) + qflx_lat_aqu_tot
                         qflx_lat_aqu_tot = 0._r8
                      end if
                    end do
                 else ! deepening water table (negative lateral flux)
                   !Remove from surface water first if available
                   if (h2osfc(c) .gt. 0) then
                     available_h2osoi_liq = min(-qflx_lat_aqu_tot, h2osfc(c))
                     h2osfc(c) = h2osfc(c) - available_h2osoi_liq
                     qflx_lat_aqu_tot = qflx_lat_aqu_tot + available_h2osoi_liq
                   end if
                   do j = jwt(c)+1, nlevbed

      ! use analytical expression for specific yield
                      s_y = eff_porosity(c,j) &
                           * ( 1. -  (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                      s_y=max(s_y,0.02_r8)
                      available_h2osoi_liq = max(h2osoi_liq(c,j) - watmin, 0._r8)
                      qflx_lat_aqu_layer(c,j)=max(qflx_lat_aqu_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                      qflx_lat_aqu_layer(c,j)=max(qflx_lat_aqu_layer(c,j),-available_h2osoi_liq)
                      qflx_lat_aqu_layer(c,j)=min(qflx_lat_aqu_layer(c,j),0._r8)
                      h2osoi_liq(c,j) = h2osoi_liq(c,j) + qflx_lat_aqu_layer(c,j)
                      qflx_lat_aqu_tot = qflx_lat_aqu_tot - qflx_lat_aqu_layer(c,j)

                      if (qflx_lat_aqu_tot >= 0.) then
                         zwt(c) = zwt(c) - qflx_lat_aqu_layer(c,j)/s_y/1000._r8
                         exit
                      else
                         zwt(c) = zi(c,j)
                      endif
                   enddo
                   if (qflx_lat_aqu_tot > 0.) zwt(c) = zwt(c) - qflx_lat_aqu_tot/1000._r8/rous
                endif
             endif
             lat_aqu_store_end = wa(c) + h2osfc(c)
             do j = 1, nlevgrnd
                lat_aqu_store_end = lat_aqu_store_end + h2osoi_liq(c,j)
             enddo
             qflx_lat_aqu_actual = (lat_aqu_store_end - lat_aqu_store_beg) / dtime
             qflx_lat_aqu(c) = qflx_lat_aqu_actual
             print_lat_aqu_diag = debug_lat_aqu_diag .and. &
                  (top_pp%is_bog(t) .and. topo_index == 3 &
                  .and. abs(top_pp%wtgcell(t) - 2.1752511754418422e-2_r8) < 1.e-10_r8 &
                  .and. abs(top_pp%peat_depth(t) - 1.8641646375629948_r8) < 1.e-8_r8 &
                  .and. nstep_mod >= 7670 .and. nstep_mod <= 7680 &
                  .and. abs(qflx_lat_aqu_requested) > 1.e-8_r8)
             if (print_lat_aqu_diag) then
                write(iulog,*)'lat_aqu diagnostic'
                write(iulog,*)'nstep                      = ',nstep_mod
                write(iulog,*)'column index               = ',c
                write(iulog,*)'gridcell index             = ',col_pp%gridcell(c)
                write(iulog,*)'topounit index             = ',t
                write(iulog,*)'topounit topo_grc_ind      = ',topo_index
                write(iulog,*)'topounit weight gridcell   = ',top_pp%wtgcell(t)
                write(iulog,*)'topounit peat_depth        = ',top_pp%peat_depth(t)
                write(iulog,*)'nlevbed                    = ',nlevbed
                write(iulog,*)'jwt before recompute       = ',jwt(c)
                write(iulog,*)'zwt after lat_aqu          = ',zwt(c)
                write(iulog,*)'wa after lat_aqu           = ',wa(c)
                write(iulog,*)'h2osfc after lat_aqu       = ',h2osfc(c)
                write(iulog,*)'qflx_lat_aqu requested     = ',qflx_lat_aqu_requested
                write(iulog,*)'qflx_lat_aqu actual        = ',qflx_lat_aqu_actual
                write(iulog,*)'qflx_lat_aqu requested dt  = ',qflx_lat_aqu_requested*dtime
                write(iulog,*)'qflx_lat_aqu actual dt     = ',qflx_lat_aqu_actual*dtime
                write(iulog,*)'qflx_lat_aqu remaining dt  = ',qflx_lat_aqu_tot
                write(iulog,*)'lat_aqu store beg          = ',lat_aqu_store_beg
                write(iulog,*)'lat_aqu store end          = ',lat_aqu_store_end
                write(iulog,*)'lat_aqu store delta        = ',lat_aqu_store_end-lat_aqu_store_beg
             endif
	      !-- recompute jwt for following calculations  ---------------------------------
      ! allow jwt to equal zero when zwt is in top layer
             jwt(c) = nlevbed
             do j = 1,nlevbed
                if(zwt(c) <= zi(c,j)) then
                   jwt(c) = j-1
                   exit
                end if
             enddo
          enddo
       end if

       if (use_humhol .and. h2osfcflag == 1) then
          natveg_col_top(:) = 0
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             l = col_pp%landunit(c)
             t = col_pp%topounit(c)
             if (lun_pp%itype(l) /= istsoil) cycle
             if (natveg_col_top(t) /= 0 .and. natveg_col_top(t) /= c) then
                call endrun(msg='HUM_HOL expects one natural vegetation hydrology column per topounit'//errMsg(__FILE__, __LINE__))
             endif
             natveg_col_top(t) = c
          enddo

          do g = bounds%begg, bounds%endg
             if (grc_pp%ntopounits(g) /= 3) cycle
             topi = grc_pp%topi(g)
             topf = grc_pp%topf(g)

             spruce_three_topounit_surface = &
                  top_pp%active(topi) .and. top_pp%active(topi+1) .and. top_pp%active(topf) .and. &
                  top_pp%topo_grc_ind(topi) == 1 .and. top_pp%topo_grc_ind(topi+1) == 2 .and. &
                  top_pp%topo_grc_ind(topf) == 3 .and. &
                  .not. top_pp%is_bog(topi) .and. top_pp%peat_depth(topi) > 0._r8 .and. &
                  top_pp%is_bog(topi+1) .and. top_pp%peat_depth(topi+1) > 0._r8 .and. &
                  top_pp%is_bog(topf) .and. top_pp%peat_depth(topf) > 0._r8
             if (.not. spruce_three_topounit_surface) cycle

             t_fen = topi
             t_hollow = topi + 1
             c_fen = natveg_col_top(t_fen)
             c_hollow = natveg_col_top(t_hollow)
             if (c_fen == 0 .or. c_hollow == 0) cycle
             if (top_pp%wtgcell(t_fen) <= 0._r8 .or. top_pp%wtgcell(t_hollow) <= 0._r8) cycle

             ! The standalone SPRUCE fen is slightly lower than the hollow.
             ! Once ponded fen water overtops the hollow surface-water level,
             ! move only the above-ground head excess directly as surface water.
             surface_head_diff = top_pp%elevation(t_fen) + h2osfc(c_fen)/1000._r8 - &
                  (top_pp%elevation(t_hollow) + h2osfc(c_hollow)/1000._r8)
             if (surface_head_diff > 0._r8 .and. h2osfc(c_fen) > 0._r8) then
                surface_area_scale = top_pp%wtgcell(t_fen) / top_pp%wtgcell(t_hollow)
                surface_spill_source = surface_head_relax_frac * surface_head_diff * 1000._r8 / &
                     (1._r8 + surface_area_scale)
                surface_spill_source = min(h2osfc(c_fen), max(0._r8, surface_spill_source))
                if (surface_spill_source > surface_spill_tol) then
                   surface_spill_receiver = surface_spill_source * surface_area_scale
                   h2osfc(c_fen) = h2osfc(c_fen) - surface_spill_source
                   h2osfc(c_hollow) = h2osfc(c_hollow) + surface_spill_receiver
                   qflx_lat_aqu(c_fen) = qflx_lat_aqu(c_fen) - surface_spill_source/dtime
                   qflx_lat_aqu(c_hollow) = qflx_lat_aqu(c_hollow) + surface_spill_receiver/dtime
                endif
             endif
          enddo
       endif

       !==  BASEFLOW ==================================================
       ! perched water table code
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          nlevbed = nlev2bed(c)
          g = col_pp%gridcell(c)
          t = col_pp%topounit(c)
          ! define frost table as first frozen layer with unfrozen layer above it
          if(t_soisno(c,1) > tfrz) then
             k_frz=nlevbed
          else
             k_frz=1
          endif

          do k=2, nlevbed
            if (t_soisno(c,k-1) > tfrz .and. t_soisno(c,k) <= tfrz) then
               k_frz=k
               exit
            endif
          enddo

          frost_table(c)=z(c,k_frz)

          zwt_perched(c)=frost_table(c)
          if (use_humhol) zwt_perched(c) = spval
          zwt_perched_k(c) = spval
          zwt_perched_kbarrier(c) = spval
          zwt_perched_sat1(c) = spval
          zwt_perched_sat2(c) = spval

          !===================  water table above frost table  =============================
          ! if water table is above frost table, do not use topmodel baseflow formulation
          if (.not. use_humhol .and. &
               zwt(c) < frost_table(c) .and. t_soisno(c,k_frz) <= tfrz &
               .and. origflag == 0) then
          else
             !===================  water table below frost table  =============================
             !--  compute possible perched water table *and* groundwater table afterwards
             ! locate perched water table from bottom up starting at frost table
             ! sat_lev is an arbitrary saturation level used to determine perched water table
             sat_lev=0.9

             k_perch=1
             do k=k_frz,1,-1
                h2osoi_vol(c,k) = h2osoi_liq(c,k)/(dz(c,k)*denh2o) &
                     + h2osoi_ice(c,k)/(dz(c,k)*denice)

                if (h2osoi_vol(c,k)/watsat(c,k) <= sat_lev) then
                   k_perch=k
                   exit
                endif
             enddo

             ! if frost_table = nlevsoi, only compute perched water table if frozen
             if (t_soisno(c,k_frz) > tfrz) k_perch=k_frz

             ! if perched water table exists
             if (.not. use_humhol .and. k_frz > k_perch) then
                ! interpolate between k_perch and k_perch+1 to find perched water table height
                s1 = (h2osoi_liq(c,k_perch)/(dz(c,k_perch)*denh2o) &
                     + h2osoi_ice(c,k_perch)/(dz(c,k_perch)*denice))/watsat(c,k_perch)
                s2 = (h2osoi_liq(c,k_perch+1)/(dz(c,k_perch+1)*denh2o) &
                     + h2osoi_ice(c,k_perch+1)/(dz(c,k_perch+1)*denice))/watsat(c,k_perch+1)
                zwt_perched_k(c) = real(k_perch, r8)
                zwt_perched_kbarrier(c) = real(k_frz, r8)
                zwt_perched_sat1(c) = s1
                zwt_perched_sat2(c) = s2

                m=(z(c,k_perch+1)-z(c,k_perch))/(s2-s1)
                b=z(c,k_perch+1)-m*s2
                zwt_perched(c)=max(0._r8,m*sat_lev+b)

             endif !k_frz > k_perch
          endif

       end do

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! Renew the ice and liquid mass due to condensation

          if (snl(c)+1 >= 1) then

             ! make consistent with how evap_grnd removed in infiltration
             if (.not.use_vsfm) then
                h2osoi_liq(c,1) = h2osoi_liq(c,1) + (1._r8 - frac_h2osfc(c))*qflx_dew_grnd(c) * dtime
                h2osoi_ice(c,1) = h2osoi_ice(c,1) + (1._r8 - frac_h2osfc(c))*qflx_dew_snow(c) * dtime
                if (qflx_sub_snow(c)*dtime > h2osoi_ice(c,1)) then
                   qflx_sub_snow(c) = h2osoi_ice(c,1)/dtime
                   h2osoi_ice(c,1) = 0._r8
                else
                   h2osoi_ice(c,1) = h2osoi_ice(c,1) - (1._r8 - frac_h2osfc(c)) * qflx_sub_snow(c) * dtime
                end if
             end if
          end if
       end do


       do fc = 1, num_urbanc
          c = filter_urbanc(fc)
          ! Renew the ice and liquid mass due to condensation for urban roof and impervious road

          if (col_pp%itype(c) == icol_roof .or. col_pp%itype(c) == icol_road_imperv) then
             if (snl(c)+1 >= 1) then
                h2osoi_liq(c,1) = h2osoi_liq(c,1) + qflx_dew_grnd(c) * dtime
                h2osoi_ice(c,1) = h2osoi_ice(c,1) + (qflx_dew_snow(c) * dtime)
                if (qflx_sub_snow(c)*dtime > h2osoi_ice(c,1)) then
                   qflx_sub_snow(c) = h2osoi_ice(c,1)/dtime
                   h2osoi_ice(c,1) = 0._r8
                else
                   h2osoi_ice(c,1) = h2osoi_ice(c,1) - (qflx_sub_snow(c) * dtime)
                end if
             end if
          end if

       end do

     end associate

   end subroutine WaterTable

   !-----------------------------------------------------------------------
   subroutine Drainage(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc,  &
        soilhydrology_vars, soilstate_vars, dtime)
     !
     ! !DESCRIPTION:
     ! Calculate subsurface drainage
     !
     ! !USES:
      !$acc routine seq
     use elm_varpar       , only : nlevsoi, nlevgrnd, nlayer, nlayert
     use elm_varcon       , only : pondmx, tfrz, watmin,rpi, secspday, nlvic, aquifer_water_baseline
     use column_varcon    , only : icol_roof, icol_road_imperv, icol_road_perv
     use elm_varctl       , only : use_vsfm, use_var_soil_thick, use_firn_percolation_and_compaction
     use SoilWaterMovementMod, only : zengdecker_2009_with_var_soil_thick
     use pftvarcon        , only : rsub_top_globalmax
     use LandunitType      , only : lun_pp
     use landunit_varcon  , only : istice_mec, istice, istsoil
     use pftvarcon        , only : humhol_ht
     use timeinfoMod      , only : nstep_mod
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds
     integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
     integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     type(soilstate_type)     , intent(in)    :: soilstate_vars
     type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
     real(r8), intent(in)  :: dtime
     !
     ! !LOCAL VARIABLES:
     !character(len=32) :: subname = 'Drainage'           ! subroutine name
     real(r8), parameter :: aquifer_water_tol = 1.e-8_r8   ! tolerance for the implicit aquifer baseline (mm)
     real(r8), parameter :: fen_rsub_top_max = 2500._r8/(365._r8*86400._r8) ! Marcell S2 total seepage routed through lagg/fen (mm/s)
     real(r8), parameter :: watmin_refill_abort_flux = 1000._r8/86400._r8 ! diagnostic threshold (mm/s)
     real(r8), parameter :: layer_water_min_abort = -100._r8  ! diagnostic threshold for layer liquid water (mm)
     real(r8), parameter :: layer_water_excess_abort = 1000._r8 ! diagnostic threshold for overfull layers (mm)
     real(r8), parameter :: h2osfc_dump_threshold = 1000._r8 ! export ponded water above this depth (mm)
     real(r8), parameter :: humhol_winter_drain_min_imped = 0.02_r8 ! minimum ice impedance for HUM/HOL ZWT drainage
     integer  :: c,j,fc,i,g,t,t2,topi,topf,t_ref      ! indices
     integer  :: nlevbed                                 ! # layers to bedrock
     real(r8) :: xs(bounds%begc:bounds%endc)             ! water needed to bring soil moisture to watmin (mm)
     real(r8) :: dzmm(bounds%begc:bounds%endc,1:nlevgrnd) ! layer thickness (mm)
     integer  :: jwt(bounds%begc:bounds%endc)            ! index of the soil layer right above the water table (-)
     real(r8) :: rsub_bot(bounds%begc:bounds%endc)       ! subsurface runoff - bottom drainage (mm/s)
     real(r8) :: rsub_top(bounds%begc:bounds%endc)       ! subsurface runoff - topographic control (mm/s)
     real(r8) :: fff(bounds%begc:bounds%endc)            ! decay factor (m-1)
     real(r8) :: xsi(bounds%begc:bounds%endc)            ! excess soil water above saturation at layer i (mm)
     real(r8) :: xsia(bounds%begc:bounds%endc)           ! available pore space at layer i (mm)
     real(r8) :: xs1(bounds%begc:bounds%endc)            ! excess soil water above saturation at layer 1 (mm)
     real(r8) :: smpfz(1:nlevsoi)                        ! matric potential of layer right above water table (mm)
     real(r8) :: wtsub                                   ! summation of hk*dzmm for layers below water table (mm**2/s)
     real(r8) :: rous                                    ! aquifer yield (-)
     real(r8) :: wh                                      ! smpfz(jwt)-z(jwt) (mm)
     real(r8) :: wh_zwt                                  ! water head at the water table depth (mm)
     real(r8) :: ws                                      ! summation of pore space of layers below water table (mm)
     real(r8) :: s_node                                  ! soil wetness (-)
     real(r8) :: dzsum                                   ! summation of dzmm of layers below water table (mm)
     real(r8) :: icefracsum                              ! summation of icefrac*dzmm of layers below water table (-)
     real(r8) :: fracice_rsub(bounds%begc:bounds%endc)   ! fractional impermeability of soil layers (-)
     real(r8) :: ka                                      ! hydraulic conductivity of the aquifer (mm/s)
     real(r8) :: dza                                     ! fff*(zwt-z(jwt)) (-)
     real(r8) :: available_h2osoi_liq                    ! available soil liquid water in a layer
     real(r8) :: bottom_liq_before_refill                 ! bottom-layer liquid before watmin refill (mm)
     real(r8) :: bottom_refill_need_initial               ! initial bottom-layer water deficit to watmin (mm)
     real(r8) :: donor_available_total                    ! water found in donor layers during watmin refill (mm)
     real(r8) :: donor_liq_total                          ! donor layer liquid before withdrawals during watmin refill (mm)
     real(r8) :: watmin_refill_flux                       ! unresolved watmin refill correction flux (mm/s)
     real(r8) :: layer_liq_min                             ! minimum layer liquid water for diagnostics (mm)
     real(r8) :: layer_excess_max                          ! maximum liquid water above capacity for diagnostics (mm)
     real(r8) :: layer_capacity                            ! layer liquid-water capacity for diagnostics (mm)
     real(r8) :: h2osfc_excess_dump                        ! ponded surface water exported above the cap (mm)
     real(r8) :: h2osfc_before_dump                        ! ponded water immediately before export cap (mm)
     real(r8) :: qflx_qrgwl_before_dump                     ! qrgwl before h2osfc cap export (mm/s)
     real(r8) :: aquifer_excess                             ! regional aquifer water above the WA cap (mm)
     real(r8) :: aquifer_bottom_add                         ! aquifer overflow retained in a soil layer (mm)
     real(r8) :: rsub_top_max
     real(r8) :: rsub_top_zwt          ! water-table depth used in the topographic drainage exponential (m)
     real(r8) :: rsub_top_ref_elev     ! reference elevation for HUM_HOL topographic drainage (m)
     real(r8) :: h2osoi_vol
     real(r8) :: imped
     real(r8) :: rsub_top_tot
     real(r8) :: rsub_top_layer
     real(r8) :: qcharge_tot
     real(r8) :: qcharge_layer
     real(r8) :: theta_unsat
     real(r8) :: f_unsat
     real(r8) :: s_y
     integer  :: k,k_frz,k_perch
     real(r8) :: sat_lev
     real(r8) :: s1
     real(r8) :: s2
     real(r8) :: m
     real(r8) :: b
     real(r8) :: q_perch
     real(r8) :: q_perch_max
     real(r8) :: vol_ice
     real(r8) :: dsmax_tmp(bounds%begc:bounds%endc)       ! temporary variable for ARNO subsurface runoff calculation
     real(r8) :: rsub_tmp                 ! temporary variable for ARNO subsurface runoff calculation
     real(r8) :: frac                     ! temporary variable for ARNO subsurface runoff calculation
     real(r8) :: rel_moist                ! relative moisture, temporary variable
     real(r8) :: wtsub_vic                ! summation of hk*dzmm for layers in the third VIC layer
     real(r8) :: peat_aquifer_excess_tot    ! peat aquifer water above the cap, returned to surface water (mm)
     real(r8) :: drainage_store_beg(bounds%begc:bounds%endc)        ! storage at start of Drainage (mm)
     real(r8) :: drainage_store_after_rsub(bounds%begc:bounds%endc) ! storage after rsub removal (mm)
     real(r8) :: drainage_store_end(bounds%begc:bounds%endc)        ! storage at end of Drainage (mm)
     real(r8) :: h2osfc_drainage_beg(bounds%begc:bounds%endc)       ! ponded water at start of Drainage (mm)
     real(r8) :: rsub_top_after_rsub(bounds%begc:bounds%endc)       ! rsub_top after direct removal (mm/s)
     real(r8) :: rsub_top_final(bounds%begc:bounds%endc)            ! rsub_top used in final qflx_drain (mm/s)
     real(r8) :: qflx_lat_aqu_beg(bounds%begc:bounds%endc)          ! qflx_lat_aqu at Drainage entry (mm/s)
     real(r8) :: qflx_lat_aqu_delta                                 ! qflx_lat_aqu created inside Drainage (mm/s)
     real(r8) :: drainage_diag_resid_no_lat                         ! Drainage residual before lateral-aquifer correction (mm)
     real(r8) :: drainage_diag_resid
     integer  :: layer_liq_min_j, layer_excess_max_j
     logical  :: print_drainage_diag
     logical, parameter :: debug_drainage_diag = .false.

     !-----------------------------------------------------------------------

     associate(                                                            &
          z                  =>    col_pp%z                                 , & ! Input:  [real(r8) (:,:) ] layer depth (m)
          zi                 =>    col_pp%zi                                , & ! Input:  [real(r8) (:,:) ] interface level below a "z" level (m)
          dz                 =>    col_pp%dz                                , & ! Input:  [real(r8) (:,:) ] layer depth (m)
          snl                =>    col_pp%snl                               , & ! Input:  [integer  (:)   ] number of snow layers
         nlev2bed            =>    col_pp%nlevbed                           , & ! Input:  [integer  (:)   ]  number of layers to bedrock

          t_soisno           =>    col_es%t_soisno         , & ! Input:  [real(r8) (:,:) ] soil temperature (Kelvin)

          h2osfc             =>    col_ws%h2osfc            , & ! Input:  [real(r8) (:)   ] surface water (mm)

          bsw                =>    soilstate_vars%bsw_col                , & ! Input:  [real(r8) (:,:) ] Clapp and Hornberger "b"
          hksat              =>    soilstate_vars%hksat_col              , & ! Input:  [real(r8) (:,:) ] hydraulic conductivity at saturation (mm H2O /s)
          sucsat             =>    soilstate_vars%sucsat_col             , & ! Input:  [real(r8) (:,:) ] minimum soil suction (mm)
          watsat             =>    soilstate_vars%watsat_col             , & ! Input:  [real(r8) (:,:) ] volumetric soil water at saturation (porosity)
          eff_porosity       =>    soilstate_vars%eff_porosity_col       , & ! Input:  [real(r8) (:,:) ] effective porosity = porosity - vol_ice
          hk_l               =>    soilstate_vars%hk_l_col               , & ! Input:  [real(r8) (:,:) ] hydraulic conductivity (mm/s)
          smp_l              =>    soilstate_vars%smp_l_col              , & ! Input:  [real(r8) (:,:) ] soil matric potential (mm)

          depth              =>    soilhydrology_vars%depth_col          , & ! Input:  [real(r8) (:,:) ] VIC soil depth
          c_param            =>    soilhydrology_vars%c_param_col        , & ! Input:  [real(r8) (:)   ] baseflow exponent (Qb)
          Dsmax              =>    soilhydrology_vars%dsmax_col          , & ! Input:  [real(r8) (:)   ] max. velocity of baseflow (mm/day)
          max_moist          =>    soilhydrology_vars%max_moist_col      , & ! Input:  [real(r8) (:,:) ] maximum soil moisture (ice + liq)
          moist              =>    soilhydrology_vars%moist_col          , & ! Input:  [real(r8) (:,:) ] soil layer moisture (mm)
          Ds                 =>    soilhydrology_vars%ds_col             , & ! Input:  [real(r8) (:)   ] fracton of Dsmax where non-linear baseflow begins
          Wsvic              =>    soilhydrology_vars%Wsvic_col          , & ! Input:  [real(r8) (:)   ] fraction of maximum soil moisutre where non-liear base flow occurs
          icefrac            =>    soilhydrology_vars%icefrac_col        , & ! Output: [real(r8) (:,:) ] fraction of ice in layer
          hkdepth            =>    soilhydrology_vars%hkdepth_col        , & ! Input:  [real(r8) (:)   ] decay factor (m)
          frost_table        =>    soilhydrology_vars%frost_table_col    , & ! Input:  [real(r8) (:)   ] frost table depth (m)
          zwt                =>    soilhydrology_vars%zwt_col            , & ! Input:  [real(r8) (:)   ] water table depth (m)
          zwt_perched        =>    soilhydrology_vars%zwt_perched_col    , & ! Input:  [real(r8) (:)   ] perched water table depth (m)
          zwt_perched_k      =>    soilhydrology_vars%zwt_perched_k_col  , & ! Output: [real(r8) (:)   ] perched water table crossing layer
          zwt_perched_kbarrier =>  soilhydrology_vars%zwt_perched_kbarrier_col, & ! Output: [real(r8) (:)   ] perched water table barrier layer
          zwt_perched_sat1   =>    soilhydrology_vars%zwt_perched_sat1_col, & ! Output: [real(r8) (:)   ] saturation ratio at crossing layer
          zwt_perched_sat2   =>    soilhydrology_vars%zwt_perched_sat2_col, & ! Output: [real(r8) (:)   ] saturation ratio below crossing layer
          wa                 =>    soilhydrology_vars%wa_col             , & ! Input:  [real(r8) (:)   ] water in the unconfined aquifer (mm)
          ice                =>    soilhydrology_vars%ice_col            , & ! Input:  [real(r8) (:,:) ] soil layer moisture (mm)
          qcharge            =>    soilhydrology_vars%qcharge_col        , & ! Input:  [real(r8) (:)   ] aquifer recharge rate (mm/s)
          origflag           =>    soilhydrology_vars%origflag           , & ! Input:  logical
          h2osfcflag         =>    soilhydrology_vars%h2osfcflag         , & ! Input:  logical

          qflx_snwcp_liq     =>    col_wf%qflx_snwcp_liq     , & ! Output: [real(r8) (:)   ] excess rainfall due to snow capping (mm H2O /s) [+]
          qflx_snwcp_ice     =>    col_wf%qflx_snwcp_ice     , & ! Output: [real(r8) (:)   ] excess snowfall due to snow capping (mm H2O /s) [+]
          qflx_ice_runoff_xs =>    col_wf%qflx_ice_runoff_xs , & ! Output: solid runoff from excess ice in soil (mm H2O /s) [+]
          !qflx_dew_grnd      =>    col_wf%qflx_dew_grnd      , & ! Output: [real(r8) (:)   ] ground surface dew formation (mm H2O /s) [+]
          !qflx_dew_snow      =>    col_wf%qflx_dew_snow      , & ! Output: [real(r8) (:)   ] surface dew added to snow pack (mm H2O /s) [+]
          !qflx_sub_snow      =>    col_wf%qflx_sub_snow      , & ! Output: [real(r8) (:)   ] sublimation rate from snow pack (mm H2O /s) [+]
          qflx_drain         =>    col_wf%qflx_drain         , & ! Output: [real(r8) (:)   ] sub-surface runoff (mm H2O /s)
          qflx_qrgwl         =>    col_wf%qflx_qrgwl         , & ! Output: [real(r8) (:)   ] qflx_surf at glaciers, wetlands, lakes (mm H2O /s)
          qflx_rsub_sat      =>    col_wf%qflx_rsub_sat      , & ! Output: [real(r8) (:)   ] soil saturation excess [mm h2o/s]
          qflx_drain_perched =>    col_wf%qflx_drain_perched , & ! Output: [real(r8) (:)   ] perched wt sub-surface runoff (mm H2O /s)
          qflx_drain_vr       =>    col_wf%qflx_drain_vr       , & ! Output: [real(r8) (:)   ] sub-surface runoff (mm H2O /time step)

          h2osoi_liq         =>    col_ws%h2osoi_liq        , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
          qflx_surf_input    =>    col_wf%qflx_surf_input      , & ! Output: [real(r8) (:,:) ] surface runoff input to hollow (mmH2O/s)
          qflx_lat_aqu       =>    col_wf%qflx_lat_aqu         , & ! Output: [real(r8) (:,:) ] total lateral flow
          qflx_lat_aqu_layer =>    col_wf%qflx_lat_aqu_layer   , & ! Output: [real(r8) (:,:) ] lateral flow for each layer
          salinity             =>    col_ws%salinity           , & ! Input:  [real(r8) (:,:)   ] salinity concentration (ppt)
          qflx_adv           =>     col_wf%qflx_adv            , & ! Input:  [real(r8) (:,:) ] vertical water flux between layers (mm/s)
          h2osoi_ice         =>    col_ws%h2osoi_ice         & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)   
          )


       ! Convert layer thicknesses from m to mm

        do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          nlevbed = nlev2bed(c)
          do j = 1,nlevbed
             dzmm(c,j) = dz(c,j)*1.e3_r8

             vol_ice = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
             icefrac(c,j) = min(1._r8,vol_ice/watsat(c,j))
          end do
       end do

       ! Initial set

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         drainage_store_beg(c) = wa(c) + h2osfc(c)
         h2osfc_drainage_beg(c) = h2osfc(c)
         do j = 1, nlev2bed(c)
            drainage_store_beg(c) = drainage_store_beg(c) + h2osoi_liq(c,j) + h2osoi_ice(c,j)
          enddo
          drainage_store_after_rsub(c) = spval
          drainage_store_end(c) = spval
          rsub_top_after_rsub(c) = spval
          rsub_top_final(c) = spval
          qflx_lat_aqu_beg(c) = qflx_lat_aqu(c)
          qflx_drain(c)    = 0._r8
          rsub_bot(c)      = 0._r8
          qflx_rsub_sat(c) = 0._r8
          rsub_top(c)      = 0._r8
          fracice_rsub(c)  = 0._r8
          qflx_qrgwl(c)    = 0._r8
          qflx_drain_vr(c,:) = 0._r8
       end do

       ! The layer index of the first unsaturated layer, i.e., the layer right above
       ! the water table

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          nlevbed = nlev2bed(c)
          jwt(c) = nlevbed
          ! allow jwt to equal zero when zwt is in top layer
          do j = 1,nlevbed
             if(zwt(c) <= zi(c,j)) then
                if (zengdecker_2009_with_var_soil_thick .and. zwt(c) == zi(c,nlevbed)) then
                   exit
                else
                   jwt(c) = j-1
                   exit
                end if
             end if
          enddo
       end do

       if (use_humhol .and. any(top_pp%peat_depth(bounds%begt:bounds%endt) > 0._r8)) then
          ! Apply the normal peat WA excess export to every peat topounit.
          ! These columns can receive implicit aquifer additions from the
          ! lower-boundary solve and layer overflow correction; without an
          ! explicit export, WA can drift upward over time.
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             t = col_pp%topounit(c)

             if (.not. (top_pp%peat_depth(t) > 0._r8)) cycle

             peat_aquifer_excess_tot = max(0._r8, wa(c) - aquifer_water_baseline)
             if (peat_aquifer_excess_tot > aquifer_water_tol) then
                wa(c) = wa(c) - peat_aquifer_excess_tot
                if (h2osfcflag == 1) then
                   h2osfc(c) = h2osfc(c) + peat_aquifer_excess_tot
                else
                   qflx_qrgwl(c) = qflx_qrgwl(c) + peat_aquifer_excess_tot/dtime
                endif
             endif
          enddo
       endif

       rous = 0.2_r8

       !==  BASEFLOW ==================================================
       ! perched water table code
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          nlevbed = nlev2bed(c)
          t = col_pp%topounit(c)
          !  specify maximum drainage rate
          q_perch_max = 1.e-5_r8 * sin(col_pp%topo_slope(c) * (rpi/180._r8))

          ! if layer containing water table is frozen, compute the following:
          !     frost table, perched water table, and drainage from perched saturated layer

          ! define frost table as first frozen layer with unfrozen layer above it
          if(t_soisno(c,1) > tfrz) then
             k_frz=nlevbed
          else
             k_frz=1
          endif

          do k=2, nlevbed
             if (t_soisno(c,k-1) > tfrz .and. t_soisno(c,k) <= tfrz) then
                k_frz=k
                exit
             endif
          enddo

          frost_table(c)=z(c,k_frz)

          zwt_perched(c)=frost_table(c)
          if (use_humhol) zwt_perched(c) = spval
          qflx_drain_perched(c) = 0._r8
          zwt_perched_k(c) = spval
          zwt_perched_kbarrier(c) = spval
          zwt_perched_sat1(c) = spval
          zwt_perched_sat2(c) = spval

          !===================  water table above frost table  =============================
          ! if water table is above frost table, do not use topmodel baseflow formulation

          if (.not. use_humhol .and. zwt(c) < frost_table(c) .and. t_soisno(c,k_frz) <= tfrz &
               .and. origflag == 0) then
             ! compute drainage from perched saturated region
             wtsub = 0._r8
             q_perch = 0._r8
             do k = jwt(c)+1, k_frz
                imped=10._r8**(-e_ice*(0.5_r8*(icefrac(c,k)+icefrac(c,min(nlevbed, k+1)))))
                q_perch = q_perch + imped*hksat(c,k)*dzmm(c,k)
                wtsub = wtsub + dzmm(c,k)
             end do
             if (wtsub > 0._r8) q_perch = q_perch/wtsub

             qflx_drain_perched(c) = q_perch_max * q_perch &
                  *(frost_table(c) - zwt(c))

             ! remove drainage from perched saturated layers
             rsub_top_tot = -  qflx_drain_perched(c) * dtime
             do k = jwt(c)+1, k_frz
                rsub_top_layer=max(rsub_top_tot,-(h2osoi_liq(c,k)-watmin))
                rsub_top_layer=min(rsub_top_layer,0._r8)
                if (use_vsfm) then
                   rsub_top_layer = 0._r8
                endif
                rsub_top_tot = rsub_top_tot - rsub_top_layer

                h2osoi_liq(c,k) = h2osoi_liq(c,k) + rsub_top_layer
                qflx_drain_vr(c,k) = qflx_drain_vr(c,k) - rsub_top_layer

                if (rsub_top_tot >= 0.) then
                   zwt(c) = zwt(c) - rsub_top_layer/eff_porosity(c,k)/1000._r8
                   exit
                else
                   zwt(c) = zi(c,k)
                endif
             enddo

             ! if rsub_top_tot is greater than available water (above frost table),
             !     then decrease qflx_drain_perched by residual amount for water balance
             qflx_drain_perched(c) = qflx_drain_perched(c) + rsub_top_tot/dtime

             !-- recompute jwt  ---------------------------------------------------------
             ! allow jwt to equal zero when zwt is in top layer
             jwt(c) = nlevbed
             do j = 1,nlevbed
                if(zwt(c) <= zi(c,j)) then
                   if (zengdecker_2009_with_var_soil_thick .and. zwt(c) == zi(c,nlevbed)) then
                      exit
                   else
                      jwt(c) = j-1
                      exit
                   end if
                end if
             enddo
          else
             !===================  water table below frost table  =============================
             !--  compute possible perched water table *and* groundwater table afterwards
             ! locate perched water table from bottom up starting at frost table
             ! sat_lev is an arbitrary saturation level used to determine perched water table
             sat_lev=0.9

             k_perch=1
             do k=k_frz,1,-1
                h2osoi_vol = h2osoi_liq(c,k)/(dz(c,k)*denh2o) &
                     + h2osoi_ice(c,k)/(dz(c,k)*denice)

                if (h2osoi_vol/watsat(c,k) <= sat_lev) then
                   k_perch=k
                   exit
                endif
             enddo

             ! if frost_table = nlevsoi, only compute perched water table if frozen
             if (t_soisno(c,k_frz) > tfrz) k_perch=k_frz

             ! if perched water table exists
             if (.not. use_humhol .and. k_frz > k_perch) then
                ! interpolate between k_perch and k_perch+1 to find perched water table height
                s1 = (h2osoi_liq(c,k_perch)/(dz(c,k_perch)*denh2o) &
                     + h2osoi_ice(c,k_perch)/(dz(c,k_perch)*denice))/watsat(c,k_perch)
                s2 = (h2osoi_liq(c,k_perch+1)/(dz(c,k_perch+1)*denh2o) &
                     + h2osoi_ice(c,k_perch+1)/(dz(c,k_perch+1)*denice))/watsat(c,k_perch+1)
                zwt_perched_k(c) = real(k_perch, r8)
                zwt_perched_kbarrier(c) = real(k_frz, r8)
                zwt_perched_sat1(c) = s1
                zwt_perched_sat2(c) = s2

                m=(z(c,k_perch+1)-z(c,k_perch))/(s2-s1)
                b=z(c,k_perch+1)-m*s2
                zwt_perched(c)=max(0._r8,m*sat_lev+b)

                ! compute drainage from perched saturated region
                wtsub = 0._r8
                q_perch = 0._r8
                do k = k_perch, k_frz
                   imped=10._r8**(-e_ice*(0.5_r8*(icefrac(c,k)+icefrac(c,min(nlevbed, k+1)))))
                   q_perch = q_perch + imped*hksat(c,k)*dzmm(c,k)
                   wtsub = wtsub + dzmm(c,k)
                end do
                if (wtsub > 0._r8) q_perch = q_perch/wtsub

                qflx_drain_perched(c) = q_perch_max * q_perch &
                     *(frost_table(c) - zwt_perched(c))

                ! no perched water table drainage if using original formulation
                if(origflag == 1) qflx_drain_perched(c) = 0._r8

                ! remove drainage from perched saturated layers
                rsub_top_tot = -  qflx_drain_perched(c) * dtime
                do k = k_perch+1, k_frz
                   rsub_top_layer=max(rsub_top_tot,-(h2osoi_liq(c,k)-watmin))
                   rsub_top_layer=min(rsub_top_layer,0._r8)
                   if (use_vsfm) rsub_top_layer = 0._r8
                   rsub_top_tot = rsub_top_tot - rsub_top_layer

                   h2osoi_liq(c,k) = h2osoi_liq(c,k) + rsub_top_layer
                   qflx_drain_vr(c,k) = qflx_drain_vr(c,k) - rsub_top_layer

                   if (rsub_top_tot >= 0.) then
                      zwt_perched(c) = zwt_perched(c) - rsub_top_layer/eff_porosity(c,k)/1000._r8
                      exit
                   else
                      if (k < k_frz) then
                         s1 = (h2osoi_liq(c,k)/(dz(c,k)*denh2o) &
                              + h2osoi_ice(c,k)/(dz(c,k)*denice))/watsat(c,k)
                         s2 = (h2osoi_liq(c,k+1)/(dz(c,k+1)*denh2o) &
                              + h2osoi_ice(c,k+1)/(dz(c,k+1)*denice))/watsat(c,k+1)
                         if (s1 > sat_lev) then
                            zwt_perched(c) = 0._r8
                         else if (s2 < sat_lev) then
                            zwt_perched(c) = frost_table(c)
                         else if (abs(s2-s1) > 1.e-12_r8) then
                            m=(z(c,k+1)-z(c,k))/(s2-s1)
                            b=z(c,k+1)-m*s2
                            zwt_perched(c)=min(frost_table(c),max(0._r8,m*sat_lev+b))
                         else
                            zwt_perched(c) = min(frost_table(c), zi(c,k))
                         endif
                      else
                         zwt_perched(c) = frost_table(c)
                      endif
                   endif

                enddo

                ! if rsub_top_tot is greater than available water (above frost table),
                !     then decrease qflx_drain_perched by residual amount for water balance
                qflx_drain_perched(c) = qflx_drain_perched(c) + rsub_top_tot/dtime

             else
                qflx_drain_perched(c) = 0._r8
             endif !k_frz > k_perch

             !-- Topographic runoff  ----------------------------------------------------------------------
             ! Apply the standard regional subsurface drainage calculation to
             ! all topounits. Frozen-soil effects enter through imped below.
             fff(c)         = 1._r8/ hkdepth(c)
             dzsum = 0._r8
             icefracsum = 0._r8
             do j = max(jwt(c),1), nlevbed
                dzsum  = dzsum + dzmm(c,j)
                icefracsum = icefracsum + icefrac(c,j) * dzmm(c,j)
             end do
             ! add ice impedance factor to baseflow
             if(origflag == 1) then
                if (use_vichydro) then
#ifndef _OPENACC
                   call endrun(msg="VICHYDRO is not available for origflag=1"//errmsg(__FILE__, __LINE__))
#endif
                else
                   fracice_rsub(c) = max(0._r8,exp(-3._r8*(1._r8-(icefracsum/dzsum))) &
                        - exp(-3._r8))/(1.0_r8-exp(-3._r8))
                   imped=(1._r8 - fracice_rsub(c))
                   rsub_top_max = 5.5e-3_r8
#if (defined MARSH)                   
                   rsub_top_max = rsub_top_globalmax
#endif
                end if
             else
                if (use_vichydro) then
                   imped=10._r8**(-e_ice*min(1.0_r8,ice(c,nlayer)/max_moist(c,nlayer)))
                   dsmax_tmp(c) = Dsmax(c) * dtime/ secspday !mm/day->mm/dtime
                   rsub_top_max = dsmax_tmp(c)
                else
                   imped=10._r8**(-e_ice*(icefracsum/dzsum))
                   rsub_top_max = min(10._r8 * sin((rpi/180.) * col_pp%topo_slope(c)), rsub_top_globalmax)
                end if
             endif
             if (use_humhol .and. top_pp%is_bog(t) .and. top_pp%peat_depth(t) > 0._r8) then
                imped = max(imped, humhol_winter_drain_min_imped)
             endif

             if (use_humhol .and. top_pp%peat_depth(t) > 0._r8 .and. .not. top_pp%is_bog(t)) then
                ! Non-bog peat units represent fen/lagg seepage rather than
                ! a freely draining upland hillslope. Use the Marcell S2
                ! watershed seepage estimate concentrated through the
                ! lagg/fen area (250 cm/yr local for 4% area) as the
                ! TOPMODEL maximum before applying water-table depth.
                rsub_top_max = fen_rsub_top_max
             endif

             rsub_top_zwt = zwt(c)
             if (use_humhol .and. grc_pp%ntopounits(col_pp%gridcell(c)) > 1) then
                g = col_pp%gridcell(c)
                topi = grc_pp%topi(g)
                topf = grc_pp%topf(g)
                rsub_top_ref_elev = top_pp%elevation(t)
                do t2 = topi, topf
                   if (top_pp%active(t2)) then
                      rsub_top_ref_elev = max(rsub_top_ref_elev, top_pp%elevation(t2))
                   endif
                enddo
                ! Standalone SPRUCE has no explicit upland topounit; assume
                ! the missing upland reference is 3 m above the represented set.
                if (grc_pp%ntopounits(g) == 3) rsub_top_ref_elev = rsub_top_ref_elev + 3._r8
                rsub_top_zwt = max(0._r8, zwt(c) + rsub_top_ref_elev - top_pp%elevation(t))
             endif

             if (use_vichydro) then
                ! ARNO model for the bottom soil layer (based on bottom soil layer
                ! moisture from previous time step
                ! use watmin instead for resid_moist to be consistent with default hydrology
                rel_moist = (moist(c,nlayer) - watmin)/(max_moist(c,nlayer)-watmin)
                frac = (Ds(c) * rsub_top_max )/Wsvic(c)
                rsub_tmp = (frac * rel_moist)/dtime
                if(rel_moist > Wsvic(c))then
                   frac = (rel_moist - Wsvic(c))/(1.0_r8 - Wsvic(c))
                   rsub_tmp = rsub_tmp + (rsub_top_max * (1.0_r8 - Ds(c)/Wsvic(c)) *frac**c_param(c))/dtime
                end if
                rsub_top(c) = imped * rsub_tmp
                ! make sure baseflow isn't negative
                rsub_top(c) = max(0._r8, rsub_top(c))
             else
                if (jwt(c) == nlevbed .and. zengdecker_2009_with_var_soil_thick) then
                   rsub_top(c)    = 0._r8
                else
                   rsub_top(c)    = imped * rsub_top_max* exp(-fff(c)*rsub_top_zwt)
                end if
             end if
             if (use_vsfm) rsub_top(c) = 0._r8

             ! use analytical expression for aquifer specific yield
             rous = watsat(c,nlevbed) &
                  * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevbed))**(-1./bsw(c,nlevbed)))
             rous=max(rous,0.02_r8)

             !--  water table is below the soil column  --------------------------------------
             if(jwt(c) == nlevbed) then
               if (zengdecker_2009_with_var_soil_thick) then
                 if (-1._r8 * smp_l(c,nlevbed) < 0.5_r8 * dzmm(c,nlevbed)) then
                     zwt(c) = z(c,nlevbed) - (smp_l(c,nlevbed) / 1000._r8)
                 end if
                 rsub_top_zwt = zwt(c)
                 if (use_humhol .and. grc_pp%ntopounits(col_pp%gridcell(c)) > 1) then
                    g = col_pp%gridcell(c)
                    topi = grc_pp%topi(g)
                    topf = grc_pp%topf(g)
                    rsub_top_ref_elev = top_pp%elevation(t)
                    do t2 = topi, topf
                       if (top_pp%active(t2)) then
                          rsub_top_ref_elev = max(rsub_top_ref_elev, top_pp%elevation(t2))
                       endif
                    enddo
                    ! Standalone SPRUCE has no explicit upland topounit; assume
                    ! the missing upland reference is 3 m above the represented set.
                    if (grc_pp%ntopounits(g) == 3) rsub_top_ref_elev = rsub_top_ref_elev + 3._r8
                    rsub_top_zwt = max(0._r8, zwt(c) + rsub_top_ref_elev - top_pp%elevation(t))
                 endif
                 rsub_top(c) = imped * rsub_top_max * exp(-fff(c) * rsub_top_zwt)
                 rsub_top_tot = - rsub_top(c) * dtime
                 s_y = watsat(c,nlevbed) &
                     * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,nlevbed))**(-1./bsw(c,nlevbed)))
                 s_y=max(s_y,0.02_r8)
                 rsub_top_layer=max(rsub_top_tot,-(s_y*(zi(c,nlevbed) - zwt(c))*1.e3_r8))
                 rsub_top_layer=min(rsub_top_layer,0._r8)
                 h2osoi_liq(c,nlevbed) = h2osoi_liq(c,nlevbed) + rsub_top_layer
                 qflx_drain_vr(c,nlevbed) = qflx_drain_vr(c,nlevbed) - rsub_top_layer
                 rsub_top_tot = rsub_top_tot - rsub_top_layer
                 if (rsub_top_tot >= 0.) then
                    zwt(c) = zwt(c) - rsub_top_layer/s_y/1000._r8
                 else
                    zwt(c) = zi(c,nlevbed)
                 end if
	         if (rsub_top_tot < 0.) then
	            rsub_top(c) = rsub_top(c) + rsub_top_tot / dtime
                    rsub_top_tot = 0.
                 end if
                else
                   wa(c)  = wa(c) - rsub_top(c) * dtime
                   zwt(c)     = zwt(c) + (rsub_top(c) * dtime)/1000._r8/rous
                   aquifer_excess = max(0._r8, wa(c) - 5000._r8)
                   if (aquifer_excess > 0._r8) then
                      do j = nlevbed, 1, -1
                         layer_capacity = max(0._r8, &
                              eff_porosity(c,j)*dz(c,j)*denh2o - h2osoi_liq(c,j))
                         aquifer_bottom_add = min(aquifer_excess, layer_capacity)
                         if (aquifer_bottom_add > 0._r8) then
                            h2osoi_liq(c,j) = h2osoi_liq(c,j) + aquifer_bottom_add
                            aquifer_excess = aquifer_excess - aquifer_bottom_add
                         endif
                         if (aquifer_excess <= 0._r8) exit
                      enddo
                      if (aquifer_excess > 0._r8) then
                         if (h2osfcflag == 1) then
                            h2osfc(c) = h2osfc(c) + aquifer_excess
                         else
                            qflx_qrgwl(c) = qflx_qrgwl(c) + aquifer_excess/dtime
                         endif
                      endif
                      call rediagnose_zwt_from_liq_profile(nlevbed, zengdecker_2009_with_var_soil_thick, &
                           h2osoi_liq(c,1:), eff_porosity(c,1:), dzmm(c,:), zi(c,0:), zwt(c), jwt(c))
                   endif
                   wa(c)  = min(wa(c), 5000._r8)
                end if
             else
                !-- water table within soil layers 1-9  -------------------------------------
                !============================== RSUB_TOP =========================================
                !--  Now remove water via rsub_top
                rsub_top_tot = - rsub_top(c) * dtime
                !should never be positive... but include for completeness
                if(rsub_top_tot > 0.) then !rising water table
#ifndef _OPENACC
                   call endrun(msg="RSUB_TOP IS POSITIVE in Drainage!"//errmsg(__FILE__, __LINE__))
#endif
                else ! deepening water table
                   if (use_vichydro) then
                      wtsub_vic = 0._r8
                      do j = (nlvic(1)+nlvic(2)+1), nlevbed
                         wtsub_vic = wtsub_vic + hk_l(c,j)*dzmm(c,j)
                      end do

                      do j = (nlvic(1)+nlvic(2)+1), nlevbed
                         rsub_top_layer=max(rsub_top_tot, rsub_top_tot*hk_l(c,j)*dzmm(c,j)/wtsub_vic)
                         rsub_top_layer=min(rsub_top_layer,0._r8)
                         if (use_vsfm) rsub_top_layer = 0._r8
                         h2osoi_liq(c,j) = h2osoi_liq(c,j) + rsub_top_layer
                         rsub_top_tot = rsub_top_tot - rsub_top_layer
                      end do
                   else
                      do j = jwt(c)+1, nlevbed
                         ! use analytical expression for specific yield
                         s_y = watsat(c,j) &
                              * ( 1. - (1.+1.e3*zwt(c)/sucsat(c,j))**(-1./bsw(c,j)))
                         s_y=max(s_y,0.02_r8)

                         rsub_top_layer=max(rsub_top_tot,-(s_y*(zi(c,j) - zwt(c))*1.e3))
                         rsub_top_layer=min(rsub_top_layer,0._r8)
                         if (use_vsfm) rsub_top_layer = 0._r8
                         h2osoi_liq(c,j) = h2osoi_liq(c,j) + rsub_top_layer
                         qflx_drain_vr(c,j) = qflx_drain_vr(c,j) - rsub_top_layer

                         rsub_top_tot = rsub_top_tot - rsub_top_layer

                         if (rsub_top_tot >= 0.) then
                            zwt(c) = zwt(c) - rsub_top_layer/s_y/1000._r8

                            exit
                         else
                            zwt(c) = zi(c,j)
                         endif
                      enddo
                   end if

                   !--  remove residual rsub_top  ---------------------------------------------
                   if (zengdecker_2009_with_var_soil_thick) then
                      if (rsub_top_tot < 0.) then
                         rsub_top(c) = rsub_top(c) + rsub_top_tot / dtime
                         rsub_top_tot = 0._r8
                      end if
                   else
                      zwt(c) = zwt(c) - rsub_top_tot/1000._r8/rous
                      wa(c) = wa(c) + rsub_top_tot
                   end if
                endif

                !-- recompute jwt  ---------------------------------------------------------
                ! allow jwt to equal zero when zwt is in top layer
                jwt(c) = nlevbed
                do j = 1,nlevbed
                   if(zwt(c) <= zi(c,j)) then
                      if (zengdecker_2009_with_var_soil_thick .and. zwt(c) == zi(c,nlevbed)) then
                         exit
                      else
                         jwt(c) = j-1
                         exit
                      end if
                   end if
                enddo
             end if! end of jwt if construct

             zwt(c) = max(0.0_r8,zwt(c))
             zwt(c) = min(80._r8,zwt(c))

             rsub_top_after_rsub(c) = rsub_top(c)
             drainage_store_after_rsub(c) = wa(c) + h2osfc(c)
             do j = 1, nlevbed
                drainage_store_after_rsub(c) = drainage_store_after_rsub(c) + h2osoi_liq(c,j) + h2osoi_ice(c,j)
             enddo

          endif

	       end do

	       do fc = 1, num_hydrologyc
	          c = filter_hydrologyc(fc)
	          nlevbed = nlev2bed(c)
	          layer_liq_min = huge(1._r8)
	          layer_excess_max = -huge(1._r8)
	          layer_liq_min_j = 0
	          layer_excess_max_j = 0
	          do j = 1, nlevbed
	             layer_capacity = eff_porosity(c,j)*dzmm(c,j)
	             if (h2osoi_liq(c,j) < layer_liq_min) then
	                layer_liq_min = h2osoi_liq(c,j)
	                layer_liq_min_j = j
	             endif
	             if (h2osoi_liq(c,j) - layer_capacity > layer_excess_max) then
	                layer_excess_max = h2osoi_liq(c,j) - layer_capacity
	                layer_excess_max_j = j
	             endif
	          enddo
	          if (layer_liq_min < layer_water_min_abort) then
	             t = col_pp%topounit(c)
	             g = col_pp%gridcell(c)
	             write(iulog,*)'drainage layer diagnostic before excess bucket'
	             write(iulog,*)'nstep                      = ',nstep_mod
	             write(iulog,*)'column index               = ',c
	             write(iulog,*)'gridcell index             = ',g
	             write(iulog,*)'Latdeg,Londeg              = ',grc_pp%latdeg(g),grc_pp%londeg(g)
	             write(iulog,*)'topounit index             = ',t
	             write(iulog,*)'topounit topo_grc_ind      = ',top_pp%topo_grc_ind(t)
	             write(iulog,*)'topounit weight gridcell   = ',top_pp%wtgcell(t)
	             write(iulog,*)'topounit peat_depth        = ',top_pp%peat_depth(t)
	             write(iulog,*)'topounit is_bog            = ',top_pp%is_bog(t)
	             write(iulog,*)'nlevbed                    = ',nlevbed
	             write(iulog,*)'jwt                        = ',jwt(c)
	             write(iulog,*)'layer liq min              = ',layer_liq_min
	             write(iulog,*)'layer liq min j            = ',layer_liq_min_j
	             write(iulog,*)'layer excess max           = ',layer_excess_max
	             write(iulog,*)'layer excess max j         = ',layer_excess_max_j
	             write(iulog,*)'rsub_top                   = ',rsub_top(c)
	             write(iulog,*)'rsub_top d                 = ',rsub_top(c)*86400._r8
	             write(iulog,*)'qflx_rsub_sat              = ',qflx_rsub_sat(c)
	             write(iulog,*)'qflx_rsub_sat d            = ',qflx_rsub_sat(c)*86400._r8
	             write(iulog,*)'qflx_lat_aqu               = ',qflx_lat_aqu(c)
	             write(iulog,*)'qflx_lat_aqu d             = ',qflx_lat_aqu(c)*86400._r8
	             write(iulog,*)'qcharge                    = ',qcharge(c)
	             write(iulog,*)'qcharge d                  = ',qcharge(c)*86400._r8
	             write(iulog,*)'wa                         = ',wa(c)
	             write(iulog,*)'zwt                        = ',zwt(c)
	             write(iulog,*)'zwt_perched                = ',zwt_perched(c)
	             write(iulog,*)'h2osfc                     = ',h2osfc(c)
	             write(iulog,*)'drainage store beg         = ',drainage_store_beg(c)
	             write(iulog,*)'drainage store after rsub  = ',drainage_store_after_rsub(c)
	             write(iulog,*)'h2osoi_liq                 = ',h2osoi_liq(c,1:nlevbed)
	             write(iulog,*)'h2osoi_ice                 = ',h2osoi_ice(c,1:nlevbed)
	             write(iulog,*)'eff_porosity               = ',eff_porosity(c,1:nlevbed)
	             write(iulog,*)'watsat                     = ',watsat(c,1:nlevbed)
	             write(iulog,*)'dzmm                       = ',dzmm(c,1:nlevbed)
	             write(iulog,*)'qflx_drain_vr              = ',qflx_drain_vr(c,1:nlevbed)
	             write(iulog,*)'qflx_adv                   = ',qflx_adv(c,1:nlevbed)
	             call endrun(msg='Drainage layer water below minimum before excess bucket'//errMsg(__FILE__, __LINE__))
	          endif
	       end do

	       !  excessive water above saturation added to the above unsaturated layer like a bucket
	       !  if column fully saturated, excess water goes to runoff

	       do fc = 1, num_hydrologyc
	          c = filter_hydrologyc(fc)
            nlevbed = nlev2bed(c)
          do j = nlevbed,2,-1
             xsi(c)            = max(h2osoi_liq(c,j)-eff_porosity(c,j)*dzmm(c,j),0._r8)
             if (use_vsfm) then
                xsi(c) = 0._r8
             else
                h2osoi_liq(c,j)   = min(eff_porosity(c,j)*dzmm(c,j), h2osoi_liq(c,j))
                h2osoi_liq(c,j-1) = h2osoi_liq(c,j-1) + xsi(c)
	             endif
	          end do
	       end do

	       do fc = 1, num_hydrologyc
	          c = filter_hydrologyc(fc)
	          nlevbed = nlev2bed(c)
	          layer_liq_min = huge(1._r8)
	          layer_excess_max = -huge(1._r8)
	          layer_liq_min_j = 0
	          layer_excess_max_j = 0
	          do j = 1, nlevbed
	             layer_capacity = eff_porosity(c,j)*dzmm(c,j)
	             if (h2osoi_liq(c,j) < layer_liq_min) then
	                layer_liq_min = h2osoi_liq(c,j)
	                layer_liq_min_j = j
	             endif
	             if (h2osoi_liq(c,j) - layer_capacity > layer_excess_max) then
	                layer_excess_max = h2osoi_liq(c,j) - layer_capacity
	                layer_excess_max_j = j
	             endif
	          enddo
	          if (layer_liq_min < layer_water_min_abort .or. &
	               (layer_excess_max > layer_water_excess_abort .and. layer_excess_max_j > 1)) then
	             t = col_pp%topounit(c)
	             g = col_pp%gridcell(c)
	             write(iulog,*)'drainage layer diagnostic after excess bucket'
	             write(iulog,*)'nstep                      = ',nstep_mod
	             write(iulog,*)'column index               = ',c
	             write(iulog,*)'gridcell index             = ',g
	             write(iulog,*)'Latdeg,Londeg              = ',grc_pp%latdeg(g),grc_pp%londeg(g)
	             write(iulog,*)'topounit index             = ',t
	             write(iulog,*)'topounit topo_grc_ind      = ',top_pp%topo_grc_ind(t)
	             write(iulog,*)'topounit weight gridcell   = ',top_pp%wtgcell(t)
	             write(iulog,*)'topounit peat_depth        = ',top_pp%peat_depth(t)
	             write(iulog,*)'topounit is_bog            = ',top_pp%is_bog(t)
	             write(iulog,*)'nlevbed                    = ',nlevbed
	             write(iulog,*)'jwt                        = ',jwt(c)
	             write(iulog,*)'layer liq min              = ',layer_liq_min
	             write(iulog,*)'layer liq min j            = ',layer_liq_min_j
	             write(iulog,*)'layer excess max           = ',layer_excess_max
	             write(iulog,*)'layer excess max j         = ',layer_excess_max_j
	             write(iulog,*)'rsub_top                   = ',rsub_top(c)
	             write(iulog,*)'rsub_top d                 = ',rsub_top(c)*86400._r8
	             write(iulog,*)'qflx_rsub_sat              = ',qflx_rsub_sat(c)
	             write(iulog,*)'qflx_rsub_sat d            = ',qflx_rsub_sat(c)*86400._r8
	             write(iulog,*)'qflx_lat_aqu               = ',qflx_lat_aqu(c)
	             write(iulog,*)'qflx_lat_aqu d             = ',qflx_lat_aqu(c)*86400._r8
	             write(iulog,*)'qcharge                    = ',qcharge(c)
	             write(iulog,*)'qcharge d                  = ',qcharge(c)*86400._r8
	             write(iulog,*)'wa                         = ',wa(c)
	             write(iulog,*)'zwt                        = ',zwt(c)
	             write(iulog,*)'zwt_perched                = ',zwt_perched(c)
	             write(iulog,*)'h2osfc                     = ',h2osfc(c)
	             write(iulog,*)'drainage store beg         = ',drainage_store_beg(c)
	             write(iulog,*)'drainage store after rsub  = ',drainage_store_after_rsub(c)
	             write(iulog,*)'h2osoi_liq                 = ',h2osoi_liq(c,1:nlevbed)
	             write(iulog,*)'h2osoi_ice                 = ',h2osoi_ice(c,1:nlevbed)
	             write(iulog,*)'eff_porosity               = ',eff_porosity(c,1:nlevbed)
	             write(iulog,*)'watsat                     = ',watsat(c,1:nlevbed)
	             write(iulog,*)'dzmm                       = ',dzmm(c,1:nlevbed)
	             write(iulog,*)'qflx_drain_vr              = ',qflx_drain_vr(c,1:nlevbed)
	             write(iulog,*)'qflx_adv                   = ',qflx_adv(c,1:nlevbed)
	             call endrun(msg='Drainage layer water out of bounds after excess bucket'//errMsg(__FILE__, __LINE__))
	          endif
	       end do

	       do fc = 1, num_hydrologyc
	          c = filter_hydrologyc(fc)

	          !scs: watmin addition to fix water balance errors
          xs1(c)          = max(max(h2osoi_liq(c,1)-watmin,0._r8)- &
               max(0._r8,(pondmx+watsat(c,1)*dzmm(c,1)-h2osoi_ice(c,1)-watmin)),0._r8)
          if (use_vsfm) xs1(c) = 0._r8
          h2osoi_liq(c,1) = h2osoi_liq(c,1) - xs1(c)

          if (lun_pp%urbpoi(col_pp%landunit(c))) then
             qflx_rsub_sat(c)     = xs1(c) / dtime
          else
             if(h2osfcflag == 1) then
                ! send this water up to h2osfc rather than sending to drainage
                h2osfc(c) = h2osfc(c) + xs1(c)
                qflx_rsub_sat(c)     = 0._r8
             else
	                ! use original code to send water to drainage (non-h2osfc case)
	                qflx_rsub_sat(c)     = xs1(c) / dtime
	             endif
	             qflx_drain_vr(c,1) = qflx_drain_vr(c,1) - xs1(c)
	          endif

	          if (h2osfc(c) > h2osfc_dump_threshold) then
	             h2osfc_before_dump = h2osfc(c)
	             qflx_qrgwl_before_dump = qflx_qrgwl(c)
	             h2osfc_excess_dump = h2osfc(c) - h2osfc_dump_threshold
	             h2osfc(c) = h2osfc_dump_threshold
	             qflx_qrgwl(c) = qflx_qrgwl(c) + h2osfc_excess_dump/dtime
	             t = col_pp%topounit(c)
	             g = col_pp%gridcell(c)
	             write(iulog,*)'WARNING: h2osfc dump nstep=',nstep_mod,' column=',c,' gridcell=',g, &
	                  ' lat=',grc_pp%latdeg(g),' lon=',grc_pp%londeg(g),' topounit=',t, &
	                  ' topo_grc_ind=',top_pp%topo_grc_ind(t),' wtgcell=',top_pp%wtgcell(t), &
	                  ' peat_depth=',top_pp%peat_depth(t),' is_bog=',top_pp%is_bog(t), &
	                  ' dumped_mm=',h2osfc_excess_dump, &
	                  ' h2osfc_start=',h2osfc_drainage_beg(c),' h2osfc_pre_dump=',h2osfc_before_dump, &
	                  ' xs1=',xs1(c),' qflx_qrgwl_pre=',qflx_qrgwl_before_dump, &
	                  ' qflx_lat_aqu_entry_dt=',qflx_lat_aqu_beg(c)*dtime, &
	                  ' qflx_lat_aqu_now_dt=',qflx_lat_aqu(c)*dtime, &
	                  ' qcharge_dt=',qcharge(c)*dtime,' wa=',wa(c),' zwt=',zwt(c), &
	                  ' qflx_qrgwl=',qflx_qrgwl(c),' h2osfc=',h2osfc(c)
	          endif

	          if (abs(qflx_drain_vr(c,1)) > layer_water_excess_abort .and. h2osfcflag /= 1) then
	             nlevbed = nlev2bed(c)
	             t = col_pp%topounit(c)
	             g = col_pp%gridcell(c)
	             write(iulog,*)'drainage top excess diagnostic'
	             write(iulog,*)'nstep                      = ',nstep_mod
	             write(iulog,*)'column index               = ',c
	             write(iulog,*)'gridcell index             = ',g
	             write(iulog,*)'Latdeg,Londeg              = ',grc_pp%latdeg(g),grc_pp%londeg(g)
	             write(iulog,*)'topounit index             = ',t
	             write(iulog,*)'topounit topo_grc_ind      = ',top_pp%topo_grc_ind(t)
	             write(iulog,*)'topounit weight gridcell   = ',top_pp%wtgcell(t)
	             write(iulog,*)'topounit peat_depth        = ',top_pp%peat_depth(t)
	             write(iulog,*)'topounit is_bog            = ',top_pp%is_bog(t)
	             write(iulog,*)'nlevbed                    = ',nlevbed
	             write(iulog,*)'jwt                        = ',jwt(c)
	             write(iulog,*)'xs1 top excess             = ',xs1(c)
	             write(iulog,*)'qflx_rsub_sat              = ',qflx_rsub_sat(c)
	             write(iulog,*)'qflx_rsub_sat d            = ',qflx_rsub_sat(c)*86400._r8
	             write(iulog,*)'qflx_qrgwl                 = ',qflx_qrgwl(c)
	             write(iulog,*)'qflx_qrgwl d               = ',qflx_qrgwl(c)*86400._r8
	             write(iulog,*)'h2osfc                     = ',h2osfc(c)
	             write(iulog,*)'wa                         = ',wa(c)
	             write(iulog,*)'zwt                        = ',zwt(c)
	             write(iulog,*)'zwt_perched                = ',zwt_perched(c)
	             write(iulog,*)'drainage store beg         = ',drainage_store_beg(c)
	             write(iulog,*)'drainage store after rsub  = ',drainage_store_after_rsub(c)
	             write(iulog,*)'h2osoi_liq                 = ',h2osoi_liq(c,1:nlevbed)
	             write(iulog,*)'h2osoi_ice                 = ',h2osoi_ice(c,1:nlevbed)
	             write(iulog,*)'eff_porosity               = ',eff_porosity(c,1:nlevbed)
	             write(iulog,*)'watsat                     = ',watsat(c,1:nlevbed)
	             write(iulog,*)'dzmm                       = ',dzmm(c,1:nlevbed)
	             write(iulog,*)'qflx_drain_vr              = ',qflx_drain_vr(c,1:nlevbed)
	             write(iulog,*)'qflx_adv                   = ',qflx_adv(c,1:nlevbed)
	             call endrun(msg='Drainage top excess created implausible surface water'//errMsg(__FILE__, __LINE__))
	          endif

	          if (use_vsfm) qflx_rsub_sat(c) = 0._r8

	          ! add in ice check
          xs1(c)          = max(max(h2osoi_ice(c,1),0._r8)-max(0._r8,(pondmx+watsat(c,1)*dzmm(c,1)-h2osoi_liq(c,1))),0._r8)
          h2osoi_ice(c,1) = min(max(0._r8,pondmx+watsat(c,1)*dzmm(c,1)-h2osoi_liq(c,1)), h2osoi_ice(c,1))
          if ( (lun_pp%itype(col_pp%landunit(c)) == istice .or. lun_pp%itype(col_pp%landunit(c)) == istice_mec) .or. (.not. use_firn_percolation_and_compaction)) then      
                qflx_snwcp_ice(c) = qflx_snwcp_ice(c) + xs1(c) / dtime
          else
                qflx_ice_runoff_xs(c) = qflx_ice_runoff_xs(c) + xs1(c) / dtime
          endif
       end do

       ! Limit h2osoi_liq to be greater than or equal to watmin.
       ! Get water needed to bring h2osoi_liq equal watmin from lower layer.
       ! If insufficient water in soil layers, get from aquifer water

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          nlevbed = nlev2bed(c)
          do j = 1, nlevbed-1
             if (h2osoi_liq(c,j) < watmin) then
                xs(c) = watmin - h2osoi_liq(c,j)
                ! deepen water table if water is passed from below zwt layer
                if(j == jwt(c)) then
                   zwt(c) = zwt(c) + xs(c)/eff_porosity(c,j)/1000._r8
                endif
             else
                xs(c) = 0._r8
             end if
             h2osoi_liq(c,j  ) = h2osoi_liq(c,j  ) + xs(c)
             h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) - xs(c)
          end do
       end do

       ! Get water for bottom layer from layers above if possible
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          nlevbed = nlev2bed(c)
          j = nlevbed
          bottom_liq_before_refill = h2osoi_liq(c,j)
          bottom_refill_need_initial = 0._r8
          donor_available_total = 0._r8
          donor_liq_total = 0._r8
          if (h2osoi_liq(c,j) < watmin) then
             xs(c) = watmin-h2osoi_liq(c,j)
             bottom_refill_need_initial = xs(c)
             searchforwater: do i = nlevbed-1, 1, -1
                donor_liq_total = donor_liq_total + h2osoi_liq(c,i)
                available_h2osoi_liq = max(h2osoi_liq(c,i)-watmin-xs(c),0._r8)
                donor_available_total = donor_available_total + available_h2osoi_liq
                if (available_h2osoi_liq >= xs(c)) then
                   h2osoi_liq(c,j) = h2osoi_liq(c,j) + xs(c)
                   h2osoi_liq(c,i) = h2osoi_liq(c,i) - xs(c)
                   xs(c) = 0._r8
                   exit searchforwater
                else
                   h2osoi_liq(c,j) = h2osoi_liq(c,j) + available_h2osoi_liq
                   h2osoi_liq(c,i) = h2osoi_liq(c,i) - available_h2osoi_liq
                   xs(c) = xs(c) - available_h2osoi_liq
                end if
             end do searchforwater
          else
             xs(c) = 0._r8
          end if

          watmin_refill_flux = xs(c)/dtime
          if (watmin_refill_flux > watmin_refill_abort_flux .or. &
               rsub_top(c) - watmin_refill_flux < -watmin_refill_abort_flux) then
             t = col_pp%topounit(c)
             g = col_pp%gridcell(c)
             write(iulog,*)'watmin refill diagnostic'
             write(iulog,*)'nstep                      = ',nstep_mod
             write(iulog,*)'column index               = ',c
             write(iulog,*)'gridcell index             = ',g
             write(iulog,*)'Latdeg,Londeg              = ',grc_pp%latdeg(g),grc_pp%londeg(g)
             write(iulog,*)'topounit index             = ',t
             write(iulog,*)'topounit topo_grc_ind      = ',top_pp%topo_grc_ind(t)
             write(iulog,*)'topounit weight gridcell   = ',top_pp%wtgcell(t)
             write(iulog,*)'topounit peat_depth        = ',top_pp%peat_depth(t)
             write(iulog,*)'topounit is_bog            = ',top_pp%is_bog(t)
             write(iulog,*)'dtime                      = ',dtime
             write(iulog,*)'nlevbed                    = ',nlevbed
             write(iulog,*)'jwt                        = ',jwt(c)
             write(iulog,*)'watmin                     = ',watmin
             write(iulog,*)'bottom layer               = ',j
             write(iulog,*)'bottom liq before refill   = ',bottom_liq_before_refill
             write(iulog,*)'bottom liq after search    = ',h2osoi_liq(c,j)
             write(iulog,*)'initial refill need        = ',bottom_refill_need_initial
             write(iulog,*)'unresolved refill xs       = ',xs(c)
             write(iulog,*)'unresolved refill flux     = ',watmin_refill_flux
             write(iulog,*)'unresolved refill flux d   = ',watmin_refill_flux*86400._r8
             write(iulog,*)'donor liquid total         = ',donor_liq_total
             write(iulog,*)'donor available total      = ',donor_available_total
             write(iulog,*)'rsub_top before refill     = ',rsub_top(c)
             write(iulog,*)'rsub_top after refill      = ',rsub_top(c)-watmin_refill_flux
             write(iulog,*)'rsub_top after refill d    = ',(rsub_top(c)-watmin_refill_flux)*86400._r8
             write(iulog,*)'qflx_rsub_sat              = ',qflx_rsub_sat(c)
             write(iulog,*)'qflx_rsub_sat d            = ',qflx_rsub_sat(c)*86400._r8
             write(iulog,*)'qflx_lat_aqu               = ',qflx_lat_aqu(c)
             write(iulog,*)'qflx_lat_aqu d             = ',qflx_lat_aqu(c)*86400._r8
             write(iulog,*)'qcharge                    = ',qcharge(c)
             write(iulog,*)'qcharge d                  = ',qcharge(c)*86400._r8
             write(iulog,*)'wa                         = ',wa(c)
             write(iulog,*)'zwt                        = ',zwt(c)
             write(iulog,*)'zwt_perched                = ',zwt_perched(c)
             write(iulog,*)'h2osfc                     = ',h2osfc(c)
             write(iulog,*)'drainage store beg         = ',drainage_store_beg(c)
             write(iulog,*)'drainage store after rsub  = ',drainage_store_after_rsub(c)
             write(iulog,*)'h2osoi_liq                 = ',h2osoi_liq(c,1:nlevbed)
             write(iulog,*)'h2osoi_ice                 = ',h2osoi_ice(c,1:nlevbed)
             write(iulog,*)'eff_porosity               = ',eff_porosity(c,1:nlevbed)
             write(iulog,*)'watsat                     = ',watsat(c,1:nlevbed)
             write(iulog,*)'dzmm                       = ',dzmm(c,1:nlevbed)
             write(iulog,*)'qflx_drain_vr              = ',qflx_drain_vr(c,1:nlevbed)
             write(iulog,*)'qflx_adv                   = ',qflx_adv(c,1:nlevbed)
             call endrun(msg='Watmin refill correction produced implausible negative drainage'//errMsg(__FILE__, __LINE__))
          end if

          ! Needed in case there is no water to be found
          h2osoi_liq(c,j) = h2osoi_liq(c,j) + xs(c)
          ! Instead of removing water from aquifer where it eventually
          ! shows up as excess drainage to the ocean, take it back out of
          ! drainage
          rsub_top(c) = rsub_top(c) - xs(c)/dtime
          rsub_top_final(c) = rsub_top(c)

       end do

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! Sub-surface runoff and drainage

          qflx_drain(c) = qflx_rsub_sat(c) + rsub_top(c)

          drainage_store_end(c) = wa(c) + h2osfc(c)
          do j = 1, nlev2bed(c)
             drainage_store_end(c) = drainage_store_end(c) + h2osoi_liq(c,j) + h2osoi_ice(c,j)
          enddo

          qflx_lat_aqu_delta = qflx_lat_aqu(c) - qflx_lat_aqu_beg(c)
          drainage_diag_resid_no_lat = drainage_store_end(c) - drainage_store_beg(c) &
               + (qflx_drain(c) + qflx_drain_perched(c) + qflx_qrgwl(c) &
               + qflx_snwcp_ice(c) + qflx_ice_runoff_xs(c)) * dtime
          drainage_diag_resid = drainage_diag_resid_no_lat - qflx_lat_aqu_delta*dtime
          print_drainage_diag = debug_drainage_diag .and. &
               (top_pp%is_bog(col_pp%topounit(c)) &
               .and. top_pp%topo_grc_ind(col_pp%topounit(c)) == 3 &
               .and. abs(top_pp%wtgcell(col_pp%topounit(c)) - 2.1752511754418422e-2_r8) < 1.e-10_r8 &
               .and. abs(top_pp%peat_depth(col_pp%topounit(c)) - 1.8641646375629948_r8) < 1.e-8_r8 &
               .and. nstep_mod >= 7670 .and. nstep_mod <= 7680 &
               .and. (abs(drainage_diag_resid) > 1.e-7_r8 .or. abs(drainage_diag_resid_no_lat) > 1.e-7_r8))
          if (print_drainage_diag) then
             write(iulog,*)'drainage diagnostic'
             write(iulog,*)'nstep                      = ',nstep_mod
             write(iulog,*)'column index               = ',c
             write(iulog,*)'gridcell index             = ',col_pp%gridcell(c)
             write(iulog,*)'topounit index             = ',col_pp%topounit(c)
             write(iulog,*)'topounit topo_grc_ind      = ',top_pp%topo_grc_ind(col_pp%topounit(c))
             write(iulog,*)'nlevbed                    = ',nlev2bed(c)
             write(iulog,*)'jwt                        = ',jwt(c)
             write(iulog,*)'zwt                        = ',zwt(c)
             write(iulog,*)'wa                         = ',wa(c)
             write(iulog,*)'h2osfc                     = ',h2osfc(c)
             write(iulog,*)'h2osoi_liq bottom          = ',h2osoi_liq(c,nlev2bed(c))
             write(iulog,*)'h2osoi_ice bottom          = ',h2osoi_ice(c,nlev2bed(c))
             write(iulog,*)'drainage store beg         = ',drainage_store_beg(c)
             write(iulog,*)'drainage store after rsub  = ',drainage_store_after_rsub(c)
             write(iulog,*)'drainage store end         = ',drainage_store_end(c)
             write(iulog,*)'drainage store delta       = ',drainage_store_end(c)-drainage_store_beg(c)
             write(iulog,*)'rsub_top after rsub        = ',rsub_top_after_rsub(c)
             write(iulog,*)'rsub_top final             = ',rsub_top_final(c)
             write(iulog,*)'qflx_rsub_sat              = ',qflx_rsub_sat(c)
             write(iulog,*)'qflx_drain                 = ',qflx_drain(c)
             write(iulog,*)'qflx_drain_perched         = ',qflx_drain_perched(c)
             write(iulog,*)'qflx_qrgwl                 = ',qflx_qrgwl(c)
             write(iulog,*)'qflx_lat_aqu entry         = ',qflx_lat_aqu_beg(c)
             write(iulog,*)'qflx_lat_aqu final         = ',qflx_lat_aqu(c)
             write(iulog,*)'qflx_lat_aqu delta         = ',qflx_lat_aqu_delta
             write(iulog,*)'qflx_lat_aqu delta dt      = ',qflx_lat_aqu_delta*dtime
             write(iulog,*)'qflx_snwcp_ice             = ',qflx_snwcp_ice(c)
             write(iulog,*)'qflx_ice_runoff_xs          = ',qflx_ice_runoff_xs(c)
             write(iulog,*)'drainage flux dt           = ',qflx_drain(c)*dtime
             write(iulog,*)'drainage resid no lat      = ',drainage_diag_resid_no_lat
             write(iulog,*)'drainage diagnostic resid  = ',drainage_diag_resid
          endif

          ! Set imbalance for snow capping without overwriting other qrgwl
          ! terms such as peat aquifer overflow.
          qflx_qrgwl(c) = qflx_qrgwl(c) + qflx_snwcp_liq(c)

       end do

       ! No drainage for urban columns (except for pervious road as computed above)

       do fc = 1, num_urbanc
          c = filter_urbanc(fc)
          if (col_pp%itype(c) /= icol_road_perv) then
             qflx_drain(c) = 0._r8
             ! This must be done for roofs and impervious road (walls will be zero)
             qflx_qrgwl(c) = qflx_snwcp_liq(c)
          end if
       end do

     end associate

   end subroutine Drainage

   !-----------------------------------------------------------------------
   subroutine DrainageVSFM(bounds, num_hydrologyc, filter_hydrologyc, num_urbanc, filter_urbanc,  &
         soilhydrology_vars, soilstate_vars, dtime)
     !
     ! !DESCRIPTION:
     ! Calculate subsurface drainage
     !
     ! !USES:
      !$acc routine seq
     use elm_varpar       , only : nlevsoi, nlevgrnd, nlayer, nlayert
     use elm_varcon       , only : pondmx, tfrz, watmin,rpi, secspday, nlvic
     use column_varcon    , only : icol_roof, icol_road_imperv, icol_road_perv
     use elm_varctl       , only : use_vsfm
     use pftvarcon        , only : rsub_top_globalmax
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds
     integer                  , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
     integer                  , intent(in)    :: num_urbanc           ! number of column urban points in column filter
     integer                  , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
     integer                  , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
     type(soilstate_type)     , intent(in)    :: soilstate_vars
     type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
     real(r8), intent(in)  :: dtime
     !
     ! !LOCAL VARIABLES:
     character(len=32) :: subname = 'Drainage'           ! subroutine name
     integer  :: c,j,fc,i                                ! indices
     real(r8) :: xs(bounds%begc:bounds%endc)             ! water needed to bring soil moisture to watmin (mm)
     real(r8) :: dzmm(bounds%begc:bounds%endc,1:nlevgrnd) ! layer thickness (mm)
     integer  :: jwt(bounds%begc:bounds%endc)            ! index of the soil layer right above the water table (-)
     real(r8) :: rsub_bot(bounds%begc:bounds%endc)       ! subsurface runoff - bottom drainage (mm/s)
     real(r8) :: rsub_top(bounds%begc:bounds%endc)       ! subsurface runoff - topographic control (mm/s)
     real(r8) :: fff(bounds%begc:bounds%endc)            ! decay factor (m-1)
     real(r8) :: xsi(bounds%begc:bounds%endc)            ! excess soil water above saturation at layer i (mm)
     real(r8) :: xsia(bounds%begc:bounds%endc)           ! available pore space at layer i (mm)
     real(r8) :: xs1(bounds%begc:bounds%endc)            ! excess soil water above saturation at layer 1 (mm)
     real(r8) :: smpfz(1:nlevgrnd)                        ! matric potential of layer right above water table (mm)
     real(r8) :: wtsub                                   ! summation of hk*dzmm for layers below water table (mm**2/s)
     real(r8) :: rous                                    ! aquifer yield (-)
     real(r8) :: wh                                      ! smpfz(jwt)-z(jwt) (mm)
     real(r8) :: wh_zwt                                  ! water head at the water table depth (mm)
     real(r8) :: ws                                      ! summation of pore space of layers below water table (mm)
     real(r8) :: s_node                                  ! soil wetness (-)
     real(r8) :: dzsum                                   ! summation of dzmm of layers below water table (mm)
     real(r8) :: icefracsum                              ! summation of icefrac*dzmm of layers below water table (-)
     real(r8) :: fracice_rsub(bounds%begc:bounds%endc)   ! fractional impermeability of soil layers (-)
     real(r8) :: ka                                      ! hydraulic conductivity of the aquifer (mm/s)
     real(r8) :: dza                                     ! fff*(zwt-z(jwt)) (-)
     real(r8) :: available_h2osoi_liq                    ! available soil liquid water in a layer
     real(r8) :: rsub_top_max
     real(r8) :: h2osoi_vol
     real(r8) :: imped
     real(r8) :: rsub_top_tot
     real(r8) :: rsub_top_layer
     real(r8) :: qcharge_tot
     real(r8) :: qcharge_layer
     real(r8) :: theta_unsat
     real(r8) :: f_unsat
     real(r8) :: s_y
     integer  :: k,k_frz,k_perch
     real(r8) :: sat_lev
     real(r8) :: s1
     real(r8) :: s2
     real(r8) :: m
     real(r8) :: b
     real(r8) :: q_perch
     real(r8) :: q_perch_max
     real(r8) :: vol_ice
     real(r8) :: dsmax_tmp(bounds%begc:bounds%endc)       ! temporary variable for ARNO subsurface runoff calculation
     real(r8) :: rsub_tmp                 ! temporary variable for ARNO subsurface runoff calculation
     real(r8) :: frac                     ! temporary variable for ARNO subsurface runoff calculation
     real(r8) :: rel_moist                ! relative moisture, temporary variable
     real(r8) :: wtsub_vic                ! summation of hk*dzmm for layers in the third VIC layer
     integer  :: idx                      ! 1D index for VSFM
     !-----------------------------------------------------------------------

     associate(                                                            &
          z                  =>    col_pp%z                                 , & ! Input:  [real(r8) (:,:) ] layer depth (m)
          zi                 =>    col_pp%zi                                , & ! Input:  [real(r8) (:,:) ] interface level below a "z" level (m)
          dz                 =>    col_pp%dz                                , & ! Input:  [real(r8) (:,:) ] layer depth (m)
          snl                =>    col_pp%snl                               , & ! Input:  [integer  (:)   ] number of snow layers

          t_soisno           =>    col_es%t_soisno         , & ! Input:  [real(r8) (:,:) ] soil temperature (Kelvin)

          h2osfc             =>    col_ws%h2osfc            , & ! Input:  [real(r8) (:)   ] surface water (mm)

          bsw                =>    soilstate_vars%bsw_col                , & ! Input:  [real(r8) (:,:) ] Clapp and Hornberger "b"
          hksat              =>    soilstate_vars%hksat_col              , & ! Input:  [real(r8) (:,:) ] hydraulic conductivity at saturation (mm H2O /s)
          sucsat             =>    soilstate_vars%sucsat_col             , & ! Input:  [real(r8) (:,:) ] minimum soil suction (mm)
          watsat             =>    soilstate_vars%watsat_col             , & ! Input:  [real(r8) (:,:) ] volumetric soil water at saturation (porosity)
          eff_porosity       =>    soilstate_vars%eff_porosity_col       , & ! Input:  [real(r8) (:,:) ] effective porosity = porosity - vol_ice
          hk_l               =>    soilstate_vars%hk_l_col               , & ! Input:  [real(r8) (:,:) ] hydraulic conductivity (mm/s)

          depth              =>    soilhydrology_vars%depth_col          , & ! Input:  [real(r8) (:,:) ] VIC soil depth
          c_param            =>    soilhydrology_vars%c_param_col        , & ! Input:  [real(r8) (:)   ] baseflow exponent (Qb)
          Dsmax              =>    soilhydrology_vars%dsmax_col          , & ! Input:  [real(r8) (:)   ] max. velocity of baseflow (mm/day)
          max_moist          =>    soilhydrology_vars%max_moist_col      , & ! Input:  [real(r8) (:,:) ] maximum soil moisture (ice + liq)
          moist              =>    soilhydrology_vars%moist_col          , & ! Input:  [real(r8) (:,:) ] soil layer moisture (mm)
          Ds                 =>    soilhydrology_vars%ds_col             , & ! Input:  [real(r8) (:)   ] fracton of Dsmax where non-linear baseflow begins
          Wsvic              =>    soilhydrology_vars%Wsvic_col          , & ! Input:  [real(r8) (:)   ] fraction of maximum soil moisutre where non-liear base flow occurs
          icefrac            =>    soilhydrology_vars%icefrac_col        , & ! Output: [real(r8) (:,:) ] fraction of ice in layer
          hkdepth            =>    soilhydrology_vars%hkdepth_col        , & ! Input:  [real(r8) (:)   ] decay factor (m)
          frost_table        =>    soilhydrology_vars%frost_table_col    , & ! Input:  [real(r8) (:)   ] frost table depth (m)
          zwt                =>    soilhydrology_vars%zwt_col            , & ! Input:  [real(r8) (:)   ] water table depth (m)
          zwt_perched        =>    soilhydrology_vars%zwt_perched_col    , & ! Input:  [real(r8) (:)   ] perched water table depth (m)
          wa                 =>    soilhydrology_vars%wa_col             , & ! Input:  [real(r8) (:)   ] water in the unconfined aquifer (mm)
          ice                =>    soilhydrology_vars%ice_col            , & ! Input:  [real(r8) (:,:) ] soil layer moisture (mm)
          qcharge            =>    soilhydrology_vars%qcharge_col        , & ! Input:  [real(r8) (:)   ] aquifer recharge rate (mm/s)
          origflag           =>    soilhydrology_vars%origflag           , & ! Input:  logical
          h2osfcflag         =>    soilhydrology_vars%h2osfcflag         , & ! Input:  logical

          qflx_snwcp_liq     =>    col_wf%qflx_snwcp_liq     , & ! Output: [real(r8) (:)   ] excess rainfall due to snow capping (mm H2O /s) [+]
          qflx_snwcp_ice     =>    col_wf%qflx_snwcp_ice     , & ! Output: [real(r8) (:)   ] excess snowfall due to snow capping (mm H2O /s) [+]
          qflx_dew_grnd      =>    col_wf%qflx_dew_grnd      , & ! Output: [real(r8) (:)   ] ground surface dew formation (mm H2O /s) [+]
          qflx_dew_snow      =>    col_wf%qflx_dew_snow      , & ! Output: [real(r8) (:)   ] surface dew added to snow pack (mm H2O /s) [+]
          qflx_sub_snow      =>    col_wf%qflx_sub_snow      , & ! Output: [real(r8) (:)   ] sublimation rate from snow pack (mm H2O /s) [+]
          qflx_drain         =>    col_wf%qflx_drain         , & ! Output: [real(r8) (:)   ] sub-surface runoff (mm H2O /s)
          qflx_qrgwl         =>    col_wf%qflx_qrgwl         , & ! Output: [real(r8) (:)   ] qflx_surf at glaciers, wetlands, lakes (mm H2O /s)
          qflx_rsub_sat      =>    col_wf%qflx_rsub_sat      , & ! Output: [real(r8) (:)   ] soil saturation excess [mm h2o/s]
          qflx_drain_perched =>    col_wf%qflx_drain_perched , & ! Output: [real(r8) (:)   ] perched wt sub-surface runoff (mm H2O /s)

          mflx_drain_perched_col_1d =>    col_wf%mflx_drain_perched_1d   , & ! Input:  [real(r8) (:)   ]  drainage from perched water table (kg H2O /s)

          h2osoi_liq         =>    col_ws%h2osoi_liq        , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
          h2osoi_ice         =>    col_ws%h2osoi_ice          & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)
          )

       ! Convert layer thicknesses from m to mm

       do j = 1,nlevgrnd
          do fc = 1, num_hydrologyc
             c = filter_hydrologyc(fc)
             dzmm(c,j) = dz(c,j)*1.e3_r8

             vol_ice = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
             icefrac(c,j) = min(1._r8,vol_ice/watsat(c,j))
          end do
       end do

       ! Initial set

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          qflx_drain(c)          = 0._r8
          rsub_bot(c)            = 0._r8
          qflx_rsub_sat(c)       = 0._r8
          rsub_top(c)            = 0._r8
          fracice_rsub(c)        = 0._r8
          qflx_qrgwl(c)          = 0._r8
          qflx_drain_perched(c)  = 0._r8
       end do

       mflx_drain_perched_col_1d(:) = 0._r8

       ! The layer index of the first unsaturated layer, i.e., the layer right above
       ! the water table

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          jwt(c) = nlevgrnd
          ! allow jwt to equal zero when zwt is in top layer
          do j = 1,nlevgrnd
             if(zwt(c) <= zi(c,j)) then
                jwt(c) = j-1
                exit
             end if
          enddo
       end do

       rous = 0.2_r8

       !==  BASEFLOW ==================================================
       ! perched water table code
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          !  specify maximum drainage rate
          q_perch_max = 1.e-5_r8 * sin(col_pp%topo_slope(c) * (rpi/180._r8))

          ! if layer containing water table is frozen, compute the following:
          !     frost table, perched water table, and drainage from perched saturated layer

          ! define frost table as first frozen layer with unfrozen layer above it
          if(t_soisno(c,1) > tfrz) then
             k_frz=nlevsoi
          else
             k_frz=1
          endif

          do k=2, nlevsoi
             if (t_soisno(c,k-1) > tfrz .and. t_soisno(c,k) <= tfrz) then
                k_frz=k
                exit
             endif
          enddo

          frost_table(c)=z(c,k_frz)

          ! initialize perched water table to frost table, and qflx_drain_perched(c) to zero
          zwt_perched(c)=frost_table(c)
          qflx_drain_perched(c) = 0._r8

          !===================  water table above frost table  =============================
          ! if water table is above frost table, do not use topmodel baseflow formulation
          if (zwt(c) < frost_table(c) .and. t_soisno(c,k_frz) <= tfrz &
               .and. origflag == 0) then

             ! compute drainage from perched saturated region
             wtsub = 0._r8
             q_perch = 0._r8
             do k = jwt(c)+1, k_frz
                imped=10._r8**(-e_ice*(0.5_r8*(icefrac(c,k)+icefrac(c,min(nlevsoi, k+1)))))
                q_perch = q_perch + imped*hksat(c,k)*dzmm(c,k)
                wtsub = wtsub + dzmm(c,k)
             end do
             if (wtsub > 0._r8) q_perch = q_perch/wtsub

             qflx_drain_perched(c) = q_perch_max * q_perch &
                  *(frost_table(c) - zwt(c))

             ! remove drainage from perched saturated layers
             rsub_top_tot = -  qflx_drain_perched(c) * dtime
             do k = jwt(c)+1, k_frz
                rsub_top_layer=max(rsub_top_tot,-(h2osoi_liq(c,k)-watmin))
                rsub_top_layer=min(rsub_top_layer,0._r8)
                rsub_top_tot = rsub_top_tot - rsub_top_layer

                ! save the flux for VSFM
                idx = (c-bounds%begc)*nlevgrnd + k
                mflx_drain_perched_col_1d(idx) = rsub_top_layer/dtime

                if (rsub_top_tot >= 0.) then
                   zwt(c) = zwt(c) - rsub_top_layer/eff_porosity(c,k)/1000._r8
                   exit
                else
                   zwt(c) = zi(c,k)
                endif
             enddo

             ! if rsub_top_tot is greater than available water (above frost table),
             !     then decrease qflx_drain_perched by residual amount for water balance
             qflx_drain_perched(c) = qflx_drain_perched(c) + rsub_top_tot/dtime

             !-- recompute jwt  ---------------------------------------------------------
             ! allow jwt to equal zero when zwt is in top layer
             jwt(c) = nlevsoi
             do j = 1,nlevsoi
                if(zwt(c) <= zi(c,j)) then
                   jwt(c) = j-1
                   exit
                end if
             enddo
          else
             !===================  water table below frost table  =============================
             !--  compute possible perched water table *and* groundwater table afterwards
             ! locate perched water table from bottom up starting at frost table
             ! sat_lev is an arbitrary saturation level used to determine perched water table
             sat_lev=0.9

             k_perch=1
             do k=k_frz,1,-1
                h2osoi_vol = h2osoi_liq(c,k)/(dz(c,k)*denh2o) &
                     + h2osoi_ice(c,k)/(dz(c,k)*denice)

                if (h2osoi_vol/watsat(c,k) <= sat_lev) then
                   k_perch=k
                   exit
                endif
             enddo

             ! if frost_table = nlevsoi, only compute perched water table if frozen
             if (t_soisno(c,k_frz) > tfrz) k_perch=k_frz

             ! if perched water table exists
             if (k_frz > k_perch) then

                ! interpolate between k_perch and k_perch+1 to find perched water table height
                s1 = (h2osoi_liq(c,k_perch)/(dz(c,k_perch)*denh2o) &
                     + h2osoi_ice(c,k_perch)/(dz(c,k_perch)*denice))/watsat(c,k_perch)
                s2 = (h2osoi_liq(c,k_perch+1)/(dz(c,k_perch+1)*denh2o) &
                     + h2osoi_ice(c,k_perch+1)/(dz(c,k_perch+1)*denice))/watsat(c,k_perch+1)

                m=(z(c,k_perch+1)-z(c,k_perch))/(s2-s1)
                b=z(c,k_perch+1)-m*s2
                zwt_perched(c)=max(0._r8,m*sat_lev+b)

                ! compute drainage from perched saturated region
                wtsub = 0._r8
                q_perch = 0._r8
                do k = k_perch, k_frz
                   imped=10._r8**(-e_ice*(0.5_r8*(icefrac(c,k)+icefrac(c,min(nlevsoi, k+1)))))
                   q_perch = q_perch + imped*hksat(c,k)*dzmm(c,k)
                   wtsub = wtsub + dzmm(c,k)
                end do
                if (wtsub > 0._r8) q_perch = q_perch/wtsub

                qflx_drain_perched(c) = q_perch_max * q_perch &
                     *(frost_table(c) - zwt_perched(c))

                ! no perched water table drainage if using original formulation
                if(origflag == 1) qflx_drain_perched(c) = 0._r8

                ! remove drainage from perched saturated layers
                rsub_top_tot = -  qflx_drain_perched(c) * dtime
                do k = k_perch+1, k_frz
                   rsub_top_layer=max(rsub_top_tot,-(h2osoi_liq(c,k)-watmin))
                   rsub_top_layer=min(rsub_top_layer,0._r8)
                   rsub_top_tot = rsub_top_tot - rsub_top_layer

                   ! save the flux for VSFM
                   idx = (c-bounds%begc)*nlevgrnd + k
                   mflx_drain_perched_col_1d(idx) = rsub_top_layer/dtime

                   if (rsub_top_tot >= 0.) then
                      zwt_perched(c) = zwt_perched(c) - rsub_top_layer/eff_porosity(c,k)/1000._r8
                      exit
                   else
                      if (k < k_frz) then
                         s1 = (h2osoi_liq(c,k)/(dz(c,k)*denh2o) &
                              + h2osoi_ice(c,k)/(dz(c,k)*denice))/watsat(c,k)
                         s2 = (h2osoi_liq(c,k+1)/(dz(c,k+1)*denh2o) &
                              + h2osoi_ice(c,k+1)/(dz(c,k+1)*denice))/watsat(c,k+1)
                         if (s1 > sat_lev) then
                            zwt_perched(c) = 0._r8
                         else if (s2 < sat_lev) then
                            zwt_perched(c) = frost_table(c)
                         else if (abs(s2-s1) > 1.e-12_r8) then
                            m=(z(c,k+1)-z(c,k))/(s2-s1)
                            b=z(c,k+1)-m*s2
                            zwt_perched(c)=min(frost_table(c),max(0._r8,m*sat_lev+b))
                         else
                            zwt_perched(c) = min(frost_table(c), zi(c,k))
                         endif
                      else
                         zwt_perched(c) = frost_table(c)
                      endif
                   endif

                enddo

                ! if rsub_top_tot is greater than available water (above frost table),
                !     then decrease qflx_drain_perched by residual amount for water balance
                qflx_drain_perched(c) = qflx_drain_perched(c) + rsub_top_tot/dtime

             else
                qflx_drain_perched(c) = 0._r8
             endif !k_frz > k_perch

             !-- Topographic runoff  ----------------------------------------------------------------------
             fff(c)         = 1._r8/ hkdepth(c)
             dzsum = 0._r8
             icefracsum = 0._r8
             do j = max(jwt(c),1), nlevgrnd
                dzsum  = dzsum + dzmm(c,j)
                icefracsum = icefracsum + icefrac(c,j) * dzmm(c,j)
             end do
             ! add ice impedance factor to baseflow
             if(origflag == 1) then
                if (use_vichydro) then
#ifndef _OPENACC
                   call endrun(msg="VICHYDRO is not available for origflag=1"//errmsg(__FILE__, __LINE__))
#endif
                else
                   fracice_rsub(c) = max(0._r8,exp(-3._r8*(1._r8-(icefracsum/dzsum))) &
                        - exp(-3._r8))/(1.0_r8-exp(-3._r8)) 
                   imped=(1._r8 - fracice_rsub(c))
                   rsub_top_max = 5.5e-3_r8
                end if
             else
                if (use_vichydro) then
                   imped=10._r8**(-e_ice*min(1.0_r8,ice(c,nlayer)/max_moist(c,nlayer)))
                   dsmax_tmp(c) = Dsmax(c) * dtime/ secspday !mm/day->mm/dtime
                   rsub_top_max = dsmax_tmp(c)
                else
                   imped=10._r8**(-e_ice*(icefracsum/dzsum))
                   rsub_top_max = min(10._r8 * sin((rpi/180.) * col_pp%topo_slope(c)), rsub_top_globalmax)
                end if
             endif
             if (use_vichydro) then
                ! ARNO model for the bottom soil layer (based on bottom soil layer
                ! moisture from previous time step
                ! use watmin instead for resid_moist to be consistent with default hydrology
                rel_moist = (moist(c,nlayer) - watmin)/(max_moist(c,nlayer)-watmin)
                frac = (Ds(c) * rsub_top_max )/Wsvic(c)
                rsub_tmp = (frac * rel_moist)/dtime
                if(rel_moist > Wsvic(c)) then
                   frac = (rel_moist - Wsvic(c))/(1.0_r8 - Wsvic(c))
                   rsub_tmp = rsub_tmp + (rsub_top_max * (1.0_r8 - Ds(c)/Wsvic(c)) *frac**c_param(c))/dtime
                end if
                rsub_top(c) = imped * rsub_tmp
                ! make sure baseflow isn't negative
                rsub_top(c) = max(0._r8, rsub_top(c))
             else
                rsub_top(c)    = imped * rsub_top_max* exp(-fff(c)*zwt(c))
             end if
          endif

       end do

       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)

          ! Sub-surface runoff and drainage

          qflx_drain(c) = qflx_rsub_sat(c) + rsub_top(c)

          ! Set imbalance for snow capping

          qflx_qrgwl(c) = qflx_qrgwl(c) + qflx_snwcp_liq(c)

       end do

       ! No drainage for urban columns (except for pervious road as computed above)

       do fc = 1, num_urbanc
          c = filter_urbanc(fc)
          if (col_pp%itype(c) /= icol_road_perv) then
             qflx_drain(c) = 0._r8
             ! This must be done for roofs and impervious road (walls will be zero)
             qflx_qrgwl(c) = qflx_snwcp_liq(c)
          end if
       end do

     end associate

   end subroutine DrainageVSFM

  !-----------------------------------------------------------------------
  subroutine ELMVICMap(bounds, numf, filter, &
       soilhydrology_vars)
     !
     ! !DESCRIPTION:
     ! Performs  the mapping from CLM layers to VIC layers
     ! Specifically, 10 (or 23 when more_vertlayers == .true.)
     ! CLM hydrologically active soil layers are mapped to three VIC layers
     ! by assigning the first nlvic(1) layers to VIC layer 1
     !              the next nlvic(2) layers  to VIC alyer 2
     !              and the remaining to VIC layer 3
     ! mapping from VIC to CLM layers, M.Huang
     !
     ! !USES:
      !$acc routine seq
     use elm_varcon  , only : denh2o, denice, watmin
     use elm_varpar  , only : nlevsoi, nlayer, nlayert, nlevgrnd
     use decompMod   , only : bounds_type
     !
     ! !REVISION HISTORY:
     ! Created by Maoyi Huang
     ! 11/13/2012, Maoyi Huang: rewrite the mapping modules in CLM4VIC
     !
     ! !ARGUMENTS:
     type(bounds_type)        , intent(in)    :: bounds
     integer                  , intent(in)    :: numf      ! number of column soil points in column filter
     integer                  , intent(in)    :: filter(:) ! column filter for soil points
     type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
     !
     ! !LOCAL VARIABLES
     real(r8) :: ice0(1:nlayer)            ! last step ice lens (mm)  (new)
     real(r8) :: moist0(1:nlayer)          ! last step soil water (mm)  (new)
     integer  :: i, j, c, fc
     ! note: in CLM4 h2osoil_liq unit is kg/m2, in VIC moist is mm
     ! h2osoi_ice is actually water equivalent ice content.
     !-----------------------------------------------------------------------

     associate(                                                   &
          dz            => col_pp%dz                               , & ! Input:  [real(r8) (:,:)   ] layer depth (m)
          zi            => col_pp%zi                               , & ! Input:  [real(r8) (:,:)   ] interface level below a "z" level (m)
          z             => col_pp%z                                , & ! Input:  [real(r8) (:,:)   ] layer thickness (m)

          h2osoi_liq    => col_ws%h2osoi_liq       , & ! Input:  [real(r8) (:,:)   ] liquid water (kg/m2)
          h2osoi_ice    => col_ws%h2osoi_ice       , & ! Input:  [real(r8) (:,:)   ] ice lens (kg/m2)
          h2osoi_vol    => col_ws%h2osoi_vol       , & ! Input:  [real(r8) (:,:)   ] volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)

          depth         => soilhydrology_vars%depth_col         , & ! Input:  [real(r8) (:,:)   ] layer depth of upper layer (m)
          porosity      => soilhydrology_vars%porosity_col      , & ! Input:  [real(r8) (:,:)   ] soil porisity (1-bulk_density/soil_density)
          max_moist     => soilhydrology_vars%max_moist_col     , & ! Input:  [real(r8) (:,:)   ] max layer moist + ice (mm)
          vic_elm_fract => soilhydrology_vars%vic_elm_fract_col , & ! Input:  [real(r8) (:,:,:) ] fraction of VIC layers in each CLM layer
          moist         => soilhydrology_vars%moist_col         , & ! Output: [real(r8) (:,:)   ] liquid water (mm)
          ice           => soilhydrology_vars%ice_col           , & ! Output: [real(r8) (:,:)   ] ice lens (mm)
          moist_vol     => soilhydrology_vars%moist_vol_col       & ! Output: [real(r8) (:,:)   ] volumetric soil moisture for VIC soil layers
          )

       ! map CLM to VIC
       do fc = 1, numf
          c = filter(fc)
          do i = 1, nlayer
             ice0(i) = ice(c,i)
             moist0(i) = moist(c,i)
             ice(c,i) = 0._r8
             moist(c,i) = 0._r8
             do j = 1, nlevsoi
                ice(c,i) = ice(c,i) + h2osoi_ice(c,j) * vic_elm_fract(c,i,j)
                moist(c,i) = moist(c,i) + h2osoi_liq(c,j) * vic_elm_fract(c,i,j)
             end do
             ice(c,i)       = min((moist0(i) + ice0(i)), ice(c,i))
             ice(c,i)       = max(0._r8, ice(c,i))
             moist(c,i)     = max(watmin, moist(c,i))
             moist(c,i)     = min(max_moist(c,i)-ice(c,i), moist(c,i))
             moist_vol(c,i) = moist(c,i)/(depth(c,i)*denice) + ice(c,i)/(depth(c,i)*denh2o)
             moist_vol(c,i) = min(porosity(c,i), moist_vol(c,i))
             moist_vol(c,i) = max(0.01_r8, moist_vol(c,i))
          end do

          ! hydrologic inactive layers
          ice(c, nlayer+1:nlayert)       = h2osoi_ice(c, nlevsoi+1:nlevgrnd)
          moist(c, nlayer+1:nlayert)     = h2osoi_liq(c, nlevsoi+1:nlevgrnd)
          moist_vol(c, nlayer+1:nlayert) = h2osoi_vol(c, nlevsoi+1:nlevgrnd)
       end do

     end associate

   end subroutine ELMVICMap

end module SoilHydrologyMod
