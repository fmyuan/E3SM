#ifdef USE_ATS_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:

  ! The ELM Interface to ATS (Fortran interface),
  ! corresponding to 'ATS elm_ats_api.h' (extern "C" interface)
  !
  ! 2020-04-26: hydrological process
  !
  interface
    !-----------------------------
    function ats_create_f(comm, input_filename) bind(c, name="ats_create_c")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      integer(C_INT), intent(in) :: comm
      character(C_CHAR), intent(in) :: input_filename(*)
      type(C_PTR) :: ats_create_f
    end function ats_create_f

    subroutine ats_delete_f(ats) bind(c, name="ats_delete_c")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR) :: ats
    end subroutine ats_delete_f

    ! the following is for ats-mesh passing into ELM.
    subroutine ats_get_mesh_info_f(ats, ncols_local, ncols_global, lat, lon, elev, surf_area, nlevgrnd, depth) &
        bind(c, name="ats_get_mesh_info_c")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      integer(C_INT), intent(out) :: ncols_local        ! number of columns on this rank
      integer(C_INT), intent(out) :: ncols_global       ! number of total columns
      real(C_DOUBLE), dimension(*), intent(out) :: lat  ! latitude of each grid cell, units: decimal degree [col]
      real(C_DOUBLE), dimension(*), intent(out) :: lon  ! longitude of each grid cell, units: decimal degree [col]
      real(C_DOUBLE), dimension(*), intent(out) :: elev ! elevation of each surface grid cell, units: m [col]
      real(C_DOUBLE), dimension(*), intent(out) :: surf_area ! surface area of each surface grid cell, units: m^2 [col]

      ! NOTE: these parameters are a check to confirm that ATS and ELM
      ! have the same mesh structure, and cannot be changed in either
      ! ATS or ELM.
      integer(C_INT), intent(out) :: nlevgrnd             ! number of cells in the vertical column
      real(C_DOUBLE), dimension(*), intent(out) :: depth  ! depth from the surface of each cell, units: m [nlevgrnd]
    end subroutine ats_get_mesh_info_f
    
    subroutine ats_setup_f(ats) bind(c, name="ats_setup_c")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
    end subroutine ats_setup_f
    
    subroutine ats_initialize_f(ats, start_t, patm, pres) bind(c, name="ats_initialize_c")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), intent(in) :: start_t                       ! unit: second
      real(C_DOUBLE), dimension(*), intent(in) :: patm     ! atm pressure,                 unit: Pa, [col]
      real(C_DOUBLE), dimension(*), intent(in) :: pres     ! water pressure,               unit: -, [col,nlevgrnd]
    end subroutine ats_initialize_f

    subroutine ats_advance_f(ats, dt, checkpoint, visualout) &
        bind(c, name="ats_advance_c")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), intent(in) :: dt                       ! unit: seconds
      logical(C_BOOL), intent(in) :: checkpoint              ! write ATS checkpoint file after advancing
      logical(C_BOOL), intent(in) :: visualout               ! write ATS vis file after advancing
    end subroutine ats_advance_f

    !----------------------------

    subroutine ats_set_soil_hydrologic_parameters_f(ats, &
      base_porosity, hksat, CH_bsw, CH_smpsat, CH_sr) &
           bind(c, name="ats_set_soil_hydrologic_parameters_c")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(in) :: base_porosity! soil porosity, not including compressibility or ice, unit: -, [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(in) :: hksat        ! sat. hydraulic conductivity, unit: mm/s, [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(in) :: CH_bsw       ! Clapp-Hornberger "b",        unit: -, [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(in) :: CH_smpsat    ! Clapp-Hornberger "smpsat",   unit: Pa,[col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(in) :: CH_sr        ! Clapp-Hornberger res. sat. if any, unit: -, [col,nlevgrnd]
    end subroutine ats_set_soil_hydrologic_parameters_f


    subroutine ats_set_veg_parameters_f(ats, &
      mafic_pot_full_turgor, mafic_pot_wilt_point) &
           bind(c, name="ats_set_veg_parameters_c")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats

      real(C_DOUBLE), dimension(*), intent(in) :: mafic_pot_full_turgor ! mafic potential at full turgor, unit: [Pa], [pft]
      real(C_DOUBLE), dimension(*), intent(in) :: mafic_pot_wilt_point  ! mafic potential at the wilting point, unit: [Pa], [pft]
    end subroutine ats_set_veg_parameters_f
      
    !----------------------------

    subroutine ats_set_veg_properties_f(ats, &
         rooting_fraction) &
           bind(c, name="ats_set_veg_properties_c")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats

      real(C_DOUBLE), dimension(*), intent(in) :: rooting_fraction   ! unit: [-], [col, nlevgrnd, pft]
    end subroutine ats_set_veg_properties_f
    
    subroutine ats_set_potential_sources_f(ats, ss_surface_sources, ss_soilevap, ss_vegtran) &
         bind(c, name="ats_set_potential_sources_c")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(in) :: ss_surface_sources ! water source/sink into top of the column, unit: m/s, [col]
      real(C_DOUBLE), dimension(*), intent(in) :: ss_soilevap        ! potential evaporation into the top of the column (positive = condensation)
                                                                     !   unit: m/s, [col]
      real(C_DOUBLE), dimension(*), intent(in) :: ss_vegtran         ! potential water source/sink of transpiration, (negative = water loss)
                                                                     !   unit: m/s, [col, pft]
    end subroutine ats_set_potential_sources_f

    subroutine ats_set_soil_hydrologic_properties_f(ats, eff_porosity) &
         bind(c, name="ats_set_soil_hydrologic_properties_c")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(in) :: eff_porosity ! soil porosity, with compressibility & excluding ice, unit: -, [col,nlevgrnd]
    end subroutine ats_set_soil_hydrologic_properties_f
    
    
    !----------------------------
    subroutine ats_get_waterstate_f(ats, surf_pd, soil_pres, soil_psi, h2oliq, h2oice) &
         bind(c, name="ats_get_waterstate_c")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats

      real(C_DOUBLE), dimension(*), intent(out) :: surf_pd           ! surface water ponding depth, unit: [m] [col]
      real(C_DOUBLE), dimension(*), intent(out) :: soil_pres         ! soil water pressure, unit: [Pa] [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(out) :: soil_psi          ! soil matric potential, unit: [+Pa] [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(out) :: h2oliq            ! soil liq. water content, unit: [kgH2O / m3] [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(out) :: h2oice            ! soil ice water content, unit: [kgH2O / m3] [col,nlevgrnd]
    end subroutine ats_get_waterstate_f

    subroutine ats_get_water_fluxes_f(ats, soil_infl, evap, vegtran, net_subsurface_flux, net_runon) &
         bind(c, name="ats_get_water_fluxes_c")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(out) :: soil_infl   ! water source/sink into top soil - gross infiltration, unit: m/s, [col]
      real(C_DOUBLE), dimension(*), intent(out) :: evap        ! actual evaporation into top of col (positive = condensation), unit: m/s, [col]
      real(C_DOUBLE), dimension(*), intent(out) :: vegtran     ! actual transpiration, (negative = water loss), unit: m/s, [col, nlevgrnd, pft]
      real(C_DOUBLE), dimension(*), intent(out) :: net_subsurface_flux ! net flux of water across all column boundaries (positive = inward), unit: m/s, [col]
      real(C_DOUBLE), dimension(*), intent(out) :: net_runon   ! net flux of water across all surface boundaries (positive = inward), unit: m/s, [col]
    end subroutine ats_get_water_fluxes_f

  end interface

  ! c++-fortran interface
  !------------------------------------------------------------------------


#endif

