#ifdef USE_ATS_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:

  ! Fortran bindings to the ATS library C interface 

  ! This should probably live within ATS
  interface
    !-----------------------------
    function ats_create_c(comm, input_filename) bind(c, name="ats_create")
      use, intrinsic :: iso_c_binding
      implicit none
      integer, intent(in) :: comm
      character(C_CHAR), intent(in) :: input_filename(*)
      type(C_PTR) :: ats_create_c
    end function ats_create_c

    subroutine ats_delete_c(ats) bind(c, name="ats_delete")
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR) :: ats
    end subroutine ats_delete_c

    subroutine ats_get_mesh_info_c(ats, ncols_local, ncols_global, lat, lon, elev, surf_area, pft, nlevgrnd, depth) &
        bind(c, name="ats_get_mesh_info")
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value :: ats
      integer(C_INT), intent(out) :: ncols_local        ! number of columns on this rank
      integer(C_INT), intent(out) :: ncols_global       ! number of total columns
      real(C_DOUBLE), dimension(*), intent(out) :: lat  ! latitude of each grid cell, units: decimal degree [col]
      real(C_DOUBLE), dimension(*), intent(out) :: lon  ! longitude of each grid cell, units: decimal degree [col]
      real(C_DOUBLE), dimension(*), intent(out) :: elev ! elevation of each surface grid cell, units: m [col]
      real(C_DOUBLE), dimension(*), intent(out) :: surf_area ! surface area of each surface grid cell, units: m^2 [col]
      integer(C_INT), dimension(*), intent(out) :: pft ! pft type of each grid cell [col]
      integer(C_INT), intent(out) :: nlevgrnd             ! number of cells in the vertical column
      real(C_DOUBLE), dimension(*), intent(out) :: depth  ! depth from the surface of each cell, units: m [nlevgrnd]
    end subroutine ats_get_mesh_info_c
    
    subroutine ats_setup_c(ats) bind(c, name="ats_setup")
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value :: ats
    end subroutine ats_setup_c
    
    subroutine ats_initialize_c(ats, start_time, soil_water_content, soil_pressure) bind(c, name="ats_initialize")
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), intent(in) :: start_time                         ! unit: second
      real(C_DOUBLE), dimension(*), intent(in) :: soil_water_content   ! unit: kgH2O/m2, [col, nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(in) :: soil_pressure        ! unit: Pa, [col, nlevgrnd]
    end subroutine ats_initialize_c

    subroutine ats_advance_c(ats, dt, checkpoint, visualout) &
        bind(c, name="ats_advance")
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), intent(in) :: dt                       ! unit: seconds
      logical(C_BOOL), intent(in) :: checkpoint              ! write ATS checkpoint file after advancing
      logical(C_BOOL), intent(in) :: visualout               ! write ATS vis file after advancing
    end subroutine ats_advance_c

    !----------------------------

    subroutine ats_set_soil_hydrologic_parameters_c(ats, &
      base_porosity, hydraulic_conductivity, clapp_horn_b, clapp_horn_smpsat, clapp_horn_sr) &
           bind(c, name="ats_set_soil_hydrologic_parameters")
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(in) :: base_porosity          ! soil porosity, not including compressibility or ice, unit: -, [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(in) :: hydraulic_conductivity ! sat. hydraulic conductivity, unit: mm/s, [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(in) :: clapp_horn_b           ! Clapp-Hornberger "b",        unit: -, [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(in) :: clapp_horn_smpsat      ! Clapp-Hornberger "smpsat",   unit: Pa, [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(in) :: clapp_horn_sr          ! Clapp-Hornberger res. sat. if any, unit: -, [col,nlevgrnd]
    end subroutine ats_set_soil_hydrologic_parameters_c
      
    !----------------------------

    subroutine ats_set_veg_properties_c(ats, rooting_fraction) &
           bind(c, name="ats_set_veg_properties")
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(in) :: rooting_fraction   ! unit: [-], [col, nlevgrnd]
    end subroutine ats_set_veg_properties_c
    
    subroutine ats_set_potential_sources_c(ats, surface_source, potential_evaporation, potential_transpiration) &
         bind(c, name="ats_set_sources")
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(in) :: surface_source           ! surface water source (precip, snowmelt, throughfall)  unit: mm/s, [col]
      real(C_DOUBLE), dimension(*), intent(in) :: potential_evaporation    ! potential evaporation (positive = condensation)  unit: mm/s, [col]
      real(C_DOUBLE), dimension(*), intent(in) :: potential_transpiration  ! potential column-summed transpiration (negative = water loss)  unit: mm/s, [col]
    end subroutine ats_set_potential_sources_c

    subroutine ats_set_soil_hydrologic_properties_c(ats, eff_porosity) &
         bind(c, name="ats_set_soil_hydrologic_properties")
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(in) :: eff_porosity ! soil porosity, with compressibility & excluding ice, unit: -, [col,nlevgrnd]
    end subroutine ats_set_soil_hydrologic_properties_c
    
    
    !----------------------------
    subroutine ats_get_waterstate_c(ats, ponded_depth, water_table_depth, mass_water_content, mass_ice_content, soil_psi) &
         bind(c, name="ats_get_waterstate")
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(out) :: ponded_depth       ! surface water ponding depth, unit: [m] [col]
      real(C_DOUBLE), dimension(*), intent(out) :: water_table_depth  ! water table depth, unit: [m] [col]
      real(C_DOUBLE), dimension(*), intent(out) :: mass_water_content ! soil liq. water content, unit: [kgH2O / m2] [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(out) :: mass_ice_content   ! soil ice water content, unit: [kgH2O / m2] [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(out) :: soil_psi           ! soil matric potential, unit: [Pa] [col,nlevgrnd]
    end subroutine ats_get_waterstate_c

    subroutine ats_get_water_fluxes_c(ats, infiltration, evaporation, transpiration, root_fluxes, net_subsurface_flux, net_runon) &
         bind(c, name="ats_get_water_fluxes")
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(out) :: infiltration        ! actaual liq. water source/sink into top soil - net infiltration, unit: mm/s, [col]
      real(C_DOUBLE), dimension(*), intent(out) :: evaporation         ! actual evaporation into top of col (positive = condensation), unit: mm/s, [col]
      real(C_DOUBLE), dimension(*), intent(out) :: transpiration       ! actual transpiration, (negative = water loss), unit: mm/s, [col]
      real(C_DOUBLE), dimension(*), intent(out) :: root_fluxes         ! actual transpiration distributed to root layers, (negative = water loss), unit: mm/s, [col, nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(out) :: net_subsurface_flux ! net flux of water across all column boundaries (positive = inward), unit: mm/s, [col]
      real(C_DOUBLE), dimension(*), intent(out) :: net_runon           ! net flux of water across all surface boundaries (positive = inward), unit: mm/s, [col]
    end subroutine ats_get_water_fluxes_c

  end interface

  ! c++-fortran interface
  !------------------------------------------------------------------------


#endif

