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

    subroutine ats_advance_f(ats, dt) &
        bind(c, name="ats_advance_c")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), intent(in) :: dt                       ! unit: seconds
    end subroutine ats_advance_f

    !----------------------------

    subroutine ats_set_soil_hydrologic_properties_f(ats, &
      porosity, hksat, CH_bsw, CH_smpsat, CH_sr) &
           bind(c, name="ats_set_soil_hydrologic_properties_c")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(in) :: porosity     ! soil porosity,               unit: -, [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(in) :: hksat        ! sat. hydraulic conductivity, unit: -, [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(in) :: CH_bsw       ! Clapp-Hornberger "b",        unit: -, [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(in) :: CH_smpsat    ! Clapp-Hornberger "smpsat",   unit: Pa,[col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(in) :: CH_sr        ! Clapp-Hornberger res. sat. if any, unit: -, [col,nlevgrnd]
    end subroutine ats_set_soil_hydrologic_properties_f

    subroutine ats_set_potential_sources_f(ats, ss_soilinfl, ss_soilevap, ss_vegtran) &
         bind(c, name="ats_set_potential_sources_c")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(in) :: ss_soilinfl        ! water source/sink into/on top soil - gross infiltration, unit: m/s, [col]
      real(C_DOUBLE), dimension(*), intent(in) :: ss_soilevap        ! water source/sink into/on top soil - potential evaporation, unit: m/s, [col]
      real(C_DOUBLE), dimension(*), intent(in) :: ss_vegtran         ! water source/sink of total potential vegetation transpiration, unit: m/s, [col, pft]
    end subroutine ats_set_potential_sources_f

    subroutine ats_get_actual_sources_f(ats, ss_soilinfl, ss_soilevap, ss_vegtran) &
         bind(c, name="ats_get_actual_sources_c")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(out) :: ss_soilinfl        ! water source/sink into/on top soil - gross infiltration, unit: m/s, [col]
      real(C_DOUBLE), dimension(*), intent(out) :: ss_soilevap        ! water source/sink into/on top soil - potential evaporation, unit: m/s, [col]
      real(C_DOUBLE), dimension(*), intent(out) :: ss_vegtran         ! water source/sink of total potential vegetation transpiration, unit: m/s, [col, pft]
    end subroutine ats_get_actual_sources_f
    
    !----------------------------
    subroutine ats_get_waterstate_f(ats, surf_pd, soil_pres, soil_psi, h2oliq, h2oice) &
         bind(c, name="ats_get_waterstate_c")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(out) :: surf_pd           ! surface water ponding depth (m)
      real(C_DOUBLE), dimension(*), intent(out) :: soil_pres         ! soil water pressure (Pa)
      real(C_DOUBLE), dimension(*), intent(out) :: soil_psi          ! soil matric potential (+Pa)
      real(C_DOUBLE), dimension(*), intent(out) :: h2oliq            ! soil liq. water content [kgH2O/m3]
      real(C_DOUBLE), dimension(*), intent(out) :: h2oice            ! soil ice water content [kgH2O/m3]
    end subroutine ats_get_waterstate_f

    ! the following is for ats-mesh passing into ELM.
    subroutine ats_get_mesh_info_f(ats, ncols_local, ncols_global, ncells_per_col, lat, lon, elev, surf_area, dz, depth) &
        bind(c, name="ats_get_mesh_info_c")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      integer(C_INT), intent(out) :: ncols_local
      integer(C_INT), intent(out) :: ncols_global
      integer(C_INT), intent(out) :: ncells_per_col
      real(C_DOUBLE), dimension(*), intent(out) :: lat
      real(C_DOUBLE), dimension(*), intent(out) :: lon
      real(C_DOUBLE), dimension(*), intent(out) :: elev
      real(C_DOUBLE), dimension(*), intent(out) :: surf_area
      real(C_DOUBLE), dimension(*), intent(out) :: dz
      real(C_DOUBLE), dimension(*), intent(out) :: depth
    end subroutine ats_get_mesh_info_f

  end interface

  ! c++-fortran interface
  !------------------------------------------------------------------------


#endif

