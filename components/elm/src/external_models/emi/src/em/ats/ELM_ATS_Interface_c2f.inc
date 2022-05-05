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
    function ats_elm_create() bind(c, name="ats_create")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR) :: ats_elm_create
    end function ats_elm_create

    subroutine ats_elm_delete(ats) bind(c, name="ats_delete")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR) :: ats
    end subroutine ats_elm_delete

    subroutine ats_elm_setup(ats, comm, input_filename) bind(c, name="ats_setup")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      integer(C_INT), intent(in) :: comm
      character(C_CHAR), intent(in) :: input_filename(*)
    end subroutine ats_elm_setup

    subroutine ats_elm_initialize(ats) bind(c, name="ats_initialize")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
    end subroutine

    subroutine ats_elm_advance(ats, start_t, dt, visout, chkout) &
        bind(c, name="ats_advance")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), intent(in) :: start_t          ! unit: second
      real(C_DOUBLE), intent(in) :: dt               ! unit: seconds
      logical(C_BOOL),intent(in), value :: visout            ! instruct ATS output data
      logical(C_BOOL),intent(in), value :: chkout            ! instruct ATS output checkpoint
    end subroutine

    !----------------------------

    subroutine ats_elm_setmesh(ats, gridsX, gridsY, gridsZ, soilcol_nodes, lx, ly, lz) bind(c, name="ats_set_mesh")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(in) :: gridsX         ! unit: meters (grid nodes, not mid-cells)
      real(C_DOUBLE), dimension(*), intent(in) :: gridsY         ! unit: meters
      real(C_DOUBLE), dimension(*), intent(in) :: gridsZ         ! unit: meters
      real(C_DOUBLE), dimension(*), intent(in) :: soilcol_nodes  ! unit: meters
      integer(C_INT), value, intent(in) :: lx                    ! length of gridsX
      integer(C_INT), value, intent(in) :: ly                    ! length of gridsY
      integer(C_INT), value, intent(in) :: lz                    ! length of soilcol_nodes
    end subroutine ats_elm_setmesh

    subroutine ats_elm_setmats(ats, porosity, hksat, CH_bsw, CH_smpsat, CH_sr, &
             eff_porosity) bind(c, name="ats_set_materials")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(in) :: porosity     ! soil porosity,               unit: -, [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(in) :: hksat        ! sat. hydraulic conductivity, unit: -, [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(in) :: CH_bsw       ! Clapp-Hornberger "b",        unit: -, [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(in) :: CH_smpsat    ! Clapp-Hornberger "smpsat",   unit: Pa,[col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(in) :: CH_sr        ! Clapp-Hornberger res. sat. if any, unit: -, [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(in) :: eff_porosity ! effective soil porosity (adjusted, e.g. by ice etc), unit: -, [col,nlevgrnd]
    end subroutine ats_elm_setmats

    subroutine ats_elm_setIC(ats, patm, soilpressure, wtd) bind(c, name="ats_set_initialconditions")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(in) :: patm            ! atm. air pressure, unit: Pa, [col]
      real(C_DOUBLE), dimension(*), intent(in) :: soilpressure    ! soil hydraulic pressure, unit: Pa, [col,nlevgrnd]
      real(C_DOUBLE), dimension(*), intent(in) :: wtd             ! soil water table depth, unit: m, [col]
    end subroutine ats_elm_setIC

    subroutine ats_elm_setBC(ats) bind(c, name="ats_set_boundaryconditions")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
    end subroutine ats_elm_setBC

    subroutine ats_elm_setSS(ats, ss_soilinfl, ss_soilevap, ss_roottran, ncols, ncells) bind(c, name="ats_set_sources")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(in) :: ss_soilinfl        ! water source/sink into/on top soil - gross infiltration, unit: m/s, [col]
      real(C_DOUBLE), dimension(*), intent(in) :: ss_soilevap        ! water source/sink into/on top soil - potential evaporation, unit: m/s, [col]
      real(C_DOUBLE), dimension(*), intent(in) :: ss_roottran        ! water root extranction potential, unit: m/s, [col,nlevgrnd]
      integer(C_INT), intent(in) :: ncols                            ! col (for data checking)
      integer(C_INT), intent(in) :: ncells                           ! col*nlevgrnd (for data checking)
    end subroutine ats_elm_setSS

    !----------------------------
    subroutine ats_elm_getwaterstate(ats, surf_pres, soil_pres, satliq, ncols, ncells) bind(c, name="ats_get_waterstate")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(out) :: surf_pres
      real(C_DOUBLE), dimension(*), intent(out) :: soil_pres
      real(C_DOUBLE), dimension(*), intent(out) :: satliq
      integer(C_INT), intent(in) :: ncols                            ! col (for data checking)
      integer(C_INT), intent(in) :: ncells                           ! col*nlevgrnd (for data checking)
     end subroutine ats_elm_getwaterstate

    ! the following is for ats-mesh passing into ELM.
    subroutine ats_elm_getmesh(ats, ncols_local, ncols_global, ncells_per_col, dz, depth, elev, surf_area_m2, lat, lon) &
        bind(c, name="ats_get_mesh_info")
      use, intrinsic :: iso_c_binding
      import
      implicit none
      type(C_PTR), value :: ats
      real(C_DOUBLE), dimension(*), intent(out) :: dz
      real(C_DOUBLE), dimension(*), intent(out) :: depth
      real(C_DOUBLE), dimension(*), intent(out) :: elev
      real(C_DOUBLE), dimension(*), intent(out) :: surf_area_m2
      real(C_DOUBLE), dimension(*), intent(out) :: lat
      real(C_DOUBLE), dimension(*), intent(out) :: lon
      integer(C_INT), intent(out) :: ncols_local
      integer(C_INT), intent(out) :: ncols_global
      integer(C_INT), intent(out) :: ncells_per_col
    end subroutine

  end interface

  ! c++-fortran interface
  !------------------------------------------------------------------------


#endif

