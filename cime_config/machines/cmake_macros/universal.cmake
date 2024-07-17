# Default guess at KOKKOS_OPTIONS. These can be overridden by
# appending -D${OPTION_NAME}=Off. The CMAKE_CXX_COMPILER will
# be handled automatically by CIME unless explicitly set in
# KOKKOS_OPTIONS. CMAKE_INSTALL_PREFIX will be handled automatically
# by CIME always.
set(KOKKOS_OPTIONS "-DKokkos_ENABLE_SERIAL=On")
if (compile_threaded)
  string(APPEND KOKKOS_OPTIONS " -DKokkos_ENABLE_OPENMP=On")
endif()

# Unless told otherwise, set has_contiguous to FALSE
set(HAS_F2008_CONTIGUOUS "FALSE")

# By default, link with CXX
set(E3SM_LINK_WITH_FORTRAN "FALSE")

# Do not use any automatic settings for these
set(CMAKE_Fortran_FLAGS "")
set(CMAKE_Fortran_FLAGS_DEBUG "")
set(CMAKE_Fortran_FLAGS_RELEASE "")
set(CMAKE_C_FLAGS "")
set(CMAKE_C_FLAGS_DEBUG "")
set(CMAKE_C_FLAGS_RELEASE "")
set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_DEBUG "")
set(CMAKE_CXX_FLAGS_RELEASE "")
set(CMAKE_Fortran_FORMAT_FIXED_FLAG "")
set(CMAKE_Fortran_FORMAT_FREE_FLAG "")

# ATS libraries, previously in userdefined.cmake, but which excluded in cime now.
set(AMANZI_TPLS_DIR "$ENV{AMANZI_TPLS_DIR}")
set(ATS_DIR "$ENV{ATS_DIR}")
if (COMP_NAME STREQUAL elm)
  if (NOT ${AMANZI_TPLS_DIR} STREQUAL "")

    string(APPEND CMAKE_Fortran_FLAGS " -I${AMANZI_TPLS_DIR}/trilinos-15-1-0/include")
    string(APPEND CMAKE_Fortran_FLAGS " -I${AMANZI_TPLS_DIR}/SEACAS/include ")
    string(APPEND CMAKE_Fortran_FLAGS " -I${AMANZI_TPLS_DIR}/petsc-3.20/include -I${AMANZI_TPLS_DIR}/pflotran/src ")
    if (NOT ${ATS_DIR} STREQUAL "")
      string(APPEND CPPDEFS " -DUSE_ATS_LIB ")
      string(APPEND CMAKE_Fortran_FLAGS " -I${ATS_DIR}/include ")
    endif()
  endif()
endif()
if (COMP_NAME STREQUAL cpl)
  string(APPEND CMAKE_EXE_LINKER_FLAGS " -lstdc++")
  if (NOT ${AMANZI_TPLS_DIR} STREQUAL "")
    string(APPEND CMAKE_EXE_LINKER_FLAGS " -L${AMANZI_TPLS_DIR}/lib")
    string(APPEND CMAKE_EXE_LINKER_FLAGS " -L${AMANZI_TPLS_DIR}/trilinos-15-1-0/lib")
    string(APPEND CMAKE_EXE_LINKER_FLAGS " -L${AMANZI_TPLS_DIR}/SEACAS/lib ")
    string(APPEND CMAKE_EXE_LINKER_FLAGS " -L${AMANZI_TPLS_DIR}/petsc-3.20/lib -L${AMANZI_TPLS_DIR}/pflotran/src ")
    if (NOT ${ATS_DIR} STREQUAL "")
      string(APPEND CMAKE_EXE_LINKER_FLAGS " -L${ATS_DIR}/lib -lerror_handling -latk -lfunctions -lgeometry -lgeochemutil -lgeochemsolvers -lgeochembase -lgeochemrxns -lgeochemistry -lmesh -lmesh_simple -lmesh_mstk -lmesh_extracted -lmesh_logical -lmesh_factory -ldbg -lwhetstone -ldata_structures -lmesh_functions -loutput -lstate -lsolvers -ltime_integration -loperators -lpks -lchemistry_pk -ltransport -lshallow_water -lats_operators -lats_eos -lats_surf_subsurf -lats_generic_evals -lats_column_integrator -lats_pks -lats_energy_relations -lats_energy -lats_flow_relations -lats_flow -lats_transport -lats_deform -lats_surface_balance -lats_bgc -lats_mpc_relations -lats_mpc -lats_executable -lelm_ats")
    endif()
  endif()
endif()


string(APPEND CPPDEFS " -DCPL_BYPASS")
