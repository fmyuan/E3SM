set(CXX_LINKER "CXX")
set(NETCDF_PATH "$ENV{NETCDF_PATH}")
set(NETCDF_C_PATH "$ENV{NETCDF_C_PATH}")
set(NETCDF_FORTRAN_PATH "$ENV{NETCDF_FORTRAN_PATH}")
set(PNETCDF_PATH "$ENV{PNETCDF_PATH}")
set(HDF5_PATH "$ENV{HDF5_PATH}")
set(LAPACK_LIBDIR "$ENV{BLASLAPACK_DIR}")
string(APPEND CFLAGS " -mcmodel=small")
string(APPEND CFLAGS " -I${NETCDF_PATH}/include")
string(APPEND FFLAGS " -mcmodel=small -fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none")
string(APPEND FFLAGS " -I${NETCDF_PATH}/include")
string(APPEND LDFLAGS " -framework Accelerate")
if (COMP_CLASS STREQUAL cpl)
  string(APPEND LDFLAGS " -L${LAPACK_LIBDIR} -llapack -lblas")
endif()
string(APPEND LDFLAGS " -L$ENV{CC_ROOT}/lib -lgomp")
execute_process(COMMAND $ENV{NETCDF_PATH}/bin/nc-config --flibs OUTPUT_VARIABLE SHELL_CMD_OUTPUT_BUILD_INTERNAL_IGNORE0 OUTPUT_STRIP_TRAILING_WHITESPACE)
string(APPEND SLIBS " ${SHELL_CMD_OUTPUT_BUILD_INTERNAL_IGNORE0} -lnetcdf")
execute_process(COMMAND $ENV{NETCDF_PATH}/bin/nf-config --flibs OUTPUT_VARIABLE SHELL_CMD_OUTPUT_BUILD_INTERNAL_IGNORE0 OUTPUT_STRIP_TRAILING_WHITESPACE)
string(APPEND SLIBS " ${SHELL_CMD_OUTPUT_BUILD_INTERNAL_IGNORE0} -lnetcdff")
set(SFC "$ENV{FC_ROOT}/bin/gfortran")
set(SCC "$ENV{CC_ROOT}/bin/gcc")
set(SCXX "$ENV{CC_ROOT}/bin/g++")
set(MPICC "$ENV{MPI_ROOT}/bin/mpicc")
set(MPICXX "$ENV{MPI_ROOT}/bin/mpicxx")
set(MPIFC "$ENV{MPI_ROOT}/bin/mpifort")