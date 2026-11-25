#!/usr/bin/env bash

# exit on error
set -e

# set up local variables
export RUN_NAME="tempest_serial"
export INPUTDATA_NAME="110x1pt_US-GC_TransTEMPEST"
export COMPSET="ICB20TRCNPRDCTCBC"
export CASE_DIR="${E3SM_WORK_DIR}/output/cases/${RUN_NAME}"
export GITHUB_ACTIONS=TRUE
export E3SM_SRC_DIR="${ELM_ATS_SRC_DIR}/E3SM"

# create the case
echo "Creating case"
echo "----------------------"
echo "${E3SM_SRC_DIR}/cime/scripts/create_newcase --case ${CASE_DIR} --res ELM_USRDAT --mach ${MACHINE_NAME} --compiler ${COMPILER_NAME} --compset ${COMPSET}"
echo "----------------------"
${E3SM_SRC_DIR}/cime/scripts/create_newcase --case ${CASE_DIR} --res ELM_USRDAT --mach ${MACHINE_NAME} --compiler ${COMPILER_NAME} --compset ${COMPSET}

# cp over example files to the case directory
echo ""
echo "Setting up case input"
echo "----------------------"
if [ $GITHUB_ACTIONS ]; then
  cp $(find / -path '*/examples/tempest_transect')/* ${CASE_DIR}
else
  cp ./* ${CASE_DIR}
fi
cd ${CASE_DIR}
sed -i "s^MESH_FILENAME^${CASE_DIR}/${RUN_NAME}.exo^g" ${CASE_DIR}/${RUN_NAME}.xml

./xmlchange MOSART_MODE=NULL,DOUT_S=FALSE,DIN_LOC_ROOT=${E3SM_WORK_DIR}/inputdata
./xmlchange DIN_LOC_ROOT_CLMFORC=\$DIN_LOC_ROOT/atm/datm7
./xmlchange ELM_USRDAT_NAME=${INPUTDATA_NAME}

./xmlchange ATM_DOMAIN_PATH=\$DIN_LOC_ROOT/share/domains/domain.clm
./xmlchange LND_DOMAIN_PATH=\$DIN_LOC_ROOT/share/domains/domain.clm

# these are now example specific
./xmlchange ATM_DOMAIN_FILE=domain.lnd.${INPUTDATA_NAME}_c20230901.nc
./xmlchange LND_DOMAIN_FILE=domain.lnd.${INPUTDATA_NAME}_c20230901.nc

./xmlchange NTASKS=1
./xmlchange NTASKS_PER_INST=1
./xmlchange PIO_TYPENAME=netcdf
./xmlchange RUN_STARTDATE=2000-07-15
./xmlchange STOP_OPTION=nmonths,STOP_N=2
./xmlchange BATCH_SYSTEM=none
./xmlchange DEBUG=TRUE

# set up the user_nl_elm file prior to setup
# generic part!
echo " ats_inputdir = '${CASE_DIR}'" >> user_nl_elm
echo " ats_inputfile = '${RUN_NAME}.xml'" >> user_nl_elm

cat user_nl_elm

# setup the case
echo ""
echo "Running case.setup"
echo "----------------------"
./case.setup
echo -e '\nstring(APPEND CPPDEFS " -DCPL_BYPASS")' >> cmake_macros/universal.cmake

# build
echo ""
echo "Running case.build"
echo "----------------------"
./case.build


echo ""
echo "Run the case yourself:"
echo "----------------------"
echo "cd ${CASE_DIR} && ./case.submit"
if [ $GITHUB_ACTIONS ]; then
  cd ${CASE_DIR} && ./case.submit --no-batch 
fi 
