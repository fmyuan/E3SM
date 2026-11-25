#!/usr/bin/env bash

set -e
shopt -s nullglob

if [ -z "${USE_ATS}" ]; then
    echo "Set USE_ATS prior to running the build_example script"
    exit 1
fi

if [ ! -v NTASKS ]; then
    echo "Setting NTASKS = 1"
    NTASKS=1
fi

# ----
# process directory structure
if [ -z "${ELM_ATS_SRC_DIR}" ]; then
    echo "Set ELM_ATS_SRC_DIR before running."
    exit 1
fi
if [ -z "${E3SM_WORK_DIR}" ]; then
    echo "Set E3SM_WORK_DIR before running."
    exit 1
fi

NAME=coweeta
export CASE_NAME=coweeta.np${NTASKS}
export DOMAIN_FILE=coweeta_domain.nc
export INPUTDATA_FILE=coweeta_surfdata.nc
export COMPSET="ICB20TRCNPRDCTCBC"


# use ats?
CASE_SUFFIX=
if [ -z "${USE_ATS}" ]; then
    echo "Set USE_ATS before running."
    exit 1
elif [ "${USE_ATS}" == "TRUE" ]; then
    CASE_SUFFIX=elm-ats
elif [ "${USE_ATS}" == "IC_ONLY" ]; then
    CASE_SUFFIX=ic_only
elif [ "${USE_ATS}" == "FALSE" ]; then
    CASE_SUFFIX=elm
else
    echo "USE_ATS must be one of TRUE, IC_ONLY, or FALSE"
fi



CASE_DIR="${E3SM_WORK_DIR}/cases/${CASE_NAME}.${CASE_SUFFIX}"
E3SM_SRC_DIR="${ELM_ATS_SRC_DIR}/E3SM"
# Use GNU sed if available (gsed on macOS), otherwise plain sed
if command -v gsed &> /dev/null; then
    SED=gsed
else
    SED=sed
fi


# create the case
echo "Creating case"
echo "----------------------"
echo "${E3SM_SRC_DIR}/cime/scripts/create_newcase --case ${CASE_DIR} --res ELM_USRDAT --mach ${MACHINE_NAME} --compiler ${COMPILER_NAME} --compset ${COMPSET}"
echo "----------------------"
${E3SM_SRC_DIR}/cime/scripts/create_newcase --case ${CASE_DIR} --res ELM_USRDAT --mach ${MACHINE_NAME} --compiler ${COMPILER_NAME} --compset ${COMPSET}

# cp over example files to the case directory
cp -r ./* ${CASE_DIR}/
cd ${CASE_DIR}

echo " fsurdat = '${CASE_DIR}/${INPUTDATA_FILE}'" >> user_nl_elm

# ATS-specific
if [ "${USE_ATS}" != "FALSE" ]; then
    if [ "${USE_ATS}" == "TRUE" ]; then
	echo " use_ats = .true." >> user_nl_elm
    else
	echo " use_ats_ic = .true." >> user_nl_elm
    fi
    echo " domain_decomp_type = 'ats'" >> user_nl_elm
fi

# make sure there is a clean endline -- an extra doesn't hurt
echo "" >> user_nl_elm



# ELM
./xmlchange MOSART_MODE=NULL,DOUT_S=FALSE,DIN_LOC_ROOT=${ELM_ATS_SRC_DIR}/inputdata
./xmlchange DIN_LOC_ROOT_CLMFORC=${ELM_ATS_SRC_DIR}/inputdata/atm/datm7
./xmlchange ELM_USRDAT_NAME=${NAME}

# try to find the domain.nc file: these are now example specific
./xmlchange ATM_DOMAIN_PATH=${CASE_DIR}
./xmlchange LND_DOMAIN_PATH=${CASE_DIR}
./xmlchange ATM_DOMAIN_FILE=${DOMAIN_FILE}
./xmlchange LND_DOMAIN_FILE=${DOMAIN_FILE}

# set the number of tasks
./xmlchange NTASKS=${NTASKS}
./xmlchange NTASKS_PER_INST=${NTASKS}

./xmlchange PIO_TYPENAME=netcdf
./xmlchange RUN_STARTDATE=2020-10-01
./xmlchange STOP_OPTION=nyears,STOP_N=2
./xmlchange BATCH_SYSTEM=none
./xmlchange DEBUG=TRUE

# setup the case
echo ""
echo "Running case.setup"
echo "----------------------"
./case.setup

# Get run directory and configure ATS paths
if [ "${USE_ATS}" != "FALSE" ]; then
    RUN_DIR=$(./xmlquery RUNDIR --value)
    echo "Run directory: ${RUN_DIR}"

    # Fix XML to point to run directory
    ${SED} -i "s^elm_output_data^${RUN_DIR}^g" ${NAME}.xml

    # Set ATS input paths in user_nl_elm
    echo " ats_inputdir = '${RUN_DIR}'" >> user_nl_elm
    echo " ats_inputfile = '${NAME}.xml'" >> user_nl_elm
fi

echo -e '\nstring(APPEND CPPDEFS " -DCPL_BYPASS -DUSE_ATS_LIB")' >> cmake_macros/universal.cmake

# build
echo ""
echo "Running case.build"
echo "----------------------"
./case.build

# Copy mesh files to run directory
if [ "${USE_ATS}" != "FALSE" ]; then
    echo ""
    echo "Copying mesh files to run directory"
    echo "----------------------"
    RUN_DIR=$(./xmlquery RUNDIR --value)
    if [ -e "${CASE_DIR}/${NAME}.exo" ]; then
        cp ${CASE_DIR}/${NAME}.exo ${RUN_DIR}/
        echo "Copied ${NAME}.exo to ${RUN_DIR}"
    fi
    if [ -e "${CASE_DIR}/${NAME}.h5" ]; then
        cp ${CASE_DIR}/${NAME}.h5 ${RUN_DIR}/
        echo "Copied ${NAME}.h5 to ${RUN_DIR}"
    fi
    if [ -e "${CASE_DIR}/${NAME}.xml" ]; then
        cp ${CASE_DIR}/${NAME}.xml ${RUN_DIR}/
        echo "Copied ${NAME}.xml to ${RUN_DIR}"
    fi
fi


echo ""
echo "Run the case yourself:"
echo "----------------------"
echo "pushd ${CASE_DIR} && ./case.submit"
if [ "$GITHUB_ACTIONS" = "TRUE" ]; then
  ./case.submit --no-batch 
fi 




