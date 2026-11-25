#!/usr/bin/env bash
if [ -z "${USE_ATS}" ]; then
    echo "Set USE_ATS prior to running the build_example script"
    exit 1
fi

if [ -z "${NTASKS}" ]; then
    echo "Setting NTASKS = 1"
    NTASKS=1
fi

NTASKS=${NTASKS} \
      EXAMPLE_DIR=oakharbor_transect \
      CASE_NAME=oakharbor_transect.np${NTASKS} \
      ATS_CASE_NAME=oakharbor_transect \
      INPUTDATA_NAME=5x1pt_Oakharbor \
      DOMAIN_NAME=oakharbor_transect \
      DOMAIN_FILE=domain.lnd.oakharbor_transect.nc \
      SURF_DATA_FILE=surfdata_oakharbor_transect.nc \
      ../build_example.sh

