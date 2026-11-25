#!/usr/bin/env bash
if [ -z "${USE_ATS}" ]; then
    echo "Set USE_ATS prior to running the build_example script"
    exit 1
fi


if [ ! -v NTASKS ]; then
    echo "Setting NTASKS = 1"
    NTASKS=1
fi

NTASKS=${NTASKS} \
      CASE_NAME=tempest_transect.np${NTASKS} \
      ATS_CASE_NAME=tempest_transect \
      INPUTDATA_NAME="110x1pt_US-GC_TransTEMPEST" \
      DOMAIN_NAME=tempest_transect \
      ../build_example.sh
