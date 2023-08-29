#!/usr/bin/env bash

# set NCOLS=1 to run single-column example
# set NCOLS=5 to run 5-column hillslope example
export NCOLS=5

# set compset name
export COMPSET=ICB1850CNPRDCTCBC
export COMPDIR="${PROJECT_E3SM}/cases/${COMPSET}"
export RUNDIR="${PROJECT_E3SM}/scratch/${COMPSET}/run"

cd ${E3SM_ROOT}/cime/scripts && 
./create_newcase --case ${COMPDIR} --res ELM_USRDAT --mach ${MACH_NAME} --compiler gnu --compset ${COMPSET} --walltime 06:00:00

cd ${COMPDIR} && 
./xmlchange --id ELM_USRDAT_NAME --val "${NCOLS}x1pt_Oakharbor-GRID"
./xmlchange NTASKS=1
./xmlchange NTASKS_PER_INST=1
./xmlchange PIO_TYPENAME=netcdf
./xmlchange RUN_STARTDATE=2000-07-15
./xmlchange STOP_N=100
./xmlchange HIST_N=1

# expect to error here
# part of pt-mode ELM hack
cd ${COMPDIR}; ./case.setup

# exit on error beyond this point
set -e

# write to user_nl_elm
if [[ $NCOLS -eq 1 ]]
then
## single-column
  cd ${COMPDIR} &&
  echo "metdata_type = 'gswp3'
 metdata_bypass = '/Users/80x/Software/ats_newstate/repos/E3SM/pt-e3sm-inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v2.c180716_Oakharbor-Grid/cpl_bypass_full'
 fsurdat = '/Users/80x/Software/ats_newstate/repos/E3SM/pt-e3sm-inputdata/lnd/clm2/surfdata_map/surfdata_1x1pt_Oakharbor-GRID_simyr1850_c360x720_c20230522.nc'
 aero_file = '/Users/80x/Software/ats_newstate/repos/E3SM/pt-e3sm-inputdata/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_monthly_1850_mean_1.9x2.5_c090803.nc'
 CO2_file = '/Users/80x/Software/ats_newstate/repos/E3SM/pt-e3sm-inputdata/atm/datm7/CO2/fco2_datm_1765-2007_c100614.nc'
 nyears_ad_carbon_only = 0
 spinup_mortality_factor = 10
 hist_empty_htapes = .true.
 hist_nhtfrq = -24
 hist_fincl1 = 'TBOT', 'PBOT','RH','RAIN','SNOW', 'TLAI', 'ZWT', 'SMP', 'SOILLIQ', 'SOILICE','SOIL_PRESSURE'
 hist_mfilt = 1
 use_ats = .true.
 ats_inputdir = '/Users/80x/Software/ats_newstate/repos/E3SM/pt-e3sm-inputdata/lnd/clm2/ats'
 ats_inputfile = 'column_jb.xml'
&dynamic_subgrid
 do_harvest = .false.
 do_transient_pfts = .false.
 flanduse_timeseries = ''
/
&finidat_consistency_checks
 check_finidat_fsurdat_consistency = .false.
 check_finidat_year_consistency = .false.
/
 " >> user_nl_elm;
fi

if [[ $NCOLS -eq 5 ]]
then
## 5-column 2D hillslope
  cd ${COMPDIR} &&
  echo "metdata_type = 'gswp3'
 fsurdat = '/Users/80x/Software/ats_newstate/repos/E3SM/pt-e3sm-inputdata/lnd/clm2/surfdata_map/surfdata_5x1pt_Oakharbor-GRID_simyr1850_c360x720_c20230522.nc'
 flanduse_timeseries = '/Users/80x/Software/ats_newstate/repos/E3SM/pt-e3sm-inputdata/lnd/clm2/surfdata_map/surfdata_5x1pt_Oakharbor-GRID_simyr1850_c360x720_c20230522.nc'
 nyears_ad_carbon_only = 25
 spinup_mortality_factor = 10
 hist_empty_htapes = .true.
 hist_nhtfrq = -24
 hist_fincl1 = 'TBOT', 'PBOT','RH','RAIN','SNOW', 'TLAI', 'ZWT', 'SMP', 'SOILLIQ', 'SOILICE','SOIL_PRESSURE'
 hist_mfilt = 365
 metdata_bypass = '/Users/80x/Software/ats_newstate/repos/E3SM/pt-e3sm-inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v2.c180716_Oakharbor-Grid/cpl_bypass_full'
 aero_file = '/Users/80x/Software/ats_newstate/repos/E3SM/pt-e3sm-inputdata/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_monthly_1850_mean_1.9x2.5_c090421.nc'
 CO2_file = '/Users/80x/Software/ats_newstate/repos/E3SM/pt-e3sm-inputdata/atm/datm7/CO2/fco2_datm_1765-2007_c100614.nc'
 use_ats = .true.
 ats_inputdir = '/Users/80x/Software/ats_newstate/repos/E3SM/pt-e3sm-inputdata/lnd/clm2/ats'
 ats_inputfile = 'hillslope_jb.xml'
&dynamic_subgrid
 do_harvest = .false.
 do_transient_pfts = .false.
 flanduse_timeseries = ''
/
&finidat_consistency_checks
 check_finidat_fsurdat_consistency = .false.
 check_finidat_year_consistency = .false.
/
 " >> user_nl_elm;
fi

# setup again - this time with correct nl variables in place
cd ${COMPDIR}; ./case.setup

# build
cd ${COMPDIR}; ./case.build
