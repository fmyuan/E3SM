#!/usr/bin/env python

import netcdf4_functions as nffun
import socket, os, sys, csv, time, math, numpy
import re, subprocess
from optparse import OptionParser
import pandas
#from Numeric import *


#runcaseCol3rd_JD.py does the following:
#   this script is modified by Junyan based on runcase.py in the OLMT: 
#       https://github.com/dmricciuto/OLMT/blob/whuang/coastal_3col/runcase.py
# to configure col3rd ELM version
# 1. modify parameter file to include tide pars
# 2. add tide forcing file in ser_nl_elm 
# 3. use ./xmlchange to add definition of col3rd and cp_bypass in pre-process compiling to turn 3d version on
#  Note: this script is written for modifying parameter file 
#------------------------------------------------------------------
#Pre-setup

mylsm='ELM'
model_name='elm'
parm_file = '/clm_params_c180301.nc'
tmpdir = '/compyfs/ding567/COMPASS/InputData/LakeErie_Wei3col/lnd/clm2/paramdata' 
myncap='ncap2'  # define NetCDF file modification command
tide_components_file = '/compyfs/ding567/COMPASS/InputData/LakeErie_Wei3col/lnd/clm2/paramdata/harmonic_Annapolis.csv'
# 1. make a copy of elm default parameter file as clm_params.nc
# 2. modify the parameter file to add tide pars

# os.system('nccopy -3 '+tmpdir+parm_file+' '+tmpdir+'/clm_params.nc')


flnr = nffun.getvar(tmpdir+'/clm_params.nc','flnr')
os.system(myncap+' -O -s "humhol_ht = br_mr*0+0.15" '+tmpdir+'/clm_params.nc '+tmpdir+'/clm_params.nc')
os.system(myncap+' -O -s "humhol_ht_frac = br_mr*0+1" '+tmpdir+'/clm_params.nc '+tmpdir+'/clm_params.nc')
os.system(myncap+' -O -s "hum_frac = br_mr*0+0.50" '+tmpdir+'/clm_params.nc '+tmpdir+'/clm_params.nc')
os.system(myncap+' -O -s "qflx_h2osfc_surfrate = br_mr*0+0.0" '+tmpdir+'/clm_params.nc '+tmpdir+'/clm_params.nc')
  
os.system(myncap+' -O -s "rsub_top_globalmax = br_mr*0+1.2e-5" '+tmpdir+'/clm_params.nc '+tmpdir+'/clm_params.nc')
os.system(myncap+' -O -s "h2osoi_offset = br_mr*0" '+tmpdir+'/clm_params.nc '+tmpdir+'/clm_params.nc')
#   flnr = nffun.getvar(tmpdir+'/clm_params.nc','flnr')
#   os.system(myncap+' -O -s "br_mr = flnr" '+tmpdir+'/clm_params.nc '+tmpdir+'/clm_params.nc')
#   ierr = nffun.putvar(tmpdir+'/clm_params.nc','br_mr', flnr*0.0+2.52e-6)
#if ((options.marsh or options.col3rd) and options.tide_components_file != ''):


tidecomps=pandas.read_csv(tide_components_file)
for comp in range(len(tidecomps)):
    os.system(myncap+' -O -s "tide_coeff_amp_%d = humhol_ht*0+%1.4e" '%(comp+1,tidecomps['Amplitude'].iloc[comp]*1000)+tmpdir+'/clm_params.nc '+tmpdir+'/clm_params.nc')
    #os.system(myncap+' -O -s "tide_coeff_period_%d = humhol_ht*0+%1.4e" '%(comp+1,360*3600/tidecomps['Speed'].iloc[comp])+tmpdir+'/clm_params.nc '+tmpdir+'/clm_params.nc')
    #editted by Wei Huang for the 1 line below to correct the tidal period
    os.system(myncap+' -O -s "tide_coeff_period_%d = humhol_ht*0+%1.4e" '%(comp+1,3600/tidecomps['Speed'].iloc[comp])+tmpdir+'/clm_params.nc '+tmpdir+'/clm_params.nc')
    os.system(myncap+' -O -s "tide_coeff_phase_%d = humhol_ht*0+%1.4e" '%(comp+1,tidecomps['Phase'].iloc[comp]*math.pi/180)+tmpdir+'/clm_params.nc '+tmpdir+'/clm_params.nc')
os.system(myncap+' -O -s "tide_baseline = humhol_ht*0+0.0" '+tmpdir+'/clm_params.nc '+tmpdir+'/clm_params.nc')

os.system(myncap+' -O -s "crit_gdd1 = flnr" '+tmpdir+'/clm_params.nc '+tmpdir+'/clm_params.nc')
os.system(myncap+' -O -s "crit_gdd2 = flnr" '+tmpdir+'/clm_params.nc '+tmpdir+'/clm_params.nc')
ierr = nffun.putvar(tmpdir+'/clm_params.nc','crit_gdd1', flnr*0.0+4.8)
ierr = nffun.putvar(tmpdir+'/clm_params.nc','crit_gdd2', flnr*0.0+0.13)

# BSulman: These Nfix constants can break the model if they don't have the right length.
#os.system(myncap+' -O -s "Nfix_NPP_c1 = br_mr*+1.8" '+tmpdir+'/clm_params.nc '+tmpdir+'/clm_params.nc')
#os.system(myncap+' -O -s "Nfix_NPP_c2 = br_mr*0.003" '+tmpdir+'/clm_params.nc '+tmpdir+'/clm_params.nc')
#os.system(myncap+' -O -s "nfix_timeconst = br_mr*10.0" '+tmpdir+'/clm_params.nc '+tmpdir+'/clm_params.nc')

# JD1 end
os.system('chmod u+w ' +tmpdir+'/clm_params.nc')

#Added option for , 3rd column [Wei Huang 2022-07-11]
# This is done , Junyan added the ./xmlchange -id ELM_CONFIG_OPTS --append --val 'cppdefs -DCOL3RD' in the run script

# print("Turning on  modification\n")
# os.system("./xmlchange -id "+mylsm+"_CONFIG_OPTS --append --val '-cppdefs -DCOL3RD'")


# END
