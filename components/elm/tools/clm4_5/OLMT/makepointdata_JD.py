#!/usr/bin/env python3

# Junyan changed the code to create domain and surface data for 3col ELM
import os, sys, csv, time, math
from optparse import OptionParser
import numpy
import netcdf4_functions as nffun

#---------------------Create domain data --------------------------------------------------
n_grids = 3
issite = True

# ===== the following codes might not be correct ======
resx=0.5
resy=0.5
xgrid_min = [0, 0, 0]
xgrid_max = [0, 0, 0]
ygrid_min = [0, 0, 0]
ygrid_max = [0, 0, 0]
# define lat and long for the three column
lon = [276.936, 276.9361,  276.9362 ]
lat = [41.485, 41.485, 41.485]
for n in range(0,n_grids):
   for i in range(0,1):
      xgrid_min[n] = i
      xgrid_max[n] = i
      ygrid_min[n] = i
      ygrid_max[n] = i
      
# ================================== end ==========

domainfile_orig = '/compyfs/ding567/COMPASS/InputData/LakeErie_Wei3col/domain_LEptr_from_CAGouldenCZ2a.nc'

domainfile_list=''
for n in range(0,n_grids):
    nst = str(100000+n)[1:]
    domainfile_new = '/compyfs/ding567/COMPASS/InputData/LakeErie_Wei3col/domain_LEptr_from_CAGouldenCZ2'+nst+'.nc'
#    if (not os.path.exists(domainfile_orig)):
#        print('Error:  '+domainfile_orig+' does not exist.  Aborting')
#        sys.exit(1)

#    if (isglobal):
#        os.system('cp '+domainfile_orig+' '+domainfile_new)
#    else:
    os.system('ncks -d ni,'+str(xgrid_min[n])+','+str(xgrid_max[n])+' -d nj,'+str(ygrid_min[n])+ \
          ','+str(ygrid_max[n])+' '+domainfile_orig+' '+domainfile_new)

#    if (issite):
    frac = nffun.getvar(domainfile_new, 'frac')
    mask = nffun.getvar(domainfile_new, 'mask')
    xc = nffun.getvar(domainfile_new, 'xc')
    yc = nffun.getvar(domainfile_new, 'yc')
    xv = nffun.getvar(domainfile_new, 'xv')
    yv = nffun.getvar(domainfile_new, 'yv')
    area = nffun.getvar(domainfile_new, 'area')
    frac[0] = 1.0
    mask[0] = 1
#    if (options.site != ''):
    xc[0] = lon[n]
    yc[0] = lat[n]
    xv[0][0][0] = lon[n]-resx/2
    xv[0][0][1] = lon[n]+resx/2
    xv[0][0][2] = lon[n]-resx/2
    xv[0][0][3] = lon[n]+resx/2
    yv[0][0][0] = lat[n]-resy/2
    yv[0][0][1] = lat[n]-resy/2
    yv[0][0][2] = lat[n]+resy/2
    yv[0][0][3] = lat[n]+resy/2
    area[0] = resx*resy*math.pi/180*math.pi/180
    ierr = nffun.putvar(domainfile_new, 'xc', xc)
    ierr = nffun.putvar(domainfile_new, 'yc', yc)
    ierr = nffun.putvar(domainfile_new, 'xv', xv)
    ierr = nffun.putvar(domainfile_new, 'yv', yv)
    ierr = nffun.putvar(domainfile_new, 'area', area)
   
    ierr = nffun.putvar(domainfile_new, 'frac', frac)
    ierr = nffun.putvar(domainfile_new, 'mask', mask)
    os.system('ncks -O --mk_rec_dim nj '+domainfile_new+' '+domainfile_new)
#    elif (options.mymask != ''):
#       print 'Applying mask from '+options.mymask
#       os.system('ncks -d lon,'+str(xgrid_min[n])+','+str(xgrid_max[n])+' -d lat,'+str(ygrid_min[n])+ \
#              ','+str(ygrid_max[n])+' '+options.mymask+' mask_temp.nc')
#       newmask = nffun.getvar('mask_temp.nc', 'PNW_mask')
#       ierr = nffun.putvar(domainfile_new, 'mask', newmask)
#       os.system('rm mask_temp.nc')

    domainfile_list = domainfile_list+' '+domainfile_new

domainfile_new = '/compyfs/ding567/COMPASS/InputData/LakeErie_Wei3col/domain_LEptr_from_CAGouldenCZ23col.nc'

if (n_grids > 1):
    os.system('ncrcat '+domainfile_list+' '+domainfile_new)
    os.system('nccopy -u  '+domainfile_new+' '+domainfile_new+'.tmp')
    os.system('ncpdq -O -a ni,nj '+domainfile_new+'.tmp '+domainfile_new)
    #os.system('ncwa -O -a ni -d ni,0,0 '+domainfile_new+'.tmp1 '+domainfile_new+'.tmp2')
    os.system('ncrename -h -O -d ni,ni_temp '+domainfile_new+' '+domainfile_new+' ')
    os.system('ncrename -h -O -d nj,ni '+domainfile_new+' '+domainfile_new+' ')
    os.system('ncrename -h -O -d ni_temp,nj '+domainfile_new+' '+domainfile_new+' ')
    #os.system('rm /compyfs/ding567/COMPASS/InputData/LakeErie_Wei3col/domain_LEptr_from_CAGouldenCZ2?????.nc*')
    #os.system('mv '+domainfile_new+'.tmp3 '+domainfile_new)
    #os.system('rm '+domainfile_new+'.tmp*')
else:
    os.system('mv '+domainfile_list+' '+domainfile_new)

#-------------------- create surface data ----------------------------------
issite = False
mypft = -1
print('Creating surface data')
surffile_orig = '/compyfs/ding567/COMPASS/InputData/LakeErie_Wei3col/surfdata_LEptr_enax4v1a.nc'
surffile_list = ''
for n in range(0,n_grids):
    nst = str(100000+n)[1:]
    surffile_new =  '/compyfs/ding567/COMPASS/InputData/LakeErie_Wei3col/surfdata_LEptr_enax4v1'+nst+'.nc'

    #os.system('ncks --fix_rec_dmn time -d lsmlon,'+str(xgrid_min[n])+','+str(xgrid_max[n])+ \
    #   ' -d lsmlat,'+str(ygrid_min[n])+','+str(ygrid_max[n])+' '+surffile_orig+' '+surffile_new)

    if (issite):
        landfrac_pft = nffun.getvar(surffile_new, 'LANDFRAC_PFT')
        pftdata_mask = nffun.getvar(surffile_new, 'PFTDATA_MASK')
        longxy       = nffun.getvar(surffile_new, 'LONGXY')
        latixy       = nffun.getvar(surffile_new, 'LATIXY')
        area         = nffun.getvar(surffile_new, 'AREA')
        pct_wetland  = nffun.getvar(surffile_new, 'PCT_WETLAND')
        pct_lake     = nffun.getvar(surffile_new, 'PCT_LAKE')
        pct_glacier  = nffun.getvar(surffile_new, 'PCT_GLACIER')
        pct_urban    = nffun.getvar(surffile_new, 'PCT_URBAN')

        soil_order   = nffun.getvar(surffile_new, 'SOIL_ORDER')
        labilep      = nffun.getvar(surffile_new, 'LABILE_P')
        primp        = nffun.getvar(surffile_new, 'APATITE_P')
        secondp      = nffun.getvar(surffile_new, 'SECONDARY_P')
        occlp        = nffun.getvar(surffile_new, 'OCCLUDED_P')

        #input from site-specific information
        soil_color   = nffun.getvar(surffile_new, 'SOIL_COLOR')
        pct_sand     = nffun.getvar(surffile_new, 'PCT_SAND')
        pct_clay     = nffun.getvar(surffile_new, 'PCT_CLAY')
        organic      = nffun.getvar(surffile_new, 'ORGANIC')
        fmax         = nffun.getvar(surffile_new, 'FMAX')
        pct_nat_veg  = nffun.getvar(surffile_new, 'PCT_NATVEG')
        pct_pft      = nffun.getvar(surffile_new, 'PCT_NAT_PFT') 
        monthly_lai  = nffun.getvar(surffile_new, 'MONTHLY_LAI')
        monthly_sai  = nffun.getvar(surffile_new, 'MONTHLY_SAI')
        monthly_height_top = nffun.getvar(surffile_new, 'MONTHLY_HEIGHT_TOP')
        monthly_height_bot = nffun.getvar(surffile_new, 'MONTHLY_HEIGHT_BOT')

        npft = 17

        #read file for site-specific PFT information
        mypft_frac=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        mypct_sand = 0.0 
        mypct_clay = 0.0
 
        #if (options.surfdata_grid == False and options.site != ''):
        AFdatareader = csv.reader(open(ccsm_input+'/lnd/clm2/PTCLM/'+options.sitegroup+'_pftdata.txt','rb'))
        for row in AFdatareader:
            if row[0] == options.site:
                for thispft in range(0,5):
                    mypft_frac[int(row[2+2*thispft])]=float(row[1+2*thispft])
                    if (n==0):
                        for thispft in range(0,5):
                            mypft_frac[int(row[2+2*thispft])]=float(row[1+2*thispft])
                    else:
                        mypft_frac=[100.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

        if (sum(mypft_frac[0:npft]) == 0.0):
            print('*** Warning:  PFT data NOT found.  Using gridded data ***')
     #read file for site-specific soil information
        AFdatareader = csv.reader(open(ccsm_input+'/lnd/clm2/PTCLM/'+options.sitegroup+'_soildata.txt','rb'))
        for row in AFdatareader:
            if row[0] == options.site:
                mypct_sand = row[4]
                mypct_clay = row[5]
        if (mypct_sand == 0.0 and mypct_clay == 0.0):
            print('*** Warning:  Soil data NOT found.  Using gridded data ***')
        #else:
        #  try:
        #    mypft_frac[point_pfts[n]] = 100.0
        #  except NameError:
        #    print('using PFT information from surface data')

        #landfrac_pft[0][0] = 1.0
        #pftdata_mask[0][0] = 1

        #if (options.site != ''):
        # run following codes as iste !='' 
        longxy[0][0] = lon[n]
        latixy[0][0] = lat[n]
        area[0] = 111.2*resy*111.321*math.cos((lon[n]*resx)*math.pi/180)*resx

        #if (not options.surfdata_grid or sum(mypft_frac[0:npft]) > 0.0):
        pct_wetland[0][0] = 0.0
        pct_lake[0][0]    = 0.0
        pct_glacier[0][0] = 0.0

          # if ('US-GC4' or 'US-GC3' or 'US-SPR' in options.site and mysimyr !=2000): #TAO replace SPR with GC4/GC3
               #SPRUCE P initial data
               # this is for 3col 
        soil_order[0][0] = 3
        labilep[0][0]    = 4.0
        primp[0][0]      = 1.0
        secondp[0][0]    = 10.0
        occlp[0][0]      = 5.0

        pct_nat_veg[0][0] = 100.0
        for k in range(0,3):
            pct_urban[k][0][0] = 0.0
        for k in range(0,10):
            if (mypct_sand > 0.0 or mypct_clay > 0.0):
                if (k == 0):
                   print('Setting %sand to '+str(mypct_sand))
                   print('Setting %clay to '+str(mypct_clay))
                pct_sand[k][0][0]   = mypct_sand
                pct_clay[k][0][0]   = mypct_clay
          # if ('US-GC4' or 'US-GC3' or 'US-SPR' in options.site):
        if (k < 8):
            organic[k][0][0] = 130.0
        elif (k == 8):
            organic[k][0][0] = 65.0
        pft_names=['Bare ground','ENF Temperate','ENF Boreal','DNF Boreal','EBF Tropical', \
                   'EBF Temperate', 'DBF Tropical', 'DBF Temperate', 'DBF Boreal', 'EB Shrub' \
                   , 'DB Shrub Temperate', 'BD Shrub Boreal', 'C3 arctic grass', \
                   'C3 non-arctic grass', 'C4 grass', 'Crop','xxx','xxx']
           # muypft = -1 in Wei's setting
           #if (options.mypft >= 0):
           #  print 'Setting PFT '+str(options.mypft)+'('+pft_names[int(options.mypft)]+') to 100%'
           #  pct_pft[:,0,0] = 0.0
           #  pct_pft[int(options.mypft),0,0] = 100.0
           #else:
        for p in range(0,npft):
          if (sum(mypft_frac[0:npft]) > 0.0):
              if (mypft_frac[p] > 0.0):
                  if (p < 16):
                     print( 'Setting PFT '+str(p)+'('+pft_names[p]+') to '+ \
                     str(mypft_frac[p])+'%')
                  else:
                     print( 'Setting PFT '+str(p)+' to '+str(mypft_frac[p])+'%')
              pct_pft[p][0][0] = mypft_frac[p]
          #maxlai = (monthly_lai).max(axis=0)
          for t in range(0,12):
              if (float(options.lai) > 0):
                monthly_lai[t][p][0][0] = float(options.lai)
                #monthly_lai[t][p][j][i] = monthly_lai[t][p][0][0] 
                #monthly_sai[t][p][j][i] = monthly_sai[t][p][0][0]
                #monthly_height_top[t][p][j][i] = monthly_height_top[t][p][0][0]
                #monthly_height_bot[t][p][j][i] = monthly_height_bot[t][p][0][0]

        ierr = nffun.putvar(surffile_new, 'LANDFRAC_PFT', landfrac_pft)
        ierr = nffun.putvar(surffile_new, 'PFTDATA_MASK', pftdata_mask)
        ierr = nffun.putvar(surffile_new, 'LONGXY', longxy)
        ierr = nffun.putvar(surffile_new, 'LATIXY', latixy)
        ierr = nffun.putvar(surffile_new, 'AREA', area)
        ierr = nffun.putvar(surffile_new, 'PCT_WETLAND', pct_wetland)
        ierr = nffun.putvar(surffile_new, 'PCT_LAKE', pct_lake)
        ierr = nffun.putvar(surffile_new, 'PCT_GLACIER',pct_glacier)
        ierr = nffun.putvar(surffile_new, 'PCT_URBAN', pct_urban)

        ierr = nffun.putvar(surffile_new, 'SOIL_ORDER', soil_order)
        ierr = nffun.putvar(surffile_new, 'LABILE_P', labilep)
        ierr = nffun.putvar(surffile_new, 'APATITE_P', primp)
        ierr = nffun.putvar(surffile_new, 'SECONDARY_P', secondp)
        ierr = nffun.putvar(surffile_new, 'OCCLUDED_P', occlp)
        ierr = nffun.putvar(surffile_new, 'SOIL_COLOR', soil_color)
        ierr = nffun.putvar(surffile_new, 'FMAX', fmax)
        ierr = nffun.putvar(surffile_new, 'ORGANIC', organic)
        ierr = nffun.putvar(surffile_new, 'PCT_SAND', pct_sand)
        ierr = nffun.putvar(surffile_new, 'PCT_CLAY', pct_clay)
        ierr = nffun.putvar(surffile_new, 'PCT_NATVEG', pct_nat_veg)
        ierr = nffun.putvar(surffile_new, 'PCT_NAT_PFT', pct_pft)
        ierr = nffun.putvar(surffile_new, 'MONTHLY_HEIGHT_TOP', monthly_height_top)
        ierr = nffun.putvar(surffile_new, 'MONTHLY_HEIGHT_BOT', monthly_height_bot)
        ierr = nffun.putvar(surffile_new, 'MONTHLY_LAI', monthly_lai)
    else:
        if (int(mypft) >= 0):
          pct_pft      = nffun.getvar(surffile_new, 'PCT_NAT_PFT')
          pct_pft[:,:,:] = 0.0
          pct_pft[int(options.mypft),:,:] = 100.0
          ierr = nffun.putvar(surffile_new, 'PCT_NAT_PFT', pct_pft)

    surffile_list = surffile_list+' '+surffile_new

surffile_new = '/compyfs/ding567/COMPASS/InputData/LakeErie_Wei3col/surfdata_LEptr_enax4v13colx.nc'
surffile_final =  '/compyfs/ding567/COMPASS/InputData/LakeErie_Wei3col/surfdata_LEptr_enax4v13col.nc'

if (n_grids > 1):
  os.system('ncecat '+surffile_list+' '+surffile_new)
  #os.system('rm ./temp/surfdata?????.nc*')
  #remove ni dimension
  # clapse gridcell dimension to 1
  os.system('ncwa -O -a gridcell -d gridcell,0,0 '+surffile_new+' '+surffile_new+'.tmp')
  # remove gridcell dimension
  os.system('ncks -O -C -x -v gridcell '+surffile_new+'.tmp '+surffile_new+'.tmpa')
  
  # rearrange dimensions so that record to be the last dimension
  os.system('nccopy -3 -u '+surffile_new+'.tmpa'+' '+surffile_new+'.tmp1')
  os.system('ncpdq -a numurbl,record '+surffile_new+'.tmp1 '+surffile_new+'.tmp2')  
  os.system('ncpdq -a natpft,record '+surffile_new+'.tmp2 '+surffile_new+'.tmp3')
  os.system('ncpdq -a nlevurb,record '+surffile_new+'.tmp3 '+surffile_new+'.tmp4')
  os.system('ncpdq -a lsmpft,record '+surffile_new+'.tmp4 '+surffile_new+'.tmp5')
  os.system('ncpdq -a nlevsoi,record '+surffile_new+'.tmp5 '+surffile_new+'.tmp6')  
  
  # rename record to be gridcell
  #os.system('ncrename -h -O -d record,gridcell '+surffile_new+'.tmp5 '+surffile_new+'.tmp6')
  os.system('ncrename -h -O -d record,gridcell '+surffile_new+'.tmp6 '+surffile_new+'.tmp7')

  os.system('mv '+surffile_new+'.tmp7 '+surffile_final)
 # os.system('rm '+surffile_new+'.tmp*')
else:
  os.system('mv '+surffile_list+' '+surffile_new)


