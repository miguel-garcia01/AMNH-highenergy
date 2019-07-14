import time
import math
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd
from astropy.io import fits 
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, get_body 
from gt_apps import filter, evtbin, maketime, gtexpcube2, expCube, expMap
from coordinates import equatorial
from GtApp import GtApp
start_time = time.time()

expMap = GtApp('gtexpmap','Likelihood')
expCube2 = GtApp('gtexpcube2','Likelihood')

#files
data_file = 'L1907101448550182BCD945_PH00.fits'
sc_file = 'L1907101448550182BCD945_SC00.fits'

#get parameters from Fermi photon event file
hdulist = fits.open(data_file)
hdu = hdulist[1]

roi_data = (hdu.header['DSVAL3']).split(',')
roi = int(float(((roi_data[2])[:-1]))) #degrees
time_start = Time(hdu.header['DATE-OBS'], format='fits')
time_end = Time(hdu.header['DATE-END'], format='fits') 
time_step = 0.5 #in units of days
time_start_MET = hdu.header['TSTART'] #MET (s)
time_end_MET = hdu.header['TSTOP'] #MET (s)

times = np.arange(time_start_MET, time_end_MET, time_step*(24*60*60), dtype=float)
#times = np.append(times, time_end_MET)

energyrange = (hdu.header['DSVAL5']).split(':')
emin = float(energyrange[0]) #MeV
emax = float(energyrange[1]) #MeV
num_bins = int(math.log10(emax/emin)*10)
zmax = 90 #deg
bin_size = 0.2 #deg/pix
dim = (2*roi)/bin_size #gives dimension of cmap diameter/bin size = num pixels

hdulist.close()


#make array that will hold stacked counts map data
total_cmap_data = np.zeros((int(dim),int(dim)), dtype = float)
#total_ccube_data = np.zeros((20,int(dim),int(dim)), dtype = float)
#total_expcube2_data = np.zeros((int(dim),int(dim)), dtype = float)
total_expcube_data = np.zeros((21, int(dim), int(dim)), dtype = float)
#make list of SkyCoord objects with coordinate data
skycoords = equatorial(obj='sun', tstart=time_start, tend=time_end, tstep=time_step)


#loop through times to create stacked counts map
for val in range(1,len(skycoords)):
     ra = (skycoords[val-1].ra).degree
     dec = (skycoords[val-1].dec).degree
     
     #gtselect
     filter['evclass'] = 128
     filter['evtype'] = 3
     filter['rad'] = roi
     filter['emin'] = emin
     filter['emax'] = emax
     filter['zmax'] = zmax
     filter['infile'] = data_file
     filter['ra'] = ra
     filter['dec'] = dec
     filter['tmin'] = times[val-1] #must be in MET (s)
     filter['tmax'] = times[val]
     filtered_file = 'SUN_filtered_' + str(val) + '.fits'
     filter['outfile'] = filtered_file
     filter.run()

     #gtmktime
     maketime['scfile'] = sc_file
     maketime['filter'] = '(DATA_QUAL==1)&&(LAT_CONFIG==1)&&(ROCK_ANGLE)<52'
     maketime['apply_filter'] = 'yes'
     maketime['roicut'] = 'no'
     maketime['evfile'] = filtered_file
     mktime_file = 'SUN_mktime_' + str(val) + '.fits'
     maketime['outfile'] = mktime_file
     maketime.run()
    
     #counts cube
     evtbin['evfile'] = mktime_file
     ccube_file = 'Sun_ccube_' + str(val) + '.fits'
     evtbin['outfile'] = ccube_file
     evtbin['algorithm'] = 'CCUBE'
     evtbin['nxpix'] = int(dim)
     evtbin['nypix'] = int(dim)
     evtbin['binsz'] = bin_size
     evtbin['coordsys'] = 'CEL'
     evtbin['xref'] = ra
     evtbin['yref'] = dec
     evtbin['axisrot'] = 0.0
     evtbin['proj'] = 'AIT'
     evtbin['ebinalg'] = 'LOG'
     evtbin['emin'] = emin
     evtbin['emax'] = emax
     evtbin['enumbins'] = num_bins
     evtbin.run()
     
     #livetime cube
     expCube['evfile'] = mktime_file
     expCube['scfile'] = sc_file
     ltCube_file = 'Sun_ltCube_' + str(val) + '.fits'
     expCube['outfile'] = ltCube_file
     expCube['zmax'] = zmax
     expCube['dcostheta'] = 0.025
     expCube['binsz'] = 1.0
     expCube.run()


     """#exposure map
     expMap['evfile'] = mktime_file 
     expMap['scfile'] = sc_file
     expMap['expcube'] = ltCube_file
     expMap_file = 'Sun_expmap_' + str(val) + '.fits'
     expMap['outfile'] = expMap_file
     expMap['irfs'] = 'CALDB'
     expMap['srcrad'] = roi
     expMap['nlong'] = 120
     expMap['nlat'] = 120
     expMap['nenergies'] = 20
     expMap.run()"""


     #exposure cube creates exposure map
     expCube2['infile'] = ltCube_file 
     expCube2['cmap'] = ccube_file
     expmap_file = 'Sun_expmap_' + str(val) + '.fits'
     expCube2['outfile'] = expmap_file
     expCube2['irfs'] = 'P8R3_SOURCE_V2'
     expCube2['nxpix'] = int(dim)
     expCube2['nypix'] = int(dim)
     expCube2['binsz'] = bin_size
     expCube2['coordsys'] = 'CEL'
     expCube2['xref'] = ra
     expCube2['yref'] = dec
     expCube2['axisrot'] = 0.0
     expCube2['proj'] = 'CAR'
     expCube2['ebinalg'] = 'LOG'
     expCube2['emin'] = emin
     expCube2['emax'] = emax
     expCube2['enumbins'] = num_bins
     expCube2.run()

     
     #stack each exposure map                                                  
     total_expmap_data += fits.getdata(expmap_file)

#cmap_outfile = 'stacked_cmap.fits'
#hdu = fits.PrimaryHDU(total_cmap_data)

expmap_outfile = 'stacked_expmap.fits'
hdu1 = fits.PrimaryHDU(total_expmap_data)
hdu1.writeto(expmap_outfile, overwrite=True)


"""plt.imshow(total_cmap_data, cmap='gray')
plt.colorbar()
plt.show()"""

print('Runtime:',time.time() - start_time)
