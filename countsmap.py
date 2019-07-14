import time
import math
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd
from astropy.io import fits 
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, get_body 
from gt_apps import filter, evtbin, maketime
from coordinates import equatorial

#files
data_file = 'L1907101406300182BCD940_PH00.fits'
sc_file = 'L1907101406300182BCD940_SC00.fits'

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

#make list of SkyCoord objects with coordinate data
skycoords = equatorial(obj='sun', tstart=time_start, tend=time_end, tstep=time_step)


#loop through times to create stacked counts map
for val in range(1,len(skycoords)):
     ra = (skycoords[val-1].ra).degree
     dec = (skycoords[val-1].dec).degree
     
     #gtselect
     filter['evclass'] = 128
     filter['evtype'] = 3
     filter['rad'] = 'INDEF'
     filter['emin'] = emin
     filter['emax'] = emax
     filter['zmax'] = zmax
     filter['infile'] = data_file
     filter['ra'] = ra
     filter['dec'] = dec
     filter['tmin'] = times[val-1] #must be in MET (s)
     filter['tmax'] = times[val]
     outfile = 'SUN_binned_filtered_' + str(val) + '.fits'
     filter['outfile'] = outfile
     filter.run()

     #gtmktime
     maketime['scfile'] = sc_file
     maketime['filter'] = '(DATA_QUAL==1)&&(LAT_CONFIG==1)&&(ROCK_ANGLE)<52'
     #'(DATA_QUAL==1)&&(LAT_CONFIG==1)&&(ABS(ROCK_ANGLE)<52)'
     maketime['roicut'] = 'yes'
     maketime['evfile'] = outfile
     gti_file = 'SUN_binned_gti_' + str(val) + '.fits'
     maketime['outfile'] = gti_file
     maketime.run()

     #counts map
     evtbin['scfile'] = sc_file
     evtbin['nxpix'] = int(dim)
     evtbin['nypix'] = int(dim)
     evtbin['binsz'] = bin_size
     evtbin['axisrot'] = 0
     evtbin['coordsys'] = 'CEL'
     evtbin['proj'] = 'AIT'
     evtbin['algorithm'] = 'CMAP'
     evtbin['evfile'] = gti_file
     cmap_file = 'SUN_binned_cmap_'+ str(val) + '.fits'
     evtbin['outfile'] = cmap_file
     evtbin['xref'] = ra
     evtbin['yref'] = dec
     evtbin.run() 
    
     #stack each cmap 
     total_cmap_data += fits.getdata(cmap_file)

cmap_outfile = 'stacked_cmap.fits'
hdu = fits.PrimaryHDU(total_cmap_data)
hdu.writeto(cmap_outfile, overwrite=True)

plt.imshow(total_cmap_data, cmap='gray')
plt.colorbar()
plt.show()
