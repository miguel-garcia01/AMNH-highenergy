"""
To use JPL ephemeris data you must pip install jplephem
Start times, end times and individual times must be inputted as astropy Time objects
Time steps must be entered as floats in units of days

Returns SkyCoord objects containing coordinates, either equitorial or galactic, of Jupiter at the given time or over the given timerange.

"""

from astropy.coordinates import solar_system_ephemeris, get_body, SkyCoord
from astropy.time import Time
import numpy as np

def equitorial(tstart=None, tend=None, tstep=None, time=None):
    if time is not None and (tstart is not None or tend is not None or tstep is not None):
        raise AssertionError('Must enter one time, or a start time, end time, and time step')

    elif time is not None:
        with solar_system_ephemeris.set('jpl'):
            jup_ephem= get_body('sun',time)

    elif time is None and (tstart is None or tend is None or tstep is None):
        raise AssertionError('Must enter start time, end time, and time step')

    elif time is None and (tstart is not None and tend is not None and tstep is not None):
#        tstart = tstart.mjd + (tstep/2) #use this line if you want times in the center of timeranges, dont use if you want coordinates at the time intervals themselves
        tstart = tstart.mjd # use if you want coordinates at time intervals not at center of intervals
        tend = tend.mjd
        times = np.arange(tstart, tend, tstep, dtype=float)
        times_list = Time(times, format='mjd').fits
        times_list = Time(times_list, format='fits')
    
        jup_ephem = list()

        for date in times_list:
            with solar_system_ephemeris.set('jpl'):
                jup = get_body('sun',date)
                jup_ephem.append(jup)

    return jup_ephem

def galactic(tstart=None, tend=None, tstep=None, time=None):
    if time is not None and (tstart is not None or tend is not None or tstep is not None):
        raise AssertionError('Must enter one time, or a start time, end time, and time step')

    elif time is not None:
        with solar_system_ephemeris.set('jpl'):
            jup_ephem= get_body('sun',time).galactic

    elif time is None and (tstart is None or tend is None or tstep is None):
        raise AssertionError('Must enter start time, end time, and time step')

    elif time is None and (tstart is not None and tend is not None and tstep is not None):
#        tstart = tstart.mjd + (tstep/2)  #use this line if you want times in the center of timeranges, dont use if you want coordinates at the time intervals themselves 
        tstart= tstart.mjd # use if you want coordinates at time intervals not at center of intervals
        tend = tend.mjd
        times = np.arange(tstart, tend, tstep, dtype=float)
        times_list = Time(times, format='mjd').fits
        times_list = Time(times_list, format='fits')

        jup_ephem = list()
    
        for date in times_list:
            with solar_system_ephemeris.set('jpl'):
                jup = get_body('sun',date).galactic
                jup_ephem.append(jup)

    return jup_ephem
