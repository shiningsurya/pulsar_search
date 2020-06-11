# coding: utf-8
"""
This script reads csv outputted by detect_blind and for a given 
stationid and observing time, finds source rise time and set time
using alt/az constraints.
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astropy.time as at
import astropy.coordinates as asc
import astropy.units as au
from   stations_earthlocations import getEarthLoc
############# OPTIONS
WRITECSV   = True
stationID  = 3
OBSTIME    = "2020-05-20T00:00:31.000"
DETECTCSV  = "./blind_detect_psr.csv"
ALPHA      = 1.25
DM_DELAY   = 4.15E-3 * (0.020**-2 - 0.080**-2)
TSTEP      = 10 * au.minute
BEAM_W     = 23.0 * au.degree
############# STEP 0 
## read what blind search gave
"""
this csv is made using detect_blind
"""
df           = pd.read_csv (DETECTCSV,)
N,D          = df.shape
############# STEP 0.1 
## prepare skycoords
df['sc']     = asc.SkyCoord (asc.Angle(df['RAJ'], unit=au.hour), asc.Angle(df['DECJ'], unit=au.degree), frame='icrs')
############# STEP 1
## setup horizon coordinates 
tstart       = at.Time (OBSTIME, format='isot')
nstep        = 1.0 / ( TSTEP.to (au.day) )
trange       = tstart + (TSTEP*np.arange(nstep.value, dtype=np.uint32))
station      = getEarthLoc (stationID)
saltaz       = asc.AltAz (location=station, obstime=trange)
zenith       = asc.AltAz (az=np.zeros(trange.size)*au.degree, alt=90.0*np.ones(trange.size)*au.degree, location=station, obstime=trange)
############# STEP 2
## find time ranges the psr is in primary beam 
"""
Write output in pd.DataFrame
PSRJ DM P0 RISE SET
"""
odf          = pd.DataFrame (columns=['PSRJ', 'DM', 'P0', 'RISE', 'SET'])
## time is `nt` * TSTEP
def resolver (w, tstart, tstep):
    """finds continuous intervals"""
    dw,      = np.where (np.diff (w) != 1)
    # find points of discontinuity
    ll       = w[0]
    hh       = w[-1]
    # list of continuous intervals
    # rr for rise times
    # ss for set  times
    rr       = []
    ss       = []
    for idw in dw:
        if idw == 0:
            ll = 1
            continue
        hh   = idw + 1
        rr.append ( tstart + ( ll * tstep ) )
        ss.append ( tstart + ( hh * tstep ) )
        ll   = hh
    rr.append ( tstart + ( ll * tstep ) )
    hh       = w[w.size -1] 
    ss.append ( tstart + ( hh * tstep ) )
    return rr,ss
# JJ is loop variable
JJ = 0 
for ii in range (N):
    aird     = df.loc[ii,'sc'].transform_to (saltaz)
    sep      = []
    for a,z in zip (aird, zenith):
        sep.append (a.separation (z) <= BEAM_W)
    ww,      = np.where ( sep )
    uu,vv    = resolver (ww, tstart, TSTEP)
    for iu,iv in zip(uu,vv):
        odf.loc[JJ] = df.loc[ii, ['PSRJ', 'DM', 'P0']].tolist() + [iu, iv]
        JJ   = JJ + 1
if WRITECSV:
    odf.to_csv ("target_detect_psr.csv", index=False, )
else:
    print (odf)
