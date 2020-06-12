# coding: utf-8
"""
This script applies three selection logic for every pulsar to determine if it is detectable 
(1) if it's flux density at LoFASM frequency range is atleast `MINMJY`
(2) if the source is observed for atleast `NSUB` pulses accounting for the DM delay in a continous stretch.
(3) if the source's altitude/azimuth  is in the LoFASM's primary beam width.
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astropy.time as at
import astropy.coordinates as asc
import astropy.units as au
from stations_earthlocations import getEarthLoc
############# OPTIONS
WRITECSV   = True
stationID  = 3
OBSTIME    = "2020-05-20T00:00:31.000"
MINMJY     = 1000
MINBIN     = 8
NSUB       = 100
ALPHA      = 1.25
DM_DELAY   = 4.15E-3 * (0.020**-2 - 0.080**-2)
TSTEP      = 10 * au.minute
BEAM_W     = 23.0 * au.degree
TSAMP      = 83.33E-3
############# STEP 0 
## Read psrs
"""
this csv is made using psrcat
$ psrcat -c "JName RAJ DECJ DM P0 S400" -s S400 -nonumber -o short_csv | grep -v "*" > all_psr_dm_p0_s400.csv
Grep is needed to remove entries with "*" in S400
Have to remove the units line manually.
"""
df           = pd.read_csv ("all_psr_dm_p0_s400.csv", sep=";",)
N,D          = df.shape
############# STEP 1 
## prepare stats
df['nbin']   = df['P0'] / TSAMP
df['tP0']    = NSUB * df['P0'] * au.second
df['dmd']    = DM_DELAY * df['DM'] * au.second
df['treq']   = df['tP0'] + df['dmd']
df['s100']   = df['S400'] * 4**ALPHA
df['sc']     = asc.SkyCoord (asc.Angle(df['RAJ'], unit=au.hour), asc.Angle(df['DECJ'], unit=au.degree), frame='icrs')
"""
nbin is number of bins possible
tP0  is the time length of `nsub` periods
dmd  is the DM delay 
treq is the total observing time required for a pulsar assuming we would need `nsub` folding.
s100 is the flux den scaled using power law with alpha.
sc   is skycoords
"""
############# STEP 2
## setup horizon coordinates 
tstart       = at.Time (OBSTIME, format='isot')
nstep        = 1.0 / ( TSTEP.to (au.day) )
trange       = tstart + (TSTEP*np.arange(nstep.value, dtype=np.uint32))
station      = getEarthLoc (stationID)
saltaz       = asc.AltAz (location=station, obstime=trange)
zenith       = asc.AltAz (az=np.zeros(trange.size)*au.degree, alt=90.0*np.ones(trange.size)*au.degree, location=station, obstime=trange)
############# STEP 3
## find total time the psr is in primary beam 
## time is `nt` * TSTEP
def resolver (w):
    """finds the size of smallest and continuous interval"""
    dw,      = np.where (np.diff (w) != 1)
    # find points of discontinuity
    ll       = 0 
    hh       = w.size - 1
    # sum continuous intervals
    rr       = []
    for idw in dw:
        if idw == 0:
            ll = 1
            continue
        hh   = idw + 1
        rr.append (np.sum(w[ll:hh]))
        ll   = hh
    rr.append (np.sum(w[ll:hh]))
    return np.min (rr)

df['nt']     = 0 * TSTEP
for ii in range (N):
    aird     = df['sc'].iloc[ii].transform_to (saltaz)
    sep      = []
    for a,z in zip (aird, zenith):
        sep.append (a.separation (z) <= BEAM_W)
    ww,      = np.where ( sep )
    uu       = resolver (ww) * TSTEP
    df.loc[ii,'nt'] = uu.to(au.second).value 
df['tyn']    =  df['nt'] >= df['treq']
############# STEP 4
## minimum s100 cut off
## minimum bins
df['syn']    = df['s100'] >= MINMJY 
df['byn']    = df['nbin'] >= MINBIN
############# STEP 5
yn           = df['tyn'] & df['syn'] & df['byn']
qdf          = df.loc[yn,['PSRJ', 'RAJ', 'DECJ', 'DM', 'P0', 's100']]
if WRITECSV:
    qdf.to_csv ("blind_detect_psr.csv", index=False, )
else:
    print (qdf)
