# coding: utf-8
"""
This script reads csv outputted by detect_blind and for a given list of time continuous bbx files, 

This script assumes the input file list is:
(a) all the auto-cors have been coadded
(b) frequency dimension has been sliced

This script performs the following:
(1) finds source rise time and set time using alt/az constraints,
(2) concatenation with time continuous files spaning over rise and set times for each pulsar
(3) Prepares `lf2fil` commands to convert freshly cat'd bbx files to filterbank format.
(3) Prepare prepfold commands suited for each pulsar

We would add the functionality to generate scripts for (a,b) but it's too effort.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime as dt
import astropy.time as at
import astropy.coordinates as asc
import astropy.units as au
from   stations_earthlocations import getEarthLoc
############# OPTIONS
stationID  = 3
ILIST      = "./20200520_XX.list"
IDIR       = "XX"
IPOL       = None
ROOT       = None
DETECTCSV  = "./blind_detect_psr.csv"
BEAM_W     = 23.0 * au.degree
DATEFMT    = "%Y%m%d_%H%M%S_{pol}.bbx.gz\n"
### switches
DO_CONCAT  = True
DO_TOFIL   = True
DO_PF      = True
### Job file
CONCAT_JOB = "{root}_concat.job"
TOFIL_JOB  = "{root}_tofil.job"
PF_JOB     = "{root}_pf.job"
### CMDs
CONCAT_CMD = "lfcat -i"
TOFIL_CMD  = "lf2fil {0}/{root}_pf{pf}.{ix} {1}/{root}_pf{pf}.{jx}\n"
PF_CMD     = "prepfold -ncpus 3 -filterbank -zerodm -noxwin -coarse -dm {dm:6.3f} -p {p:6.3f} -nsub 256 -n 64 -npart 100 -o {odir}/{root}_pf{pf}_pfd  {idir}/{root}_pf{pf}.fil\n"
### fixed options
CCDIR      = "PF"
CCNAME     = "{0}/{root}_pf{pf}.{ext}"
CCLIST     = []
############# STEP 0 
## read what is given
"""
this csv is made using detect_blind
"""
df           = pd.read_csv (DETECTCSV,)
N,D          = df.shape
############# STEP 0.1 
## prepare skycoords
df['sc']     = asc.SkyCoord (asc.Angle(df['RAJ'], unit=au.hour), asc.Angle(df['DECJ'], unit=au.degree), frame='icrs')
############# STEP 0.2 
## prepare trange, root and ipol
with open (ILIST, "r") as f:
    IL       = f.readlines ()
### resolve input pol
if IPOL is None:
    IPOL     = IL[0].strip().split("_")[2].split(".")[0]
### resolve root 
if ROOT is None:
    ROOT = IL[0].strip().split("_")[0]
### prepare
DTs          = []
for i in IL:
    DTs.append (dt.datetime.strptime (i, DATEFMT.format(pol = IPOL)))
trange       = at.Time (DTs)
############# STEP 1
## setup horizon coordinates 
station      = getEarthLoc (stationID)
saltaz       = asc.AltAz (location=station, obstime=trange)
zenith       = asc.AltAz (az=np.zeros(trange.size)*au.degree, alt=90.0*np.ones(trange.size)*au.degree, location=station, obstime=trange)
############# STEP 2
## concat
## find time ranges the psr is in primary beam 
def resolver (w):
    """finds continuous intervals"""
    dw,      = np.where (np.diff (w) != 1)
    # find points of discontinuity
    ll       = 0
    hh       = w.size -1
    # list of continuous intervals
    # rr for rise indices
    # ss for set  set indices
    rr       = []
    ss       = []
    for idw in dw:
        if idw == 0:
            ll = 1
            continue
        hh   = idw + 1
        rr.append ( ll )
        ss.append ( hh )
        ll   = hh
    rr.append ( ll )
    hh       = w.size - 1
    ss.append ( hh )
    return rr,ss
# JJ is loop variable
JJ = 0 
with open (CONCAT_JOB.format (root = ROOT), "w") as f:
    for ii in range (N):
        aird     = df.loc[ii,'sc'].transform_to (saltaz)
        sep      = []
        # find when in beam
        for a,z in zip (aird, zenith):
            sep.append (a.separation (z) <= BEAM_W)
        ww,      = np.where ( sep )
        # compute rise index and set index
        uu,vv    = resolver (ww)
        # run on zip
        for iu, iv  in zip (uu,vv):
            tc_cmd   = CONCAT_CMD + " "
            for ii in range (iu, iv+1):
                tc_cmd = tc_cmd + " " + os.path.join (IDIR, IL[ii].strip())
            # output name
            occ =  CCNAME.format (CCDIR, root = ROOT, pf = JJ, ext = "bbx.gz") + "\n"
            CCLIST.append (occ)
            JJ = JJ + 1
            f.write (tc_cmd + " " + occ)
############# STEP 2.1
## to fil?
if DO_TOFIL:
    with open (TOFIL_JOB.format (root = ROOT), "w") as f:
        for ii in range (JJ):
            f.write ( 
                TOFIL_CMD.format (CCDIR, CCDIR, root = ROOT, pf = ii, ix = "bbx.gz", jx = "fil")
                )
############# STEP 3
## prepfold
if DO_PF:
    with open (PF_JOB.format (root = ROOT), "w") as f:
        for ii in range (JJ):
            f.write (
            PF_CMD.format (
                root = ROOT, pf = ii, 
                idir = "PF", odir = "PF",
                dm = df.loc[ii, 'DM'], p = df.loc[ii, 'P0'],
            ))
