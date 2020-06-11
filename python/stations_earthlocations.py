import astropy.coordinates as asc
import astropy.units as au
import numpy as np
"""
Used 
https://viewer.nationalmap.gov/theme/elevation/
for the elevations
"""
### 
BEAM_W    = 23.0 * au.degree
SIDE_DAY  = 23.9344696 * au.hour
###
LF1_LAT   = asc.Latitude ( "26:33:19.676", unit=au.degree)
LF1_LONG  = asc.Longitude( "-97:26:31.174",unit=au.degree)
LF1_ELVF  = ( 10.36 * au.imperial.foot ).to (au.meter)
LF2_LAT   = asc.Latitude ( "34:4:43.497" , unit=au.degree)
LF2_LONG  = asc.Longitude( "-107:37:5.819",unit=au.degree)
LF2_ELVF  = ( 6967.08 * au.imperial.foot ).to (au.meter)
LF3_LAT   = asc.Latitude ( "38:25:59.0" ,  unit=au.degree)
LF3_LONG  = asc.Longitude( "-79:50:23.0",  unit=au.degree)
LF3_ELVF  = ( 2464.89 * au.imperial.foot ).to (au.meter)
LF4_LAT   = asc.Latitude ( "34:12:3.0" ,   unit=au.degree)
LF4_LONG  = asc.Longitude( "-118:10:18.0", unit=au.degree)
LF4_ELVF  = ( 1167.89 * au.imperial.foot ).to (au.meter)
# Earth locations
LF1       = asc.EarthLocation (LF1_LONG, LF1_LAT, LF1_ELVF)
LF2       = asc.EarthLocation (LF2_LONG, LF2_LAT, LF2_ELVF)
LF3       = asc.EarthLocation (LF3_LONG, LF3_LAT, LF3_ELVF)
LF4       = asc.EarthLocation (LF4_LONG, LF4_LAT, LF4_ELVF)
#
def getEarthLoc (stationid=1):
    """Returns Earth Location"""
    if stationid == 1:
        return LF1
    elif stationid == 2:
        return LF2
    elif stationid == 3:
        return LF3
    elif stationid == 4:
        return LF4
    else:
        raise ValueError ("Station ID not recognized")
