from __future__ import print_function
# Script to calculate magnetic apex coordinates at a given set of latitudes

from apexpy import Apex
import numpy as np
import pdb

date = 2014.0
alt = 300
A = Apex(date=date)
glat = np.arange(-90, 90, 1)
glon = np.arange(0, 360, 1)
glat_2d, glon_2d = np.meshgrid(glat, glon)
with open('apex_%i_%ikm.txt' % (date, alt), 'w') as f:
    for lat, lon in zip(glat_2d.flatten(), glon_2d.flatten()):
        mlat, mlon = A.convert(lat, lon, 'geo', 'apex', height=alt) 
        f.write("%2.2f, %2.2f, %2.2f, %2.2f\n" % (lat, lon, mlat, mlon))
