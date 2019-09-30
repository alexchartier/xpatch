import nvector as nv
import numpy as np


wgs84 = nv.FrameE(name='WGS84')
mcm = wgs84.GeoPoint(latitude=-77.8564, longitude=166.6881, z=0, degrees=True).to_ecef_vector().pvector
zsp = wgs84.GeoPoint(latitude=-90.0, longitude=166.6881, z=0, degrees=True).to_ecef_vector().pvector
mcm_zsp_strt = np.sqrt(np.sum((mcm - zsp) ** 2)) / 1E3 

print("MCM-ZSP dist: %1.1f km" % mcm_zsp_strt)
