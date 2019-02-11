import matplotlib.pyplot as plt
import numpy as np
import pdb

inputdir = '/Volumes/seagate_8/data/sami3/jan_2014/'
dims = 'glat', 'glon'

# jan:set no
tec_shape = 672, 96, 100
# jul
# tec_shape = 100, 96, 670A


with open(inputdir + 'time.dat', 'rb') as f: 
    time = np.fromfile(f,dtype='float32')[1:-1]
with open(inputdir + 'zalt0B.dat', 'rb') as f: 
    alt = np.fromfile(f,dtype='float32')[1:-1]
dim_shape = (tec_shape[1], len(alt), tec_shape[2])  # lon, alt, lat

with open(inputdir + 'glon0B.dat', 'rb') as f: 
    lon = np.fromfile(f, dtype='float32')[1:-1]
    lon.shape = dim_shape 
with open(inputdir + 'glat0B.dat', 'rb') as f: 
    lat = np.fromfile(f, dtype='float32')[1:-1]
    lat.shape = dim_shape
with open(inputdir + 'tecuB.dat', 'rb') as f: 
    tec = np.fromfile(f, dtype='float32')[1:-1]
    tec.shape = tec_shape


"""
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
colat_r = 90 - lat[:, 0, :].flatten()
lon_r = lon[:, 0, :].flatten()

im = ax.tricontourf(lon_r, colat_r, tec[100, :, :].flatten(), 50, vmin=0, vmax=90)
ax.set_rlim(0, 45)
fig.colorbar(im)
plt.show()
"""
"""
phi = {}
for dim in dims:
    with open(dim + fname_tec, 'rb') as f: 
        phi[dim] = np.fromfile(f, dtype='float32')[1:-1]
        phi[dim].shape = (96, 100, 100)
with open('phi' + fname_tec, 'rb') as f: 
    phi['vals'] = np.fromfile(f, dtype='float32')[1:-1]
    phi['vals'].shape = (len(time), 100, 96)


"""


