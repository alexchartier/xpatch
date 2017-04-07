"""
count_passes.py
Script to count the number of times the swarm satellites pass through each lat/lon bin each day, then sum them over the period. 
Script works on both geomagnetic and geographic bins
"""

from spacepy import pycdf
import pdb 
import numpy as np
import scipy as sp
import datetime as dt
import matplotlib.pyplot as plt 
import glob
import pickle
import sys 
import collections
import proc_swarm_lp
sys.path.insert(0, '/users/chartat1/fusionpp/fusion/')
import physics


def main(ipath='./data/swarm_lp/',
         opath='./data/pass_ct/', 
         time=dt.datetime(2016, 1, 1), 
         step=dt.timedelta(days=1),
         endtime=dt.datetime(2017, 1, 1), 
         sats = ['A', 'B', 'C'],
         save=True):

    while time <= endtime:
        timestr = time.strftime('%Y-%m-%d')
        print(timestr)
        vals = {}
        pass_ct = {}

        for sat in sats:
            print('\nSatellite %s' % sat)
            fname_format = ipath + 'SW_OPER_EFI%s' % sat + '_PL*%Y%m%d*.cdf'
            try:
                fname = glob.glob(time.strftime(fname_format))[0]
                vals[sat] = proc_swarm_lp.load_lp(fname)
                pass_ct[sat] = count_passes(vals[sat])
            except:
                print('Could not count passes for satellite %s on %s' % (sat, timestr))

        if save:
            fout = opath + time.strftime('/pass_%Y%m%d.pkl')
            with open(fout, 'wb') as f:
                pickle.dump(pass_ct, f)
            print('Saving %s' % fout)
        time += dt.timedelta(days=1)

def count_passes(vals, crdtype=['geo', 'mag']):
    # Transform lats/lons to magnetic 
    vals['lat_geo'] *= np.pi / 180
    vals['lon_geo'][vals['lon_geo'] < 0] += 360
    vals['lon_geo'] *= np.pi / 180
    alts, vals['lat_mag'], vals['lon_mag'] = physics.transform(vals['rad'], vals['lat_geo'], \
                          vals['lon_geo'], from_=['GEO', 'sph'], to=['MAG', 'sph'])
    vals['lon_mag'][vals['lon_mag'] < 0] += 2 * np.pi
    
    latbins = np.deg2rad(np.arange(-91, 91.1, 2)) 
    lonbins = np.deg2rad(np.arange(-5, 365.1, 10))
    latcentre = (latbins[:-1] + latbins[1:]) / 2
    loncentre = (lonbins[:-1] + lonbins[1:]) / 2
    passes = {}
    for crd in crdtype:
        latind = find_closest(latcentre, vals['lat_' + crd])
        lonind = find_closest(loncentre, vals['lon_' + crd])
        combind = np.zeros(len(latind) + 1)
        combind[1:] += latind * 1E5 + lonind  # Leading zero ensures first block gets counted
        blockind = (combind[:-1] - combind[1:]) != 0
        lats = latcentre[latind[blockind].tolist()]
        lons = loncentre[lonind[blockind].tolist()]
        counts = np.histogram2d(lats, lons, np.array((latbins, lonbins)))
        passes[crd] = counts[0]

    passes['latbins'] = latbins
    passes['lonbins'] = lonbins
    return passes


def find_closest(A, target):
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx

if __name__ == '__main__':
    main()
