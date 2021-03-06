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
import plot_patch_ct

def main(ipath='/Volumes/Seagate/data/swarm/lp/',
         opath='/Volumes/Seagate/data/swarm/pass_ct/', 
         time=dt.datetime(2014, 1, 1), 
         endtime=dt.datetime(2017, 7, 1), 
         sats=['A', 'B', 'C'],
         lat_cutoff=70,
         save=True):

    while time <= endtime:
        timestr = time.strftime('%Y-%m-%d')
        print(timestr)
        vals = {}
        pass_ct = {}
        pass_norm = {}

        for sat in sats:
            print('\nSatellite %s' % sat)
            fname_format = ipath + 'SW_*EFI%s' % sat + '*%Y%m%d*.cdf'
            try:
                fname = glob.glob(time.strftime(fname_format))[0]
                vals[sat] = proc_swarm_lp.load_lp(fname)
                # pass_ct[sat] = count_passes(vals[sat])
                pass_norm[sat] = norm_passes(vals[sat], lat_cutoff=lat_cutoff)  
            except:
                print('Could not count passes for satellite %s on %s' % (sat, timestr))

        if save:
            # fout = opath + time.strftime('/pass_%Y%m%d.pkl')
            # with open(fout, 'wb') as f:
            #     pickle.dump(pass_ct, f)
            # print('Saving %s' % fout)
            fnorm_out = opath + time.strftime('/pass_norm_%Y%m%d') + '_%ideg.pkl' % lat_cutoff
            with open(fnorm_out, 'wb') as f:
                pickle.dump(pass_norm, f)
            print('Saving %s' % fnorm_out)
        time += dt.timedelta(days=1)


def get_ct(ipath='./data/pass_ct/pass_%Y%m%d.pkl', 
                 time=dt.datetime(2016, 1, 1), 
              endtime=dt.datetime(2017, 1, 1), 
                 sats=['A', 'B', 'C'], 
                  crd='mag'):
    pass_count = {}
    for sat in sats:
        pass_count[sat] = {}
    while time <= endtime:
        fin = time.strftime(ipath)
        with open(fin, 'rb') as f:
            pass_ct = pickle.load(f)
        for key, val in pass_ct.items():
            for k in ['times', 'hem']:
                try:
                    pass_count[key][k] = np.append(pass_count[key][k], val[crd][k])
                except:
                    pdb.set_trace()
                    pass_count[key][k] = val[crd][k]
        time += dt.timedelta(days=1)
    return pass_count


def get_norm_ct(ipath='./data/pass_ct/pass_norm_%Y%m%d.pkl', 
            starttime=dt.datetime(2016, 1, 1), 
              endtime=dt.datetime(2017, 1, 1), 
                 sats=['A', 'B', 'C']):
    pass_count = {}
    for sat in sats:
        pass_count[sat] = {}

    n_months = 12
    n_hours = 24
    norm_2dvars = 'ut', 'lt', 'mlt' 
    hems = 'nh_', 'sh_'

    for sat in sats:
        for hem in hems:
            for v in norm_2dvars:
                pass_count[sat][hem + v + '_2d'] = np.zeros((n_months, n_hours))

    t = starttime
    while t <= endtime:
        fin = t.strftime(ipath)
        with open(fin, 'rb') as f:
            pass_ct = pickle.load(f)
        for sat in sats:
            # make a 2D array of hour vs cal. month
            try:
                for key, val in pass_ct[sat].items():
                    if t == starttime:
                        pass_count[sat][key] = np.array(val)
                    else:
                        pass_count[sat][key][0] += val[0]
                month = t.month
                for hem in hems:
                    for var in norm_2dvars:
                        pass_count[sat][hem + var + '_2d'][month - 1, :] += pass_ct[sat][hem + var][0]

            except:
                print('No entry for satellite %s on %s' % (sat, t))
        t += dt.timedelta(days=1)
    return pass_count


def norm_passes(vals, lat_cutoff=70,
                ):
    # Transform lats/lons to magnetic 
    vals['lat_geo'] *= np.pi / 180 
    vals['lon_geo'][vals['lon_geo'] < 0] += 360 
    vals['lon_geo'] *= np.pi / 180 
    alts, vals['lat_mag'], vals['lon_mag'] = proc_swarm_lp.transform(vals['rad'], vals['lat_geo'], \
                          vals['lon_geo'], from_=['GEO', 'sph'], to=['MAG', 'sph'])
    vals['lon_mag'][vals['lon_mag'] < 0] += 2 * np.pi
    vals['lon_geo'][vals['lon_geo'] < 0] += 2 * np.pi
    ut = np.array([t.hour + t.minute / 60 for t in vals['times']]) 
    mlt = plot_patch_ct.calc_mlt(ut, np.rad2deg(vals['lon_mag']))
    lt = plot_patch_ct.calc_mlt(ut, np.rad2deg(vals['lon_geo']))

    latind = np.abs(vals['lat_mag']) >= np.deg2rad(lat_cutoff) # Filter out low latitudes
    nhind = vals['lat_mag'] > 0
    shind = vals['lat_mag'] < 0

    norm_params = {}
    norm_params['hem'] = np.histogram(vals['lat_mag'][latind], [-6, 0, 6])
    norm_params['nh_mlt'] = np.histogram(mlt[np.logical_and(latind, nhind)], np.arange(0, 24.1, 1))
    norm_params['sh_mlt'] = np.histogram(mlt[np.logical_and(latind, shind)], np.arange(0, 24.1, 1))
    norm_params['nh_lt'] = np.histogram(lt[np.logical_and(latind, nhind)], np.arange(0, 24.1, 1))
    norm_params['sh_lt'] = np.histogram(lt[np.logical_and(latind, shind)], np.arange(0, 24.1, 1))
    norm_params['nh_ut'] = np.histogram(ut[np.logical_and(latind, nhind)], np.arange(0, 24.1, 1))
    norm_params['sh_ut'] = np.histogram(ut[np.logical_and(latind, shind)], np.arange(0, 24.1, 1))
    norm_params['nh_mlon'] = np.histogram(vals['lon_mag'][np.logical_and(latind, nhind)], np.deg2rad(np.arange(-5, 365.1, 10)))
    norm_params['sh_mlon'] = np.histogram(vals['lon_mag'][np.logical_and(latind, shind)], np.deg2rad(np.arange(-5, 365.1, 10)))

    return norm_params 


def count_passes(vals, crdtype='mag'):
    # Transform lats/lons to magnetic 
    pdb.set_trace()  # NOTE: Not currently working. Look for an earlier commit to use
    vals['lat_geo'] *= np.pi / 180
    vals['lon_geo'][vals['lon_geo'] < 0] += 360
    vals['lon_geo'] *= np.pi / 180
    alts, vals['lat_mag'], vals['lon_mag'] = proc_swarm_lp.transform(vals['rad'], vals['lat_geo'], \
                          vals['lon_geo'], from_=['GEO', 'sph'], to=['MAG', 'sph'])
    vals['lon_mag'][vals['lon_mag'] < 0] += 2 * np.pi
    
    passes = {}
    utbins = np.arange(0, 24.1)
    hems = {'nh': vals['lat_mag'] >= np.deg2rad(70),
            'sh': vals['lat_mag'] <= np.deg2rad(-70)}
    vals['ut'] = np.array([t.hour + t.minute / 60 for t in vals['times']])
    for hem, hemind in hems.items():
        pdb.set_trace()
        passes[hem] = np.histogram(vals['ut'][hemind], binedges=utbins) 

    return passes


def find_closest(A, target):
    # A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx

if __name__ == '__main__':
    main()
