#!/usr/local/bin/python3
"""
proc_swarm_lp.py
Script to process the SWARM langmuir probe data and analyse for patches. 
"""
import pdb 
import spacepy.coordinates as crd 
from spacepy.time import Ticktock
import numpy as np
import scipy as sp
import datetime as dt
import matplotlib.pyplot as plt 
import matplotlib
import glob
import pickle
import sys, os
import collections
from proc_swarm_lp import load_lp
sys.path.append('../sounder/')
from plot_rtd import concat_files, calc_vht


def main():
    swarm, hf = load_data()
   
    sats = {'A': 'rx', 'B': 'gx'}
    for sat, clr in sats.items():
        # Figure out swarm chunks
        tot_sec = [(t - swarm[sat]['times'][0]).total_seconds() for t in swarm[sat]['times']]
        clumps = np.argwhere(np.diff(tot_sec) > 10)
        mof = np.zeros(clumps.shape)
        ne = np.zeros(clumps.shape)
        rmse = np.zeros(clumps.shape)

        start_ind = 0 
        for clumpind, end_ind in enumerate(clumps):
            tlim = swarm[sat]['times'][start_ind], swarm[sat]['times'][end_ind]
            # Calculate mean and RMSE Ne from Swarm
            swarm_tind = np.logical_and(swarm[sat]['times'] >= tlim[0], swarm[sat]['times'] <= tlim[1])
            swarm_ne = swarm[sat]['ne'][swarm_tind]
            ne[clumpind] = np.mean(swarm_ne)
            rmse[clumpind] = np.sqrt(np.mean(swarm_ne ** 2))

            # Calculate MOF
            for freq, spec in hf.items():
                tind = np.logical_and(spec['time'] > tlim[0], spec['time'] < tlim[1])
                if np.sum(tind) > 0:
                    if (spec['max_pwr_db'][tind].max() > 0) and freq > mof[clumpind]:
                        mof[clumpind] = freq
            start_ind = end_ind + 1
        from scipy.stats import pearsonr
        rs = pearsonr(mof[mof > 0], ne[mof > 0])
        print('Ne    Sat %s: R: %1.2f, P: %1.2f (<0.05 is significant)' % (sat, rs[0], rs[1]))
        rs = pearsonr(mof[mof > 0], rmse[mof > 0])
        print('RMS   Sat %s: R: %1.2f, P: %1.2f (<0.05 is significant)' % (sat, rs[0], rs[1]))
    
            
    plt.errorbar(mof[mof > 0], ne[mof > 0], yerr=rmse[mof > 0], fmt='.')
    plt.xlabel('Max. Obs. Freq. (MHz)')
    plt.ylabel('Mean & RMS Swarm Ne (el. m3)')
    plt.show()


def load_data( 
    swarmpath='./data/swarm/lp/',
    radiopath='/Users/chartat1/sounder/data/prc_analysis/no_badfrq/spectra/%Y%m%d',
    hf_fname_fmt='/Users/chartat1/sounder/data/prc_analysis/no_badfrq/daily/data/%Y%b%d_analysis.pkl',
    hf_out_fname_fmt=['./data/hf/hf_%Y%b%d_', '%Y%b%d.pkl'],
    swarm_out_fname_fmt=['./data/swarm_agg/swarm_%Y%b%d_', '%Y%b%d.pkl'],
    time=dt.datetime(2019, 2, 28),
    step=dt.timedelta(days=1),
    endtime=dt.datetime(2019, 3, 14), 
    sats=['A', 'B'],
    lat=-83.93,
    lon=166.69,
    lat_cutoff=3,
    lon_cutoff=5,
):

    hf = load_hf(time, step, endtime, radiopath, hf_fname_fmt, hf_out_fname_fmt) 
    swarm = load_swarm(time, swarmpath, swarm_out_fname_fmt, sats, lat, lon, lat_cutoff, lon_cutoff)
    return swarm, hf


def load_hf(time, step, endtime, radiopath, hf_fname_fmt, hf_out_fname_fmt): 
    # Load HF data
    hf_out_fname = time.strftime(hf_out_fname_fmt[0]) + endtime.strftime(hf_out_fname_fmt[1])
    try:
        with open(hf_out_fname, 'rb') as f:
            hf = pickle.load(f)    
    except:
        print('Could not load %s\n Recalculating...' % hf_out_fname)
        dirn = []
        tm = time 
        while tm <= endtime: 
            dirn.append(tm.strftime(radiopath))
            tm += dt.timedelta(days=1)
        hf = concat_files(dirn, hf_fname_fmt)
        for freq, vals in hf.items():
            for k, v in vals.items():
                hf[freq][k] = np.array(v)
            hf[freq]['alt'] = calc_vht(vals['range'])
        with open(hf_out_fname, 'wb') as f:
            pickle.dump(hf, f)
    return hf


def load_swarm(time, swarmpath, swarm_out_fname_fmt, sats, lat, lon, lat_cutoff, lon_cutoff):
    # Load Swarm
    swarm_out_fname = time.strftime(swarm_out_fname_fmt[0]) + endtime.strftime(swarm_out_fname_fmt[1])
    try:
        with open(swarm_out_fname, 'rb') as f:
            swarm = pickle.load(f)    
    except:
        swarm = {}
        for sat in sats:
            print('\nSatellite %s' % sat)
            swarm[sat] = {}
            tm = time
            while tm <= endtime: 
                print(tm.strftime('%Y/%b/%d'))
                fname_format = swarmpath + 'SW_*_EFI%s' % sat + '*%Y%m%d*.cdf' 
                fname_str = glob.glob(tm.strftime(fname_format))
                if len(fname_str) == 0:
                    print('No CDF file matching %s' % tm.strftime(fname_format))
                fname = fname_str[0]
                swarm_ts = load_lp(fname)
                latind = np.abs(swarm_ts['lat_geo'] - lat) < lat_cutoff
                lonind = np.abs(swarm_ts['lon_geo'] - lon) < lon_cutoff
                locind = np.logical_and(latind, lonind)
                for k, v in swarm_ts.items():
                    if k in swarm[sat].keys():
                        swarm[sat][k] = np.concatenate([swarm[sat][k], v[locind]])
                    else:
                        swarm[sat][k] = v[locind]
                tm += dt.timedelta(days=1)
        with open(swarm_out_fname, 'wb') as f:
            pickle.dump(swarm, f)
    return swarm

if __name__ == '__main__':
    main()


