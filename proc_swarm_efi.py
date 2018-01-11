"""
proc_swarm_efi.py
Script to perform basic processing of the Swarm EFI data to determine relationship between density patches and drifts
Information from AGU is that quality flags are important, validation with SuperDARN gives good correlation, maybe 2X bias
Questions to answer:
    1. What is the dependence of drift velocity on DOY in each of the polar caps
    2. What is the relationship between patch drift velocity and background polar cap drift velocity
        Background means:
            i) daily-averaged
            ii) pass-averaged
We only have 2016 data, and only ~4 passes/day
"""

import pdb 
import numpy as np
import scipy as sp
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt 
import matplotlib
import glob
import pickle
import sys 
import collections
import nvector as nv
from spacepy import pycdf
from plot_patch_ct import get_patch_ct
import seaborn as sns


def main(fname_format=['./data/swarm_efi/SW_EXPT_EFI%s', '_TIICT_%Y%m%d*.cdf'], \
   drift_fname_format='./data/drifts_%s_to_%s.pkl', \
                 sats=['A', 'B'], \
            starttime=dt.datetime(2016, 1, 1), 
              endtime=dt.datetime(2016, 12, 31), \
           lat_cutoff=70, \
        ):

    drift_fname = drift_fname_format % (starttime.strftime('%Y%m%d'), endtime.strftime('%Y%m%d'))
    try:
        efi = pd.read_pickle(drift_fname)
        print("Loading preprocessed drift file %s" % drift_fname)
    except:
        print("No preprocessed file found - loading drifts")
        efi = {}
        for sat in sats:
            efi[sat] = load_efi([fname_format[0] % sat, fname_format[1]], starttime, endtime)
        with open(drift_fname, 'wb') as f:
            pickle.dump(efi, f)  

    efi_df = {}
    for sat in sats:  
        efi_df[sat] = efi_to_dataframe(efi[sat], lat_cutoff)
    plot_drifts(efi_df)
 

def plot_drifts(efi):
    hemnames = {'nh': 'Northern', 'sh': 'Southern'}
    ct = 1
    for hem, hemname in hemnames.items():
        for sat, df in efi.items():
            plt.subplot(len(efi), len(hemnames), ct)
            hemind = df['latitude'] > 0 if hem == 'nh' else df['latitude'] < 0
            # df.loc[hemind, 'abs_viy'].plot(marker='.')
            plt.plot(df.loc[hemind, 'abs_viy'].resample('5min').median(), '.')
            plt.ylim([0, 2000])
            plt.grid()
            plt.title('Swarm %s %s hemisphere' % (sat, hemname))
            ct += 1
    plt.show()

"""
            hemind = efi[sat]['latitude'] > 0 if hem == 'nh' else efi[sat]['latitude'] < 0
            vels = efi[sat]['viy'][hemind]
            times = efi[sat]['times'][hemind]

            plt.plot_date(times, vels, '.k')

            ax = plt.gca()
            ax.set_xticks(ax.get_xticks()[::2])
            if ct <= len(sats):
                plt.title('Swarm %s %s hemisphere' % (sat, hemname))
            # plt.ylim([-1000, 1000])
            plt.grid()
            # plt.xlim(0, doymax)

            frame = plt.gca()
            if np.mod(ct, 2) == 0:
                frame.axes.yaxis.set_ticklabels([])
            else:
                plt.ylabel('Five-day median ion drift magnitude ($m s^{-1}$)')
            ct += 1

    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 15}
    matplotlib.rc('font', **font)

    plt.show()
"""


def load_efi(fname_format, starttime, endtime):
    """
    Variable definitions:
        'timestamp' seconds from 1 Jan 2000 00:00:00 UT
        'latitude', 'longitude' in degrees geographic
        'radius' in metres
        'viy' cross-track horizontal flow speed in satellite frame (m/s)
        'qy' quality of viy (2 is good)  
        'mlt' magnetic local time (in hours)
    """
    vars = 'timestamp', 'latitude', 'longitude', 'radius', 'viy', 'qy', 'mlt'
    efi = {v: [] for v in vars}
    time = starttime
    while time <= endtime:
        try:
            fin = pycdf.CDF(glob.glob(fname_format[0] + time.strftime(fname_format[1]))[0])
            data = {v: fin[v][...] for v in vars}
            for k, v in data.items():
                efi[k].append(v) 
        except:
            print(time.strftime('No file on %Y/%m/%d'))
        time += dt.timedelta(days=1)
   
    efi = {k: np.hstack(v) for k, v in efi.items()}
    return efi


def efi_to_dataframe(efi, lat_cutoff):
    good_data = efi['qy'] == 2
    above_lat_cutoff = np.abs(efi['latitude']) >= lat_cutoff
    sane_timestamp = (np.abs(efi['timestamp']) < 1E10) & (efi['timestamp'] > 10)
    efi = {k: v[good_data & above_lat_cutoff & sane_timestamp] for k, v in efi.items()}
    efi['times'] = np.array([dt.datetime(2000, 1, 1) + dt.timedelta(seconds=t) for t in efi['timestamp']])
    efi['abs_viy'] = np.abs(efi['viy'])
    efi_df = pd.DataFrame(efi).set_index('times')
    
    return efi_df


if __name__ == '__main__':
    main()
