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
import proc_swarm_lp


def main(efi_fname_fmt=['./data/swarm_efi/SW_EXPT_EFI%s', '_TIICT_%Y%m%d*.cdf'], \
             drift_fin='./data/drifts_%s_to_%s.pkl', \
          ne_fname_fmt=['./data/swarm_lp/SW_*EFI%s', '*%Y%m%d*.cdf'], \
                ne_fin='./data/ne_%s_to_%s.pkl', \
                ct_fin='/Volumes/Seagate/data/swarm/proc_lp/alex/55_deg/lp_%Y%m%d_55deg.pkl', \
                  sats=['A'], \
             starttime=dt.datetime(2016, 1, 1), 
               endtime=dt.datetime(2016, 1, 31), \
            lat_cutoff=70, \
         ):

    # Load electric fields
    drift_fname = drift_fin % (starttime.strftime('%Y%m%d'), endtime.strftime('%Y%m%d'))
    try:
        efi = pd.read_pickle(drift_fname)
        print("Loading preprocessed drift file %s" % drift_fname)
    except:
        print("No preprocessed file found - loading drifts")
        efi = {}
        for sat in sats:
            efi[sat] = load_efi([efi_fname_fmt[0] % sat, efi_fname_fmt[1]], starttime, endtime)
        with open(drift_fname, 'wb') as f:
            pickle.dump(efi, f)  

    efi_df = {}
    for sat in sats:  
        efi_df[sat] = efi_to_dataframe(efi[sat], lat_cutoff)

    # Load electron densities
    ne_fname = ne_fin % (starttime.strftime('%Y%m%d'), endtime.strftime('%Y%m%d'))
    try:
        ne = pd.read_pickle(ne_fname)
        print("Loading preprocessed ne file %s" % ne_fname)
    except:
        print("No preprocessed file found - loading ne")
        ne = {}
        for sat in sats:
            ne[sat] = load_ne([ne_fname_fmt[0] % sat, ne_fname_fmt[1]], starttime, endtime)
        with open(ne_fname, 'wb') as f:
            pickle.dump(ne, f)  

    ne_df = {}
    for sat in sats:  
        ne_df[sat] = ne_to_dataframe(ne[sat], lat_cutoff)

    # Load patch counts
    patch_ct = get_patch_ct(starttime, endtime, sats, ct_fin)

    # Plot density timeseries and EFI timeseries
    month = 1
    start_looktimes = [dt.datetime(2016, month, 1), dt.datetime(2016, month, 5), dt.datetime(2016, month, 10), \
                       dt.datetime(2016, month, 15), dt.datetime(2016, month, 20), dt.datetime(2016, month, 25)]
    for t in start_looktimes:
        plot_dens_efi(efi_df['A'], ne_df['A'], patch_ct['A'], start_looktime=t)

    # Plot histogram of EFI values
    if False:
        plot_drifts(efi_df)


def plot_dens_efi(efi, ne, patch_ct, start_looktime=dt.datetime(2016, 1, 1)):
    # find a patch with EFI data
    ct = 0
    good_patch = False
    try:
        while good_patch == False:
            time = patch_ct['times'][ct][0]
            window = dt.timedelta(minutes=5)
            starttime = time - window
            endtime = time + window
            efi_tind = (efi.index > starttime) & (efi.index < endtime)
            good_patch = (efi[efi_tind].index.min() < time) & (efi[efi_tind].index.max() > time) & (np.sum(efi_tind) > 0) & (time > start_looktime)
            ct += 1

        plt.subplot(3, 1, 1)
        ax = efi[efi_tind]['viy'].plot()
        ax.set_ylabel('X-track vel (m/s)')
        ax.set_title(time.strftime('%Y/%m/%d %H:%M'))
        ax.set_ylim(-4000, 4000)

        plt.subplot(3, 1, 2)
        ne_tind = (ne.index >= efi[efi_tind].index.min()) & (ne.index <= efi[efi_tind].index.max())
        ax2 = ne[ne_tind]['n'].plot()
        ax2.set_ylabel('Electron dens. (e-/cm3)')
        ax2.set_ylim(0, 500000)

        plt.subplot(3, 1, 3)
        ne_tind = (ne.index >= efi[efi_tind].index.min()) & (ne.index <= efi[efi_tind].index.max())
        ax3 = ne[ne_tind]['T_elec'].plot()
        ax3.set_ylabel('Electron Temp. (K)')
        ax3.set_ylim(0, 6500)
        
        plt.show()
    except:
        None


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


def load_ne(fname_format, starttime, endtime):
    """
    Variable definitions:
        'timestamp' seconds from 1 Jan 2000 00:00:00 UT
        'latitude', 'longitude' in degrees geographic
        'radius' in metres
        'mlt' magnetic local time (in hours)
    """
    vars = 'Timestamp', 'Latitude', 'Longitude', 'n', 'T_elec'
    ne = {v: [] for v in vars}
    time = starttime
    while time <= endtime:
        try:
            fin = pycdf.CDF(glob.glob(fname_format[0] + time.strftime(fname_format[1]))[0])
            data = {v: fin[v][...] for v in vars}
            for k, v in data.items():
                ne[k].append(v) 
        except:
            print(time.strftime('No file on %Y/%m/%d'))
        time += dt.timedelta(days=1)
   
    ne = {k: np.hstack(v) for k, v in ne.items()}
    return ne


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


def ne_to_dataframe(ne, lat_cutoff):
    above_lat_cutoff = np.abs(ne['Latitude']) >= lat_cutoff
    ne = {k: v[above_lat_cutoff] for k, v in ne.items()}
    ne_df = pd.DataFrame(ne).set_index('Timestamp')
    
    return ne_df


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
