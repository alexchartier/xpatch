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
from spacepy import pycdf
from plot_patch_ct import get_patch_ct
import proc_swarm_lp
import physics


def main(efi_fname_fmt=['./data/swarm_efi/SW_EXPT_EFI%s', '_TIICT_%Y%m%d*.cdf'], \
             drift_fin='./data/drifts_%s_to_%s.pkl', \
          ne_fname_fmt=['./data/swarm_lp/SW_*EFI%s', '*%Y%m%d*.cdf'], \
                ne_fin='./data/ne_%s_to_%s.pkl', \
                ct_fin='./data/swarm/proc_lp/alex/lp_%Y%m%d_55deg.pkl', \
                  sats=['A', 'B'], \
             starttime=dt.datetime(2016, 1, 1), 
               endtime=dt.datetime(2016, 2, 1), \
            lat_cutoff=0, \
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
    # patch_ct = get_patch_ct(starttime, endtime, sats, ct_fin)


    # Spaghetti plots
    plot_spaghetti(efi_df['A'], ne_df['A'])

    # Plot histogram of patch quantities vs. average quantities
    if False:
        plot_hists(efi_df['A'], ne_df['A'], patch_ct['A'])

    if False:
        # Plot density timeseries and EFI timeseries
        month = 1
        start_looktimes = [dt.datetime(2016, month, 1), dt.datetime(2016, month, 5), dt.datetime(2016, month, 10), \
                           dt.datetime(2016, month, 15), dt.datetime(2016, month, 20), dt.datetime(2016, month, 25)]
        for t in start_looktimes:
            plot_dens_efi(efi_df['A'], ne_df['A'], patch_ct['A'], start_looktime=t)

    if False:
        # Plot histogram of EFI values
        plot_drifts(efi_df)


def plot_spaghetti(efi_df, ne_df):
    xvar = 'lat_mag'
    ct = 0
    nrows = 2
    ncols = 2
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols)
    pltvars = {'n':[0, 1500000], 'T_elec':[0, 14000]}
    for k, v in pltvars.items():
        c = 'k' if k == 'n' else 'r'
        ne_df[ne_df['lat_mag'] > 0].plot(x=xvar, y=k, marker='.', markersize=0.1, c=c, legend=False, \
                                        linewidth=0, ax=axes[ct, 0], ylim=v,)# xlim=[0, 24], )
        ct += 1
    ct = 0
    for k, v in pltvars.items():
        c = 'k' if k == 'n' else 'r'
        ne_df[ne_df['lat_mag'] < 0].plot(x=xvar, y=k, marker='.', markersize=0.1, c=c, legend=False, \
                                        linewidth=0, ax=axes[ct, 1], ylim=v)# xlim=[0, 24], )
        ct += 1
    """
    efi_df['B'] = np.sqrt(efi_df['bx'] ** 2 + efi_df['by'] ** 2 + efi_df['bz'] ** 2)
    efi_vars = {'abs_viy':[0, 5000], 'B':[30000, 55000]}
    for k, v in efi_vars.items():
        efi_df[efi_df['lat_mag'] > 0].plot(x=xvar, y=k, marker='.', markersize=0.1, legend=False, \
                                            linewidth=0, ax=axes[ct, 0], ylim=v)#, xlim=[0, 24])
        efi_df[efi_df['lat_mag'] < 0].plot(x=xvar, y=k, marker='.', markersize=0.1, legend=False, \
                                            linewidth=0, ax=axes[ct, 1], ylim=v)#, xlim=[0, 24])
        if ct == 2:
            ct += 1
    
    """
    plt.suptitle(ne_df.index[0].strftime('%b %Y'))
    ylabs = 'electrons / cm3', 'elec. temp (K)', 'abs. vel. (m/s)', 'B (nT)'
    hem = 'North', 'South'
    for r in np.arange(nrows):
        for c in np.arange(ncols):
            if r < nrows - 1:
                axes[r, c].tick_params(labelbottom='off')
                axes[r, c].set_xlabel('')
            if r == 0:
                axes[r, c].set_title(hem[c])
            if c > 0:
                axes[r, c].tick_params(labelleft='off')
            else:
                axes[r, c].set_ylabel(ylabs[r])
            axes[r, c].grid()
    plt.show()


def plot_hists(efi_df, ne_df, patch_ct):

    # Get the times of patches in each file
    looktime = patch_ct['times'][0][0]
    ne_times = []
    efi_ts = []
    while looktime < patch_ct['times'][-1][0]:
        try:
            ne_time, efi_t, efi_tind_window = find_good_patch(efi_df, patch_ct, start_looktime=looktime) 
            ne_times.append(ne_time)
            efi_ts.append(efi_t)
            looktime = ne_time + dt.timedelta(seconds=15)
        except:
            looktime = patch_ct['times'][-1][0] 
    ne_times = np.array(ne_times)

    patch_ct.pop('params')
    # Store variables
    efi_vars = {} 
    ne_vars = {} 
    for k in efi_df.keys():
        efi_vars[k] = []
    for k in patch_ct.keys():
        ne_vars[k] = []
    for tind, t in np.ndenumerate(ne_times):
        ne_tind = patch_ct['times'] == np.array(t)
        for k in ne_vars.keys():
            ne_vars[k].append(np.array(patch_ct[k])[np.squeeze(ne_tind)])
        for k in efi_vars.keys():
            efi_vars[k].append(efi_df.loc[efi_ts[tind[0]], [k]])
    for k in ne_vars.keys():
        ne_vars[k] = np.array(ne_vars[k])
    for k in efi_vars.keys():
        efi_vars[k] = np.array(efi_vars[k])
   
    efi_plots = 'viy', 'abs_viy', 'mlt',  
    ne_plots = 'ne_rm', 'ne_bg', 'T_elec', 'lat_mag'
    for p in np.arange(len(efi_plots)):
        plt.subplot(len(efi_plots), 1, p + 1)
        plt.hist(efi_vars[efi_plots[p]])
        plt.title(efi_plots[p])
    plt.show()
    for p in np.arange(len(ne_plots)):
        plt.subplot(len(ne_plots), 1, p + 1)
        plt.hist(np.squeeze(ne_vars[ne_plots[p]]))
        plt.title(ne_plots[p])
    plt.show()
        

def find_good_patch(efi, patch_ct, start_looktime=dt.datetime(2016, 1, 1)):
    # find a patch with EFI data
    ct = 0
    good_patch = False
    try:
        while good_patch == False:
            time = patch_ct['times'][ct][0]
            window = dt.timedelta(minutes=5)
            starttime = time - window
            endtime = time + window
            efi_tind_window = (efi.index > starttime) & (efi.index < endtime)
            good_patch = (efi[efi_tind_window].index.min() < time) & (efi[efi_tind_window].index.max() > time) \
                       & (np.sum(efi_tind_window) > 0) & (time > start_looktime)
            ct += 1
    except:
        None
    if good_patch != False:
        efi_tind = np.abs(efi[efi_tind_window].index - time) == np.abs(efi[efi_tind_window].index - time).min()
        efi_t = efi[efi_tind_window][efi_tind].index._data[0]
        return time, efi_t, efi_tind_window
    else:
        return None


def plot_dens_efi(efi, ne, patch_ct, start_looktime=dt.datetime(2016, 1, 1)):
    # Plot electron density and EFI parameters 
    time, efi_tind, efi_tind_window = find_good_patch(efi, patch_ct, start_looktime=dt.datetime(2016, 1, 1))
    plt.subplot(3, 1, 1)
    ax = efi[efi_tind_window]['viy'].plot()
    ax.set_ylabel('X-track vel (m/s)')
    ax.set_title(time.strftime('%Y/%m/%d %H:%M'))
    ax.set_ylim(-4000, 4000)

    plt.subplot(3, 1, 2)
    ne_tind = (ne.index >= efi[efi_tind_window].index.min()) & (ne.index <= efi[efi_tind_window].index.max())
    ax2 = ne[ne_tind]['n'].plot()
    ax2.set_ylabel('Electron dens. (e-/cm3)')
    ax2.set_ylim(0, 500000)

    plt.subplot(3, 1, 3)
    ne_tind = (ne.index >= efi[efi_tind_window].index.min()) & (ne.index <= efi[efi_tind_window].index.max())
    ax3 = ne[ne_tind]['T_elec'].plot()
    ax3.set_ylabel('Electron Temp. (K)')
    ax3.set_ylim(0, 6500)
    
    plt.show()


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
    vars = 'Timestamp', 'Latitude', 'Longitude', 'n', 'T_elec', 'Flags_LP', 'Flags_LP_T_elec', 'Flags_LP_n'
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
    vars = 'timestamp', 'latitude', 'longitude', 'radius', 'qdlat', 'viy', 'qy', 'mlt', 'bx', 'by', 'bz'
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
    try:
        ne['Latitude'] = ne['lat_geo']
    except:
        None
    ne['lat_mag'], ne['lon_mag'], ne['mlt'] = calc_mlt(ne['Latitude'], ne['Longitude'], ne['Timestamp'])
    above_lat_cutoff = np.abs(ne['lat_mag']) >= lat_cutoff
    not_crazy = (ne['n'] < 1E8) & (ne['T_elec'] < 1E6)
    good_flags = (ne['Flags_LP_n'] == 20) & (ne['Flags_LP_T_elec'] == 20)  # & (ne['Flags_LP'] == 1) 
    ne = {k: v[above_lat_cutoff & not_crazy & good_flags] for k, v in ne.items()}
    ne_df = pd.DataFrame(ne).set_index('Timestamp')
    return ne_df


def efi_to_dataframe(efi, lat_cutoff):
    efi['lat_mag'] = efi['qdlat']
    good_data = efi['qy'] == 2
    above_lat_cutoff = np.abs(efi['lat_mag']) >= lat_cutoff
    sane_timestamp = (np.abs(efi['timestamp']) < 1E10) & (efi['timestamp'] > 10)
    sane_mlt = efi['mlt'] <= 24
    efi = {k: v[good_data & above_lat_cutoff & sane_timestamp & sane_mlt] for k, v in efi.items()}
    efi['times'] = np.array([dt.datetime(2000, 1, 1) + dt.timedelta(seconds=t) for t in efi['timestamp']])
    efi['abs_viy'] = np.abs(efi['viy'])
    efi_df = pd.DataFrame(efi).set_index('times')
    return efi_df
    

def find_closest(alist, target):
    return min(alist, key=lambda x:abs(x-target))


def list_matching(list1, list2):
    list1_copy = list1[:]
    pairs = []
    for i, e in enumerate(list2):
        elem = find_closest(list1_copy, e)
        pairs.append([i, list1.index(elem)])
        list1_copy.remove(elem)
    return pairs


def calc_mlt(lat, lon, time):
    alts, mlat, mlon = physics.transform(lat * 0 + 1, \
                                         np.deg2rad(lat), \
                                         np.deg2rad(lon), \
                                         from_=['GEO', 'sph'], to=['MAG', 'sph']) 
    mlat, mlon = np.rad2deg(mlat), np.rad2deg(mlon)
    dh = np.array([t.hour + t.minute / 60 for t in time])
    mlt = mlon / 360 * 24 + dh
    mlt[mlt < 0] += 24
    mlt[mlt >= 24] -= 24
    return mlat, mlon, mlt



if __name__ == '__main__':
    main()
