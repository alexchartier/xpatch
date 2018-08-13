#!/usr/local/bin/python3
"""
proc_swarm_bkgd_dens.py
Script to process the SWARM electron density data and get the average background value each day

"""
import pdb
import numpy as np
import scipy as sp
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib
import glob
import pickle
import sys 
import collections
sys.path.insert(0, '/users/chartat1/fusionpp/fusion/')
import socket
import proc_swarm_lp
from scipy import stats


def main(
        save=False, 
        plot=True,
        time=dt.datetime(2014, 8, 1),
        endtime=dt.datetime(2018, 8, 1),
        datapath='data/swarm/lp/',
        out_type='lat_bin',
        procpath=['data/swarm/proc_latbin_dens/%s', '_dens_%Y%m%d.pkl'],
        ):
    procpath = procpath[0] % out_type + procpath[1]
    lat_bins = np.arange(-90, 91, 5)
    if save:
        save_bkgd_dens(ipath=datapath, opath=procpath, 
                time=time, endtime=endtime, lat_bins=lat_bins, out_type=out_type)

    if plot:
        if out_type == 'lat_bin':
            plot_avg_dens(procpath, time, endtime, lat_bins, sats=['A', 'B'])
        else:
            plot_bkgd_dens(ipath=procpath, time=time, endtime=endtime)


def plot_avg_dens(
                  ipath, starttime, endtime, lat_bins, sats=['A', 'B'], 
                ):
    # Calculate timeseries 
    times = []
    time = starttime
    while time <= endtime:
        times.append(time)
        time += dt.timedelta(days=1)

    # set up holder
    lat_c = (lat_bins[:-1] + lat_bins[1:]) / 2
    dens_ts = {}
    for sat in sats: 
       dens_ts[sat] = np.zeros((len(times), len(lat_bins) - 1)) * np.nan
   
    # Load files 
    for tind, time in enumerate(times):
        print(time)
        with open(time.strftime(ipath), 'rb') as f: 
            bkgd_dens = pickle.load(f)
        for sat in sats: 
            try:
                dens_ts[sat][tind, :] = bkgd_dens[sat][0]
            except:
                print('Failed to load %s' % sat)
    
    for satind, sat in enumerate(sats):
        plt.subplot(len(sats), 1, satind + 1)
        pdb.set_trace()
        plt.pcolor(times, lat_c, dens_ts[sat].T)
        plt.clim([0, 2E6])
    plt.show()
    


def plot_bkgd_dens(
                  ipath='/Volumes/Seagate/data/swarm/proc_bkgd_dens/bkgd_dens_%Y%m%d.pkl',
                  time=dt.datetime(2016, 1, 1),
                  step=dt.timedelta(days=1),
                  endtime=dt.datetime(2016, 12, 31),
                  sats=['A', 'B'],
                  hems=['nh', 'sh'],
                  ):
 
    # set up holder
    dens_ts = {}
    for sat in sats: 
       dens_ts[sat] = {}
       for hem in hems: 
           dens_ts[sat][hem] = []
   
    # Load files 
    t = time
    times = []
    while t <= endtime:
        with open(t.strftime(ipath), 'rb') as f: 
            bkgd_dens = pickle.load(f)
        for sat in sats: 
            for hem in hems: 
                try:
                    dens_ts[sat][hem].append(bkgd_dens[sat][hem])
                except:
                    dens_ts[sat][hem].append(np.nan)
        times.append(t)
        t += step
 
    hemnames = {'nh': 'Northern', 'sh': 'Southern'} 
    ct = 1 
    for hem in hems: 
        plt.subplot(1, len(hems), ct)
        # Get avg. dens across sats
        dens = np.zeros(len(times))
        for sat in sats: 
            dens += np.array(dens_ts[sat][hem])
        dens /= len(sats) 
        dens_med = sp.signal.medfilt(dens, 5)
        plt.plot_date(times, dens_med, '.k')

        ax = plt.gca()
        ax.set_xticks(ax.get_xticks()[::2])
        plt.title('%s hemisphere' % hemnames[hem])
        ymax = 250000 
        plt.ylim([0, ymax])
        plt.grid()
        # plt.xlim(0, doymax)

        yr = time.year
        dec_sols = []
        jun_sols = []

        while yr < endtime.year:
            dec_sols.append(dt.datetime(yr, 12, 21))
            jun_sols.append(dt.datetime(yr, 6, 21))
            yr += 1

        cnt = 1
        for d in jun_sols:
            if cnt == 1:
                plt.plot_date([d, d], [0, ymax], 'b--', label='June Solstice')
                cnt += 1
            else:
                plt.plot_date([d, d], [0, ymax], 'b--')

        for d in dec_sols:
            if cnt == 2:
                plt.plot_date([d, d], [0, ymax], 'r--', label='December Solstice')
                cnt += 1
            else:
                plt.plot_date([d, d], [0, ymax], 'r--')

        if ct == 2:
            plt.legend()

        frame = plt.gca()
        if np.mod(ct, 2) == 0:
            frame.axes.yaxis.set_ticklabels([])
        else:
            plt.ylabel('Five-day median electon density ($cm^{-3}$)')
        ct += 1
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 15}
    matplotlib.rc('font', **font) 

    plt.show()
 
 
def save_bkgd_dens(
                  time=dt.datetime(2016, 1, 1),
                  step=dt.timedelta(days=1),
                  endtime=dt.datetime(2016, 2, 1),
                  sats=['A', 'B'],
                  ipath='/Volumes/Seagate/data/swarm/lp/',
                  opath='/Volumes/Seagate/data/swarm/proc_bkgd_dens/',
                  out_type='lat_bin',
                  lat_bins=None,
                  ):
    
    while time <= endtime: 
        timestr = time.strftime('%Y-%m-%d')
        print(timestr)
        vals = {}
        bkgd_dens = {}

        for sat in sats:
            print('Satellite %s' % sat)
            fname_format = ipath + 'SW_*_EFI%s' % sat + '*%Y%m%d*.cdf'
            try:
                fname = glob.glob(time.strftime(fname_format))[0]
            except:
                print('No file for satellite %s on %s' % (sat, timestr))
                continue
        
            if out_type == 'bkgd':
                bkgd_dens[sat] = calc_bkgd_dens(fname)
            elif out_type == 'lat_bin':
                try:
                    bkgd_dens[sat] = calc_avg_dens(fname, lat_bins, lt_lim=[12, 18])
                except:
                    pdb.set_trace()
                    print('Could not process satellite %s on %s' % (sat, timestr))
                    bkgd_dens[sat] = np.ones((2, len(lat_bins) - 1)) * np.nan

        fout = time.strftime(opath)
        with open(fout, 'wb') as f:
            pickle.dump(bkgd_dens, f) 
        print('Saving %s' % fout)

        time += dt.timedelta(days=1)


def calc_avg_dens(fname, latbins, lt_lim=[12, 18]):
    # Calculate the latitude-binned mean density  each day
    vals = get_swarm_vals(fname)
    try:
        vals['Diplat'] = vals['lat_mag']
    except:
        None
    ut = np.array([(t.hour + t.minute / 60) for t in vals['times']])
    lt = ut + vals['lon_geo'] / 360 * 24 
    lt[lt < 0] += 24
    lt[lt >= 24] -= 24
    lt_ind = (lt > lt_lim[0]) & (lt < lt_lim[1])
    bin_avg = stats.binned_statistic(vals['Diplat'][lt_ind], vals['ne'][lt_ind], bins=latbins)
    return bin_avg


def calc_bkgd_dens(fname, lat_cutoff=55):
    """
    Calculates the average electron density in each SWARM file

    Inputs: 
        fname = '~/Downloads/****.cdf'
        lat_cutoff  # degrees magnetic
        elev_cutoff  # degrees
    
    Returns: 
       bkgd_dens - daily average background dens within the limits specified
    """ 
    vals = get_swarm_vals(fname)

    # Take out values below latitude cutoff
    index = np.abs(vals['lat_mag']) > lat_cutoff
    for key, val in vals.items():  
        vals[key] = val[index]

    # Remove flagged values
    pdb.set_trace()

    bkgd_dens = {}
    bkgd_dens['nh'] = np.mean(vals['ne'][vals['lat_geo'] > 0])
    bkgd_dens['sh'] = np.mean(vals['ne'][vals['lat_geo'] < 0])

    return bkgd_dens 


def get_swarm_vals(fname):
    vals = proc_swarm_lp.load_lp(fname)
    # Transform lats/lons to magnetic 
    alts, vals['lat_mag'], vals['lon_mag'] = proc_swarm_lp.transform(vals['rad'], np.deg2rad(vals['lat_geo']), \
                          np.deg2rad(vals['lon_geo']), from_=['GEO', 'sph'], to=['MAG', 'sph'])
    vals['lat_mag'] *= 180 / np.pi
    vals['lon_mag'] *= 180 / np.pi
    return vals


if __name__ == '__main__':
    main()






