#!/usr/local/bin/python3
"""
proc_swarm_bkgd_dens.py
Script to process the SWARM electron density data and get the average background value each day

"""
from spacepy import pycdf
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


def main(save=True, plot=True):
    time=dt.datetime(2014, 8, 1)
    endtime=dt.datetime(2017, 6, 29)

    if save:
        save_bkgd_dens(time=time, endtime=endtime)

    if plot:
        plot_bkgd_dens(time=time, endtime=endtime)


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
                  ):

    ipath = '/Volumes/Seagate/data/swarm/lp/'
    opath = '/Volumes/Seagate/data/swarm/proc_bkgd_dens/'
    
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
            bkgd_dens[sat] = calc_bkgd_dens(fname)

        fout = opath + time.strftime('bkgd_dens_%Y%m%d.pkl')
        with open(fout, 'wb') as f:
            pickle.dump(bkgd_dens, f) 
        print('Saving %s' % fout)

        time += dt.timedelta(days=1)


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






