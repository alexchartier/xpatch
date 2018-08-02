#!/usr/local/bin/python3
"""
proc_swarm_bkgd_tec.py
Script to process the SWARM upward-looking TEC data and get the average background value each day

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
import proc_swarm_lp
import socket


def main(save=False, plot=True):
    time=dt.datetime(2014, 8, 1)
    endtime=dt.datetime(2017, 6, 29)

    if save:
        save_bkgd_tec(time=time, endtime=endtime)

    if plot:
        plot_bkgd_tec(time=time, endtime=endtime)


def plot_bkgd_tec(
                  time=dt.datetime(2016, 1, 1),
                  step=dt.timedelta(days=1),
                  endtime=dt.datetime(2016, 12, 31),
                  sats=['A', 'B'],
                  hems=['nh', 'sh'],
                 ):
 
    if socket.gethostname() == 'chartat1-ml2':
        ipath = '/Volumes/Seagate/data/swarm/proc_bkgd_tec/bkgd_tec_%Y%m%d.pkl'
    elif socket.gethostname() == 'chartat1-ml1':
        ipath = 'data/proc_bkgd_tec/bkgd_tec_%Y%m%d.pkl'   

    # set up holder
    tec_ts = {}
    for sat in sats: 
       tec_ts[sat] = {}
       for hem in hems: 
           tec_ts[sat][hem] = []
    
    t = time
    times = []
    while t <= endtime:
        with open(t.strftime(ipath), 'rb') as f: 
            bkgd_tec = pickle.load(f)
        for sat in sats: 
            for hem in hems: 
                tec_ts[sat][hem].append(bkgd_tec[sat][hem])
        times.append(t)
        t += step
  
    ct = 1 
    for hem in hems: 
        for sat in sats: 
            plt.subplot(len(sats), len(hems), ct)
            tec_ts[sat][hem] = np.array(tec_ts[sat][hem])
            
            doy = [(t - time).days for t in times]
            plt.plot(doy, tec_ts[sat][hem], 'k')
    
            ax = plt.gca()
            ax.set_xticks(ax.get_xticks()[::2])
            plt.title('sat %s %s hemisphere' % (sat, hem))
            ymax = 18 
            plt.ylim([0, ymax])
            plt.grid()
            doymin = min(doy)
            doymax = max(doy)
            plt.ylim(0, ymax)
            plt.xlim(0, doymax)


            yr = time.year
            dec_sols = []
            jun_sols = []

            while yr < endtime.year:
                dec_sols.append((dt.datetime(yr, 12, 21) - time).days)
                jun_sols.append((dt.datetime(yr, 6, 21) - time).days)
                yr += 1

            cnt = 1
            for d in jun_sols:
                if cnt == 1:
                    plt.plot([d, d], [0, doymax], 'b--', label='June Solstice')
                    cnt += 1
                else:
                    plt.plot([d, d], [0, doymax], 'b--')

            for d in dec_sols:
                if cnt == 2:
                    plt.plot([d, d], [0, doymax], 'r--', label='December Solstice')
                    cnt += 1
                else:
                    plt.plot([d, d], [0, doymax], 'r--')

            if ct == 2:
                plt.legend()

            frame = plt.gca()
            if np.mod(ct, 2) == 0:
                frame.axes.yaxis.set_ticklabels([])
            else:
                plt.ylabel('average slant TEC (TECU)')
            if ct < 3:
                frame.axes.xaxis.set_ticklabels([])
            else:
                plt.xlabel(time.strftime('Days after %d %b %Y'))
            ct += 1
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 15}
    matplotlib.rc('font', **font) 

    plt.show()
 
 
def save_bkgd_tec(
                  time=dt.datetime(2016, 1, 1),
                  step=dt.timedelta(days=1),
                  endtime=dt.datetime(2016, 2, 1),
                  sats=['A', 'B'],
                  ):

    if socket.gethostname() == 'chartat1-ml2':
        ipath = '/Volumes/Seagate/data/swarm/gps_tec/'
        opath = '/Volumes/Seagate/data/swarm/proc_bkgd_tec/'
    elif socket.gethostname() == 'chartat1-ml1':
        ipath = 'data/swarm_tec/'
        opath = 'data/proc_bkgd_tec/'
    
    while time <= endtime: 
        timestr = time.strftime('%Y-%m-%d')
        print(timestr)
        vals = {}
        bkgd_tec = {}

        for sat in sats:
            print('Satellite %s' % sat)
            fname_format = time.strftime(ipath) + 'SW_OPER_TEC%s' % sat + '*%Y%m%d*.cdf'
            try:
                fname = glob.glob(time.strftime(fname_format))[0]
            except:
                print('No file for satellite %s on %s' % (sat, timestr))
                continue
            bkgd_tec[sat] = calc_bkgd_tec(fname)

        fout = opath + time.strftime('bkgd_tec_%Y%m%d.pkl')
        with open(fout, 'wb') as f:
            pickle.dump(bkgd_tec, f) 
        print('Saving %s' % fout)

        time += dt.timedelta(days=1)


def calc_bkgd_tec(fname, lat_cutoff=55, elev_cutoff=25):
    """
    Calculates the average VTEC in each SWARM file

    Inputs: 
        fname = '~/Downloads/SW_OPER_TECATMS_2F_20131125T105953_20131125T235953_0201.DBL'
        lat_cutoff  # degrees magnetic
        elev_cutoff  # degrees
    
    Returns: 
       bkgd_tec - daily average background tec within the limits specified
    """ 
    vals, vars = get_swarm_vals(fname)

    pdb.set_trace()
    # Take out values below elevation and latitude cutoffs 
    index = np.logical_and(vals['elev'] >= elev_cutoff, np.abs(vals['lat_mag']) > lat_cutoff)
    for key, val in vars.items():  
        vals[val] = vals[val][index]

    bkgd_tec = {}
    bkgd_tec['nh'] = np.mean(vals['tec'][vals['lat_geo'] > 0])
    bkgd_tec['sh'] = np.mean(vals['tec'][vals['lat_geo'] < 0])

    return bkgd_tec 


def get_swarm_vals(fname):
    vals, vars = load_swarm(fname)
    # Preliminary calculations
    rad = np.sqrt(np.sum(vals['leo_pos'] ** 2, axis=1))
    unused_alts, vals['lat_mag'], vals['lon_mag'] = proc_swarm_lp.transform(rad, vals['lat_geo'] * np.pi / 180, \
                          vals['lon_geo'] * np.pi / 180, from_=['GEO', 'sph'], to=['MAG', 'sph'])
    vals['lat_mag'] *= 180 / np.pi
    vals['lon_mag'] *= 180 / np.pi
    vals['elev'] = proc_swarm_lp.elevation(vals['gps_pos'], vals['leo_pos']) * 180 / np.pi
    new_vars = 'lat_mag', 'lon_mag', 'elev'
    vars.update(dict(zip(new_vars, new_vars)))
    return vals, vars


def load_swarm(fname):
    cdf = pycdf.CDF(fname)
    vars = {'Latitude': 'lat_geo',   # Geographic latitude
            'Longitude': 'lon_geo',  # geographic longitude
            'Absolute_STEC': 'tec',  # TEC in TECU
            'GPS_Position' : 'gps_pos',  # XYZ position of GPS (m)
            'LEO_Position' : 'leo_pos',  # XYZ position of SWARM (m)
            'Timestamp': 'times',  # Datetime times
            'PRN': 'prn'}  # GPS pseudo-random ID
    vals = {}
    for key, val in vars.items():
       vals[val] = cdf[key][...]

    return vals, vars


def localtime(fname):
    """
    Calculate local time of the Swarm satellites
    """
    vals = load_swarm(fname)
    outvals = {}
    utsec = np.array([(t - dt.datetime(t.year, t.month, t.day)).total_seconds() for t in vals['times']])
    lt = utsec / 3600 + vals['lon_geo'] / 360 * 24 
    lt[lt > 24] -= 24
    lt[lt < 0] += 24
    outvals['lt'] = lt
    outvals['times'] = np.array(vals['times'])
    return outvals 


if __name__ == '__main__':
    main()






