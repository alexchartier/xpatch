#!/Library/Frameworks/Python.framework/Versions/3.6/bin/python3
"""
plot_patch_ct.py
Script to plot the output of the SWARM patch counter (either TEC or LP)
"""
import pdb
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pickle
import collections
import datetime as dt
import numpy as np
import socket
import sys
sys.path.insert(0, '/users/chartat1/fusionpp/fusion')
import physics
import count_passes

starttime = dt.datetime(2016, 1, 1)
endtime = dt.datetime(2016, 12, 31)
timestep = dt.timedelta(days=5)
satellites = 'A', 'B', 'C'
cutoff_crd = 'mag'
instrument = 'Langmuir Probe'  # or 'GPS'

# Langmuir Probe
if instrument == 'Langmuir Probe':
    fin = './data/proc_lp/%s/' % cutoff_crd + 'lp_%Y%m%d_70deg.pkl'

elif instrument == 'GPS':
    if socket.gethostname() == 'chartat1-ml2':
        # Work GPS
        fin = '/Volumes/Seagate/data/swarm/proc/patch_ct_%Y%m%d.pkl'

    elif socket.gethostname() == 'chartat1-ml2':
        # Home GPS
        fin = './gps_proc/patch_ct_%Y%m%d.pkl'

def main():
    patch_ct = get_patch_ct(starttime, endtime, satellites, fin)
    norm_patch_ct = normalize(patch_ct)
    # plot_oneday_timeseries(patch_ct)
    # plot_utmlt(patch_ct, plot_type='MLT')
    # plot_utmlt(patch_ct, plot_type='UT')
    # plot_magnitudes(patch_ct)
    plot_polar(patch_ct, crd='geo')
    # plot_hist(patch_ct)

def normalize(patch_ct, ipath='./data/pass_ct/pass_norm_%Y%m%d.pkl', 
                    starttime=dt.datetime(2016, 1, 1), 
                         step=dt.timedelta(days=1),
                      endtime=dt.datetime(2016, 12, 31),
                         sats=['A', 'B', 'C']):
    times = []
    time = starttime
    while time <= endtime:
        times.append(time)
        time += step
    
    day_count = {}
    hems = 'nh', 'sh'
    for sat in sats:
        day_count[sat] = {}
        for hem in hems:
            day_count[sat][hem + '_doy'] = np.zeros(len(times)) 
    for tind, t in enumerate(times):
        fin = t.strftime(ipath)
        with open(fin, 'rb') as f:
            ct = pickle.load(f)
        for sat in sats:
            for hem in hems:
                try:
                    day_count[sat][hem + '_doy'][tind] = ct[sat][hem + '_t'] 
                except:
                    print('No count for %s sat %s' % (t.strftime(ipath), sat))
    pdb.set_trace()
    return nh_ct, sh_ct 

def plot_magnitudes(patch_ct, magtype='absolute'):
    nbins = 50
    if magtype == 'relative':
        bins = np.linspace(2, 9, nbins)
        ylimit = [0, 3000]
        xlabel = 'Relative magnitude'
    else:
        bins = np.linspace(1E3, 2E5, nbins)
        ylimit = [0, 300]
    hems = 'north', 'south'
    ct = 0
    for sat in satellites:
        if magtype == 'relative':
            mag = np.array(patch_ct[sat]['ne']).flatten() / np.array(patch_ct[sat]['ne_bg']).flatten()
        elif magtype == 'absolute':
            mag = np.array(patch_ct[sat]['ne']).flatten() - np.array(patch_ct[sat]['ne_bg']).flatten()
        
        sat_lats = np.array([x[0] for x in patch_ct[sat]['lat_geo']])
        nh_ind = sat_lats > 0
        sh_ind = sat_lats < 0
        for hem in hems:
            ct += 1
            plt.subplot(len(satellites), 2, ct)
            mag_h = mag[nh_ind] if hem == 'north' else mag[sh_ind]
            print('Sat %s, %s hemisphere, median: %2.2g, mean: %2.2g, max: %2.2g' % (sat, hem, np.median(mag_h), np.mean(mag_h), np.max(mag_h)))
            pdb.set_trace()
            n, bins, patches = plt.hist(mag_h, bins=bins)
            plt.ylim(ylimit)
            if np.ceil(ct / 2) == len(satellites):
                plt.xlabel(r'Absolute magnitude (electrons / $cm^{-3}$)')
            else:
                plt.tick_params(
                                axis='x',          # changes apply to the x-axis
                                which='both',      # both major and minor ticks are affected
                                bottom='off',      # ticks along the bottom edge are off
                                top='off',         # ticks along the top edge are off
                                labelbottom='off') # labels along the bottom edge are off
            if np.mod(ct, 2) != 0:
                plt.ylabel('Patch count')
            
            plt.grid()
            plt.title('Satellite %s %s hemisphere' % (sat, hem))

    plt.suptitle(instrument, fontweight='bold')
    plt.show() 

def plot_oneday_timeseries(patch_ct, sat='B', start=dt.datetime(2015, 12, 20, 16, 37, 30), \
                                               stop=dt.datetime(2015, 12, 20, 16, 55)):
    ut = np.array([t[0] for t in patch_ct[sat]['times']])
    timeind = (ut > start) and (ut < stop)
    ne = patch_ct[sat]['n'][timeind]
    ne_err = patch_ct[sat]['n_error'][timeind]
    ut = ut[timeind]

def plot_utmlt(patch_ct, plot_type='UT'):
    sats = [s for s in patch_ct.keys()]
    sats.sort()
    nbins = 24
    ct = 0  
    for sat in sats:
        mlon = np.squeeze(np.array(patch_ct[sat]['lon_mag']))
        mlon[mlon < 0] += 2 * np.pi
        time = {}
        ut = np.array([t[0].hour + t[0].minute / 60 for t in patch_ct[sat]['times']])
        mlt = calc_mlt(ut, mlon)

        sat_lats = np.array([x[0] for x in patch_ct[sat]['lat_geo']])
        nh_ind = sat_lats > 0
        sh_ind = sat_lats < 0
        hems = 'north', 'south'
        for hem in hems:
            ct += 1
            plt.subplot(len(satellites), 2, ct)
            plt.xlim(0, 24)
            if instrument == 'GPS':
                plt.ylim(0, 800)
            else:
                plt.ylim(0, 400)
            if plot_type is 'UT':
                ut_h = ut[nh_ind] if hem == 'north' else ut[sh_ind]
                plt.hist(ut_h, bins=nbins)
            else:
                mlt_h = mlt[nh_ind] if hem == 'north' else mlt[sh_ind]
                n, bins, patches = plt.hist(mlt_h, bins=nbins)
            if np.ceil(ct / 2) == len(sats):
                plt.xlabel('%s Hour' % plot_type)
            else:
                plt.tick_params(
                                axis='x',          # changes apply to the x-axis
                                which='both',      # both major and minor ticks are affected
                                bottom='off',      # ticks along the bottom edge are off
                                top='off',         # ticks along the top edge are off
                                labelbottom='off') # labels along the bottom edge are off
            if np.mod(ct, 2) != 0:
                plt.ylabel('Patch count')
            plt.grid()
            plt.title('Satellite %s %s hemisphere' % (sat, hem))

    plt.suptitle(instrument, fontweight='bold')
    plt.show()
                


def plot_hist(patch_ct):
    ct = 0
    for sat in satellites:
        doy = np.array([time[0].timetuple().tm_yday for time in patch_ct[sat]['times']])
        nbins = round((endtime - starttime + dt.timedelta(days=1)) / timestep)
        sat_lats = np.array([x[0] for x in patch_ct[sat]['lat_geo']])
        nh_ind = sat_lats > 0
        sh_ind = sat_lats < 0
        hems = 'north', 'south'
        for hem in hems:
            ct += 1
            plt.subplot(len(satellites), 2, ct)
            doy_h = doy[nh_ind] if hem == 'north' else doy[sh_ind]
            plt.hist(doy_h, bins=nbins)
            plt.title('Satellite %s, %s hemisphere' % (sat, hem))
            plt.ylabel('Patch count')
            frame = plt.gca()
            plt.ylim(0, 150)
            plt.xlim(min(doy), max(doy))
            if ct >= 5:
                plt.xlabel('Day of year')
            else:
                frame.axes.xaxis.set_ticklabels([])
    plt.show()


def plot_polar(patch_ct, crd='mag'):
    passes = sum_passes('./data/pass_ct/pass_%Y%m%d.pkl', crd=crd)
    sats = [s for s in patch_ct.keys()]
    sats.sort()

    hems = 'north', 'south'
    latbins = np.deg2rad(np.arange(-91, 91.1, 2)) 
    lonbins = np.deg2rad(np.arange(-5, 365.1, 10))
    if crd == 'mag':
        latlim = 70
    else:
        latlim = 60

    ct = 0  
    for sat in sats:
        if crd == 'mag':
            lats = np.deg2rad(np.squeeze(np.array(patch_ct[sat]['lat_mag'])))
            lons = np.deg2rad(np.squeeze(np.array(patch_ct[sat]['lon_mag'])))
            unused_alts, pole_lat, pole_lon = physics.transform([1, 1], np.deg2rad([90, -90]), [0, 0], from_=['GEO', 'sph'], to=['MAG', 'sph'])
        elif crd == 'geo':
            lats = np.deg2rad(np.squeeze(np.array(patch_ct[sat]['lat_geo'])))
            lons = np.deg2rad(np.squeeze(np.array(patch_ct[sat]['lon_geo'])))
            unused_alts, pole_lat, pole_lon = physics.transform([1, 1], np.deg2rad([90, -90]), [0, 0], from_=['MAG', 'sph'], to=['GEO', 'sph'])
        pole_lat[1] = - pole_lat[1]
        lons[lons < 0] += 2 * np.pi

        counts = np.histogram2d(lats, lons, np.array((latbins, lonbins)))
        vals = counts[0] / passes[sat][crd]
        vals[vals == 0] = np.nan
        latvec, lonvec = np.meshgrid((latbins[:-1] + latbins[1:]) / 2, (lonbins[:-1] + lonbins[1:]) / 2, indexing='ij')
        
        for hem in hems: 
            ct += 1
            ax = plt.subplot(len(sats), len(hems), ct, polar=True)
            hemlat = latvec
            if hem == 'south':
                hemlat = -hemlat
            plt.ylim(0, np.deg2rad(90 - latlim))
            pdb.set_trace()
            sc = plt.pcolor(lonvec, np.pi / 2 - hemlat, vals)
            hemind = 1 if hem == 'south' else 0
            plt.plot(pole_lon[hemind], np.pi / 2 - pole_lat[hemind], '.m', markersize=15)
            labels = ['%2.0f' % (90 - val) for val in np.linspace(0, 90 - latlim, 7)]
            labels = labels[1:]
            ax.set_yticklabels(labels)
            sc.cmap.set_under('white')
            """
            m = Basemap(projection='npstere',boundinglat=70,lon_0=270,resolution='l')
            # draw parallels and meridians.
            #m.drawparallels(np.arange(-70.,81.,20.))
            #m.drawmeridians(np.arange(-180.,181.,20.))
            sc = m.pcolor(lonvec, hemlat, vals)
            """
            plt.clim(0, 0.1)
            plt.colorbar(sc)
            if ct < 3:
                plt.title('Sat %s: %s hemisphere' % (sat, hem))
            else:
                plt.title('Sat %s                                 ' % sat)

    plt.show()


def sum_passes(fname_format, crd='mag'):
    time = starttime

    while time <= endtime:
        with open(time.strftime(fname_format), 'rb') as f:
            pass_ct = pickle.load(f)
        if time == starttime:
            pass_ct_full = pass_ct
        else:
            for s in satellites:
                try:
                    pass_ct_full[s][crd] += pass_ct[s][crd]
                except:
                    print('Missing file for satellite %s on %s' % (s, time.strftime('%Y %m %d')))
        time += dt.timedelta(days=1)

    for sat in satellites:
        pass_ct_full[sat][crd][pass_ct_full[sat][crd] == 0] = np.nan
    return pass_ct_full


def get_patch_ct(starttime, endtime, satellites, fin):
    patch_ct = {}
    time = starttime
    while time < endtime:
        with open(time.strftime(fin), 'rb') as f:
            count = pickle.load(f)
        if patch_ct == {}:
            patch_ct = count
        else:
            for sat in satellites:
                try:
                    for key, val in count[sat].items():
                        if key != 'params':
                            patch_ct[sat][key] = patch_ct[sat][key] + val
                except:
                    print('%s Missing file on satellite %s' % (time.strftime('%Y-%m-%d'), sat))
        time += dt.timedelta(days=1)
    return patch_ct

def calc_mlt(ut, mlon):
    mlt = ut + mlon * 24 / 360
    mlt[mlt > 24] -= 24
    mlt[mlt < 0] += 24
    return mlt

if __name__ == '__main__':
    main()
