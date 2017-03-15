#!/Library/Frameworks/Python.framework/Versions/3.6/bin/python3
"""
plot_patch_ct.py
Script to plot the output of the SWARM patch counter (either TEC or LP)
"""
import pdb
import matplotlib.pyplot as plt
import pickle
import collections
import datetime as dt
import numpy as np
import socket

starttime = dt.datetime(2016, 1, 1)
endtime = dt.datetime(2016, 12, 31)
timestep = dt.timedelta(days=5)
satellites = 'A', 'B', 'C'
cutoff_crd = 'mag'
instrument = 'GPS'  # or 'GPS'

# Langmuir Probe
if instrument is 'Langmuir Probe':
    fin = './data/proc_lp/%s/' % cutoff_crd + 'lp_%Y%m%d.pkl'
    colour = 'b'

elif instrument is 'GPS':
    colour = 'r'
    if socket.gethostname() == 'chartat1-ml2':
        # Work GPS
        fin = '/Volumes/Seagate/data/swarm/proc/patch_ct_%Y%m%d.pkl'

    elif socket.gethostname() == 'chartat1-ml2':
        # Home GPS
        fin = './gps_proc/patch_ct_%Y%m%d.pkl'


def main():
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

    # plot_oneday_timeseries(patch_ct)
    # plot_utmlt(patch_ct, plot_type='MLT')
    # plot_utmlt(patch_ct, plot_type='UT')
    #plot_polar(patch_ct)
    plot_hist(patch_ct)


def plot_oneday_timeseries(patch_ct, sat='B', start=dt.datetime(2015, 12, 20, 16, 37, 30), \
                                               stop=dt.datetime(2015, 12, 20, 16, 55)):
    ut = np.array([t[0] for t in patch_ct[sat]['times']])
    pdb.set_trace()
    timeind = (ut > start) and (ut < stop)
    ne = patch_ct[sat]['n'][timeind]
    ne_err = patch_ct[sat]['n_error'][timeind]
    ut = ut[timeind]
    pdb.set_trace() 

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
        mlt = ut + mlon * 24 / 360
        mlt[mlt > 24] -= 24
        mlt[mlt < 0] += 24

        sat_lats = np.array([x[0] for x in patch_ct[sat]['lat_geo']])
        nh_ind = sat_lats > 0
        sh_ind = sat_lats < 0
        hems = 'north', 'south'
        for hem in hems:
            ct += 1
            plt.subplot(len(satellites), 2, ct)
            plt.xlim(0, 24)
            if instrument is 'GPS':
                plt.ylim(0, 800)
            else:
                plt.ylim(0, 600)
            if plot_type is 'UT':
                ut_h = ut[nh_ind] if hem is 'north' else ut[sh_ind]
                plt.hist(ut_h, bins=nbins)
            else:
                mlt_h = mlt[nh_ind] if hem is 'north' else mlt[sh_ind]
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
            doy_h = doy[nh_ind] if hem is 'north' else doy[sh_ind]
            plt.hist(doy_h, color=colour, bins=nbins)
            plt.title('Satellite %s, %s hemisphere' % (sat, hem))
            if np.mod(ct, 2) != 0:
                plt.ylabel('Patch count / 5 days')
            frame = plt.gca()
            plt.ylim(0, 300)
            plt.xlim(min(doy), max(doy))
            if ct >= 5:
                plt.xlabel('Day of year')
            else:
                frame.axes.xaxis.set_ticklabels([])
    plt.show()


def plot_polar(patch_ct):
    sats = [s for s in patch_ct.keys()]
    sats.sort()
    ct = 0  
    for sat in sats:
        lats = np.squeeze(np.array(patch_ct[sat]['lat_mag'])) * np.pi / 180
        lons = np.squeeze(np.array(patch_ct[sat]['lon_mag'])) * np.pi / 180
        lons[lons < 0] += 2 * np.pi

        latbins = np.arange(-91, 93, 2) * np.pi / 180
        lonbins = np.arange(-5, 375, 10) * np.pi / 180
        counts = np.histogram2d(lats, lons, np.array((latbins, lonbins)))
        vals = counts[0]
        latvec, lonvec = np.meshgrid((latbins[:-1] + latbins[1:]) / 2, (lonbins[:-1] + lonbins[1:]) / 2, indexing='ij')
        hems = collections.OrderedDict({'north': latvec > np.deg2rad(0), 'south': latvec < np.deg2rad(0)})
        
        for hem, ind in hems.items(): 
            ct += 1
            ax = plt.subplot(len(sats), len(hems), ct, polar=True)
            # labels = ['%2.0f' % (90 - np.rad2deg(val)) for val in np.arange(0.1, 0.6, 0.1)]
            # ax.set_yticklabels(labels)
            plt.ylim(0, 35 * np.pi / 180)
            hemlat = latvec
            if hem is 'south':
                hemlat = -hemlat
            sc = plt.pcolor(lonvec, np.pi / 2 - hemlat, vals)
            # sc.cmap.set_under('white')
            
            plt.clim(0, 150)
            plt.colorbar(sc)
            if ct < 3:
                plt.title('Sat %s: %s hemisphere' % (sat, hem))
            else:
                plt.title('Sat %s                                 ' % sat)

    plt.show()


if __name__ == '__main__':
    main()
