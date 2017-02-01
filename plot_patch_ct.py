#!/Library/Frameworks/Python.framework/Versions/3.4/bin/python3
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

starttime = dt.datetime(2016, 1, 1)
endtime = dt.datetime(2016, 12, 31)
timestep = dt.timedelta(days=5)
satellites = 'A', 'B', 'C'
cutoff_crd = 'geo'
fin = './data/proc_lp/%s/' % cutoff_crd + 'lp_%Y%m%d.pkl'


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

    # plot_polar(patch_ct)
    plot_hist(patch_ct)


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
            plt.hist(doy_h, bins=nbins)
            plt.title('Satellite %s, %s hemisphere' % (sat, hem))
            plt.ylabel('Patch count')
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
