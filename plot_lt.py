#!/Library/Frameworks/Python.framework/Versions/3.4/bin/python3
"""
plot_lt.py
Script to plot the output of the Swarm langmuir probe
"""
import pdb
import matplotlib.pyplot as plt
import pickle
import datetime as dt
import numpy as np

starttime = dt.datetime(2016, 1, 1)
endtime = dt.datetime(2016, 2, 1)
timestep = dt.timedelta(days=1)
satellites = 'A', 'B', 'C'
fin = '/Volumes/Seagate/data/swarm/lt/lt_%Y%m%d.pkl'

vals = {}
time = starttime
while time < endtime:
    with open(time.strftime(fin), 'rb') as f:
        count = pickle.load(f)
    if vals == {}:
        vals = count
    else:
        for sat in satellites:
            try:
                vals[sat] = {key: np.append(vals[sat][key], val) for key, val in count[sat].items()}
            except:
                print('%s Missing file on satellite %s' % (time.strftime('%Y-%m-%d'), sat))
    time += dt.timedelta(days=1)

def main():
    # plot_polar()
    plot_hist()

def plot_hist():
    ct = 0
    for sat in satellites:
        doy = np.array([time.timetuple().tm_yday for time in vals[sat]['times']])
        nbins = round((endtime - starttime + dt.timedelta(days=1)) / timestep)
        ind = np.abs(vals[sat]['lt'] - 20) < 15/60
        ct += 1
        plt.subplot(len(satellites), 1, ct)
        plt.hist(doy[ind], bins=nbins)
        plt.title('Satellite %s' % sat)
        plt.ylabel('# Obs. < 15 mins from 20 LT')
        frame = plt.gca()
        plt.xlim(min(doy), max(doy))
        if ct > 2:
            plt.xlabel('Day of year')
        else:
            frame.axes.xaxis.set_ticklabels([])

    fig = plt.gcf()
    fig.suptitle("# of 2Hz Ne obs. within 15 mins of 20LT in Jan. 2016", fontsize=14)
    plt.show()

def plot_polar():
    sats = [s for s in vals.keys()]
    sats.sort()
    ct = 0  
    for sat in sats:
        lats = np.squeeze(np.array(vals[sat]['lat_mag']))
        lons = np.squeeze(np.array(vals[sat]['lon_mag'])) * np.pi / 180
        lons[lons < 0] += 2 * np.pi

        latbins = np.arange(-90, 92, 2)
        lonbins = np.arange(0, 2 * np.pi, 10)
        counts = np.histogram2d(lats, lons, np.array((latbins, lonbins)))
        xlat = (latbins[:-1] + latbins[1:]) / 2 
        ylon = (lonbins[:-1] + lonbins[1:]) / 2 
        theta, r = np.mgrid[ylon, xlat]
        vals = counts[0]
        hems = {'north': xlat > 60, 'south': xlat < -60}
        
        for hem, ind in hems.items(): 
            ct += 1
            ax = plt.subplot(len(sats), len(hems), ct, polar=True)
            plt.ylim(0, 30 * np.pi / 180)
            ax.set_yticklabels(['', '80', '', '70', '', '60'])
            colats = 90 - xlat
            if hem is 'south':
                colats = 180 - colats
            sc = plt.pcolor(ylon * np.pi / 180, colats[ind] * np.pi / 180, vals[ind, :], vmin=0.000001)
            sc.cmap.set_under('white')
            
            pdb.set_trace()
            plt.clim(0, 5)
            plt.colorbar(sc)
            if ct < 3:
                plt.title('Sat %s: %s hemisphere' % (sat, hem))
            else:
                plt.title('Sat %s                                 ' % sat)
    plt.show()

# def plot_dens_time():
    

if __name__ == '__main__':
    main()
