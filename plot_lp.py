#!/Library/Frameworks/Python.framework/Versions/3.4/bin/python3
"""
plot_lp.py
Script to plot the output of the Swarm langmuir probe
"""
import pdb
import matplotlib.pyplot as plt
import pickle
import datetime as dt
import numpy as np
import glob
import proc_swarm_lp
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.dates import DateFormatter, MinuteLocator
import numpy as np


def main():
    ipath = './data/swarm_lp/'
    opath = './data/proc_lp/'
    time = dt.datetime(2015, 12, 20, 18, 45)
    window = dt.timedelta(minutes=30)
    cutoff_crd = 'geo'
    sats = 'A', 'B', 'C' 
    vals = {}
    patch_ct = {}    

    print('Loading data')
    for sat in sats:
        print('\nSatellite %s' % sat)
        fname_format = ipath + 'SW_OPER_EFI%s' % sat + '_PL*%Y%m%d*.cdf' 
        try:
            fname = glob.glob(time.strftime(fname_format))[0]
            vals[sat] = proc_swarm_lp.load_lp(fname)
            vals[sat]['lt'] = proc_swarm_lp.localtime(vals[sat])
        except:
            print('No file for satellite %s on %s' % (sat, timestr))
        
        # Time filter
        timeind = np.logical_and((vals[sat]['times'] >= time), (vals[sat]['times'] <= time + window))
        for key, val in vals[sat].items(): 
            vals[sat][key] = val[timeind]
    
    #plt.subplot(2, 1, 1) 
    #plot_time_series(vals, key='A')
    vals.pop('B')
    vals.pop('C')
    # plt.subplot(2, 1, 2) 
    plot_sat_pos(vals, plot_type='polar')
    plt.show()

        
def plot_time_series(vals, key='A', xaxis='time'):
    """
    Plot the satellite output vs. time. 
    Give latitude and longitude as additional x axes
    """
    val = vals[key]   
 
    # Plot the data
    dates = matplotlib.dates.date2num(val['times'])
    ne = val['ne'] / 1E5
    if xaxis is 'lat':
        plt.plot(val['lat_geo'], ne, '.')
        plt.xlabel(r'Latitude $(degrees)$')
    else:  # time
        plt.plot_date(dates, ne, '.')
        ax = plt.gca()
        ax.xaxis.set_major_locator(MinuteLocator(byminute=range(0, 60, 10)))
        ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
        plt.xlabel(r'Time $(UT)$')

    centre_time = val['times'][np.abs(dates - np.median(dates)).argmin()]
    # plt.title('Satellite %s centred on %s' % (key, centre_time.strftime('%H:%M:%S %d %B %Y')))
    plt.ylabel(r'Electron density $(10^5 m^-3)$')
    plt.grid()


def plot_sat_pos(vals, key='A', plot_type='3dglobe'):
    """
    Plot the satellite position on a map
    Set up orthographic map projection with perspective of satellite looking down at 50N, 100W.
    Use low resolution coastlines.
    Don't plot features that are smaller than 1000 square km.
    """
    # For 3D globe:
    if plot_type is '3dglobe':
        map = Basemap(projection='ortho', lat_0=70, lon_0=0,
                      resolution='l', area_thresh=1000.)

    elif plot_type is 'miller':
        # Miller cylindrical example (suitable for global proj.)
        lat_bound=[-85, 85]
        lon_bound=[-180, 180]
        lon_0 = 90
        map = Basemap(projection='mill', llcrnrlat=lat_bound[0], urcrnrlat=lat_bound[1],
               llcrnrlon=lon_bound[0], urcrnrlon=lon_bound[1], resolution='l', lon_0=lon_0)
    elif plot_type is 'polar':
        map = Basemap(projection='npstere', boundinglat=28, lon_0=270, resolution='l')
        # draw parallels and meridians.

    # draw coastlines, country boundaries, fill continents.
    map.drawcoastlines()
    map.drawcountries()
    map.fillcontinents(color='coral',lake_color='aqua')
    # draw the edge of the map projection region (the projection limb)
    map.drawmapboundary(fill_color='aqua')
    # draw lat/lon grid lines every 30 degrees.
    map.drawparallels(np.arange(-90, 90, 30), labels=[False, True, True, False])
    map.drawmeridians(np.arange(-180.,181.,20.), labels=[True, False, False, True])

    # Plot the data
    colours = 'r', 'g', 'b'
    count = 0
    for key, val in vals.items():
        x, y = map(val['lon_geo'], val['lat_geo'])
        map.plot(x, y, colours[count] + '.')
        count += 1







   

if __name__ == '__main__':
    main()
