"""
plot_eclipse_tec.py
Script to plot TEC for the eclipse
"""
import glob
import pdb 
import matplotlib.pyplot as plt 
import pickle
import collections
import datetime as dt
import numpy as np
import matplotlib.dates as mdates
import matplotlib
import socket
import proc_swarm_tec


def main(
        starttime = dt.datetime(2017, 8, 16, 17, 0), 
          endtime = dt.datetime(2017, 8, 16, 20, 0),
             sats = ['A', 'B'],
         ):

    vals = {}
    for sat in sats:
        fname_format = 'data/eclipse/gps/' + 'SW_OPER_TEC%s' % sat + '*%Y%m%d*.cdf'
        fname = glob.glob(starttime.strftime(fname_format))[0]
        vals[sat], vars = proc_swarm_tec.get_swarm_vals(fname)
    plot_tec_timeseries(vals, sat=sats[0], start=starttime, stop=endtime) 


def plot_tec_timeseries(vals, sat='B', \
                                      start=dt.datetime(2015, 12, 20, 16, 35), \
                                       stop=dt.datetime(2015, 12, 20, 16, 59, 32), \
                                 lat_bounds=[15, 70], \
                                 lon_bounds=[-135, -65], \
                                    ):
    ut = np.array([t for t in vals[sat]['times']])
    lat = vals[sat]['lat_geo']
    lon = vals[sat]['lon_geo']
    timeind = np.logical_and(ut > start, ut < stop)
    latind = np.logical_and(lat > lat_bounds[0], lat < lat_bounds[1])
    lonind = np.logical_and(lon > lon_bounds[0], lon < lon_bounds[1])
    ind = np.logical_and(np.logical_and(latind, timeind), lonind)
    tec = vals[sat]['vtec'][ind]
    ut = ut[ind]
    prn = vals[sat]['prn'][ind]
    lat = lat[ind]
    lon = lon[ind]

    fig, ax1 = plt.subplots()
    # loop over PRNs
    unique_prns = np.unique(prn)
    vals_p = {}
    utd = mdates.date2num(ut)

    # TEC timeseries
    ctr = 0
    for p in unique_prns:
        vals_prnind = prn == p
        plt.plot_date(utd[vals_prnind][::5], tec[vals_prnind][::5], '.')

        if ctr == 0:
            handles, labels = ax1.get_legend_handles_labels()
            ax1.legend(handles, labels)
            ctr += 1

    maxval = 15
    # maxval = tec.max() * 1.1
    fig.autofmt_xdate()
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.xlabel(r'Time ($UT$)')
    plt.xlim(mdates.date2num(start), mdates.date2num(stop)) 
    plt.ylim(0, maxval) 
    plt.ylabel(r'Total Electron Content ($10^{16} m^{-2}$)')
    plt.grid(which='both')
    plt.title('Vertical TEC %s to %s UT' % (start.strftime('%Y-%m-%d %H:%M'), stop.strftime('%H:%M')), y=1.08)

    # plot magnetic latitude ticks
    nticks = 6
    ax2 = ax1.twiny()
    new_tick_locations = utd[0:-1:int(len(utd) / nticks)]
    new_ticks = lon[0:-1:int(len(utd) / nticks)]
    new_ticklabels = ['%2.0f' % t for t in new_ticks]
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(new_ticklabels)
    ax2.set_xlabel(r"Geo. lon. ($deg$)")
    matplotlib.rcParams.update({'font.size': 18})
    plt.show()

if __name__ == '__main__':
    main()
