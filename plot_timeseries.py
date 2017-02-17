"""
plot_timeseries.py
Script to plot timeseries for comparison against the MIDAS contours
"""

import pdb 
import matplotlib.pyplot as plt 
import pickle
import collections
import datetime as dt
import numpy as np
import matplotlib.dates as mdates
import matplotlib
import socket


def main():
    starttime = dt.datetime(2015, 12, 20)
    endtime = dt.datetime(2015, 12, 21)
    sat = 'B'
    instrument = 'Langmuir Probe'  # or 'GPS'
    cutoff_crd = 'mag'

    # Langmuir Probe
    if instrument is 'Langmuir Probe':
        import proc_swarm_lp
        patch_ct, vals = proc_swarm_lp.main(time=starttime, endtime=endtime, sats=sat, save=False)
        plot_oneday_timeseries(patch_ct, vals) 

    elif instrument is 'TEC':        
        import proc_swarm_lp
        patch_ct, vals = proc_swarm_lp.main(time=starttime, endtime=endtime, sats=sat, save=False)
        plot_oneday_timeseries(patch_ct, vals) 

def plot_oneday_timeseries(patch_ct, vals, sat='B', \
                                           start=dt.datetime(2015, 12, 20, 16, 37, 30), \
                                           stop=dt.datetime(2015, 12, 20, 16, 55)):
    ut = np.array([t for t in vals[sat]['times']])
    timeind = np.logical_and(ut > start, ut < stop)
    ne = vals[sat]['ne'][timeind]
    ut = ut[timeind]
    mlat = vals[sat]['lat_mag'][timeind]

    count_ut = np.array([t[0] for t in patch_ct[sat]['times']])
    timeind = np.logical_and(count_ut > start, count_ut < stop)
    count = {}
    patch_ct[sat].pop('params')
    for key, val in patch_ct[sat].items():
        count[key] = np.array(val)[timeind]

    fig, ax1 = plt.subplots()
    utd = mdates.date2num(ut)
    # Ne timeseries
    plt.plot_date(utd, ne, 'b--')
    plt.plot_date(utd, ne, 'b.', markersize=3)

    # plot peak Ne
    plt.plot(mdates.date2num(count['times'][0]), count['ne'], 'kx', markersize=10, mew=4)

    # plot vertical lines at start and end of window
    maxval = ne.max() * 1.1
    plt.plot([mdates.date2num(count['t_start'][0]), mdates.date2num(count['t_start'][0])], [0, maxval], 'k--', mew=2)
    plt.plot([mdates.date2num(count['t_end'][0]), mdates.date2num(count['t_end'][0])], [0, maxval], 'k--', mew=2)

    # plot b1, b2 and bg levels of Ne
    plt.plot(mdates.date2num(count['t1'][0]), count['ne_b1'], 'gx', markersize=10, mew=4)
    plt.plot(mdates.date2num(count['t2'][0]), count['ne_b2'], 'gx', markersize=10, mew=4)
    plt.plot(mdates.date2num(count['times'][0]), count['ne_bg'], 'rx', markersize=10, mew=4)
    plt.plot([mdates.date2num(count['t1'][0]), mdates.date2num(count['t2'][0])], [count['ne_b1'], count['ne_b2']],'g--')

    fig.autofmt_xdate()
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    plt.xlabel(r'Time ($UT$)')
    plt.ylim(0, maxval) 
    plt.ylabel(r'Electron density ($cm^{-3}$)')
    plt.grid(which='both')
    t = ut.min()
    major_x = []
    while t <= ut.max():
        major_x.append(t)
        t += dt.timedelta(minutes=5)
    major_x = np.array(major_x)
    major_y = np.arange(0, maxval, 50000)                                              
    
    ax1.set_xticks(major_x)                                                       
    ax1.set_yticks(major_y)                                                       

    # plot magnetic latitude ticks
    ax2 = ax1.twiny()
    new_tick_locations = utd[0:-1:int(len(utd) / 5)]
    new_ticks = mlat[0:-1:int(len(utd) / 5)]
    new_ticklabels = ['%2.1f' % t for t in new_ticks]
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(new_ticklabels)
    ax2.set_xlabel(r"Mag. lat. ($deg$)")
    matplotlib.rcParams.update({'font.size': 24})
    plt.show()
 

if __name__ == '__main__':
    main()
