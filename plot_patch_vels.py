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


def main(
        starttime = dt.datetime(2016, 12, 21, 0,  0),
          endtime = dt.datetime(2016, 12, 21, 23, 0),
              sats=['A'],
          approach='alex',
         ):

    # pdb.set_trace()
    # '/Volumes/Seagate/data/swarm_efi/SW_EXPT_EFIA_TIICT_20161231T131014_20161231T223104_0101.cdf
    # Langmuir Probe
    lat_cutoff = 55
    import proc_swarm_lp
    patch_ct, vals = proc_swarm_lp.main(time=starttime, endtime=endtime, approach=approach, sats=sats, lat_cutoff=55, save=False)
    plot_timeseries(patch_ct, vals, sat=sats[0], start=starttime, stop=endtime)

    
def plot_timeseries(patch_ct, vals, sat='B', \
                             lat_cutoff=55, \
                                  start=dt.datetime(2015, 12, 20, 16, 35), \
                                   stop=dt.datetime(2015, 12, 20, 17, 5)):
    ut = np.array([t for t in vals[sat]['times']])
    mlat = vals[sat]['lat_mag']
    timeind = np.logical_and(ut > start, ut < stop)
    latind = np.abs(mlat) > lat_cutoff
    ind = np.logical_and(latind, timeind)
    Te = vals[sat]['T_elec'][ind]
    ut = ut[ind]
    mlat = mlat[ind]
    count_ut = np.array([t[0] for t in patch_ct[sat]['times']])
    timeind = np.logical_and(count_ut > start, count_ut < stop)
    count = {}
    patch_ct[sat].pop('params')
    for key, val in patch_ct[sat].items():
        count[key] = np.array(val)[timeind]

    fig, ax1 = plt.subplots()
    utd = mdates.date2num(ut)

    # Te timeseries
    obs = plt.plot_date(utd, Te, 'b.', label='Observed Te')

    # plot peak Te
    pk = plt.plot(mdates.date2num(count['times'].tolist()), np.squeeze(count['T_elec']), 'rx', markersize=20, mew=4, label='Patch')
    # bg = plt.plot(mdates.date2num(count['times'].tolist()), np.squeeze(count['ne_bg']), 'gx', label='background',  markersize=10, mew=4)

    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels)

    # maxval = ne.max() * 1.1
    maxval = 1E4
    """
    # plot vertical lines at start and end of window
    plt.plot([mdates.date2num(count['t_start'][0]), mdates.date2num(count['t_start'][0])], [0, maxval], 'k--', mew=2)
    plt.plot([mdates.date2num(count['t_end'][0]), mdates.date2num(count['t_end'][0])], [0, maxval], 'k--', mew=2)

    # horizontal line between the two points
    plt.plot([mdates.date2num(count['t_start'][0]), mdates.date2num(count['t_end'][0])], [count['ne_bg'], count['ne_bg']], 'g--')
    """

    fig.autofmt_xdate()
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.xlabel(r'Time ($UT$)')
    plt.ylim(0, maxval)
    plt.ylabel(r'Electron Temperature ($K$)')
    plt.grid(which='both')
    t = ut.min()
    """
    major_x = []
    while t <= ut.max():
        major_x.append(t)
        t += dt.timedelta(minutes=5)
    major_x = np.array(major_x)
    major_y = np.arange(0, maxval, 50000)                                              
    
    ax1.set_xticks(major_x)                                                       
    ax1.set_yticks(major_y)                                                       
    """

    # plot magnetic latitude ticks
    nticks = 15
    ax2 = ax1.twiny()
    new_tick_locations = utd[0:-1:int(len(utd) / nticks)]
    new_ticks = mlat[0:-1:int(len(utd) / nticks)]
    new_ticklabels = ['%2.1f' % t for t in new_ticks]
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(new_ticklabels)
    ax2.set_xlabel(r"Mag. lat. ($deg$)")
    matplotlib.rcParams.update({'font.size': 24})
    # plt.title(start.strftime('%Y/%b/%d'))
    plt.show()


if __name__ == '__main__':
    main()
