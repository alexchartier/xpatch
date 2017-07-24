"""
plot_timeseries.py
Script to plot timeseries for comparison against the MIDAS contours
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


def main(
        starttime = dt.datetime(2014, 12, 21),
          endtime = dt.datetime(2014, 12, 21, 23, 59),
        # starttime = dt.datetime(2016, 5, 8, 16, 37, 10),
        #   endtime = dt.datetime(2016, 5, 8, 17, 0, 0),
              sats=['A'],
        instrument='Langmuir Probe',
          approach='coley',
           procgps=False,
         ):

    # Langmuir Probe
    if instrument == 'Langmuir Probe':
        import proc_swarm_lp
        patch_ct, vals = proc_swarm_lp.main(time=starttime, endtime=endtime, approach=approach, sats=sats, save=False)
        plot_ne_timeseries(patch_ct, vals, sat=sats[0], start=starttime, stop=endtime) 

    elif instrument == 'GPS':        
        import proc_swarm_tec
        if procgps:
            patch_ct, vals = proc_swarm_tec.main(time=starttime, endtime=endtime, sats=sats, save=True)
        else:
            vals = {}
            for sat in sats:
                fname_format = '/Volumes/Seagate/data/swarm/gps_tec/' + 'SW_OPER_TEC%s' % sat + '*%Y%m%d*.cdf'
                fname = glob.glob(starttime.strftime(fname_format))[0]
                vals[sat], vars = proc_swarm_tec.get_swarm_vals(fname)
            with open(starttime.strftime('/Volumes/Seagate/data/swarm/proc_gps/patch_ct_%Y%m%d.pkl'), 'rb') as f:
                patch_ct = pickle.load(f)
        plot_tec_timeseries(patch_ct, vals, sat=sats[0], start=starttime, stop=endtime) 
        if socket.gethostname() == 'chartat1-ml2':
            # Work GPS
            fin = '/Volumes/Seagate/data/swarm/proc_gps/patch_ct_%Y%m%d.pkl'

        elif socket.gethostname() == 'chartat1-ml1':
            # Home GPS
            fin = './proc/gps/patch_ct_%Y%m%d.pkl'


def plot_tec_timeseries(patch_ct, vals, sat='B', \
                                           start=dt.datetime(2015, 12, 20, 16, 35), \
                                           stop=dt.datetime(2015, 12, 20, 16, 59, 32)):
    ut = np.array([t for t in vals[sat]['times']])
    mlat = vals[sat]['lat_mag']
    timeind = np.logical_and(ut > start, ut < stop)
    latind = np.abs(mlat) > 45
    ind = np.logical_and(latind, timeind)
    tec = vals[sat]['tec'][ind]
    ut = ut[ind]
    prn = vals[sat]['prn'][ind]
    mlat = mlat[ind]

    count_ut = np.array([t[0] for t in patch_ct[sat]['times']])
    timeind = np.logical_and(count_ut > start, count_ut < stop)
    count = {}
    patch_ct[sat].pop('params')
    for key, val in patch_ct[sat].items():
        count[key] = np.array(val)[timeind]

    fig, ax1 = plt.subplots()
    # loop over PRNs
    unique_prns = np.unique(prn)
    vals_p = {}
    utd = mdates.date2num(ut)

    # TEC timeseries
    for p in unique_prns:
        vals_prnind = prn == p
        count_prnind = count['prn'].flatten() == p
        plt.plot_date(utd[vals_prnind], tec[vals_prnind], '--')
        # plot peak TEC
        plt.plot(mdates.date2num(count['times'].flatten()[count_prnind]), count['tec'].flatten()[count_prnind], 'kx', markersize=10, mew=4)

        """     
        # plot b1, b2 and bg levels of Ne
        plt.plot(mdates.date2num(count['t1'].flatten()[count_prnind]), count['tec_b1'][count_prnind], 'gx', markersize=10, mew=4)
        plt.plot(mdates.date2num(count['t2'].flatten()[count_prnind]), count['tec_b2'][count_prnind], 'gx', markersize=10, mew=4)
        plt.plot([mdates.date2num(count['t1'].flatten()[count_prnind]), \
            mdates.date2num(count['t2'].flatten()[count_prnind])], [count['tec_b1'][count_prnind], count['tec_b2'][count_prnind]],'g--')
        """

    maxval = tec.max() * 1.1
    fig.autofmt_xdate()
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    plt.xlabel(r'Time ($UT$)')
    plt.ylim(0, maxval) 
    plt.ylabel(r'Total Electron Content ($10^{11} m^{-2}$)')
    plt.grid(which='both')
    """
    t = ut.min()
    major_x = []
    while t <= ut.max():
        major_x.append(t)
        t += dt.timedelta(minutes=5)
    major_x = np.array(major_x)
    major_y = np.arange(0, maxval, 2)                                              
    
    ax1.set_xticks(major_x)                                                       
    ax1.set_yticks(major_y)                                                       
    """

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

 
def plot_ne_timeseries(patch_ct, vals, sat='B', \
                                       lat_cutoff=-90, \
                                       start=dt.datetime(2015, 12, 20, 16, 35), \
                                       stop=dt.datetime(2015, 12, 20, 17, 5)):
    ut = np.array([t for t in vals[sat]['times']])
    mlat = vals[sat]['lat_mag']
    timeind = np.logical_and(ut > start, ut < stop)
    latind = mlat > lat_cutoff
    ind = np.logical_and(latind, timeind)
    ne = vals[sat]['ne'][ind]
    ne_rm = vals[sat]['ne_rm'][ind]
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

    # Ne timeseries
    obs = plt.plot_date(utd, ne, 'b--', label='Observed Ne')
    mf = plt.plot_date(utd, ne_rm, 'k', label='Median-smoothed Ne', linewidth=2)

    # plot peak and background Ne
    pk = plt.plot(mdates.date2num(np.squeeze(count['times'])), np.squeeze(count['ne_rm']), 'rx', markersize=10, mew=4, label='Peak')
    bg = plt.plot(mdates.date2num(np.squeeze(count['times'])), np.squeeze(count['ne_bg']), 'gx', label='background',  markersize=10, mew=4)

    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels)

    maxval = ne.max() * 1.1
    """
    # plot vertical lines at start and end of window
    plt.plot([mdates.date2num(count['t_start'][0]), mdates.date2num(count['t_start'][0])], [0, maxval], 'k--', mew=2)
    plt.plot([mdates.date2num(count['t_end'][0]), mdates.date2num(count['t_end'][0])], [0, maxval], 'k--', mew=2)

    # horizontal line between the two points
    plt.plot([mdates.date2num(count['t_start'][0]), mdates.date2num(count['t_end'][0])], [count['ne_bg'], count['ne_bg']], 'g--')
    """

    fig.autofmt_xdate()
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    plt.xlabel(r'Time ($UT$)')
    plt.ylim(0, maxval) 
    plt.ylabel(r'Electron density ($cm^{-3}$)')
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
    ax2 = ax1.twiny()
    new_tick_locations = utd[0:-1:int(len(utd) / 4)]
    new_ticks = mlat[0:-1:int(len(utd) / 4)]
    new_ticklabels = ['%2.1f' % t for t in new_ticks]
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(new_ticklabels)
    ax2.set_xlabel(r"Mag. lat. ($deg$)")
    matplotlib.rcParams.update({'font.size': 24})
    plt.show()
 

if __name__ == '__main__':
    main()
