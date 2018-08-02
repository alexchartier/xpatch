# Find pass counts that are different across Langmuir Probe and GPS
# Could also find 'unique' pass counts in each GPS - removing patches within X seconds of other patches

import pdb 
import matplotlib.pyplot as plt 
import matplotlib
import pickle
import collections
import datetime as dt
import numpy as np
import socket
import sys 
sys.path.insert(0, '/users/chartat1/fusionpp/fusion')
import count_passes
import plot_patch_ct

def main():
    starttime = dt.datetime(2014, 8, 1)
    endtime = dt.datetime(2017, 7, 1)
    satellites = 'A', 'B' 

    approach = 'coley'
    lat_cutoff = 70
    lp_fin = '/Volumes/Seagate/data/swarm/proc_lp_comb/%s/' % approach + 'lp_%Y%m%d_' + '%ideg.pkl' % lat_cutoff
    gps_fin = '/Volumes/Seagate/data/swarm/proc_gps/patch_ct_%Y%m%d.pkl'

    lp_patch_ct = plot_patch_ct.get_patch_ct(starttime, endtime, satellites, lp_fin)
    gps_patch_ct = plot_patch_ct.get_patch_ct(starttime, endtime, satellites, gps_fin)

    gps_pop = 't1', 't2', 'tec_b1', 'tec_b2', 'params'
    lp_pop = ['params']

    for sat in satellites:
        # remove unnecessary misshapen variables
        for p in gps_pop:
            gps_patch_ct[sat].pop(p)
        for p in lp_pop:
            lp_patch_ct[sat].pop(p)

        # convert to arrays
        for key, val in lp_patch_ct[sat].items():
            lp_patch_ct[sat][key] = np.squeeze(np.array(val))
        for key, val in gps_patch_ct[sat].items():
            gps_patch_ct[sat][key] = np.squeeze(np.array(val))

        # filter out low-lat GPS
        lat_ind = np.abs(gps_patch_ct[sat]['lat_mag']) >= lat_cutoff
        for key, val in gps_patch_ct[sat].items():
            gps_patch_ct[sat][key] = val[lat_ind]

        # Find GPS patches within 10-s of each LP patch
            # 1. convert datetimes to numbers and round off to desired precision.
            # 2. Do a set comparison to get indices

        gps_ts = np.round(np.array([t.timestamp() for t in gps_patch_ct[sat]['times']]) / 10)
        lp_ts = np.round(np.array([t.timestamp() for t in lp_patch_ct[sat]['times']]) / 10)

        gps_ts_set = set(gps_ts)
        lp_ts_set = set(lp_ts)
        both_sets = np.array(list(gps_ts_set & lp_ts_set))
        just_gps = np.array(list(gps_ts_set - lp_ts_set))
        just_lp = np.array(list(lp_ts_set - gps_ts_set))

        both_gps_ind = np.in1d(gps_ts, both_sets)  # index of GPS times where both GPS and LP saw a patch
        just_gps_ind = np.in1d(gps_ts, just_gps)   # index of GPS times where both GPS and LP saw a patch
        nh_ind = gps_patch_ct[sat]['lat_geo'] > 0
        sh_ind = gps_patch_ct[sat]['lat_geo'] < 0
        
        both_gps_times_sh = [t for t in gps_patch_ct[sat]['times'][np.logical_and(both_gps_ind, sh_ind)]]
        just_gps_times_sh = [t for t in gps_patch_ct[sat]['times'][np.logical_and(just_gps_ind, sh_ind)]]
        start = dt.datetime(2014, 12, 15)
        stop = dt.datetime(2014, 12, 30)
        pdb.set_trace()
        print(np.array(just_gps_times_sh)[np.logical_and(np.array(just_gps_times_sh) > start, np.array(just_gps_times_sh) < stop)])

        # plot
        both_gps_timestamps = np.array([t.timestamp() for t in both_gps_times]) / 86400
        plt.hist(both_gps_timestamps, 50)
        plt.show()

if __name__ == '__main__':
    main()

