#!/usr/local/bin/python3
"""
proc_swarm_tec.py
Script to process the SWARM upward-looking TEC data and analyse for patches. 
This is intended to follow the method of Noja et al. [2013]
No slant-to-vertical conversion (justified in the paper)
Elevation cutoff: 25 degrees
Latitude cutoff: 55 degrees geomagnetic
Sliding window: 200s
Check that there is a positive slope followed by a negative slope. No minimum slope magnitude is required
TECp: Patch peak magnitude is the largest value between the two slopes
TECbg: Background level is calculated by symmetric linear interpolation of the 'boundary' values to the location of the peak
TECstart: First boundary. Find smallest value to the left and right of peak. Take the greatest of those two. 
TECend: Second boundary. On the other side of the peak, find value closest in magnitude to the first boundary. 
Relative magnitude criterion: (TECp - TECbg) ** 2 / TECbg >= 1.2
Absolute magnitude criterion: TECp - TECstart >= 4 TECU and TECp - TECend >= 4
"""
from spacepy import pycdf
import pdb
import numpy as np
import scipy as sp
import datetime as dt
import matplotlib.pyplot as plt
import glob
import pickle
import sys 
import collections
sys.path.insert(0, '/users/chartat1/fusionpp/fusion/')
import physics
import socket

def main(
        time=dt.datetime(2014, 1, 1),
        step=dt.timedelta(days=1),
        endtime=dt.datetime(2014, 8, 1),
        sats=['A', 'B'],
        save=True
        ):

    if socket.gethostname() == 'chartat1-ml2':
        ipath = '/Volumes/Seagate/data/swarm/gps_tec/'
        opath = '/Volumes/Seagate/data/swarm/proc_gps/'
    elif socket.gethostname() == 'chartat1-ml1':
        ipath = 'data/swarm_tec/'
        opath = 'data/proc_gps/'
    
    while time <= endtime: 
        timestr = time.strftime('%Y-%m-%d')
        print(timestr)
        patch_ct = {}
        vals = {}

        for sat in sats:
            print('Satellite %s' % sat)
            fname_format = time.strftime(ipath) + 'SW_OPER_TEC%s' % sat + '*%Y%m%d*.cdf'
            try:
                fname = glob.glob(time.strftime(fname_format))[0]
            except:
                print('No file for satellite %s on %s' % (sat, timestr))
                continue
            patch_ct[sat], vals[sat] = count_patches(fname)

        if save:
            fout = opath + time.strftime('patch_ct_%Y%m%d.pkl')
            with open(fout, 'wb') as f:
                pickle.dump(patch_ct, f) 
            print('Saving %s' % fout)

        time += dt.timedelta(days=1)

    return patch_ct, vals


def count_patches(fname, lat_cutoff=55, elev_cutoff=25, TEC_abs_cutoff=4, TEC_rel_cutoff=1.2, window_sec=200, cadence_sec=1):
    """
    Counts the patches in each SWARM file

    Inputs: 
        fname = '~/Downloads/SW_OPER_TECATMS_2F_20131125T105953_20131125T235953_0201.DBL'
        lat_cutoff  # degrees magnetic
        elev_cutoff  # degrees
        TEC_abs_cutoff  # TECU
        TEC_rel_cutoff  #  Relative cutoff
        window  # seconds
        cadence  # seconds
    
    Returns: 
        patch_ct - A dictionary of the attributes of each patch
    """ 
    window = dt.timedelta(seconds=window_sec)  
    cadence = dt.timedelta(seconds=cadence_sec) 
    vals, vars = get_swarm_vals(fname)
    # Take out values below elevation and latitude cutoffs 
    index = np.logical_and(vals['elev'] >= elev_cutoff, np.abs(vals['lat_mag']) > lat_cutoff)
    for key, val in vars.items():  
        vals[val] = vals[val][index]

    # loop over PRNs
    unique_prns = np.unique(vals['prn'])
    vals_p = {}
    patch_ct = {}
    for key, var in vars.items():
        patch_ct[var] = []
    new_varnames = 'tec_bg', 'tec_b1', 'tec_b2', 't1', 't2'
    for v in new_varnames:
        patch_ct[v] = []

    for p in unique_prns:
        index = vals['prn'] == p
        p = str(p)
        vals_p[p] = {}
        for key, val in vars.items():
            vals_p[p][val] = vals[val][index]

        # Problem with duplicate values in SWARM data. This will fix the problems it causes.
        if len(set(vals_p[p]['tec'])) < len(vals_p[p]['tec']):
            print('Found duplicates in TEC data')
            vals_p[p]['tec'] += np.random.randn(len(vals_p[p]['tec'])) * 1E-9

        # Sliding window filtering
        tind = 0
        while tind < len(vals_p[p]['times']) - window / cadence:
            # NOTES FOR PUBLICATION: 
            #   Not clear if Noja et al. skipped forward if they found a patch within the window. Assume they did
            #   What happens if there are less than max. points available in a segment? We throw the segment out.
            #   There seem to be a lot of repeated TEC values
            t = vals_p[p]['times'][tind]
            tind += 1
            ind = np.logical_and(vals_p[p]['times'] >= t, vals_p[p]['times'] <= t + window)
            if sum(ind) < window / cadence:  # Require full window
                continue

            times = vals_p[p]['times'][ind]
            tec_vals = vals_p[p]['tec'][ind]
            grads = np.diff(tec_vals)
            mag_lats = vals_p[p]['lat_mag'][ind]

            # Check window does not cut across a hemispheric boundary
            lat_steps = np.diff(mag_lats)
            if lat_steps.max() > lat_cutoff:
                print('Found hemispheric jump at %s' % t)
                # vals['times'].tolist().index(times[lat_steps == lat_steps.max()][0])
                # tind += 
                continue

            # Algorithm requires a positive gradient ...
            if grads.max() <= 0:
                continue
            pos_ind = np.argmax(grads >= 0)
            # ... followed by a negative gradient
            if grads[pos_ind:].min() >= 0:
                continue

            TECp = tec_vals.max()  # Patch maximum is the highest value in the window
            ind_TECp = np.where(tec_vals == TECp)[0][0]
            
            # The next part won't work if the first/last value is the largest
            if (ind_TECp == 0) or (ind_TECp == len(tec_vals) - 1):
                continue
            
            # Define the two boundary values: b1 - greater of the two minimum values either side
            #                                 b2 - closest value on other side to b1     
            lhs = tec_vals[:ind_TECp]
            rhs = tec_vals[ind_TECp + 1:]
            TEC_b1 = max([np.min(lhs), np.min(rhs)])
            if TEC_b1 in lhs:
                TEC_b2 = rhs[(np.abs(rhs - TEC_b1)).argmin()]
            else:
                TEC_b2 = lhs[(np.abs(lhs - TEC_b1)).argmin()]
        
            # Determine background TEC: Symmetric linear interpolation of values to location of the peak
            ind_TEC_b1 = np.where(tec_vals == TEC_b1)[0][0]
            ind_TEC_b2 = np.where(tec_vals == TEC_b2)[0][0]  
            time_b1 = (times[ind_TEC_b1] - times.min()).seconds
            time_b2 = (times[ind_TEC_b2] - times.min()).seconds
            time_p = (times[ind_TECp] - times.min()).seconds
            try:
                TECbg = sp.interpolate.interp1d([time_b1, time_b2], [TEC_b1, TEC_b2])(time_p).tolist()
            except:
                pdb.set_trace()

            # Relative magnitude test
            if (TECp - TECbg) ** 2 / TECbg < TEC_rel_cutoff:
                continue
            
            # Absolute magnitude test
            if (TECp - TEC_b1 < TEC_abs_cutoff) or (TECp - TEC_b2 < TEC_abs_cutoff):
                continue

            # If we're still going at this point, we have found a patch. Store the details and skip forward to the next window
            patch_index = vals_p[p]['tec'] == TECp
            assert patch_index.sum() == 1, 'There should be exactly one patch index for each patch'
            for key, var in vars.items():
                patch_ct[var].append(vals_p[p][var][patch_index])
            patch_ct['tec_bg'].append(TECbg)
            patch_ct['t1'].append(times[ind_TEC_b1])
            patch_ct['t2'].append(times[ind_TEC_b2])
            patch_ct['tec_b1'].append(TEC_b1)
            patch_ct['tec_b2'].append(TEC_b2)
            tind += sum(ind)

    patch_ct['params'] = {'lat_cutoff': lat_cutoff,
                         'elev_cutoff': elev_cutoff,
                      'TEC_abs_cutoff': TEC_abs_cutoff,
                      'TEC_rel_cutoff': TEC_rel_cutoff,
                          'window_sec': window_sec,
                         'cadence_sec': cadence_sec}
    return patch_ct, vals


def get_swarm_vals(fname):
    vals, vars = load_swarm(fname)
    # Preliminary calculations
    rad = np.sqrt(np.sum(vals['leo_pos'] ** 2, axis=1))
    unused_alts, vals['lat_mag'], vals['lon_mag'] = physics.transform(rad, vals['lat_geo'] * np.pi / 180, \
                          vals['lon_geo'] * np.pi / 180, from_=['GEO', 'sph'], to=['MAG', 'sph'])
    vals['lat_mag'] *= 180 / np.pi
    vals['lon_mag'] *= 180 / np.pi
    new_vars = 'lat_mag', 'lon_mag'
    vars.update(dict(zip(new_vars, new_vars)))
    return vals, vars


def load_swarm(fname):
    cdf = pycdf.CDF(fname)
    vars = {'Latitude': 'lat_geo',   # Geographic latitude
            'Longitude': 'lon_geo',  # geographic longitude
            'Absolute_STEC': 'tec',  # TEC in TECU
            'Absolute_VTEC': 'vtec',  # Vertical TEC in TECU
            'GPS_Position' : 'gps_pos',  # XYZ position of GPS (m)
            'LEO_Position' : 'leo_pos',  # XYZ position of SWARM (m)
            'Elevation_Angle' : 'elev',  # Elevation angle reported by Swarm
            'Timestamp': 'times',  # Datetime times
            'PRN': 'prn'}  # GPS pseudo-random ID
    vals = {}
    for key, val in vars.items():
       vals[val] = cdf[key][...]

    return vals, vars


def localtime(fname):
    """
    Calculate local time of the Swarm satellites
    """
    vals = load_swarm(fname)
    outvals = {}
    utsec = np.array([(t - dt.datetime(t.year, t.month, t.day)).total_seconds() for t in vals['times']])
    lt = utsec / 3600 + vals['lon_geo'] / 360 * 24 
    lt[lt > 24] -= 24
    lt[lt < 0] += 24
    outvals['lt'] = lt
    outvals['times'] = np.array(vals['times'])
    return outvals 


if __name__ == '__main__':
    main()






