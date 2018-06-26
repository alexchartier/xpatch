#!/usr/local/bin/python3
"""
proc_swarm_lp.py
Script to process the SWARM langmuir probe data and analyse for patches. 
"""
import pdb
from spacepy import pycdf
import spacepy.coordinates as crd
import numpy as np
import scipy as sp
import datetime as dt
import matplotlib.pyplot as plt
import glob
import pickle
import sys 
import collections
import aacgmv2

RAD_EARTH = 6371E3


def main(ipath='./data/swarm_lp/',
         opath='./data/swarm/proc_lp/', 
         time=dt.datetime(2017, 1, 1),
         step=dt.timedelta(days=1),
         endtime=dt.datetime(2018, 1, 1),
         lat_cutoff=70,
         sats = ['A', 'B'],
         save=True,
         approach='alex'):
   
    print(approach) 
    while time <= endtime: 
        timestr = time.strftime('%Y-%m-%d')
        print('\n' + timestr)
        vals = {}
        patch_ct = {}    

        for sat in sats:
            print('\nSatellite %s' % sat)
            fname_format = ipath + 'SW_*_EFI%s' % sat + '*%Y%m%d*.cdf' 
            #try:
            fname_str = glob.glob(time.strftime(fname_format))
            if len(fname_str) == 0:
                print('No CDF file matching %s' % time.strftime(fname_format))
            fname = fname_str[0]
            vals[sat] = load_lp(fname)
            vals[sat]['lt'] = localtime(vals[sat])
            if approach == 'coley':
                patch_ct[sat] = coley_patches(vals[sat], lat_cutoff=lat_cutoff)
            elif approach == 'alex':
                patch_ct[sat] = alex_patches(vals[sat], lat_cutoff=lat_cutoff)
            else:
                patch_ct[sat] = count_patches(vals[sat], cutoff_crd=cutoff_crd, lat_cutoff=lat_cutoff)
            #except:
            #    print('Could not count patches for satellite %s on %s' % (sat, timestr))

        if save:
            fout = opath + approach + time.strftime('/lp_%Y%m%d_') + '%ideg.pkl' % lat_cutoff
            with open(fout, 'wb') as f:
                pickle.dump(patch_ct, f) 
            print('\nSaving %s' % fout)
        time += dt.timedelta(days=1)

    return patch_ct, vals


def alex_patches(vals, lat_cutoff=55, window_sec=200, cadence_sec=0.5, filter_pts=30, \
                            edge_mag=1.4, peak_f107_mult=1000, edge_pts=36, cutoff_crd='mag'):
    from pyglow import pyglow
    # Specify cutoff value according to F10.7
    pt = pyglow.Point(vals['times'][0], 100, 0, 0)
    if np.isnan(pt.f107a):
        print('F107a is NaN - stopping')
    assert not np.isnan(pt.f107a), 'F107a is NaN - stopping'
    peak_mag = pt.f107a * peak_f107_mult

    # Count the patches from Langmuir probe data 
    window = dt.timedelta(seconds=window_sec)  
    cadence = dt.timedelta(seconds=cadence_sec) 

    # Transform lats/lons to magnetic 
    alts, vals['lat_mag'], vals['lon_mag'] = transform(vals['rad'], np.deg2rad(vals['lat_geo']), \
                          np.deg2rad(vals['lon_geo']), from_=['GEO', 'sph'], to=['MAG', 'sph'])
    vals['lat_mag'] *= 180 / np.pi
    vals['lon_mag'] *= 180 / np.pi

    # Median-filter the data to remove high frequency noise
    idx = np.arange(filter_pts) + np.arange(len(vals['ne']) - filter_pts + 1)[:, None]
    ne_rm = np.median(vals['ne'][idx], axis=1)
        
    # add random lowlevel noise to the ne_vals to compensate for their low numerical precision for the filter
    ne_rm += (np.random.rand(len(ne_rm)) - 0.5) * 1E-5

    # Shorten all the other datasets to match ne_rm
    for key, val in vals.items():
        vals[key] = val[int(filter_pts / 2): - int(filter_pts / 2) + 1]
    vals['ne_rm'] = ne_rm

    # Reject low-latitude data
    index = (np.abs(vals['lat_' + cutoff_crd]) > lat_cutoff)
    vals_ind = {}
    for key, val in vals.items():  
        vals_ind[key] = vals[key][index]

    # Initialise data storage dictionary
    patch_ct = {}
    for key in vals_ind.keys():
        patch_ct[key] = []

    new_vars = 'ne_rm', 'ne_bg', 't_start', 't_end'
    for v in new_vars:
        patch_ct[v] = []

    # Convert times to integers for faster execution (datetime comparisons are slow)
    times_sec = np.array([(t - vals_ind['times'][0]).total_seconds() for t in vals_ind['times']])
    # Sliding window filtering
    tind = -1
    window_pts = window / cadence
    while tind < len(vals_ind['times']) - window_pts:
        tind += 1  # Has to happen at the top because we use 'continue' to escape the loop at various points
        t = times_sec[tind]
        sys.stdout.write("Time in seconds %s \r" % t)
        ind = np.logical_and(times_sec >= t, times_sec <= t + window_sec)
        sumind = ind.sum()
        # Require full window
        if sumind < window / cadence:  
            # print('Incomplete window found at %s' % t)
            tind += sumind
            continue
        
        times = vals_ind['times'][ind]
        ne_rm_vals = vals_ind['ne_rm'][ind]
        ne_vals = vals_ind['ne'][ind]
        mag_lats = vals_ind['lat_mag'][ind]

        # Check window does not cut across a hemispheric boundary
        lat_steps = np.diff(mag_lats)
        if lat_steps.max() > lat_cutoff:
            print('Found hemispheric jump at %s' % t)
            continue
            
        # Check for 40% increase over 140 km...
        upgrad = 0
        ptind = 0
        while (upgrad < edge_mag) and (ptind < window_pts - edge_pts):
            upgrad = ne_rm_vals[ptind + edge_pts] / ne_rm_vals[ptind] 
            ptind += 1

        if upgrad < edge_mag:
            continue
        upgrad_ind = ptind
        ptind += edge_pts

        # ...followed by 40% decrease over 140 km
        downgrad = 1
        while (downgrad > 1 / edge_mag) and (ptind < window_pts - edge_pts):
            downgrad = ne_rm_vals[ptind + edge_pts] / ne_rm_vals[ptind] 
            ptind += 1
        if downgrad > 1 / edge_mag:
            continue
        downgrad_ind = ptind

        # Specify background density (median over window according to Coley and Heelis)
        NEbg = np.median(ne_vals)

        # Check peak is above threshold
        NEp = ne_rm_vals[upgrad_ind:downgrad_ind].max()  # Patch maximum is the highest value in the window

        # Perform abs magnitude test
        if NEp - NEbg < peak_mag:
            continue

        # If we're still going at this point, we have found a patch. Store the details and skip forward to the next window
        patch_index = vals_ind['ne_rm'] == NEp
        assert patch_index.sum() == 1, 'There should be exactly one patch index for each patch'
        for key, var in vals_ind.items():
            patch_ct[key].append(vals_ind[key][patch_index])
        patch_ct['ne_bg'].append(NEbg)
        patch_ct['t_start'].append(times.min())
        patch_ct['t_end'].append(times.max())
        print('\nFound a patch')
        tind += sumind

    print('4')
    # Count the patches f
    patch_ct['params'] = {'lat_cutoff': lat_cutoff,
                            'peak_mag': peak_mag,
                            'edge_mag': edge_mag,
                            'edge_pts': edge_pts,
                          'filter_pts': filter_pts,
                          'window_sec': window_sec,
                         'cadence_sec': cadence_sec}
    return patch_ct


def coley_patches(vals, lat_cutoff=70, window_sec=165, cadence_sec=0.5, filter_pts=30, \
                            edge_mag=1.4, edge_pts=36, peak_mag=2, cutoff_crd='mag'):
    # Count the patches from Langmuir probe data using Coley and Heelis (1995) approach
    window = dt.timedelta(seconds=window_sec)  
    cadence = dt.timedelta(seconds=cadence_sec) 

    # Transform lats/lons to magnetic 
    alts, vals['lat_mag'], vals['lon_mag'] = transform(vals['rad'], np.deg2rad(vals['lat_geo']), \
                          np.deg2rad(vals['lon_geo']), from_=['GEO', 'sph'], to=['MAG', 'sph'])
    vals['lat_mag'] *= 180 / np.pi
    vals['lon_mag'] *= 180 / np.pi

    # add random lowlevel noise to the ne_vals to compensate for their low numerical precision for the filter
    vals['ne'] += (np.random.rand(len(vals['ne'])) - 0.5) * 1E-5

    # Median-filter the data to remove high frequency noise
    idx = np.arange(filter_pts) + np.arange(len(vals['ne']) - filter_pts + 1)[:, None]
    ne_rm = np.mean(vals['ne'][idx], axis=1)

    # Shorten all the other datasets to match ne_rm
    for key, val in vals.items():
        vals[key] = val[int(filter_pts / 2): - int(filter_pts / 2) + 1]
    vals['ne_rm'] = ne_rm

    # Reject low-latitude data
    index = (np.abs(vals['lat_' + cutoff_crd]) > lat_cutoff)
    vals_ind = {}
    for key, val in vals.items():  
        vals_ind[key] = vals[key][index]

    # Initialise data storage dictionary
    patch_ct = {}
    for key, val in vals_ind.items():
        patch_ct[key] = []

    new_vars = 'ne_rm', 'ne_bg', 't_start', 't_end'
    for v in new_vars:
        patch_ct[v] = []

    # Convert times to integers for faster execution (datetime comparisons are slow)
    times_sec = np.array([(t - vals_ind['times'][0]).total_seconds() for t in vals_ind['times']])
    # Sliding window filtering
    tind = -1
    window_pts = window / cadence
    while tind < len(vals_ind['times']) - window_pts:
        tind += 1  # Has to happen at the top because we use 'continue' to escape the loop at various points
        t = times_sec[tind]
        sys.stdout.write("Time in seconds %s \r" % t)
        ind = np.logical_and(times_sec >= t, times_sec <= t + window_sec)
        sumind = ind.sum()
        # Require full window
        if sumind < window / cadence:  
            # print('Incomplete window found at %s' % t)
            tind += sumind
            continue
        
        times = vals_ind['times'][ind]
        ne_rm_vals = vals_ind['ne_rm'][ind]
        ne_vals = vals_ind['ne'][ind]
        mag_lats = vals_ind['lat_mag'][ind]

        # Check window does not cut across a hemispheric boundary
        lat_steps = np.diff(mag_lats)
        if lat_steps.max() > lat_cutoff:
            print('Found hemispheric jump at %s' % t)
            continue
            
        # Check for 40% increase over 140 km...
        upgrad = 0
        ptind = 0
        while (upgrad < edge_mag) and (ptind < window_pts - edge_pts):
            upgrad = ne_rm_vals[ptind + edge_pts] / ne_rm_vals[ptind] 
            ptind += 1

        if upgrad < edge_mag:
            continue
        upgrad_ind = ptind
        ptind += edge_pts

        # ...followed by 40% decrease over 140 km
        downgrad = 1
        while (downgrad > 1 / edge_mag) and (ptind < window_pts - edge_pts):
            downgrad = ne_rm_vals[ptind + edge_pts] / ne_rm_vals[ptind] 
            ptind += 1
        if downgrad > 1 / edge_mag:
            continue
        downgrad_ind = ptind

        # Specify background density (median over window according to Coley and Heelis)
        NEbg = np.median(ne_vals)

        # Check peak is 2x background
        NEp = ne_rm_vals[upgrad_ind:downgrad_ind].max()  # Patch maximum is the highest value in the window

        # Perform relative magnitude test
        if NEp / NEbg < peak_mag:
            continue

        # If we're still going at this point, we have found a patch. Store the details and skip forward to the next window
        patch_index = vals_ind['ne_rm'] == NEp
        assert patch_index.sum() == 1, 'There should be exactly one patch index for each patch'
        for key, var in vals_ind.items():
            patch_ct[key].append(vals_ind[key][patch_index])
        patch_ct['ne_bg'].append(NEbg)
        patch_ct['t_start'].append(times.min())
        patch_ct['t_end'].append(times.max())
        # print('\nFound a patch')
        tind += sumind

    # Store the count parameters
    patch_ct['params'] = {'lat_cutoff': lat_cutoff,
                            'peak_mag': peak_mag,
                            'edge_mag': edge_mag,
                            'edge_pts': edge_pts,
                          'filter_pts': filter_pts,
                          'window_sec': window_sec,
                         'cadence_sec': cadence_sec}
    return patch_ct


def count_patches(vals, lat_cutoff=70, window_sec=200, min_time_sec=10, cadence_sec=0.5, rel_mag_cutoff=2, cutoff_crd='mag'):
    # Count the patches from Langmuir probe data
    # Patch = 2x background density (defined by ne_fac) over 78 < d < 1560 km. Translates to 10 < t < 200s
    window = dt.timedelta(seconds=window_sec)  
    cadence = dt.timedelta(seconds=cadence_sec) 

    # Transform lats/lons to magnetic 
    alts, vals['lat_mag'], vals['lon_mag'] = physics.transform(vals['rad'], np.deg2rad(vals['lat_geo']), \
                          np.deg2rad(vals['lon_geo']), from_=['GEO', 'sph'], to=['MAG', 'sph'])
    vals['lat_mag'] *= 180 / np.pi
    vals['lon_mag'] *= 180 / np.pi
        
    # add random lowlevel noise to the ne_vals to compensate for their low numerical precision for the filter
    vals['ne'] += (np.random.rand(len(vals['ne'])) - 0.5) * 1E-5

    # Reject low-latitude data
    index = (np.abs(vals['lat_' + cutoff_crd]) > lat_cutoff)
    vals_ind = {}
    for key, val in vals.items():  
        vals_ind[key] = vals[key][index]

    # Initialise data storage dictionary
    patch_ct = {}
    for key, val in vals_ind.items():
        patch_ct[key] = []

    new_vars = 'ne_bg', 'ne_b1', 'ne_b2', 't1', 't2', 't_start', 't_end'
    for v in new_vars:
        patch_ct[v] = []

    # Convert times to integers for faster execution (datetime comparisons are slow)
    times_sec = np.array([(t - vals_ind['times'][0]).total_seconds() for t in vals_ind['times']])
    # Sliding window filtering
    tind = -1
    while tind < len(vals_ind['times']) - window / cadence:
        # NOTES FOR PUBLICATION: 
        #   Not clear if Noja et al. skipped forward if they found a patch within the window. Assume they did
        #   What happens if there are less than max. points available in a segment? We throw the segment out.
        tind += 1  # Has to happen at the top because we use 'continue' to escape the loop at various points
        t = times_sec[tind]
        sys.stdout.write("Time in seconds %s \r" % t)
        ind = np.logical_and(times_sec >= t, times_sec <= t + window_sec)
        sumind = ind.sum()
        # Require full window
        if sumind < window / cadence:  
            # print('Incomplete window found at %s' % t)
            tind += sumind
            continue
        
        times = vals_ind['times'][ind]
        ne_vals = vals_ind['ne'][ind]
        mag_lats = vals_ind['lat_mag'][ind]
        grads = np.diff(ne_vals)

        # Check window does not cut across a hemispheric boundary
        lat_steps = np.diff(mag_lats)
        if lat_steps.max() > lat_cutoff:
            print('Found hemispheric jump at %s' % t)
            continue
            
        # Algorithm requires a positive gradient ...
        if grads.max() <= 0:
            continue
        pos_ind = np.argmax(grads >= 0)
        # ... followed by a negative gradient
        if grads[pos_ind:].min() >= 0:
            continue

        NEp = ne_vals.max()  # Patch maximum is the highest value in the window
        ind_NEp = np.where(ne_vals == NEp)[0][0]
        
        # The next part won't work if the first/last value is the largest
        if (ind_NEp == 0) or (ind_NEp == len(ne_vals) - 1):
            continue

        # Define the two boundary values: b1 - greater of the two minimum values either side
        #                                 b2 - closest value on other side to b1     
        lhs = ne_vals[:ind_NEp]
        rhs = ne_vals[ind_NEp + 1:]
        NE_b1 = max([np.min(lhs), np.min(rhs)])
        if NE_b1 in lhs:
            NE_b2 = rhs[(np.abs(rhs - NE_b1)).argmin()]
        else:
            NE_b2 = lhs[(np.abs(lhs - NE_b1)).argmin()]

        # Determine background NE: Symmetric linear interpolation of values to location of the peak
        ind_NE_b1 = np.where(ne_vals == NE_b1)[0][0]
        ind_NE_b2 = np.where(ne_vals == NE_b2)[0][0]
        assert (ind_NE_b1 > ind_NEp) ^ (ind_NE_b2 > ind_NEp), 'Background indices should be either side of peak'
        time_b1 = (times[ind_NE_b1] - times.min()).seconds
        time_b2 = (times[ind_NE_b2] - times.min()).seconds
        time_p = (times[ind_NEp] - times.min()).seconds
    
        # Interpolate background values to location of patch    
        try:
            NEbg = sp.interpolate.interp1d([time_b1, time_b2], [NE_b1, NE_b2])(time_p).tolist()
        except:
            pdb.set_trace()
    
        # Perform relative magnitude test
        if NEp / NEbg < rel_mag_cutoff:
            continue
        
        if abs(time_b1 - time_b2) < min_time_sec:
            print('found short isolated spike')
            continue

        # If we're still going at this point, we have found a patch. Store the details and skip forward to the next window
        patch_index = vals_ind['ne'] == NEp
        assert patch_index.sum() == 1, 'There should be exactly one patch index for each patch'
        for key, var in vals_ind.items():
            patch_ct[key].append(vals_ind[key][patch_index])
        patch_ct['ne_bg'].append(NEbg)
        patch_ct['ne_b1'].append(NE_b1)
        patch_ct['ne_b2'].append(NE_b2)
        patch_ct['t1'].append(times[ind_NE_b1])
        patch_ct['t2'].append(times[ind_NE_b2])
        patch_ct['t_start'].append(times.min())
        patch_ct['t_end'].append(times.max())
        # print('\nFound a patch')
        tind += sumind

    # newvars = 'ne_bg', 'ne_b1', 'ne_b2', 't1', 't2', 't_start', 't_end'
    # for v in newvars:
    #     patch_ct[v] = np.array(patch_ct[v])
    patch_ct['params'] = {'lat_cutoff': lat_cutoff,
                      'rel_mag_cutoff': rel_mag_cutoff,
                          'window_sec': window_sec,
                        'min_time_sec': min_time_sec,
                         'cadence_sec': cadence_sec}
    return patch_ct


def load_lp(fname):
    """
    Load the Swarm langmuir probe data
    """ 
    cdf = pycdf.CDF(fname)
    vars = {'Latitude': 'lat_geo',   # Geographic latitude
            'Longitude': 'lon_geo',  # geographic longitude
            'Diplat': 'Diplat', 
            'Diplon': 'Diplon', 
            'MLat': 'AACGMLat', 
            'MLon': 'AACGMLon', 
            'Radius': 'rad',  # Radial distance
            'n': 'ne',  # Electron density
            'Timestamp': 'times',  # Datetime times
            'T_elec': 'T_elec',  # Electron temp.
            }   
    vals = {}
        
    for key, val in vars.items():
        try:
            vals[val] = cdf[key][...]
        except:
            None

    return vals


def localtime(vals):
    """
    Calculate local time of the Swarm satellites
    """
    outvals = {}
    utsec = np.array([(t - dt.datetime(t.year, t.month, t.day)).total_seconds() for t in vals['times']])
    lt = utsec / 3600 + vals['lon_geo'] / 360 * 24 
    lt[lt > 24] -= 24
    lt[lt < 0] += 24
    return lt 


def transform(rads, lats, lons, from_=['GEO', 'sph'], to=['MAG', 'sph']):
    # Convert 3D arrays (or vectors) from geographic to geomagnetic coordinates
    # Function takes radians, meters; returns lats, lons
    rads = np.array(rads)
    lats = np.array(lats) 
    lons = np.array(lons) 
    assert ((rads.size == lats.size) and (rads.size == lons.size)), 'Inputs must have the same number of elements'
    in_crd_array = np.array([np.ravel(rads / RAD_EARTH).tolist(), \
                           np.ravel(np.rad2deg(lats)).tolist(), \
                           np.ravel(np.rad2deg(lons)).tolist()]).T.tolist()
   
    in_crd_array = np.array(in_crd_array)
    pdb.set_trace() 
    aacgmv2.convert_latlon_arr(in_crd_array[:, 1], in_crd_array[:, 2], 100, dt.datetime(2015, 1, 1))
    # in_crd = crd.Coords(in_crd_array, from_[0], from_[1])
    in_crd.ticks = Ticktock(np.tile(['2014-12-31T12:00:00'], len(in_crd_array)), 'ISO')  # Needs a time to specify mag. field
    out_crd = in_crd.convert(to[0], to[1]) 
    out_alts = np.reshape(np.array(out_crd.data[:, 0]), rads.shape)
    out_lats = np.reshape(np.array(out_crd.data[:, 1]), lats.shape)
    out_lons = np.reshape(np.array(out_crd.data[:, 2]), lons.shape)
    if to[1] == 'sph':
        out_lats = np.deg2rad(out_lats)
        out_lons = np.deg2rad(out_lons)
    return out_alts, out_lats, out_lons


if __name__ == '__main__':
    main()






