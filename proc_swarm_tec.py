#!/Library/Frameworks/Python.framework/Versions/3.4/bin/python3
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

def main():
    ipath = '/Volumes/Seagate/data/swarm/gps_tec/'
    opath = '/Volumes/Seagate/data/swarm/lt/'
    time = dt.datetime(2014, 10, 20)
    step = dt.timedelta(days=1)
    endtime = dt.datetime(2014, 12, 31)
    
    while time < endtime: 
        timestr = time.strftime('%Y-%m-%d')
        print(timestr)
        # sys.stdout.write("%s\r" % timestr) 
        sats = 'A', 'B', 'C'
        patch_ct = {}
        vals = {}

        for sat in sats:
            print('Satellite %s' % sat)
            fname_format = ipath + 'SW_OPER_TEC%s' % sat + '*%Y%m%d*.DBL'
            try:
                fname = glob.glob(time.strftime(fname_format))[0]
                patch_ct[sat] = count_patches(fname)
            except:
                print('No file for satellite %s on %s' % (sat, timestr))
        fout = opath + time.strftime('lt_%Y%m%d.pkl')
        with open(fout, 'wb') as f:
            pickle.dump(patch_ct, f) 
        print('Saving %s' % fout)
        time += dt.timedelta(days=1)

def count_patches(fname, lat_cutoff=55, elev_cutoff=25, TEC_abs_cutoff=4, TEC_rel_cutoff=1.2, window_sec=200, cadence_sec=10):
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
    vals = load_swarm(fname)
    # Preliminary calculations
    rad = np.sqrt(np.sum(vals['leo_pos'] ** 2, axis=1))
    vals['lat_mag'], vals['lon_mag'] = physics.transform(rad, vals['lat_geo'] * np.pi / 180, \
                          vals['lon_geo'] * np.pi / 180, from_=['GEO', 'sph'], to=['MAG', 'sph'])
    vals['lat_mag'] *= 180 / np.pi
    vals['lon_mag'] *= 180 / np.pi
    vals['elev'] = physics.elevation(vals['gps_pos'], vals['leo_pos']) * 180 / np.pi
    new_vars = 'lat_mag', 'lon_mag', 'elev'
    vars.update(dict(zip(new_vars, new_vars)))

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
    patch_ct['tec_bg'] = []

    for p in unique_prns:
        index = vals['prn'] == p
        p = str(p)
        vals_p[p] = {}
        for key, val in vars.items():
            vals_p[p][val] = vals[val][index]

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
            mag_lats = vals['lat_mag'][ind]

            # Problem with duplicate values in SWARM data. This will fix the problems it causes.
            duplicates = [tec_val for tec_val, count in collections.Counter(tec_vals).items() if count > 1]
            for duplicate in duplicates:
                print('Found duplicate TEC value in PRN %s on %s' % (p, t))
                tec_vals[np.argmax(tec_vals == duplicate)] += 1E-9
                vals_p[p]['tec'][np.argmax(vals_p[p]['tec'] == duplicate)] += 1E-9

            # Check window does not cut across a hemispheric boundary
            lat_steps = np.diff(mag_lats)
            if lat_steps.max() > lat_cutoff:
                print('Found hemispheric jump at %s' % t)
                vals['times'].index(times[lat_steps == lat_steps.max()])
                pdb.set_trace()
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
            tind += sum(ind)

    patch_ct['params'] = {'lat_cutoff': lat_cutoff,
                         'elev_cutoff': elev_cutoff,
                      'TEC_abs_cutoff': TEC_abs_cutoff,
                      'TEC_rel_cutoff': TEC_rel_cutoff,
                          'window_sec': window_sec,
                         'cadence_sec': cadence_sec}
    return patch_ct

def load_swarm(fname):
    cdf = pycdf.CDF(fname)
    vars = {'Latitude': 'lat_geo',   # Geographic latitude
            'Longitude': 'lon_geo',  # geographic longitude
            'Absolute_STEC': 'tec',  # TEC in TECU
            'GPS_Position' : 'gps_pos',  # XYZ position of GPS (m)
            'LEO_Position' : 'leo_pos',  # XYZ position of SWARM (m)
            'Timestamp': 'times',  # Datetime times
            'PRN': 'prn'}  # GPS pseudo-random ID
    vals = {}
    for key, val in vars.items():
       vals[val] = cdf[key][...]

    return vals


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






