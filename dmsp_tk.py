"""
dmsp_tk.py
Script to perform data analysis of DMSP in situ electron densities

Tailored data-download command for madrigal Python downloader:
globalIsprint.py --verbose --url=http://cedar.openmadrigal.org --parms=YEAR,MONTH,DAY,HOUR,MIN,SEC,MLT,GDALT,GDLAT,GLON,AACGM_LAT,MLAT,AACGM_LONG,MLONG,NE,VERT_ION_V,HOR_ION_V --output=/tmp --user_fullname="Alex+T+Chartier" --user_email=alex.chartier@outlook.com --user_affiliation="APL" --startDate="05/01/2014" --endDate="08/01/2018" --inst=8100 --format=Hdf5 

Bulk download command:
globalDownload.py --verbose --url=http://cedar.openmadrigal.org --outputDir=/home/alex/Downloads/ --user_fullname="Alex+T+Chartier" --user_email=alex.chartier@outlook.com --user_affiliation="APL" --format="hdf5" --startDate="01/01/2014" --endDate="02/01/2014" --inst=8100


"""

import h5py
import numpy as np
import pdb
import datetime as dt
import pickle
import os
import matplotlib.pyplot as plt
import scipy.stats


def main(
        full_pkl_fname=['./data/dmsp/dmsp_%s_', '%Y%m%d_to_', '%Y%m%d.pkl'],
        starttime=dt.datetime(2014, 8, 1),
        endtime=dt.datetime(2018, 8, 1),
        sat=17,
    ):
    full_pkl_fname = full_pkl_fname[0] % sat + starttime.strftime(full_pkl_fname[1]) + endtime.strftime(full_pkl_fname[2])
    """
    try:
        with open(full_pkl_fname, 'rb') as f:
            vals = pickle.load(f)
    except:
    """
    vals = preproc_data()
    pdb.set_trace()
"""
 Things we need: 
    0. Data coverage plot - lat vs MLT/LT 
    1. MLAT/date distribution of velocities - where does the convection pattern extend to? Only look at evening sector
    2. Time series of density at the lowest MLAT where convection happens. Maybe select only the euro-american sector 
"""




def preproc_data(
        sat='17',
        in_fname=['./data/dmsp/dms_%Y%m%d',  '_%ss1.001.hdf5'],
        pkl_fname=['./data/dmsp/dmsp_%s_', '%Y%m%d.pkl'],
        full_pkl_fname=['./data/dmsp/dmsp_%s_', '%Y%m%d_to_', '%Y%m%d.pkl'],
        proc_pkl_fname=['./data/dmsp_proc/dmsp_%s_', '%Y%m%d.pkl'],
        starttime=dt.datetime(2014, 8, 1),
        endtime=dt.datetime(2018, 8, 1),
        timestep=dt.timedelta(days=1),
        needed_vars=['NE', 'VERT_ION_V', 'HOR_ION_V', 'GLON', 'GDALT', 'MLONG', 'MLAT', 'MLT', 'time'],
    ):

    in_fname = in_fname[0] + in_fname[1] % sat
    pkl_fname = pkl_fname[0] % sat + pkl_fname[1]
    full_pkl_fname = full_pkl_fname[0] % sat + starttime.strftime(full_pkl_fname[1]) + endtime.strftime(full_pkl_fname[2])
    
    all_vals = {}
    time = starttime - timestep
    while time < endtime:
        time += timestep
        print(time)
        in_fname_t = time.strftime(in_fname)
        pkl_fname_t = time.strftime(pkl_fname)
        try:
            try: 
                vals = load_pkl(pkl_fname_t)
            except:
                vals = preproc_dmsp(in_fname_t, pkl_fname_t)
        except:
            print(time.strftime('Could not load data for %Y %b %d'))
            continue
    
        pdb.set_trace()

        bn_vals = binned_vals(vals)

        if time == starttime:
            for k in needed_vars:
                all_vals[k] = vals[k]
        else:
            for k in needed_vars:
                all_vals[k] = np.append(all_vals[k], vals[k])



    with open(full_pkl_fname, 'wb') as f:
        print('writing file to full_pklk_fname')
        pickle.dump(all_vals, f) 
    return all_vals 


def binned_vals(vals):
    if type(vals['HOUR']) == np.int64:
        return None
    lt = vals['HOUR'] + vals['MIN'] / 60 + vals['SEC'] / 3600 + vals['GLON'] * 24 / 360
    lt[lt >= 24] -= 24
    lt[lt < 0] += 24

    ev_ind = lt > 12
    morn_ind = lt < 12
    v2 = {}
    v2['evening'] = {k: v[ev_ind] for k, v in vals.items()}
    v2['morning'] = {k: v[morn_ind] for k, v in vals.items()}
    binned_vel = {}
   
    for k, v in v2.items():
        finind_v = np.isfinite(v['HOR_ION_V'])
        #finind_ne = np.isfinite(v['NE'])
        #scipy.stats.binned_statistic(v['GDLAT'][finind_v], v['NE'][finind_ne])
        bin_meds, bin_edges, binnumber = scipy.stats.binned_statistic(v['MLAT'][finind_v], \
                  np.abs(v['HOR_ION_V'][finind_v]), statistic='median', bins=np.linspace(-90, 90, 91))
        """
        plt.hlines(bin_meds, bin_edges[:-1], bin_edges[1:], colors='g', lw=5)
        plt.title(k)
        plt.show()
        """
        pdb.set_trace()
        binned_vel[k] = bin_meds
    binned_vel['bin_edges'] = bin_edges
    return binned_vel
    pdb.set_trace()
    


def load_pkl(in_fname):
    with open(in_fname, 'rb') as f:
        vals = pickle.load(f)
    return vals


def preproc_dmsp(in_fname, out_fname):
    hf = h5py.File(in_fname, 'r')
    vals = {}
    headerdata = hf['Metadata']['Data Parameters'][...]
    header = [h[0].decode('UTF-8') for h in headerdata]
    assert 'NE' in header, 'NE not in header'
    data = hf['Data']['Table Layout'][...]  # 1Hz DMSP data
   
    for h in header:
        vals[h] = np.array([])
    vals['time'] = np.array([])

    for dind, d in enumerate(data):
        for hind, h in enumerate(header):
            vals[h] = np.append(vals[h], d[hind])
        time = dt.datetime(int(vals['YEAR'][dind]), int(vals['MONTH'][dind]), int(vals['DAY'][dind]), \
                           int(vals['HOUR'][dind]), int(vals['MIN'][dind]), int(vals['SEC'][dind]))

        vals['time'] = np.append(vals['time'], time)
    
    with open(out_fname, 'wb') as f:
        pickle.dump(vals, f)
        print('Dumping to %s' % out_fname)
    return vals
        
        
if __name__ == '__main__':
    main() 






















 
