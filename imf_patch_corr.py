"""
imf_patch_corr.py
Correlation between IMF and patch occurrence rates
"""
import pdb
import os
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import scipy.stats.stats as st
import sys 
sys.path.append('/Users/chartat1/fusionpp/glimpse/')
import gl_datetime
import plot_patch_ct

starttime = dt.datetime(2016, 1, 1)
endtime = dt.datetime(2017, 1, 1) 
timestep = dt.timedelta(hours=1)
satellites = 'A', 'B', 'C' 
cutoff_crd = 'mag'
fin = './data/proc_lp/%s/' % cutoff_crd + 'lp_%Y%m%d_70deg.pkl'

def main():
    omni_data = read_omni_ascii()
    patch_ct = plot_patch_ct.get_patch_ct(starttime, endtime, satellites, fin)
    correl_sw(omni_data, patch_ct)

def correl_sw(omni_data, patch_ct):
    omni_time = np.array(omni_data['time'])
    timeind = np.logical_and(omni_time >= starttime, omni_time < endtime)
    omni_varnames = 'Bz', 'By', 'B', 'Ey'
    omni = {}
    for var in omni_varnames:
        omni[var] = omni_data[var][timeind]

    bin_edges = []
    t = 0
    while t <= (endtime - starttime).total_seconds():
        bin_edges.append(t) 
        t = t + timestep.total_seconds()
    
    binned_patch_ct = {}
    
    for sat in satellites:
        latarr = np.array(patch_ct[sat]['lat_mag'])
        hems = {'north':  latarr > 0,
                'south': latarr < 0,
                'full': latarr > -95}
        time_ts = np.array([(t[0] - starttime).total_seconds() for t in patch_ct[sat]['times']])
        binned_patch_ct[sat] = {}
        for hem, hemind in hems.items():
            binned_patch_ct[sat][hem] = np.histogram(time_ts[hemind.flatten()], bins=bin_edges)
            for var in omni_varnames:
                corr, p = st.pearsonr(omni[var], binned_patch_ct[sat][hem][0]) 
                if p <= 0.05:
                    print('Sat %s, %s hemisphere, %s correlation: %2.3f, p-factor: %2.3f' % (sat, hem, var, corr, p))
    
    

def read_omni_ascii(fname='./data/omni2_392.lst'):
    with open(os.path.abspath(fname), 'r') as f:
        lines = f.readlines()
        floats = []
        for l in lines:
            floats.append([float(li) for li in l.split()])
        floatarr = np.array(floats)
        varnames = 'year', 'doy', 'hour', 'B', 'By', 'Bz', 'Ey'
        ct = 0
        data = {}
        for v in varnames:
            if (v == 'year') or (v == 'doy'):
                data[v] = floatarr[:, ct].astype(int)
            else:
                data[v] = floatarr[:, ct]
            ct += 1
        data['time'] = []
        for t in range(len(data['year'])):
            data['time'].append(gl_datetime.idadate2event('%i%03d:%02d0000' % (data['year'][t], data['doy'][t], data['hour'][t])))

        return data


if __name__ == '__main__':
    main()
