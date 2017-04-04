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
import operator
import sys 
sys.path.append('/Users/chartat1/fusionpp/glimpse/')
import gl_datetime
import plot_patch_ct
import count_passes

starttime = dt.datetime(2016, 1, 1)
endtime = dt.datetime(2017, 1, 1) 
timestep = dt.timedelta(hours=1)
satellites = 'A', 'B', 'C' 
cutoff_crd = 'mag'
fin = './data/proc_lp/%s/' % cutoff_crd + 'lp_%Y%m%d_70deg.pkl'

def main():
    pass_ct = count_passes.get_pass_ct()
    omni_data = read_omni_ascii()
    patch_ct = plot_patch_ct.get_patch_ct(starttime, endtime, satellites, fin)
    correl_sw(omni_data, patch_ct, pass_ct)

def correl_sw(omni_data, patch_ct, pass_ct):
    omni_time = np.array(omni_data['time'])
    timeind = np.logical_and(omni_time >= starttime, omni_time < endtime)
    omni_varnames = set(omni_data.keys()) - set(['DOY', 'time', 'Hour'])
    omni = {}
    for var in omni_varnames:
        omni[var] = omni_data[var][timeind]

    bin_edges = []
    t = 0
    while t <= (endtime - starttime).total_seconds():
        bin_edges.append(t) 
        t = t + timestep.total_seconds()
    
    binned_patch_ct = {}
    norm_patch_ct = {}
    binned_pass_ct = {}
   
    # Sum up all the satellites 
    for sat in satellites:
        latarr = np.array(patch_ct[sat]['lat_mag'])
        hems = {'north': latarr > 0,
                'south': latarr < 0,
                 'full': latarr > -95}
        pass_hems = {'north': pass_ct[sat]['hem'] > 0,
                     'south': pass_ct[sat]['hem'] < 0,
                      'full': pass_ct[sat]['hem'] > -95}
        time_ts = np.array([(t[0] - starttime).total_seconds() for t in patch_ct[sat]['times']])
        pass_time_ts = np.array([(t - starttime).total_seconds() for t in pass_ct[sat]['times']])
        binned_patch_ct[sat] = {}
        norm_patch_ct[sat] = {}
        binned_pass_ct[sat] = {}
        for hem, hemind in hems.items():
            pass_hemind = pass_hems[hem].flatten()
            binned_patch_ct[sat][hem] = np.histogram(time_ts[hemind.flatten()], bins=bin_edges)
            binned_pass_ct[sat][hem] = np.histogram(pass_time_ts[pass_hemind], bins=bin_edges)
            norm_patch_ct[sat][hem] = binned_patch_ct[sat][hem][0] / binned_pass_ct[sat][hem][0]
            try:
                norm_patch_ct[hem] += norm_patch_ct[sat][hem]
            except:
                norm_patch_ct[hem] = norm_patch_ct[sat][hem]

    for hem in hems.keys():
        print('\n\n%s hemisphere\n\n' % hem)
        corrs = {}
        finind = np.isfinite(norm_patch_ct[hem])
        for var in omni_varnames:
            corr, p = st.pearsonr(omni[var][finind], norm_patch_ct[hem][finind]) 
            if p <= 0.05:
                corrstr = '%s correlation: %2.2f, R-squared: %2.2f, p-factor: %2.3f' % (var, corr, corr ** 2, p) 
                corrs[corrstr] = corr ** 2
        sorted_corr = sorted(corrs.items(), key=operator.itemgetter(1))
        sorted_corr.reverse()
        for sc in sorted_corr:
            print(sc[0])
    pdb.set_trace()


def read_omni_ascii(fname='./data/omni2_full.lst', format='./data/omni2_full.fmt'):
    with open(os.path.abspath(format), 'r') as f:
        lines = f.readlines()[4:-2]
        varnames = []
        types = []
        for line in lines:
            varnames.append(' '.join(line.split()[1:-1]))
            types.append(line.split()[-1])
    with open(os.path.abspath(fname), 'r') as f:
        lines = f.readlines()
        floats = []
        for l in lines:
            floats.append([float(li) for li in l.split()])
        floatarr = np.array(floats)
        ct = 0
        data = {}
        for vind, v in enumerate(varnames):
            if 'I' in types[vind]:
                data[v] = floatarr[:, ct].astype(int)
            else:
                data[v] = floatarr[:, ct]
            ct += 1
        data['time'] = []
        for t in range(len(data['YEAR'])):
            data['time'].append(gl_datetime.idadate2event('%i%03d:%02d0000' % (data['YEAR'][t], data['DOY'][t], data['Hour'][t])))

        for varname in varnames:
            if len(np.unique(data[varname])) < 20:
                data.pop(varname)
                print('Removed %s due to lack of unique values' % varname)

        return data


if __name__ == '__main__':
    main()
