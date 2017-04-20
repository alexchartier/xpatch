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
import pickle
import sys 
sys.path.append('/Users/chartat1/fusionpp/glimpse/')
import gl_datetime
import plot_patch_ct
import count_passes


def main(starttime=dt.datetime(2016, 1, 1),
           endtime=dt.datetime(2016, 12, 31),
          timestep=dt.timedelta(hours=1),
        satellites=['A', 'B'],
               fin='./data/proc_lp/coley/lp_%Y%m%d_70deg.pkl',
          norm_fin='./data/pass_ct/pass_norm_%Y%m%d.pkl'):

    omni_data = read_omni_ascii()
    omni_filtered, norm_patch_ct = norm_patch_omni(omni_data, fin, norm_fin, starttime, endtime, satellites)
    correl_sw(omni_filtered, norm_patch_ct)


def correl_sw(omni_filtered, norm_patch_ct):
    for sat, v1 in omni_filtered.items():    
        for hem, v2 in v1.items():
            if sat == 'both':
                print('\n\n%s %s hemisphere\n\n' % (sat, hem))
                corrs = {}
                for var in v2.keys():
                    corr, p = st.pearsonr(omni_filtered[sat][hem][var], norm_patch_ct[sat][hem]) 
                    if p <= 0.05:
                        corrstr = '%s correlation: %2.2f, p-factor: %2.3f' % (var, corr, p) 
                        corrs[corrstr] = corr ** 2
                sorted_corr = sorted(corrs.items(), key=operator.itemgetter(1))
                sorted_corr.reverse()
                for sc in sorted_corr:
                    print(sc[0])


def norm_patch_omni(omni, fin, norm_fin, starttime, endtime, satellites):
    # Normalize the patch count according to UT bins. Remove OMNI data where there are no patches
    time = starttime
    timestep = dt.timedelta(days=1)
    omni_filtered = {'A': {'nh': {}, 'sh': {}, 'full': {}},
                     'B': {'nh': {}, 'sh': {}, 'full': {}}}
    norm_patch_ct = {'A': {'nh': [], 'sh': [], 'full': []},
                     'B': {'nh': [], 'sh': [], 'full': []}}
            
    omni_time = np.array(omni['time'])
    omni.pop('time')
    while time < endtime:
        print(time)
        with open(time.strftime(fin), 'rb') as f:
            patch_ct = pickle.load(f)
        with open(time.strftime(norm_fin), 'rb') as f:
            norm_ct = pickle.load(f)
        utbins = np.arange(0, 24.1)
        omni_time_ind = np.logical_and(omni_time >= time, omni_time < time + timestep)
        for sat in satellites:
            try:
                for key, val in patch_ct[sat].items():
                    patch_ct[sat][key] = np.squeeze(np.array(val))
                hems = {'nh': [patch_ct[sat]['lat_mag'] > 0],
                        'sh': [patch_ct[sat]['lat_mag'] < 0],
                      'full': [patch_ct[sat]['lat_mag'] > -100]}
                norm_ct[sat]['full_ut'] = []
                for n in [0, 1]:  # Append the histogram counts, then bins. Can't do it all at once
                    norm_ct[sat]['full_ut'].append(norm_ct[sat]['nh_ut'][n] + norm_ct[sat]['sh_ut'][n])

                for hem, hemind in hems.items():
                    # Normalize patches according to time spent in polar caps
                    ct_ut = np.array([t.hour + t.minute / 60 for t in patch_ct[sat]['times'][hemind]])
                    norm_ut_hist = np.histogram(ct_ut, utbins)[0] / norm_ct[sat][hem + '_ut'][0]
                    assert np.sum(norm_ut_hist == np.inf) == 0, 'If you have infs here, you counted patches in UT bins where the satellite spent 0 time'
                    norm_patch_ct[sat][hem].extend(norm_ut_hist.tolist())

                    # NAN out the omni data where there is no SWARM pass through the polar cap
                    for key, val in omni.items():
                        v_fil = np.asfarray(val[omni_time_ind])
                        v_fil[np.isnan(norm_ut_hist)] = np.nan
                        if time == starttime:
                            omni_filtered[sat][hem][key] = v_fil.tolist()
                        else:
                            omni_filtered[sat][hem][key].extend(v_fil.tolist())
            except:
                # Put in 24 NaNs when you have no Swarm file so the two sats match in length
                for hem, hemind in hems.items():
                    norm_patch_ct[sat][hem].extend((np.ones(24) * np.nan).tolist())  
                    for key, val in omni.items():
                        v_fil = np.asfarray(val[omni_time_ind]) * np.nan
                        omni_filtered[sat][hem][key].extend(v_fil.tolist())
                print('Could not load file for sat %s' % sat)
        time += timestep

    # Sum over both satellites and hemispheres
    norm_patch_ct['both'] = {'full': {}}
    satnorm_ind = np.zeros(len(norm_patch_ct[sat]['full']))
    for sat in satellites:
        ct_sat = np.array(norm_patch_ct[sat]['full'])
        satnorm_ind[np.isfinite(ct_sat)] += 1
        ct_sat_nonan = ct_sat.copy()
        ct_sat_nonan[np.isnan(ct_sat_nonan)] = 0
        if sat == 'A':
            norm_patch_ct['both']['full'] = ct_sat_nonan
        else:
            norm_patch_ct['both']['full'] += ct_sat_nonan
    norm_patch_ct['both']['full'] /= satnorm_ind  # Normalize for the presence of one or both satellites during a given hour
    finite_ind = np.isfinite(norm_patch_ct['both']['full'])

    # remove NaNs from list entries
    norm_patch_ct['both']['full'] = norm_patch_ct['both']['full'][finite_ind]
    omni_filtered['both'] = {'full': {}}
    omni_time_ind = np.logical_and(omni_time >= starttime, omni_time < endtime)
    for key, val in omni.items(): 
        omni_filtered['both']['full'][key] = val[omni_time_ind][finite_ind]

    for sat in satellites:
        for hem in hems.keys():
            fin_ind = np.isfinite(norm_patch_ct[sat][hem])
            norm_patch_ct[sat][hem] = np.array(norm_patch_ct[sat][hem])[fin_ind]

            # Remove the omni data where there is no SWARM pass through the polar cap
            for key, val in omni_filtered[sat][hem].items():
                omni_filtered[sat][hem][key] = np.array(val)[fin_ind]

    return omni_filtered, norm_patch_ct
        

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
