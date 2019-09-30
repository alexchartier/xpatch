"""
calc_geophys_params.py
Get NmF2 out of the sounder data
""" 
import pdb
import datetime as dt
import numpy as np
from plot_lp_radioice import load_hf
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import nvector as nv
import sys
sys.path.append('/Users/chartat1/fusionpp/src/nimo')
sys.path.append('/Users/chartat1/sounder')
import nc_utils
import pickle

def main(
    radiopath = '/Users/chartat1/sounder/data/prc_analysis/no_badfrq/spectra/%Y%m%d',
    starttime = dt.datetime(2019, 2, 28) ,
    step = dt.timedelta(minutes=1),
    endtime = dt.datetime(2019, 3, 14) ,
    hf_fname_fmt = '/Users/chartat1/sounder/data/prc_analysis/no_badfrq/daily/data/%Y%b%d_analysis.pkl',
    hf_out_fname_fmt = ['./data/hf/hf_%Y%b%d_', '%Y%b%d.pkl'],
    vipir_fname_fmt = '/Users/chartat1/xpatch/data/vipir/Result/result_%Y_%m_%d.nc',
    param_out_fname_fmt = ['./data/hf/hf_params_%Y%b%d_', '%Y%b%d.pkl'],
):

    hf = load_hf(starttime, dt.timedelta(days=1), endtime, radiopath, hf_fname_fmt, hf_out_fname_fmt)

    times = []
    time = starttime 
    while time <= endtime:
       times.append(time)
       time += step 
    times = np.array(times)

    #out_fname = starttime.strftime(param_out_fname_fmt[0]) + endtime.strftime(param_out_fname_fmt[1]) 
    #try: 
    #    with open(out_fname, 'rb') as f:
    #        hf = pickle.load(f)
    #except:
    hf = get_hf_params(times, hf)
    #    with open(out_fname, 'wb') as f:
    #        pickle.dump(hf, f)

    vipir = load_vipir_nmf2(starttime, endtime, dt.timedelta(days=1), vipir_fname_fmt)
    
    plts = [
        ['hmE', 'km',  180],
        ['NmE', 'el./m3', 2E11],
    ]
    plot(hf, vipir, plts, starttime, endtime)
    
    plts = [
        ['hmF2', 'km', 600],
        ['NmF2', 'el./m3', 5E11],
    ]
    plot(hf, vipir, plts, starttime, endtime)

def load_vipir_nmf2(starttime, endtime, step, in_fname_fmt):
    time = starttime
    vipir_out = {}
    while time <= endtime:
        vipir_fname =  time.strftime(in_fname_fmt)
        print('Trying to load %s' % vipir_fname)
        vipir = nc_utils.ncread_vars(vipir_fname)
        if 'time' not in vipir_out.keys():
            vipir_out['time'] = []
        for ind, yr in enumerate(vipir['yr']):
            vipir_out['time'].append(dt.datetime(
                yr, vipir['mon'][ind], vipir['day'][ind],
                vipir['hr'][ind], vipir['min'][ind], vipir['sec'][ind],
            ))
        kv = {
            'foF2': 'foF2',
            'hmF2': 'HMXF',
            'foE': 'foE',
            'hmE': 'HMXE',
        }
            
        if 'foF2' not in vipir_out.keys():
            for k, v in kv.items():
                vipir_out[k] = np.ma.filled(vipir[v])
        else:
            for k, v in kv.items():
                vipir_out[k] = np.concatenate([vipir_out[k], np.ma.filled(vipir[v])], axis=0)
        time += step

    for k, v in vipir_out.items():
        vipir_out[k] = np.array(v)
    vipir_out['NmF2'] = (vipir_out['foF2'] * 1E6 / 9) ** 2
    vipir_out['NmE'] = (vipir_out['foE'] * 1E6 / 9) ** 2
   
    e_badind = np.logical_or(vipir_out['hmE'] < 50, vipir_out['NmE'] < 1E3)
    f_badind = np.logical_or(vipir_out['hmF2'] < 50, vipir_out['NmF2'] < 1E3) 
    vipir_out['hmF2'][f_badind] *= np.nan
    vipir_out['NmF2'][f_badind] *= np.nan
    vipir_out['hmE'][e_badind] *= np.nan
    vipir_out['NmE'][e_badind] *= np.nan

    return vipir_out


def get_hf_params(times, hf):
    # Get record of MOF and virt. ht
    wgs84 = nv.FrameE(name='WGS84')
    mcm = wgs84.GeoPoint(latitude=-77.8564, longitude=166.6881, z=0, degrees=True).to_ecef_vector().pvector
    zsp = wgs84.GeoPoint(latitude=-90.0, longitude=166.6881, z=0, degrees=True).to_ecef_vector().pvector
    mcm_zsp_strt = np.sqrt(np.sum((mcm - zsp) ** 2)) / 1E3
    
    hmf2, rg_hmf2, MOF = get_hm_rg_mof(hf, times, 180, 600)
    hmE, rg_hmE, MOF_E = get_hm_rg_mof(hf, times, 80, 160)
    
    # Get NmF2 from MUF and virt. ht
    incid = np.arcsin(mcm_zsp_strt / rg_hmf2)  * 0.9
    crit_frq = MOF * np.cos(incid)
    NmF2 = (crit_frq * 1E6 / 9) ** 2

    # Get NmF2 from MUF and virt. ht
    incid_E = np.arcsin(mcm_zsp_strt / rg_hmE) * 0.9
    crit_frq_E = MOF_E * np.cos(incid_E)
    NmE = (crit_frq_E * 1E6 / 9) ** 2

    out = {
        'hmE': hmE,
        'NmE': NmE,
        'hmF2': hmf2,
        'NmF2': NmF2,
        'f0f2': crit_frq,
        'MOF': MOF,
        'time': times,
    }
    return out


def get_hm_rg_mof(hf, times, alt_min, alt_max):
    MOF = np.zeros(times.shape)
    hm = np.zeros(times.shape)
    rg_hm = np.zeros(times.shape)
    for freq, vals in hf.items():
        hf_times = vals['time']  
        alts = vals['alt']
        ranges = vals['range']
        altind = np.logical_and(alts > alt_min, alts < alt_max)
        
        for hf_tind, hf_t in enumerate(hf_times):
            td = times - hf_t
            td = np.abs(np.array([t.total_seconds() for t in td]))
            t_ind = td == np.min(td)

            if MOF[t_ind][0] < freq:
                if np.sum(vals['max_pwr_db'][hf_tind, altind]) > 0:
                    MOF[t_ind] = freq 
                    pwr = vals['max_pwr_db'][hf_tind, :]
                    pwr_ind = np.logical_and(pwr > 0, altind)
                    hm[t_ind] = alts[pwr_ind].min()
                    rg_hm[t_ind] = ranges[pwr_ind].min()

    lowind = hm < 5 
    hm[lowind] *= np.nan
    rg_hm[lowind] *= np.nan
    MOF[lowind] *= np.nan
    return hm, rg_hm, MOF


def plot(hf, vipir, plts, starttime, endtime):
    plt.rcParams.update({'font.size': 14})
    fig, ax = plt.subplots(len(plts), 1, sharex=True, sharey=False, figsize=(8, 6))
    for ind, prm in enumerate(plts):
        ax[ind].plot(vipir['time'], vipir[prm[0]], 'b.', label='VIPIR Jang Bogo')
        ax[ind].plot(hf['time'], hf[prm[0]], 'r.', label='Oblique MCM-ZSP')
        ax[ind].set_xlim(starttime, endtime)
        ax[ind].set_ylim(0, prm[2])
        ax[ind].grid(which='both')   
        ax[ind].set_ylabel('%s  (%s)' % (prm[0], prm[1]))
        if ind == len(plts) - 1:
            ax[ind].legend(loc='best')
            ax[ind].set_xlabel('Time  (UT)')
            ax[ind].xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    plt.show()




if __name__ == '__main__':
    main()

