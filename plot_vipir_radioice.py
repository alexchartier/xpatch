import sys
import numpy as np
import pdb
import datetime as dt
import matplotlib.pyplot as plt
from plot_lp_radioice import load_hf
sys.path.append('/Users/chartat1/fusionpp/src/nimo')
sys.path.append('/Users/chartat1/sounder')
import nc_utils
from plot_rtd import plt_frq
import pickle

radiopath='/Users/chartat1/sounder/data/prc_analysis/no_badfrq/spectra/%Y%m%d'
hf_fname_fmt='/Users/chartat1/sounder/data/prc_analysis/no_badfrq/daily/data/%Y%b%d_analysis.pkl'
hf_out_fname_fmt=['./data/hf/hf_%Y%b%d_', '%Y%b%d.pkl']
time=dt.datetime(2019, 2, 28)
step=dt.timedelta(days=1)
endtime = dt.datetime(2019, 3, 14)

r_min = 100
r_max = 600
ob_ang = 25

# Load RadioICE
ri_frq = 7.2
hf = load_hf(time, step, endtime, radiopath, hf_fname_fmt, hf_out_fname_fmt)
spec = hf[ri_frq] 

# Load VIPIR
vipir_fname_fmt = './data/vipir/Amplitude/ampl_%Y_%m_%d.nc'
vipir_out_fname = './data/vipir/out_%1.1f' % ri_frq + time.strftime('%Y%b%d_to') + endtime.strftime('%Y%b%d.pkl')
try:
    with open(vipir_out_fname, 'rb') as f:
        vipir_out = pickle.load(f)
except:
    vipir_out = {}
    while time < endtime:
        vipir_fname = time.strftime(vipir_fname_fmt)
        print('Trying to load %s' % vipir_fname)
        vipir = nc_utils.ncread_vars(vipir_fname)
        if 'time' not in vipir_out.keys():
            vipir_out['time'] = []
        for ind, yr in enumerate(vipir['yr']):
            vipir_out['time'].append(dt.datetime(
                yr, vipir['mon'][ind], vipir['day'][ind], 
                vipir['hr'][ind], vipir['min'][ind], vipir['sec'][ind],
            ))
        vipir_out['alt'] = vipir['tofa']

        # Vipir
        vipir_f = ri_frq * np.cos(np.deg2rad(90 - ob_ang))
        diff = np.abs(vipir['frrq'] - vipir_f)
        f_ind = diff == np.min(diff)
        vipir_out['freq'] = np.squeeze(vipir['frrq'][f_ind])
        if 'ampl' not in vipir_out.keys():
            vipir_out['ampl'] = np.squeeze(vipir['ampl'][:, :, f_ind])
        else:
            vipir_out['ampl'] = np.concatenate([vipir_out['ampl'], np.squeeze(vipir['ampl'][:, :, f_ind])], axis=0)
        time += step
    with open(vipir_out_fname, 'wb') as f:
        pickle.dump(vipir_out, f)

vipir_ampl = np.array(vipir_out['ampl'].astype('float64'))
vipir_ampl[vipir_out['ampl'] < 2000] = np.nan
v_r_ind = np.logical_and(vipir_out['alt'] >= r_min, vipir_out['alt'] <= r_max)
vipir_ampl = vipir_ampl[:, v_r_ind]
vipir_rg = vipir_out['alt'][v_r_ind]

fig, ax = plt.subplots(2, 1, figsize=(8, 6))  

# RadioICE
tlim = time, endtime
im, colorlabel = plt_frq(spec, fig, ax[1], tlim, xlabels=True, vht=True, daylim=True)   
ax[1].set_ylim([r_min, r_max])
ax[1].set_title('%1.1f MHz (%i deg. incidence)' % (ri_frq, ob_ang))

# VIPIR
cmapn = 'gist_heat_r'
cmap = plt.get_cmap(cmapn)
ax[0].pcolormesh(np.array(vipir_out['time']), np.array(vipir_rg), vipir_ampl.T, cmap=cmap, shading='flat',)
ax[0].set_ylim([r_min, r_max])
ax[0].set_ylabel('vHt (km)')
ax[0].set_title('%1.1f MHz (Vertical)' % vipir_out['freq'])
ax[0].grid()

plt.show() 
