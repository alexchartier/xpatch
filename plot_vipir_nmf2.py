import pdb 
import datetime as dt
import numpy as np
from plot_lp_radioice import load_hf
import matplotlib.pyplot as plt 
import nvector as nv
import pickle
import sys
sys.path.append('/Users/chartat1/fusionpp/src/nimo')
import nc_utils

starttime = dt.datetime(2019, 2, 28)
step = dt.timedelta(days=1)
endtime = dt.datetime(2019, 3, 14)

in_fname_fmt = '/Users/chartat1/xpatch/data/vipir/Result/result_%Y_%m_%d.nc'
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
        if 'foF2' not in vipir_out.keys():
            vipir_out['foF2'] = np.ma.filled(vipir['foF2'])
        else:
            vipir_out['foF2'] = np.concatenate([vipir_out['foF2'], np.ma.filled(vipir['foF2'])], axis=0)
        time += step

    for k, v in vipir_out.items():
        vipir_out[k] = np.array(v)
    vipir_out['NmF2'] = (vipir_out['foF2'] * 1E6 / 9) ** 2

    return vipir_out


plt.plot(vipir_out['time'], vipir_out['NmF2'], '.')
plt.show()




 
