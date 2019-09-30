"""
plt_hf_rg_spread.py
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
import matplotlib.pyplot as plt
import matplotlib
import pickle

def main(
    starttime = dt.datetime(2019, 2, 28) ,
    step = dt.timedelta(minutes=1),
    endtime = dt.datetime(2019, 3, 14) ,
    radiopath = '/Users/chartat1/sounder/data/prc_analysis/no_badfrq/spectra/%Y%m%d',
    hf_fname_fmt = '/Users/chartat1/sounder/data/prc_analysis/no_badfrq/daily/data/%Y%b%d_analysis.pkl',
    hf_out_fname_fmt = ['./data/hf/hf_%Y%b%d_', '%Y%b%d.pkl'],
    freq=5.1,
):

    hf = load_hf(starttime, dt.timedelta(days=1), endtime, radiopath, hf_fname_fmt, hf_out_fname_fmt)
    hff =  hf[freq]
    alt = np.tile(hf[5.1]['alt'], (hff['time'].shape[0], 1))
    alt[hff['max_pwr_db'] == 0] *= np.nan
    alt[alt < 200] *= np.nan
    alt[alt > 800] *= np.nan
    rg_spread = np.nanmax(alt, 1) - np.nanmin(alt, 1)

    ax = plt.subplot(111)
    ax.bar(hff['time'], rg_spread, width=0.1, color='k')
    ax.xaxis_date()
    ax.set_xlim(starttime, endtime)
    ax.set_ylabel('F-region virt. ht spread @ %1.1f MHz (km)' % freq)
    ax.grid(which='both')
    matplotlib.rc('font', size=16)
    plt.show() 



if __name__ == "__main__":
    main()

