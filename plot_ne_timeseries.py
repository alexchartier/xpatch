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
import proc_swarm_lp
import matplotlib
import matplotlib.dates as mdates

ipath = './data/swarm_lp/'
time = dt.datetime(2016, 5, 8)
sat = 'A'
fname_format = ipath + 'SW_OPER_EFI%s' % sat + '_PL*%Y%m%d*.cdf'
fname = glob.glob(time.strftime(fname_format))[0]

vals = proc_swarm_lp.load_lp(fname)
lt = proc_swarm_lp.localtime(vals)
fig, ax1 = plt.subplots()
utd = mdates.date2num(vals['times'])
plt.plot_date(utd, vals['ne'], '-')
fig.autofmt_xdate()
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
plt.xlabel(r'Time ($UT$)')
plt.ylabel(r'Electron density ($cm^{-3}$)')
plt.grid(which='both')

# plot magnetic latitude ticks
ticks = 20
ax2 = ax1.twiny()
new_tick_locations = utd[0:-1:int(len(utd) / ticks)]
new_ticks = lt[0:-1:int(len(utd) / ticks)]
new_ticklabels = ['%2.1f' % t for t in new_ticks]
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(new_ticklabels)
ax2.set_xlabel(r"Local Time ($Hour$)")
matplotlib.rcParams.update({'font.size': 12})
plt.show()

