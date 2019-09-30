# O/N2 plotter

import pdb
from pyglow import pyglow
import datetime as dt
import pickle
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy.interpolate import griddata, interp2d

# Set inputs
times_l = []
for start_d in np.arange(1, 31):
    for start_t in [dt.datetime(2014, 1, start_d), dt.datetime(2014, 7, start_d)]:
        times_l.append(
            [start_t + dt.timedelta(hours=x) for x in range(0, 24, 3)],
        )

d_alt = 10
alts = np.arange(200, 403, d_alt)
lats = np.arange(-90, 95, 5)
lons = np.arange(0, 360, 30)
pkl_fname = 'data/on2/O_N2_%Y%m%d_%H%M.pkl'
recalc = False
lat_2d, lon_2d = np.meshgrid(lats, lons)
"""
# Coordinates calculation (requires python3)
from apexpy import Apex
date = 2014.0
alt = 300 
A = Apex(date=date)
m_crd = {
    'lat': [],
    'lon': [],
}
for lat, lon in zip(lat_2d.flatten(), lon_2d.flatten()):
    mlat, mlon = A.convert(lat, lon, 'geo', 'apex', height=alt)
    m_crd['lat'].append(mlat)
    m_crd['lon'].append(mlon)
for k, v in m_crd.items():
    m_crd[k] = np.array(v)

with open('crd.pkl', 'wb') as f:
    pickle.dump(m_crd, f, protocol=2)
"""
with open('crd.pkl', 'rb') as f:
    m_crd = pickle.load(f)

# Calc O/N2
if recalc:
    for times in times_l:
        O = np.zeros((len(alts), len(lats), len(lons))) * np.nan
        N2 = np.zeros((len(alts), len(lats), len(lons))) * np.nan
        Tn = np.zeros((len(alts), len(lats), len(lons))) * np.nan
        Ti = np.zeros((len(alts), len(lats), len(lons))) * np.nan
        for timeind, time in enumerate(times):
            print(time)
            out = {}
            for altind, alt in enumerate(alts):
                for latind, lat in enumerate(lats):
                    for lonind, lon in enumerate(lons):
                        pt = pyglow.Point(time, lat, lon, alt)
                        pt.run_msis()
                        pt.run_iri()
                        O[altind, latind, lonind] = pt.nn['O']
                        N2[altind, latind, lonind] = pt.nn['N2']
                        Tn[altind, latind, lonind] = pt.Tn_msis
                        Ti[altind, latind, lonind] = pt.Ti
            out['O'] = O
            out['N2'] = N2
            out['Tn'] = Tn
            out['Ti'] = Ti
            out['alt'] = alts
            out['lat'] = lat
            out['lon'] = lon
            out['time'] = time
            with open(time.strftime(pkl_fname), 'wb') as f:
                pickle.dump(out, f)

# Load into a holder
ct = 0
species = ['O', 'N2']
color = 'g', 'b'

fig, ax  = plt.subplots(2, 1, sharex=True)
Tn_NSm = np.zeros((len(times_l), 2))
for mind, times in enumerate(times_l):
    nden = {}
    for spec in species:
        nden[spec]  = np.zeros((len(times), len(lats), len(lons)))
    Tn = np.zeros((len(times), len(lats)))
    Ti = np.zeros((len(times), len(lats)))
    ON2 = np.zeros((len(times), len(lats)))
    N2 = np.zeros((len(times), len(lats)))
    Tn_NS = np.zeros((len(times), 2)) 
    for timeind, time in enumerate(times):
        with open(time.strftime(pkl_fname), 'rb') as f:
            out = pickle.load(f)
    
        # Column-integrated density of O and N2
        amin = 200
        altind = alts >= amin
        #for spec in ['O', 'N2']:
            # ON2[timeind, :] = np.mean(np.sum(out['O'][altind, :, :], 0) / np.sum(out['N2'][altind, :, :], 0), 1)
        #N2[timeind, :] = np.mean(np.mean(out['N2'][altind, :, :], 0), 1)
        #Ti[timeind, :] = np.mean(np.mean(out['Ti'][altind, :, :], 0), 1)
        Tn_NS[timeind, 0] = np.mean(out['Tn'][-1, :, :].flatten()[m_crd['lat'] < -70])
        Tn_NS[timeind, 1] = np.mean(out['Tn'][-1, :, :].flatten()[m_crd['lat'] > 70])
        

    month = time.strftime('%b')
    Tn_NSm[mind, :] = np.mean(Tn_NS, 0)
    print('%s\n South    North' % month)
    print('%2.1f K  %2.1f K' % (Tn_NSm[mind, 0], Tn_NSm[mind, 1]))

print('******\n\n')
tarr = np.array(times_l)
janind = (tarr < dt.datetime(2014, 6, 1))[:, 0]
julind = (tarr > dt.datetime(2014, 6, 1))[:, 0]
pdb.set_trace()
print('Jan\n%2.1f K  %2.1f K' % (np.mean(Tn_NSm[janind, 0]), np.mean(Tn_NSm[janind, 1])))
print('Jul\n%2.1f K  %2.1f K' % (np.mean(Tn_NSm[julind, 0]), np.mean(Tn_NSm[julind, 1])))
"""
    #ax[0].plot(lats, np.mean(Ti, axis=0), '-%s' % color[mind], label=month)
    #ax[0].plot(lats, np.max(Ti, axis=0),'--%s' % color[mind])
    #ax[0].plot(lats, np.min(Ti, axis=0),'--%s' % color[mind])
    #ax[0].set_ylabel('Mean Ion Temp. from %i-400km (K)' % amin)
    ax[0].plot(lats, np.mean(Tn, axis=0), '-%s' % color[mind], label=month)
    ax[0].plot(lats, np.max(Tn, axis=0),'--%s' % color[mind])
    ax[0].plot(lats, np.min(Tn, axis=0),'--%s' % color[mind])
    ax[0].set_ylabel('Neutral Temp. @ 400km (K)')
    ax[1].plot(lats, np.mean(N2, axis=0), '-%s' % color[mind], label=month)
    ax[1].plot(lats, np.max(N2, axis=0),'--%s' % color[mind])
    ax[1].plot(lats, np.min(N2, axis=0),'--%s' % color[mind])
    #ax[1].plot(lats, np.mean(ON2, axis=0), '-%s' % color[mind], label=month)
    #ax[1].plot(lats, np.max(ON2, axis=0),'--%s' % color[mind])
    #ax[1].plot(lats, np.min(ON2, axis=0),'--%s' % color[mind])
    ax[1].set_xlabel('Lat. (Degrees N)')
    ax[1].set_ylabel('Mean N2 from %i-400 km' % amin)
"""


for ind in [0,]:  # 1]:
    ax[ind].legend()
    ax[ind].grid()
    ax[ind].set_xlim(-90, 90)
plt.show()


"""
    for specind, spec in enumerate(species):
        im = ax[ct, specind].contourf(times, lats, np.mean(nden[spec], 2).T)
        im.set_clim(1E17, 5E17)
    ct += 1
cb_ax = fig.add_axes([0.83, 0.1, 0.02, 0.8])
cbar = fig.colorbar(im, cax=cb_ax)
plt.show()






#    plt.plot(times_l[0], np.mean(O_N2_ratio, (1, 2)))
plt.grid()
plt.ylabel('O/N2 120-300 km')
plt.ylim([0, 1])
from matplotlib.dates import  DateFormatter
ax = plt.gca()
ax.xaxis.set_major_formatter( DateFormatter('%Y/%m/%d') )

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 15}
matplotlib.rc('font', **font)
plt.show()
"""
