"""
calc_prod_loss.py
Calculate the production and loss rates based on MSIS
"""

import pdb 
import sys 
import numpy as np
import pickle
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates
sys.path.append('../fusionpp/src/nimo/')
from physics import production, loss_rate, EUVAC
from pyglow import pyglow
from dateutil.tz import tzoffset 


def main(calc=False, calc_crd=False, plot=True):
    # Set inputs
    pkl_fname = 'data/msis_iri/%Y%m%d_%H%M.pkl'
    out_fname = 'data/prod_loss/%Y%m%d_%H%M.pkl'

    times_l = []
    for start_t in [dt.datetime(2014, 1, 1), dt.datetime(2014, 7, 1)]:
        times_l.append(
            [start_t + dt.timedelta(hours=x) for x in range(0, 24 * 9, 3)],
        )
    times_l = np.array(times_l)

    crd_fname = './data/apex_crd.pkl'
    
    if calc_crd:  # requires python 3
        from apexpy import Apex
        crd = {}
        d_alt = 20 
        alts = np.arange(200, 403, d_alt)
        A = Apex(date=times_l.flatten()[0])
        mlats = np.concatenate([np.arange(-90, -69, 2), np.arange(70, 91, 2)])
        mlons = np.arange(0, 360, 10)
        lats, lons = [], []
        for mlat in mlats:
            for mlon in mlons:
                lat, lon = A.convert(mlat, mlon, 'apex', 'geo', height=300)
                lats.append(lat)
                lons.append(lon)
        lats = np.array(lats)
        lons = np.array(lons)
        alts.shape = (alts.shape[0], 1)
        lats.shape = (1, lats.shape[0])
        lons.shape = (1, lons.shape[0])
        crd['alt'] = np.tile(alts, (1, lats.shape[1]))
        crd['lat'] = np.tile(lats, (alts.shape[0], 1))
        crd['lon'] = np.tile(lons, (alts.shape[0], 1))
        with open(crd_fname, 'wb') as f:
            pickle.dump(crd, f, protocol=0)
    else:# requires python 2
        with open(crd_fname, 'rb') as f:
            crd = pickle.load(f)
        for k, v in crd.items():
            crd[k] = v.flatten()

    if calc:
        for time in times_l.flatten():
            msis = load_msis(time, crd['alt'], crd['lat'], crd['lon'], pkl_fname)
            calc_prod_loss(time, msis, out_fname)
    if plot:
        losses = []
        ne = []
        for time_l in times_l:
            out_f = {}
            for time in time_l:
                with open(time.strftime(out_fname), 'rb') as f:
                    out = pickle.load(f)
                for key, val in out.items():
                    val.shape = (1, val.shape[0])
                    if time == time_l[0]:
                        out_f[key] = val
                    else:
                        out_f[key] = np.concatenate([out_f[key], val])
            pdb.set_trace()
            losses.append(n2loss(out_f))
            # losses.append(out_f['loss_r'])
            ne.append(out_f['ne'])
       
        nalt = len(np.unique(crd['alt']))
        shape = (nalt, len(crd['alt']) / nalt)
        for k, v in crd.items():
            crd[k].shape = shape
        losses = np.array(losses)
        ne = np.array(ne)
        varshape = (losses.shape[0], losses.shape[1], shape[0], shape[1]) 
        losses.shape = varshape
        ne.shape = varshape
        plot_loss(times_l, crd['alt'], crd['lat'], losses, ne)


def plot_loss(times_l, alt, lat, loss_r, ne):
    colrs = {
        'jan': ['r-', 'g--'],
        'jul': ['b--', 'y-'], 
    }
    heminds = [lat[0, :] > 0, lat[0, :] < 0]
    time_0 = times_l[0]
    fig, ax = plt.subplots()
    alt_s = np.unique(np.diff(alt[:, 0])) * 1E5
    alt_min = 100
    altind = np.unique(alt[:, 0]) > alt_min
    assert len(alt_s) == 1, 'Assuming uniform height grid - must stop'
    loss_r *= alt_s  # multiply by height of each voxel ahead of integration
    ne *= alt_s

   # fig, ax = plt.subplots(2, 1)
   # im = ax[0].contourf(times_l[0], lat, (np.sum(np.mean(loss_r[0, :, :, :], 3), 1) / ne_int[0]).T, vmin=0, vmax=0.5)
   # fig.colorbar(im, ax=ax[0])
   # im2 = ax[1].contourf(times_l[1], lat, (np.sum(np.mean(loss_r[1, :, :, :], 3), 1) / ne_int[1]).T, vmin=0, vmax=0.5)
   # fig.colorbar(im2, ax=ax[1])
    
    for mon, cs in colrs.items():
        mon_ind = 0 if mon == 'jan' else 1
        for ind, hemind in enumerate(heminds):
            ne_H = ne[:, :, :, hemind]
            loss_r_H = loss_r[:, :, :, hemind]
            ne_int = np.sum(ne_H[mon_ind, :, altind, :], 0)
            loss_n = np.sum(loss_r_H[mon_ind, :, altind, :], 0) / ne_int
            hem = 'N' if min(lat[0, hemind]) > 0 else 'S'
            lifetimes = 1/np.median(loss_n, 1)
            lt_daily = np.median(np.reshape(lifetimes, (28, 8)), 1)
            times = time_0[::8]
            ax.plot(times, lt_daily, cs[ind], label='%s >70 %s MLAT' % (mon, hem))

    dateFmt = matplotlib.dates.DateFormatter('%d')
    ax.xaxis.set_major_formatter(dateFmt)
    fig.autofmt_xdate()
    ax.set_ylim([5, 140])
    ax.grid()
    ax.legend()
    ax.set_xlabel('Day of Month')
    ax.set_ylabel(r'%i-400 km Avg. Plasma Lifetime $(s^{-1})$' % alt_min)
    plt.show()
    

def calc_prod_loss(time, msis, out_fname):
    out = {}
    solar_flux = EUVAC(msis['f107'], msis['f107a'])
    sf = EUVAC(np.array(200), np.array(200))
    msis['neutrals']['temp'] = msis['temps']['Tn_msis']
    out['N2'] = msis['neutrals']['N2']
    out['O2'] = msis['neutrals']['O2']
    out['Ti'] = msis['temps']['Ti']
    out['loss_r'] = loss_rate(
        msis['neutrals']['N2'], msis['neutrals']['O2'], msis['temps']['Ti'], 
        msis['iono']['ne'],
    )
    out['ne'] = msis['iono']['ne']
#    out['prod'] = production(
#        time, msis['dims']['alt'] * 1E5, msis['dims']['lat'], msis['dims']['lon'],
#        msis['neutrals'], solar_flux,
#    )

    with open(time.strftime(out_fname), 'wb') as f:
        pickle.dump(out, f)
        print('Saved %s' % time.strftime(out_fname))


def load_msis(time, alts, lats, lons, pkl_fname, reload=True):
    # Save MSIS/IRI input data
    sf = 'f107', 'f107a'
    terms = {
        'neutrals': ['O', 'O2', 'N2'],
        'iono': ['ne'],
        'temps': ['Te', 'Ti', 'Tn_msis', 'Tn_iri'],
        'dims': ['alt', 'lat', 'lon'],
    }
    arr_shape = alts.shape
    if reload: 
        msis = {}
        for name, term in terms.items():
            msis[name] = {}
            for t in term:
                msis[name][t] = np.ones(arr_shape) * np.nan
    
        for ind, alt in enumerate(alts):
            lat = lats[ind]
            lon = lons[ind]
            pt = pyglow.Point(time, lat, lon, alt)
            pt.run_msis()
            pt.run_iri()
            for term in sf:
                msis[term] = pt.__dict__[term]
            for term in terms['neutrals']:
                msis['neutrals'][term][ind] = pt.nn[term]
            for term in terms['temps']:
                msis['temps'][term][ind] = pt.__dict__[term]
            msis['iono']['ne'][ind] = pt.ne
            msis['dims']['alt'][ind] = alt
            msis['dims']['lat'][ind] = lat
            msis['dims']['lon'][ind] = lon
        with open(time.strftime(pkl_fname), 'wb') as f:
            pickle.dump(msis, f)
            print('Saved %s' % time.strftime(pkl_fname))
    else:
        with open(time.strftime(pkl_fname), 'rb') as f:
            msis = pickle.load(f)
            print('Loaded %s' % time.strftime(pkl_fname))
    return msis


def sami_o_loss(
):
    # o+ + h --> h+ + o   (bb)

    chrate_Op['H'] = 2.5e-11                 \
               * sqrt( tn(iz,nfl,nll) )      

    # o+ + n2 --> no+ + n (bb)

    chrate_Op['N2'] = 1.533e-12 - 5.920e-13 * ti300o + 8.600e-14 * ti300o ** 2

    if ti > 1700:
      chrate_Op['N2'] = 2.730e-12 -                \
                      1.155e-12 * ti300o +       \
                      1.483e-13 * ti300o ** 2

    # o+ + o2 --> o2+ + o
    chrate_Op['O2'] = 2.820e-11 -                    \
               7.740e-12 * ti300o +              \
               1.073e-12 * ti300o ** 2 -         \
               5.170e-14 * ti300o ** 3 +         \
               9.650e-16 * ti300o ** 4

    # o+ + no --> no+ + o
    chrate_Op['NO'] = 1.0e-12   



if __name__ == '__main__':
    main()
