#!/Library/Frameworks/Python.framework/Versions/3.6/bin/python3
"""
plot_patch_ct.py
Script to plot the output of the SWARM patch counter (either TEC or LP)
"""
import pdb
import matplotlib.pyplot as plt
import matplotlib
import pickle
import collections
import datetime as dt
import numpy as np
import socket
import sys
sys.path.insert(0, '/users/chartat1/fusionpp/fusion')
import physics
import count_passes

starttime = dt.datetime(2014, 8, 1)
endtime = dt.datetime(2017, 6, 1)
satellites = 'A', 'B'

lat_cutoff = 55
plot_lat_cutoff = 70
approach = 'alex'
instrument = 'Langmuir Probe'  # or 'GPS'

# Langmuir Probe
if instrument == 'Langmuir Probe':
    fin = '/Volumes/Seagate/data/swarm/proc_lp/%s/%i_deg/' % (approach, lat_cutoff) \
          + 'lp_%Y%m%d_' + '%ideg.pkl' % lat_cutoff    
    norm_fin = '/Volumes/Seagate/data/swarm/pass_ct/pass_norm_%Y%m%d' + '_%ideg.pkl' % lat_cutoff
    colour = 'g'
    freq = 2  # hz

elif instrument == 'GPS':
    elev_cutoff = 25
    colour = 'y'
    freq = 1  # hz
    if socket.gethostname() == 'chartat1-ml2':
        # Work GPS
        fin = '/Volumes/Seagate/data/swarm/proc_gps/patch_ct_%Y%m%d.pkl'
        norm_fin = '/Volumes/Seagate/data/swarm/pass_ct/pass_norm_%Y%m%d_' + '%ideg.pkl' % lat_cutoff

    elif socket.gethostname() == 'chartat1-ml1':
        # Home GPS
        fin = './gps_proc/patch_ct_%Y%m%d.pkl'
        norm_fin = '/Volumes/Seagate/data/swarm/pass_ct/pass_norm_%Y%m%d.pkl'


def main():
    patch_ct = get_patch_ct(starttime, endtime, satellites, fin)
    if instrument == 'GPS':
        for sat in satellites:
            elev_ind = np.array(patch_ct[sat]['elev']) > elev_cutoff
            for key, val in patch_ct[sat].items():
                if key not in ('t1', 't2', 'tec_b1', 'tec_b2', 'params'):
                    if len(np.array(val).shape) == 1:
                        val = np.expand_dims(np.array(val), 1)  # expand 1D variables
                    patch_ct[sat][key] = np.array(val)[elev_ind]
    else:
        for sat in satellites:
            lat_ind = np.array(np.abs(patch_ct[sat]['lat_mag']) > plot_lat_cutoff)
            for key, val in patch_ct[sat].items():
                print(key)
                pdb.set_trace()
                patch_ct[sat][key] = np.array(val)[lat_ind]
            pdb.set_trace()
        
            
            
    # norm_ct = count_passes.get_norm_ct(norm_fin, starttime=starttime, endtime=endtime, sats=satellites)
    # plot_t_doy(patch_ct, norm_ct, vartype='lt')
    plot_hist(patch_ct)
    # plot_magnitudes(patch_ct)  # Determine the relative magnitude of all the patches counted in each hemisphere
    # oneday_timeseries()
    # plot_ut(patch_ct, norm_ct)
    # plot_mlt(patch_ct, norm_ct)
    # plot_polar(patch_ct, crd='mag')


def plot_magnitudes(patch_ct):
    nbins = 50
    bins = np.linspace(2, 15, nbins)
    ylimit = [0, 300]
    xlabel = 'Relative magnitude'
    hems = 'north', 'south'
    ct = 0
    for sat in satellites:
        if instrument == 'Langmuir Probe':
            mag = np.array(patch_ct[sat]['ne_rm']).flatten() / np.array(patch_ct[sat]['ne_bg']).flatten()
        else:
            mag = np.array(patch_ct[sat]['tec']).flatten() / np.array(patch_ct[sat]['tec_bg']).flatten()
        sat_lats = np.array([x[0] for x in patch_ct[sat]['lat_geo']])
        nh_ind = sat_lats > 0
        sh_ind = sat_lats < 0
        for hem in hems:
            ct += 1
            plt.subplot(len(satellites), 2, ct)
            mag_h = mag[nh_ind] if hem == 'north' else mag[sh_ind]
            print('Sat %s, %s hemisphere, median: %2.2g, mean: %2.2g, max: %2.2g' % (sat, hem, np.median(mag_h), np.mean(mag_h), np.max(mag_h)))
            n, bins, patches = plt.hist(mag_h, bins=bins)
            plt.ylim(ylimit)
            if np.ceil(ct / 2) == len(satellites):
                plt.xlabel(r'Patch magnitude')
            else:
                plt.tick_params(
                                axis='x',          # changes apply to the x-axis
                                which='both',      # both major and minor ticks are affected
                                bottom='off',      # ticks along the bottom edge are off
                                top='off',         # ticks along the top edge are off
                                labelbottom='off') # labels along the bottom edge are off
            if np.mod(ct, 2) != 0:
                plt.ylabel('Patch count')
            
            plt.grid()
            plt.title('Satellite %s %s hemisphere' % (sat, hem))

    plt.suptitle(instrument, fontweight='bold')
    plt.show() 


def oneday_timeseries():
    import plot_timeseries
    plot_timeseries.main(instrument='Langmuir Probe')


def plot_t_doy(patch_ct, norm_ct, sats=['A', 'B'], vartype='mlt'):
    ct = 0  
    hems = 'north', 'south'
    hemdict = {'north': 'nh', 'south': 'sh'}
    for hem in hems:
        for sat in sats:
            sat_lats = np.squeeze(np.array(patch_ct[sat]['lat_mag']))
            hem_ind = sat_lats > 0 if hem == 'north' else sat_lats < 0
            month = np.array([t[0].month for t in patch_ct[sat]['times']])
            month_h = month[hem_ind]
            ut = np.array([t[0].hour + t[0].minute / 60 for t in patch_ct[sat]['times']])
             
            if vartype == 'ut':
                t_h = ut[hem_ind]
            elif vartype == 'mlt':
                mlon = np.squeeze(np.array(patch_ct[sat]['lon_mag']))
                mlon[mlon < 0] += 360
                mlt = calc_mlt(ut, mlon)
                t_h = mlt[hem_ind]
            elif vartype == 'lt':
                lon = np.squeeze(np.array(patch_ct[sat]['lon_geo']))
                lon[lon < 0] += 360
                lt = calc_mlt(ut, lon)
                t_h = lt[hem_ind]

            tbins = np.arange(0, 24.1)
            monthbins = np.arange(1, 13.1)
            H, xedges, yedges = np.histogram2d(month_h, t_h, bins=(monthbins, tbins))
            X, Y = np.meshgrid(monthbins, tbins)
            norm_h = norm_ct[sat]['%s_%s_2d' % (hemdict[hem], vartype)]
            norm_h[norm_h < (3600 * freq)] = np.nan  # mask out if less than an hour's data available
            H /= norm_h 
            H *= 3600 * freq
            H_masked = np.ma.masked_invalid(H).T
            
            ct += 1
            plt.subplot(len(hems), len(sats), ct)
            my_cmap = matplotlib.cm.get_cmap('viridis')
            my_cmap.set_under('w')
            plt.pcolormesh(X, Y, H_masked, vmin=0, vmax=15, cmap=my_cmap)
            plt.title('%s hemisphere, satellite %s' % (hem, sat))
            plt.xlim((monthbins[0], monthbins[-1]))
            plt.ylim((tbins[0], tbins[-1]))
            frame = plt.gca()
            if ct > 2:
                plt.xlabel('month')
            else:
                frame.axes.xaxis.set_ticklabels([])
            if np.mod(ct, 2) != 0:
                plt.ylabel(vartype.upper())
            else:
                frame.axes.yaxis.set_ticklabels([])
                plt.colorbar(label='normalized patch count per hour')
    plt.suptitle('%s %s to %s' % (instrument, starttime.strftime('%Y/%m/%d'), endtime.strftime('%Y/%m/%d')))

    plt.show()


def plot_ut(patch_ct, norm_ct, sats=['A', 'B']):
    ct = 0  
    hems = 'north', 'south'
    for hem in hems:
        ut_hist = np.zeros(norm_ct[sats[0]]['nh_mlt'][0].shape)
        for sat in sats:
            time = {}
            ut = np.array([t[0].hour + t[0].minute / 60 for t in patch_ct[sat]['times']])
            sat_lats = np.squeeze(np.array(patch_ct[sat]['lat_mag']))
            nh_ind = sat_lats > 0
            sh_ind = sat_lats < 0
            ut_h = ut[nh_ind] if hem == 'north' else ut[sh_ind]
            norm = norm_ct[sat]['nh_ut'] if hem == 'north' else norm_ct[sat]['sh_ut']
            binedges = norm[1]
            ut_hist_sat = np.histogram(ut_h, binedges)
            ut_hist += ut_hist_sat[0] / (norm[0] / 7200)  # Patches observed per satellite hour spent in each MLT bin

        ct += 1
        plt.subplot(1, 2, ct)
        plt.xlim(0, 24)
        if instrument == 'GPS':
            plt.ylim(0, 800)
        else:
            plt.ylim(0, 12)
       
        plt.bar(binedges[:-1], ut_hist, width=np.diff(binedges), edgecolor=colour, linewidth=0)
        plt.xlabel('UT Hour')
        if np.mod(ct, 2) != 0:
            plt.ylabel('Patch count / hour')
        plt.grid()
        plt.title('%s hemisphere' % (hem))

    plt.suptitle(instrument, fontweight='bold')
    plt.show()
                

def plot_mlt(patch_ct, norm_ct, sats=['A', 'B']):
    ct = 0  
    hems = 'north', 'south'
    for hem in hems:
        mlt_hist = np.zeros(norm_ct[sats[0]]['nh_mlt'][0].shape)
        for sat in sats:
            mlon = np.squeeze(np.array(patch_ct[sat]['lon_mag']))
            mlon[mlon < 0] += 360
            time = {}
            ut = np.array([t[0].hour + t[0].minute / 60 for t in patch_ct[sat]['times']])
            mlt = calc_mlt(ut, mlon)

            sat_lats = np.squeeze(np.array(patch_ct[sat]['lat_mag']))
            nh_ind = sat_lats > 0
            sh_ind = sat_lats < 0
            if (hem == 'south') and (sat == 'A'):  # Store a couple of things for the 15-16 MLT plots coming up
                times = np.squeeze(np.array(patch_ct[sat]['times']))
                glon = np.squeeze(np.array(patch_ct[sat]['lon_geo']))
                sh_15ind = np.logical_and(np.logical_and(mlt > 15, mlt < 16), sh_ind)
                times_15mlt = times[sh_15ind]
                glon_15mlt = glon[sh_15ind]

            mlt_h = mlt[nh_ind] if hem == 'north' else mlt[sh_ind]
            norm = norm_ct[sat]['nh_mlt'] if hem == 'north' else norm_ct[sat]['sh_mlt']
            binedges = norm[1]
            mlt_hist_sat = np.histogram(mlt_h, binedges)
            mlt_hist += mlt_hist_sat[0] / (norm[0] / 7200)  # Patches observed per satellite hour spent in each MLT bin

        ct += 1
        plt.subplot(1, 2, ct)
        plt.xlim(0, 24)
        if instrument == 'GPS':
            plt.ylim(0, 800)
        else:
            plt.ylim(0, 12)
       
        plt.bar(binedges[:-1], mlt_hist, width=np.diff(binedges))
        plt.xlabel('MLT Hour')
        if np.mod(ct, 2) != 0:
            plt.ylabel('Patch count / hour')
        plt.grid()
        plt.title('%s hemisphere' % (hem))

    plt.suptitle(instrument, fontweight='bold')
    plt.show()

    # 15 MLT spike plot 
    plt.subplot(1, 2, 1)
    plt.hist([t.hour for t in times_15mlt], np.arange(0, 24.1), color='g'); 
    plt.xlim(0, 24)
    plt.ylim(0, 36)
    plt.xlabel('UT hour')
    plt.ylabel('Patches detected')
    plt.title('(a)')
    plt.grid()
    
    plt.subplot(1, 2, 2)
    plt.hist(glon_15mlt, np.arange(-180, 180.1, 15), color='y')
    plt.xlim(-180, 180)
    plt.ylim(0, 36)
    plt.xlabel('geo. lon. (deg.)')
    plt.title('(b)')
    plt.grid()

    plt.show()               

def plot_hist(patch_ct, timestep=dt.timedelta(days=5)):
    ct = 0
    hems = 'north', 'south'
    for hem in hems:
        for sat in satellites:
            ct += 1
            ax = plt.subplot(len(satellites), 2, ct)
            times = np.squeeze(patch_ct[sat]['times'])
            day_ct = np.array([(time - starttime).days for time in times])
            nbins = round((endtime - starttime + dt.timedelta(days=1)) / timestep)
            sat_lats = np.squeeze(patch_ct[sat]['lat_geo'])
            nh_ind = sat_lats > 0
            sh_ind = sat_lats < 0
            day_ct_h = day_ct[nh_ind] if hem == 'north' else day_ct[sh_ind]
            day = dt.datetime(starttime.year, starttime.month, starttime.day)
            days = []            
            while day < endtime:
                days.append(day)
                day += timestep
        
            cts, bins = np.histogram(day_ct_h, len(days))
            ax.bar(days, cts, width=5, color=colour, edgecolor=colour, linewidth=0)
            ax.xaxis_date()
            ax.set_xticks(ax.get_xticks()[::2])

            plt.title('Satellite %s, %s hemisphere' % (sat, hem))
            if np.mod(ct, 2) != 0:
                plt.ylabel('Patch count / 5 days')
            frame = plt.gca()
            ymax = 80
            day_ctmin = min(day_ct)
            day_ctmax = max(day_ct)
            plt.ylim(0, ymax)
            # plt.xlim(0, day_ctmax)
  
            yr = starttime.year 
            dec_sols = []
            jun_sols = []
            while yr < endtime.year: 
                dec_sols.append(dt.datetime(yr, 12, 21))
                jun_sols.append(dt.datetime(yr, 6, 21))
                yr += 1

            cnt = 1
            for d in jun_sols:
                if cnt == 1:
                    plt.plot_date([d, d], [0, ymax], 'r--', label='June Solstice')
                    cnt += 1
                else:
                    plt.plot([d, d], [0, day_ctmax], 'r--')
             
            for d in dec_sols:
                if cnt == 2:
                    plt.plot([d, d], [0, day_ctmax], 'b--', label='December Solstice')
                    cnt += 1
                else:
                    plt.plot([d, d], [0, day_ctmax], 'b--')

            if ct == 2:
                plt.legend()
            if np.mod(ct, 2) == 0: 
                frame.axes.yaxis.set_ticklabels([])
            if ct < 3:
                frame.axes.xaxis.set_ticklabels([])
            plt.grid()
    plt.suptitle(instrument, fontweight='bold')
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 15} 
    matplotlib.rc('font', **font) 
    plt.show()


def plot_polar(patch_ct, crd='mag'):
    passes = sum_passes('./data/pass_ct/pass_%Y%m%d.pkl', crd=crd)
    sats = [s for s in patch_ct.keys()]
    sats.sort()

    hems = 'north', 'south'
    if crd == 'mag':
        latlim = 70
    else:
        latlim = 60

    ct = 0  
    for sat in sats:
        if crd == 'mag':
            lats = np.deg2rad(np.squeeze(np.array(patch_ct[sat]['lat_mag'])))
            lons = np.deg2rad(np.squeeze(np.array(patch_ct[sat]['lon_mag'])))
            unused_alts, pole_lat, pole_lon = physics.transform([1, 1], np.deg2rad([90, -90]), [0, 0], from_=['GEO', 'sph'], to=['MAG', 'sph'])
        elif crd == 'geo':
            lats = np.deg2rad(np.squeeze(np.array(patch_ct[sat]['lat_geo'])))
            lons = np.deg2rad(np.squeeze(np.array(patch_ct[sat]['lon_geo'])))
            unused_alts, pole_lat, pole_lon = physics.transform([1, 1], np.deg2rad([90, -90]), [0, 0], from_=['MAG', 'sph'], to=['GEO', 'sph'])
        pole_lat[1] = - pole_lat[1]
        lons[lons < 0] += 2 * np.pi
       
        latbins, lonbins =  passes[sat][crd][1:]
        counts = np.histogram2d(lats, lons, np.array((latbins, lonbins)))
        vals = counts[0] / (passes[sat][crd][0] / 7200)
        vals[vals == 0] = np.nan
        latvec, lonvec = np.meshgrid((latbins[:-1] + latbins[1:]) / 2, (lonbins[:-1] + lonbins[1:]) / 2, indexing='ij')
        
        for hem in hems: 
            ct += 1
            ax = plt.subplot(len(sats), len(hems), ct, polar=True)
            hemlat = latvec
            if hem == 'south':
                hemlat = -hemlat
            plt.ylim(0, np.deg2rad(90 - latlim))
            # vals[np.isnan(vals)] = 0
            sc = plt.pcolor(lonvec, np.pi / 2 - hemlat, vals, vmin=np.nanmin(vals), vmax=np.nanmax(vals))
            sc.cmap.set_under('white')
            hemind = 1 if hem == 'south' else 0
            plt.plot(pole_lon[hemind], np.pi / 2 - pole_lat[hemind], '.m', markersize=15)
            labels = ['%2.0f' % (90 - val) for val in np.linspace(0, 90 - latlim, 7)]
            labels = labels[1:]
            ax.set_yticklabels(labels)
            sc.cmap.set_under('white')
            """
            from mpl_toolkits.basemap import Basemap
            m = Basemap(projection='npstere',boundinglat=70,lon_0=270,resolution='l')
            # draw parallels and meridians.
            #m.drawparallels(np.arange(-70.,81.,20.))
            #m.drawmeridians(np.arange(-180.,181.,20.))
            sc = m.pcolor(lonvec, hemlat, vals)
            """
            plt.clim(0, 15)
            plt.colorbar(sc)
            if ct < 3:
                plt.title('Sat %s: %s hemisphere' % (sat, hem))
            else:
                plt.title('Sat %s                                 ' % sat)

    plt.show()


def sum_passes(fname_format, crd='mag'):
    time = starttime
    pass_ct_full = {}
    while time <= endtime:
        with open(time.strftime(fname_format), 'rb') as f:
            pass_ct = pickle.load(f)
        if time == starttime:
            for s in satellites:
                pass_ct_full[s] = {}
                pass_ct_full[s][crd] = np.array(pass_ct[s][crd])
        else:
            for s in satellites:
                try:
                    pass_ct_full[s][crd][0] += pass_ct[s][crd][0]
                except:
                    print('No counts on satellite %s on %s' % (s, time.strftime('%Y %m %d')))
        time += dt.timedelta(days=1)

    # for sat in satellites:
    #     pass_ct_full[sat][crd][pass_ct_full[sat][crd] == 0] = np.nan
    return pass_ct_full


def get_patch_ct(starttime, endtime, satellites, fin):
    patch_ct = {}
    time = starttime
    while time < endtime:
        try:
            with open(time.strftime(fin), 'rb') as f:
                count = pickle.load(f)
            if patch_ct == {}:
                patch_ct = count
            else:
                for sat in satellites:
                    try:
                        for key, val in count[sat].items():
                            if key != 'params':
                                patch_ct[sat][key] = patch_ct[sat][key] + val
                    except:
                        print('%s Missing file on satellite %s' % (time.strftime('%Y-%m-%d'), sat))
        except:
            print('No file on %s' % (time.strftime('%Y-%m-%d')))

        time += dt.timedelta(days=1)
    return patch_ct

def calc_mlt(ut, mlon):
    mlt = ut + mlon * 24 / 360
    mlt[mlt > 24] -= 24
    mlt[mlt < 0] += 24
    return mlt

if __name__ == '__main__':
    main()
