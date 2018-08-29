# Plot the median density-latitude distribution of Swarm perpendicular velocities

from spacepy import pycdf
import glob
import pdb
import numpy as np
import pickle
import datetime as dt
from scipy import stats
from proc_swarm_efi import efi_to_dataframe
import matplotlib.pyplot as plt

def main(fname_format=['./data/swarm_efi/SW_EXPT_EFI%s', '_TIICT_%Y%m%d*.cdf'],
                  sats=['A', 'B'], 
             starttime=dt.datetime(2016, 1, 1), 
               endtime=dt.datetime(2016, 12, 31), 
            latbins=np.linspace(-90, 90, 91),
            vars=['timestamp', 'latitude', 'longitude', 'radius', 'qdlat', 'viy', 'qy', 'mlt', 'bx', 'by', 'bz'], 
            out_fname='data/swarm_vels.pkl',
    ):  
    try:
        with open(out_fname, 'rb') as f:
            median_vels = pickle.load(f)
    except:
        print("Couldn't load %s" % out_fname)
        # Load in SWARM velocity files each day and calculate the MLAT median
        median_vels = {}
        for sat in sats:
            median_vels[sat] = {}
            median_vels[sat]['vels'] = np.zeros((len(latbins) - 1, (endtime - starttime).days + 1)) * np.nan
        times = []
        time = starttime
        while time <= endtime:
            print(time)
            times.append(time)
            for sat in sats:
                try:
                    fin = pycdf.CDF(glob.glob(fname_format[0] % sat + time.strftime(fname_format[1]))[0])
                    data = {v: fin[v][...] for v in vars}
                    efi_df = efi_to_dataframe(data, lat_cutoff=0)
                    bin_meds, bin_edges, binnumber = stats.binned_statistic(efi_df['lat_mag'], \
                          np.abs(efi_df['viy']), statistic='median', bins=latbins) 
                    median_vels[sat]['vels'][:, (time - starttime).days] = bin_meds

                except:
                    print(time.strftime('No file on %Y/%m/%d'))
            time += dt.timedelta(days=1)
        for sat in sats:
            median_vels[sat]['times'] = np.array(times)
        with open(out_fname, 'wb') as f:
            pickle.dump(median_vels, f)

    # Plot 
    bin_centres = (latbins[:-1] + latbins[1:]) / 2

    fig, axes = plt.subplots(1, 1, sharey=True)
    mean_vels = np.nanmean(np.stack((median_vels['A']['vels'], median_vels['B']['vels'])), axis=0)
    im = axes.pcolor(median_vels['A']['times'], bin_centres, mean_vels, vmin=0, vmax=1000)
    axes.axhline(y=60, linewidth=2, color='r')
    axes.axhline(y=-60, linewidth=2, color='r')
    axes.set_ylabel('MLAT') 
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    clb = fig.colorbar(im, cax=cbar_ax)
    clb.set_label(r'Swarm median horz. $\perp$ ion drift speed (m/s)') 
    plt.show()


if __name__ == '__main__':
    main()

