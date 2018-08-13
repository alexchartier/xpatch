import pdb
import pickle
import numpy as np
import datetime as dt


def main(
        in_fname = '/media/alex/USB/2016_Ampere_data/%Y%m%dAmp_invert.ncdf',
        out_fname = '/media/alex/USB/ampere_pkl/ampere_%y%m%d.pkl',
        starttime=dt.datetime(2016, 1, 1),
        endtime=dt.datetime(2016, 12, 31),
        timestep=dt.timedelta(days=1),
    ):

    time = starttime
    while time <= endtime:
        ncdf4_to_pkl(time.strftime(in_fname), time.strftime(out_fname))    
        print('saved %s' % time.strftime(out_fname))
        time += timestep 
    
    with open(out_fname, 'rb') as f:
        data = pickle.load(f)

    rms_curr_dens = rmse(data['b_eci'])
    pdb.set_trace()
    plot_currents(data)
    


def rmse(var):
    return np.sqrt(np.mean(var ** 2))


def plot_currents(data):
    import matplotlib.pyplot as plt
    from matplotlib import cm
    plt.polar() 
    plt.tricontourf(data['mlt'][0, :], data['colat'][0, :], data['Jr'][0, :], cmap='bwr', vmin=-1, vmax=1)
    plt.colorbar() 
    plt.show()
    

def ncdf4_to_pkl(in_fname, out_fname):
    from netCDF4 import Dataset
    rootgrp = Dataset(in_fname)
    data = {}
    for k, v in rootgrp.variables.items():
        data[k] = v[...].filled().astype(float)
        data[k][data[k] == v[...].fill_value] = np.nan

    with open(out_fname, 'wb') as f:
        pickle.dump(data, f)


if __name__ == '__main__':
    main() 
