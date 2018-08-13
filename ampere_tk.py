"""
ampere_tk.py
Created by Alex Chartier on 8/13/2018
"""

import pdb
import pickle
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt


def main(
        in_fname = '/media/alex/USB/2016_Ampere_data/%Y%m%dAmp_invert.ncdf',
        out_fname = '/media/alex/USB/ampere_pkl/ampere_%y%m%d.pkl',
        starttime=dt.datetime(2016, 1, 1),
        endtime=dt.datetime(2016, 12, 31),
        timestep=dt.timedelta(days=1),
    ):

    time = starttime
    rms_B = {'north': [], 'south': []}
    while time <= endtime:
        try:
            with open(time.strftime(out_fname), 'rb') as f:
                data = pickle.load(f)
        except:
            data = ncdf4_to_pkl(time.strftime(in_fname), time.strftime(out_fname))    
            print('saved %s' % time.strftime(out_fname))

        for hem in rms_B.keys():
            rms_B[hem].append(rms(data['b_eci']))
        time += timestep 
    
    plt.plot(rms_B)
    plt.show() 


def rms(var):
    return np.sqrt(np.mean(var ** 2))


def plot_currents(data):
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
    return data


if __name__ == '__main__':
    main() 
