"""
dmsp_tk.py
Script to perform data analysis of DMSP in situ electron densities

Tailored data-download command for madrigal Python downloader:
globalIsprint.py --verbose --url=http://cedar.openmadrigal.org --parms=YEAR,MONTH,DAY,HOUR,MIN,SEC,MLT,GDALT,GDLAT,GLON,AACGM_LAT,MLAT,AACGM_LONG,MLONG,NE,VERT_ION_V,HOR_ION_V --output=/tmp --user_fullname="Alex+T+Chartier" --user_email=alex.chartier@outlook.com --user_affiliation="APL" --startDate="05/01/2014" --endDate="08/01/2018" --inst=8100 --format=Hdf5 

Bulk download command:
globalDownload.py --verbose --url=http://cedar.openmadrigal.org --outputDir=/home/alex/Downloads/ --user_fullname="Alex+T+Chartier" --user_email=alex.chartier@outlook.com --user_affiliation="APL" --format="hdf5" --startDate="01/01/2014" --endDate="02/01/2014" --inst=8100


"""

import h5py
import numpy as np
import pdb
import datetime as dt
import pickle

sat = '17'
fname = '/home/alex/Downloads/dms_%Y%m%d' +  '_%ss1.001.hdf5' % sat
out_fname = './data/dmsp/dmsp_%s_' % sat + '%Y%m%d.pkl'
starttime = dt.datetime(2014, 1, 1)
endtime = dt.datetime(2014, 2, 1)
timestep = dt.timedelta(days=1)

time = starttime
while time < endtime:
    hf = h5py.File(time.strftime(fname), 'r')
    vals = {}
    headerdata = hf['Metadata']['Data Parameters'][...]
    header = [h[0].decode('UTF-8') for h in headerdata]
    assert 'NE' in header, 'NE not in header'
    data = hf['Data']['Table Layout'][...]  # 1Hz DMSP data
    for d in data:
        for ind, h in enumerate(header):
            vals[h] = d[ind]
        vals['time'] = dt.datetime(vals['YEAR'], vals['MONTH'], vals['DAY'], vals['HOUR'], vals['MIN'], vals['SEC'])

    with open(time.strftime(out_fname), 'wb') as f:
        pickle.dump(vals, f)
    time += timestep
        
        
 
 
