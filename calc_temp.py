import pdb 
from pyglow import pyglow
import datetime as dt
import pickle
import matplotlib.pyplot as plt 
import matplotlib
import numpy as np


"""
calc Tn at the magnetic poles
"""

poles = [[-74.7648,124.6146], [83.4092, -84.7716]]
alt = 600
months = ['jan', 'jul']
starttime = {
    'jan': dt.datetime(2014, 1, 1),
    'jul': dt.datetime(2014, 7, 1),
}
endtime = {
    'jan': dt.datetime(2014, 2, 1), 
    'jul': dt.datetime(2014, 8, 1), 
}
step = dt.timedelta(hours=3)
times = {'jan': [], 'jul': []}
Tn = {}
t = {}
for month in months:
    t[month] = starttime[month]
    while t[month] < endtime[month]:
        times[month].append(t[month])
        t[month] += step 

    Tn[month] = np.zeros([len(times[month]), 2])
    for tind, time in enumerate(times[month]):
        for pind, pole in enumerate(poles):
            pt = pyglow.Point(time, pole[0], pole[1], alt)
            pt.run_msis()
            Tn[month][tind, pind] = pt.Tn_msis
    print('%s\n S: %1.1f  N: %1.1f' % (month, np.mean(Tn[month][:, 0]), np.mean(Tn[month][:, 1])))







