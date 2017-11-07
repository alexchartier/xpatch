# F10.7 plotter

import pdb
from pyglow import pyglow
import datetime as dt
import pickle
import matplotlib.pyplot as plt
import matplotlib

time = dt.datetime(2014, 8, 1)
endtime = dt.datetime(2017, 7, 1)
times = []
f107 = []
while time <= endtime:
    pt = pyglow.Point(time, 100, 0, 0)
    f107.append(pt.f107a)
    times.append(time)
    time += dt.timedelta(days=1)

plt.plot_date(times, f107, '-')
plt.grid()
plt.ylabel('81-day averaged F10.7')
plt.ylim([0, 180])
pdb.set_trace()
from matplotlib.dates import  DateFormatter
ax = plt.gca()
ax.xaxis.set_major_formatter( DateFormatter('%Y/%m/%d') )

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 15}
matplotlib.rc('font', **font)
plt.show()
