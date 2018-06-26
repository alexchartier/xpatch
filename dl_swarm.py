#!/usr/env/bin/python3
""" 
dl_swarm.py
Script to download and preprocess the SWARM GPS TEC data/LP data
Swarm provides a Level-2 TEC product, handily preprocessed and ready for use. 
The files come as .ZIP, which contains a .DBL that is actually a CDF file
Details here:
http://swarm-wiki.spacecenter.dk/mediawiki-1.21.1/index.php/Level_2_Product_Definitions#Intermediate_and_final_products
Swarm files:
    ftp://swarm-diss.eo.esa.int/Level2daily/Current/TEC/TMS/

unzip files: 
"""
import datetime, pdb, sys, os
sys.path.append('/Users/chartat1/fusionpp/src/data_prep/gps/')
from get_gps_data import dl_data

instrument = 'lp'
dl_times = [datetime.datetime(2017, 4, 1), datetime.datetime(2018, 1, 1)]

if instrument == 'gps':
    servername = 'data/swarm_server_names.txt'
    # dl_dir = './data/swarm_tec/'
    dl_dir = '/Volumes/Seagate/data/swarm/gps_tec/'

elif instrument == 'lp':
    servername = 'data/swarm_lp.txt'
    # dl_dir = '/Volumes/Seagate/data/swarm/lp/'
    dl_dir = './data/swarm_lp/'
    # servername = 'data/swarm_lp_adv.txt'
    # dl_dir = '/Volumes/Seagate/data/swarm/lp_adv/'

elif instrument == 'efi':
    servername = 'data/swarm_efi.txt'
    dl_dir = '/Volumes/Seagate/data/swarm_efi/'

dl_days = dl_data(dl_times[0], dl_times[1], dl_dir, servername, datatype='swarm', dirnames='not_smart')



