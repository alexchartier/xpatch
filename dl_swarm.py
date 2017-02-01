#!/Library/Frameworks/Python.framework/Versions/3.4/bin/python3
""" 
dl_swarm.py
Script to download and preprocess the SWARM GPS TEC data
Swarm provides a Level-2 TEC product, handily preprocessed and ready for use. 
The files come as .ZIP, which contains a .DBL that is actually a CDF file
Details here:
http://swarm-wiki.spacecenter.dk/mediawiki-1.21.1/index.php/Level_2_Product_Definitions#Intermediate_and_final_products
Swarm files:
    ftp://swarm-diss.eo.esa.int/Level2daily/Current/TEC/TMS/

unzip files: 
"""
import sys
sys.path.append('~/fusionpp/glimpse/')
import downloader, datetime, pdb

servername = 'swarm_lp.txt'
dl_times = [datetime.datetime(2015, 12, 20), datetime.datetime(2015, 12, 20)]
# dl_dir = '/Volumes/Seagate/data/swarm/lp/'
dl_dir = '~/xpatch/swarm_lp/'
dl_days = downloader.dl_data(dl_times, dl_dir, servername, datatype='swarm', dirnames='not_smart')






#def proc_gps_upward(infile, navfile, ofile, proc_dir):
#    """ 
#    Process upward-looking GPS RINEX data 
#    Function uses GPStk executables
#    NOTE: GPStk v1.* does not work on RINEX 3 format, so this does not work for SWARM RINEX files
#          However, Swarm has a level 2 processed TEC product that can be used instead
#    """
#    os.system(procdir + '/DiscFix --inputfile ' + infile + ' --RinexFile D.obs --dt 10 --smoothPR --smoothPH -f' + procdir + '/DiscFixTs.inp')
#    os.system(procdir + '/PRSolve --obs D.obs --nav ' + navfile + ' --outRinex prsolve.obs')
#    os.system(procdir + '/ResCor -IFprsolve.obs -OFR.obs -f' + procdir + '/ResCorTs.inp --nav ' + navfile)
#    os.system(procdir + '/gobias ' + infile)
#    os.system(procdir + '/posInterp --obs prsolve.obs  --mult 10 --outRinex PosInt10.obs')
#    os.system(procdir + '/ResCor -IFR.obs -OF' + ofile + ' -fsvbias.inp --RxRinex prsolve.obs')
#    os.system(procdir + '/Occult2Ncdf ' +  ofile)
#"""
##!/bin/tcsh
#set procdir=${IDA_INPUT_DIR}/bin/
#set infile=$argv[1]
#set ofile=$argv[2]
#set navfile=$argv[3]
#   
#
#"""
