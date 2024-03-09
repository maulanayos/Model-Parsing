"""
On development
@mamyos
"""

import os
# import pygrib
import datetime
import numpy as np
# import numexpr as ne
# import netCDF4
from time import sleep
from pathlib import Path
import urllib.request
# from tools.db_insert import mongoinsert
# import sentry_sdk
# sentry_sdk.init("https://0c2dbdaf918b4494b0f905b3a86d8eab@sentry.circlegeo.com/16")

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
maxtry = 100 # Total 8.3 Hours
waitsec = 300 # 5 minutes
tmpdiryos = Path('/home/opn/script/siam-netcdf-generator/coba/tes/envyosa/')
""" DOWNLOAD DATA """
# https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?file=gfs.t00z.pgrb2.0p25.anl&lev_1000_mb=on&lev_100_mb=on&lev_10_m_above_ground=on&lev_150_mb=on&lev_200_mb=on&lev_250_mb=on&lev_2_m_above_ground=on&lev_300_mb=on&lev_400_mb=on&lev_500_mb=on&lev_600_mb=on&lev_700_mb=on&lev_800_mb=on&lev_850_mb=on&lev_900_mb=on&lev_925_mb=on&lev_950_mb=on&lev_entire_atmosphere=on&lev_mean_sea_level=on&lev_surface=on&var_APCP=on&var_CAPE=on&var_GUST=on&var_HCDC=on&var_HGT=on&var_LCDC=on&var_LFTX=on&var_MCDC=on&var_PRATE=on&var_RH=on&var_TCDC=on&var_TMP=on&var_VIS=on&leftlon=95&rightlon=141&toplat=6&bottomlat=-11&dir=%2Fgfs.20231117%2F00%2Fatmos

# / format link gfsurl+[initialtime]+midsufgfsfile+gfsopt+(path)
#(path) = 'dir=%2Fgfs.+[tahun]+[bulan]+[tanggal]+%2F+[initialtime]+%2Fatmos
    
def dldURLCreator(iTime:str,tabutang:str):
    
    gfsurl = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?file=gfs.t'
    gfsopt = '&lev_1000_mb=on&lev_100_mb=on&lev_10_m_above_ground=on&lev_150_mb=on&lev_200_mb=on&lev_250_mb=on&lev_2_m_above_ground=on&lev_300_mb=on&lev_400_mb=on&lev_500_mb=on&lev_600_mb=on&lev_700_mb=on&lev_800_mb=on&lev_850_mb=on&lev_900_mb=on&lev_925_mb=on&lev_950_mb=on&lev_entire_atmosphere=on&lev_mean_sea_level=on&lev_surface=on&var_APCP=on&var_CAPE=on&var_GUST=on&var_HCDC=on&var_HGT=on&var_UGRD=on&var_VGRD=on&var_LCDC=on&var_LFTX=on&var_MCDC=on&var_PRATE=on&var_RH=on&var_TCDC=on&var_TMP=on&var_VIS=on&subregion=&leftlon=95&rightlon=141&toplat=6&bottomlat=-11&'
    midsufgfsfile = 'z.pgrb2.0p25.'
    extension = [] 
    for i in range(1, 240):
        tmpext = 'f'+str(i).zfill(3)
        extension.append(tmpext)
        urls = []
    for i in range(len(extension)):
        path = 'dir=%2Fgfs.' + (tabutang)+ '%2F'+(iTime)+'%2Fatmos'
        url = gfsurl + (iTime)+midsufgfsfile+extension[i]+gfsopt+path
        urls.append(url)
    return(urls)

def gfsdown(urls):
    for i in range(len(urls)):
        gfsfile = urls[i][653:661]+urls[i][65:68]+urls[i][80:85]
        gfslocalfile = tmpdiryos / 'gfsdowntemp' / gfsfile
        gfslocalfile.parent.mkdir(parents=True, exist_ok=True)
        
        print(gfslocalfile)
        print('temp file : %s' % gfslocalfile)
        for i in range(maxtry):
            try:
                response = urllib.request.urlopen(urls[i], timeout = 5)
                data = response.read()
                with open(gfslocalfile, 'wb') as f:  
                    f.write(data)
                break
            except:
                print(f'Download Failed, Attemp : {i}')
                sleep(waitsec)
                continue
        else:
            raise StopIteration('Maximum retry has been reached') 
    return(gfslocalfile)