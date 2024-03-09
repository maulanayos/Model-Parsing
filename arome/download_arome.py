import xarray as xr
import numpy as np
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from dateutil.relativedelta import relativedelta
import bz2
import ftputil
import os
import time
import glob
from datetime import datetime, timedelta
import pygrib
import pandas as pd
import re

'''
Data Needed
	Gust
	CAPE v
	Temperature (hPa) n
	ISO (0o C) x
	Ceiling x
	Visibility x
	Cloudiness Convective x
	Cloud Cover v
	Cloudiness Low v
	Cloudiness Medium v
	Cloudiness High v
	WindTem n wind chill?
	TKE (hPa)

arrunknown =  ['CAPE_INS','cape']

arrisobaric = [['cc','tcc'], ]

arrsurface = ['hcc','['vis','vis'],]


Data faulty's - Assessed from 2024/02/19_00 data
0. !!CANT BE ACESSED WITH XARRAY EVEN THO FILTERED BY ONE SHORTNAME!!
0.1 PATTERN DIDNT INDICATE FORECAST TIME, PATTERN-0 IS NOT GRIB WITH FCST TIME 0
1. variative stepType and there are variables and type of level named "unknown"
2. tke missing first 4 levels
3. vis has different stepType than the rest surface data
4. iconsistent gribmessage, some grib have several gribmessages missing (making it more difficult to do machine learning)
'''


arrisobaric = [['cc','tcc'],
               ['t','t'],
               ['tke','tke']]

arrsurface = [['hcc','hcc'],
              ['lcc','lcc'],
              ['mcc','mcc'],
              ['vis','vis']]

arrother = [['10efg','efg10',None],
          ['10nfg','nfg10',None],
          ['CAPE_INS','CAPE_INS','cape']]

inpath = '/mnt/data/input/arome'

def gendts():
    import datetime
    dtn = datetime.datetime.utcnow() - timedelta(hours=12)
    if (dtn.hour >= 0) & (dtn.hour < 12):
        cycl = 0
    elif (dtn.hour >= 12) & (dtn.hour <= 24):
        cycl = 12
    
    dtc = datetime.datetime(dtn.year,dtn.month,dtn.day,cycl)
    return dtc
    


def fileNamegen(dtn):
    pattern = f'_arome_indo_{dtn:%Y}{dtn:%m}{dtn:%d}_{dtn:%H}'
    inpathgrib = Path(f'{inpath}/{dtn:%Y}/{dtn:%m}/{dtn:%d}_{dtn:%H}')
    with ftputil.FTPHost("172.19.0.47", "transmet2023", "99xrAQbW2qN7h") as ftp_host:
            names = ftp_host.listdir('MODELAROME/')
            filename = [x for x in names if re.search(pattern,x)]
    return filename, inpathgrib

def fileDownload(grbzname, dtn, maxretry : int=10,retryfreq : int = 60):    #nanti ganti 60 aja untuk debug 3 detik  
    import time
    with ftputil.FTPHost("172.19.0.47", "transmet2023", "99xrAQbW2qN7h") as ftp_host:
        file = f'MODELAROME/{grbzname}'
        basename = os.path.basename(file)
        inpathgrib = Path(f'{inpath}/{dtn:%Y}/{dtn:%m}/{dtn:%d}_{dtn:%H}')
        fout = Path(inpathgrib, basename)
        fout.parent.mkdir(parents=True, exist_ok=True)
    
        if fout.exists():
            return fout
        else:
            for i in range(maxretry):
                if ftp_host.path.exists(file):
                    ftp_host.download(file,fout)
                    print(f'file downloaded {fout}')
                    return fout
                else:
                    print(f'File {file} not found')
                time.sleep(retryfreq)
            print(f'File {file} not found after {maxretry} retries')
            pass


def initdfandisobar(filenames):
    dscc = xr.open_mfdataset(filenames[1:],engine="cfgrib",
                       backend_kwargs=dict(filter_by_keys={'typeOfLevel': 'isobaricInhPa','shortName':'cc'}),
                       concat_dim='valid_time', combine='nested', parallel=True)
    lats,lons,lev, time = dscc['latitude'].values, dscc['longitude'].values, dscc['isobaricInhPa'].values, dscc['valid_time'].values
    cc = dscc['cc'].values
    
    dst = xr.open_mfdataset(filenames[1:],engine="cfgrib",
                       backend_kwargs=dict(filter_by_keys={'typeOfLevel': 'isobaricInhPa','shortName':'t'}),
                       concat_dim='valid_time', combine='nested', parallel=True)

    t = dst['t'].values

    dumdum = np.empty((len(time),5,len(lats),len(lons)))
    dumdum[:,:,:,:] = np.nan
    dstke = xr.open_mfdataset(filenames[1:],engine="cfgrib",
                       backend_kwargs=dict(filter_by_keys={'typeOfLevel': 'isobaricInhPa','shortName':'tke'}),
                       concat_dim='valid_time', combine='nested', parallel=True)
    
    tke = dstke['tke'].values
    tke = np.concatenate((tke,dumdum), axis=1)
    
    ds = xr.Dataset(coords=dict(
            time=(["time"],time),
            level=(["level"], lev),
            lat=(["lat"], lats),
            lon=(["lon"], lons)))
    
    ds['cc']=(('time','level','lat','lon'),cc)
    ds['t'] = (('time','level','lat','lon'),t)
    ds['tke'] = (('tkeime','level','lat','lon'),tke)

    return ds

    

def surfaceret(var,filenames):
    if var[0] != 'vis' :    
        ds = xr.open_mfdataset(filenames[1:],engine="cfgrib",
                        backend_kwargs=dict(filter_by_keys={'typeOfLevel': 'surface','stepType': 'instant','shortName':var[0]}),
                        concat_dim='valid_time', combine='nested', parallel=True)
        return ds[var[1]].values
    else :
        ds = xr.open_mfdataset(filenames[1:],engine="cfgrib",
                       backend_kwargs=dict(filter_by_keys={'typeOfLevel': 'surface','stepType': 'min','shortName':var[0]}),
                       concat_dim='valid_time', combine='nested', parallel=True)
        return ds[var[1]].values
    
def otherret(vars,filenames):
    value = []
    for var in vars:
        ds = xr.open_mfdataset(filenames[1:],engine="cfgrib",
                       backend_kwargs=dict(filter_by_keys={'shortName':var[0]}),
                       concat_dim='valid_time', combine='nested', parallel=True)
        value.append(ds[var[1]].values)

    return value
def initNC():
    dtn = gendts()
    filenames, inpathgrib = fileNamegen(dtn)

    now = datetime.now()
    for filename in filenames:
        fileDownload(filename,dtn)
        
    print(f'done downloading, elapsed time:{(datetime.now() - now)}')

    filenames = glob.glob(f'{inpathgrib}/*.grib')
    filenames.sort()

    #Dataset initiation
    ds = initdfandisobar(filenames)

    for var in arrsurface:
        value = surfaceret(var,filenames)
        ds[var[0]] = (('time','lat','lon'),value)

    others = otherret(arrother,filenames)

    for var in enumerate(arrother):
        name = var[1][2]
        if name:
            ds[name] = (('time','lat','lon'),others[var[0]])
        else :
            ds[var[1][1]] = (('time','lat','lon'),others[var[0]])

    print(f'done processing nc, elapsed time:{(datetime.now() - now)}')


    now = datetime.now()
    #EXPORT
    chunk3d = [3,6,12]
    chunk4d = [3,3,6,12]
    encoding = {}
    encoding_keys = ("_FillValue", "dtype", "scale_factor", "add_offset", "grid_mapping", "chunksizes")
    for data_var in ds.data_vars:
        if len(ds[data_var].shape) == 3:
            chunks = chunk3d
        elif len(ds[data_var].shape) == 4:
            chunks = chunk4d

        encoding[data_var] = {key: value for key, value in ds[data_var].encoding.items() if key in encoding_keys}
        encoding[data_var].update(zlib=True, complevel=5, chunksizes = chunks)
    fno = f'/mnt/data/output/arome/arome_{dtn:%Y%m%d}_{dtn:%H%M}.nc'
    fno = Path(fno)
    # save to netcdf
    fno.unlink(missing_ok=True)
    fno.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(
        fno,
        'w',
        format = "NETCDF4",
        encoding = encoding
    )
    print(f'done saving, elapsed time:{(datetime.now() - now)}')


if __name__ == '__main__':

    initNC()