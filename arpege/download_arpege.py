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

'''
Data's Faulty :

1. fct 0 cant be accessed with xarray
2. gribs validDate is not working and it refers to analdate instead
3. data name patterns and upload date are confusing
4. 
5. pattern seems in alphabhetical oder but several models are missing (v,w,x,y), 
6. fct time seems not in the correct order, with steps
    array([ 3,  6,  9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 54, 60, 66, 72, 51]
7. absolute vorticity, potential vorticity cant be read with xarray
8. height above ground should be manually extracted
9. absv lev idx 10:13 kosong, pv lev idx 3(400), 7(800), 


'''

#%% static 
inpath = '/mnt/data/input/arpege/'

#%% Functions
def gendts():
    import datetime
    dtn = datetime.datetime.utcnow() - timedelta(hours=3)
    if (dtn.hour >= 0) & (dtn.hour < 6):
        cycl = 0
    elif (dtn.hour >= 6) & (dtn.hour < 12):
        cycl = 6
    elif (dtn.hour >= 12) & (dtn.hour <= 18):
        cycl = 12
    elif (dtn.hour >= 18) & (dtn.hour <= 24):
        cycl = 18
    
    dtc = datetime.datetime(dtn.year,dtn.month,dtn.day,cycl)
    return dtc

# A_YXGR01ECMA060000_C_ECMF_20231206000000_D1S12060000120601001.bz2.bin, A_YXGR02ECMA060000_C_ECMF_20231206000000_D1S12060000120602001.bz2.bin


def fileNameGen(dtn):
    import string
    alphabet_list = [letter for letter in string.ascii_uppercase]
    filenames = [f'A_YMJ{x}89ARPM{dtn:%d}{dtn:%H}00_C_LFPW_------{dtn:%d}{dtn:%H}00--.bin' for x in alphabet_list]   

    return filenames

def fileDownload(grbzname, dtn, maxretry : int=1,retryfreq : int = 3):    #nanti ganti 60 aja untuk debug 3 detik  
    import time
    with ftputil.FTPHost("172.19.0.47", "transmet2023", "99xrAQbW2qN7h") as ftp_host:
        file = f'MODELARPEGE/{grbzname}'
        basename = os.path.basename(file)
        inpathbin = Path(f'{inpath}/{dtn:%Y}/{dtn:%m}/{dtn:%d}_{dtn:%H}')
        fout = Path(inpathbin, basename)
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

def isobaricRet (arr,filenames):
    arrval = []
    for var in arr:
        ds = xr.open_mfdataset(filenames, engine="cfgrib",
        backend_kwargs=dict(filter_by_keys={"typeOfLevel": "isobaricInhPa", "shortName":var, 'edition': 1}),
        concat_dim='valid_time', combine='nested', parallel=True)

        arrval.append(ds[var].values)
    return arrval

def varret(filename):
    dsr = xr.open_mfdataset(filename, engine="cfgrib",
        backend_kwargs=dict(filter_by_keys={"typeOfLevel": "heightAboveGround",'shortName' : 'r', 'level':2, 'edition': 1}),
        concat_dim='valid_time', combine='nested', parallel=True)
    dsr = dsr.sortby(dsr.valid_time)
    r2 = dsr['r'].values
    dst = xr.open_mfdataset(filename, engine="cfgrib",
        backend_kwargs=dict(filter_by_keys={"typeOfLevel": "heightAboveGround",'shortName' : '2t', 'level':2, 'edition': 1}),
        concat_dim='valid_time', combine='nested', parallel=True)
    dst = dst.sortby(dst.valid_time)
    t2 = dst['t2m'].values
    dsu = xr.open_mfdataset(filename, engine="cfgrib",
        backend_kwargs=dict(filter_by_keys={"typeOfLevel": "heightAboveGround",'shortName' : '10u', 'level':10, 'edition': 1}),
        concat_dim='valid_time', combine='nested', parallel=True)
    dsu = dsu.sortby(dsu.valid_time)
    u10 = dsu['u10'].values
    dsv = xr.open_mfdataset(filename, engine="cfgrib",
        backend_kwargs=dict(filter_by_keys={"typeOfLevel": "heightAboveGround",'shortName' : '10v', 'level':10, 'edition': 1}),
        concat_dim='valid_time', combine='nested', parallel=True)
    dsv = dsv.sortby(dsv.valid_time)
    v10 = dsv['v10'].values
    dsm = xr.open_mfdataset(filename, engine="cfgrib",
        backend_kwargs=dict(filter_by_keys={"shortName": "msl"}),
        concat_dim='valid_time', combine='nested', parallel=True)
    dsm = dsm.sortby(dsm.valid_time)
    msl = dsm['msl'].values
    return r2, t2, u10, v10,msl

def rvortextract(gribname):
    grbs = pygrib.open(gribname)
    avgrb = grbs.select(shortName='absv',typeOfLevel='isobaricInhPa')
    pvgrb = grbs.select(shortName='pv',typeOfLevel='isobaricInhPa')
    avgrb,pvgrb
    av,pv = ([] for i in range(2))
    i = 0
    j = 0
    k=0
    nans = np.empty((np.shape(avgrb[0].values)))
    nans[:] = np.nan
    while i < 13:
        if i<9:
            av.append(avgrb[i].values)
        else:
            av.append(nans)  
        i+=1

    while j<13:
        if k<7:
            pv.append(pvgrb[k].values)
            k+=1
        else:
            pv.append(nans)
        j+=1
    return av,pv

def windextract(u,v):
    from metpy.calc import wind_direction, wind_speed
    return wind_speed(u,v), wind_direction(u,v)


#%% initnc
def initNC():
    now = datetime.now()
    dtn = gendts()
    filenames = fileNameGen(dtn)
    
    for filename in filenames :
        fileDownload(filename, dtn)

    print(f'done downloading, elapsed time:{(datetime.now()-now)}')
    inpathbin = Path(f'{inpath}/{dtn:%Y}/{dtn:%m}/{dtn:%d}_{dtn:%H}')
    gribnames = glob.glob(f"{inpathbin}/*.bin")

    now = datetime.now()
    gribnames.sort()
    gribnames = gribnames[1:]

    #ds will also be used as the dataset frame
    ds = xr.open_mfdataset(gribnames, engine="cfgrib",
        backend_kwargs=dict(filter_by_keys={"typeOfLevel": "surface", 'edition': 1}),
        concat_dim='valid_time', combine='nested', parallel=True)
    
    ds = ds.sortby(ds.valid_time)


    isobaric = ['r', 't', 'u', 'v', 'w', 'z'] # w : vertivorti, 

    dsforlev = xr.open_mfdataset(gribnames[0], engine="cfgrib",
        backend_kwargs=dict(filter_by_keys={"typeOfLevel": "isobaricInhPa", "shortName":'r', 'edition': 1}),
        concat_dim='valid_time', combine='nested', parallel=True)
    levels = dsforlev.isobaricInhPa.values
    ds = ds.assign_coords(level=("level" ,levels))
    dsforlev.close()

    arrval = isobaricRet(isobaric,gribnames)
    for var in enumerate(isobaric):
        idx = var[0]
        ds[var[1]] = (('valid_time','level','latitude','longitude'), arrval[idx])
    
    # ds vorticity
    # levels = [ 100.,  150.,  200.,  250.,  300.,  400.,  500.,  600.,  700., 800.,  850.,  900.,  925.,  950., 1000.]
    
    dumdum = np.empty((len(ds.valid_time),len(levels),len(ds.latitude),len(ds.longitude)))
    dumdum[:,:,:,:] = np.nan
    var4 = ['av','pv']
    for var in var4:
        ds[var] =(('valid_time','level','latitude','longitude'),dumdum)

    for gribname in enumerate(gribnames):
        av,pv = rvortextract(gribname[1])
        ds['av'][gribname[0],:,:,:] = av
        ds['pv'][gribname[0],:,:,:] = pv


    r2, t2, u10, v10,msl = varret(gribnames)

    ds['rh2'] = (('valid_time','latitude','longitude'),r2)
    ds['t2'] = (('valid_time','latitude','longitude'),t2)
    ds['u10'] = (('valid_time','latitude','longitude'),u10)
    ds['v10'] = (('valid_time','latitude','longitude'),v10)
    ds['msl'] = (('valid_time','latitude','longitude'),msl)


    print(f'done processing, elapsed time:{(datetime.now()- now)}')
    ds = ds.drop_vars(['time','step','surface','p3099']).rename_dims({'valid_time':'time'}).rename_vars({'valid_time':'time', 'latitude':'lat','longitude':'lon',
                                                                                                         'CAPE_INS':'cape','r':'rh','w':'vv','msl':'mslp'})
    ds['u'].attrs['units'] = 'm s**-1'
    ds['v'].attrs['units'] = 'm s**-1'
    wspd, wdir = windextract(ds['u'],ds['v'])
    ds['wspd'] = (('time','level','lat','lon'),wspd)
    ds['wdir'] = (('time','level','lat','lon'),wdir)

    now = datetime.now()
    # EXPORTING
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

    fno = f'/mnt/data/output/arpege/arpege_{dtn:%Y%m%d}_{dtn:%H%M}.nc'
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