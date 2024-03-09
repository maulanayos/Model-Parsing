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

1. gribs unordered and cannot be accessed with cfgrib engine
2. gribs validDate is not working and it refers to analdate instead
3. different level on r,t variables with the rest
    r,t -> include lvl 600, 800, 900, 950
    skipped index level for the rest 7,9,11,13

Update :
1. sort ketinggalan, added after 13 Feb 2024 4.54 utc
2. per 19 Feb upload data maju 12 jam
'''

#%% static 
inpath = '/mnt/data/input/ecmwf_v2/zip'
outpath = '/mnt/data/input/ecmwf_v2/unzipped'
ds3demand = ['blh','cape','cbh','cp','deg0l','fg10','gh','lcc','lsp','mcc','tp','u10','v10','p3020']
ds4demand = ['u','v']
ds3demand2 = ['ceil']
niteration = []
niter = 1
# 3 days try
while niter < 241:
    niteration.append(niter)
    if niter < 90 :
        niter = niter + 1
    elif niter < 144:
        niter = niter + 3
    else:
        niter = niter + 6


def gendts():
    import datetime
    dtn = datetime.datetime.utcnow()
    if (dtn.hour >= 0) & (dtn.hour < 12):
        cycl = 0
    elif (dtn.hour >= 12) & (dtn.hour <= 24):
        cycl = 12
    
    dtc = datetime.datetime(dtn.year,dtn.month,dtn.day,cycl)
    return dtc - timedelta(hours=12)

# A_YXGR01ECMA060000_C_ECMF_20231206000000_D1S12060000120601001.bz2.bin, A_YXGR02ECMA060000_C_ECMF_20231206000000_D1S12060000120602001.bz2.bin

def fileNameGen(baserun,fct):
    year = baserun.strftime('%Y')
    mon = baserun.strftime('%m')
    day = baserun.strftime('%d')
    cycl = baserun.strftime('%H')

    fctrun = baserun + relativedelta(hours = fct)

    ftcmon = fctrun.strftime('%m')
    fctday = fctrun.strftime('%d')
    fcthour = fctrun.strftime('%H')

    if fct < 100:
        fctstr = 'R%02i' % fct
    elif fct < 200:
        fctstr = 'S%02i' % (fct - 100)
    elif fct < 300:
        fctstr = 'T%02i' % (fct - 200)
    
    if cycl == '00':
        cyclstr = 'ECMA'
    elif cycl == '12':
        cyclstr = 'ECMC'

    if fct == 00:
        A1 = "H"
        DDD = "D1S"
        D = "1"
    elif fct <= 90:
        A1 = "G"
        DDD = "D1S"
        D = "0"
    else:
        A1 = "E"
        DDD = "D1D"
        D = "0"

    grbzname = f'A_YX{A1}{fctstr}{cyclstr}{day}{cycl}00_C_ECMF_{year}{mon}{day}{cycl}0000_{DDD}{mon}{day}{cycl}00{ftcmon}{fctday}{fcthour}0{D}1.bz2.bin'
    grbname = f'A_YX{A1}{fctstr}{cyclstr}{day}{cycl}00_C_ECMF_{year}{mon}{day}{cycl}0000_{DDD}{mon}{day}{cycl}00{ftcmon}{fctday}{fcthour}0{D}1.bin'

    return grbzname, grbname, fctrun

def filedownload(grbzname, dtn, maxretry : int=10,retryfreq : int = 1800):      
    with ftputil.FTPHost("172.19.0.47", "transmet2023", "99xrAQbW2qN7h") as ftp_host:
        file = f'MODELECMWF/{grbzname}'
        basename = os.path.basename(file)
        inpathbin = Path(f'{inpath}/{dtn:%Y}/{dtn:%m}/{dtn:%d}')
        fout = Path(inpathbin, basename)
        fout.parent.mkdir(parents=True, exist_ok=True)
    
        if fout.exists():
            return fout
        else:
            for i in range(maxretry):
                if ftp_host.path.exists(file):
                    ftp_host.download(file,fout)
                    # print(f'file downloaded {fout}')
                    return fout
                else:
                    print(f'File {file} not found')
                time.sleep(retryfreq)
            print(f'File {file} not found after {maxretry} retries')
            return None

def bunzip2_file(f):
    fin = f[0]
    fout = f[1]
    fout.parent.mkdir(parents=True, exist_ok=True)
    with open(fin, 'rb') as f_in, open(fout, 'wb') as f_out:
        f_out.write(bz2.decompress(f_in.read()))

def sortFunc(e):
    return e['level'], e['shortName']

def rtextract(grms):
    r = []
    t = []
    i=0
    while i < 30:
        r.append(grms[i].values)
        t.append(grms[i+1].values)
        i+=2
    t = np.array(t) - 273.15
    return r, t

def resextract(gribnames, var, dumdum):
    dstemp = xr.open_mfdataset(gribnames, engine="cfgrib",
                            backend_kwargs=dict(filter_by_keys={"shortName":var,'typeOfLevel': 'isobaricInhPa'}),
                            concat_dim='valid_time', combine='nested', parallel=True)
    val = dstemp[var].values
    new_val = np.concatenate((val[:,0:1,:,:], dumdum, val[:,1:2,:,:], dumdum, val[:,2:3,:,:], dumdum, val[:,3:4,:,:], dumdum,val[:,4:,:,:] ), axis=1)
    return new_val
def initNC():
    now = datetime.now()
    zippedfiles = []
    bins = []
    dts = []
    dtn = gendts()
    for iter in niteration :
        zipped,_, _ = fileNameGen(dtn,iter)
        filedownload(zipped,dtn)
    for iter in niteration : 
        zipped,ncs, dtf = fileNameGen(dtn,iter)
        inpathbin = Path(f'{inpath}/{dtn:%Y}/{dtn:%m}/{dtn:%d}')
        zippedin = inpathbin / zipped
        outbin = Path(f'{outpath}/{dtn:%Y}/{dtn:%m}/{dtn:%d}/{dtn:%H}')
        binout = outbin / ncs
        zippedfiles.append(zippedin)
        bins.append(binout)
        dts.append(dtf)
    fiter = []
    for i in range(len(zippedfiles)):
            fiter.append([zippedfiles[i], bins[i]])

    with ProcessPoolExecutor(max_workers=24) as executor:
        executor.map(bunzip2_file, fiter)

    
    gribnames = glob.glob(f"{outbin}/*.bin")
    print(f'done unzipping, elapsed time:{(datetime.now()-now)}')

    xr.set_options(keep_attrs=False)

    now = datetime.now()
    # Dataset build initiation
    ds4 = xr.open_mfdataset(gribnames, engine="cfgrib",
        backend_kwargs=dict(filter_by_keys={"shortName":'r'}),
        concat_dim='valid_time', combine='nested', parallel=True)
    
    lats, lons = ds4['latitude'].values, ds4['longitude'].values
    timestep = ds4['valid_time']
    
    dst = xr.open_mfdataset(gribnames, engine="cfgrib",
        backend_kwargs=dict(filter_by_keys={"shortName":'t'}),
        concat_dim='valid_time', combine='nested', parallel=True)
    
    ds4 = xr.merge([ds4,dst], compat='override')
    dst.close()

    dumdum = np.empty((len(timestep),1,len(lats),len(lons)))
    dumdum[:,:,:,:] = np.nan

    var4 = ['gh', 'pv', 'q', 'u', 'v', 'w', 'z']
    
    for var in var4:
        values = resextract(gribnames,var,dumdum)
        ds4[var] = (('valid_time','isobaricInhPa','latitude','longitude'),values)
        

    ds4 = ds4.drop_vars(['time','step','number']).rename_dims({'valid_time':'time','isobaricInhPa':'level'}).rename_vars({'valid_time':'time','isobaricInhPa':'level', 'latitude':'lat','longitude':'lon'})

    ds3 = xr.open_mfdataset(gribnames, engine="cfgrib",
        backend_kwargs=dict(filter_by_keys={"typeOfLevel": "surface", 'edition': 1}),
        concat_dim='valid_time', combine='nested', parallel=True)
    # ds3 = ds3[ds3demand]
    ds3 = ds3.drop_vars(['time','step','surface','number']).rename_dims({'valid_time':'time'}).rename_vars({'valid_time':'time', 'latitude':'lat','longitude':'lon'})

    ds32 = xr.open_mfdataset(gribnames, engine="cfgrib",
        backend_kwargs=dict(filter_by_keys={"typeOfLevel": "surface", 'edition': 2}),
        concat_dim='valid_time', combine='nested', parallel=True)
    # ds32 = ds32[ds3demand2]
    ds32 = ds32.drop_vars(['time','step','surface']).rename_dims({'valid_time':'time'}).rename_vars({'valid_time':'time','latitude':'lat','longitude':'lon'})
    nds = xr.merge([ds4,ds3,ds32], compat='override')
    nds = nds.sortby("time")

    nds = nds.drop_vars(['u100','v100','fzra','cin','sst','q'])
    print(f'done processed ncdf, elapsed time: {(datetime.now() - now)}')


    
    now = datetime.now()
    # EXPORTING
    chunk3d = [3,6,12]
    chunk4d = [3,3,6,12]
    encoding = {}
    encoding_keys = ("_FillValue", "dtype", "scale_factor", "add_offset", "grid_mapping", "chunksizes")
    for data_var in nds.data_vars:
        if len(nds[data_var].shape) == 3:
            chunks = chunk3d
        elif len(nds[data_var].shape) == 4:
            chunks = chunk4d

        encoding[data_var] = {key: value for key, value in nds[data_var].encoding.items() if key in encoding_keys}
        encoding[data_var].update(zlib=True, complevel=5, chunksizes = chunks)

    fno = f'/mnt/data/output/ECMWF_V2/{dtn:%Y}/{dtn:%m}//ecmwf_{dtn:%Y%m%d}_{dtn:%H%M}.nc'
    fno = Path(fno)
    # save to netcdf
    fno.unlink(missing_ok=True)
    fno.parent.mkdir(parents=True, exist_ok=True)
    nds.to_netcdf(
        fno,
        'w',
        format = "NETCDF4",
        encoding = encoding
    )
    print(f'done saving, elapsed time:{(datetime.now() - now)}')
    
if __name__ == '__main__':

    initNC()