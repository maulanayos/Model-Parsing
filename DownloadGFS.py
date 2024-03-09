"""
On development
@mamyos
"""

import os
import datetime
import numpy as np
import numexpr as ne
from time import sleep
from pathlib import Path
import urllib.request
import pygrib
import pandas as pd
import xarray as xr


# from tools.db_insert import mongoinsert
# import sentry_sdk
# sentry_sdk.init("https://0c2dbdaf918b4494b0f905b3a86d8eab@sentry.circlegeo.com/16")

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
maxtry = 100 # Total 8.3 Hours
waitsec = 300 # 5 minutes
tmpdown = Path('/mnt/data/input/gfs')
optdir = Path('/mnt/data/output/gfs_indo_v2/')
chunk2d = [3,6,12]
chunk3d = [3,3,6,12]
""" DOWNLOAD DATA """
# https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?file=gfs.t00z.pgrb2.0p25.anl&lev_1000_mb=on&lev_100_mb=on&lev_10_m_above_ground=on&lev_150_mb=on&lev_200_mb=on&lev_250_mb=on&lev_2_m_above_ground=on&lev_300_mb=on&lev_400_mb=on&lev_500_mb=on&lev_600_mb=on&lev_700_mb=on&lev_800_mb=on&lev_850_mb=on&lev_900_mb=on&lev_925_mb=on&lev_950_mb=on&lev_entire_atmosphere=on&lev_mean_sea_level=on&lev_surface=on&var_APCP=on&var_CAPE=on&var_GUST=on&var_HCDC=on&var_HGT=on&var_LCDC=on&var_LFTX=on&var_MCDC=on&var_PRATE=on&var_RH=on&var_TCDC=on&var_TMP=on&var_VIS=on&leftlon=95&rightlon=141&toplat=6&bottomlat=-11&dir=%2Fgfs.20231117%2F00%2Fatmos

# / format link gfsurl+[initialtime]+midsufgfsfile+gfsopt+(path)
#(path) = 'dir=%2Fgfs.+[tahun]+[bulan]+[tanggal]+%2F+[initialtime]+%2Fatmos
    
#%% download func
def gendts():
    import datetime
    dtn = datetime.datetime.utcnow()

    if (dtn.hour >= 3) & (dtn.hour < 15):
        iTime = '00'
    elif (dtn.hour >= 15) & (dtn.hour <= 24):
        iTime = '12'
    else:
        iTime = '12'
        dtn = dtn - datetime.timedelta(hours = 12)
    
    
    tabutang = dtn.strftime('%Y%m%d')
    
    return tabutang,iTime

def dldURLCreator(iTime:str,tabutang:str):
   
    gfsurl = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?file=gfs.t'
    gfsopt = '&lev_1000_mb=on&lev_100_mb=on&lev_10_m_above_ground=on&lev_150_mb=on&lev_200_mb=on&lev_250_mb=on&lev_2_m_above_ground=on&lev_300_mb=on&lev_400_mb=on&lev_500_mb=on&lev_600_mb=on&lev_700_mb=on&lev_800_mb=on&lev_850_mb=on&lev_900_mb=on&lev_925_mb=on&lev_950_mb=on&lev_entire_atmosphere=on&lev_mean_sea_level=on&lev_surface=on&var_APCP=on&var_CAPE=on&var_GUST=on&var_HCDC=on&var_HGT=on&var_UGRD=on&var_VGRD=on&var_LCDC=on&var_LFTX=on&var_MCDC=on&var_PRATE=on&var_RH=on&var_TCDC=on&var_TMP=on&var_VIS=on&subregion=&leftlon=70&rightlon=170&toplat=40&bottomlat=-25&'
    midsufgfsfile = 'z.pgrb2.0p25.'
    extension = [] 
    for i in range(1, 121):
        tmpext = 'f'+str(i).zfill(3)
        extension.append(tmpext)
        urls = []
    for i in range(len(extension)):
        path = 'dir=%2Fgfs.' + (tabutang)+ '%2F'+(iTime)+'%2Fatmos'
        url = gfsurl + (iTime)+midsufgfsfile+extension[i]+gfsopt+path
        urls.append(url)
    return(urls)

def gfsdown(urls,itime):
    j = 0
    for j in range(len(urls)):
        gfsdir = urls[j][654:662]
        cycdir = urls[j][65:68]
        gfsfile = urls[j][80:85]
        fname = 'gfs.'+cycdir+'z.pgrb2.0p25'+gfsfile
        if itime == '00':
            gfslocalfile = tmpdown / 'cy00'/ gfsdir / fname
        else :
            gfslocalfile = tmpdown / 'cy12'/ gfsdir / fname

        gfslocalfile.parent.mkdir(parents=True, exist_ok=True)
        
        print(gfslocalfile)
        print('temp file : %s' % gfslocalfile)
        if not os.path.exists(gfslocalfile):
            for i in range(maxtry):
                try:
                    response = urllib.request.urlopen(urls[j], timeout = 5)
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
        else :
            print ("file aready downloaded"+str(gfslocalfile))
        j+=1
    # return(gfslocalfile)

def mstoknot(var):
    return ne.evaluate('var*1.9438445')

def KtodegC(var):
    return ne.evaluate('var-273.15')

def patombar(var):
    return ne.evaluate('var*0.01')

def remask(var,varmask):
    import numpy.ma as ma
    masktemp = varmask.mask
    var = ma.masked_array(var,mask=masktemp)
    return var

param2d = [
    ['2r', 'rh2', None, '%',None],                  # Relative humidity at 2m
    ['2t', 't2', KtodegC, 'degC',None],             # Temperature at 2m
    ['10u', 'u10', None, 'm s**-1',None],           # U-component of wind at 10m
    ['10v', 'v10', None, 'm s**-1',None],           # V-component of wind at 10m
    ['vis', 'vis', None, 'm',None],                 # Visibility
    ['gust', 'gust', None, 'm s**-1',None],         # Wind gust
    ['prate', 'prate', None, 'kg m**-2 s**-1','avg'],# Precipitation rate
    ['tp', 'rr', None, 'mm',None],                  # Total precipitation
#     ['acpcp', 'acpcp', None, 'kg m**-2'],      # Accumulated convective precipitation
    ['lftx', 'lftx', None, 'K',None],               # Surface lifted index
    ['cape', 'cape', None, 'J kg**-1',None],        # Convective Available Potential Energy
    ['tcc', 'tcc_atmospheric', None,'%','avg']
]

param3d = [
    ['u', 'u', None, 'm s**-1'],      # U-component of wind
    ['v', 'v', None, 'm s**-1'],      # V-component of wind
    ['gh', 'gh', None, 'gpm'],        # Geopotential Height
    ['t', 't', KtodegC, 'C'],            # Temperature
    ['r', 'r', None, '%'],            # Relative Humidity
    ['tcc', 'tcc_isobaric', None, '%'],        # Total Cloud Cover
]
#%%encode
encode = {"rh2": {'zlib' : True, "complevel":5, "chunksizes":chunk2d },
          "t2": {'zlib' : True, "complevel":5,"chunksizes":chunk2d},
          "u10": {'zlib' : True, "complevel":5,"chunksizes":chunk2d},
          "v10": {'zlib' : True, "complevel":5,"chunksizes":chunk2d},
          "vis": {'zlib' : True, "complevel":5,"chunksizes":chunk2d},
          "gust": {'zlib' : True, "complevel":5,"chunksizes":chunk2d},
          "prate": {'zlib' : True, "complevel":5,"chunksizes":chunk2d},
          "rr": {'zlib' : True, "complevel":5,"chunksizes":chunk2d},
          "lftx": {'zlib' : True, "complevel":5,"chunksizes":chunk2d},
          "cape": {'zlib' : True, "complevel":5,"chunksizes":chunk2d},
          "tcc_atmospheric": {'zlib' : True, "complevel":5,"chunksizes":chunk2d},
          "u": {'zlib' : True, "complevel":5,"chunksizes":chunk3d},
          "v": {'zlib' : True, "complevel":5,"chunksizes":chunk3d},
          "gh": {'zlib' : True, "complevel":5,"chunksizes":chunk3d},
          "t": {'zlib' : True, "complevel":5,"chunksizes":chunk3d},
          "r": {'zlib' : True, "complevel":5,"chunksizes":chunk3d},
          "tcc_isobaric": {'zlib' : True, "complevel":5,"chunksizes":chunk3d}
          } 
#%% main
def InitNc():
    tabutang,itime = gendts()
    # tabutang, itime = '20231128', '12'
    urls = dldURLCreator(itime, tabutang)
    gfsdown(urls,itime)
    nclocalfile = optdir / tabutang[0:4] / tabutang[4:6]
    if itime == '00':
        ncname = str(nclocalfile)+'/'+'gfs_'+tabutang+'_0000.nc' # masih statis
    else :
        ncname =str(nclocalfile)+'/'+'gfs_'+tabutang+'_1200.nc' # masih statis
    if not os.path.exists(ncname):
        gribdir = Path('/mnt/data/input/gfs/' +'cy'+itime+'/'+tabutang+'/')
        filenames = os.listdir(gribdir)
        grbs = pygrib.open(str(gribdir / filenames[0]))
        gsel = grbs.select(shortName = param2d[1][0]) #rh
        g = gsel[0]
        lat, lon = g.latlons()
        lat = lat[:,0]
        lon = lon[0,:]
        dtni = g.analDate + datetime.timedelta(hours = 1)
        tunits = 'minutes since %s' % dtni.strftime('%Y-%m-%d %H:%M') #set attrib later on
        levels = []

        for filename in filenames:
            grbs = pygrib.open(str(gribdir / filename))
        gselz = grbs.select(shortName = 'u')
        for lev in gselz:
            levels.append(float(lev.level))
        levels = np.unique(levels)  

        timestep = pd.date_range(start=dtni, periods=120, freq='H') # <---DIGANTI BERAPA PERIODS

        ds = xr.Dataset(coords=dict(
            lon=(["lon"], lon),
            lat=(["lat"], lat),
            level=(["level"], levels),
            time=(["time"],timestep)))
        
        for allvar in param2d:
            dummy = np.empty([len(timestep),len(lat),len(lon)])
            ds[allvar[1]] = (('time','lat','lon'),dummy)
        for allvar in param3d:
            dummy = np.empty([len(timestep),len(levels), len(lat), len(lon)])
            ds[allvar[1]] = (('time','level','lat','lon'),dummy)
        timeid = 0
        for filename in filenames:
            grbpath = str(gribdir / filename)
            grbs = pygrib.open(grbpath)
            for allvar in param2d:
                grbs.seek(0)
                if allvar[4] != None:
                    gsel = grbs.select(shortName = allvar[0],stepType=allvar[4])
                else:
                    gsel = grbs.select(shortName = allvar[0]) #rh
                g = gsel[0]
                var = g.values
                converter = allvar[2]
                if converter:
                    if type(var) is np.ma.core.MaskedArray:
                        vartemp = var
                        var = converter(var)
                        var = remask(var,vartemp)
                    else:
                        var = converter(var)
                ds[allvar[1]][timeid,:,:] = var #nanti indexingnya benerin buat for loop

            for allvar in param3d:
                grbs.seek(0)
                data = np.empty([len(levels), len(lat), len(lon)])

                g_sel = grbs.select(shortName = allvar[0],
                                    typeOfLevel = 'isobaricInhPa')
                converter = allvar[2]
                
                levid = 0
                if converter:      
                    for gselsel in g_sel:
                        data = gselsel.values
                        if type(data) is np.ma.core.MaskedArray:
                            vartemp = data
                            data = converter(data)
                            data = remask(data,vartemp)
                        else:
                            data = converter(data)
                        ds[allvar[1]][timeid,levid,:,:]=data
                        levid+=1
                else : 
                    for gselsel in g_sel:
                        data = gselsel.values
                        ds[allvar[1]][timeid,levid,:,:]=data
                        levid+=1
            timeid+=1
        
        if not os.path.exists(nclocalfile):
            os.makedirs(nclocalfile)
        if itime == '00':
            ds.to_netcdf(ncname,encoding = encode) # masih statis
            print(ncname + 'created')
        else :
            ds.to_netcdf(ncname,encoding = encode) # masih statis
            print(ncname + 'created')
    else:
        print("nc already downloaded")



# %%
if __name__=="__main__":
    InitNc()