

import xarray as xr
from netCDF4 import Dataset
import wrf
from wrf import getvar, interplevel
import numpy as np
import numexpr as ne
        
"""
relative vorticity & prate needs some work
"""
#%% Units and derived variables functions
def qketotke (var):
    return ne.evaluate('var/2')

def zdestag (var) :
    var = wrf.destagger(var,1)
    return var

def latdestag (var) :
    var = wrf.destagger(var,2)
    return var

def londestag (var) :
    var = wrf.destagger(var,3)
    return var

def KtodegC(var):
    return ne.evaluate('var-273.15')

def viscalc(ds):
    
    # Extract required variables
    temp = ds['t'][:,:,:,:] 
    dewpoint = ds['td'][:,:,:,:]
    
    # Calculate vapor pressure and relative humidity
    e = 6.112 * np.exp((17.67*dewpoint)/(dewpoint + 243.5))
    es = 6.112 * np.exp((17.67*temp)/(temp + 243.5)) 
    rh = (e/es) * 100
    
    # Constants
    diameter = 20*10**(-6) # 20 μm  
    gamma = 0.65
        
    # Visibility function 
    visibility = (3.0 / (rh**0.75)) / (1 + ((5.0*10**(-5)) / 
                (diameter**2)))**(1.0/(gamma-1))
    
    return visibility.values


def calculate_wind_gust(wrfout_file):
    import netCDF4
    # Open the NetCDF file
    nc = netCDF4.Dataset(wrfout_file)
    
    # Extract needed variables
    u10 = nc.variables['U10'][:] 
    v10 = nc.variables['V10'][:]
    
    # Calculate the 10m wind speed  
    wind10 = np.sqrt(u10**2 + v10**2)

    # Empirical equation to compute gust
    gust_factor = 0.43 + (0.097 * wind10) 

    # Compute final gust speed  
    surface_gust = gust_factor * wind10
    nc.close()
    surface_gust = surface_gust.data
    return surface_gust

def calculate_ceiling(ds, cloud_threshold=0.8):

    # Calculate cloud cover based on your specific criteria
    cloud_cover = ds['tcc']  # Cloud fraction variable from WRF output

    # Find the height where cloud cover exceeds the threshold
    ceiling_height = ds['hgt'].where(cloud_cover > cloud_threshold).max(dim='level')

    # Close the WRF output file

    return ceiling_height

def calculate_k_index(wrfnc,timeidx):

    t = getvar(wrfnc,"tc",timeidx)
    p = getvar(wrfnc, "pressure",timeidx)
    td = getvar(wrfnc,"td",timeidx)
    t700 = interplevel(t, p, 700.0)
    t850 = interplevel(t, p, 850.0)
    t500 = interplevel(t, p, 500.0) 

    # Extract 850 mb dewpoint  
    td850 = interplevel(td, p, 850.0)
    td700 = interplevel(td, p, 700.0)
    # Calculate K Index 
    k_index = (t850 - t500) + td850 - (t700 - td700)

    return k_index.values

def calculate_prate(wrfnc):
    ncdf = Dataset(wrfnc)  

    # Get accumulated precipitation
    precip_acc = ncdf.variables['RAINNC'][:]  

    # Get time variable  
    time_var = ncdf.variables['Times']  

    # Get time conversion factor 
    acc_period = time_var.units.split()[2]
    factor = float(acc_period[:-1])/float(time_var.units.split()[4])

    # Calculate instantaneous precip rate   
    precip_rate = precip_acc/factor

    return precip_rate.values


def calculate_showalter(wrfnc,timeidx):
    t = getvar(wrfnc,"tc",timeidx)
    p = getvar(wrfnc, "pressure",timeidx)
    td = getvar(wrfnc,"td",timeidx)
    t500 = interplevel(t, p, 500.0)
    t850 = interplevel(t, p, 850.0)
    td850 = interplevel(td, p, 850.0) 

    # Lifted parcel calculation 
    parcel_t = t850 + (8000 - 500)/(3000 - (td850-t850)) * (td850-t850)  

    # Calculate Showalter Index
    showalter_index = t500 - parcel_t

    return showalter_index.values

def calculate_surface_li(wrfnc,timeidx):
    # Specify the pressure levels of interest
    t = getvar(wrfnc,"tk",timeidx)  
    t2 = getvar(wrfnc,"T2",timeidx)      
    p = getvar(wrfnc, "pressure",timeidx)
    # Extract temperature and geopotential height at surface and 500 hPa
    t500 = interplevel(t, p, 500.0)

    # Calculate the potential temperature of a parcel lifted adiabatically from the surface to 500 hPa
    potential_temperature_parcel = t2 * (500 / 1000) ** (287.05 / 1005.0)

    # Calculate the surface lifted index
    surface_lifted_index = potential_temperature_parcel - t500

    return surface_lifted_index.values

def calculate_totaltotal_sweat(wrfnc,timeidx) :
    """
    TT = VT + CT
    VT = T(850 mb) - T(500 mb)
    CT = Td(850 mb) - T(500 mb)
    """
    t = getvar(wrfnc,"tc",timeidx)
    p = getvar(wrfnc, "pressure",timeidx)
    td = getvar(wrfnc,"td",timeidx)
    wspd_wdir = getvar(wrfnc, "wspd_wdir", units="kt")
    wspd = wspd_wdir[0,:,:,:]
    wdir = wspd_wdir[1,:,:,:]

    t500 = interplevel(t, p, 500.0)
    t850 = interplevel(t, p, 850.0)
    td850 = interplevel(td, p, 850.0) 
    f5 = interplevel(wspd, p, 500.0) 
    f8 = interplevel(wspd, p, 850.0)
    fd5 = interplevel(wdir, p, 500.0)
    fd8 = interplevel(wdir, p, 850.0)
    S = np.sin(fd5 - fd8) 

    vt = t850 - t500
    ct = td850 - t500
    tt = vt +ct

    sweat = 12*td850 + (20*(tt - 49)) + (2*f8) + f5 + (125*(S + 0.2))
    return tt.values, sweat.values

def calculate_relative_vorticity (ds,timeidx):
    lats = ds['lat'][:]
    long = ds['lon'][:]

    dlat = lats[1] - lats[0]
    dlon = long[1] - long[0]

    vortarr = np.zeros((len(ds['level']), len(lats),len(long)))

    for levid in range (len(ds['level'])):
        u = ds['u'][timeidx,levid,:,:].values
        v = ds['v'][timeidx,levid,:,:].values
        dudy = np.gradient(u, dlat,dlon)
        dvdx = np.gradient(v, dlat,dlon)
        vort = dvdx - dudy
        vortarr[levid,:,:] = vort

    return vortarr

def calculate_deg0l(ds):
    temperature = ds['t']
    time_dim, vert_dim, lat_dim, lon_dim = temperature.shape 
    # Find 0 C isothermal layers
    iso0_layers = np.zeros((time_dim, vert_dim-1, lat_dim, lon_dim))
    for t in range(time_dim):
        for k in range(vert_dim-1):
            iso0_cond = (temperature[t,k,:,:] == 273.15) & (temperature[t,k+1,:,:] == 273.15)  
            iso0_layers[t, k, :, :][iso0_cond] = 1

    return iso0_layers

"""
K INDEX eval
1. ZNU levels to mb
    -> p = PTOP + Π×(P_sfc - PTOP)
    Where:
    PTOP = pressure at the model top
    P_sfc = pressure at the model surface
    Π = ETA value (ZNU) 

2. 
    
"""

#%%  Demanded Variable Declaration
dmd2ds = [['AFWA_CAPE','cape'],
        ['PBLH','blh'],
        ['CBASEHT','cbh'],
        ['RAINC','cp'],
        ['RAINNC','tp'],
        ['AFWA_MSLP','mslp'],
        ['U10','u10'],
        ['V10','v10'],
        ['PSFC','sfc'],
        ['HGT','hgt']]

dmd3ds = [['P','p',None],
        ['CLDFRA','tcc',None],
        ['T','t', KtodegC],
        ['QKE','tke', qketotke],
        ['PH','g', zdestag],
        ['U', 'u', londestag],
        ['V','v', latdestag]]

wrftl = [['td',None],
        ['avo',None],
        ['pvo',None],
        ['rh',None],
        ['rh2',None],
        ['T2',KtodegC],
        ['td2',None],
        ['twb',KtodegC],
        ['low_cloudfrac',None],
        ['mid_cloudfrac',None],
        ['high_cloudfrac',None]]
#%% Basic variables Retrieval
file_nc = "coba/tes/envyosa/wrfout_d01_2024-01-08_12:00:00"
wrfnc = Dataset(file_nc)

ds = xr.open_dataset(file_nc)
lon, lat, time, levels = ds['XLONG'][0,0,:].values, ds['XLAT'][0,:,0].values, ds['XTIME'].values, ds['ZNU'][0,:].values
dsnew = xr.Dataset(coords=dict(
    lon=(["lon"], lon),
    lat=(["lat"], lat),
    level=(["level"], levels),
    time=(["time"],time)))

#wrftools bag1

for wrftld in wrftl:
    i=0
    fdummy = getvar(wrfnc, wrftld[0])

    if np.shape(fdummy) == (32,311,620):   
        dummy = np.empty((25,32,311,620))
        dsnew[wrftld[0]] = (('time','level','lat','lon'), dummy)
        while i<len(time):
            dsnew[wrftld[0]][i,:,:,:] = getvar(wrfnc,wrftld[0],timeidx=i).values
            i+=1

        converter = wrftld[1]
        if converter :
            var = dsnew[wrftld[0]]
            var = converter(var)
            dsnew[wrftld[0]] = (('time','level','lat','lon'), var)

    else :
        dummy = np.empty((25,311,620))
        dsnew[wrftld[0]] = (('time','lat','lon'), dummy)
        while i<len(time):
            dsnew[wrftld[0]][i,:,:] = getvar(wrfnc,wrftld[0],timeidx=i).values
            i+=1

        converter = wrftld[1]
        if converter :
            var = dsnew[wrftld[0]]
            var = converter(var)
            dsnew[wrftld[0]] = (('time','lat','lon'), var)

# 2D
for dmd2d in dmd2ds:
    dsnew[dmd2d[1]] = (('time','lat','lon'), ds[dmd2d[0]].values)
# 3D
for dmd3d in dmd3ds:
    converter = dmd3d[2]
    if converter :
        data = converter(ds[dmd3d[0]])
        dsnew[dmd3d[1]] = (('time','level','lat','lon'), data)
    else :
        dsnew[dmd3d[1]] = (('time','level','lat','lon'), ds[dmd3d[0]].values)

#%% Derived Variables Retrieval
dsnew['lsp'] = dsnew['tp'] - dsnew['cp']
# dsnew['prate'] = calculate_prate(file_nc)
        
dsnew['vis'] =(('time','level','lat','lon'), viscalc(dsnew)) #vis

dsnew['gust'] = (('time','lat','lon'), calculate_wind_gust(file_nc))


dummy2 = np.empty((25,311,620))
dummy3 = np.empty((25,32,311,620))
dsnew['wdir'] = (('time','level','lat','lon'), dummy3)
dsnew['wind'] = (('time','level','lat','lon'), dummy3)
dsnew['wdir10'] = (('time','level','lat','lon'), dummy3)
dsnew['wind10'] = (('time','level','lat','lon'), dummy3)
dsnew['deg0l'] = (('time','level','lat','lon'), dummy3)
# dsnew['rv'] = (('time','level','lat','lon'), dummy3)
dsnew['ki'] = (('time','lat','lon'),dummy2)
dsnew['ceil'] = (('time','lat','lon'),dummy2)
dsnew['si'] = (('time','lat','lon'),dummy2)
dsnew['sli'] = (('time','lat','lon'),dummy2)
dsnew['tt'] = (('time','lat','lon'),dummy2)
dsnew['sweat'] = (('time','lat','lon'),dummy2)
i=0 
while i<len(time):
        # uvmet = getvar(wrfnc,'uvmet10',timeidx=i)
        dsnew['ki'][i,:,:]=calculate_k_index(wrfnc,i)
        dsnew['si'][i,:,:]=calculate_showalter(wrfnc,i)
        dsnew['sli'][i,:,:]=calculate_surface_li(wrfnc,i)
        dsnew['tt'][i,:,:], dsnew['sweat'][i,:,:] = calculate_totaltotal_sweat(wrfnc,i)    
        wspd_wdir = getvar(wrfnc,'wspd_wdir',timeidx=i)
        wspd_wdir10 = getvar(wrfnc,'wspd_wdir',timeidx=i)
        # dsnew['rv'][i,:,:,:] = calculate_relative_vorticity(dsnew,i)
        dsnew['wind10'][i,:,:,:] = wspd_wdir10[0,:,:,:]
        dsnew['wdir10'][i,:,:,:] = wspd_wdir10[1,:,:,:]
        dsnew['wind'][i,:,:,:] = wspd_wdir[0,:,:,:]
        dsnew['wdir'][i,:,:,:] = wspd_wdir[1,:,:,:]
        i+=1

dsnew['ceil'] = calculate_ceiling(dsnew)
deg0l = calculate_deg0l(dsnew)
dsnew['deg0l'][:,:-1,:,:] = deg0l

dsnew.to_netcdf("tesinanwp.nc")