# -*- coding: utf-8 -*-
'''
Created on Jul 13, 2022
Written by Meng Gao, revised by K
@author: 31417
'''
from pylab import *
import netCDF4 as nc
import scipy.io
import numpy as np
import pandas
from shapely.geometry import Point
import geopandas as gpd
import shapely
import matplotlib.pyplot as plt
import matplotlib 
import cartopy
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
import matplotlib.cm as cm
from sf import *
ds_chl = nc.Dataset('C:\\Users\\mathp\\Downloads\\erdMH1chlamday_Lon0360_6ff8_fcd0_97cc_U1712153979623.nc')

chl_tmp = ds_chl.variables["chlorophyll"][:]
chl = np.nanmean(chl_tmp,0)


## getting the k490...
# ds_k490 = nc.Dataset('C:\\Users\\mathp\\Downloads\\erdMH1kd490mday_Lon0360_d7f9_721c_abcb_U1712155772765.nc')
# k490_tmp = ds_chl.variables["chlorophyll"][:]
# k490 = np.nanmean(k490_tmp,0)
# "C:\Users\mathp\Downloads\erdVH2018parmday_Lon0360_e8ab_da17_06c7.nc"
# "C:\Users\mathp\Downloads\erdVH2018parmday_6679_a71c_226a.nc"
ds_par = nc.Dataset('C:\\Users\\mathp\\Downloads\\erdVH2018parmday_Lon0360_e8ab_da17_06c7.nc')
par = np.flip(np.nanmean(ds_par.variables['par'][:],0),0)
# Photosynthetically Available Radiation, R. Frouin, einstein m^-2 day^-1
def Ical():
    I0 = genfromtxt("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization\\LightData.csv",delimiter=',').T
    new_light = np.zeros([180,90])
    for i in np.arange(180):
        for j in np.arange(90):
            new_light[i,j] = I0[i*2,j*2]+I0[i*2,j*2+1]+I0[i*2+1,j*2]+I0[i*2+1,j*2+1]/4
    new_light2=np.flip(new_light,1)
    zn=12
    xn=180
    yn=90
    I = zeros([zn,xn,yn])*nan
    z = genfromtxt("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization\\Depth.csv",delimiter=',')
    C = np.flip(np.transpose(chl),1)
    k = 0.121*C**0.428
    # k =0.0166+ 0.077298*C**0.67155
    # k = np.flip(np.transpose(k490),1)
    for i in range(size(z)):
        I[i] = new_light* exp(-k*z[i]) #multiply by k instead of dividing by z0
        # I[i] = np.transpose(par) * exp(-z[i]*k) #multiply by k instead of dividing by z0
    return I

I = Ical()


# cmap1 = cm.get_cmap("RdYlBu_r",lut=20)
# fig, ax = plt.subplots(figsize = [12,8],subplot_kw={'projection':ccrs.PlateCarree()})
# x,y = np.meshgrid(lon_c, lat_c)
# im = ax.pcolormesh(lon_c, lat_c, np.transpose(I[0,:,:]),cmap=cmap1)
# cbar = plt.colorbar(im, ax = ax, orientation = 'vertical', fraction = 0.02, pad = 0.02)
# cbar.set_label('Attenuation coefficient  ($m^{-1}$)')
# ax.set_xticks(np.arange(-180,181,60), crs=ccrs.PlateCarree())
# lon_formatter = cticker.LongitudeFormatter()
# ax.xaxis.set_major_formatter(lon_formatter)
# ax.set_yticks(np.arange(-90,91,30), crs=ccrs.PlateCarree())
# lat_formatter = cticker.LatitudeFormatter()
# ax.yaxis.set_major_formatter(lat_formatter)
# ax.coastlines()
# sf("k_map_without_chla")