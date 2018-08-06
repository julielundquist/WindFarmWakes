#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 14:39:05 2017

@author: Jessica Tomaszewski
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap

plt.close('all')

# U turbine data
datafile = open('TX_turbines_U.csv', 'r')
data = datafile.readlines()[1:]
data = [x.split(',') for x in data]

lats_U = np.array([x[12] for x in data])
lons_U = np.array([x[13] for x in data])

# D turbine data
datafile = open('TX_turbines_D.csv', 'r')
data = datafile.readlines()[1:]
data = [x.split(',') for x in data]

lats_D = np.array([x[12] for x in data])
lons_D = np.array([x[13] for x in data])

# C turbine data
datafile = open('TX_turbines_C.csv', 'r')
data = datafile.readlines()[1:]
data = [x.split(',') for x in data]

lats_C = np.array([x[12] for x in data])
lons_C = np.array([x[13] for x in data])

# basic simulation data
nc_f = 'wrfout_d04_2012-12-31_12:00:00'
path = ''
nc_fid = Dataset(path+nc_f, 'r')  
            
lat = nc_fid.variables['XLAT'][:]
lon = nc_fid.variables['XLONG'][:]
time = nc_fid.variables['Times'][:]

lat_stag = 0.5*(lat[:,:,:-1]+lat[:,:,1:])
lon_stag = 0.5*(lon[:,:-1,:]+lon[:,1:,:])

# simulation data for 1km
power = nc_fid.variables['POWER'][:]

# elevation data
url = 'http://ferret.pmel.noaa.gov/pmel/thredds/dodsC/data/PMEL/etopo2.nc'
etopodata = Dataset(url)

topoin = etopodata.variables['rose'][:]
print 'done with topoin'
e_lon = etopodata.variables['etopo2_x'][:]
e_lat = etopodata.variables['etopo2_y'][:]
print 'done with lat/lon'

# font sizes
cb_fs = 18
lb_fs = 22
tk_fs = 14
bx_fs = 20

lat_min = np.min(lat)
lat_max = np.max(lat)
lon_min = np.min(lon)
lon_max = np.max(lon)

lllon=-100.8
lllat=32.35
urlon=-100.52
urlat=32.55


fig = plt.figure(figsize=[14,8])
ax = plt.subplot(111)

m = Basemap(projection='lcc', lat_0=(lat_min+lat_max)/2, lon_0=(lon_min+lon_max)/2,
    resolution = 'h', area_thresh = 10, epsg = 32039,
    llcrnrlon=lllon, llcrnrlat=lllat,  
    urcrnrlon=urlon, urcrnrlat=urlat)

m.drawcounties
m.drawmapboundary()

paral_values = np.arange(32.371-(0.027/2),32.8,0.027)
merid_values = np.arange(-100.794-(0.0325/2),-100.4,0.0325)

parallels = m.drawparallels(paral_values,labels=[1,0,0,0],fontsize=tk_fs+1,fmt='%.2f',linewidth=0.0)#,rotation=90) 
meridians = m.drawmeridians(merid_values,labels=[0,0,0,1],fontsize=tk_fs+1,fmt='%.2f',linewidth=0.0)

nx = int((m.xmax-m.xmin)/10.)+1; ny = int((m.ymax-m.ymin)/10.)+1
topodat = m.transform_scalar(topoin,e_lon,e_lat,nx,ny)

# plot image over map with imshow.
im = m.imshow(topodat, cmap='gist_earth')#, alpha=0.8)

# for pcolormesh plotting        
lon_shift = 0.5*(lon[0,:,1:] + lon[0,:,:-1])
lat_shift = 0.5*(lat[0,1:,:] + lat[0,:-1,:])
x, y = m(lon_shift[:-1,:], lat_shift[:,:-1])

CS = m.pcolormesh(x,y,power[0,1:,1:]*0, vmax=20, vmin=0, cmap='binary', alpha=0.3,
                edgecolors='k',linestyle=':',lw=0.5)   

# plot wind turbines
x_farm, y_farm = m(lons_U, lats_U)
m.plot(x_farm, y_farm, '.', ms=6, c='#f7fcb9', label='"Upwind": Loraine')

# plot wind turbines
x_farm, y_farm = m(lons_D, lats_D)
m.plot(x_farm, y_farm, '.', ms=6, c='#a50f15', label='"Downwind": Roscoe')

# plot wind turbines
x_farm, y_farm = m(lons_C, lats_C)
m.plot(x_farm, y_farm, '.', ms=6, c='k', label='"Control": Champion')   

m.drawmapscale(-100.739, 32.362, (lllon+urlon)/2., (lllat+urlat)/2., 10, barstyle='fancy', 
               yoffset=300, units='km', fontsize=12, fillcolor1='w', fillcolor2='0.1',
               zorder=10)

cb = plt.colorbar(im, format='%.0f', pad=0.04)
cb.set_label(r'Elevation  (m)', fontsize=cb_fs, labelpad=-78)
cb.ax.tick_params(labelsize=tk_fs) 
plt.tick_params(axis='both', which='major', labelsize=tk_fs)

leg = plt.legend(loc='upper right', numpoints=1, fontsize=14, 
            handletextpad=0.1, borderaxespad=0.3, markerfirst=False)
frame = leg.get_frame()
frame.set_alpha(0.6)

x0 = 1500
y0 = 2200

ax.annotate('', xy=(x0+3112,y0+1500), xytext=(x0,y0),
            arrowprops=dict(facecolor='k', width=2))

plt.tick_params(direction='inout', length=6)

for m in meridians:
    try:
        meridians[m][1][0].set_rotation(25)
    except:
        pass

#plt.savefig('zoom_elevation_map_1km_arrow.pdf', bbox_inches='tight', dpi=1200)

