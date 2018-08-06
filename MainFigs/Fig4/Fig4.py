#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 19:42:45 2017

@author: Jessica Tomaszewski
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Rectangle
import netCDF4 as nc4



#################################
# Read NETCDF of simulation data
path = ''
f = 'fig4_spatial_deficits.nc'
nc_fid = nc4.Dataset(path+f, 'r')  
#################################


#################################
# Extract variables
u_UDC = nc_fid.variables['u_udc'][:]          #['days', 'hours', 'lat', 'lon']
v_UDC = nc_fid.variables['v_udc'][:]          #['days', 'hours', 'lat', 'lon']
u_DC = nc_fid.variables['u_dc'][:]            #['days', 'hours', 'lat', 'lon']
v_DC = nc_fid.variables['v_dc'][:]            #['days', 'hours', 'lat', 'lon']
power_UDC = nc_fid.variables['power_udc'][:]  #['days', 'hours', 'lat', 'lon']
power_DC = nc_fid.variables['power_dc'][:]    #['days', 'hours', 'lat', 'lon']
#################################


#################################
# Read turbine data

# Upwind turbine data
datafile = open('/Users/jeto6273/Documents/CU/Research/CNH/TX_turbines_U.csv', 'r')
data = datafile.readlines()[1:]
data = [x.split(',') for x in data]

lats_U = np.array([x[12] for x in data])
lons_U = np.array([x[13] for x in data])

# Downwind turbine data
datafile = open('/Users/jeto6273/Documents/CU/Research/CNH/TX_turbines_D.csv', 'r')
data = datafile.readlines()[1:]
data = [x.split(',') for x in data]

lats_D = np.array([x[12] for x in data])
lons_D = np.array([x[13] for x in data])

# Control turbine data
datafile = open('/Users/jeto6273/Documents/CU/Research/CNH/TX_turbines_C.csv', 'r')
data = datafile.readlines()[1:]
data = [x.split(',') for x in data]

lats_C = np.array([x[12] for x in data])
lons_C = np.array([x[13] for x in data])
#################################


#################################
# Read one full WRF file to get coordinates

nc_f = 'wrfout_d04_2012-12-31_12:00:00_UDC'
path = '/Users/jeto6273/Documents/CU/Research/CNH/'
nc_fid = Dataset(path+nc_f, 'r')  
            
lat = nc_fid.variables['XLAT'][:]    # [Time, south_north, west_east]
lon = nc_fid.variables['XLONG'][:]   # [Time, south_north, west_east]
time = nc_fid.variables['Times'][:]  # [Time, DateStrLen]
#################################


#################################
# Prep data
wspd_UDC = np.sqrt(u_UDC**2+v_UDC**2)
wspd_DC = np.sqrt(u_DC**2+v_DC**2)

wspd_deficit = wspd_UDC - wspd_DC
power_deficit = power_UDC - power_DC

sum_deficit = []
#################################


#################################
# Make figures

# font sizes
cb_fs = 18
lb_fs = 22
tk_fs = 14
bx_fs = 20

lat_min = np.min(lat)
lat_max = np.max(lat)
lon_min = np.min(lon)
lon_max = np.max(lon)

lllon=-100.99
lllat=32.2
urlon=-100.25
urlat=32.74

d = 23
h = 2


if h < 10:
    hourstring = '0'
else:
    hourstring = ''
if d < 10:
    daystring = '0'
else:
    daystring = ''


fig1 = plt.figure(figsize=[14,8])
ax = fig1.add_subplot(111)

m = Basemap(projection='lcc', lat_0=(lat_min+lat_max)/2, lon_0=(lon_min+lon_max)/2,
    resolution = 'h', area_thresh = 100,
    llcrnrlon=lllon, llcrnrlat=lllat,  
    urcrnrlon=urlon, urcrnrlat=urlat)
m.drawcounties()      
m.drawmapboundary()
parallels = np.arange(lat_min.astype(int),lat_max.astype(int)+1,0.1)
meridians = np.arange(lon_min.astype(int),lon_max.astype(int)+3,0.1)          
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=tk_fs,rotation=90)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=tk_fs)

# compute map x,y coordinates of grid
x, y = m(lon[0,:,:], lat[0,:,:])

stag = 5

# plot wspd contour
CS = m.scatter(x,y,c=wspd_deficit[d,h,:,:],marker='s',s=150,cmap='magma',edgecolor='none',vmax=0,vmin=-2) 

# plot UDC wind barbs
Urot, Vrot = m.rotate_vector(u_UDC[d,h,:,:],v_UDC[d,h,:,:],lon[0],lat[0])

# convert m/s to kts
u_kts = Urot * 1.94384449
v_kts = Vrot * 1.94384449

# thin grid for plotting
u_plot = u_kts[::stag,::stag]
v_plot = v_kts[::stag,::stag]
x_plot = x[::stag,::stag]
y_plot = y[::stag,::stag]

# plot wind barbs over map
m.barbs(x_plot,y_plot,u_plot,v_plot,length=7,barbcolor='k',flagcolor='k',linewidth=0.5)

# plot wind turbines
x_farm, y_farm = m(lons_U, lats_U)
m.plot(x_farm, y_farm, '.', ms=2, c='k', label='"Upwind": Loraine' )#f7fcb9

# plot wind turbines
x_farm, y_farm = m(lons_D, lats_D)
m.plot(x_farm, y_farm, '.', ms=2, c='k', label='"Downwind": Roscoe')

# plot wind turbines
x_farm, y_farm = m(lons_C, lats_C)
m.plot(x_farm, y_farm, '.', ms=2, c='k', label='"Control": Champion')   

plt.text(0.05,0.95, '(b)', ha='center', va='center', transform=ax.transAxes,
        fontsize=lb_fs, bbox=dict(facecolor='w', edgecolor='k', alpha=0.7))

ax.add_artist(Rectangle([18830, 19050], 2.3e4, 1.9e4, ec='k', fc='none', lw=1.2))

m.drawmapscale(-100.7, 32.24, (lllon+urlon)/2., (lllat+urlat)/2., 50, barstyle='fancy', 
               yoffset=1200, units='km', fontsize=12, fillcolor1='w', fillcolor2='0.1',
               zorder=10)
    
cb = plt.colorbar(CS, format='%.1f', pad=0.05, extend='both')#, cax=cbaxes)
cb.set_label(r'Wind Speed Deficit (m s$^{-1}$)', fontsize=cb_fs, labelpad=-75)
cb.ax.tick_params(labelsize=tk_fs) 
plt.tick_params(axis='both', which='major', labelsize=tk_fs)

#plt.savefig('/Users/jeto6273/Documents/CU/Research/CNH/Final_Figs/wspd_deficit_Jan_'+daystring+str(d+1)+'_'+hourstring+str(h)+'.pdf', bbox_inches='tight')


lllon=-100.8
lllat=32.37
urlon=-100.52
urlat=32.55

fig2 = plt.figure(figsize=[14,8])
ax = fig2.add_subplot(111)

m = Basemap(projection='lcc', lat_0=(lat_min+lat_max)/2, lon_0=(lon_min+lon_max)/2,
    resolution = 'h', area_thresh = 100,
    llcrnrlon=lllon, llcrnrlat=lllat,  
    urcrnrlon=urlon, urcrnrlat=urlat)
m.drawcounties
m.drawmapboundary()

paral_values = np.arange(32.38,32.55,0.04)
merid_values = np.arange(-100.78,-100.52,0.04)

parallels = m.drawparallels(paral_values,labels=[1,0,0,0],fontsize=tk_fs+1,fmt='%.2f',linewidth=0.0,rotation=90) 
meridians = m.drawmeridians(merid_values,labels=[0,0,0,1],fontsize=tk_fs+1,fmt='%.2f',linewidth=0.0)

# for pcolormesh plotting        
lon_shift = (lon[0,:,1:] + lon[0,:,:-1])/2
lat_shift = (lat[0,1:,:] + lat[0,:-1,:])/2
x, y = m(lon_shift[:-1,:], lat_shift[:,:-1])
  
# thin grid for plotting
u_plot = u_kts
v_plot = v_kts
x_plot = 0.5*(x[:,1:] + x[:,:-1])
y_plot = 0.5*(y[1:,:] + y[:-1,:])

sum_deficit = np.nansum(((power_UDC[d,h,43,24] - power_DC[d,h,43,24]),
(power_UDC[d,h,46,25] - power_DC[d,h,46,25]),
(power_UDC[d,h,45,25] - power_DC[d,h,45,25]),
(power_UDC[d,h,44,25] - power_DC[d,h,44,25]),
(power_UDC[d,h,43,25] - power_DC[d,h,43,25]),
(power_UDC[d,h,42,25] - power_DC[d,h,42,25]),
(power_UDC[d,h,41,25] - power_DC[d,h,41,25]),

(power_UDC[d,h,45,26] - power_DC[d,h,45,26]),
(power_UDC[d,h,44,26] - power_DC[d,h,44,26]),
(power_UDC[d,h,43,26] - power_DC[d,h,43,26]),
(power_UDC[d,h,42,26] - power_DC[d,h,42,26]),
(power_UDC[d,h,41,26] - power_DC[d,h,41,26]),

(power_UDC[d,h,44,27] - power_DC[d,h,44,27]),
(power_UDC[d,h,43,27] - power_DC[d,h,43,27]),
(power_UDC[d,h,42,27] - power_DC[d,h,42,27]),
(power_UDC[d,h,41,27] - power_DC[d,h,41,27]),
(power_UDC[d,h,40,27] - power_DC[d,h,40,27]),

(power_UDC[d,h,43,28] - power_DC[d,h,43,28]),
(power_UDC[d,h,42,28] - power_DC[d,h,42,28]),
(power_UDC[d,h,41,28] - power_DC[d,h,41,28]),
(power_UDC[d,h,40,28] - power_DC[d,h,40,28]),

(power_UDC[d,h,43,29] - power_DC[d,h,43,29]),
(power_UDC[d,h,42,29] - power_DC[d,h,42,29]),
(power_UDC[d,h,41,29] - power_DC[d,h,41,29]),

(power_UDC[d,h,42,30] - power_DC[d,h,42,30]),
(power_UDC[d,h,41,30] - power_DC[d,h,41,30]),
(power_UDC[d,h,40,30] - power_DC[d,h,40,30]),

(power_UDC[d,h,42,31] - power_DC[d,h,42,31]),
(power_UDC[d,h,41,31] - power_DC[d,h,41,31]),
(power_UDC[d,h,40,31] - power_DC[d,h,40,31]),

(power_UDC[d,h,40,32] - power_DC[d,h,40,32]),
(power_UDC[d,h,40,33] - power_DC[d,h,40,33]),

(power_UDC[d,h,39,32] - power_DC[d,h,39,32]),
(power_UDC[d,h,39,33] - power_DC[d,h,39,33]),
(power_UDC[d,h,39,34] - power_DC[d,h,39,34]),
(power_UDC[d,h,39,35] - power_DC[d,h,39,35]),

(power_UDC[d,h,38,32] - power_DC[d,h,38,32]),
(power_UDC[d,h,38,33] - power_DC[d,h,38,33]),
(power_UDC[d,h,38,34] - power_DC[d,h,38,34]),
(power_UDC[d,h,38,35] - power_DC[d,h,38,35]),
(power_UDC[d,h,38,36] - power_DC[d,h,38,36]),
(power_UDC[d,h,38,37] - power_DC[d,h,38,37]),
(power_UDC[d,h,38,38] - power_DC[d,h,38,38]),
(power_UDC[d,h,38,39] - power_DC[d,h,38,39]),
(power_UDC[d,h,38,40] - power_DC[d,h,38,40]),
(power_UDC[d,h,38,41] - power_DC[d,h,38,41]),        

(power_UDC[d,h,37,31] - power_DC[d,h,37,31]),
(power_UDC[d,h,37,32] - power_DC[d,h,37,32]),
(power_UDC[d,h,37,33] - power_DC[d,h,37,33]),
(power_UDC[d,h,37,34] - power_DC[d,h,37,34]),
(power_UDC[d,h,37,35] - power_DC[d,h,37,35]),
(power_UDC[d,h,37,36] - power_DC[d,h,37,36]),
(power_UDC[d,h,37,37] - power_DC[d,h,37,37]),
(power_UDC[d,h,37,38] - power_DC[d,h,37,38]),
(power_UDC[d,h,37,39] - power_DC[d,h,37,39]),
(power_UDC[d,h,37,40] - power_DC[d,h,37,40]),
(power_UDC[d,h,37,41] - power_DC[d,h,37,41]),

(power_UDC[d,h,36,32] - power_DC[d,h,36,32]),
(power_UDC[d,h,36,33] - power_DC[d,h,36,33]),
(power_UDC[d,h,36,34] - power_DC[d,h,36,34]),
(power_UDC[d,h,36,35] - power_DC[d,h,36,35]),
(power_UDC[d,h,36,36] - power_DC[d,h,36,36]),
(power_UDC[d,h,36,37] - power_DC[d,h,36,37]),
(power_UDC[d,h,36,38] - power_DC[d,h,36,38]),
(power_UDC[d,h,36,39] - power_DC[d,h,36,39]),
(power_UDC[d,h,36,40] - power_DC[d,h,36,40]),

(power_UDC[d,h,35,32] - power_DC[d,h,35,32]),
(power_UDC[d,h,35,33] - power_DC[d,h,35,33]),
(power_UDC[d,h,35,34] - power_DC[d,h,35,34]),
(power_UDC[d,h,35,35] - power_DC[d,h,35,35]),
(power_UDC[d,h,35,36] - power_DC[d,h,35,36]),

(power_UDC[d,h,36,31]*4/5. - power_DC[d,h,36,31]),
(power_UDC[d,h,35,30]*1/3. - power_DC[d,h,35,30]),
(power_UDC[d,h,38,31]*1/3. - power_DC[d,h,38,31]),
(power_UDC[d,h,39,30]*1/4. - power_DC[d,h,39,30]),
(power_UDC[d,h,39,31]*3/4. - power_DC[d,h,39,31])))

pd1 = power_UDC[d,h,36,31]*4/5. - power_DC[d,h,36,31] 
pd2 = power_UDC[d,h,35,30]*1/3. - power_DC[d,h,35,30]
pd3 = power_UDC[d,h,38,31]*1/3. - power_DC[d,h,38,31]
pd4 = power_UDC[d,h,39,30]*1/4. - power_DC[d,h,39,30]
pd5 = power_UDC[d,h,39,31]*3/4. - power_DC[d,h,39,31] 

power_deficit[d,h,36,31] = pd1
power_deficit[d,h,35,30] = pd2
power_deficit[d,h,38,31] = pd3
power_deficit[d,h,39,30] = pd4
power_deficit[d,h,39,31] = pd5

# plot p contour
CS = m.pcolormesh(x,y,power_deficit[d,h,1:,1:]/1e6, vmax=0, vmin=-2, cmap='magma',
        edgecolors='k',linestyle=':',lw=0.5) 

# plot wind turbines
x_farm, y_farm = m(lons_U, lats_U)
m.plot(x_farm, y_farm, '.', ms=5, c='k', label='"Upwind": Loraine') #navy

# plot wind turbines
x_farm, y_farm = m(lons_D, lats_D)
m.plot(x_farm, y_farm, '^', ms=4, c='k', label='"Downwind": Roscoe') #a50f15

# plot wind turbines
x_farm, y_farm = m(lons_C, lats_C)
m.plot(x_farm, y_farm, '.', ms=5, c='k', label='"Control": Champion')   

plt.text(19183,18800, 'Total power deficit = ' + str(np.round(sum_deficit/1e6,1)) + ' MW', ha='center', va='center',
        fontsize=bx_fs, bbox=dict(facecolor='w', edgecolor='k', alpha=0.7))

plt.text(0.05,0.95, '(a)', ha='center', va='center', transform=ax.transAxes,
        fontsize=26, bbox=dict(facecolor='w', edgecolor='k', alpha=0.7))

m.drawmapscale(-100.739, 32.382, (lllon+urlon)/2., (lllat+urlat)/2., 10, barstyle='fancy', 
       yoffset=350, units='km', fontsize=12, fillcolor1='w', fillcolor2='0.1',
       zorder=10)
       
cb = plt.colorbar(CS, format='%.1f', pad=0.05, extend='both')#, cax=cbaxes)
cb.set_label(r'Power Deficit (MW)', fontsize=cb_fs+1, labelpad=-78)
cb.ax.tick_params(labelsize=tk_fs+1) 
plt.tick_params(axis='both', which='major', labelsize=tk_fs)
 
#plt.savefig(‘power_deficit_Jan_'+daystring+str(d+1)+'_'+hourstring+str(h)+'_nobarbs.pdf', bbox_inches='tight')
#################################