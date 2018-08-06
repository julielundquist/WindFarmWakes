#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 15:16:58 2017

@author: Jessica Tomaszewski
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4


#################################
# Read NETCDF of simulation data
path = ''
f = 'fig3_power_timeseries.nc'
nc_fid = nc4.Dataset(path+f, 'r')  
#################################


#################################
# Extract variables
u = nc_fid.variables['hourly_u'][:]
v = nc_fid.variables['hourly_v'][:]
hourly_total_power_downwind = nc_fid.variables['hourly_total_power_downwind'][:]
hourly_total_power_downwind_waked = nc_fid.variables['hourly_total_power_downwind_waked'][:]
hourly_total_power_upwind = nc_fid.variables['hourly_total_power_upwind'][:]
#################################


#################################
# Read turbine data
datafile = open('TX_turbines_UDC.csv', 'r')
data = datafile.readlines()[1:]
data = [x.split(',') for x in data]

lats_farm = np.array([x[12] for x in data])
lons_farm = np.array([x[13] for x in data])
#################################


#################################
# Make figure

# font sizes
cb_fs = 18
lb_fs = 22
tk_fs = 14
bx_fs = 20

# prep plot bottom panel
lw = 0.5
stag = 6
time = np.arange(744)/24.+1
barb_line = np.ones(len(time))*-50
bottom_line = np.ones(len(time))*-95


fig = plt.figure(figsize=[17.5,8.25])
ax = plt.subplot(111)

# plot power time series
plt.fill_between(time, 0, hourly_total_power_downwind, alpha=1,
                 facecolor='#c57c3c',lw=lw,label='Only Downwind') 
plt.fill_between(time, 0, hourly_total_power_downwind_waked, alpha=0.9,
                 facecolor='#66c2a5',lw=lw,label='Downwind waked\n by Upwind')  
plt.fill_between(time, 0, hourly_total_power_upwind, alpha=0.8,
                 facecolor='#7570b3',lw=lw,label='Only Upwind')


plt.fill_between(time, 0, bottom_line, alpha=0.3, facecolor='gray',lw=lw)
bb = plt.barbs(time[::stag],barb_line[::stag],u[::stag]*1.94,v[::stag]*1.94,length=8.25)
bb.set_clip_box(None) 
plt.plot(time,barb_line,'k')

plt.ylim(-95,350)
plt.xlim(1,31.9)
plt.xticks(np.arange(min(time), max(time), 1))
ax.set_yticklabels(['','',0,50,100,150,200,250,300,350])

plt.tick_params(axis='both', which='major', labelsize=tk_fs)
plt.xlabel('Time (UTC day of month, Jan)', size=lb_fs)
plt.ylabel(r'Total farm power and winds (MW, kts)', size=lb_fs, labelpad=15)
plt.grid()
plt.legend(loc='upper left')
    
plt.subplots_adjust(left=0.065, bottom=None, right=0.98, top=None, wspace=None, hspace=None)
#################################