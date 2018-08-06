#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 15:36:03 2017

@author: jeto6273
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4


#################################
# Read NETCDF of simulation data
path = ''
f = 'fig5_total_power_deficits.nc'
nc_fid = nc4.Dataset(path+f, 'r')  
#################################


#################################
# Extract variables
hourly_total_power_deficit = nc_fid.variables['hourly_total_power_deficit'][:]
hourly_wd = nc_fid.variables['hourly_wd'][:]
hourly_ws = nc_fid.variables['hourly_ws'][:]
#################################


#################################
# Make figures

# font sizes
cb_fs = 18
lb_fs = 18
tk_fs = 14
bx_fs = 20

mat = np.vstack([hourly_total_power_deficit, hourly_wd*(np.pi/180), hourly_ws])
sort = np.argsort(mat[0])
sorted_power = mat[0,sort]
sorted_wdir = mat[1,sort]
sorted_wspd = mat[2,sort]

fig1 = plt.figure(figsize=[9,8])

ax = plt.subplot(111, projection='polar')
s = ax.scatter(sorted_wdir, sorted_wspd, c=sorted_power, vmin=0, vmax=60,
               s=sorted_power*30, cmap='viridis', alpha=0.6)

ax.plot(np.arange(360)*(np.pi/180), np.ones(360)*12, c='darkred', lw=2, zorder=0)
ax.set_theta_offset(np.pi/2)
ax.set_theta_direction(-1)
ax.set_rmax(19.7)
ax.set_rlabel_position(113)  # get radial labels away from plotted line, 20
plt.text(-0.035, 1.01, '(b)', ha='center', va='center', transform=ax.transAxes, fontsize=lb_fs+8)#,
         #bbox=dict(facecolor='w', edgecolor='k', alpha=0.7))
x0 = 0.13
y0 = 0.37
ax.annotate('', xy=(4.206,13), xytext=(4.206,22), #textcoords='figure fraction',
            arrowprops=dict(facecolor='k', width=2))
plt.tick_params(axis='x', which='major', labelsize=tk_fs+2)

#cbaxes = plt.axes([0.1, -0.05, 0.83, 0.04]) # left, bottom, width, height
cbaxes = plt.axes([1, 0.1, 0.04, 0.83]) # left, bottom, width, height
cb = plt.colorbar(s, format='%.1f', orientation='vertical', cax=cbaxes)
cb.ax.tick_params(labelsize=tk_fs)
cb.set_alpha(1)
cb.draw_all()

ax.set_yticks(range(0, 20, 2))
plt.ylabel(r'Hourly total farm instantaneous power deficit (MW)', fontsize=lb_fs, labelpad=-90)
ax.grid(True)
plt.text(1.01, 0.315, 'm s$^{-1}$', ha='center', va='center', transform=ax.transAxes,
                fontsize=12, bbox=dict(facecolor='w', edgecolor='w', alpha=0))

#plt.savefig(‘Jan_power-deficit_v2.pdf', bbox_inches='tight')#, dpi=1200)
#################################