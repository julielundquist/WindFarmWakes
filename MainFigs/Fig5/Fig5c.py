#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 19:18:07 2017

@author: Jessica Tomaszewski
"""

import numpy as np
import matplotlib.pyplot as plt
import cmocean
import matplotlib.colors as mcolors
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
hourly_hfx = nc_fid.variables['hourly_hfx'][:] 
#################################


#################################
# Make figures

# font sizes
cb_fs = 18
lb_fs = 18
tk_fs = 14
bx_fs = 20

# circular colormap
old1 = cmocean.cm.turbid
old2 = cmocean.cm.speed
indices = np.linspace(0,1,61)
colors1 = old1(indices[4:])
colors2 = old2(indices[::-1][:])
colors = np.vstack((colors1, colors2))
mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

mat = np.vstack([hourly_total_power_deficit, hourly_hfx, hourly_wd, hourly_ws])
sort = np.argsort(mat[0])
sorted_power = mat[0,sort]
sorted_HFX = mat[1,sort]
sorted_wdir = mat[2,sort]
sorted_wspd = mat[3,sort]

fig1 = plt.figure(figsize=[18,8])

ax = plt.subplot(111)
s = ax.scatter(sorted_HFX, sorted_wspd, c=sorted_wdir, vmin=0, vmax=360,
               s=sorted_power*30, cmap=mymap, alpha=0.6) #cmocean.cm.delta

ax.plot(np.zeros(100), np.linspace(-10,80,100), 'k-', lw=2, zorder=0)
ax.plot(np.linspace(-100,350,100), np.ones(100)*12, '-', c='darkred', lw=2, zorder=0)

a = plt.cm.viridis
indices = np.linspace(0,1,61)
cmap_indices = a(indices)

ax.scatter(-200, -200, facecolor='none', s=50*30, label='50', alpha=0.8)
ax.scatter(-200, -200, facecolor='none', s=40*30, label='40', alpha=0.8)
ax.scatter(-200, -200, facecolor='none', s=30*30, label='30', alpha=0.8)
ax.scatter(-200, -200, facecolor='none', s=20*30, label='20', alpha=0.8)
ax.scatter(-200, -200, facecolor='none', s=10*30, label='10', alpha=0.8)

plt.text(242, 15.32, 'Hourly Total Power Deficit (MW)', ha='center', va='center', 
         rotation=90, fontsize=tk_fs, bbox=dict(facecolor='w', edgecolor='none', alpha=0.8))

plt.text(345, 0, 'N', ha='center', va='center', fontweight='medium',
         fontsize=tk_fs+4, bbox=dict(facecolor='none', edgecolor='none'))
plt.text(345, 10, 'S', ha='center', va='center', fontweight='medium', 
         fontsize=tk_fs+4, bbox=dict(facecolor='none', edgecolor='none'))
plt.text(345, 5, 'E', ha='center', va='center', fontweight='medium',
         fontsize=tk_fs+4, bbox=dict(facecolor='none', edgecolor='none'))
plt.text(345, 15, 'W', ha='center', va='center', fontweight='medium',
         fontsize=tk_fs+4, bbox=dict(facecolor='none', edgecolor='none'))
plt.text(345, 20, 'N', ha='center', va='center', fontweight='medium',
         fontsize=tk_fs+4, bbox=dict(facecolor='none', edgecolor='none'))

cb = plt.colorbar(s, pad=0.035)
cb.ax.tick_params(labelsize=tk_fs)
cb.set_label(r'Wind Direction ($\circ$)', fontsize=cb_fs, labelpad=-78)
cb.set_ticks(np.arange(0,370,45)) 
cb.set_alpha(1)
cb.draw_all()

ax.set_xlim(-75,300)
ax.set_ylim(0,20)
ax.set_yticks(np.arange(0,22,2))

plt.xlabel(r'Sensible Heat Flux (W m$^{-2}$)', fontsize=lb_fs)
plt.ylabel(r'Wind Speed (m s$^{-1}$)', fontsize=lb_fs)
plt.tick_params(axis='both', which='major', labelsize=tk_fs)
ax.grid(True)

leg = plt.legend(loc='upper right', scatterpoints=1, borderpad=1.3, handletextpad=1.3,
           fontsize=13.5, labelspacing=2, bbox_to_anchor=(0.96, 0.985)),#, framealpha=0.8,

plt.text(0.022, 0.96, '(c)', ha='center', va='center', transform=ax.transAxes, fontsize=lb_fs+8)#,

plt.tight_layout()
#plt.savefig('HFX-vs-winds_power.pdf', bbox_inches='tight') #, dpi=1200)
#################################