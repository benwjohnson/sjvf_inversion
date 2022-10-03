#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Driver script used in Johnson et al., 2022 to calculate alteration fluid d18O values in the San Juan Mountains CO

@author: bwj
"""
# import sys
# sys.path.insert(0, '/Users/bwj/Google_Drive/python_scripts/')
# sys.path.insert(0, '/Users/bwj/Google_Drive/current_projects/geochemical_inverse_modeling')
# sys.path.insert(0, '/Users/bwj/Google_Drive/current_projects/geochemical_inverse_modeling/meteoric_volcanic/san_juan_volcanic_field')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


plt.close('all')

#%% Calling individual functions for each location 
from plot_lc import plot_lc
from plot_bonanza import plot_bonanza
from plot_silverton import plot_silverton
from plot_stony import plot_stony

[LC_water,LC_X,LC_Z,del18O_LC,temp_LC,
 LC_x_grid,LC_z_grid,LC_iso_grid,temp_grid_LC,LC_W_R,
 x_center_points,y_center_points,summed_horz_vectors,summed_vert_vectors] = plot_lc()

[bon_water,bon_X,bon_Z,del18O_bon,temp_bon,
  bonz_x_grid,bonz_z_grid,bonz_iso_grid,temp_grid_bonz,bonz_W_R,
  bx_center_points,by_center_points,bsummed_horz_vectors,bsummed_vert_vectors] = plot_bonanza()

[silverton_water,silverton_X,silverton_Z,del18O_silverton,temp_silverton,
  silverton_x_grid,silverton_z_grid,silverton_iso_grid,temp_grid_silverton,silverton_W_R,
  six_center_points,siy_center_points,sisummed_horz_vectors,sisummed_vert_vectors] = plot_silverton()

[stony_water,stony_X,stony_Z,del18O_stony,temp_stony,
  SM_x_grid,SM_z_grid,SM_iso_grid,temp_grid_SM,SM_W_R,
  stx_center_points,sty_center_points,stsummed_horz_vectors,stsummed_vert_vectors] = plot_stony()

#%% Package into dataframe
calderaNames = ['Bonanza','Lake City','Silverton','Stony Mountain']
initial_water = [bon_water,LC_water,silverton_water,stony_water]
max_ages = [33,30,20,25]
min_ages = [32,28,18,14]
caldera_info = pd.DataFrame([initial_water,max_ages,min_ages],columns=calderaNames)

#%% Plots
## In manuscript plots

#age ranges, based on Lipman and Bachmann, 2015. Plotting total range of activity
#for each volcanic center, since we can't at the moment say when hydrothermal activity
#was at its peak



plt.close('all')
iso_colormap = plt.cm.cividis
temp_colormap = plt.cm.magma

bonz_O_contours = [5,6,7,8,9]
LC_O_contours = [-1,1,3,5,7]
silv_O_contours = [-5,-3,-1,1,3,5,7]
SM_O_contours = LC_O_contours

bonz_temp_contours = [90,95,100,105,110]
LC_temp_contours = [60,90,120,150,180,210]
silv_temp_contours = [100,120,150,160,180,200,220]
SM_temp_contours = [20,60,100,140,200]

fig2, (ax1,ax2) = plt.subplots(1,2,figsize=(12, 5),num=2); plt.clf() 
plt.subplot(1,2,1)   
cplot=plt.contourf(LC_x_grid/1e5,LC_z_grid,LC_iso_grid,LC_O_contours,extend ='both',cmap=plt.cm.cividis)
cbar=plt.colorbar(fraction=0.046, pad=0.04)
cbar.ax.set_ylabel('$\delta^{18}$O (‰ V-SMOW)',rotation=270)
plt.scatter(LC_X/1e5,LC_Z,s=60,marker='^',c='w',edgecolors='k')
plt.quiver(x_center_points/1e5,y_center_points,summed_horz_vectors,-summed_vert_vectors,
           units = 'x',color='w',edgecolor='k',width=.03,linewidth = .5)
plt.xlabel('Horizontal (km)');plt.ylabel('Vertical (m)');
plt.xlim([LC_x_grid.min()/1e5,LC_x_grid.max()/1e5]);plt.ylim([LC_z_grid.min(),LC_z_grid.max()]);
# cm = plt.cm.get_cmap('RdYlBu')

plt.title('Lake City')
plt.subplot(1,2,2)   
plt.scatter(LC_X/1e5,LC_Z,s=60,marker='^',c='w',edgecolors='k')

ax = plt.gca()
plt.xlabel('Horizontal (km)');plt.ylabel(''); ax.set_yticklabels([''])
plt.xlim([LC_x_grid.min()/1e5,LC_x_grid.max()/1e5]);plt.ylim([LC_z_grid.min(),LC_z_grid.max()]);



#%%
import scipy.ndimage
from scipy.interpolate import griddata

numsteps= 15 #default = 15
z=LC_Z;x=LC_X
del18O=del18O_LC
Zmax = max(z); Zmin=min(z)
Xmax = max(x); Xmin=min(x)

xi = np.arange(Xmin,Xmax+1,(Xmax-Xmin)/numsteps)
zi = np.arange(Zmin,Zmax+1,(Zmax-Zmin)/numsteps)
xi,zi = np.meshgrid(xi,zi)  

 
iso_grid = griddata((x,z),del18O,(xi,zi))
plt.figure()
plt.contourf(xi/1e5,zi,iso_grid)

fig2new, (ax1,ax2) = plt.subplots(1,2,figsize=(12, 5),num=2); plt.clf() 
plt.subplot(1,2,1)   
cplot=plt.contourf(xi/1e5,zi,iso_grid,LC_O_contours,extend ='both',cmap=plt.cm.cividis)
cbar=plt.colorbar(fraction=0.046, pad=0.04)
cbar.ax.set_ylabel('$\delta^{18}$O (‰ V-SMOW)',rotation=270)
plt.scatter(LC_X/1e5,LC_Z,s=60,marker='^',c='w',edgecolors='k')
plt.xlabel('Horizontal (km)');plt.ylabel('Vertical (m)');
plt.xlim([LC_x_grid.min()/1e5,LC_x_grid.max()/1e5]);plt.ylim([LC_z_grid.min(),LC_z_grid.max()]);


plt.title('Lake City')
plt.subplot(1,2,2)   
plt.scatter(LC_X/1e5,LC_Z,s=60,marker='^',c='w',edgecolors='k')
ax = plt.gca()
plt.xlabel('Horizontal (km)');plt.ylabel(''); ax.set_yticklabels([''])
plt.xlim([LC_x_grid.min()/1e5,LC_x_grid.max()/1e5]);plt.ylim([LC_z_grid.min(),LC_z_grid.max()]);

#%% Supplemental plots 


sfig1 = plt.figure();plt.clf()
x_grid_list=[bonz_x_grid/1e5,LC_x_grid/1e5,silverton_x_grid,SM_x_grid/1e4]
z_grid_list = [bonz_z_grid,LC_z_grid,silverton_z_grid,SM_z_grid]
iso_grid_list = [bonz_iso_grid,LC_iso_grid,silverton_iso_grid,SM_iso_grid]
temp_grid_list = [temp_grid_bonz,temp_grid_LC,temp_grid_silverton,temp_grid_SM]

O_contour_list = [bonz_O_contours,LC_O_contours,silv_O_contours,SM_O_contours]
temp_contour_list = [bonz_temp_contours,LC_temp_contours,silv_temp_contours,SM_temp_contours]
iso_plot_id =-1
temp_plot_id = 0
title_list = ['Bonanza','Lake City','Silverton','Stony Mountain']
label_list = ['a','b','c','d']
for iplot in range(0,len(x_grid_list)):
    iso_plot_id = iso_plot_id+2
    temp_plot_id =temp_plot_id+2

    plt.subplot(4,2,iso_plot_id)
    plt.contourf(x_grid_list[iplot],z_grid_list[iplot],
                  iso_grid_list[iplot],O_contour_list[iplot],cmap=iso_colormap)
    plt.xlim([x_grid_list[iplot].min(),x_grid_list[iplot].max()])
    plt.ylim([z_grid_list[iplot].min(),z_grid_list[iplot].max()])
    plt.ylabel('Height (m)')
    plt.title(title_list[iplot],x=1.4,y=1)

  
    if iso_plot_id == 7:
         plt.xlabel('Horizontal (km)')

    
    cbar=plt.colorbar()
    cbar.ax.set_ylabel('$\delta^{18}$O',rotation=270)
    cbar.ax.get_yaxis().labelpad = 15
    
    plt.subplot(4,2,temp_plot_id)
    plt.contourf(x_grid_list[iplot],z_grid_list[iplot],
                  temp_grid_list[iplot],temp_contour_list[iplot],cmap=temp_colormap)
    cbar=plt.colorbar()
    ax=plt.gca()
    ax.set_yticklabels([])
    cbar.ax.set_ylabel('Temp ($^\circ$C)',rotation=270)
    cbar.ax.get_yaxis().labelpad = 15
    
    if temp_plot_id == 8:
        plt.xlabel('Horizontal (km)')
                   
sfig1.tight_layout()

sfig2  = plt.figure();plt.clf() 
plt.subplot(2,2,1)
cplot=plt.contourf(bonz_x_grid/1e5,bonz_z_grid,bonz_iso_grid,bonz_O_contours,extend ='both',cmap=plt.cm.cividis)
cbar=plt.colorbar()
cbar.ax.set_ylabel('$\delta^{18}$O (‰ V-SMOW)',rotation=270)
cbar.ax.get_yaxis().labelpad = 15
plt.quiver(bx_center_points/1e5,by_center_points,bsummed_horz_vectors,-bsummed_vert_vectors,units = 'y',color='w')
plt.ylabel('Vertical (m)'); ax.set_yticklabels([''])
plt.xlim([bonz_x_grid.min()/1e5,bonz_x_grid.max()/1e5]);plt.ylim([bonz_z_grid.min(),bonz_z_grid.max()]);
plt.title('Bonanza')

    
 
plt.subplot(2,2,2)
cplot=plt.contourf(LC_x_grid/1e5,LC_z_grid,LC_iso_grid,LC_O_contours,extend ='both',cmap=plt.cm.cividis)
cbar=plt.colorbar()
cbar.ax.set_ylabel('$\delta^{18}$O (‰ V-SMOW)',rotation=270)
cbar.ax.get_yaxis().labelpad = 15
plt.quiver(x_center_points/1e5,y_center_points,summed_horz_vectors,-summed_vert_vectors,units = 'y',color='w')
plt.ylabel(''); ax.set_yticklabels([''])
plt.xlim([LC_x_grid.min()/1e5,LC_x_grid.max()/1e5]);plt.ylim([LC_z_grid.min(),LC_z_grid.max()]);
plt.title('Lake City')
 
plt.subplot(2,2,3)
cplot=plt.contourf(silverton_x_grid/1e5,silverton_z_grid,silverton_iso_grid,silv_O_contours,extend ='both',cmap=plt.cm.cividis)
ax=plt.gca()
ax.set_yticklabels([])
cbar=plt.colorbar()
cbar.ax.set_ylabel('$\delta^{18}$O (‰ V-SMOW)',rotation=270)
cbar.ax.get_yaxis().labelpad = 15
plt.quiver(six_center_points/1e5,siy_center_points,sisummed_horz_vectors,-sisummed_vert_vectors,units = 'y',color='w')
plt.ylabel('Vertical (m)'); plt.xlabel('Horizontal (km)');
plt.xlim([silverton_x_grid.min()/1e5,silverton_x_grid.max()/1e5]);plt.ylim([silverton_z_grid.min(),silverton_z_grid.max()]);
plt.title('Silverton')
 
plt.subplot(2,2,4)
cplot=plt.contourf(SM_x_grid/1e4,SM_z_grid,SM_iso_grid,SM_O_contours,extend ='both',cmap=plt.cm.cividis)
ax=plt.gca()
ax.set_yticklabels([])
cbar=plt.colorbar()
cbar.ax.set_ylabel('$\delta^{18}$O (‰ V-SMOW)',rotation=270)
cbar.ax.get_yaxis().labelpad = 15
plt.quiver(stx_center_points/1e4,sty_center_points,stsummed_horz_vectors,-stsummed_vert_vectors,units = 'y',color='w')
plt.xlabel('Horizontal (km)'); ax.set_yticklabels([''])
plt.xlim([SM_x_grid.min()/1e4,SM_x_grid.max()/1e4]);plt.ylim([SM_z_grid.min(),SM_z_grid.max()]);
plt.title('Stony Mountain')
 

#%% Output table with isotope data, location information, elevation, and temperature used in inversion

col_names = ['delta18O','Temperature (deg. C)','X-coordinate','Z-coordinate']

import math

full_d18O = del18O_bon+del18O_LC+del18O_silverton.tolist()+del18O_stony
full_temp = temp_bon.tolist()+temp_LC.tolist()+temp_silverton.tolist()+temp_stony.tolist()
full_X = bon_X+ LC_X.tolist() + silverton_X.tolist() + stony_X.tolist()
full_Z = bon_Z + LC_Z.tolist() +silverton_Z.tolist() + stony_Z.tolist()

table_d18O = np.zeros(len(full_d18O))
table_temp = np.zeros(len(full_d18O))
table_X = np.zeros(len(full_d18O))
table_Z = np.zeros(len(full_d18O))

for inum in range(0,len(full_d18O)):
    table_d18O[inum] = round(full_d18O[inum],2)
    table_temp[inum] = int(round(full_temp[inum]))
    table_X[inum] = int(round(full_X[inum]))
    table_Z[inum] = round(full_Z[inum])
    table_Z[inum]= int(table_Z[inum])
    



full_data = pd.DataFrame({col_names[0]: table_d18O,
                          col_names[1]: table_temp,
                          col_names[2]: table_X,
                          col_names[3]: table_Z})


output_supp=full_data.to_latex(index=False,float_format=lambda x: '%.0f' % x)

si_saved=full_data.to_csv('input_data.csv')