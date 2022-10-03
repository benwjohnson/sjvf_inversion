#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 16:41:43 2019

@author: benjohnson
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 16:06:55 2019

@author: benjohnson
"""

#####-------Bonanza Caldera Inversion ------###
# Based on data from Varga and Smith, 1984. 
def plot_bonanza():
    # import sys
    # sys.path.insert(0, '/Users/bwj/Google_Drive/python_scripts/')
    # sys.path.insert(0, '/Users/bwj/Google_Drive/current_projects/geochemical_inverse_modeling/')
    
    import pandas as pd 
    import numpy as np
    import matplotlib.pyplot as plt
    
    from o_isotope_invert import o_isotope_invert
    plt.close('all')
    #%% ----create regular dataset ------------- %%# 
    data_raw = pd.read_csv('Bonanza_data.csv')
    
    del18O = pd.Series.tolist(data_raw['d18O'])
    x = pd.Series.tolist(data_raw['X'])
    z = pd.Series.tolist(data_raw['Z'])
    # temp = pd.Series.tolist(data_raw['Temp'])
    max_temp = 100
    min_temp = 50
    temp_slope = np.divide((max_temp-min_temp),(np.max(del18O)- np.min(del18O)))
    temp = np.multiply(temp_slope,del18O) +min_temp
    
    #%% Inversion
    d18Oinit_bonz=7.8
    numsamples = len(x)
    bonz_water_initial = np.zeros(numsamples)
    bonz_water_outgoing = np.zeros(numsamples)
    bonz_moles_fluid = np.zeros(numsamples)
    bonz_W_R = np.zeros(numsamples)
    bonz_final_rock = np.zeros(numsamples)
    bonz_moles_Orock = np.zeros(numsamples)
    for irun in range(0,numsamples):
        del18O_in_bonz = del18O
        Temp_in_bonz = temp
        x_bonz_in = x
        z_bonz_in = z
        
        num_takeout = np.random.randint(0,len(del18O_in_bonz))
        
        del18O_in_bonz = np.delete(del18O_in_bonz,num_takeout)
        Temp_in_bonz =np.delete(Temp_in_bonz,num_takeout)
        x_bonz_in =  np.delete(x,num_takeout)
        z_bonz_in =np.delete(z,num_takeout)
        
        [bonz_iso_grid,bonz_x_grid,bonz_z_grid,conc_matrix_bonz,change_matrix_bonz,qsolved_bonz,\
         moles_fluid_bonz,W_R_taper,water_initial_bonz, water_outgoing_bonz, temp_test_bonz, temp_grid_bonz,moles_Orock_bonz,bonz_epsilon_w_r,\
         x_center_points,y_center_points,summed_horz_vectors,summed_vert_vectors] = \
        o_isotope_invert(del18O_in_bonz,Temp_in_bonz,x_bonz_in,z_bonz_in,d18Oinit_bonz)
        bonz_water_initial[irun] = water_initial_bonz
        bonz_water_outgoing[irun] = water_outgoing_bonz
        bonz_moles_fluid[irun] = moles_fluid_bonz
        bonz_W_R[irun] = W_R_taper
        bonz_final_rock[irun] = del18O_in_bonz.mean()
        bonz_moles_Orock[irun] = moles_Orock_bonz
    #%% Plots
    # fig1=plt.figure(num=1); plt.clf()    
    # cplot=plt.contourf(bonz_x_grid,bonz_z_grid,bonz_iso_grid)
    # plt.plot(x,z,'w^')
    # plt.xlabel('Horizontal (m)');plt.ylabel('Vertical (m)');
    # cm = plt.cm.get_cmap('RdYlBu')
    # cbar=plt.colorbar()
    # cbar.ax.set_ylabel('$\delta^{18}$O',rotation=270)
    
    # fig2 = plt.figure(num=2); plt.clf()
    # plt.hist(bonz_water_initial)
    
    return [bonz_water_initial,x,z,del18O,temp,bonz_x_grid,bonz_z_grid,
            bonz_iso_grid,temp_grid_bonz,bonz_W_R,
            x_center_points,y_center_points,summed_horz_vectors,summed_vert_vectors]