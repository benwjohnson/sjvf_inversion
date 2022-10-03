#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 16:06:55 2019

@author: benjohnson
"""

#####-------Silverton Caldera Inversion ------###
# Based on data from Eberl et al., 1987. Their samples RM and LF formed 
# in a 21 Ma event, and I go forward here assuming these represent a crustal
# section that is 3500 m x 300 m. They measured a number of mineral specific
#d18O values, mostly in sericite, which according to Elsinger and Savin, 1973
#are usually about 3permil depleted compared to the whole rock 
def plot_silverton():
    # import sys
    # sys.path.insert(0, '/Users/bwj/Google_Drive/python_scripts/')
    # sys.path.insert(0, '/Users/bwj/Google_Drive/current_projects/geochemical_inverse_modeling/')
    
    import pandas as pd 
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import random 
    from o_isotope_invert import o_isotope_invert
    plt.close('all')
    #%% ----create regular dataset ------------- %%# 
    
    #spatial grid, in m
    Zmax = 300; Zmin=0 #pile of crust 7 km thick
    Xmax = 3500; Xmin=0 #and 10 km wide
    Omax = 5.4
    Omin =  3.3
    tempmin = 100
    tempmax = 200
    
    numsteps = 50 #50x50 grid
    
    xi = np.arange(Xmin,Xmax+1,(Xmax-Xmin)/numsteps)
    zi = np.arange(Zmin,Zmax+1,(Zmax-Zmin)/numsteps)
    xgrid,zgrid = np.meshgrid(xi,zi)  
    
    ##create regular d18O and temp pattern
    del18O = np.zeros([numsteps+1,numsteps+1])
    numrows,numcols = del18O.shape
    temperature =  np.zeros([numsteps+1,numsteps+1])
    
    #interpolate
    for icol in range(0,numcols):
        del18O[:,icol] = np.linspace(Omin,Omax,numrows)
        temperature[:,icol] = np.flip(np.linspace(tempmin,tempmax,numrows),0)
        
    numsamples = 125
    d18Oinit =7.8
    
    rand_sample_row = np.random.randint(0,high=numrows,size=numsamples)
    rand_sample_col = np.random.randint(0,high=numcols,size=numsamples)
    
    randd18O = del18O[rand_sample_row,rand_sample_col]
    randtemp = temperature[rand_sample_row,rand_sample_col]
    randx = xgrid[rand_sample_row,rand_sample_col]
    randz = zgrid[rand_sample_row,rand_sample_col]
    randxi, randzi = np.meshgrid(randx,randz)
    
    #%% Read in data 
    silv_data = pd.read_csv('silverton_data.csv')
    #Estimate WR values from quartz and sericite
    #WR is usually 3permil higher than sericite, and qtz ~1 permil higher than WR
    d18O_raw =  silv_data['d18O']
    wr_d18O = np.zeros(len(silv_data['Sample']))
    ser_idx = silv_data['phase']=='sericite'        
    q_idx = silv_data['phase']=='q'  
       
    for isamp in range(0,len(wr_d18O)):
        if ser_idx[isamp]==True:
            wr_d18O[isamp] = d18O_raw[isamp]+3
        elif q_idx[isamp]==True:
            wr_d18O[isamp] = d18O_raw[isamp] -3
        else:
            wr_d18O[isamp] = d18O_raw[isamp]
            
    silv_x = np.asarray(silv_data['UTM E'])
    silv_z = np.asarray(silv_data['Elevation'])
    
    max_temp = 300
    min_temp = 200
    temp_slope = np.divide((max_temp-min_temp),(np.max(wr_d18O)- np.min(wr_d18O)))
    temp = np.multiply(temp_slope,wr_d18O) +min_temp
        
     
    
    d18Oinit_silv=7.8
    numsamples=len(d18O_raw)
    silv_water_initial = np.zeros(numsamples)
    silv_water_outgoing = np.zeros(numsamples)
    silv_moles_fluid = np.zeros(numsamples)
    silv_W_R = np.zeros(numsamples)
    silv_final_rock = np.zeros(numsamples)
    silv_moles_Orock = np.zeros(numsamples)
    for irun in range(0,numsamples):
        del18O_in_silv = wr_d18O
        Temp_in_silv = temp
        x_silv_in = silv_x
        z_silv_in = silv_z
        
        num_takeout = np.random.randint(0,len(del18O_in_silv))
        
        del18O_in_silv = np.delete(del18O_in_silv,num_takeout)
        Temp_in_silv =np.delete(Temp_in_silv,num_takeout)
        x_silv_in =  np.delete(x_silv_in,num_takeout)
        z_silv_in =np.delete(z_silv_in,num_takeout)
        
        [silv_iso_grid,silv_x_grid,silv_z_grid,conc_matrix_silv,change_matrix_silv,qsolved_silv,\
         moles_fluid_silv,W_R_taper,water_initial_silv, water_outgoing_silv, temp_test_silv, temp_grid_silv,moles_Orock_silv,silv_epsilon_w_r,
         x_center_points,y_center_points,summed_horz_vectors,summed_vert_vectors] = \
        o_isotope_invert(del18O_in_silv,Temp_in_silv,x_silv_in,z_silv_in,d18Oinit_silv)
        silv_water_initial[irun] = water_initial_silv
        silv_water_outgoing[irun] = water_outgoing_silv
        silv_moles_fluid[irun] = moles_fluid_silv
        silv_W_R[irun] = W_R_taper
        silv_final_rock[irun] = del18O_in_silv.mean()
        silv_moles_Orock[irun] = moles_Orock_silv
        fig1=plt.figure(num=1); plt.clf()    
        cplot=plt.contourf(xgrid,zgrid,del18O)
        plt.plot(randx,randz,'w^')
        plt.xlabel('Horizontal (m)');plt.ylabel('Vertical (m)');
        cm = plt.cm.get_cmap('RdYlBu')
        cbar=plt.colorbar()
        cbar.ax.set_ylabel('$\delta^{18}$O',rotation=270)
        
        fig2 = plt.figure(num=2); plt.clf()
        plt.hist(silv_water_initial)
    
    return [silv_water_initial,randx,randz,randd18O,randtemp,
            silv_x_grid,silv_z_grid,silv_iso_grid,temp_grid_silv,silv_W_R,
            x_center_points,y_center_points,summed_horz_vectors,summed_vert_vectors]