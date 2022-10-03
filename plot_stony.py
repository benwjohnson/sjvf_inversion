#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  2 15:44:37 2018

@author: benjohnson
"""

##----- Inversion for Stony mountain area ---------##
#add path to python scripts 
def plot_stony():
    import sys
    # sys.path.insert(0, '/Users/bwj/Google_Drive/python_scripts/')
    # sys.path.insert(0, '/Users/bwj/Google_Drive/current_projects/geochemical_inverse_modeling/')
    
    import pandas as pd 
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.tri as tri
    import scipy as sp
    from scipy.interpolate import griddata
    import utm
    from matrix_average import matrix_average 
    from o_isotope_invert import o_isotope_invert
    
    font = {'family' : 'sans-serif',
            'weight' : 'bold',
            'size'   : 12}
    
    matplotlib.rc('font', **font)
    plt.close('all')
    data_raw = pd.read_csv('stony_mountain.csv')
    
    del18O_SM = pd.Series.tolist(data_raw['d18O'])
    elevation = pd.Series.tolist(data_raw['Elevation'])
    elevation = np.asarray(elevation)
    d18Of = np.nanmean(del18O_SM)
    latitude =  pd.Series.tolist(data_raw['Lat'])
    longitude =  pd.Series.tolist(data_raw['Long'])
    UTM_local = np.zeros([len(data_raw.Long),2])
    for i in range(len(data_raw.Long)): 
        coord = utm.from_latlon(latitude[i],longitude[i],force_zone_number=13)
        UTM_local[i,:] = coord[0:2] #Easting first, Northing second
    
    x = UTM_local[:,0]; y = UTM_local[:,1]
    
       
    plt.figure(num=1); plt.clf()
    plt.plot(longitude,latitude,'o')
    plt.xlabel('Long'); plt.ylabel('Lat')
    
    plt.figure(num=2); plt.clf()
    plt.plot(UTM_local[:,0],UTM_local[:,1],'o')
    plt.xlabel('Easting'); plt.ylabel('Northing')
    
    
    midline_points = np.array([[256049,258338],[4208016, 4207876]])  #first row are x points, second row y poitns
    midline_slope =  (midline_points[1,1] - midline_points[1,0])/(midline_points[0,1]-midline_points[0,0])
    midline_intercept = midline_points[1,0]-midline_slope*midline_points[0,0]
    
    numpoints = len(x)
    
    #calculate slope for each point, assuming intercepts midline at 90 degree angle
    points_slope = np.zeros([numpoints,1])
    point_angle = np.zeros([numpoints,1])
    point_distance = np.zeros([numpoints,1])
    along_line_distance = np.zeros([numpoints,1])
    for ipoint in range(len(x)):
        points_slope[ipoint] = (y[ipoint]-midline_intercept)/x[ipoint]
        point_angle[ipoint] = np.arctan((midline_slope - points_slope[ipoint])/(1+points_slope[ipoint]*midline_slope))
        point_angle[ipoint] = np.rad2deg(point_angle[ipoint])
        point_distance[ipoint] = np.sqrt(x[ipoint]**2+(y[ipoint] - midline_intercept)**2)
        along_line_distance[ipoint] = np.sin(point_angle[ipoint])*point_distance[ipoint]
    
    along_line_distance = along_line_distance.reshape(len(along_line_distance))
    
    depth = np.subtract(max(elevation),elevation)  
    
    
    
    ktemp = 100
    aa = np.add(depth,ktemp); bb = np.divide(depth,aa)
    temp_SM_raw = np.multiply(200,bb)
    #temp_SM = np.add(depth,20)
    #Inversion
    num_runs = depth.size
    d18Oinit_SM = 7.5 #Crowley and Giletti, 1983
    SM_water_initial = np.zeros(num_runs)
    SM_water_outgoing = np.zeros(num_runs)
    SM_W_R = np.zeros(num_runs)
    SM_final_rock = np.zeros(num_runs)
    
    for irun in range(0,num_runs):
        del18O_in_SM = del18O_SM
        Temp_redux_SM_in = temp_SM_raw 
        X_SM_in = along_line_distance
        Z_SM_in = elevation
        
        num_takeout = np.random.randint(0,len(del18O_in_SM))
    
        
        del18O_in_SM = np.delete(del18O_in_SM,num_takeout)
        Temp_redux_SM_in =np.delete(Temp_redux_SM_in,num_takeout)
        X_SM_in =  np.delete(X_SM_in,num_takeout)
        Z_SM_in =np.delete(Z_SM_in,num_takeout)
       
        [SM_iso_grid,SM_x_grid,SM_z_grid,\
                conc_matrix_SM,change_SM,qsolved_SM,moles_fluid_SM,W_R_taper_SM,\
                water_initial_SM, water_outgoing_SM,temp_SM,temp_grid_SM,SM_moles_O_rock,epsilon_SM,
                x_center_points,y_center_points,summed_horz_vectors,summed_vert_vectors] = \
        o_isotope_invert(del18O_in_SM,Temp_redux_SM_in,X_SM_in,Z_SM_in,d18Oinit_SM)
    
    
    
        SM_water_initial[irun] = water_initial_SM
        SM_water_outgoing[irun] = water_outgoing_SM
        SM_W_R[irun] = W_R_taper_SM
        SM_final_rock[irun] = del18O_in_SM.mean()
    
    
    # fignum=3    
    # fig3=plt.figure(num=fignum); plt.clf(); #plt.title('grid')
    # plt.contourf(SM_x_grid,SM_z_grid,SM_iso_grid,cmap=plt.cm.BrBG)
    # plt.title('Stony Mountain')
    # plt.plot(along_line_distance,elevation,'w^',markeredgecolor='k',markersize=12)
    # cm = plt.cm.get_cmap('RdYlBu')
    # cbar=plt.colorbar()
    # cbar.ax.get_yaxis().labelpad = 15
    # cbar.ax.set_ylabel('$\delta^{18}$O',rotation=270)    
        
    
    
    #%% Inverting same data, but using cross section from Taylor, 1974
    #
    #data_taylor = pd.read_csv('stony_taylor_xsection.csv')
    #x = pd.Series.tolist(data_taylor['x'])
    #z = pd.Series.tolist(data_taylor['z'])
    #dO_taylor = pd.Series.tolist(data_taylor['d18O'])
    #temp_taylor = pd.Series.tolist(data_taylor['temp'])
    #
    #num_runsT = len(x)
    #
    #SMT_water_initial = np.zeros(len(dO_taylor))
    #SMT_water_outgoing = np.zeros(len(dO_taylor))
    #SMT_W_R = np.zeros(len(dO_taylor))
    #SMT_final_rock = np.zeros(len(dO_taylor))
    #
    #for irun in range(0,num_runsT):
    #    del18O_in_taylor = dO_taylor
    #    Temp_taylor_in = temp_taylor 
    #    X_taylor_in = x
    #    Z_taylor_in = z
    #    
    #    num_takeout = np.random.randint(0,len(dO_taylor))
    #
    #    
    #    del18O_in_taylor = np.delete(del18O_in_taylor,num_takeout)
    #    Temp_taylor_in = np.delete(Temp_taylor_in,num_takeout)
    #    X_taylor_in =  np.delete(X_taylor_in,num_takeout)
    #    Z_taylor_in = np.delete(Z_taylor_in,num_takeout)
    #   
    #    [SMT_iso_grid,SMT_x_grid,SMT_z_grid,\
    #            conc_matrix_SMT,change_SMT,qsolved_SMT,moles_fluid_SMT,W_R_taper_SMT,\
    #            water_initial_SMT, water_outgoing_SMT,temp_SMT,temp_grid_SMT,SMT_moles_O_rock,epsilon_SMT] = \
    #    o_isotope_invert(del18O_in_taylor,Temp_taylor_in,X_taylor_in,Z_taylor_in,d18Oinit_SM)
    #
    #
    #
    #    SMT_water_initial[irun] = water_initial_SMT
    #    SMT_water_outgoing[irun] = water_outgoing_SMT
    #    SMT_W_R[irun] = W_R_taper_SMT
    #    SMT_final_rock[irun] = del18O_in_taylor.mean()
        
    # fig4=plt.figure(num=4);
    # plt.hist(SM_water_initial)
    # #plt.hist(SMT_water_initial)
    # plt.legend(['F/T','Taylor'])
    # plt.xlabel('$\delta^{18}$O - Incoming')
    # plt.ylabel('Number')
    # plt.title('Stony Mountain') 
    
    return[SM_water_initial,along_line_distance,elevation,del18O_SM,temp_SM_raw,
           SM_x_grid,SM_z_grid,SM_iso_grid,temp_grid_SM,SM_W_R,
           x_center_points,y_center_points,summed_horz_vectors,summed_vert_vectors]