 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 17:10:24 2019

@author: benjohnson
"""

#%% Double checking method with Lake City Caldera, which should spit out a depleted (-8 to -14) initial water value
def plot_lc():
    # import sys
    # sys.path.insert(0, '/Users/bwj/Google_Drive/python_scripts/')
    # sys.path.insert(0, '/Users/bwj/Google_Drive/current_projects/geochemical_inverse_modeling')
    
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
    d18Oinit_LC=7.2
    
    
    data_raw = pd.read_csv('Larson_and_Taylor_Lake_City_caldera.csv')
    del18O_LC = pd.Series.tolist(data_raw['d18O'])
    elevation = pd.Series.tolist(data_raw['Elevation_m'])
    elevation = np.asarray(elevation)
    d18Of = np.nanmean(del18O_LC)
    latitude =  np.zeros([len(data_raw.deg_N_lat_37),1])
    longitude =  np.zeros([len(data_raw.deg_N_lat_37),1])
    UTM_local = np.zeros([len(data_raw.deg_N_lat_37),2])
    for i in range(len(data_raw.deg_N_lat_37)): 
        latitude[i] = (37+data_raw.deg_N_lat_37[i]/60)
        longitude[i] = -(107+data_raw.deg_W_long_107[i]/60)
        coord = utm.from_latlon(latitude[i],longitude[i],force_zone_number=13)
        UTM_local[i,:] = coord[0:2] #Easting first, Northing second
    
    x = UTM_local[:,0]; y = UTM_local[:,1]
    
    midline_points = np.array([[282000,294000],[4202000, 4204000]])  #first row are x points, second row y poitns
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
    ptemp = 200 #275
    temp_LC = np.multiply(ptemp,np.divide(depth,np.add(depth,ktemp)))
    num_runs=len(del18O_LC)
    LC_water_initial = np.zeros(num_runs)
    LC_water_outgoing = np.zeros(num_runs)
    LC_W_R = np.zeros(num_runs)
    LC_final_rock = np.zeros(num_runs)
    
    for irun in range(0,num_runs):
        del18O_in_LC = del18O_LC
        Temp_redux_LC_in = temp_LC 
        X_LC_in = along_line_distance
        Z_LC_in = elevation
        
        num_takeout = np.random.randint(0,len(del18O_in_LC))
    #    num_replace = np.random.randint(0,len(del18O_in_MW)) 
    #    del18O_in_MW[num_takeout] = del18O_in_MW[num_replace]
        
        del18O_in_LC = np.delete(del18O_in_LC,num_takeout)
        Temp_redux_LC_in =np.delete(Temp_redux_LC_in,num_takeout)
        X_LC_in =  np.delete(X_LC_in,num_takeout)
        Z_LC_in =np.delete(Z_LC_in,num_takeout)
        [LC_iso_grid,LC_x_grid,LC_z_grid,conc_matrix_LC,change_matrix_LC,qsolved_LC,\
         moles_fluid_LC,W_R_taper,water_initial_LC, water_outgoing_LC, temp_test_LC, temp_grid_LC,moles_Orock_LC,LC_epsilon_w_r,
         x_center_points,y_center_points,summed_horz_vectors,summed_vert_vectors]=\
        o_isotope_invert(del18O_in_LC,Temp_redux_LC_in,X_LC_in,Z_LC_in,d18Oinit_LC)
        LC_water_initial[irun] = water_initial_LC
        LC_water_outgoing[irun] = water_outgoing_LC
        LC_W_R[irun] = W_R_taper
        LC_final_rock[irun] = del18O_in_LC.mean()
        
    return[LC_water_initial,along_line_distance,elevation,del18O_LC,temp_LC,
           LC_x_grid,LC_z_grid,LC_iso_grid,temp_grid_LC,LC_W_R,
           x_center_points,y_center_points,summed_horz_vectors,summed_vert_vectors]