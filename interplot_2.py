import matplotlib.pylab as plt                                                 
import pandas as pd                                                            
from scipy.interpolate import griddata                                         
import numpy as np                                                             
import math                                                                    
from colour_bar import code_TUD_cbar as cbar                                   
                                                                               
plt.style.use("seddon_TUD")                                                    

def twod_plot(f_1, f_name, axis_1, axis_2, axis_3_limits, output_folder):
    Q_coords = f_1[0]
    print(Q_coords)

#    x = []                                                                         
#    y = []                                                                         
#    z = []                                                                         
#    points = []                                                                    
#                                                                                   
#                                                                                   
#    ylim = [-0.01,0.01]                                                            
#                                                                                   
#                                                                                   
#    for _ in range(len(df[0])):                                                    
#        if ylim[0] <= df[1][_] <= ylim[1]:                                         
#            points.append((df[0][_],df[1][_]))                                     
#            x.append(df[0][_])                                                     
#            y.append(df[1][_])                                                     
#            z.append(df[2][_])                                                     
#                                                                                   
#                                                                                   
#    grid_x, grid_y = np.mgrid[min(x):max(x):128j, min(y):max(y):128j]              
#                                                                                   
#    x_min = min(x)                                                                 
#    x_max = max(x)                                                                 
#    y_min = min(y)                                                                 
#    y_max = max(y)                                                                 
#                                                                                   
#    XY = griddata(points,np.log(z),(grid_x, grid_y))                               
#    for i in range(len(grid_x)):                                                   
#        for j in range(len(grid_y)):                                               
#            if math.isnan(XY[j,i]) or math.isinf(XY[j,i]):                         
#                XY[j,i] = 0                                                        
#    fig, ax = plt.subplots(figsize = (5,4), dpi = 128)                             
#    ax.contourf(XY.T, cmap = cbar, vmin=1, levels = 64)                            
#    ax.contour(XY.T, levels = 16, colors = 'C3', linewidths = 0.1)                 
#    xtick_pos = [16,64,112]                                                        
#    ax.set_xticks(xtick_pos)                                                       
#    ax.set_xticklabels([round(x_min+(x_max-x_min)*(xtick_pos[0]/len(grid_x)),3),   
#                        round(x_min+(x_max-x_min)*(xtick_pos[1]/len(grid_x)),3),   
#                        round(x_min+(x_max-x_min)*(xtick_pos[2]/len(grid_x)),3)])  
#    ax.set_xlabel(r"\rm Q$_z"+r"$ $(\rm\AA^{-1})$")                                
#    ytick_pos = [16,64,112]                                                        
#    ax.set_yticks(ytick_pos)                                                       
#    ax.set_yticklabels([round(y_min+(y_max-y_min)*(yt0ck_pos[0]/len(grid_y)),3),   
#                        round(y_min+(y_max-y_min)*(ytick_pos[1]/len(grid_y)),3),   
#                        round(y_min+(y_max-y_min)*(ytick_pos[2]/len(grid_y)),3)])  
#    ax.set_ylabel(r"\rm Q$_x"+r"$ $(\rm\AA^{-1})$")                                
#    #plt.contourf(np.log10(z), extent = [np.min(x),np.max(x),np.min(y),np.max(y)], levels=64)
#    plt.show() 
