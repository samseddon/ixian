import glob
import inspect
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from colour_bar import code_TUD_cbar as cbar
from scipy import interpolate


def file_checker(s_num, s_ind, input_path):                                    
    """                                                                        
    """                                                                        
#    print("Function" \                                                        
#            + str(inspect.currentframe()).split(",")[-1][5:-1] \              
#            + " called from"\                                                 
#            + str(inspect.currentframe()).split(",")[1])                      
    files = glob.glob(os.path.join(input_path, str(s_num) + '*'))              
    file_name = max(files, key=os.path.getctime)                               
    file_name = file_name[len(input_path):]                                    
    if s_ind == -1:                                                            
        with open(os.path.join(input_path,file_name), 'rb') as handle:         
             read_file = pickle.load(handle)                                   
                                                                               
    elif s_ind == 0:                                                           
        file_name = file_name[:len(file_name)-9]+file_name[-7:]                
        with open(os.path.join(input_path,file_name), 'rb') as handle:         
             read_file = pickle.load(handle)                                   
    else:                                                                      
        file_name = file_name[:(len(file_name)-8)]+str(s_ind)+file_name[-7:]   
        with open(os.path.join(input_path,file_name), 'rb') as handle:         
             read_file = pickle.load(handle)                                   
    return read_file, file_name 



def plot(scan_num):
    data_location = "local/processed_files/"                                   
    f_1, f_name = file_checker(scan_num, -1, data_location)         
    q_1, q_2 = "qx", "qz"

    if q_1 == 'qx':
        values_1 = f_1.q_x
    if q_1 == 'qy':
        values_1 = f_1.q_y
    if q_1 == 'qz':
        values_1 = f_1.q_z
    if q_2 == 'qx':
        values_2 = f_1.q_x
    if q_2 == 'qy':
        values_2 = f_1.q_y
    if q_2 == 'qz':
        values_2 = f_1.q_z
    print(values_2)
    print(scan_num)
    volume = np.array(f_1.data)


    qxqz = np.mean(f_1.data, axis = 1)
    qyqz = np.mean(f_1.data, axis = 0)
    h, w = np.shape(qxqz)
    data = qxqz.ravel()
    fq_1, fq_2 = np.meshgrid(values_1, values_2)

    fq_2 = np.array(fq_2).flatten()  # this is y in y
    fq_1 = np.array(fq_1).flatten()



    #  FIGURE                                                                  
    grid_q = np.log10(data.reshape(h,w))
    fig= plt.figure(figsize = (4,4))                                           
    ax = plt.subplot(111)                                                      
    vmin = 1                                                                   
    vmax = 6                                                                   
    
    img = ax.contourf(grid_q, origin='lower',                                  
                       # extent=[mesh_q_2_min, mesh_q_2_max, mesh_q_1_min,mesh_q_1_max],
            levels = 256, vmin = vmin, vmax = vmax,                            
                         cmap=cbar)                                            
    ax.contour(grid_q, origin='lower',                               
                    #extent=[mesh_q_2_min, mesh_q_2_max, mesh_q_1_min,mesh_q_1_max],
            levels = 32, colors = "#00305d", linewidths = 0.1, vmin = vmax)    
    plt.show()

def twoD_Gaussian(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):     
    x, y = xy                                                                  
    xo = float(xo)                                                             
    yo = float(yo)                                                             
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)  
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)   
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)  
    g = offset + np.abs(amplitude)*np.exp( - (a*((x-np.abs(xo))**2) + 2*b*(x-np.abs(xo))*(y-np.abs(yo))
                            + c*((y-np.abs(yo))**2)))                          
    return g.ravel()                                                           
                                                                               
                                                                               
def multi_Gauss(xy,  *parameter):                                              
    x, y = xy                                                                  
    J = np.zeros_like(x.ravel())                                               
    for i in range(0, len(parameter), 7):                                      
        new_params = (parameter[0+i],                                          
                     parameter[1+i],                                           
                     parameter[2+i],                                           
                     parameter[3+i],                                           
                     parameter[4+i],                                           
                     parameter[5+i],                                           
                     parameter[6+i])                                           
        J += twoD_Gaussian((x,y), *new_params)                                 
    return J                                                                   


def fit(scan_num):                                                            
#    fig = plt.figure(figsize = (12,5))                                        
#    ax_xz = plt.subplot(131, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
#    ax_yz = plt.subplot(132, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
#    ax_xy = plt.subplot(133, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
#                                                                              
    file_index = [-1] #for most recent processed data just put -1, otherwise 0 will do the first or its index
#    axis_1 = 'qx'                                                             
#    axis_2 = 'qz'                                                             
#    axis_3_limits = [1,240]                                                   
#    ##                                                                        
    data_location = "local/processed_files/"                                   
    f_1, f_name = file_checker(scan_num, file_index[0], data_location)         
    qyqz = np.mean(f_1.data, axis = 1)                                         
    h, w = np.shape(qyqz)                                                      
    data = qyqz.ravel()                                                        
    x = np.linspace(0, w, w)                                                   
    y = np.linspace(0, h, h)                                                   
    x, y = np.meshgrid(x, y)                                                   
                                                                               
    guess_1 = (max(data), 38, 42, 0.7, 1.1 , 0, 0)                                
    guess_2 = (max(data), 38, 28, 0.7, 1.1 , 0, 0)                                
    guess_3 = (0.5*max(data), 34, 35, 0.5, 0.5,  0, 0)                                
    guess_4 = (0.5*max(data), 35, 39, 0.5, 0.5,  0, 0)                                
    #guess_5 = (max(data)/8, 39, 35, 0.5, 0.5,  0, 0)                                
    #guess_6 = (max(data)/8, 34, 28, 0.5, 0.5,  0, 0)                                
    initial_guess = guess_1 + guess_2 + guess_3 + guess_4# + guess_5 + guess_6#+ guess_3 + guess_4  
    #print(initial_guess)                                                      
    #initial_guess = (max(data),25, 30, 20, 20,0, 50)                          
    popt, pcov = curve_fit(multi_Gauss, (x, y), data, p0=initial_guess)        
                                                                               
    #popt, pcov = curve_fit(multi_Gauss, (x, y), data, p0=initial_guess)       
    #print(popt)                                                               
    data_fitted = multi_Gauss((x, y), *popt)                                   
                                                                               
                                                                               
                                                                               
                                                                               
    fig= plt.figure(figsize = (4,4))                                           
    ax = plt.subplot(111)                                                      
    vmin, vmax = 1, 6                                                          
    ax.contourf(np.log10(qyqz), origin = "lower", levels = 256, cmap = cbar)#, vmin = vmin, vmax = vmax)
    for i in range(0, len(initial_guess), 7):

   ##     ax.contour(qyqz, origin = "lower", levels = 64, linewidths = 0.1, colors="#00305d")
    #    print(popt[i],popt[i+1], popt[i+2])
        new_params = (popt[0+i],
                      popt[1+i],
                      popt[2+i],
                      popt[3+i],
                      popt[4+i],
                      popt[5+i],
                      popt[6+i])
        print(new_params)
    #    #plt.imshow(twoD_Gaussian((x,y), *new_params).reshape(h,w))
    #    #plt.show()
        ax.contour(x, y, twoD_Gaussian((x,y), *new_params).reshape(h,w),17, colors = "#00305d", linewidths = 0.2)
    #    _+=1
    #ax.contourf(np.log10(data_fitted.reshape(h,w)), origin = "lower", levels = 256, cmap = cbar)#, vmin = vmin, vmax = vmax)
    ax.contour(x, y, np.log10(data_fitted.reshape(h, w)), 18, colors = "#00305d", linewidths = 0.1, vmin=1, vmax=6)
    plt.show()
#    print('plotting', f_name)
#    vmin = 1
#    vmax = 6
#    fig, ax_xz = rsm_plot(f_1, f_name, "qx", "qz",\
#            axis_3_limits, data_location + "images/", scan_num, fig, ax_xz)
#    fig, ax_yz = rsm_plot(f_1, f_name, "qy", "qz",\
#            axis_3_limits, data_location + "images/",  scan_num, fig, ax_yz)
#    fig, ax_xy = rsm_plot(f_1, f_name, "qx", "qy",\
#            axis_3_limits, data_location + "images/", scan_num, fig, ax_xy)
#    plt.tight_layout()
#    plt.show()
