import numpy as np
from colorama import Fore
import time
import os
import inspect
import pickle
import math
from numba import jit


def pickle_jar(filename, pickled_obj):
    with open(filename,'wb') as handle:
        pickle.dump(pickled_obj, handle, protocol=pickle.HIGHEST_PROTOCOL)


def pickle_unjar(filename):
    with open(filename, 'rb') as handle:
        return pickle.load(handle)


def progress_bar(progress, total,start_t): #here progress is an integer for a loop say, and total is a length of the loop
    percent = 100 * (progress/float(total))
    bar = '█' * int(percent/2) + '-' *int((100-percent)/2)
    time_n = (time.time()-start_t)*(100/percent)
    time_elapsed = (time.time()-start_t)
    if progress == 1:
        print(Fore.RESET + f"\n|{bar}| {percent:.2f}%", end = ", eta = " + str(int((time_n-time_elapsed))) +" seconds "+"\r")
    print(Fore.RESET + f"\r|{bar}| {percent:.2f}%", end = ", eta = " + str(int((time_n-time_elapsed))) +" seconds "+"\r")
    if progress==total:
        print(Fore.GREEN + f"\r|{bar}| {percent:.2f}%" , end = ", took "+str(int(time_elapsed))+ " seconds "+ Fore.RESET + "\n\n") 

def existential_check(o_f, f_s, folder): #file name, file suffix, destination folder
    print("Function" \
        + str(inspect.currentframe()).split(",")[-1][5:-1] \
        + " called from"\
        + str(inspect.currentframe()).split(",")[1])
    f_n = str(o_f + f_s)
    while os.path.exists(folder + f_n):
        if len(f_n) > len(o_f + f_s):
            filenum = int(f_n[(len(o_f + f_s) - (len(f_s)-1)):-len(f_s)]) + 1 
            f_n = f_n[:(len(o_f + f_s) - (len(f_s)-1))] + str(filenum) + f_n[-len(f_s):]
        else:
            f_n = f_n[:-len(f_s)] + '_1' + f_n[-len(f_s):]
    return folder + f_n 


@jit(nopython=True, cache=True)                                                                       
def calc_qx(WL, O, TT_hor, TT_ver): # c1 = i  c2 = j
    qx = ((2 * np.pi / WL)                                                    
            * (np.cos(O / 180 * np.pi)                                     
            - np.cos(TT_hor                                 
                    / 180                                                  
                    * np.pi                                                
                    - O                                                    
                    / 180                                                  
                    * np.pi)                                               
             * np.cos(TT_ver                                 
                     / 180                                                 
                     * np.pi)))                                            
    return(qx) #0.02748846124221295           

@jit(nopython=True, cache=True)                                                                       
def calc_qy(WL, O, TT_hor, TT_ver):
    qy =  ((2 * np.pi / WL)
            * np.sin(TT_ver
                    / 180
                    * np.pi)
            * np.cos(TT_hor
                    / 180
                    * np.pi
                    - O
                    / 180
                    * np.pi))
    return qy

@jit(nopython=True, cache=True)                                                                       
def calc_qz(WL, O, TT_hor):
    qz =    ((2 * np.pi / WL)
            * (np.sin(O
                     / 180
                     * np.pi)
            + np.sin(TT_hor
                    / 180
                    * np.pi
                    - O
                    / 180
                    * np.pi)))
    return qz


def calc_Q_coordinates(filename):
    dec_image = pickle_unjar(filename)
    for coordinate_1 in range(np.shape(dec_image.data)[1]):
        row = []
        row_x = []
        row_y = []
        row_z = []

        for coordinate_2 in range(np.shape(dec_image.data)[0]):
            dec_image.pixel_list.append([calc_qx(dec_image.WAVELENGTH, dec_image.OMEGA, dec_image.two_theta_horizontal_range[coordinate_2], dec_image.two_theta_vertical_range[coordinate_1]),
                             calc_qy(dec_image.WAVELENGTH, dec_image.OMEGA, dec_image.two_theta_horizontal_range[coordinate_2], dec_image.two_theta_vertical_range[coordinate_1]),
                             calc_qz(dec_image.WAVELENGTH, dec_image.OMEGA, dec_image.two_theta_horizontal_range[coordinate_2]),
                             dec_image.data[coordinate_2, coordinate_1]])

            row_x.append(dec_image.pixel_list[-1][0])
            row_y.append(dec_image.pixel_list[-1][1])
            row_z.append(dec_image.pixel_list[-1][2])
            row.append([row_x[-1], row_y[-1], row_z[-1]])
        dec_image.Q_x.append(row_x)
        dec_image.Q_y.append(row_y)
        dec_image.Q_z.append(row_z)
        dec_image.Q_coords.append(row)
    pickle_jar(filename, dec_image)
    return 

def nr_pts_finder(filename):
    set_x = set()                                                              
    set_y = set()
    set_z = set()
    list_x = []                                                            
    list_y = []                                                            
    list_z = []                                                            
    dec_image = pickle_unjar(filename)                                    
    temp_array = np.array(dec_image.pixel_list)                            
                                                                           
#    for _ in range(len(dec_image.pixel_list)):                             
#        list_x.append(dec_image.pixel_list[_][0])                         
#        list_y.append(dec_image.pixel_list[_][1])                         
#        list_z.append(dec_image.pixel_list[_][2])                         
#        total_list.append(dec_image.pixel_list[_])                        
                                                                           
    list_x = np.round(temp_array[:,0], 3)                                  
    list_y = np.round(temp_array[:,1], 3)                                  
    list_z = np.round(temp_array[:,2], 3)                                  
    for number in np.unique(list_x):                                       
        set_x.add(number)                                                  
    for number in np.unique(list_y):                                       
        set_y.add(number)                                                  
    for number in np.unique(list_z):                                       
        set_z.add(number)                                                  
#   set_x = set.add(list_x)                                                
#   set_y = set.add(list_y)                                                
#   set_z = set.add(list_z)                                                
    del dec_image
    return [list(set_x), list(set_y), list(set_z)]

@jit(cache = True)
def find_q_linear_fast(
    qx,
    qy,
    qz,
    QX_MAX, 
    QY_MAX, 
    QZ_MAX, 
    QX_MIN, 
    QY_MIN, 
    QZ_MIN, 
    QX_GRAD,
    QY_GRAD,
    QZ_GRAD):
    if qx < QX_MIN \
            or qx >= QX_MAX \
            or qy <  QY_MIN \
            or qy >= QY_MAX \
            or qz <  QZ_MIN \
            or qz >= QZ_MAX:
        return (False,[0,0,0])
    else:
        return(True, [math.floor((qx-QX_MIN)/QX_GRAD)+1,
                      math.floor((qy-QY_MIN)/QY_GRAD)+1,
                      math.floor((qz-QZ_MIN)/QZ_GRAD)+1])
