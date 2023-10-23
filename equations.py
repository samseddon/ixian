import numpy as np
from colorama import Fore
import time
import os
import inspect
from numba import jit


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


def calc_Q_coordinates(dec_image):
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
    return dec_image
