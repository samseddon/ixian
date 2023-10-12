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
