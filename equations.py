import numpy as np
from colorama import Fore
import time
import os
import inspect



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


