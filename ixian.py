import os
import sys
import inspect
#from pixel_selection import data_fill
import numpy as np
from matplotlib.cm import ScalarMappable
from rsm_object_plotter import file_checker,  rsm_plot, slicer_and_dicer_3000, test_slicer, threeD_rsm_plot
from object_approach import omega_scan
import matplotlib.pyplot as plt

# NOTE TO DO 
'''
Make command line accesible 
check if optimal Qspace already exists 
'''

def main(scan_num):
    if "-h" in sys.argv and os.path.exists("setup/history.txt"):
        scan_num = load_scan_num_history()
    if os.path.exists("local/processed_files/") ==  False:
       os.mkdir("local/processed_files/")
    if os.path.exists("local/temp/") ==  False:
        os.mkdir("local/")
        os.mkdir("local/temp")  
    if os.path.exists("local/qlim/") ==  False:
        os.mkdir("local/qlim")  
    
    if os.path.exists("setup/data_info.txt") == False \
            or "-dataset" in sys.argv:
                data_setup()

    if "-oscan" in sys.argv:
        if "-h" not in sys.argv:
            scan_num = scan_num_input()
        with open("setup/data_info.txt","r") as inf:   
            dict1 = eval(inf.read())

        directory = dict1["directory"]
        file_reference = dict1["file_reference"]
        



        omega_scan(directory, file_reference, scan_num, create_files = True)
         
    if "-plot" in sys.argv:
        if len(scan_num) > 0:
            plot(scan_num)

        else:
            scan_num = scan_num_input()
            plot(scan_num)
    return scan_num
            

def scan_num_input():
    print("please input scan number/s as integers separated by commas if "+\
            "appropriate")
    return [int(scan_num) for scan_num in input().split(",")]

def data_setup():
    print("please input the absolute path to data directory")
    while True:
        directory = input()
        print(directory)
        if os.path.exists(directory) == True:
            break
        print("unfortunately that directory does not exist"\
                +" please try another")
    
    data_setup = dict()
    data_setup["directory"] = directory
    
    print("Now please enter unique experiment identifier, if not known"+\
            " please immediately contact your beamline scientist ;)")
    
    
    data_setup["file_reference"] = input()
    print("Recording these in file, to edit these values run the code"+\
            "with -dataset flag, or delete file 'data_info.txt' in setup"+\
            "folder")

    with open("setup/data_info.txt",'w') as inf:
        inf.write(str(data_setup))

def plot(scan_num):
    fig = plt.figure(figsize = (12,5))
    ax_xz = plt.subplot(131, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
    ax_yz = plt.subplot(132, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
    ax_xy = plt.subplot(133, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')

    file_index = [-1] #for most recent processed data just put -1, otherwise 0 will do the first or its index
    axis_1 = 'qx'
    axis_2 = 'qz'
    axis_3_limits = [1,240]
    ##
    data_location = "local/processed_files/"
    f_1, f_name = file_checker(scan_num, file_index[0],data_location)
    print('plotting', f_name)
    vmin = 1
    vmax = 6
    fig, ax_xz = rsm_plot(f_1, f_name, "qx", "qz",\
            axis_3_limits, data_location + "images/", scan_num, fig, ax_xz)
    fig, ax_yz = rsm_plot(f_1, f_name, "qy", "qz",\
            axis_3_limits, data_location + "images/",  scan_num, fig, ax_yz)
    fig, ax_xy = rsm_plot(f_1, f_name, "qx", "qy",\
            axis_3_limits, data_location + "images/", scan_num, fig, ax_xy)
    plt.tight_layout()
    plt.show()

def add_history(scan_num):
    with open("setup/history.txt",'w') as inf:
        inf.write(str(scan_num))

def load_scan_num_history():
    with open("setup/history.txt","r") as inf:   
        return eval(inf.read())

    

if __name__ == "__main__":
    scan_num = []
    scan_num = main(scan_num)
    add_history(scan_num)
