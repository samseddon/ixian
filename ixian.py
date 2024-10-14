import os
import sys
import inspect
import time
#from pixel_selection import data_fill
import numpy as np
from matplotlib.cm import ScalarMappable
#from rsm_object_plotter import file_checker,  rsm_plot, slicer_and_dicer_3000, test_slicer, threeD_rsm_plot
from object_approach import omega_scan
import glob
import pickle
from scipy import interpolate
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from plotting_functions import plot, fit
from object_approach import image_read, calc_Q_of_dectris_objects,\
        Q_space_optimisation, multi_image_populate
from q_space_class import Q_Space
from equations import pickle_jar, pickle_unjar




# NOTE TO DO 
'''
check if optimal Qspace already exists 
simplify file checker 
'''
def file_checker(s_num, s_ind, input_path):                                    
    """                                                                        
    """                                                                        
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

def local_setup():
    if os.path.exists("local/") == False:
        os.mkdir("local/")
    if os.path.exists("local/processed_files/") ==  False:
       os.mkdir("local/processed_files/")
    if os.path.exists("local/temp/") ==  False:
        os.mkdir("local/temp")  
    if os.path.exists("local/qlim/") ==  False:
        os.mkdir("local/qlim")  
    if os.path.exists("local/images/") ==  False:
        os.mkdir("local/images")  
    if os.path.exists("local/gif/") == False:
        os.mkdir("local/gif/")
    if os.path.exists("local/gif/images/") ==  False:
        os.mkdir("local/gif/images")  
    if os.path.exists("setup/data_info.txt") == False \
            or "-dataset" in sys.argv:
                data_setup()

def oscan(scan_num):
    dectris_objects_filenames = image_read(scan_num) 
    master_sets = calc_Q_of_dectris_objects(dectris_objects_filenames)         
    Q_space_optimisation(scan_num, master_sets, dectris_objects_filenames)        
    
    q_space = Q_Space(scan_num, symmetric = True, SPACE_3D = True)             
    q_space.multi_image_populate(dectris_objects_filenames)         
    q_space.normalise_3D()       
    q_space.save()
    
    data_location = "local/processed_files/"
    loaded_qspace, f_name = file_checker(scan_num, -1, data_location)
    loaded_qspace.plot()
        

    if "-singplot" in sys.argv:
        if len(scan_num) == 0:
            scan_num = scan_num_input()
        data_location = "local/processed_files/"
        loaded_qspace, f_name = file_checker(scan_num, -1, data_location)
        loaded_qspace.plot_sing()
            
def tscan(scan_num):
    dectris_object_filenames = image_read(scan_num) 
    master_sets = calc_Q_of_dectris_objects(dectris_object_filenames)         
    Q_space_optimisation(scan_num, master_sets, dectris_object_filenames)        
    q_space_list = []    
    for filename in dectris_object_filenames:
        q_space = Q_Space(scan_num, symmetric = False, SPACE_3D = True)
        q_space.sing_image_populate(filename)
        q_space.normalise_3D()
        q_space_list.append(q_space)
    filename ="local/processed_files" + str(scan_num[0]) + "_1D_tscan"+ ".pickle"
    pickle_jar(filename, q_space_list)

def mess_around(scan_num):
    filename ="local/processed_files" + str(scan_num[0]) + "_1D_tscan"+ ".pickle"
    q_space_list = pickle_unjar(filename)
    int_area = []
    temps = []
    times = []
    mag_fields = []
    for q_space in q_space_list: 
        times.append(q_space.TIME) 
        temps.append(q_space.TEMP)
        mag_fields.append(q_space.MAG_FIELD)
        int_area.append(sum(sum(sum(q_space.data))))
    min_times = min(times)
    norm_times = []
    for time in times:
        norm_times.append(time-min_times)
    gif_filenames = []
    for q_space in q_space_list:
        gif_filenames.append(q_space.plot_for_gif())
     
    import imageio
    images = []
    print(gif_filenames)
    with imageio.get_writer('test.gif', mode='I') as writer:
        for filename in gif_filenames:
            image = imageio.imread(filename)
            writer.append_data(image)



    #plt.figure()
    #plt.plot(mag_fields, int_area)
    #plt.show()
    #q_space = Q_Space(scan_num, symmetric = True, SPACE_3D = True)             
    #q_space.multi_image_populate(dectris_objects_filenames)         
    #q_space.normalise_3D()       
    #q_space.save()
    #

    #data_location = "local/processed_files/"
    #loaded_qspace, f_name = file_checker(scan_num, -1, data_location)
    #print(loaded_qspace)
    #print(loaded_qspace.time)
    #loaded_qspace.plot()

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

def add_history(scan_num):
    with open("setup/history.txt",'w') as inf:
        inf.write(str(scan_num))

def load_scan_num_history():
    if "-r" in sys.argv:
        try:
            with open("setup/history.txt","r") as inf:   
                return eval(inf.read())
        except FileNotFoundError:
                return scan_num_input()
    else:
        return scan_num_input()

if __name__ == "__main__":
    start_t = time.time()

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--rerun-last", action="store_true", help="runs last scan number")
    parser.add_argument("-Q", "--Q-Space-threeD", action="store_true", help="creates a Q-Space object from a given scan number (or multiple scan numbers) and populates it with the scans")
    parser.add_argument("-s", "--scratch", action="store_true", help="Used in the instance that you want to run the code compltely again, overwriting parameter files etc, that may have been manually changed. When not used, if a scan has been run before the existing paramter files will be used to save time")
    parser.add_argument("-T", "--timescan", action="store_true", help="creates a list of q_spaces, which can be then later opened and plotted as a function of time or motor position")
    parser.add_argument("-m", "--mess", action="store_true", help="opens a list of timescan generated qspaces and allows the plotting of these")
    args = parser.parse_args()
    
    local_setup()
    scan_num = load_scan_num_history()

    if args.Q_Space_threeD:
        oscan(scan_num) 
    if args.timescan:
        tscan(scan_num)
    if args.mess:
        mess_around(scan_num)
    add_history(scan_num)
    print("Complete, total duration: {}".format(int(time.time()                
                                                    - start_t))                
                                                    + ' seconds')

#def twoD_Gaussian(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
#    x, y = xy
#    xo = float(xo)
#    yo = float(yo)
#    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
#    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
#    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
#    g = offset + np.abs(amplitude)*np.exp( - (a*((x-np.abs(xo))**2) + 2*b*(x-np.abs(xo))*(y-np.abs(yo))
#                            + c*((y-np.abs(yo))**2)))
#    return g.ravel()
#
#
#def multi_Gauss(xy,  *parameter):
#    x, y = xy
#    J = np.zeros_like(x.ravel())
#    for i in range(0, len(parameter), 7):
#        new_params = (parameter[0+i],
#                     parameter[1+i],
#                     parameter[2+i],
#                     parameter[3+i],
#                     parameter[4+i],
#                     parameter[5+i],
#                     parameter[6+i])
#        J += twoD_Gaussian((x,y), *new_params)
#    return J
#
#def plot(scan_num):
##    fig = plt.figure(figsize = (12,5))
##    ax_xz = plt.subplot(131, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
##    ax_yz = plt.subplot(132, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
##    ax_xy = plt.subplot(133, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
##
#    file_index = [-1] #for most recent processed data just put -1, otherwise 0 will do the first or its index
##    axis_1 = 'qx'
##    axis_2 = 'qz'
##    axis_3_limits = [1,240]
##    ##
#    data_location = "local/processed_files/"
#    f_1, f_name = file_checker(scan_num, file_index[0], data_location)
#    qyqz = np.mean(f_1.data, axis = 1)
#    h, w = np.shape(qyqz)
#    data = qyqz.ravel()
#    x = np.linspace(0, w, w)
#    y = np.linspace(0, h, h)
#    x, y = np.meshgrid(x, y)
#    
#    guess_1 = (max(data)/4, 38, 30, 5, 5, 0, 0)
#    guess_2 = (max(data)/4, 37, 36, 5, 5, 0, 0)
#    guess_3 = (max(data)/4, 40, 34, 5, 5, 0, 0)
#    guess_4 = (max(data)/4, 38, 42, 5, 5, 0, 0)
#    initial_guess = guess_1 + guess_2 + guess_3+ guess_4 #+ guess_3 + guess_4
#    #print(initial_guess)
#    #initial_guess = (max(data),25, 30, 20, 20,0, 50)
#    popt, pcov = curve_fit(multi_Gauss, (x, y), data, p0=initial_guess)
#  
#    #popt, pcov = curve_fit(multi_Gauss, (x, y), data, p0=initial_guess)
#    #print(popt)
#    data_fitted = multi_Gauss((x, y), *popt)
#
#
#
#
#    fig= plt.figure(figsize = (4,4))
#    ax = plt.subplot(111)
#    vmin, vmax = 1, 6
#    fin = qyqz.flatten("F")
#    fq_1 = []
#    fq_2 = []
#    for i in range(len(y)):                                           
#        fq_2.append(x)                                                  
#                                                                               
#    fq_2 = np.array(fq_2).flatten()  # this is y in y                          
#                                                                               
#    c = 0  # incrementing variable                                             
#    nr = 0  # incrementing variable                                            
#    for j in range(len(y)):                                             
#        fq_1.append([y[c]] * len(x))                             
#        c = c + 1                                                              
#                                                                               
#    fq_1 = np.array(fq_1).flatten()  
#    mesh_q_1_min = np.min(fq_2)  # sets the minimum qx value for new mesh (np.min(fq_2) is the minimum of the input data)
#    mesh_q_1_max = np.max(fq_2)  # sets the maximum qx value for new mesh (np.max(fq_2) is the maximum of the input data)
#    mesh_q_2_min = np.min(fq_1)
#    mesh_q_2_max = np.max(fq_1)
#    nr_fqz = int(len(np.unique(fq_1)) * 1)  # sets the number of mesh points in qx  *1.0 keeps as original
#    nr_fqx = int(len(np.unique(fq_2)) * 1)  # sets the number of mesh points in qz  *1.0 keeps as original
#
#    # this is y for y, #grid qx
#    grid_qx, grid_qz = np.mgrid[mesh_q_1_min:mesh_q_1_max:(nr_fqx * 1j), mesh_q_2_min:mesh_q_2_max:(nr_fqz * 1j)]
#    print(fq_2, fq_1, grid_qx, grid_qz)
#    #grid_q = interpolate.griddata((fq_2, fq_1), fin, (grid_qx, grid_qz), method='nearest')
#    #plt.imshow(multi_Gauss((x,y), *popt).reshape(h,w))
#    _ = 0
#    ax.contourf(qyqz, origin = "lower", levels = 256, cmap = cmap)#, vmin = vmin, vmax = vmax)
#    #for i in range(0, len(initial_guess), 7):
#   
#   ##     ax.contour(qyqz, origin = "lower", levels = 64, linewidths = 0.1, colors="#00305d")
#    #    print(popt[i],popt[i+1], popt[i+2])
#    #    new_params = (popt[0+i],
#    #                  popt[1+i],
#    #                  popt[2+i],
#    #                  popt[3+i],
#    #                  popt[4+i],
#    #                  popt[5+i],
#    #                  popt[6+i])
#    #    #plt.imshow(twoD_Gaussian((x,y), *new_params).reshape(h,w))
#    #    #plt.show()
#    #    ax.contour(x, y, twoD_Gaussian((x,y), *new_params).reshape(h,w),17, colors = "#00305d", linewidths = 0.2)
#    #    _+=1
#    ax.contour(x, y, data_fitted.reshape(h, w), 18, colors = "#00305d", linewidths = 0.1)
#    plt.show()
##    print('plotting', f_name)
##    vmin = 1
##    vmax = 6
##    fig, ax_xz = rsm_plot(f_1, f_name, "qx", "qz",\
##            axis_3_limits, data_location + "images/", scan_num, fig, ax_xz)
##    fig, ax_yz = rsm_plot(f_1, f_name, "qy", "qz",\
##            axis_3_limits, data_location + "images/",  scan_num, fig, ax_yz)
##    fig, ax_xy = rsm_plot(f_1, f_name, "qx", "qy",\
##            axis_3_limits, data_location + "images/", scan_num, fig, ax_xy)
##    plt.tight_layout()
##    plt.show()
