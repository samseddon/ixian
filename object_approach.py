import os
import sys
import inspect
import fabio
import time
from equations import *
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from dectris_class import Dectris_Image
from q_space_class import Q_Space
from colour_bar import code_TUD_cbar as cbar
from multiprocessing import Pool
"""
Created on Wed Jan 11 14:49:21 2023

@author: samseddon
"""
def Q_space_optimisation(scan_num, master_sets, dectris_objects_filenames):
    qlim_dir = "local/qlim/"
    potential_qlim_filename = qlim_dir\
                         + "qlim_"\
                         + str(scan_num[0])\
                         + ".txt"                                                        
    if "-s" in sys.argv or os.path.exists(potential_qlim_filename) == False:
        set_x, set_y, set_z = set(), set(), set()
        for element in master_sets:
            unique_elements_x = np.unique(element[0])
            for number in element[0]:
                set_x.add(number)
            unique_elements_y = np.unique(element[1])
            for number in element[1]:
                set_y.add(number)
            unique_elements_z = np.unique(element[2])
            for number in element[2]:
                set_z.add(number)
         
        NX_PTS, NY_PTS, NZ_PTS = len(set_x), len(set_y), len(set_z)
        qx_min, qy_min, qz_min, qx_max, qy_max, qz_max, q_min =\
                [], [], [], [], [], [], []
        tim_check = time.time()
        for filename in dectris_objects_filenames:
            dec_image = pickle_unjar(filename)
            qx, qy, qz = dec_image.q_lim()
            qx_min.append(qx[0])
            qy_min.append(qy[0])
            qz_min.append(qz[0])
            qx_max.append(qx[1])
            qy_max.append(qy[1])
            qz_max.append(qz[1])
            del dec_image
        
        Q_max = [np.amax(qx_max), np.amax(qy_max), np.amax(qz_max)]
        Q_min = [np.amin(qx_min), np.amin(qy_min), np.amin(qz_min)]
        
        Q_limit_dict_maker(
                           scan_num, 
                           Q_max, 
                           Q_min, 
                           NX_PTS,
                           NY_PTS,
                           NZ_PTS)


def calc_Q_of_dectris_objects(dectris_objects_filenames):
    pool = Pool()
    multiprocessing_result = pool.imap_unordered(calc_Q_coordinates, dectris_objects_filenames)
    all_images = list(multiprocessing_result)
    possible_NR_PTS = []
    multiprocessing_result = pool.imap_unordered(nr_pts_finder, dectris_objects_filenames)
    master_sets = list(multiprocessing_result)
    pool.close()
    pool.join()
    return master_sets

def image_read(scan_num, whole_image=False):
    with open("setup/data_info.txt","r") as inf:                               
        dict1 = eval(inf.read())                                               
    directory = dict1["directory"]                                             
    file_reference = dict1["file_reference"]
    final_file_list, param = XMaS_parameter_setup(scan_num)
    major_list = []
    pickle_names = []
    for image_number in range(len(final_file_list)):
        major_list.append([file_reference, image_number, \
                directory, final_file_list, param])
    for image_number in range(len(final_file_list)):
        temp_file_name = "local/temp/" + str(image_number) + "_dec_class" + ".pickle"
        temp_file = Dectris_Image(major_list[image_number], whole_image)
        print(temp_file.time)
        pickle_jar(temp_file_name, temp_file)
        pickle_names.append(temp_file_name)
    return pickle_names
    
def multi_image_populate(q_space, dectris_objects_filenames):
    new_start_t = time.time()
    for index, filename in enumerate(dectris_objects_filenames):
        dec_image = pickle_unjar(filename)
        q_space.populate_3D(dec_image)
        progress_bar(index + 1, len(dectris_objects_filenames), new_start_t)
    return q_space

def omega_scan(directory, file_reference, scan_num, whole_image=False):
    """ The main function to call, this finds the relevant data files, writes
    (if create_files = True) and reads the q_limits, initiles an array in 
    q-space, populates it with pixels and normalises the result. 
           Parameters: 
               directory(string)      : string pointing to data directory 
                                        structure
               file_reference(string) : unique file identifier for a given 
                                        experiment
               scan_num(array)        : array of scan numbers being read
               create_files(bool)     : should be True for first use, allows
                                        functions for given 
    """
    
    print('\nCreating parameter files')
    tim_check = time.time()
    dectris_objects_filenames = image_read(scan_num, file_reference, directory) 
    master_sets = calc_Q_of_dectris_objects(dectris_objects_filenames)
    Q_space_optimisation(scan_num, master_sets, dectris_objects_filenames)
    q_space = Q_Space(scan_num, symmetric = True, SPACE_3D = True)
    q_space = multi_image_populate(q_space, dectris_objects_filenames)
    q_space.normalise_3D()
    q_space.save()
    
    #new_filename = existential_check(str(scan_num[0]) + "_new_3d_fill",
    #                                 ".pickle",
    #                                 "local/processed_files/")
    #pickle_jar(new_filename, q_space)

def Q_limit_dict_maker(scan_num,
                       Q_max,
                       Q_min,
                       NX_PTS,
                       NY_PTS, 
                       NZ_PTS):

    qlim_dict = {'qz_min' : Q_min[2],                                          
                  'qz_max' : Q_max[2],                                         
                  'qx_min' : Q_min[0],                                         
                  'qx_max' : Q_max[0],                                         
                  'qy_min' : Q_min[1],                                         
                  'qy_max' : Q_max[1]}                                         
    
    qlim_dict['nx_pts'] = int(0.75*NX_PTS)                                     
    qlim_dict['ny_pts'] = int(0.75*NY_PTS)                                     
    qlim_dict['nz_pts'] = int(0.75*NZ_PTS)                                     
    qlim_dict["nr_pts"] = int(0.75 * min(NX_PTS, NY_PTS, NZ_PTS))
    qlim_dir = "local/qlim/"
   
    if os.path.exists(qlim_dir) == False:
        os.makedirs(qlim_dir)
                                                                           
    filename = qlim_dir\
               + "qlim_"\
               + str(scan_num[0])\
               + ".txt"                                                        
    with open(filename,'w') as inf:                                            
        inf.write(str(qlim_dict))                                              
    print('Created file',filename)   

def XMaS_parameter_setup(scan_num,
                         window_help = False):
                               
    """ When called, checks if you want to overwrite an existing parameter 
    file, and if so, calls a standard parameter file and writes a spot 
    specific parameter file after plotting a red square indicating the real
    space crop that will take place. To change this standard, central red spot
    open the now written spot parameters and tweak the last 4 parameters
            Parameters: 
               directory(string)      : string pointing to data directory 
               final_file_list(list)     : list of files for given scan number
               file_reference(string) : unique file identifier for a given 
                                        experiment
               spot_dict(dict)        : dictionary of scans and spots
               scan_num(list)        : list of scan numbers being read
    """
    with open("setup/data_info.txt","r") as inf:                               
        dict1 = eval(inf.read())                                               
    directory = dict1["directory"]                                             
    file_reference = dict1["file_reference"]

    filename = directory \
               + '/user_defined_parameters/param/param_' \
               + str(scan_num[0]) \
               + '.txt' 

    list_of_all_files = os.listdir(directory)
     
    experiment_files = [file for file in list_of_all_files \
                    if file.startswith(file_reference) \
                    and file.endswith(".edf")]
    
    file_list = [file for file in experiment_files \
                    if int(file.split("_")[-2]) \
                    in scan_num]
    
    order = [int(file.split("_")[-1][:-4]) for file in file_list]
    final_file_list = [file[0] for file in sorted(zip(file_list, order), key = lambda x: x[1])]

    # final_file_list a list of all files from relevant scan
    # NOTE -2 is a magic number, based on the XMaS file saving format
    
    with open('setup/XMaS_standard_param.txt', 'r') as inf:
        gen_param = eval(inf.read())
    REAL_HOR_LIM_LOW = gen_param['REAL_HOR_LIM_LOW'] 
    REAL_HOR_LIM_HIG = gen_param['REAL_HOR_LIM_HIG']
    REAL_VER_LIM_LOW = gen_param['REAL_VER_LIM_LOW']
    REAL_VER_LIM_HIG = gen_param['REAL_VER_LIM_HIG']

    #window_help = True#
    if window_help == True:
        final_file_list = sorted(final_file_list)
        k = int(len(final_file_list)/4)
        f = fabio.open(os.path.join(directory, final_file_list[k]))
        f0 = f.data
        f.close()
        
        k = int(len(final_file_list)/2)
        f = fabio.open(os.path.join(directory, final_file_list[k]))
        f1 = f.data    
        f.close()
        
        k = int(3*len(final_file_list)/4)
        f = fabio.open(os.path.join(directory, final_file_list[k]))
        f2 = f.data    
        f.close()
        fig = plt.figure(figsize = (12,5))
        ax_0 = plt.subplot(131, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
        ax_1 = plt.subplot(132, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
        ax_2 = plt.subplot(133, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
        ax_0.imshow(f0, origin = "lower")
        ax_1.imshow(f1, origin = "lower")
        ax_2.imshow(f2, origin = "lower")
        rect0 = patches.Rectangle((REAL_HOR_LIM_LOW, \
                                   REAL_VER_LIM_LOW),\
                                  REAL_HOR_LIM_HIG-REAL_HOR_LIM_LOW,\
                                  REAL_VER_LIM_HIG-REAL_VER_LIM_LOW,\
                                  linewidth=1,\
                                  edgecolor='r',\
                                  facecolor='none')
        rect1 = patches.Rectangle((REAL_HOR_LIM_LOW, \
                                   REAL_VER_LIM_LOW),\
                                  REAL_HOR_LIM_HIG-REAL_HOR_LIM_LOW,\
                                  REAL_VER_LIM_HIG-REAL_VER_LIM_LOW,\
                                  linewidth=1,\
                                  edgecolor='r',\
                                  facecolor='none')
        rect2 = patches.Rectangle((REAL_HOR_LIM_LOW, \
                                   REAL_VER_LIM_LOW),\
                                  REAL_HOR_LIM_HIG-REAL_HOR_LIM_LOW,\
                                  REAL_VER_LIM_HIG-REAL_VER_LIM_LOW,\
                                  linewidth=1,\
                                  edgecolor='r',\
                                  facecolor='none')
        ax_0.add_patch(rect0) 
        ax_1.add_patch(rect1) 
        ax_2.add_patch(rect2) 
        
        plt.show()
    return final_file_list, gen_param
