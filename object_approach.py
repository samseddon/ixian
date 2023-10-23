import os
import sys
import inspect
import fabio
import time
from equations import *
import pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches
from dectris_class import *
from colour_bar import code_TUD_cbar as cbar
from multiprocessing import Pool
"""
Created on Wed Jan 11 14:49:21 2023

@author: samseddon
"""

def data_fill(directory,output_folder,file_reference,scan_num,create_files):
    print("Function" \
        + str(inspect.currentframe()).split(",")[-1][5:-1] \
        + " called from"\
        + str(inspect.currentframe()).split(",")[1])
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
    files_location = os.listdir(directory+'data/')
    master_files = [m for m in files_location \
                    if m.startswith(file_reference) \
                    and m.endswith(".edf")]
    master_files = [c for c in master_files \
                    if int(c.split("_")[-2]) \
                    in scan_num]
    
    temp= []
    mag = []
    
    with open(directory+'user_defined_parameters/spot_dict.txt','r') as inf:
        spot_dict = eval(inf.read())
    
    print('\nCreating parameter files')
    
    if create_files == True:
        filename = directory\
                   + 'user_defined_parameters/qlim/qlim_'\
                   + spot_dict[str(scan_num[0])]\
                   + '.txt'
        
        parameter_setup(directory,
                        master_files, 
                        file_reference,
                        spot_dict,
                        scan_num)


    param = param_read(spot_dict,
                       scan_num[0],
                       directory)
        
    print('\nSlicing images and calulating pixel Q values..')
    
    q_unsorted = []
    limit_dict = []
    all_images = []
    start_t = time.time()
    major_list = []
    pickle_names = []
    for image_number in range(len(master_files)):
        major_list.append([file_reference, image_number, \
                directory, master_files, param])
    for image_number in range(len(master_files)):
        temp_file_name = "local/temp/" + str(image_number) + "_dec_class" + ".pickle"
        temp_file = Dectris_Image(major_list[image_number])
        pickle_jar(temp_file_name, temp_file)
        pickle_names.append(temp_file_name)
    print("slicing and images took ",time.time() - start_t)                                    
    tim_check = time.time()
#    for filename in pickle_names:
#        calc_Q_coordinates(filename)
    pool = Pool()
    multiprocessing_result = pool.imap_unordered(calc_Q_coordinates, pickle_names)
    all_images = list(multiprocessing_result)
    pool.close()
    pool.join()
    print("calculating Q values took ", time.time()- tim_check)
    tim_check = time.time()
    possible_NR_PTS = []
    c = 0
    #for filename in enumerate(pickle_names):
    #    dec_image = pickle_unjar(filename[1])
    #    if filename[0] % 1 == 0:
    #        print(filename[0])
    #        for _ in range(0,len(dec_image.pixel_list),1):#, int(len(dec_image.pixel_list)/1000)):
    #            #print(_)
    #            total_x.append(dec_image.pixel_list[_][0])
    #            total_y.append(dec_image.pixel_list[_][1])
    #            total_z.append(dec_image.pixel_list[_][2])
    #        possible_NR_PTS.append([len(np.unique(np.round(total_x,3))),
    #                           len(np.unique(np.round(total_y,3))),
    #                           len(np.unique(np.round(total_z,3)))])
    #        
    #    else:
    #        pass
    total_list = []
    master_sets = []
    
#    for file_name in pickle_names:
#        master_sets.append(nr_pts_finder(file_name))

    pool = Pool()
    multiprocessing_result = pool.imap_unordered(nr_pts_finder, pickle_names)
    master_sets = list(multiprocessing_result)
    pool.close()
    pool.join()




    set_x = set() 
    set_y = set()
    set_z = set()
    for element in master_sets:
        unique_elements_x = np.unique(element[0])
        for number in unique_elements_x:
            set_x.add(number)
        unique_elements_y = np.unique(element[1])
        for number in unique_elements_y:
            set_y.add(number)
        unique_elements_z = np.unique(element[2])
        for number in unique_elements_z:
            set_z.add(number)
    
    NR_PTS = min(len(set_x), len(set_y), len(set_z))
#   for _ in range(len(total_list)):                                                     
#           x.append(total_list[_][0])                                                       
#           y.append(total_list[_][1])                                                       
#           z.append(total_list[_][2])                                                       
#    print(len(np.unique(np.round(x,3))))
#    print(len(np.unique(np.round(y,3))))  
#    print(len(set_x), len(set_y), len(set_z))
#    print(len(np.unique(np.round(x,3))),

    #              len(np.unique(np.round(y,3))),
    #              len(np.unique(np.round(z,3))))
    #NR_PTS = (min(len(np.unique(np.round(x,3))),
    #              len(np.unique(np.round(y,3))),
    #              len(np.unique(np.round(z,3)))))
#    NR_PTS = min(possible_NR_PTS)
    print(NR_PTS)
    print("NR_PTS Calculation took ", time.time() - tim_check)
    print('\nFinding q limits from sliced data and optimising Q_space mesh')
    qx_min = []
    qy_min = []
    qz_min = []
    qx_max = []
    qy_max = []
    qz_max = []
    q_min = []
    tim_check = time.time()
    for filename in pickle_names:
        dec_image = pickle_unjar(filename)
        qx, qy, qz = dec_image.q_lim()
        qx_min.append(qx[0])
        qy_min.append(qy[0])
        qz_min.append(qz[0])
        qx_max.append(qx[1])
        qy_max.append(qy[1])
        qz_max.append(qz[1])
        del dec_image
    
    # NOTE here the qmax/qlim are found, needs to be put into a function.
    Q_max = [np.amax(qx_max), np.amax(qy_max), np.amax(qz_max)]
    Q_min = [np.amin(qx_min), np.amin(qy_min), np.amin(qz_min)]
    print(Q_max, Q_min)

    if create_files == True:
              # or os.path.exists("local/qlim"\
              #                               + "qlim_"\
              #                               + str(scan_num[0])\
              #                               + ".txt")\
              #                               == False:                 
        Q_limit_dict_maker(directory, spot_dict, scan_num, Q_max, Q_min, NR_PTS)
    print("found q limits in this time", time.time() - tim_check)

    q_space = Q_Space(scan_num, spot_dict, directory)
    new_start_t = time.time()
    print('Populating Q_space with pixels')
    for filename in enumerate(pickle_names):
        dec_image = pickle_unjar(filename[1])
        q_space.populate_3D(dec_image)
        progress_bar(filename[0] + 1,len(master_files),new_start_t)


    q_space.normalise_3D()


    orig_filename = str(scan_num[0])+'_new_3d_fill'
    suffix = '.pickle'
    new_filename = existential_check(orig_filename,
                                     suffix, 
                                     directory + 'processed_files/')
    with open(new_filename,'wb') as handle:
        pickle.dump(q_space, handle, protocol=pickle.HIGHEST_PROTOCOL)
#    
    print("Complete, total duration: {}".format(int(time.time() 
                                                    - start_t))
                                                    + ' seconds')
#


def Q_limit_dict_maker(directory, spot_dict, scan_num, Q_max, Q_min, NR_PTS):  
    qlim_dict = {'qz_min' : Q_min[2],                                          
                  'qz_max' : Q_max[2],                                         
                  'qx_min' : Q_min[0],                                         
                  'qx_max' : Q_max[0],                                         
                  'qy_min' : Q_min[1],                                         
                  'qy_max' : Q_max[1]}                                         
    
    qlim_dict['nr_pts'] = int(0.75*NR_PTS)                                     
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

def param_read(spot_dict,scan_num,directory):
    print("Function" \
        + str(inspect.currentframe()).split(",")[-1][5:-1] \
        + " called from"\
        + str(inspect.currentframe()).split(",")[1])
    """ Takes the spot string defined in the spot dictionary and opens the 
    corresponding parameter file, returning it.
           Parameters: 
               spot_dict(dict)        : dictionary of scans and spots
               scan_num(array)        : array of scan numbers being read
               directory(string)      : string pointing to data directory 
                                        structure
           Returns:
               newly_read_param_file  : as it says on the tin
    """
    import_var = directory\
                 + "user_defined_parameters/param/param_" \
                 + spot_dict[str(scan_num)]\
                 + '.txt'
    if os.path.exists(import_var)==False:
        raise(ValueError("Param_dict for spot non-existant"))
    
    with open(import_var,'r') as inf:
        newly_read_param_dict = eval(inf.read())
    return newly_read_param_dict


def parameter_setup(directory,
                    master_files,
                    file_reference,
                    spot_dict,
                    scan_num):
    print("Function" \
        + str(inspect.currentframe()).split(",")[-1][5:-1] \
        + " called from"\
        + str(inspect.currentframe()).split(",")[1])
    """ When called, checks if you want to overwrite an existing parameter 
    file, and if so, calls a standard parameter file and writes a spot 
    specific parameter file after plotting a red square indicating the real
    space crop that will take place. To change this standard, central red spot
    open the now written spot parameters and tweak the last 4 parameters
            Parameters: 
               directory(string)      : string pointing to data directory 
               master_files(list)     : list of files for given scan number
               file_reference(string) : unique file identifier for a given 
                                        experiment
               spot_dict(dict)        : dictionary of scans and spots
               scan_num(list)        : list of scan numbers being read
    """
    filename = directory \
               + '/user_defined_parameters/param/param_' \
               + spot_dict[str(scan_num[0])] \
               + '.txt' 
    
    with open('setup/standard_param.txt', 'r') as inf:
        gen_param = eval(inf.read())
    
    REALX_LIM_LOW = gen_param['realx_lim_low'] 
    REALX_LIM_HIG = gen_param['realx_lim_hig']
    REALY_LIM_LOW = gen_param['realy_lim_low']
    REALY_LIM_HIG = gen_param['realy_lim_hig']

    k = int(len(master_files)/4)
    f = fabio.open(os.path.join(directory+'data/', master_files[k]))
    f1 = f.data
    f.close()
    
    k = int(len(master_files)/2)
    f = fabio.open(os.path.join(directory+'data/', master_files[k]))
    f2 = f.data    
    f.close()
    
    k = int(3*len(master_files)/4)
    f = fabio.open(os.path.join(directory+'data/', master_files[k]))
    f3 = f.data    
    f.close()
    

    with open(filename,'w') as inf:
        inf.write(str(gen_param))
    print('Created file',
          filename,
          '\n Change slice region in parameter file')



