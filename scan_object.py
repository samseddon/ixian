import os                           
import fabio                        
import time                         
from equations import *             
import pickle                       
import numpy as np                  
import matplotlib.pyplot as plt     
import seaborn as sns               
import matplotlib.patches as patches

"""
Created on Wed Jan 11 14:49:21 2023

@author: samseddon
"""

class pilatus_image():










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

    start_t = time.time()
    temp= []
    mag = []

    with open(directory+'user_defined_parameters/spot_dict.txt','r') as inf:
        spot_dict = eval(inf.read())

    print('\nCreating parameter files')

    if create_files == True:
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
    for _ in range(len(master_files)):
        :WQ

        q_unsorted_temp,limit_dict_temp = pixel_segmenter(_, 
                                                          directory,
                                                          master_files,
                                                          param,       
                                                          start_t,     
                                                          temp,        
                                                          mag,         
                                                          file_reference)
        
        q_unsorted.append(q_unsorted_temp)
        limit_dict.append(limit_dict_temp)
        progress_bar(k+1,len(master_files),start_t)

    print('\nFinding q limits from sliced data and optimising Q_space mesh')
                                                                            
    if create_files == True:                                                
        find_q_lim(q_unsorted,directory,spot_dict,scan_num,limit_dict)      
                                                                            
    q_final,q_x_fin,q_y_fin,q_z_fin,qlim_dict = q_array_init(param,         
                                                             spot_dict,     
                                                             scan_num,      
                                                             directory)     
    new_start_t = time.time()                                               
    print('Populating Q_space with pixels')                                 
                                                                            
    # NOTE These remain ugly as they are going to be put into a function soon
                                                                               
    for k in range(len(q_unsorted)): 
        for i in range(np.shape(q_unsorted[k])[1]):
            for j in range(np.shape(q_unsorted[k])[0]):
                q_index_dict = find_q_index(q_unsorted[k][j][i]['qx'],\
                        q_unsorted[k][j][i]['qy'],q_unsorted[k][j][i]['qz'],\
                        q_x_fin,q_y_fin,q_z_fin,qlim_dict)
                if q_index_dict == 1:
                    pass
                else:
                    q_final['q_data'][q_index_dict['i_x'],q_index_dict['i_y'],q_index_dict['i_z']] =\
                            q_unsorted[k][j][i]['data'] +\
                            q_final['q_data'][q_index_dict['i_x'],q_index_dict['i_y'],q_index_dict['i_z']]
                    q_final['q_idx'][q_index_dict['i_x'],q_index_dict['i_y'],q_index_dict['i_z']] =\
                            1 + q_final['q_idx'][q_index_dict['i_x'],q_index_dict['i_y'],q_index_dict['i_z']]
        progress_bar(k+1,len(master_files),new_start_t)


    print('\n Normalising Q Space')

    for s in range(q_final['q_idx'].shape[0]):
        for t in range(q_final['q_idx'].shape[1]):
            for u in range(q_final['q_idx'].shape[2]):
                if q_final['q_idx'][s, t, u] == 0:
                    pass
                else:
                    q_final['q_data'][s, t, u] = q_final['q_data'][s, t, u] / q_final['q_idx'][s, t, u]


    orig_filename = str(scan_num[0])+'_new_3d_fill'

    suffix = '.pickle'
    print(orig_filename)
    export = {  'qx':   q_final['q_x_axis'],\
                'qy':   q_final['q_y_axis'],\
                'qz':   q_final['q_z_axis'],\
                'data': q_final['q_data']}
    new_filename = existential_check(orig_filename,
                                     suffix,
                                     directory + 'processed_files/')

    with open(new_filename,'wb') as handle:
        pickle.dump(export, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print("Complete, total duration: {}".format(int(time.time()
                                                    - start_t))
                                                    + ' seconds')
