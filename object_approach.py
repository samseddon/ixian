import os
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
"""
Created on Wed Jan 11 14:49:21 2023

@author: samseddon
"""

def q_lim(spot_dict,
          scan_num,
          directory):
    """ Uses scan number to read a file to open relevant scan dictionary 
    and returns a dictionary of the Q_space limits
           Parameters:
               spot_dict(dict)   : dictionary of scans and spots
               scan_num(array)   : array of scan numbers being read
               directory(string) : string pointing to data directory structure
           Returns
               dict1(dict)       : dictionary of q_limits and n_pts
    """
    import_var = directory \
                 + "user_defined_parameters/qlim/qlim_" \
                 + spot_dict[str(scan_num[0])] \
                 + ".txt" \

    if os.path.exists(import_var)==False:
       raise(ValueError("qlim file does not exist, consider running"\
                         "again with create_files = True"))
    
    with open(import_var,"r") as inf:
        dict1 = eval(inf.read())
    return dict1
            

def q_array_init(param,
                 spot_dict,
                 scan_num,
                 directory):
    '''Here the Q array is created from the dictionary of Q_limits, sliced
    into a number of points n_pts. 
           Parameters: 
               param(dict)       : dictionary of general spot parameters 
               spot_dict(dict)   : dictionary of scans and spots
               scan_num(array)   : array of scan numbers being read
               directory(string) : string pointing to data directory structure
           Returns:
               dict1(dict)       : dictionary containing a two arrrays, q_data
                                   to be populated with pixels, q_idx to keep
                                   track of how many data points have been put
                                   into q_data, and the q axes of these arrays
               3 x q-axes(array) : required later so exported simultaneously
                                   to save time
               q_lim_dict(dict)  : dictionary of qlims and number of points q
                                   space was split into
    '''
    qlim_dict = q_lim(spot_dict,scan_num,directory)
    
    dict1 = {}
    dict1["q_data"] = np.zeros((qlim_dict["nr_pts"],
                                qlim_dict["nr_pts"],
                                qlim_dict["nr_pts"]))
    
    dict1['q_idx']  = np.zeros((qlim_dict["nr_pts"],
                                qlim_dict["nr_pts"],
                                qlim_dict["nr_pts"]))  
    
    dict1["q_x_axis"] = np.linspace(qlim_dict["qx_min"],
                                    qlim_dict["qx_max"], 
                                    qlim_dict["nr_pts"])
    
    dict1["q_y_axis"] = np.linspace(qlim_dict["qy_min"], 
                                    qlim_dict["qy_max"], 
                                    qlim_dict["nr_pts"])
    
    dict1["q_z_axis"] = np.linspace(qlim_dict["qz_min"], 
                                    qlim_dict["qz_max"], 
                                    qlim_dict["nr_pts"])
    
    return dict1,\
           dict1["q_x_axis"],\
           dict1["q_y_axis"],\
           dict1["q_z_axis"],\
           qlim_dict


def find_q_index(qx,
                 qy,
                 qz,
                 q_x_fin,
                 q_y_fin,
                 q_z_fin,
                 qlim_dict):
    """ Checks if the given Q coordinates for a pixel in q-space are in 
    the new q-space array, and returns the index of the points in the 
    q_array ready for population of that space
           Parameters: 
               qx(int)           : Pixel qx coordinate
               qy(int)           : Pixel qy coordinate
               qz(int)           : Pixel qz coordinate
               q_x_fin(array)    : qx coordinates of initilised Q array 
               q_y_fin(array)    : qy coordinates of initilised Q array
               q_z_fin(array)    : qz coordinates of initilised Q array
               q_lim_dict(dict)  : dictionary of qlims and number of points q
                                   space was split into
           Returns:
               1             : When pixel is out of range
               q_index_dict(dict) : dictionary of indexes of pixels position
                                    in the initilised Q_space
    """
    q_index_dict = {}
    if qx < qlim_dict['qx_min'] \
       or qx > qlim_dict['qx_max'] \
       or qy < qlim_dict['qy_min'] \
       or qy > qlim_dict['qy_max'] \
       or qz < qlim_dict['qz_min'] \
       or qz > qlim_dict['qz_max']:    
        return(1)
    else: 
        for i_x in range(len(q_x_fin)):
            if(qx > q_x_fin[i_x]):
                 pass
            else:        
                q_index_dict['i_x'] = i_x
                break

        for i_y in range(len(q_y_fin)):
            if(qy > q_y_fin[i_y]):
                pass
            else:        
                q_index_dict['i_y'] = i_y
                break
        
        for i_z in range(len(q_z_fin)):
            if(qz > q_z_fin[i_z]):
                pass
            else:        
                q_index_dict['i_z'] = i_z
                break
        return q_index_dict


def find_n_lim(directory,spot_dict,scan_num,limit_dict,qlim_dict):
    """ This program takes the qlimits, and the average pixel size in q-space,
    and maximises n_pts, the number of points our new q_array is split up 
    into, thus achieving the best possible resolution without underpopulating 
    the new q_array.
           Parameters: 
               directory(string) : string pointing to data directory structure
               spot_dict(dict)   : dictionary of scans and spots
               scan_num(array)   : array of scan numbers being read
               limit_dict        : dictionary of size of average size of
                                   pixels in qspace to avoid striping
               qlim_dict(dict)   : dictionary of qlims 
                                    
           Returns:
               qlim_dict(dict)  : dictionary of qlims with the newly appended
                                  n_pts, so q_array can be properly initilised
    """
    q_x_lim_min = []
    q_y_lim_min = []
    q_z_lim_min = []
    
    for k in range(len(limit_dict)):
        q_y_lim_min.append(limit_dict[k]['pixel_qy'])
        q_x_lim_min.append(limit_dict[k]['pixel_qx'])
        q_z_lim_min.append(limit_dict[k]['pixel_qz'])
    
    q_x_lim_min = np.average(q_x_lim_min)
    q_z_lim_min = np.average(q_z_lim_min)
    q_y_lim_min = abs(max(q_y_lim_min)-min(q_y_lim_min))/len(q_y_lim_min)
    
    nr_pts = 1
    while abs(qlim_dict['qx_max']-qlim_dict['qx_min'])/nr_pts > q_x_lim_min \
      and abs(qlim_dict['qy_max']-qlim_dict['qy_min'])/nr_pts > q_y_lim_min \
      and abs(qlim_dict['qz_max']-qlim_dict['qz_min'])/nr_pts > q_z_lim_min:
              nr_pts = nr_pts + 1 
    
    qlim_dict['nr_pts'] = nr_pts
    
    return qlim_dict


def find_q_lim(q,directory,spot_dict,scan_num,limit_dict):
    """ Given all of the data from the pixels, the minimum and maximum of 
    their q-space coordinates are found and the absolute limits of the 
    detector in q-space is determined, to create a q_array of perfect size. 
    It passes these to function find_n_lim, and writes the qlim dictionary 
    for future reference.

            Parameters: 
               q(dict(array))    : all pixel data across k scans
               directory(string) : string pointing to data directory structure
               spot_dict(dict)   : dictionary of scans and spots
               scan_num(array)   : array of scan numbers being read
               limit_dict        : dictionary of size of average size of
                                   pixels in qspace to avoid striping
    """
    qx = []
    qy = []
    qz = []
    
    for k in range(len(q)):
        for i in range(np.shape(q[k])[1]):
            for j in range(np.shape(q[k])[0]):
                qx.append(q[k][j][i]['qx'])
                qy.append(q[k][j][i]['qy'])
                qz.append(q[k][j][i]['qz'])
    
    q_limits_without_n_pts = {'qz_min' : min(qz),
                              'qz_max' : max(qz),
                              'qx_min' : min(qx),
                              'qx_max' : max(qx),
                              'qy_min' : min(qy),
                              'qy_max' : max(qy)}
    
    q_limits_with_n_pts = find_n_lim(directory,
                        spot_dict,
                        scan_num,
                        limit_dict,
                        q_limits_without_n_pts)
    filename = directory\
               + 'user_defined_parameters/qlim/qlim_'\
               + spot_dict[str(scan_num[0])]\
               + '.txt'
    
    if os.path.exists(filename) == True \
        and input('Overwrite existing file, [y] or n?\n') != 'y':
            pass
    else:
       with open(filename,'w') as inf:
          inf.write(str(q_limits_with_n_pts))
       print('Created file',filename)


def data_fill(directory,output_folder,file_reference,scan_num,create_files):
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
    for image_number in range(len(master_files)):
        all_images.append(Dectris_Image(file_reference, 
                                        image_number, 
                                        directory, 
                                        master_files, 
                                        param))
        #fig, ax = plt.subplots()
        #ax.contourf(all_images[-1].data,cmap = cbar)
        #plt.show()
        # print(all_images[-1].Q_coords)
#        q_unsorted_temp,limit_dict_temp = pixel_segmenter(image_number,
#                                                          directory, 
#                                                          master_files, 
#                                                          param, 
#                                                          start_t,
#                                                          temp,
#                                                          mag,
#                                                          file_reference)
#        q_unsorted.append(q_unsorted_temp)
#        limit_dict.append(limit_dict_temp)
        #print(limit_dict_temp['pixel_qz'])
        progress_bar(image_number+1,len(master_files),start_t)
    
    print('\nFinding q limits from sliced data and optimising Q_space mesh')
    qx_min = []
    qy_min = []
    qz_min = []
    qx_max = []
    qy_max = []
    qz_max = []
    q_min = []
    for dec_image in all_images:
        qx, qy, qz = dec_image.q_lim()
        qx_min.append(qx[0])
        qy_min.append(qy[0])
        qz_min.append(qz[0])
        qx_max.append(qx[1])
        qy_max.append(qy[1])
        qz_max.append(qz[1])
    
    # NOTE here the qmax/qlim are found, needs to be put into a function.
    print(np.amax(qx_max), np.amax(qy_max), np.amax(qz_max))
    print(np.amin(qx_min), np.amin(qy_min), np.amin(qz_min))
#    print(all_images[-1].file_reference)
    #print(q_min[:][:][0])

    if create_files == True:
        find_q_lim(q_unsorted,directory,spot_dict,scan_num,limit_dict)

    q_space = Q_Space(scan_num, spot_dict, directory)
#    q_final,q_x_fin,q_y_fin,q_z_fin,qlim_dict = q_array_init(param, 
#                                                             spot_dict, 
#                                                             scan_num, 
#                                                             directory)
    new_start_t = time.time()
    print('Populating Q_space with pixels')
#    for i in range(len(q_x_fin)):
#        print(q_x_fin[i], q_space.q_x[i]) 

    for _ in range(len(all_images)):
        q_space.populate_3D(all_images[_])
        progress_bar(_ + 1,len(master_files),new_start_t)
    q_space.normalise_3D()
    """    
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


    print('\n Normalising Q Space')
    
    for s in range(q_final['q_idx'].shape[0]):
        for t in range(q_final['q_idx'].shape[1]):
            for u in range(q_final['q_idx'].shape[2]):
                if q_final['q_idx'][s, t, u] == 0:
                    pass
                else:
                    q_final['q_data'][s, t, u] = q_final['q_data'][s, t, u] / q_final['q_idx'][s, t, u]
    """
    orig_filename = str(scan_num[0])+'_new_3d_fill'
    suffix = '.pickle'
    print(orig_filename) 
#    export = {  'qx':   q_final['q_x_axis'],\
#                'qy':   q_final['q_y_axis'],\
#                'qz':   q_final['q_z_axis'],\
#                'data': q_final['q_data']}
    new_filename = existential_check(orig_filename,
                                     suffix, 
                                     directory + 'processed_files/')
                                     
    with open(new_filename,'wb') as handle:
        pickle.dump(q_space, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    print("Complete, total duration: {}".format(int(time.time() 
                                                    - start_t))
                                                    + ' seconds')

def param_read(spot_dict,scan_num,directory):
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


def header_strip(f): # Incorporated into INIT 
    """ Takes each file when called and returns the header motors, counter and
    ub matrix for each scan. 
           Parameters: 
               f(fabio.edfimage)      : raw image file
           Returns:
               dict_count(dict)       : dictionary of counter paramters 
                                        with the keys as their mnemonics
               dict_motor(dict)       : dictionary of motor paramters 
                                        with the keys as their mnemonics
               ub(array)              : ub matrix outputted as an array
               WAVELENGTH             : WAVELENGTH
    """
    dict_count  = {}
    dict_motor  = {} 

    ub = np.array(f.header.get('UB_pos').split(' '))
    count_mne   = f.header.get('counter_mne').split(' ')
    count_pos   = f.header.get('counter_pos').split(' ')
    motor_mne   = f.header.get('motor_mne').split(' ')
    motor_pos   = f.header.get('motor_pos').split(' ')
    ub          = ub.astype(float) 
    WAVELENGTH  = f.header.get('source_wavelength').split(' ')
    WAVELENGTH  = float(WAVELENGTH[0])
    for key in count_mne:
        for value in count_pos:
            dict_count[key] = value
            count_pos.remove(value)
            break
    
    for key in motor_mne:
        for value in motor_pos:
            dict_motor[key] = value
            motor_pos.remove(value)
            break
    
    return dict_count, dict_motor, ub, WAVELENGTH


def omega_offsetter(omega,
                    file_reference,
                    phi):
    """ At least temporarily removed the presence of magic numbers in the code
    which with a bit of magic from didier and I later in the experiment should
    result in not needing omega offsets at all if we just use the ub (or 
    something) 
           Parameters: 
               omega(float)           : omega values of the chosen scan image
               file_reference(string) : unique file identifier for a given     
                                        experiment
               phi(float)             : phi value of the chosen scan image
           Returns:
               float                  : omega+omega_offset, an outdated and 
                                        soon fixed
    """
    if file_reference == 'MAG001':
        if phi > - 50:
            omega_offset = 1.6923
        else:
            omega_offset = 3.8765
    else: 
        omega_offset = 0
    return omega + omega_offset


def parameter_setup(directory,
                    master_files,
                    file_reference,
                    spot_dict,
                    scan_num):
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
    if os.path.exists(filename) == True \
        and input('Overwrite existing '
                  + spot_dict[str(scan_num[0])]
                  + ' file, [y] or n?\n')\
            != 'y':
        pass
    else:
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
        
        fig = plt.figure(figsize = (12,3))
        fig.add_subplot(131)
        
        ax = sns.heatmap(f1)
        ax.add_patch(patches.Rectangle(
            (REALX_LIM_LOW, REALY_LIM_LOW),
            REALX_LIM_HIG-REALX_LIM_LOW,
            REALX_LIM_HIG-REALX_LIM_LOW,
            edgecolor='red',
            fill=False,
            lw=2))
        fig.add_subplot(132)
        ax = sns.heatmap(f2)
        ax.add_patch(patches.Rectangle(
            (REALX_LIM_LOW, REALY_LIM_LOW),
            REALX_LIM_HIG-REALX_LIM_LOW,
            REALX_LIM_HIG-REALX_LIM_LOW,
            edgecolor='red',
            fill=False,
            lw=2))
        fig.add_subplot(133)
        ax = sns.heatmap(f3)
        ax.add_patch(patches.Rectangle(
            (REALX_LIM_LOW, REALY_LIM_LOW),
            REALX_LIM_HIG-REALX_LIM_LOW,
            REALX_LIM_HIG-REALX_LIM_LOW,
            edgecolor='red',
            fill=False,
            lw=2))
        plt.show()

        with open(filename,'w') as inf:
            inf.write(str(gen_param))
        print('Created file',
              filename,
              '\n Change slice region in parameter file')



def pixel_segmenter(k,
                    directory, 
                    master_files,
                    param,
                    start_t,
                    temp,
                    mag,
                    file_reference):
    '''This function gets given the relevant files, one at time, opens them 
    and strips the header, it slices out  the data as defined in parameter 
    files, and outputs back an array item in scan numbers  of dictionaries 
    of the pixels data points and their q_x,q_y and q_z coordinates. It also 
    returns any relevant motor/counter positions for any scans that aren't 
    just an RSM/eta scan, such as an ascan field or temperature
           Parameters: 
               k(int)                 : looping variable function is looped 
                                        over corresponding to each scan
               directory(string)      : string pointing to data directory 
               master_files(list)     : list of files for given scan number
               param(dict)            : dictionary of general spot parameters 
               start_t(float)         : time function first called for 
                                        progress bar
               temp(list)             : list of temperature of scan
               mag(list)              : list of magnetic field of scan
               file_reference(string) : unique file identifier for a given 
           Returns:
               dict_of_data_and_coordinates : list of dictionaries containing 
                                              pixel data and q space 
                                              coordinate
               limit_dict                   : dictionary of size of average 
                                              size of pixels in qspace to 
                                              avoid striping
    '''
    f = fabio.open(os.path.join(directory+'data/', master_files[k]))

    count_head, motor_head, ub, wavelength = header_strip(f) 
 
    attn        =   float(count_head['attn'])
    trans       =   float(count_head['trans'])
    io          =   float(count_head['Io'])
    temp.append(    float(count_head['temp_B']))

    two_theta   =   float(motor_head['del'])
    omega       =   float(motor_head['eta'])
    chi         =   float(motor_head['chi'])   
    phi         =   float(motor_head['phi'])   
   
    if motor_head['ami_mag'] in motor_head:
        mag.append(float(motor_head['ami_mag']))
        
    two_theta_range = np.array(generate_two_thetas(two_theta,param))
    two_theta_h_range = np.array(qy_angle(param))#
    
    omega = omega_offsetter(omega,file_reference,phi)
    
    realx_lim_low = param['realx_lim_low']
    realx_lim_hig = param['realx_lim_hig']
    realy_lim_low = param['realy_lim_low']
    realy_lim_hig = param['realy_lim_hig']

    f_cut,\
    two_theta_range_new,\
    two_theta_h_range = dead_pixel(f,
                                   param,
                                   two_theta_range,
                                   two_theta_h_range,
                                   realx_lim_low,
                                   realx_lim_hig,
                                   realy_lim_low,
                                   realy_lim_hig)
    
    f.close()
    dict_of_data_and_coordinates = []
    qz = []
    qx = []
    qy = []
    for i in range(np.shape(f_cut)[1]):
        row = []
        qz_row = []
        qy_row = []
        qx_row = []
        for j in range(np.shape(f_cut)[0]):
            dict1 = {}
            dict1['data']=f_cut[j,i]
            dict1['qx'] = calc_qx(two_theta_range_new[j], 
                                  omega,
                                  wavelength,
                                  two_theta_h_range[i])
            dict1['qy'] = calc_qy(two_theta_range_new[j],
                                  omega, 
                                  wavelength, 
                                  two_theta_h_range[i])
            dict1['qz'] = calc_qz(two_theta_range_new[j],
                                  omega,
                                  wavelength)

            row.append(dict1)
            qz_row.append(dict1['qz'])
            qx_row.append(dict1['qx'])
            qy_row.append(dict1['qy'])

        dict_of_data_and_coordinates.append(row)
        qz.append(abs(max(qz_row)))
        qx.append(abs(max(qx_row)))
        qy.append(np.average(qy_row))

    limit_dict = {}

    limit_dict['pixel_qy'] = np.average(qy)
    limit_dict['pixel_qz'] = qz
    limit_dict['pixel_qx'] = np.average(qx)
    
    
    return dict_of_data_and_coordinates, limit_dict


directory="/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/"
output_folder = "/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/processed_files/"
file_reference = "MAG001"
scan_num = [180]

#data_fill(directory, output_folder, file_reference, scan_num, create_files = False)
#q_check(directory,scan_num)
