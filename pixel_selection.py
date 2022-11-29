import os
import fabio
import time
from equations import *
import pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def q_lim(spot_dict,scan_num,directory):
    import_var = directory[:-4] + "user_defined_parameters/qlim/qlim_" + spot_dict[str(scan_num[0])]+'.txt'
    if os.path.exists(import_var)==False:
       raise(ValueError("qlim file does not exist"))
    with open(import_var,'r') as inf:
        dict1 = eval(inf.read())
    return dict1
            


def q_array_init(param, spot_dict,scan_num, directory):
    nr_pts_x = param['nr_pts_x']
    nr_pts_y = param['nr_pts_y']
    nr_pts_z = param['nr_pts_z']

    qlim_dict = q_lim(spot_dict,scan_num,directory)
    
    dict1 = {}
    dict1['q_data'] = np.zeros((nr_pts_x, nr_pts_y, nr_pts_z))
    dict1['q_idx']  = np.zeros((nr_pts_x, nr_pts_y, nr_pts_z))  
    dict1['q_x_axis'] = np.linspace(qlim_dict['qx_min'], qlim_dict['qx_max'], nr_pts_x)
    dict1['q_y_axis'] = np.linspace(qlim_dict['qy_min'], qlim_dict['qy_max'], nr_pts_y)
    dict1['q_z_axis'] = np.linspace(qlim_dict['qz_min'], qlim_dict['qz_max'], nr_pts_z)
    
    return dict1, dict1['q_x_axis'],dict1['q_y_axis'],dict1['q_z_axis'],qlim_dict

def find_q_index(qx,qy,qz,q_x_fin,q_y_fin,q_z_fin,qlim_dict):
    q_index_dict = {}
    
    if qx < qlim_dict['qx_min'] or qx > qlim_dict['qx_max'] \
            or qy < qlim_dict['qy_min'] or qy > qlim_dict['qy_max'] \
            or qz < qlim_dict['qz_min'] or qz > qlim_dict['qz_max']:
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


def data_fill(directory,output_folder,file_reference,scan_num):
    files_location = os.listdir(directory)
    master_files = [m for m in files_location if m.startswith(file_reference) and m.endswith(".edf")]
    master_files = [c for c in master_files if int(c.split("_")[-2]) in scan_num]
    start_t = time.time()
    temp= []
    mag = []
    with open(directory[:-4]+'user_defined_parameters/spot_dict.txt','r') as inf:
        spot_dict = eval(inf.read())
    param = param_read(spot_dict,scan_num[0],directory)
    q_final,q_x_fin,q_y_fin,q_z_fin,qlim_dict = q_array_init(param, spot_dict, scan_num, directory)

    print('\nSlicing images and calulating pixel Q values..')
    q_unsorted = []
    for k in range(len(master_files)):
        q_unsorted.append(pixel_segmenter(k,directory, master_files, param, start_t,temp,mag))
        progress_bar(k+1,len(master_files),start_t)
    new_start_t = time.time()
    print('Inserting pixels into Q_voxels')
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

    for s in range(q_final['q_idx'].shape[0]):
        for t in range(q_final['q_idx'].shape[1]):
            for u in range(q_final['q_idx'].shape[2]):
                if q_final['q_idx'][s, t, u] == 0:
                    pass
                else:
                    q_final['q_data'][s, t, u] = q_final['q_data'][s, t, u] / q_final['q_idx'][s, t, u]


    orig_filename = str(scan_num[0])+'_new_3d_fill'
    suffix = '.pickle'
    export = {  'qx':   q_final['q_x_axis'],\
                'qy':   q_final['q_y_axis'],\
                'qz':   q_final['q_z_axis'],\
                'data': q_final['q_data']}
    
    with open(existential_check(orig_filename, suffix, output_folder), 'wb') as handle:
        pickle.dump(export, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print("Complete, total duration: {}".format(int(time.time() - start_t))+ ' seconds')

def param_read(spot_dict,scan_num,directory):
    import_var = directory[:-4] + "user_defined_parameters/param/param_" + spot_dict[str(scan_num)]+'.txt'
    if os.path.exists(import_var)==False:
       raise(ValueError("Param_dict for spot non-existant"))
    with open(import_var,'r') as inf:
        dict1 = eval(inf.read())
    return dict1


'''This function strips the motor and counter positions, and spits them back as a dictionary 
with the pnemonics '''
def header_strip(f):
    dict_count  = {}
    dict_motor  = {} 
    ub = np.array(f.header.get('UB_pos').split(' '))
    count_mne   = f.header.get('counter_mne').split(' ')
    count_pos   = f.header.get('counter_pos').split(' ')
    motor_mne   = f.header.get('motor_mne').split(' ')
    motor_pos   = f.header.get('motor_pos').split(' ')
    ub          = ub.astype(float) 
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
    return dict_count, dict_motor, ub


'''This function gets given the relevant files, one at time and will eventually return an array of '''
def pixel_segmenter(k,directory, master_files,param,start_t,temp,mag):
    f1 = fabio.open(os.path.join(directory, master_files[k]))
    count_head, motor_head, ub = header_strip(f1) 
 
    attn        =   float(count_head['attn'])
    trans       =   float(count_head['trans'])
    io          =   float(count_head['Io'])
    temp.append(    float(count_head['temp_B']))

    two_theta   =   float(motor_head['del'])
    omega       =   float(motor_head['eta'])
    chi         =   float(motor_head['chi'])   
    phi         =   float(motor_head['phi'])   
    mag.append(     float(motor_head['ami_mag']))

    wavelength  =   param['wavelength']
    number_x    =   param['nr_pts_x']-1
    number_y    =   param['nr_pts_y']-1
    number_z    =   param['nr_pts_z']-1
    
    two_theta_range = np.array(generate_two_thetas(two_theta,param))
    two_theta_h_range = np.array(qy_angle(param))#

    if phi > - 50:
        omega_offset = 1.6923
    else:
        omega_offset = 3.8765
    omega = omega + omega_offset
    
    '''These parameters need to be put into a param_004 file'''
    realx_lim_low = 200
    realx_lim_hig = 300
    realy_lim_low = 260
    realy_lim_hig = 330

    f_cut,two_theta_range_new, two_theta_h_range = dead_pixel(f1,param,two_theta_range,\
            two_theta_h_range, realx_lim_low,realx_lim_hig,realy_lim_low,realy_lim_hig)
    
    #f_cut = np.array(f_image)
    #f_cut = f_cut[realy_lim_low:realy_lim_hig,realx_lim_low:realx_lim_hig]
    array = []
    for i in range(np.shape(f_cut)[1]):
        row = []
        for j in range(np.shape(f_cut)[0]):
            dict1 = {}
            dict1['data']=f_cut[j,i]
            dict1['qx'] = calc_qx(two_theta_range_new[j], omega, wavelength, two_theta_h_range[i])
            dict1['qy'] = calc_qy(two_theta_range_new[j], omega, wavelength, two_theta_h_range[i])
            dict1['qz'] = calc_qz(two_theta_range_new[j], omega, wavelength)
           # print(i,j,dict1['qx'],dicty['qy'])
            row.append(dict1)
        array.append(row)
    #print(array)
#    fig = plt.figure(figsize = (5,5))
#    ax = sns.heatmap(array['data'])
#    plt.show()
    f1.close()
    return array


directory="/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/data"
output_folder = "/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/processed_files/"
file_reference = "MAG001"
scan_num = [152]

data_fill(directory, output_folder, file_reference, scan_num)
#q_check(directory,scan_num)
