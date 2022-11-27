import os
import fabio
import time
from equations import *
import pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


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
#for k in range(len(master_files)):
    for k in range(1):
        pixel_segmenter(k,directory, master_files, param, start_t,temp,mag)


def param_read(spot_dict,scan_num,directory):
    import_var = directory[:-4] + "user_defined_parameters/param/param_" + spot_dict[str(scan_num)]+'.txt'
    if os.path.exists(import_var)==False:
       raise(ValueError("Param_dict for spot non-existant"))
    with open(import_var,'r') as inf:
        dict1 = eval(inf.read())
    return dict1


def pixel_segmenter(k,directory, master_files,param,start_t,temp,mag):
    '''This function gets given the relevant files, one at time and will eventually return an array of '''
    f1 = fabio.open(os.path.join(directory, master_files[k]))

    attenuator = float(f1.header.get("counter_pos").split(" ")[24])
    trans = float(f1.header.get("counter_pos").split(" ")[23])
    io = float(f1.header.get("counter_pos").split(" ")[6])
    two_theta = float(f1.header.get("motor_pos").split(" ")[0])  # TWO THETA VALUE
    omega = float(f1.header.get("motor_pos").split(" ")[1])  # OMEGA VALUE
    chi = float(f1.header.get("motor_pos").split(" ")[2])  # CHI VALUE
    phi = float(f1.header.get("motor_pos").split(" ")[3])  # PHI VALUE
    temp.append(float(f1.header.get("counter_pos").split(" ")[35]))
    mag.append(float(f1.header.get("motor_pos").split(" ")[88]))
    wavelength = param['wavelength']
    sat_pix =  param['sat_pix']
    number_x = param['nr_pts_x']-1
    number_y = param['nr_pts_y']-1
    number_z = param['nr_pts_z']-1
    two_theta_range = np.array(generate_two_thetas(two_theta,param))
    two_theta_h_range = np.array(qy_angle(param))#

    if phi > - 50:
        omega_offset = 1.6923
    else:
        omega_offset = 3.8765
    omega = omega + omega_offset
    
    f_image,two_theta_range_new = dead_pixel(f1,param,two_theta_range)

    m,n = np.shape(f_image)
#    f_image = []
#    for i in range(n):
#        f_image.append([])
#    for i in range(m):
#        for j in range(n):
    realx_lim_low = 200
    realx_lim_hig = 300
    realy_lim_low = 260
    realy_lim_hig = 330
    f_cut = np.array(f_image)
    f_cut = f_cut[realy_lim_low:realy_lim_hig,realx_lim_low:realx_lim_hig]
    f_cut_qx = np.zeros((np.shape(f_cut)[0],np.shape(f_cut)[1]))
    f_cut_qy = np.zeros((np.shape(f_cut)[0],np.shape(f_cut)[1]))
    f_cut_qz = np.zeros((np.shape(f_cut)[0],np.shape(f_cut)[1]))

    for i in range(np.shape(f_cut)[1]):
        for j in range(np.shape(f_cut)[0]):
            f_cut_qx[j,i] = calc_qx(two_theta_range_new[j], omega, wavelength, two_theta_h_range[i])
            f_cut_qy[j,i] = calc_qy(two_theta_range_new[j], omega, wavelength, two_theta_h_range[i])
            f_cut_qz[j,i] = calc_qz(two_theta_range_new[j], omega, wavelength)
    
#    fig = plt.figure(figsize = (5,5))
#    ax = sns.heatmap(f_cut)
#    plt.show()
    print(f_cut_qx)    
    print(f_cut_qy)    
    print(f_cut_qz)    
directory="/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/data"
output_folder = "/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/processed_files/"
file_reference = "MAG001"
scan_num = [152]

data_fill(directory, output_folder, file_reference, scan_num)

