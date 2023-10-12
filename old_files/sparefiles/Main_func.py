import os
import fabio
import time
import numpy as np
from nexusformat import nexus
from useful_functions import progress_bar,existential_check
#from Param_dict import param
from equations import *
from spot_dict import spot_dict
from joblib import Parallel, delayed
import multiprocessing as mp 
import pickle

def sanity_check(f_new,omega, wavelength, two_theta_h_range, two_theta_range_new,
        i,j,attenuator,x_slope, x_intercept,y_slope,y_intercept,z_slope,
        z_intercept, number_x,number_y, number_z,io,sat_pix,average_qz):
    average_qz.append(calc_qz(two_theta_range_new[j], omega, wavelength)-calc_qz(two_theta_range_new[j-1], omega, wavelength))
    return average_qz

def voxel_fill(f_new,omega, wavelength, two_theta_h_range, two_theta_range_new, 
        i,j,attenuator,x_slope, x_intercept,y_slope,y_intercept,z_slope,
        z_intercept, q_3d, idx_grid, number_x,number_y, number_z,io,sat_pix):
    fqx = calc_qx(two_theta_range_new[j], omega, wavelength, two_theta_h_range[i])
    fqy = calc_qy(two_theta_range_new[j], omega, wavelength, two_theta_h_range[i])
    fqz = calc_qz(two_theta_range_new[j], omega, wavelength)
    if attenuator == 0:

        if int(fqx * x_slope + x_intercept) > number_x \
                or int(fqy * y_slope + y_intercept) > number_y \
                or int(fqz * z_slope + z_intercept) > number_z \
                or int(fqx * x_slope + x_intercept) < 0 \
                or int(fqy * y_slope + y_intercept) < 0 \
                or int(fqz * z_slope + z_intercept) < 0:

            # skipped = skipped+1
            pass
        else:
            if f_new[j, i] * io * 10 ** -6 > (sat_pix):  # 1.08*10**6
                print('hotspot')
                # hotspot=hotspot+1
                pass
            else:
                q_3d[int(fqx * x_slope + x_intercept),
                     int(fqy * y_slope + y_intercept),
                     int(fqz * z_slope + z_intercept)] \
                    = q_3d[int(fqx * x_slope + x_intercept),
                           int(fqy * y_slope + y_intercept),
                           int(fqz * z_slope + z_intercept)] + f_new[j, i]
                idx_grid[int(fqx * x_slope + x_intercept),
                         int(fqy * y_slope + y_intercept),
                         int(fqz * z_slope + z_intercept)] \
                    = idx_grid[int(fqx * x_slope + x_intercept),
                               int(fqy * y_slope + y_intercept),
                               int(fqz * z_slope + z_intercept)] + 1


    if attenuator != 0 and f_new[j, i] < 1000:  # T.H. Made 1000 from 100
        # low_count = low_count+1
        pass
    if attenuator != 0 and f_new[j, i] > 1000:
        if f_new[j, i] * io * 10 ** -6 > (sat_pix):
            print('hotspot')
            # hotspot= hotspot + 1
            pass
        else:
            if int(fqx * x_slope + x_intercept) > number_x \
                    or int(fqy * y_slope + y_intercept) > number_y \
                    or int(fqz * z_slope + z_intercept) > number_z \
                    or int(fqx * x_slope + x_intercept) < 0 \
                    or int(fqy * y_slope + y_intercept) < 0 \
                    or int(fqz * z_slope + z_intercept) < 0:

                # skipped = skipped + 1
                pass
            else:
                q_3d[int(fqx * x_slope + x_intercept),
                     int(fqy * y_slope + y_intercept),
                     int(fqz * z_slope + z_intercept)] \
                    = q_3d[int(fqx * x_slope + x_intercept),
                           int(fqy * y_slope + y_intercept),
                           int(fqz * z_slope + z_intercept)] + (f_new[j, i] / trans)
                idx_grid[int(fqx * x_slope + x_intercept),
                         int(fqy * y_slope + y_intercept),
                         int(fqz * z_slope + z_intercept)] \
                    = idx_grid[int(fqx * x_slope + x_intercept),
                               int(fqy * y_slope + y_intercept),
                               int(fqz * z_slope + z_intercept)] + 1

    return q_3d, idx_grid,i,j

def q_lim(spot_dict,scan_num):
    import_var = "qlim/var_" + spot_dict[str(scan_num)]+'.txt'
    if os.path.exists(import_var)==False:
       raise(ValueError("you fucked it son")) 
    with open(import_var,'r') as inf:
        dict1 = eval(inf.read())
    delta_qx_range = dict1['qx_max'] - dict1['qx_min']
    delta_qy_range = dict1['qy_max'] - dict1['qy_min']
    delta_qz_range = dict1['qz_max'] - dict1['qz_min']
    return dict1['qz_min'],dict1['qz_max'],delta_qz_range,dict1['qx_min'],\
            dict1['qx_max'],delta_qx_range,dict1['qy_min'],dict1['qy_max'],delta_qy_range
    
def param_read(spot_dict,scan_num):
    import_var = "param/param_" + spot_dict[str(scan_num)]+'.txt'
    if os.path.exists(import_var)==False:
       raise(ValueError("Param_dict for spot non existant"))
    with open(import_var,'r') as inf:
        dict1 = eval(inf.read())
    return dict1




def q_array_3d(param_dict,spot_dict,scan_num):
    nr_pts_x = param_dict['nr_pts_x']
    nr_pts_y = param_dict['nr_pts_y']
    nr_pts_z = param_dict['nr_pts_z']

    q_3d = np.zeros((nr_pts_x, nr_pts_y, nr_pts_z))
    ###################################################################
    qz_min,qz_max,delta_qz_range,qx_min,qx_max,delta_qx_range,qy_min,qy_max,delta_qy_range=q_lim(spot_dict,scan_num)
    z_slope = nr_pts_z / delta_qz_range  # has intercept
    z_intercept = -z_slope * qz_min

    x_slope = nr_pts_x / delta_qx_range  # has intercept
    x_intercept = -x_slope * qx_min

    y_slope = nr_pts_y / delta_qy_range  # has intercept
    y_intercept = -y_slope * qy_min
    fqz = np.linspace(qz_min, qz_min + delta_qz_range, nr_pts_z)
    #print(fqz)
    fqx = np.linspace(qx_min, qx_min + delta_qx_range, nr_pts_x)
    fqy = np.linspace(qy_min, qy_min + delta_qy_range, nr_pts_y)

    return q_3d, x_slope, y_slope, z_slope, z_intercept, x_intercept, y_intercept, fqx, fqy, fqz

def ixian(k,directory, master_files,param,q_3d, x_slope, y_slope, 
        z_slope, z_intercept, x_intercept, y_intercept, idx_grid,average_qz,start_t):
    f1 = fabio.open(os.path.join(directory, master_files[k]))
    f2 = f1
    #Stripping the header 
    attenuator = float(f2.header.get("counter_pos").split(" ")[24])
    trans = float(f2.header.get("counter_pos").split(" ")[23])
    io = float(f2.header.get("counter_pos").split(" ")[6])
    two_theta = float(f2.header.get("motor_pos").split(" ")[0])  # TWO THETA VALUE
    omega = float(f2.header.get("motor_pos").split(" ")[1])  # OMEGA VALUE
    # chi = float(f.header.get("motor_pos").split(" ")[2])  # CHI VALUE
    phi = float(f2.header.get("motor_pos").split(" ")[3])  # PHI VALUE
    wavelength = param['wavelength']
    sat_pix =  param['sat_pix']
    number_x = param['nr_pts_x']-1
    number_y = param['nr_pts_y']-1
    number_z = param['nr_pts_z']-1
    two_theta_range = np.array(generate_two_thetas(two_theta,param))
    two_theta_h_range = np.array(qy_angle(param))

    if phi > - 50: 
        omega_offset = 1.6923 
    else: 
        omega_offset = 3.8765 
    omega = omega + omega_offset 
    f2,two_theta_range_new = dead_pixel(f2,param,two_theta_range)


    for i in range(f2.shape[1]):
        for j in range(f2.shape[0]):  # FIRST ITERATE THE ROWS A.K.A TWO THETA
            q_3d,idx_grid,i,j = voxel_fill(f2,omega, wavelength, two_theta_h_range,
                    two_theta_range_new,i,j,attenuator,x_slope, 
                    x_intercept,y_slope,y_intercept,z_slope,z_intercept,q_3d, idx_grid,number_x,number_y,number_z,io,sat_pix)
#            average_qz = sanity_check(f2,omega, wavelength, two_theta_h_range,
#                    two_theta_range_new,i,j,attenuator,x_slope,
#                    x_intercept,y_slope,y_intercept,z_slope,z_intercept,number_x,number_y,number_z,io,sat_pix,
#                    average_qz)
    f1.close()

    ##########################################
    progress_bar(k+1,len(master_files),start_t)
    ##########################################
    return idx_grid, q_3d, average_qz

def data_fill(directory,output_folder,file_reference,scan_num):
    files_location = os.listdir(directory)
    master_files = [m for m in files_location if m.startswith(file_reference) and m.endswith(".edf")]
    master_files = [c for c in master_files if int(c.split("_")[-2]) in scan_num]
    param = param_read(spot_dict,scan_num[0])
    q_3d, x_slope, y_slope, z_slope, z_intercept, x_intercept, y_intercept, qx, qy, qz = q_array_3d(param,
            spot_dict,scan_num[0])
    idx_grid = index_grid(param)
    offset = 0
    start_t = time.time()
    average_qz = []
    
    for k in range(len(master_files)):
        idx_grid, q_3d, average_qz = ixian(k,directory, master_files,param,q_3d, x_slope, y_slope,
                z_slope, z_intercept, x_intercept, y_intercept,idx_grid,average_qz,start_t)
    
    for i in range(idx_grid.shape[0]):
        for j in range(idx_grid.shape[1]):
            for k in range(idx_grid.shape[2]):
                if idx_grid[i, j, k] == 0:
                    pass
                else:
                    q_3d[i, j, k] = q_3d[i, j, k] / idx_grid[i, j, k]
    
    
    #q_3d = nexus.NXfield(q_3d,'Counts')
    #q_3d.NXaxes = [qx, qy, qz]
    #a = nexus.NXdata(signal=q_3d, axes=(qx, qy, qz), axes_names={'x': qx, 'y': qy, 'z': qz})
    orig_filename = str(scan_num[0])+'_3d_fill'
    #suffix = '.nxs'
    #a.save(existential_check(orig_filename, suffix, output_folder))
    suffix = '.pickle'
    axis_export = []
    axis_export.append(qx)
    axis_export.append(qy)
    axis_export.append(qz)
    with open(existential_check(orig_filename, suffix, output_folder), 'wb') as handle:
        pickle.dump(q_3d, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open(existential_check(orig_filename, suffix, output_folder), 'wb') as handle:
        pickle.dump(q_3d, handle, protocol=pickle.HIGHEST_PROTOCOL)

directory="/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/data"
output_folder = "/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/processed_files/"
file_reference = "MAG001"
scan_num = [180]

data_fill(directory,output_folder,file_reference,scan_num)

