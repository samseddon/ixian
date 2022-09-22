import os
import fabio
import time
import numpy as np
from nexusformat import nexus
from useful_functions import progress_bar,existential_check
from Param_dict import param

def index_grid(Param_dict):  # This is the grid that keeps track of how mnay elements are being put into each Q voxel.
    # Identical size to Q volume
    idx_grid = np.zeros((Param_dict['nr_pts_x'], Param_dict['nr_pts_y'], Param_dict['nr_pts_z']))
    return idx_grid

# Makes a list of 2thetas corresponding to each row of pixels (top down) in an image
def generate_two_thetas(TwoTheta0, Param_dict):
    two_thetas = []
    for x in range(Param_dict['ShotHeight']):
        two_theta = TwoTheta0 + (Param_dict['Zero_Pixel_Ver'] - float(x)) * 0.01126  # degrees per row determined from alignment calibration (pixel # vs. two theta)
        two_thetas.append(two_theta)
    return two_thetas

# Calculate qx from 2theta, eta, and wavelength
def calc_qx(twoTheta, eta, wvl, angle_y):
    qx = (2 * np.pi / wvl) * (np.cos(eta / 180 * np.pi) - np.cos(twoTheta / 180 * np.pi - eta / 180 * np.pi)
                              * np.cos(angle_y / 180 * np.pi))
    return qx


# Calculate qz from 2theta, eta, and wavelength
def calc_qz(twoTheta, eta, wvl):
    qz = (2 * np.pi / wvl) * (np.sin(eta / 180 * np.pi) + np.sin(twoTheta / 180 * np.pi - eta / 180 * np.pi))
    return qz


# Calculate qy from 2theta, eta, and wavelength
def calc_qy(twoTheta, eta, wvl, angle_y):
    qy = (2 * np.pi / wvl) * np.sin(angle_y / 180 * np.pi) \
         * np.cos(twoTheta / 180 * np.pi - eta / 180 * np.pi)
    return qy


def qy_angle(offset, Param_dict):
    qy_angles = []
    for x in range(Param_dict['ShotHeight_y']):
        angle = (Param_dict['Zero_Pixel_Hor'] - offset - float(x)) * 0.01126  
        # degrees per column determined from alignment calibration (pixel # vs. nu)
        qy_angles.append(angle)
    return qy_angles


def q_array_3d(Param_dict):

    nr_pts_x = Param_dict['nr_pts_x']
    nr_pts_y = Param_dict['nr_pts_y']
    nr_pts_z = Param_dict['nr_pts_z']
    #if spot == '004':  # Symmetric reflection
    q_3d = np.zeros((nr_pts_x, nr_pts_y, nr_pts_z))

    qz_min = 2.9

    qz_max = 3.1
    delta_qz_range = qz_max - qz_min

    qx_min = -0.05
    qx_max = 0.05

    delta_qx_range = qx_max - qx_min

    qy_min = -0.08

    qy_max = 0.08
    delta_qy_range = qy_max - qy_min

    ###################################################################

    z_slope = nr_pts_z / delta_qz_range  # has intercept
    z_intercept = -z_slope * qz_min

    x_slope = nr_pts_x / delta_qx_range  # has intercept
    x_intercept = -x_slope * qx_min

    y_slope = nr_pts_y / delta_qy_range  # has intercept
    y_intercept = -y_slope * qy_min

    fqz = np.linspace(qz_min, qz_min + delta_qz_range, nr_pts_z)
    fqx = np.linspace(qx_min, qx_min + delta_qx_range, nr_pts_x)
    fqy = np.linspace(qy_min, qy_min + delta_qy_range, nr_pts_y)

    return q_3d, x_slope, y_slope, z_slope, z_intercept, x_intercept, y_intercept, fqx, fqy, fqz




#print(param)
directory="/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/data"
output_folder = "/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/processed_files/"
files_location = os.listdir(directory)
file_reference = "MAG001"
scan_num = [152]
# Selecting the .edf files from a messy folder full of all Data scans, and then with the target scan numbers
master_files = [m for m in files_location if m.startswith(file_reference) and m.endswith(".edf")]
master_files = [c for c in master_files if int(c.split("_")[-2]) in scan_num]
q_3d, x_slope, y_slope, z_slope, z_intercept, x_intercept, y_intercept, qx, qy, qz = q_array_3d(param)
idx_grid = index_grid(param)
offset = 0
start_t = time.time()



for k in range(len(master_files)):
    f = fabio.open(os.path.join(directory, master_files[k]))
    #Stripping the header 
    attenuator = float(f.header.get("counter_pos").split(" ")[24])
    trans = float(f.header.get("counter_pos").split(" ")[23])
    io = float(f.header.get("counter_pos").split(" ")[6])
    two_theta = float(f.header.get("motor_pos").split(" ")[0])  # TWO THETA VALUE
    omega = float(f.header.get("motor_pos").split(" ")[1])  # OMEGA VALUE
    # chi = float(f.header.get("motor_pos").split(" ")[2])  # CHI VALUE
    phi = float(f.header.get("motor_pos").split(" ")[3])  # PHI VALUE
    
    two_theta_range = np.array(generate_two_thetas(two_theta,param))
    two_theta_h_range = np.array(qy_angle(offset,param))

    if phi > - 50: 
        omega_offset = 1.6923 
    else: 
        omega_offset = 3.8765 
    omega = omega + omega_offset 
    
    #  REMOVE THE DEAD PIXELS
    DeadRow1 = param['DeadRow1']   
    # Defines the Dead regions on the Pilatus where the chips are connected - 300k so only in rows!   #######
    DeadRow2 = param['DeadRow2']
    DeadRow3 = param['DeadRow3']
    DeadRow4 = param['DeadRow4']
    wavelength = param['wavelength']
    sat_pix =  param['sat_pix']
    number_x = param['number_x']
    number_y = param['number_y']
    number_z = param['number_z']
    f_new = np.concatenate((f.data[:DeadRow1 - 1, :],
                        f.data[DeadRow2 + 1:DeadRow3 - 1, :],
                        f.data[DeadRow4 + 1:, :]), axis=0)
    two_theta_range_new = np.concatenate((two_theta_range[:DeadRow1 - 1],
                                          two_theta_range[DeadRow2 + 1:DeadRow3 - 1],
                                          two_theta_range[DeadRow4 + 1:]))

    f_new = (f_new / io) * 10 ** 6
    for i in range(f_new.shape[1]):
        for j in range(f_new.shape[0]):  # FIRST ITERATE THE ROWS A.K.A TWO THETA
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

    f.close()

    ##########################################
    progress_bar(k+1,len(master_files),start_t)
    ##########################################


for i in range(idx_grid.shape[0]):
    for j in range(idx_grid.shape[1]):
        for k in range(idx_grid.shape[2]):
            if idx_grid[i, j, k] == 0:
                pass
            else:
                q_3d[i, j, k] = q_3d[i, j, k] / idx_grid[i, j, k]


q_3d = nexus.NXfield(q_3d,'Counts')
q_3d.NXaxes = [qx, qy, qz]
a = nexus.NXdata(signal=q_3d, axes=(qx, qy, qz), axes_names={'x': qx, 'y': qy, 'z': qz})
orig_filename = str(scan_num[0])+'_3d_fill'
suffix = '.nxs'
a.save(existential_check(orig_filename, suffix, output_folder))





