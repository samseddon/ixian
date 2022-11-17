import pickle   
import os
import glob
import numpy as np
from nexusformat import nexus
import matplotlib.pyplot as plt
from scipy import interpolate


with open('182_3d_fill_5.pickle', 'rb') as handle:
   file = pickle.load(handle)
print(file['qx'])



def file_checker(s_num, s_ind, input_path):
    print(os.path.join(input_path, str(s_num) + '*'))
    files = glob.glob(os.path.join(input_path, str(s_num) + '*'))
    print(files)
    file_name = max(files, key=os.path.getctime)
    file_name = file_name[len(input_path):]
    if s_ind == -1:
        with open('182_3d_fill_5.pickle', 'rb') as handle:
             read_file = pickle.load(handle)

#    elif s_ind == 0:
#        short_file_name = min(files, key=os.path.getctime)
#        file_name = short_file_name[len(input_path):]
#        read_file = nexus.nxload(os.path.join(input_path, file_name))
#    else:
#        short_file_name = min(files, key=os.path.getctime)
#        file_name = short_file_name[len(input_path):-4] + '_' + str(s_ind) + '.nxs'
#        read_file = nexus.nxload(os.path.join(input_path, file_name))
    return read_file, file_name

def existential_check(o_f, f_s, folder): #file name, file suffix, destination folder
    f_n = str(o_f + f_s)
    while os.path.exists(folder + f_n):
        if len(f_n) > len(o_f + f_s):
            filenum = int(f_n[(len(o_f + f_s) - (len(f_s)-1)):-len(f_s)]) + 1
            f_n = f_n[:(len(o_f + f_s) - (len(f_s)-1))] + str(filenum) + f_n[-len(f_s):]
        else:
            f_n = f_n[:-len(f_s)] + '_1' + f_n[-len(f_s):]
    return folder + f_n

def input_check(q_1, q_2):
    axis_dict = {'x': 0, 'X': 0, 'y': 1, 'Y': 1, 'z': 2, 'Z': 2}
    if q_1 not in axis_dict or q_2 not in axis_dict:
        print('Axis not recognised, please use x,y or z')
        return
    if axis_dict[q_1] + axis_dict[q_2] == 1:  # int z
        integrate_axis = 2
        print('Integrating along z axis')
    if axis_dict[q_1] + axis_dict[q_2] == 2:  # int y
        integrate_axis = 1
        print('Integrating along y axis')
    if axis_dict[q_1] + axis_dict[q_2] == 3:  # int x
        integrate_axis = 0
        print('Integrating along x axis')
    return axis_dict[q_1], axis_dict[q_2], integrate_axis

def rsm_plot(f_1, file_name, q_1, q_2, q_3_lim, output_folder):
    q_1_assigned, q_2_assigned, q_3_assigned = input_check(q_1, q_2)
    values_1 = f_1['qx']
    values_2 = f_1['qz']

    volume = np.array(f_1['data'])
    print(volume)
    q1_q2_matrix = np.zeros((volume.shape[0], volume.shape[2])) # THIS IS ALWAYS 0 and 2
    print(volume.shape[2])
    if q_3_assigned == 0:
        for i in range(volume.shape[2]):  # this axis 0 for x
            q1_q2_matrix[:, i] = np.copy(volume[q_3_lim[0]:q_3_lim[1],:, i].mean(axis=0) + 10 ** (-8))
            #print(q1_q2_matrix)
    elif q_3_assigned == 1: # int y
        for i in range(volume.shape[2]):  # this axis 0 for y
            q1_q2_matrix[:, i] = np.copy(volume[:, q_3_lim[0]:q_3_lim[1], i].mean(axis=1) + 10 ** (-8))
            #print(q1_q2_matrix)
    elif q_3_assigned == 2: #int z

        for i in range(volume.shape[2]):  # this axis 0 for z
            #print(q1_q2_matrix)
            q1_q2_matrix[:, i] = np.copy(volume[:, i, q_3_lim[0]:q_3_lim[1]].mean(axis=1) + 10 ** (-8))

    fin = q1_q2_matrix.flatten('F')

    fq_1 = []
    fq_2 = []
    # INTERPOLATE TO MAKE A NEW GRID
    # Option to Interpolate to a different mesh size

    for i in range(len(values_2)):
        fq_2.append(values_1)

    fq_2 = np.array(fq_2).flatten()  # this is y in y

    c = 0  # incrementing variable
    nr = 0  # incrementing variable
    for j in range(len(values_2)):
        fq_1.append([values_2[c]] * len(values_1))
        c = c + 1

    fq_1 = np.array(fq_1).flatten()

    for i in range(len(fin)):
        if fin[i] < 0:
            fin[i] = 0.000000001

    q1_q2_matrix = np.reshape(fin, (volume.shape[0], volume.shape[2]))

    mesh_q_1_min = np.min(fq_2)  # sets the minimum qx value for new mesh (np.min(fq_2) is the minimum of the input data)
    mesh_q_1_max = np.max(fq_2)  # sets the maximum qx value for new mesh (np.max(fq_2) is the maximum of the input data)
    mesh_q_2_min = np.min(fq_1)
    mesh_q_2_max = np.max(fq_1)
    nr_fqz = int(len(np.unique(fq_1)) * 1)  # sets the number of mesh points in qx  *1.0 keeps as original
    nr_fqx = int(len(np.unique(fq_2)) * 1)  # sets the number of mesh points in qz  *1.0 keeps as original

    # this is y for y, #grid qx
    grid_qx, grid_qz = np.mgrid[mesh_q_1_min:mesh_q_1_max:(nr_fqx * 1j), mesh_q_2_min:mesh_q_2_max:(nr_fqz * 1j)]
    grid_q = interpolate.griddata((fq_2, fq_1), fin, (grid_qx, grid_qz), method='nearest')

    #  GENERATE LINE SCANS BY INTEGRATING FROM THE NEW Q GRID

    n_qx = np.linspace(mesh_q_1_min, mesh_q_1_max, nr_fqx)  # Determining the areal size of the new mesh
    n_qz = np.linspace(mesh_q_2_min, mesh_q_2_max, nr_fqz)

    #  FIGURE
    fig, ax = plt.subplots(1, 1, figsize=(1, 1), dpi=500)
    cb = img = ax.imshow(np.log10(grid_q), origin='lower', aspect=1,
                         extent=[mesh_q_2_min, mesh_q_2_max, mesh_q_1_min,
                             mesh_q_1_max], vmin=0, vmax=6,cmap='inferno')
    plt.subplots_adjust(bottom=0.15, left=0.2)
    ax.set_xlabel(r"$\rm Q_"+q_2+r"$ $(\rm\AA^{-1})$", fontsize=3)
    ax.set_ylabel(r"$\rm Q_"+q_1+r"$ $(\rm\AA^{-1})$", fontsize=3)
    ax.tick_params(axis='both', labelsize=3)
    ax.xaxis.set_tick_params(width=0.3)
    ax.yaxis.set_tick_params(width=0.3)
    output_file_name = file_name[:-4] + '_Q' + q_1 + '_v_Q' + q_2
    output_file_type = '.jpg'
    # plt.savefig(existential_check(output_file_name, output_file_type, output_folder))

    plt.show()


scan = [182]
file_index = [-1]
axis_1 = 'x'
axis_2 = 'z'
axis_3_limits= [0,-1]
axis_3_limits = [1,35]

f_1, f_name = file_checker(scan[0], file_index[0],'')
rsm_plot(f_1, f_name, axis_1, axis_2, axis_3_limits,'')


