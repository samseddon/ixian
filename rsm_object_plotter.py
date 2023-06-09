import pickle   
import os
import glob
import numpy as np
from nexusformat import nexus
from matplotlib.cm import ScalarMappable
import matplotlib.pyplot as plt
from scipy import interpolate
from equations import existential_check
from colour_bar import code_TUD_cbar as cbar
plt.style.use("seddon_TUD") 

def file_checker(s_num, s_ind, input_path):
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

#def existential_check(o_f, f_s, folder): #file name, file suffix, destination folder
#    f_n = str(o_f + f_s)
#    while os.path.exists(folder + f_n):
#        if len(f_n) > len(o_f + f_s):
#            filenum = int(f_n[(len(o_f + f_s) - (len(f_s)-1)):-len(f_s)]) + 1
#            f_n = f_n[:(len(o_f + f_s) - (len(f_s)-1))] + str(filenum) + f_n[-len(f_s):]
#        else:
#            f_n = f_n[:-len(f_s)] + '_1' + f_n[-len(f_s):]
#    return folder + f_n

def input_check(q_1, q_2):
    axis_dict = {'qx': 0, 'qX': 0, 'qy': 1, 'qY': 1, 'qz': 2, 'qZ': 2}
    if q_1 not in axis_dict or q_2 not in axis_dict:
        print('Axis not recognised, please use qx, qy, qz, T, or M')
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

def slicer_and_dicer_3000(f_1, file_name, q_1, q_2, q_3, output_folder):
    volume = np.array(f_1.data)
    q_x = np.array(f_1.q_x)
    q_y = np.array(f_1.q_y)
    q_z = np.array(f_1.q_z)


    # Test array init
    axis0 = []

    for i in range(100):
        axis1 = []
        for j in range(100):
            axis2 = []
            for k in range(100):
                axis2.append(1)
            axis1.append(axis2)
        axis0.append(axis1)
    

    Q_vec = np.array([1,1,1])

    Q_vec_perp1 = np.random.randn(3)
    Q_vec_perp1 -= Q_vec_perp1.dot(Q_vec) * Q_vec / np.linalg.norm(Q_vec)**2
    Q_vec_perp1 /= np.linalg.norm(Q_vec_perp1)
    Q_vec_perp2 = np.cross(Q_vec, Q_vec_perp1)
    
    d_1 = 0
    d_2 = 0 

    Q_init_1 = [0,0,0]
    Q_init_2 = [1,1,1]
    pass_num = 0 
    c = 0
    test = []
    while c<100:
        pass_num = 0 
        for i in range(5):
            for j in range(5): 
                for k in range(5): 
                    #print(i, j, k, "Planes", plane(Q_vec, Q_init_1, i, j, k) ,  plane(Q_vec, Q_init_2, i, j, k))
                    low_plane = plane(Q_vec, Q_init_1, i, j, k)
                    hig_plane = plane(Q_vec, Q_init_2, i, j, k)
                    if low_plane >= 0 and low_plane <= (sum(Q_vec)) and hig_plane <= 0:
                        pass_num +=1
        
        test.append([c, pass_num])
        c += 1
        Q_init_1 = Q_init_1 + Q_vec
        print(Q_init_1)
        Q_init_2 = Q_init_2 + Q_vec
        print(Q_init_2)
    print(test)

def plane(Q_vec, Q_init, x, y, z): 
    return Q_vec[0] * x + Q_vec[1]* y + Q_vec[2] * z - (Q_vec[0] * Q_init[0] + Q_vec[1] * Q_init[1]+ Q_vec[2] * Q_init[2]) 

    



def rsm_plot(f_1, file_name, q_1, q_2, q_3_lim, output_folder):
    q_1_assigned, q_2_assigned, q_3_assigned = input_check(q_1, q_2)
    if q_1 == 'qx':
        values_1 = f_1.q_x
    if q_1 == 'qy':
        values_1 = f_1.q_y
    if q_1 == 'qz':
        values_1 = f_1.q_z
    if q_2 == 'qx':
        values_2 = f_1.q_x
    if q_2 == 'qy':
        values_2 = f_1.q_y
    if q_2 == 'qz':
        values_2 = f_1.q_z

    
    volume = np.array(f_1.data)
#    for i in range(volume.shape[0]):
#        for j in range(volume.shape[1]):
#            for k in range(volume.shape[2]):
#                if volume[i,j,k] >0:
#                    print('good')
    q1_q2_matrix = np.zeros((volume.shape[0], volume.shape[2])) # THIS IS ALWAYS 0 and 2
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
    fig, ax = plt.subplots(figsize = (5,4), dpi = 128)
    #cb = img = ax.imshow(np.log10(grid_q), origin='lower', aspect=1,
    #                     extent=[mesh_q_2_min, mesh_q_2_max, mesh_q_1_min,
    #                         mesh_q_1_max], vmin=0, vmax=6,cmap='inferno')
    grid_q = np.log10(grid_q)
    vmin = 1
    vmax = 6
    img = ax.contourf(grid_q, origin='lower',           
                        extent=[mesh_q_2_min, mesh_q_2_max, mesh_q_1_min,mesh_q_1_max],
            levels = 256, vmin=vmin, vmax=vmax,
                         cmap=cbar) 
    ax.contour(np.log10(grid_q), origin='lower',         
                    extent=[mesh_q_2_min, mesh_q_2_max, mesh_q_1_min,mesh_q_1_max],
            levels = 64, colors = "#00305d", linewidths = 0.1, vmin = vmin)
    cobar = fig.colorbar(
                ScalarMappable(norm=img.norm, cmap=img.cmap),
                ticks=range(vmin, vmax+1, 1))
    cobar.ax.tick_params(size=0)
    cobar.set_label('Log(Intensity) (A. U.)', rotation=270, labelpad = 30)
    #plt.axis('square')
    #plt.subplots_adjust(bottom=0.15, left=0.2)
    ax.set_xlabel(r"\rm Q$_"+q_2[-1]+r"$ $(\rm\AA^{-1})$")
    ax.set_ylabel(r"\rm Q$_"+q_1[-1]+r"$ $(\rm\AA^{-1})$")
    output_file_name = file_name[:-4] + '_Q' + q_1 + '_v_Q' + q_2
    output_file_type = '.png'
    plt.savefig(existential_check(output_file_name, output_file_type, output_folder), bbox_inches='tight')

    plt.show()




#scan = [180]
#file_index = [-1] #for most recent processed data just put -1, otherwise 0 will do the first or its index
#axis_1 = 'qx'
#axis_2 = 'qy'
#axis_3_limits= [0,-1]
#axis_3_limits = [1,35]
#data_directory = '/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/processed_files/'
#output_directory = '/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/processed_files/images'
#
#
#f_1, f_name = file_checker(scan[0], file_index[0],data_directory)
#rsm_plot(f_1, f_name, axis_1, axis_2, axis_3_limits,output_directory)


