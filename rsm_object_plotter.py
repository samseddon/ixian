import pickle   
import inspect
import pprint as pp
import os
import glob
import time
import numpy as np
from nexusformat import nexus
from matplotlib.cm import ScalarMappable
import matplotlib.pyplot as plt
from scipy import interpolate
from equations import existential_check
from colour_bar import code_TUD_cbar as cbar
from mpl_toolkits import mplot3d
plt.style.use("seddon_TUD") 

def file_checker(s_num, s_ind, input_path):
    """
    """
    print("Function" \
            + str(inspect.currentframe()).split(",")[-1][5:-1] \
            + " called from"\
            + str(inspect.currentframe()).split(",")[1])
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

def input_check(q_1, q_2):
    """
    """
    print("Function" \
            + str(inspect.currentframe()).split(",")[-1][5:-1] \
            + " called from"\
            + str(inspect.currentframe()).split(",")[1])
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



def slicer_and_dicer_3000(f_1, file_name, Q_vec, output_folder, scan_num):
    """
    """
    print("Function" \
            + str(inspect.currentframe()).split(",")[-1][5:-1] \
            + " called from"\
            + str(inspect.currentframe()).split(",")[1])
    
    volume = np.array(f_1.data)
    q_x = np.array(f_1.q_x)
    q_y = np.array(f_1.q_y)
    q_z = np.array(f_1.q_z)
    limits = np.shape(volume)
    idx_test = []
    for i in range(limits[0]):
        idx1 = []
        for j in range(limits[1]):
            idx2 = []
            for k in range(limits[2]):
                idx2.append(0)
            idx1.append(idx2)
        idx_test.append(idx1)

    xrange = (max(q_x) - min(q_x)) 
    yrange = (max(q_y) - min(q_y)) 
    zrange = (max(q_z) - min(q_z)) 
    Q_norm = [0,0,0]                                                                           
                                                                                                                     
                                              
    Q_norm[0] = float(Q_vec[0])/xrange
    Q_norm[1] = float(Q_vec[1])/yrange
    Q_norm[2] = float(Q_vec[2])/zrange
    max_val = max(Q_norm)
    Q_real = []
    for _ in range(len(Q_norm)):
        Q_norm[_] = Q_norm[_]/max_val
        Q_real.append(Q_norm[_])
    
    alpha = np.dot([min(q_x),min(q_y),min(q_z)], Q_vec) / (np.sqrt((Q_vec[0]**2 + Q_vec[1] **2 + Q_vec[2]**2)) **2)
    print(Q_norm)
    Q_init_1 = [0,0,0]
    Q_init_2 = [0,0,0]
    print(Q_init_1)
    for _ in range(len(Q_init_1)):
        Q_init_2[_] = Q_init_1[_] + Q_norm[_]
    print(Q_init_2)
    c = 0
    test = []
    integral = []
    size = max(np.shape(volume))
    start_time = time.time()
    while c<size:
#    while c<size:
        if np.average(idx_test) == 1:
            break
        else:
            Q_real[0] =  alpha * Q_vec[0] + c * (Q_norm[0] * xrange/size)
            Q_real[1] =  alpha * Q_vec[1] + c * (Q_norm[1] * yrange/size)
            Q_real[2] =  alpha * Q_vec[2] + c * (Q_norm[2] * zrange/size)

            #print((time.time() - start_time)/60)
            checker = []
            pass_num = 0 
            plane_count = 0
            for i in range(size):
           #     print("i = ", i)
                for j in range(size): 
                    for k in range(size): 
                        if idx_test[i][j][k] == 1:
                            checker.append(1)
                            pass
                        else:
                           # print(i, j, k, "Planes", plane_check(Q_norm, Q_init_1, [i, j, k]) ,  plane_check(Q_norm, Q_init_2, [i, j, k]))
                            low_plane = plane_check(Q_norm, Q_init_1, ([i-Q_norm[0],j-Q_norm[1],k-Q_norm[2]]))
                            hig_plane = plane_check(Q_norm, Q_init_2, ([i+Q_norm[0],j+Q_norm[1],k+Q_norm[2]]))
                  #          print(i,j,k,[i,j,k]-Q_vec,low_plane, hig_plane)
                            if low_plane < 0 and hig_plane >= 0:
    #                            print('here')
                                pass_num +=1
                                idx_test[i][j][k] = 1
                                plane_count += volume[i][j][k]
            #                pp.pprint(idx_test)
            # This is extra
            #pass_num = 1
            #plane_count = 1
            # This is normal again
            print("pass_num = ", pass_num)
            print("plane_count =", plane_count)
            integral.append([c, plane_count])
            test.append([np.sqrt(Q_real[0]**2 + Q_real[1]**2 + Q_real[2]**2), pass_num])
            c += 1
            print(np.average(idx_test))
            #print(c)
            #print(len(checker))
            for _ in range(len(Q_init_1)):
                Q_init_1[_] = Q_init_1[_] + Q_norm[_]
                Q_init_2[_] = Q_init_2[_] + Q_norm[_]
            
            print(c, Q_init_1)
            print(c, Q_init_2)
    x = []
    y = []
    print(integral, test)
    for i in range(len(test)):
        x.append(test[i][0])
        if test[i][1] != 0:
            y.append(integral[i][1]/test[i][1])
        else:
            y.append(0)
    print(x, y)
    np.savetxt("krys_scans/" + str(scan_num) + 'x.txt', x)
    np.savetxt("krys_scans/" + str(scan_num) + 'y.txt', y)
    #plt.plot(x,y)
    #plt.close()
    

def test_slicer(Q_vec):
    """
    """
    print("Function" \
            + str(inspect.currentframe()).split(",")[-1][5:-1] \
            + " called from"\
            + str(inspect.currentframe()).split(",")[1])
    idx_test = []
    size = 40

    for i in range(size):
        axis1 = []
        idx1 = []
        for j in range(size):
            idx2 = []
            axis2 = []
            for k in range(size):
                idx2.append(0)
                axis2.append(10)
            axis1.append(axis2)
            idx1.append(idx2)
        axis0.append(axis1)
        idx_test.append(idx1)
     
    #Q_vec = np.array([1.4,1.2,0.3])

    Q_vec_perp1 = np.random.randn(3)
    Q_vec_perp1 -= Q_vec_perp1.dot(Q_vec) * Q_vec / np.linalg.norm(Q_vec)**2
    Q_vec_perp1 /= np.linalg.norm(Q_vec_perp1)
    Q_vec_perp2 = np.cross(Q_vec, Q_vec_perp1)
    
    d_1 = 0
    d_2 = 0 

    Q_init_1 = [0,0,0]
    Q_init_2 = Q_init_1 + Q_vec
    pass_num = 0 
    c = 0
    test = []
    while c<size:
        checker = []
        pass_num = 0 
        for i in range(size):
            for j in range(size): 
                for k in range(size): 
                    if idx_test[i][j][k] == 1:
                        checker.append(1)
                        pass
                    else:
                    #print(i, j, k, "Planes", plane(Q_vec, Q_init_1, i, j, k) ,  plane(Q_vec, Q_init_2, i, j, k))
                        low_plane = plane_check(Q_vec, Q_init_1, ([i,j,k]-Q_vec))
                        hig_plane = plane_check(Q_vec, Q_init_2, ([i,j,k]+Q_vec))
              #          print(i,j,k,[i,j,k]-Q_vec,low_plane, hig_plane)
                        if low_plane < 0 and hig_plane >= 0:
                            pass_num +=1
                            idx_test[i][j][k] = 1
            #                pp.pprint(idx_test)
        
        test.append([c, pass_num])
        c += 1
        print(np.average(idx_test))
        #print(c)
        #print(len(checker))
        Q_init_1 = Q_init_1 + Q_vec
        #print(c, Q_init_1)
        Q_init_2 = Q_init_2 + Q_vec
        #print(c, Q_init_2)
    x = []
    y = []
    for point in test:
        x.append(point[0])
        y.append(point[1])
    plt.plot(x,y)
    plt.show()

def plane_check(Q_vec, Q_init, Q_test): 
    """
    """
    print("Function" \
            + str(inspect.currentframe()).split(",")[-1][5:-1] \
            + " called from"\
            + str(inspect.currentframe()).split(",")[1])

    return (Q_vec[0] * Q_test[0] + Q_vec[1]* Q_test[1] + Q_vec[2] * Q_test[2] - (Q_vec[0] * Q_init[0] + Q_vec[1] * Q_init[1]+ Q_vec[2] * Q_init[2]))/(Q_vec[0] + Q_vec[1] + Q_vec[2])


def threeD_rsm_plot(f_1):
    volume = f_1.data
    print(volume)
    q_x = f_1.q_x
    q_y = f_1.q_y
    q_z = f_1.q_z
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter(q_x, q_y, q_z, volume)#, cmap=cbar)
    plt.show()







def rsm_plot(f_1, 
             file_name, 
             q_1, 
             q_2, 
             q_3_lim, 
             output_folder, 
             directory, 
             scan_num, 
             fig,
             ax,
             show_plot = True):
             
    
    """
    """
    print("Function" \
            + str(inspect.currentframe()).split(",")[-1][5:-1] \
            + " called from"\
            + str(inspect.currentframe()).split(",")[1])

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
    with open(directory+'user_defined_parameters/spot_dict.txt','r') as inf:
        spot_dict = eval(inf.read())

    spot = spot_dict[str(scan_num[0])]
    volume = np.array(f_1.data)
#    for i in range(volume.shape[0]):
#        for j in range(volume.shape[1]):
#            for k in range(volume.shape[2]):
#                if volume[i,j,k] >0:
#                    print('good')
    if q_3_assigned == 0: # int x
        q1_q2_matrix = np.zeros((volume.shape[1], volume.shape[2])) # THIS IS ALWAYS 0 and 2
        for i in range(volume.shape[2]):  # this axis 0 for x
            q1_q2_matrix[:, i] = np.copy(volume[q_3_lim[0]:q_3_lim[1],:, i].mean(axis=0) + 10 ** (-8))


    
    elif q_3_assigned == 1: # int y
        q1_q2_matrix = np.zeros((volume.shape[0], volume.shape[2])) # THIS IS ALWAYS 0 and 2
        for i in range(volume.shape[2]):  # this axis 0 for y
            q1_q2_matrix[:, i] = np.copy(volume[:, q_3_lim[0]:q_3_lim[1], i].mean(axis=1) + 10 ** (-8))
            #print(q1_q2_matrix)
    
    elif q_3_assigned == 2: #int z
        q1_q2_matrix = np.zeros((volume.shape[0], volume.shape[1])) # THIS IS ALWAYS 0 and 2
        for i in range(volume.shape[1]):  # this axis 0 for z
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
    if q_3_assigned == 0:
        q1_q2_matrix = np.reshape(fin, (volume.shape[1], volume.shape[2]))
    elif q_3_assigned == 1:
        q1_q2_matrix = np.reshape(fin, (volume.shape[0], volume.shape[2]))
    elif q_3_assigned == 1:
        q1_q2_matrix = np.reshape(fin, (volume.shape[0], volume.shape[1]))

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
    #fig, ax = plt.figure(figsize = (5,4), dpi = 128)
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
    #cobar = fig.colorbar(
    #            ScalarMappable(norm=img.norm, cmap=img.cmap),
    #            ticks=range(vmin, vmax+1, 1))
    #cobar.ax.tick_params(size=0)
    #cobar.set_label('Log(Intensity) (A. U.)', rotation=270, labelpad = 30)
    #plt.axis('square')
    #plt.subplots_adjust(bottom=0.15, left=0.2)
    ax.set_xlabel(r"\rm Q$_"+q_2[-1]+r"$ $(\rm\AA^{-1})$")
    ax.set_ylabel(r"\rm Q$_"+q_1[-1]+r"$ $(\rm\AA^{-1})$")
    #plt.axis("square")
    #ax.set_title(spot)
    output_file_name = file_name[:-7] + q_1 + '_v_' + q_2
    output_file_type = '.png'
    return fig, ax

#plt.savefig(existential_check(output_file_name, output_file_type, output_folder), bbox_inches='tight')
    show_plot = False
    if show_plot == True:
        plt.show()
    else:
        plt.close()


def pickle_jar(filename, pickle_object):
    with open(filename,'wb') as handle:
        pickle.dump(pickle_object, handle, protocol=pickle.HIGHEST_PROTOCOL)

def pickle_unjar(filename):
    with (filename, 'rb') as handle:
        return pickle.load(handle)
