import numpy as np
import os

def trident(all_images):
    delQ1 = pixel_estimator(all_images[0], all_images[1])
    #delQ2 = pixel_estimator(all_images[-2], all_images[-1])
    #delQ3 = pixel_estimator(all_images[int(len(all_images)/2)],
    #                        all_images[int(len(all_images)/2)+1])
    return [delQ1[0],#delQ2[0]),#, delQ3[0]),
            delQ1[1],#delQ2[1]),#, delQ3[1]),
            delQ1[2]]#delQ2[2])]#, delQ3[2])]
    

def pixel_estimator(im1, im2):
    del_Qx = union_jack(im1.Q_x, im2.Q_x)
    del_Qy = union_jack(im1.Q_y, im2.Q_y)
    del_Qz = union_jack(im1.Q_z, im2.Q_z)
    return [del_Qx, del_Qy, del_Qz]

def union_jack(Q_1, Q_2):
    delQ = []
    #delQ.append(max_spacing_middle(Q_1,Q_2))
    delQ.append(max_spacing_corner(Q_1,Q_2))
    delQ.append(max_spacing_edge(Q_1,Q_2))
    Q_1 = rotate(Q_1)
    Q_2 = rotate(Q_2)
    delQ.append(max_spacing_corner(Q_1,Q_2))
    delQ.append(max_spacing_edge(Q_1,Q_2))
    Q_1 = rotate(Q_1)
    Q_2 = rotate(Q_2)
    delQ.append(max_spacing_corner(Q_1,Q_2))
    delQ.append(max_spacing_edge(Q_1,Q_2))
    Q_1 = rotate(Q_1)
    Q_2 = rotate(Q_2)
    delQ.append(max_spacing_corner(Q_1,Q_2))
    delQ.append(max_spacing_edge(Q_1,Q_2))
    return max(delQ)

def rotate(array_2d):
    list_of_tuples = zip(*array_2d[::-1])
    return [list(elem) for elem in list_of_tuples]

def max_spacing_corner(Q_1, Q_2):
    return max([abs(Q_1[0][0]-Q_1[0][1])/len(Q_1[0]),
               abs(Q_1[0][0]-Q_1[1][0])/len(Q_1),
               abs(Q_2[0][0]-Q_1[0][0])/200])

def max_spacing_edge(Q_1, Q_2):
    mid_coord = int(len(Q_1[0])/2)
    return max([abs(Q_1[0][mid_coord]-Q_1[0][mid_coord+1])/len(Q_1[0]),
               abs(Q_1[0][mid_coord]-Q_1[1][mid_coord])/len(Q_1),
               abs(Q_2[0][mid_coord]-Q_1[0][mid_coord])/200])


def max_spacing_middle(Q_1, Q_2):
    mid_coord1 = int(len(Q_1[0])/2)
    mid_coord2 = int(len(Q_1)/2)
    return max([abs(Q_1[mid_coord2][mid_coord1]-Q_1[mid_coord2][mid_coord1+1]),
               abs(Q_1[mid_coord2][mid_coord1]-Q_1[mid_coord2+1][mid_coord1]),
               abs(Q_2[mid_coord2][mid_coord1]-Q_1[mid_coord2][mid_coord1])])


def Q_limit_dict_maker(directory, spot_dict, scan_num, Q_max, Q_min, NR_PTS):
    qlim_dict = {'qz_min' : Q_min[2],
                  'qz_max' : Q_max[2],
                  'qx_min' : Q_min[0],
                  'qx_max' : Q_max[0],
                  'qy_min' : Q_min[1],
                  'qy_max' : Q_max[1]}
#    q_x_lim_min = del_Q_all[0]
#    q_y_lim_min = del_Q_all[1]
#    q_z_lim_min = del_Q_all[2]
#    nr_pts = 1 
#    while abs(qlim_dict['qx_max']-qlim_dict['qx_min'])/nr_pts > q_x_lim_min \
#      and abs(qlim_dict['qy_max']-qlim_dict['qy_min'])/nr_pts > q_y_lim_min \
#      and abs(qlim_dict['qz_max']-qlim_dict['qz_min'])/nr_pts > q_z_lim_min:
#              nr_pts = nr_pts + 1                                              
#                                                                               
    qlim_dict['nr_pts'] = int(0.75*NR_PTS)  

    filename = directory\
               + 'user_defined_parameters/qlim/qlim_'\
               + spot_dict[str(scan_num[0])]\
               + '.txt'                                                        
                                                                               
    if os.path.exists(filename) == True \
        and input('Overwrite existing file, [y] or n?\n') != 'y':              
            pass                                                               
    else:                                                                      
       with open(filename,'w') as inf:                                         
          inf.write(str(qlim_dict))                                  
       print('Created file',filename)                                          

