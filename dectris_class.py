import os
import fabio
import numpy as np
"""
Created on Wed Jan 11 14:51:58 2023

@author: samseddon
"""

#class Q_Space():
#    def __init__(self):
#        self.data = []
#        self.Q_x = []
#        self.Q_y = []
#        self.Q_z = []
#
#
#    def add_image(self, image):
#
#        self.data.append(image.data)
#        self.Q_x.append(image.Q_x)
#        self.Q_y.append(image.Q_y)
#        self.Q_z.append(image.Q_z)




class Q_Space():
    def __init__(self, scan_num, spot_dict, directory):
        self.scan_num = scan_num
        self.spot_dict = spot_dict
        self.directory = directory
        
        self.qlim_dict = self.qlim()
        
        self.QX_MIN = self.qlim_dict["qx_min"]
        self.QX_MAX = self.qlim_dict["qx_max"]
        self.QY_MIN = self.qlim_dict["qy_min"]
        self.QY_MAX = self.qlim_dict["qy_max"]
        self.QZ_MIN = self.qlim_dict["qz_min"]
        self.QZ_MAX = self.qlim_dict["qz_max"]
        self.NR_PTS = self.qlim_dict["nr_pts"]

#        self.NR_PTS = 84
        self.data = np.zeros((self.NR_PTS,                           
                              self.NR_PTS,                           
                              self.NR_PTS))                          
        
        self.q_idx = np.zeros((self.NR_PTS,                           
                               self.NR_PTS,                           
                               self.NR_PTS))                          


        self.q_x = np.linspace(self.QX_MIN,                       
                               self.QX_MAX,                       
                               self.NR_PTS)                       
                                                                               
        self.q_y = np.linspace(self.QY_MIN,                       
                               self.QY_MAX,                       
                               self.NR_PTS)                       
                                                                               
        self.q_z = np.linspace(self.QZ_MIN,                       
                               self.QZ_MAX,                       
                               self.NR_PTS)


    def normalise_3D(self):
        for s in range(self.data.shape[0]):                                 
            for t in range(self.data.shape[1]):                             
                for u in range(self.data.shape[2]):                         

                    if self.q_idx[s, t, u] == 0:
                        pass
                    else:
                        self.data[s, t, u] = self.data[s, t, u] / self.q_idx[s, t, u]
# X_Z WORKS

    def fix_qxz(self):
        fix_idx = 1
        fix_val = 0
        for s in range(self.data.shape[0]):                                 
            for t in range(self.data.shape[1]):                             
                for u in range(self.data.shape[2]):                         
                    if self.q_idx[s, t, u] != 0:
                        fix_idx = self.q_idx[s, t, u]
                        fix_val = self.data[s, t, u]
                    if self.q_idx[s, t, u] == 0:                             
                        self.q_idx[s, t, u] = 0
                        self.data[s, t, u] =fix_val

    def fix_qxy(self): 
        fix_idx = 1
        fix_val = 0
        for u in range(self.data.shape[2]):                         
            for s in range(self.data.shape[0]):                                 
                for t in range(self.data.shape[1]):                             
                    if self.q_idx[s, t, u] != 0:
                        fix_idx = self.q_idx[s, t, u]
                        fix_val = self.data[s, t, u]
                    if self.q_idx[s, t, u] == 0:                             
                        self.q_idx[s, t, u] = 0
                        self.data[s, t, u] =fix_val



# YZ doesn't work :( 

#        fix_idx = 1
#        fix_val = 0
#        for u in range(self.data.shape[1]):                             
#            for t in range(self.data.shape[2]):                         
#                for s in range(self.data.shape[0]):                                 
#                    if            self.q_idx[s, t, u] != 0:
#                        fix_idx = self.q_idx[s, t, u]
#                        fix_val = self.data[s, t, u]
#                    if  self.q_idx[s, t, u] == 0:                             
#                        self.q_idx[s, t, u] = 0
#                        self.data[s, t, u] =fix_val
                                                                              
                                                                               


    def qlim(self):
        import_var = self.directory \
                + "user_defined_parameters/qlim/qlim_" \
                + self.spot_dict[str(self.scan_num[0])] \
                + ".txt"

        if os.path.exists(import_var) == False:
            raise(ValueError("qlim file does not exist, consider running"\
                         "again with create_files = True"))
        with open(import_var,"r") as inf:
            dict1 = eval(inf.read())
        return dict1
    
    def populate_3D(self, dec_image):
        for coordinate_1 in range(np.shape(dec_image.data)[1]):
            for coordinate_2 in range(np.shape(dec_image.data)[0]):
                Q_index = self.find_q_index(dec_image.Q_coords[coordinate_1][coordinate_2])
                if Q_index[0] == False:
                    pass
                else:
                    self.add_point_3D(Q_index[1], dec_image.data[coordinate_2][coordinate_1])

    
    def add_point_3D(self, indexes, data_point):
        self.data[indexes[0], indexes[1], indexes[2]] = \
                self.data[indexes[0], indexes[1], indexes[2]] + data_point
        self.q_idx[indexes[0], indexes[1], indexes[2]] += 1

    
    def find_q_index(self, coord):
        qx = coord[0]
        qy = coord[1]
        qz = coord[2]
        if qx < self.QX_MIN \
                or qx > self.QX_MAX \
                or qy < self.QY_MIN \
                or qy > self.QY_MAX \
                or qz < self.QZ_MIN \
                or qz > self.QZ_MAX:
            return (False,[0,0,0])
        else:
            for i_x in range(len(self.q_x)):
                if(qx > self.q_x[i_x]):
                     pass
                else:
                    break

            for i_y in range(len(self.q_y)):
                if(qy > self.q_y[i_y]):
                    pass
                else:
                    break

            for i_z in range(len(self.q_z)):
                if(qz > self.q_z[i_z]):
                    pass
                else:
                    break
            return (True, [i_x, i_y, i_z])



class Dectris_Image():
    def __init__(self, file_reference, image_number, directory, master_files, param):
        self.file_reference = file_reference
        self.image_number = image_number
        self.file = fabio.open(os.path.join(directory + "data/", 
                                            master_files[image_number]))
        self.param = param
        self.DEADROW1 = param['DeadRow1']                                               
        self.DEADROW2 = param['DeadRow2']                                               
        self.DEADROW3 = param['DeadRow3']                                               
        self.DEADROW4 = param['DeadRow4']   
        self.REALX_LIM_LOW = param['realx_lim_low']                                     
        self.REALX_LIM_HIG = param['realx_lim_hig']                                     
        self.REALY_LIM_LOW = param['realy_lim_low']                                     
        self.REALY_LIM_HIG = param['realy_lim_hig']  

        ub = np.array(self.file.header.get('UB_pos').split(' '))
        self.ub = ub.astype(float)
        self.count_mne   = self.file.header.get('counter_mne').split(' ')
        self.count_pos   = self.file.header.get('counter_pos').split(' ')
        self.motor_mne   = self.file.header.get('motor_mne').split(' ')
        self.motor_pos   = self.file.header.get('motor_pos').split(' ')
        self.WAVELENGTH  = float(self.file.header.get('source_wavelength').split(' ')[0])
        
        self.dict_count = {}
        self.dict_motor = {}
        for key in self.count_mne:
            for value in self.count_pos:
                self.dict_count[key] = value
                self.count_pos.remove(value)
                break
    
        for key in self.motor_mne:
            for value in self.motor_pos:
                self.dict_motor[key] = value
                self.motor_pos.remove(value)
                break

        self.ATTN        =   float(self.dict_count['attn']) 
        self.TRANS       =   float(self.dict_count['trans'])
        self.IO          =   float(self.dict_count['Io'])   
        self.TEMP        =   float(self.dict_count['temp_B'])
        
        self.TWO_THETA   =   float(self.dict_motor['del'])
        self.OMEGA       =   float(self.dict_motor['eta'])
        self.CHI         =   float(self.dict_motor['chi'])
        self.PHI         =   float(self.dict_motor['phi'])
        if "ami_mag" in self.dict_motor:
            self.MAG_FIELD = float(self.dict_motor['ami_mag'])
        else:
            self.MAG_FIELD = 0
        
        self.data = []
        self.Q_x = []
        self.Q_y = []
        self.Q_z = []
        self.Q_coords = []
        self.pixel_list = []
        self.omega_offsetter() # Tempory fix- Didier will help

        self.two_theta_horizontal_range = []
        self.two_theta_vertical_range = []
        self.generate_two_thetas() # Calculates range of two_theta
        self.slice_and_dice() # Chops out the dead rows and selects window
        self.file.close()
        self.calc_Q_coordinates()


    def q_lim(self):
        return (np.amin(self.Q_x), np.amax(self.Q_x)),\
               (np.amin(self.Q_y), np.amax(self.Q_y)),\
               (np.amin(self.Q_z), np.amax(self.Q_z))


    def calc_Q_coordinates(self):
        for coordinate_1 in range(np.shape(self.data)[1]):
            row = []
            row_x = []
            row_y = []
            row_z = []
            for coordinate_2 in range(np.shape(self.data)[0]):
                self.pixel_list.append([self.calc_qx(coordinate_1, coordinate_2),
                                 self.calc_qy(coordinate_1, coordinate_2),
                                 self.calc_qz(coordinate_2),
                                 self.data[coordinate_2, coordinate_1]])
#                self.data = [self.calc_qx(coordinate_1, coordinate_2),
#                             self.calc_qy(coordinate_1, coordinate_2),
#                             self.calc_qz(coordinate_2),
#                             self.data[coordinate_2][coordinate_1]]
                row_x.append(self.calc_qx(coordinate_1, coordinate_2))
                row_y.append(self.calc_qy(coordinate_1, coordinate_2))
                row_z.append(self.calc_qz(coordinate_2))
                row.append([row_x[-1], row_y[-1], row_z[-1]])
            self.Q_x.append(row_x)
            self.Q_y.append(row_y)
            self.Q_z.append(row_z)
            self.Q_coords.append(row)

    def omega_offsetter(self):
        if self.file_reference == "MAG001":
            if self.PHI > -50:
                self.OMEGA = self.OMEGA + 1.6923
            else:
                self.OMEGA = self.OMEGA + 3.8765


    def slice_and_dice(self):
        data = np.concatenate((self.file.data[:self.DEADROW1 - 1, :],
                    self.file.data[self.DEADROW2 + 1:self.DEADROW3 - 1, :],
                    self.file.data[self.DEADROW4 + 1:, :]), axis=0)
        two_theta_range_new = np.concatenate((
                      self.two_theta_horizontal_range[:self.DEADROW1 - 1],
                      self.two_theta_horizontal_range[self.DEADROW2 
                                                     + 1:self.DEADROW3 - 1],
                      self.two_theta_horizontal_range[self.DEADROW4 + 1:]))
        data = np.array(data)
        self.data = data[self.REALY_LIM_LOW:self.REALY_LIM_HIG
                        ,self.REALX_LIM_LOW:self.REALX_LIM_HIG]
        self.two_theta_horizontal_range =\
                two_theta_range_new[self.REALY_LIM_LOW:self.REALY_LIM_HIG]
        self.two_theta_vertical_range =\
          self.two_theta_vertical_range[self.REALX_LIM_LOW:self.REALX_LIM_HIG]


    def generate_two_thetas(self):
        for pixel in range(self.param["ShotHeight"]):
            self.two_theta_horizontal_range.append(self.TWO_THETA 
                                            + (self.param["Zero_Pixel_Ver"] 
                                            - float(pixel)) 
                                            * self.param["det_ang"])
        for pixel in range(self.param["ShotHeight_y"]):
            self.two_theta_vertical_range.append((self.param['Zero_Pixel_Hor'] 
                                                  - self.param['chi_offset'] 
                                                  - float(pixel)) 
                                                 * self.param['det_ang'])


    def calc_qx(self, coordinate_1, coordinate_2): # c1 = i  c2 = j
        qx = ((2 * np.pi / self.WAVELENGTH) 
                * (np.cos(self.OMEGA / 180 * np.pi) 
                - np.cos(self.two_theta_horizontal_range[coordinate_2] 
                        / 180 
                        * np.pi 
                        - self.OMEGA 
                        / 180 
                        * np.pi)
                 * np.cos(self.two_theta_vertical_range[coordinate_1] 
                         / 180 
                         * np.pi)))
        return(qx) #0.02748846124221295


    def calc_qy(self, coordinate_1, coordinate_2):
        qy =  ((2 * np.pi / self.WAVELENGTH) 
                * np.sin(self.two_theta_vertical_range[coordinate_1] 
                        / 180 
                        * np.pi)
                * np.cos(self.two_theta_horizontal_range[coordinate_2] 
                        / 180 
                        * np.pi 
                        - self.OMEGA 
                        / 180 
                        * np.pi))
        return qy 
    
    def calc_qz(self, coordinate_2):
        qz =    ((2 * np.pi / self.WAVELENGTH) 
                * (np.sin(self.OMEGA 
                         / 180 
                         * np.pi) 
                + np.sin(self.two_theta_horizontal_range[coordinate_2] 
                        / 180 
                        * np.pi 
                        - self.OMEGA 
                        / 180 
                        * np.pi)))
        return qz 
    
    
