import os
import fabio
import numpy as np

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

        self.two_theta_horizontal_range = []
        self.two_theta_vertical_range = []
        self.generate_two_thetas() # Calculates range of two_theta
        self.omega_offsetter() # Tempory fix- Didier will help
        self.slice_and_dice() # Chops out the dead rows and selects window
        self.file.close()
        self.calc_Q_coordinates()


    def calc_Q_coordinates(self):
        for coordinate_1 in range(np.shape(self.data)[1]):
            row = []
            row_x = []
            row_y = []
            row_z = []
            for coordinate_2 in range(np.shape(self.data)[0]):
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
                                            + self.param["Zero_Pixel_Ver"] 
                                            - float(pixel) 
                                            * self.param["det_ang"])
        for pixel in range(self.param["ShotHeight_y"]):
            self.two_theta_vertical_range.append((self.param['Zero_Pixel_Hor'] 
                                                  - self.param['chi_offset'] 
                                                  - float(pixel)) 
                                                 * self.param['det_ang'])


    def calc_qx(self, coordinate_1, coordinate_2): # c1 = i  c2 = j
        return ((2 * np.pi / self.WAVELENGTH) 
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
   

    def calc_qy(self, coordinate_1, coordinate_2):
        return ((2 * np.pi / self.WAVELENGTH) 
                * np.sin(self.two_theta_vertical_range[coordinate_1] 
                        / 180 
                        * np.pi)
                * np.cos(self.two_theta_horizontal_range[coordinate_2] 
                        / 180 
                        * np.pi 
                        - self.OMEGA 
                        / 180 
                        * np.pi))

    
    def calc_qz(self, coordinate_2):
        return ((2 * np.pi / self.WAVELENGTH) 
                * (np.sin(self.OMEGA 
                         / 180 
                         * np.pi)) 
                + np.sin(self.two_theta_horizontal_range[coordinate_2] 
                        / 180 
                        * np.pi 
                        - self.OMEGA 
                        / 180 
                        * np.pi))
    
    
