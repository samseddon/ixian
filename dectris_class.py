import os
import matplotlib.pyplot as plt
import inspect
import fabio
import math
import numpy as np
from numba import jit
from equations import calc_qx, calc_qy, calc_qz, find_q_linear_fast
from colour_bar import code_TUD_cbar as cbar
"""
Created on Wed Jan 11 14:51:58 2023

@author: samseddon
"""
class Dectris_Image():
    def __init__(self, major_list, whole_image):
        self.file_reference = major_list[0]
        self.image_number = major_list[1]
        self.directory = major_list[2]
        self.master_files = major_list[3]
        file = fabio.open(os.path.join(self.directory, 
                                self.master_files[self.image_number]))
        self.param = major_list[-1]
#        self.DEADROW1 =      self.param['DeadRow1']                                               
#        self.DEADROW2 =      self.param['DeadRow2']                                               
#        self.DEADROW3 =      self.param['DeadRow3']                                               
#        self.DEADROW4 =      self.param['DeadRow4']   
        if self.param["WHOLE_IMAGE"] == True or whole_image  == True:
            self.REAL_HOR_LIM_LOW = 0
            self.REAL_HOR_LIM_HIG = self.param["TOTAL_PIXELS_HOR"] 
            self.REAL_VER_LIM_LOW = 0
            self.REAL_VER_LIM_HIG = self.param["TOTAL_PIXELS_VER"]
        else:
            self.REAL_HOR_LIM_LOW = self.param['REAL_HOR_LIM_LOW']
            self.REAL_HOR_LIM_HIG = self.param['REAL_HOR_LIM_HIG']                                     
            self.REAL_VER_LIM_LOW = self.param['REAL_VER_LIM_LOW']                                     
            self.REAL_VER_LIM_HIG = self.param['REAL_VER_LIM_HIG']  

        ub = np.array(file.header.get('UB_pos').split(' '))
        self.ub = ub.astype(float)
        self.count_mne   = file.header.get('counter_mne').split(' ')
        self.count_pos   = file.header.get('counter_pos').split(' ')
        self.motor_mne   = file.header.get('motor_mne').split(' ')
        self.motor_pos   = file.header.get('motor_pos').split(' ')
        self.WAVELENGTH  = float(file.header.get('source_wavelength').split(' ')[0])
        self.TIME = float(file.header.get('time').split(' ')[0])
        
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
#        for key in self.dict_count:
#            print(key)
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
        self.slice_and_dice(file) # Chops out the dead rows and selects window
        file.close()


    def q_lim(self):
        return (np.amin(self.Q_x), np.amax(self.Q_x)),\
               (np.amin(self.Q_y), np.amax(self.Q_y)),\
               (np.amin(self.Q_z), np.amax(self.Q_z))



    def omega_offsetter(self):
        if self.file_reference == "MAG001":
            if self.PHI > -50:
                self.OMEGA = self.OMEGA + 1.6923
            else:
                self.OMEGA = self.OMEGA + 3.8765


    def slice_and_dice(self, file):
        data = file.data
        two_theta_range_new = np.array(self.two_theta_horizontal_range)
        data = np.array(data)
        self.data = data[self.REAL_VER_LIM_LOW:self.REAL_VER_LIM_HIG
                        ,self.REAL_HOR_LIM_LOW:self.REAL_HOR_LIM_HIG]
        self.two_theta_horizontal_range =\
                two_theta_range_new[self.REAL_VER_LIM_LOW:self.REAL_VER_LIM_HIG]
        self.two_theta_vertical_range =\
          self.two_theta_vertical_range[self.REAL_HOR_LIM_LOW:self.REAL_HOR_LIM_HIG]


    def generate_two_thetas(self):
        for pixel in range(self.param["TOTAL_PIXELS_VER"]):
            self.two_theta_horizontal_range.append(self.TWO_THETA 
                                            + (self.param["ZERO_PIXEL_VER"] 
                                            - float(pixel)) 
                                            * self.param["DET_ANG"])
        for pixel in range(self.param["TOTAL_PIXELS_HOR"]):
            self.two_theta_vertical_range.append((self.param['ZERO_PIXEL_HOR'] 
                                                  - self.param['CHI_OFFSET'] 
                                                  - float(pixel)) 
                                                 * self.param['DET_ANG'])

    
    
