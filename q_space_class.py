import os
import matplotlib.pyplot as plt
import math
import inspect
import time
import numpy as np
from equations import find_q_linear_fast, existential_check, pickle_jar,\
        pickle_unjar, progress_bar
from colour_bar import code_TUD_cbar as cbar

class Q_Space():                                                               
    def __init__(self, scan_num, symmetric, SPACE_3D):              
        self.scan_num = scan_num                                             
                                                                               
        self.qlim_dict = self.qlim()                                           
        self.QX_MIN = self.qlim_dict["qx_min"]                                 
        self.QX_MAX = self.qlim_dict["qx_max"]                                 
        self.QY_MIN = self.qlim_dict["qy_min"]                                 
        self.QY_MAX = self.qlim_dict["qy_max"]                                 
        self.QZ_MIN = self.qlim_dict["qz_min"]                                 
        self.QZ_MAX = self.qlim_dict["qz_max"]                                 
        self.NX_PTS = self.qlim_dict["nx_pts"]                                 
        self.NY_PTS = self.qlim_dict["ny_pts"]                                 
        self.NZ_PTS = self.qlim_dict["nz_pts"]                                 
        self.NR_PTS = self.qlim_dict["nr_pts"]                                 
        self.time = 0

        self.SPACE_3D = SPACE_3D                                               
                                                                               
                                                                               
        if symmetric == True:                                                  
            self.QX_GRAD = (self.QX_MAX-self.QX_MIN) / (self.NR_PTS - 1)       
            self.QY_GRAD = (self.QY_MAX-self.QY_MIN) / (self.NR_PTS - 1)       
            self.QZ_GRAD = (self.QZ_MAX-self.QZ_MIN) / (self.NR_PTS - 1)       
                                                                               
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

        else:                                                                  
                                                                               
            self.QX_GRAD = (self.QX_MAX-self.QX_MIN) / (self.NX_PTS - 1)       
            self.QY_GRAD = (self.QY_MAX-self.QY_MIN) / (self.NY_PTS - 1)       
            self.QZ_GRAD = (self.QZ_MAX-self.QZ_MIN) / (self.NZ_PTS - 1)       
                                                                               
    #        self.NR_PTS = 84                                                  
            self.data = np.zeros((self.NX_PTS,                                 
                                  self.NY_PTS,                                 
                                  self.NZ_PTS))                                
                                                                               
            self.q_idx = np.zeros((self.NX_PTS,                                
                                   self.NY_PTS,                                
                                   self.NZ_PTS))                               
                                                                               
                                                                               
            self.q_x = np.linspace(self.QX_MIN,                                
                                   self.QX_MAX,                                
                                   self.NX_PTS)                                
                                                                               
            self.q_y = np.linspace(self.QY_MIN,                                
                                   self.QY_MAX,                                
                                   self.NY_PTS)                                
                                                                               
            self.q_z = np.linspace(self.QZ_MIN,                                
                                   self.QZ_MAX,                                
                                   self.NZ_PTS)                                
        self.find_q_lin = [self.QX_MAX,                                        
                           self.QY_MAX,                                        
                           self.QZ_MAX,                                        
                           self.QX_MIN,                                        
                           self.QY_MIN,                                        
                           self.QZ_MIN,                                        
                           self.QX_GRAD,                                       
                           self.QY_GRAD,                                       
                           self.QZ_GRAD]                                       


    def normalise_3D(self):                                                    
        for s in range(self.data.shape[0]):                                    
            for t in range(self.data.shape[1]):                                
                for u in range(self.data.shape[2]):                            
                                                                               
                    if self.q_idx[s, t, u] == 0:                               
                        pass                                                   
                    else:                                                      
                        self.data[s, t, u] = self.data[s, t, u] / self.q_idx[s, t, u]
    def qlim(self):                                                            
        import_var = "local/"\
                + "qlim/qlim_" \
                + str(self.scan_num[0]) \
                + ".txt"                                                       
                                                                               
        if os.path.exists(import_var) == False:                                
            raise(ValueError("qlim file does not exist, consider running"\
                         "again with create_files = True"))                    
        with open(import_var,"r") as inf:                                      
            dict1 = eval(inf.read())                                           
        return dict1                                                           

    def sing_image_populate(self, dectris_object_filename):
        dec_image = pickle_unjar(dectris_object_filename)
        self.populate_3D(dec_image)
        self.TIME = dec_image.TIME
        self.TEMP = dec_image.TEMP
        try:
            self.MAG_FIELD = dec_image.MAG_FIELD
        except:
            KeyError
        #self.temp = dec_image.TEMP

    def multi_image_populate(self, dectris_objects_filenames):                  
        self.time = []
        new_start_t = time.time()                                                  
        for index, filename in enumerate(dectris_objects_filenames):               
            dec_image = pickle_unjar(filename)                                     
            self.populate_3D(dec_image)                                         
            self.time.append(dec_image.time)
            progress_bar(index + 1, len(dectris_objects_filenames), new_start_t)   
                                                                               
    def populate_3D(self, dec_image):                                          
        for coordinate_1 in range(np.shape(dec_image.data)[1]):                
            for coordinate_2 in range(np.shape(dec_image.data)[0]):            
#                Q_index_1 = self.find_q_index(dec_image.Q_coords[coordinate_1][coordinate_2])
                #Q_index = self.find_q_linear(dec_image.Q_coords[coordinate_1][coordinate_2])
                Q_index = find_q_linear_fast(*dec_image.Q_coords[coordinate_1][coordinate_2], *self.find_q_lin)
                                                                               
                if Q_index[0] == False:                                        
                    pass                                                       
                else:                                                          
                    self.add_point_3D(Q_index[1], dec_image.data[coordinate_2][coordinate_1])
                                                                               
    def find_q_linear(self, coord):                                            
        qx = coord[0]                                                          
        qy = coord[1]                                                          
        qz = coord[2]                                                          
        if qx < self.QX_MIN \
                or qx >= self.QX_MAX \
                or qy < self.QY_MIN \
                or qy >= self.QY_MAX \
                or qz < self.QZ_MIN \
                or qz >= self.QZ_MAX:                                          
            return (False,[0,0,0])                                             
        else:                                                                  
            return(True, [math.floor((qx-self.QX_MIN)/self.QX_GRAD)+1,         
                          math.floor((qy-self.QY_MIN)/self.QY_GRAD)+1,         
                          math.floor((qz-self.QZ_MIN)/self.QZ_GRAD)+1])        
                                                                               
                                                                               
    def add_point_3D(self, indexes, data_point):                               
        self.data[indexes[0], indexes[1], indexes[2]] = \
                self.data[indexes[0], indexes[1], indexes[2]] + data_point     
        self.q_idx[indexes[0], indexes[1], indexes[2]] += 1                    
                                                                               
    def plot(self):                                                            
        if self.SPACE_3D == True:                                              
            self.plot_3d()                                                     

    def plot_sing(self):
        print("input min value for plot, bare in mind its already a log plot")
        vmin = input()
        print("input max values for plot")
        vmax = input()
        vmin = int(vmin)
        vmax = int(vmax)
        fig = plt.figure(figsize = (5,5))                                     
        ax = plt.subplot(111, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
        fig, ax= self.plot_2d(fig, ax, np.mean(self.data, axis=1),\
                self.q_x, self.q_z,\
                xlabel = r"\rm Q$_z \ (\rm\AA^{-1})$",                         
                ylabel = r"\rm Q$_x \ (\rm\AA^{-1})$", vmin = vmin, vmax = vmax)                         
        plt.tight_layout()                                                     
        plt.show()
        fig = plt.figure(figsize = (5,5))                                     
        ax = plt.subplot(111, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
        fig, ax = self.plot_2d(fig, ax, np.mean(self.data, axis=0),\
                self.q_y, self.q_z,\
                xlabel = r"\rm Q$_z \ (\rm\AA^{-1})$",                         
                ylabel = r"\rm Q$_y \ (\rm\AA^{-1})$", vmin = vmin, vmax = vmax)                         
        plt.tight_layout()                                                     
        plt.show()
        fig = plt.figure(figsize = (5,5))                                     
        ax = plt.subplot(111, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
        fig, ax = self.plot_2d(fig, ax, np.mean(self.data, axis=2),\
                self.q_x, self.q_y,\
                xlabel = r"\rm Q$_y \ (\rm\AA^{-1})$",                         
                ylabel = r"\rm Q$_x \ (\rm\AA^{-1})$", vmin = vmin, vmax = vmax)                         

        plt.tight_layout()                                                     
        plt.show()
                                                                               
    def plot_2d(self, fig, ax, data, xaxis, yaxis, xlabel, ylabel,\
            vmin = 1, vmax = 6):                                               
        ax.contourf(np.log10(data), origin='lower',                            
                       extent=[min(yaxis), max(yaxis), min(xaxis),max(xaxis)], 
            levels = 256, vmin = vmin, vmax = vmax,                            
                         cmap=cbar)                                            
        ax.contour(np.log10(data), origin='lower',                             
                       extent=[min(yaxis), max(yaxis), min(xaxis),max(xaxis)], 
            levels = 32, colors = "#00305d", linewidths = 0.1, vmin = vmax)    
        ax.set_xlabel(xlabel)                                                  
        ax.set_ylabel(ylabel)                                                  
        return fig, ax                                                         

    def plot_3d(self):                                                         
        fig = plt.figure(figsize = (12,5))                                     
        ax_xz = plt.subplot(131, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
        ax_yz = plt.subplot(132, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
        ax_xy = plt.subplot(133, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
                                                                               
#        fig = plt.figure(figsize = (4,4))                                     
#        ax = plt.subplot(111)                                                 
        fig, ax_xz = self.plot_2d(fig, ax_xz, np.mean(self.data, axis=1),\
                self.q_x, self.q_z,\
                xlabel = r"\rm Q$_z \ (\rm\AA^{-1})$",                         
                ylabel = r"\rm Q$_x \ (\rm\AA^{-1})$")                         
        fig, ax_yz = self.plot_2d(fig, ax_yz, np.mean(self.data, axis=0),\
                self.q_y, self.q_z,\
                xlabel = r"\rm Q$_z \ (\rm\AA^{-1})$",                         
                ylabel = r"\rm Q$_y \ (\rm\AA^{-1})$")                         
        fig, ax_xy = self.plot_2d(fig, ax_xy, np.mean(self.data, axis=2),\
                self.q_x, self.q_y,\
                xlabel = r"\rm Q$_y \ (\rm\AA^{-1})$",                         
                ylabel = r"\rm Q$_x \ (\rm\AA^{-1})$")                         
        plt.tight_layout()                                                     
        plt.savefig("local/images/"+str(self.scan_num) + "_3d.png", bbox_inches = "tight")
        plt.show()                                                             

    def plot_for_gif(self):
        fig = plt.figure(figsize = (12,6))                                     
        ax_xz = plt.subplot(131, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
        ax_yz = plt.subplot(132, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
        ax_xy = plt.subplot(133, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
                                                                               
#        fig = plt.figure(figsize = (4,4))                                     
#        ax = plt.subplot(111)                                                 
        fig, ax_xz = self.plot_2d(fig, ax_xz, np.mean(self.data, axis=1),\
                self.q_x, self.q_z,\
                xlabel = r"\rm Q$_z \ (\rm\AA^{-1})$",                         
                ylabel = r"\rm Q$_x \ (\rm\AA^{-1})$")                         
        fig, ax_yz = self.plot_2d(fig, ax_yz, np.mean(self.data, axis=0),\
                self.q_y, self.q_z,\
                xlabel = r"\rm Q$_z \ (\rm\AA^{-1})$",                         
                ylabel = r"\rm Q$_y \ (\rm\AA^{-1})$")                         
        fig, ax_xy = self.plot_2d(fig, ax_xy, np.mean(self.data, axis=2),\
                self.q_x, self.q_y,\
                xlabel = r"\rm Q$_y \ (\rm\AA^{-1})$",                         
                ylabel = r"\rm Q$_x \ (\rm\AA^{-1})$")                         
        plt.tight_layout()                                                     
        filename = "local/gif/images/"+str(self.scan_num) + "_" + str(self.TIME)+ "_gif.png"
        plt.savefig(filename, bbox_inches = "tight")
        plt.close()
        return filename



    def save(self):
        new_filename = existential_check(str(self.scan_num[0]) + "_new_3d_fill",
                                     ".pickle",
                                     "local/processed_files/")
        pickle_jar(new_filename, self)
                                                                               
                                                                               

