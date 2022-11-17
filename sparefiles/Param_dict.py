# Global input parameters
param = {'ShotHeight': 619, # Total number of rows in camera image (including dead pixels)
         'ShotHeight_y': 487,# Total number of Columns in camera image
         'DeadRow1' : 195, # Defines the Dead regions on the Pilatus where the chips are connected - 300k so only in rows!   #######
         'DeadRow2' : 211,
         'DeadRow3' : 407,
         'DeadRow4' : 423,
         'sat_pix' : 1 * 10 ** 7,
         'Zero_Pixel_Hor' : 252,  # Horizontal Pixel corresponding to "zero" for the experiment
         'Zero_Pixel_Ver' : 307, # Vertical Pixel corresponding to "zero" for the experiment
         'wavelength' : 1.0,
         'nr_pts_x' : 80,  # Number of points along each Qx
         'nr_pts_y' : 80,  # Number of points along each Qy
         'nr_pts_z' : 80,  # Number of points along each QZ
         'chi_offset' : 0   #neccesary for the equations but typically zero
         }
