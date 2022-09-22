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
         'nr_pts_x' : 150,  # Number of points along each Qx
         'nr_pts_y' : 150,  # Number of points along each Qy
         'nr_pts_z' : 150,  # Number of points along each QZ
         'number_x' : 149,  #
         'number_y' : 149,
         'number_z' : 149
         }
 # Total number of rows in camera image (including dead pixels)
  # Total number of Columns in camera image

# Defines the Dead regions on the Pilatus where the chips are connected - 300k so only in rows!   #######

sat_pix = 1 * 10 ** 7 # depends on the counting mode (this is not for accumulate, silly didier)


Zero_Pixel_Hor = 252  # Horizontal Pixel corresponding to "zero" for the experiment
Zero_Pixel_Ver = 307  # Vertical Pixel corresponding to "zero" for the experiment

wavelength = 1.0  # Experimetnal wavelength in Angstroms

# Defining the reciprocal space volume to populate from the raw data (200x200x200 is about max)
nr_same = 150
nr_pts_x = nr_same  # Number of points along each Qx
nr_pts_y = nr_same  # Number of points along each Qy
nr_pts_z = nr_same  # Number of points along each Qz

number_x = nr_pts_x - 1  #
number_y = nr_pts_y - 1
number_z = nr_pts_z - 1
