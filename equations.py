import numpy as np

# Calculate qx from 2theta, eta, and wavelength
def calc_qx(twoTheta, eta, wvl, angle_y):
    qx = (2 * np.pi / wvl) * (np.cos(eta / 180 * np.pi) - np.cos(twoTheta / 180 * np.pi - eta / 180 * np.pi)
                              * np.cos(angle_y / 180 * np.pi))
    return qx


# Calculate qz from 2theta, eta, and wavelength
def calc_qz(twoTheta, eta, wvl):
    qz = (2 * np.pi / wvl) * (np.sin(eta / 180 * np.pi) + np.sin(twoTheta / 180 * np.pi - eta / 180 * np.pi))
    return qz


# Calculate qy from 2theta, eta, and wavelength
def calc_qy(twoTheta, eta, wvl, angle_y):
    qy = (2 * np.pi / wvl) * np.sin(angle_y / 180 * np.pi) \
         * np.cos(twoTheta / 180 * np.pi - eta / 180 * np.pi)
    return qy

def index_grid(Param_dict):  # This is the grid that keeps track of how mnay elements are being put into each Q voxel.    # Identical size to Q volume
    global idx_grid
    idx_grid = np.zeros((Param_dict['nr_pts_x'], Param_dict['nr_pts_y'], Param_dict['nr_pts_z']))
    return idx_grid

# Makes a list of 2thetas corresponding to each row of pixels (top down) in an image
def generate_two_thetas(TwoTheta0, Param_dict):
    two_thetas = []
    for x in range(Param_dict['ShotHeight']):
        two_theta = TwoTheta0 + (Param_dict['Zero_Pixel_Ver'] - float(x)) * 0.01126  # degrees per row determined from alignment calibration (pixel # vs. two theta)
        two_thetas.append(two_theta)
    return two_thetas



def qy_angle(Param_dict):
    qy_angles = []
    for x in range(Param_dict['ShotHeight_y']):
        angle = (Param_dict['Zero_Pixel_Hor'] - Param_dict['chi_offset'] - float(x)) * 0.01126
        # degrees per column determined from alignment calibration (pixel # vs. nu)
        qy_angles.append(angle)
    return qy_angles

def dead_pixel(f,param,two_theta_range):
    #  REMOVE THE DEAD PIXELS
    DeadRow1 = param['DeadRow1']
    # Defines the Dead regions on the Pilatus where the chips are connected - 300k so only in rows!   #######
    DeadRow2 = param['DeadRow2']
    DeadRow3 = param['DeadRow3']
    DeadRow4 = param['DeadRow4']
    f_new = np.concatenate((f.data[:DeadRow1 - 1, :],
                        f.data[DeadRow2 + 1:DeadRow3 - 1, :],
                        f.data[DeadRow4 + 1:, :]), axis=0)
    two_theta_range_new = np.concatenate((two_theta_range[:DeadRow1 - 1],
                                          two_theta_range[DeadRow2 + 1:DeadRow3 - 1],
                                          two_theta_range[DeadRow4 + 1:]))
    return f_new,two_theta_range_new


