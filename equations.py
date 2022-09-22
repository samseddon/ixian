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

