from equations import calc_qx, calc_qy, calc_qz
import os
import fabio
import numpy as np
from spot_dict import spot_dict

# THESE PARAMETERS ARE IN THE FUNCTION AS THE NAMES LISTED
Directory = "/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/data"
FileStart = "MAG001"
scan_num = [152]
# -----------------------------------------------------------


# Make a list of all files in the target directory
files = os.listdir(Directory)

# Filter it down to only master files
master_files = [m for m in files if m.startswith(FileStart) and m.endswith(".edf")]

# Filter list down to only master files of a certain sample with the target scan number
scan_masters = [c for c in master_files if int(c.split("_")[-2]) in scan_num]
two_theta = []
omega = []
phi = []

for f in scan_masters[:]:
    f = fabio.open(os.path.join(Directory, f)) 
    two_theta.append(float(f.header.get("motor_pos").split(" ")[0]))# TWO THETA VALUE
    omega.append(float(f.header.get("motor_pos").split(" ")[1])) # OMEGA VALUE
    # chi = float(f.header.get("motor_pos").split(" ")[2])  # CHI VALUE
    phi.append(float(f.header.get("motor_pos").split(" ")[3]))

Wavelength = 1 
two_theta = np.average(two_theta)
theta1 = min(omega)
theta2 = max(omega)
#theta1 = 20.2
#theta2 = 22.2
y_angle = 0 
phi = 0 

if phi == 0:
    omega_offset = 1.6923
else:
    omega_offset = 3.8765
theta1 = theta1+omega_offset
theta2 = theta2+omega_offset

qx1 = calc_qx(two_theta, theta1, Wavelength, y_angle)
qx2 = calc_qx(two_theta, theta2, Wavelength, y_angle)
qy1 = calc_qy(two_theta, theta1, Wavelength, y_angle)
qy2 = calc_qy(two_theta, theta2, Wavelength, y_angle)
qz1 = calc_qz(two_theta, theta1, Wavelength)
qz2 = calc_qz(two_theta, theta2, Wavelength)

print('Scan', scan_num[0], 'is of spot', spot_dict[str(scan_num[0])])
print('Qx mean is', (qx1+qx2)/2, 'with range', abs(qx2-qx1))
print('Qy mean is', (qy1+qy2)/2, 'with range', abs(qy2-qy1))
print('Qz mean is', (qz1+qz2)/2, 'with range', abs(qz2-qz1))

