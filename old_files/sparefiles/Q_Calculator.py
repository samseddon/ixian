from equations import calc_qx, calc_qy, calc_qz
import os
import fabio
import numpy as np

# THESE PARAMETERS ARE IN THE FUNCTION AS THE NAMES LISTED
directory = "/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/data"
FileStart = "MAG001"
scan_num = []
scan_num.append(int(input('Enter scan number\n'))) #try 306 and 291
# -----------------------------------------------------------
with open(directory[:-4]+'user_defined_parameters/spot_dict.txt','r') as inf:
        spot_dict = eval(inf.read())

# Make a list of all files in the target directory
files = os.listdir(directory)

# Filter it down to only master files
master_files = [m for m in files if m.startswith(FileStart) and m.endswith(".edf")]

# Filter list down to only master files of a certain sample with the target scan number
scan_masters = [c for c in master_files if int(c.split("_")[-2]) in scan_num]
two_theta = []
omega = []
phi = []

for f in scan_masters[:]:
    f = fabio.open(os.path.join(directory, f)) 
    two_theta.append(float(f.header.get("motor_pos").split(" ")[0]))# TWO THETA VALUE
    omega.append(float(f.header.get("motor_pos").split(" ")[1])) # OMEGA VALUE
    # chi = float(f.header.get("motor_pos").split(" ")[2])  # CHI VALUE
    phi.append(float(f.header.get("motor_pos").split(" ")[3]))
print(phi)
Wavelength = 1 
two_theta = np.average(two_theta)
theta1 = min(omega)
theta2 = max(omega)
#theta1 = 20.2
#theta2 = 22.2
y_angle = 0 
if phi[-1] > - 50: 
    omega_offset = 1.6923
else:
        omega_offset = 3.8765

#if phi == 0:
#    omega_offset = 1.6923
#else:
#    omega_offset = 3.8765
theta1 = theta1+omega_offset
theta2 = theta2+omega_offset

qx1 = calc_qx(two_theta, theta1, Wavelength, y_angle)
qx2 = calc_qx(two_theta, theta2, Wavelength, y_angle)
qy1 = calc_qy(two_theta, theta1, Wavelength, y_angle)
qy2 = calc_qy(two_theta, theta2, Wavelength, y_angle)
qz1 = calc_qz(two_theta, theta1, Wavelength)
qz2 = calc_qz(two_theta, theta2, Wavelength)

print('Scan', scan_num[0], 'is of spot', spot_dict[str(scan_num[0])])
print('Qx lowlim is', (qx1+qx2)/2-abs(qx2-qx1), 'Qx highlim is', (qx1+qx2)/2+abs(qx2-qx1))
print('Qy lowlim is', (qx1+qx2)/2-abs(qx2-qx1), 'Qy highlim is', (qy1+qy2)/2+abs(qx2-qx1))
print('Qz lowlim is', (qz1+qz2)/2-abs(qz2-qz1), 'Qz highlim is', (qz1+qz2)/2+abs(qz2-qz1))

output = {
 'qz_min' : (qz1+qz2)/2-abs(qz2-qz1),
 'qz_max' : (qz1+qz2)/2+abs(qz2-qz1), 
 'qx_min' : (qx1+qx2)/2-abs(qx2-qx1),
 'qx_max' : (qx1+qx2)/2+abs(qx2-qx1),
 'qy_min' : (qx1+qx2)/2-abs(qx2-qx1),
 'qy_max' : (qx1+qx2)/2+abs(qx2-qx1)}

filename = directory[:-4]+'user_defined_parameters/qlim/qlim_'+spot_dict[str(scan_num[0])]+'.txt' 
print('Created file',filename)
if os.path.exists(filename) == True and input('Overwrite existing file, [y] or n?\n') != 'y':
        pass
else:
   with open(filename,'w') as inf:
      inf.write(str(output)) 





