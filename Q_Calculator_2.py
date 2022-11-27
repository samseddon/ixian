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
Wavelength = 1 
#two_theta = np.average(two_theta)
two_theta1 = min(two_theta)-3
two_theta2 = max(two_theta)+3
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
qzmaxleft = calc_qz(two_theta2,theta1,Wavelength)
qzminleft = calc_qz(two_theta1,theta1,Wavelength)
qzmaxright = calc_qz(two_theta2,theta2,Wavelength)
qzminright = calc_qz(two_theta1,theta1,Wavelength)
qzcentre = calc_qz(two_theta2, (theta1+theta2)/2,Wavelength)

qxmaxleft = calc_qx(two_theta2,theta1,Wavelength, y_angle+2.5)
qxminleft = calc_qx(two_theta1,theta1,Wavelength, y_angle+2.5)
qxmaxright = calc_qx(two_theta2,theta2,Wavelength, y_angle+2.5)
qxminright = calc_qx(two_theta1,theta2,Wavelength, y_angle+2.5)

qymin = calc_qy(two_theta2,theta1,Wavelength,y_angle-2.5)
qymax = calc_qy(two_theta2,theta1,Wavelength,y_angle+2.5)


qzmax = max(qzmaxleft,qzmaxright,qzcentre)
qzmin = min(qzminleft,qzminright)
qxmax = max(qxmaxright,qxmaxleft)
qxmin = min(qxminright,qxminleft)

qxave = (qxmax+qxmin)/2
qyave = (qymax+qymin)/2
qzave = (qzmax+qzmin)/2
qxrange = abs(qxmax-qxmin)/2
qyrange = abs(qymax-qymin)/2
qzrange = abs(qzmax-qzmin)/2

output = {
 'qz_min' : qzave-qzrange,
 'qz_max' : qzave+qzrange,
 'qx_min' : qxave-qxrange,
 'qx_max' : qxave+qxrange,
 'qy_min' : qyave-qyrange,
 'qy_max' : qyave+qyrange}

filename = directory[:-4]+'user_defined_parameters/qlim/qlim_'+spot_dict[str(scan_num[0])]+'.txt'
if os.path.exists(filename) == True and input('Overwrite existing file, [y] or n?\n') != 'y':
        pass
else:
   with open(filename,'w') as inf:
      inf.write(str(output))
   print('Created file',filename)

