import os
import fabio
import time
import Experimental_Details
from useful_functions import progress_bar


with open('paths.txt') as f:
    lines = f.readlines()
print(lines[0])

directory="/home/sseddon/Desktop/500GB/Data/XMaS/magnetite"
files_location = os.listdir("/home/sseddon/Desktop/500GB/Data/XMaS/magnetite")
file_reference = "MAG001"
scan_num = [152]
# Selecting the .edf files from a messy folder full of all Data scans
master_files = [m for m in files_location if m.startswith(file_reference) and m.endswith(".edf")]

# Filter list down to only master files of a certain sample with the target scan number
master_files = [c for c in master_files if int(c.split("_")[-2]) in scan_num]
i=0
for f in master_files:
    f = fabio.open(os.path.join(directory, f))
    #Stripping the header 
    attenuator = float(f.header.get("counter_pos").split(" ")[24])
    trans = float(f.header.get("counter_pos").split(" ")[23])
    io = float(f.header.get("counter_pos").split(" ")[6])
    two_theta = float(f.header.get("motor_pos").split(" ")[0])  # TWO THETA VALUE
    omega = float(f.header.get("motor_pos").split(" ")[1])  # OMEGA VALUE
    # chi = float(f.header.get("motor_pos").split(" ")[2])  # CHI VALUE
    phi = float(f.header.get("motor_pos").split(" ")[3])  # PHI VALUE
    time.sleep(0.02)
    progress_bar(i+1,len(master_files))
    i=i+1
