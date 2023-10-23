import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
plt.style.use("seddon_TUD")

"""                                                                            
Created on Tue Jan 10 10:52:19 2023                                             
                                                                               
@author: samseddon                                                                
"""  

class Scan():
    def __init__(self):
        self.scan_num = 0
        self.scan_type = ""
        self.motor = ""
        self.scan_parameters = ""
        self.time = ""
        self.temp_A = 0.
        self.temp_B = 0.
        self.motor_nme = []
        self.motor_value = []
        self.iroi = []
        self.iroi2 = []
        self.cdet = []


    def set_scan_type(self, scan_type):
        self.scan_type = scan_type

    def set_motor(self, motor):
        self.motor = motor

    def set_scan_parameters(self, scan_parameters):
        self.scan_parameters = scan_parameters 

    def set_scan_num(self, scan_num):
        self.scan_num = int(scan_num)

    def set_motor_nme(self, array):
        self.motor_nme = array

    def set_motor_value(self, motor_value):
        self.motor_value.append(float(motor_value))
    
    def set_cdet(self, cdet):
        self.cdet.append(float(cdet))
    
    def set_iroi(self, iroi):
        self.iroi.append(float(iroi))
    
    def set_iroi2(self, iroi2):
        self.iroi2.append(float(iroi2))

def generate_spectre_list(spec_file_location):

    file1 = open(spec_file_location, 'r')
    Lines = file1.readlines()
    scans = []
    first_scan = True
    
    for line in Lines:
        if not line.strip():
            if first_scan == True:
                scan = Scan()
                first_scan = False
                pass
            else:
                scans.append(scan)
                scan = Scan()
                pass
        else:    
            line_array = line.strip().split(' ')
            line_id = line_array[0]
                
            if line_id == '#S':
                scan.set_scan_type(line_array[3])
                scan.set_scan_num(line_array[1])
                scan.set_motor(line_array[5])
            if line_id == '#L':
                scan.set_motor_nme(line[2:].strip().split('  '))
            if line_id[0] != "#":
                scan.set_motor_value(line_array[0])
                for i in range(1, len(scan.motor_nme)):
                        nme = scan.motor_nme[i]
                        if nme == 'cdet':
                            scan.set_cdet(line_array[i])
                        if nme == 'iroi':
                            scan.set_iroi(line_array[i])
                        if nme == 'iroi2':
                            scan.set_iroi2(line_array[i])
    return scans

if __name__ == "__main__":
    spec_file_location="/home/sseddon/Desktop/500GB/Data/XMaS/pyrrhotite/"\
              + "spec/Fe7S8_010_strain01.spec"
    scans = generate_spectre_list(spec_file_location)
#    file1 = open(spec_file, 'r')
#    Lines = file1.readlines()
#    scans = []
#    first_scan = True
#    
#    for line in Lines:
#        if not line.strip():
#            if first_scan == True:
#                scan = Scan()
#                first_scan = False
#                pass
#            else:
#                scans.append(scan)
#                scan = Scan()
#                pass
#        else:    
#            line_array = line.strip().split(' ')
#            line_id = line_array[0]
#                
#            if line_id == '#S':
#                scan.set_scan_type(line_array[3])
#                scan.set_scan_num(line_array[1])
#                scan.set_motor(line_array[5])
#            if line_id == '#L':
#                scan.set_motor_nme(line[2:].strip().split('  '))
#            if line_id[0] != "#":
#                scan.set_motor_value(line_array[0])
#                for i in range(1, len(scan.motor_nme)):
#                        nme = scan.motor_nme[i]
#                        if nme == 'cdet':
#                            scan.set_cdet(line_array[i])
#                        if nme == 'iroi':
#                            scan.set_iroi(line_array[i])
#                        if nme == 'iroi2':
#                            scan.set_iroi2(line_array[i])

    
    scan_num = 137
    
    x = []
    y = []
    key = []
    motor  = ""
    xticks = []
    scale_x = 1
    for scan in scans:
        if scan.scan_num == scan_num:
            motor = scan.motor
            x = scan.motor_value
            y.append(scan.cdet)
            key.append("c-det")
            y.append(scan.iroi2)
            key.append("roi-1")
    if motor == "amp1":
        scale_x = 0.02 # 0.02 if strain cell scan
        xlabel = "Applied Voltage (V)"
        xticks = [0, 10, 20, 30, 40]
        for i in range(len(xticks)):
            xticks[i] = xticks[i] * scale_x
    if motor == "eta":
        scale_x = 1
        xlabel = "Eta angle (deg)"

    #_________________________PLOT PARAMETERS_________________________________#
    
    title = 'Integrated Intensities'
    scale_y = 10**8
    ylabel = r"Intensity (counts $\times 10 ^" + str(int(np.log10(scale_y))) +"$)"
    yticks = [1, 2, 3, 4]

    #_________________________________________________________________________#
    for i in range(len(yticks)):
        yticks[i] = yticks[i] * scale_y
    
    fig, ax = plt.subplots()
    for _ in range(len(y)):
        ax.plot(x,y[_])
    
    plt.legend(key)
    if len(xticks) > 0:
        plt.xticks(xticks)
    if len(yticks) > 0:
        plt.yticks(yticks) 
    ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_x))
    ax.xaxis.set_major_formatter(ticks_x)
    ax.set_xlabel(xlabel)
    ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_y))
    ax.yaxis.set_major_formatter(ticks_y)
    ax.set_ylabel(ylabel)
    plt.show()

        
            
