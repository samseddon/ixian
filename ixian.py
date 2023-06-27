#from rsm_plot_pickle import file_checker, rsm_plot
#from pixel_selection import data_fill
import numpy as np
from rsm_object_plotter import file_checker,  rsm_plot, slicer_and_dicer_3000, test_slicer
from object_approach import data_fill
from interplot_2 import twod_plot
##

## MAKE SURE DATA, WHEREVER IT IS, IS IN A FOLDER CALLED MATERIAL/DATA (IE magnetite/data) and also in magnetite is a folder called user_defined_parameters, as before

##

directory="/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/"    
output_folder = "/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/processed_files/"
file_reference = "MAG001"


## NOTE test data entered here 
#file_reference = "aV2O3R1"
#directory = "/home/sseddon/Desktop/500GB/Data/XMaS/test_data/"
#output_folder = "/home/sseddon/Desktop/500GB/Data/XMaS/test_data/processed_files/"
#scan_num = [418,419]


#data_fill(directory, output_folder, file_reference, scan_num, create_files = True)




file_index = [-1] #for most recent processed data just put -1, otherwise 0 will do the first or its index
#axis_1 = 'qx'
#axis_2 = 'qz'
#axis_3_limits = [1,240]
##
#f_1, f_name = file_checker(scan_num[0], file_index[0],directory + 'processed_files/')
#print('plotting', f_name)
#rsm_plot(f_1, f_name, axis_1, axis_2, axis_3_limits,output_folder+'images')
#
#axis_1 = 'qy'
#axis_2 = 'qz'
#axis_3_limits = [1,240]
##
#f_1, f_name = file_checker(scan_num[0], file_index[0],directory + 'processed_files/')
#print('plotting', f_name)
#rsm_plot(f_1, f_name, axis_1, axis_2, axis_3_limits,output_folder+'images')
#
#axis_1 = 'qx'
#axis_2 = 'qy'
#axis_3_limits = [1,240]
##
#f_1, f_name = file_checker(scan_num[0], file_index[0],directory + 'processed_files/')
#print('plotting', f_name)
#rsm_plot(f_1, f_name, axis_1, axis_2, axis_3_limits,output_folder+'images')


#_____________SAM adds plotting code here __________#
Q_vec = np.array([2.66,0.0,5.41])
#test_slicer(Q_vec)
scan_num = [52,182,102]
for i in range(len(scan_num)):
    if i == 6: 
        pass
    else:
        print(scan_num[i])
        f_1, f_name = file_checker(scan_num[i], file_index[0],directory + 'processed_files/')
        slicer_and_dicer_3000(f_1, f_name, Q_vec, output_folder, scan_num[i])

#f_1, f_name = file_checker(scan_num[1], file_index[0],directory + 'processed_files/')
#slicer_and_dicer_3000(f_1, f_name, Q_vec, output_folder, scan_num[1])
