from reciprocal_converter import data_fill
from rsm_plot_pickle import file_checker,  rsm_plot
##

## MAKE SURE DATA, WHEREVER IT IS, IS IN A FOLDER CALLED MATERIAL/DATA (IE magnetite/data) and also in magnetite is a folder called user_defined_parameters, as before

##

directory="/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/data"    
output_folder = "/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/processed_files/"
file_reference = "MAG001"
scan_num = [239]

data_fill(directory,output_folder,file_reference,scan_num)
#data_fill(directory,output_folder,file_reference,[199])

file_index = [-1] #for most recent processed data just put -1, otherwise 0 will do the first or its index
axis_1 = 'qx'
axis_2 = 'qz'
axis_3_limits= [0,-1]
axis_3_limits = [1,140]


f_1, f_name = file_checker(scan_num[0], file_index[0],output_folder)
rsm_plot(f_1, f_name, axis_1, axis_2, axis_3_limits,output_folder+'images')
#f_1, f_name = file_checker(sc
