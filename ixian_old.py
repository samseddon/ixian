#from rsm_plot_pickle import file_checker, rsm_plot                            
#from pixel_selection import data_fill                                         
                                                                               
from rsm_plot_pickle import file_checker,  rsm_plot                         
from pixel_selection import data_fill                                          
##                                                                             
                                                                               
## MAKE SURE DATA, WHEREVER IT IS, IS IN A FOLDER CALLED MATERIAL/DATA (IE magnetite/data) and also in magnetite is a folder called user_defined_parameters, as before
                                                                               
##                                                                             
                                                                               
                                                                               
# KRISTINA - NOTE THAT DIRECTORY IS JUST TO MAGNETITE FOLDER IE YOUR DATA FILE. 
directory="/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/"                   
output_folder = "/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/processed_files/"
file_reference = "MAG001"                                                      
scan_num = [[195]]                                                             
                                                                               
                                                                               
                                                                               
for i in range(len(scan_num)):                                                 
    data_fill(directory, output_folder, file_reference, scan_num[i], create_files = True)
                                                                               
                                                                               
                                                                               
file_index = [-1] #for most recent processed data just put -1, otherwise 0 will do the first or its index
axis_1 = 'qy'                                                                  
axis_2 = 'qz'                                                                  
axis_3_limits = [1,240]                                                        
                                                                               
f_1, f_name = file_checker(scan_num[0], file_index[0],directory + 'processed_files/')
print('plotting', f_name)                                                      
rsm_plot(f_1, f_name, axis_1, axis_2, axis_3_limits,output_folder+'images')       
                                                                               
file_index = [-1] #for most recent processed data just put -1, otherwise 0 will do the first or its index
axis_1 = 'qx'                                                                  
axis_2 = 'qz'                                                                  
axis_3_limits = [1,240]                                                        
                                                                               
f_1, f_name = file_checker(scan_num[0], file_index[0],directory + 'processed_files/')
print('plotting', f_name)                                                      
rsm_plot(f_1, f_name, axis_1, axis_2, axis_3_limits,output_folder+'images')       
                                                                               
file_index = [-1] #for most recent processed data just put -1, otherwise 0 will do the first or its index
axis_1 = 'qx'                                                                  
axis_2 = 'qy'                                                                  
axis_3_limits = [1,240]                                                        
                                                                               
f_1, f_name = file_checker(scan_num[0], file_index[0],directory + 'processed_files/')
print('plotting', f_name)                                                      
rsm_plot(f_1, f_name, axis_1, axis_2, axis_3_limits,output_folder+'images')       

