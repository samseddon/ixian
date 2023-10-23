import os
import inspect
#from pixel_selection import data_fill
import numpy as np
from rsm_object_plotter import file_checker,  rsm_plot, slicer_and_dicer_3000, test_slicer
from object_approach import data_fill
##

## MAKE SURE DATA, WHEREVER IT IS, IS IN A FOLDER CALLED MATERIAL/DATA (IE magnetite/data) and also in magnetite is a folder called user_defined_parameters, as before

##
## NOTE test data entered here 
#file_reference = "aV2O3R1"
#directory = "/home/sseddon/Desktop/500GB/Data/XMaS/test_data/"
#output_folder = "/home/sseddon/Desktop/500GB/Data/XMaS/test_data/processed_files/"
#scan_num = [418,419]
## NOTE K test run data ends here

def main():
    # NOTE test data entered here 
    file_reference = "aV2O3R1"
    directory = "/home/sseddon/Documents/Local-data/test_data/"
    output_folder = directory + "processed_files/"
    scan_num = [418,419]
    # NOTE K test run data ends here

    directory="/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/"    
    output_folder = "/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/processed_files_july23/"
    file_reference = "MAG001"
   # #scan_nums = [[153],[154]]
    scan_num = [152]
    data_fill(directory, output_folder, file_reference, scan_num, create_files = True)
    file_index = [-1] #for most recent processed data just put -1, otherwise 0 will do the first or its index
    axis_1 = 'qx'
    axis_2 = 'qz'
    axis_3_limits = [1,240]
    ##
    f_1, f_name = file_checker(scan_num, file_index[0],directory + 'processed_files/')
    print('plotting', f_name)
    rsm_plot(f_1, f_name, axis_1, axis_2, axis_3_limits,output_folder+'images/', directory, scan_num)
    
    
    
    """ 
    All these steps are for running the code through every file that is there
    """
    """
    Uncomment me first! 
    
    file_check = os.listdir(directory+"data/")
    scan_nums = []
    for file in file_check:
        scan_nums.append(int(file.split("_")[-2]))
    unique_scans_master = list(np.unique(scan_nums))
    
    
    completed_files_all = os.listdir(output_folder)
    completed_files_actual = []
    for file in completed_files_all:
        if file.endswith(".pickle") == True:
            completed_files_actual.append(int(file.split("_")[0]))
    completed_files_actual = sorted(completed_files_actual)
    scans_left_to_run = [x for x in unique_scans_master if x not in completed_files_actual]
    scans_final = []
    for scan_num in scans_left_to_run:
        scans_final.append([scan_num])
    
    print(scans_final)
    
    
    
    
    
    for scan_num in scans_final:
        try:
            print(scan_num)
    #        data_fill(directory, output_folder, file_reference, scan_num, create_files = True)
            file_index = [-1] #for most recent processed data just put -1, otherwise 0 will do the first or its index
            axis_1 = 'qx'
            axis_2 = 'qz'
            axis_3_limits = [1,240]
            ##
            f_1, f_name = file_checker(scan_num, file_index[0],directory + 'processed_files_july23/')
            print('plotting', f_name)
            rsm_plot(f_1, f_name, axis_1, axis_2, axis_3_limits,output_folder+'images/', directory, scan_num)
            
            axis_1 = 'qy'
            axis_2 = 'qz'
            axis_3_limits = [1,240]
            #
            f_1, f_name = file_checker(scan_num, file_index[0],directory + 'processed_files_july23/')
            print('plotting', f_name)
            rsm_plot(f_1, f_name, axis_1, axis_2, axis_3_limits,output_folder+'images/', directory, scan_num)
            
            axis_1 = 'qx'
            axis_2 = 'qy'
            axis_3_limits = [1,240]
            #
            f_1, f_name = file_checker(scan_num, file_index[0],directory + 'processed_files_july23/')
            print('plotting', f_name)
            rsm_plot(f_1, f_name, axis_1, axis_2, axis_3_limits,output_folder+'images/', directory, scan_num)
    
        except KeyboardInterrupt:
            break
    
        except:
            continue
    """
if __name__ == "__main__":
    main()
    
    
    
    #_____________SAM adds plotting code here __________#
    #Q_vec = np.array([2.66,0.0,5.41])
    #test_slicer(Q_vec)
    #scan_num = [52,182,102]
    #for i in range(len(scan_num)):
    #    if i == 6: 
    #        pass
    #    else:
    #        print(scan_num[i])
    #        f_1, f_name = file_checker(scan_num[i], file_index[0],directory + 'processed_files/')
    #        slicer_and_dicer_3000(f_1, f_name, Q_vec, output_folder, scan_num[i])
    
    #f_1, f_name = file_checker(scan_num[1], file_index[0],directory + 'processed_files/')
    #slicer_and_dicer_3000(f_1, f_name, Q_vec, output_folder, scan_num[1])
