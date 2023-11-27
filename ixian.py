import os
import inspect
#from pixel_selection import data_fill
import numpy as np
from matplotlib.cm import ScalarMappable
from rsm_object_plotter import file_checker,  rsm_plot, slicer_and_dicer_3000, test_slicer, threeD_rsm_plot
from object_approach import omega_scan
import matplotlib.pyplot as plt

def main():
    
    directory="/home/sseddon/Downloads/ixian_testdata/"
    file_reference = "MAG001"
    scan_num = [194]

    if os.path.exists(directory + "processed_files/") ==  False:
        os.mkdir(directory + "processed_files/")
    if os.path.exists("local/temp/") ==  False:
        os.mkdir("local/")
        os.mkdir("local/temp")  
    if os.path.exists("local/qlim/") ==  False:
        os.mkdir("local/qlim")  


    omega_scan(directory, file_reference, scan_num, create_files = True)
    

#    fig, (ax_xz, ax_yz, ax_xy) = plt.subplots(nrows = 1, ncols = 3, aspect = "equal", figsize = (12,5))
    fig = plt.figure(figsize = (12,5))
    ax_xz = plt.subplot(131, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
    ax_yz = plt.subplot(132, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')
    ax_xy = plt.subplot(133, aspect = "equal")#, autoscale_on=False)#, adjustable='box-forced')

    


    file_index = [-1] #for most recent processed data just put -1, otherwise 0 will do the first or its index
    axis_1 = 'qx'
    axis_2 = 'qz'
    axis_3_limits = [1,240]
    ##
    f_1, f_name = file_checker(scan_num, file_index[0],directory + 'processed_files/')
    print('plotting', f_name)
    vmin = 1
    vmax = 6
    fig, ax_xz = rsm_plot(f_1, f_name, "qx", "qz", axis_3_limits, directory + "processed_files/"+'images/', directory, scan_num, fig, ax_xz)
    fig, ax_yz = rsm_plot(f_1, f_name, "qy", "qz", axis_3_limits, directory + "processed_files/"+'images/', directory, scan_num, fig, ax_yz)
    fig, ax_xy = rsm_plot(f_1, f_name, "qx", "qy", axis_3_limits, directory + "processed_files/"+'images/', directory, scan_num, fig, ax_xy)
#    ax_xz.axes.set_aspect('equal', adjustable='box')
#    ax_yz.axes.set_aspect("equal", adjustable = "box")
#    ax_xy.axes.set_aspect("equal", adjustable = "box")
    plt.tight_layout()

   # 
   # cobar = fig.colorbar(
   #             ScalarMappable(norm=img.norm, cmap=img.cmap),
   #             ticks=range(vmin, vmax+1, 1))
   # #cobar.ax_xy.tick_params(size=0)
   # cobar.set_label('Log(Intensity) (A. U.)', rotation=270, labelpad = 30)
    plt.show()
    #threeD_rsm_plot(f_1)
#    axis_1 = 'qy'
#    axis_2 = 'qz'
#    axis_3_limits = [1,240]
#    ##
#    f_1, f_name = file_checker(scan_num, file_index[0],directory + 'processed_files/')
#    print('plotting', f_name)
#    rsm_plot(f_1, f_name, axis_1, axis_2, axis_3_limits,output_folder+'images/', directory, scan_num)
#    ##
#    axis_1 = 'qx'
#    axis_2 = 'qy'
#    axis_3_limits = [1,240]
#    f_1, f_name = file_checker(scan_num, file_index[0],directory + 'processed_files/')
#    print('plotting', f_name)
#    rsm_plot(f_1, f_name, axis_1, axis_2, axis_3_limits,output_folder+'images/', directory, scan_num)
    
    
    
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
