import os
import fabio
import matplotlib.pyplot as plt
import time
from colour_bar import code_TUD_cbar as cbar
from equations import progress_bar


def plot_and_save(directory,output_folder,file_reference,scan_num,create_files):  
    start_t = time.time()
    if os.path.exists(output_folder) == False:
        os.makedirs(output_folder)

    files_location = os.listdir(directory+'data/')                             
    master_files = [m for m in files_location \
                    if m.startswith(file_reference) \
                    and m.endswith(".edf")]                                    
    master_files = [c for c in master_files \
                    if int(c.split("_")[-2]) \
                    in scan_num]                                               
    for image_number in range(len(master_files)):                              
        file = fabio.open(os.path.join(directory + "data/",               
                                       master_files[image_number]))

        data = file.data
        fig, ax = plt.subplots(figsize = (5,4))
        ax.imshow(data, cmap = cbar)
        plt.savefig(output_folder+str(image_number)+".png")
        plt.close()
        progress_bar(image_number+1,len(master_files),start_t)




if __name__ == "__main__":
    directory="/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/"
    output_folder = "/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/crazy_plotr" 
    file_reference = "MAG001"
    scan_num = [167]
    output_folder += "/" + str(scan_num) + "/"
    plot_and_save(directory, output_folder, file_reference, scan_num, create_files = True)
