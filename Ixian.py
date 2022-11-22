from reciprocal_converter import data_fill
##

## MAKE SURE DATA, WHEREVER IT IS, IS IN A FOLDER CALLED MATERIAL/DATA (IE magnetite/data) and also in magnetite is a folder called user_defined_parameters, as before

##

directory="/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/data"    
output_folder = "/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/processed_files/"
file_reference = "MAG001"
scan_num = [180]

data_fill(directory,output_folder,file_reference,scan_num)

