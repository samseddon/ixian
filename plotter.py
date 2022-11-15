from Projection_plot import file_checker, rsm_plot
import os


DataFolder = "\\home\\sseddon\\Desktop\\500GB\\Data\\XMaS\\magnetite\\data" 
             
# Path to output files:
#output_folder = "\\home\\sseddon\\Desktop\\500GB\\Data\\XMaS\\magnetite\\processed_files"
output_folder = "/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/processed_files"

# Enter start of scan name:
FileNameStart = "MAG001"




scan_num = [152]
file_index = [-1]
# plot x vs z, y vs z, or x vs y
axis_1 = 'x'
axis_2 = 'y'
#axis_3_limits= [0,-1]
axis_3_limits = [1,150]


f_1, f_name = file_checker(scan_num[0], file_index[0], output_folder+'/')
rsm_plot(f_1, f_name, 'x','z', axis_3_limits, output_folder+'/Images/')

