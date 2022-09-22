from Projection_plot import file_checker,rsm_plot
output_folder =  "/home/sseddon/Desktop/500GB/Data/XMaS/magnetite/processed_files/"
axis_3_limits = [1,150]
file_index = [-1]
f_1, f_name = file_checker([152], file_index[0], output_folder)
rsm_plot(f_1, f_name, 'x','z', axis_3_limits, output_folder+'Images\\')

