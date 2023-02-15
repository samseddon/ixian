import matplotlib.pylab as plt
from scipy.interpolate import griddata
import numpy as np
import math
from colour_bar import code_TUD_cbar as cbar
plt.style.use("seddon_TUD")

def twodplot(file, filename, axis_1, axis_2, axis_3_limits, output_folder):
    x = []
    y = []
    z = []
    cts = []
    points = []
    
    for _ in range(len(file)):

        points.append((file[_][0],file[_][2]))
        x.append(file[_][0])
        y.append(file[_][1])
        z.append(file[_][2])
        cts.append(file[_][3])
    
    grid_x, grid_y = np.mgrid[min(x):max(x):128j, min(z):max(z):128j]
    
    x_min = min(x)
    x_max = max(x)
    y_min = min(z)
    y_max = max(z)
    
    XY = griddata(points,np.log(cts),(grid_x, grid_y), method = 'nearest')
    for i in range(len(grid_x)):
        for j in range(len(grid_y)):
            if math.isnan(XY[j,i]) or math.isinf(XY[j,i]):
                XY[j,i] = 0
    fig, ax = plt.subplots(figsize = (5,4), dpi = 128)
    ax.contourf(XY.T, cmap = cbar, levels = 64)
    ax.contour(XY.T, levels = 16, colors = 'C3', linewidths = 0.1)
   # xtick_pos = [16,64,112]
   # ax.set_xticks(xtick_pos)
   # ax.set_xticklabels([round(x_min+(x_max-x_min)*(xtick_pos[0]/len(grid_x)),3),
   #                     round(x_min+(x_max-x_min)*(xtick_pos[1]/len(grid_x)),3),
   #                     round(x_min+(x_max-x_min)*(xtick_pos[2]/len(grid_x)),3)])
   # ax.set_xlabel(r"\rm Q$_z"+r"$ $(\rm\AA^{-1})$")
   # ytick_pos = [16,64,112]
   # ax.set_yticks(ytick_pos)
   # ax.set_yticklabels([round(y_min+(y_max-y_min)*(ytick_pos[0]/len(grid_y)),3),
   #                     round(y_min+(y_max-y_min)*(ytick_pos[1]/len(grid_y)),3),
   #                     round(y_min+(y_max-y_min)*(ytick_pos[2]/len(grid_y)),3)])
   # ax.set_ylabel(r"\rm Q$_x"+r"$ $(\rm\AA^{-1})$")
   # #plt.contourf(np.log10(z), extent = [np.min(x),np.max(x),np.min(y),np.max(y)], levels=64)
    plt.show()       
