import pickle
import numpy as np
import matplotlib.pyplot as plt
from colour_bar import code_TUD_cbar as cbar                                   
plt.style.use("seddon_TUD")        


with open('test_file.pickle', 'rb') as handle:
    read_file = pickle.load(handle)
print(read_file.Q_x[100])



fig, ax = plt.subplots(figsize = (5,4), dpi = 128)
img = ax.imshow(read_file.data[100], cmap = cbar)
#img = ax.contourf(read_file.data[100], )
plt.show()

