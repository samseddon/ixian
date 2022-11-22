import pickle
import numpy as np
import os

def existential_check(o_f, f_s, folder): #file name, file suffix, destination folder
    f_n = str(o_f + f_s)
    while os.path.exists(folder + f_n):
        if len(f_n) > len(o_f + f_s):
            filenum = int(f_n[(len(o_f + f_s) - (len(f_s)-1)):-len(f_s)]) + 1
            f_n = f_n[:(len(o_f + f_s) - (len(f_s)-1))] + str(filenum) + f_n[-len(f_s):]
        else:
            f_n = f_n[:-len(f_s)] + '_1' + f_n[-len(f_s):]
    return folder + f_n

# Store data (serialize)
#with open('filename.pickle', 'wb') as handle:
#    pickle.dump(your_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Load data (deserialize)
with open('Comparison/182_3d_fill_5.pickle', 'rb') as handle:
    data_work = pickle.load(handle)
    
with open('Comparison/180_3d_fill_3.pickle', 'rb') as handle:
    data_not_work = pickle.load(handle)

print(data_not_work)
print(data_work)
#qx = []
#qy = []
#qz = []
#
#for i in range(40):
#    qx.append(i)
#    qy.append(i)
#    qz.append(i)
#qaxes = []
#qaxes.append(qx)
#qaxes.append(qy)
#qaxes.append(qz)
##
##with open('file.txt', 'w') as file:
##     file.write(str(your_data))
#export = {'qx':qx, 'qy':qy,'qz':qz,'data':data}
#orig_filename = '182_3d_fill'
#suffix = '.pickle'
#output_folder = ''
#with open(existential_check(orig_filename,suffix,output_folder), 'wb') as handle:
#    pickle.dump(export, handle, protocol=pickle.HIGHEST_PROTOCOL)

