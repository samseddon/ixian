import os

directory="/home/sseddon/Desktop/500GB/Data/XMaS/magnetite"

with open(directory + '/user_defined_parameters/spot_dict.txt','r') as inf:
        dict1 = eval(inf.read())
flipped = {}
with open(directory+'/user_defined_parameters/param/standard_param.txt', 'r') as inf: 
    output = eval(inf.read())

for key, value in dict1.items():
    if value not in flipped:
        flipped[value] = [key]
    else:
        flipped[value].append(key)

all_spots = list(flipped.keys())

for i in range(len(all_spots)):
    filename = directory+'/user_defined_parameters/param/param_'+all_spots[i]+'.txt'
    if os.path.exists(filename) == True and input('Overwrite existing '+all_spots[i]+' file, [y] or n?\n') != 'y':
        pass
    else: 
        with open(filename,'w') as inf:
            inf.write(str(output))
        print('Created file',filename)

