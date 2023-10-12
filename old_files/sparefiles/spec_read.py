spec_file="/home/sseddon/Desktop/500GB/Data/XMaS/pyrrhotite/"\
          + "spec/Fe7S8_010_strain01.spec"

file1 = open(spec_file, 'r')
Lines = file1.readlines()
count = 0
spot_dict = {}

for line in Lines:
    if not line.strip():
        if 'loopscan' in Lines[count+1]:
            pass  
        elif 'timescan' in Lines[count+1]:
            pass  
        else:
            dict_count = {} 
            scan_num_and_type  = Lines[count+1].split('  ') 
            scan_num = scan_num_and_type[0][3:]
            scan_type = scan_num_and_type[1] + ' '\
                        + scan_num_and_type[2].split(' ')[0]
            count_mne = Lines[count+27].split('  ')
            count_pos = Lines[count+28].split(' ')
            for key in count_mne:                                                      
                for value in count_pos:                                                
                    dict_count[key] = value                                            
                    count_pos.remove(value)        
                    break  
            h = str(int(round(float(dict_count['H'])*2)/2))
            k = str(int(round(float(dict_count['K'])*2)/2))
            l = str(int(round(float(dict_count['L'])*2)/2))
            hkl = h + k + l
            
            spot_dict[scan_num] = [hkl,scan_type]
    count += 1
with open("/home/sseddon/Desktop/500GB/Data/XMaS/pyrrhotite/"
          + "user_defined_parameters/spot_dict.txt",'w') as inf:
      inf.writelines(str(spot_dict))
