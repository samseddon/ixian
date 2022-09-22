from colorama import Fore

def progress_bar(progress, total): #here progress is an integer for a loop say, and total is a length of the loop
    percent = 100 * (progress/float(total))
    bar = '█' * int(percent) + '-' *(100-int(percent))
    print(f"\r|{bar}| {percent:.2f}%", end = "\r")
    if progress==total:
        print(Fore.GREEN + f"\r|{bar}| {percent:.2f}%", end = "\n\n")

def existential_check(o_f, f_s, folder): #file name, file suffix, destination folder
    f_n = str(o_f + f_s)
    while os.path.exists(folder + f_n):

        if len(f_n) > len(o_f + f_s):
            filenum = int(f_n[(len(o_f + f_s) - 3):-len(f_s)]) + 1 
            f_n = f_n[:(len(o_f + f_s) - 3)] + str(filenum) + f_n[-len(f_s):]
        else:
            f_n = f_n[:-len(f_s)] + '_1' + f_n[-len(f_s):]
    return folder + f_n 



