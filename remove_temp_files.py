import os

temp_files = []
for file in os.listdir("local/temp"):
    temp_files.append(file)


print("Deleting " + str(len(os.listdir("local/temp"))) + \
        " files, are you sure? [y]/n")
prompt = input()
if prompt == "y" or prompt == "Y":
    for file in temp_files:
        os.remove("local/temp/" + file)
    print("Files deleted")
else:
    print("tempory files kept")
