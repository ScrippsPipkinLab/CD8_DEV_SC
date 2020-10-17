#!/usr/bin/env python3
import os
import glob
import sys

print("----------")
print("Check if files are up-to-date and move old versions into a 'old' folder")
print("----------")
print("Usage:")
print("python check_up-to-date.py path_to_folder_for_checking")
print("----------")
print("Naming convention for files:")
print("filex_YYYYMMDD")
print("filey_YYYYMMDDnew")
print("")
print("")
print("")
print("File:\tMost up to date version")
print("----------")


wkdir = sys.argv[1]

os.chdir(wkdir)

#--- Check if 'old' folder exist, create if not
if not os.path.exists("old"):
    os.mkdir("old")

#--- Find all files & file names without dates
files = glob.glob("*")
files_simp = [x.split("_")[0] for x in files]
files_simp = list(set(files_simp))


for i in files_simp:
    files_i = [x for x in files if "%s_"%i in x]
    if ("20" in "".join(files_i)): # Disregard files that do not have dates
        
        # Get date list & find the most recent date
        files_i_dates = [x.split(".")[-2].split("_")[-1].replace("new", "") for x in files_i]
        files_i_dates = [int(x) for x in list(set(files_i_dates))]
        files_i_dates.sort()
        use_date = str(files_i_dates[-1])
        if ("%snew"%use_date in "".join(files_i)):
            use_date = "%snew"%use_date
        print("%s:\t%s"%(i,use_date))
        
        # Move files that are not most recent dates into 'old' folder
        for j in files_i:
            if use_date not in j:
                os.rename(j, 'old/%s'%j)  
                #print('old/%s'%j)