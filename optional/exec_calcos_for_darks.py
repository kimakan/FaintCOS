"""This module uses CALCOS to reduce all rawtag files in the folder.

The reduced files (corrtags) are stored in the folder "reduced". 
All obsolete files (flt_a.fits, counts_a.fits etc.) will be 
automatically removed.

"""

import subprocess, sys, os

all_files = os.listdir(".")
files = []
for f in all_files:
    filename = f.split("_")
    if (len(filename) > 2):
        if (f.split("_")[1] == "rawtag" and f.split("_")[2] == "a.fits"):
            files.append(f)
for f in files:
    process = subprocess.Popen("calcos -o reduced " + f + " > log.txt", shell=True)
    print("Working on " + f)
    process.wait()
    os.system("rm -r reduced/*counts_a.fits")
    os.system("rm -r reduced/*counts_b.fits")
    os.system("rm -r reduced/*flt_a.fits")
    os.system("rm -r reduced/*flt_b.fits")
    os.system("rm -r reduced/*.tra")
