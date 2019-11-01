import glob, sys, os

prefix="/eos/uscms/store/user/zhangj/HGC/HDBSCAN/"

files=glob.glob(prefix+"*v2*")

for fil in files:
    cmd="mv "+fil+" "+fil.replace("nonMIP", "")
    os.system(cmd)
