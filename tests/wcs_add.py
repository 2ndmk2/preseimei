
from preseimei import darkflat
from preseimei import install_indexfiles
from preseimei import seimei_ana
import argparse
import pandas as pd
import os
parser = argparse.ArgumentParser()
parser.add_argument('--date', type=str)
args = parser.parse_args()

data_dir = '/alps/north/'
out_dir = "/alps/south/"
obs_date = args.date

df = pd.read_csv("../dark_flat_info/dark_flat_info_all.csv")
date = df["date"].values.astype("str")
flat = df["flat"].values
dark = df["dark"].values
cadence = df["cadence"].values.astype("int")
dirs= df["dir"].values.astype("str")

for i in range(len(date)):
    data_dir = dirs[i]
    obs_date_now =date[i]
    gain_exp_for_obj =dark[i]
    gain_exp_for_flat=flat[i]
    
    if  obs_date_now == obs_date:

        ## download index files for astrometry.net
        ## dark & flat
        out_dir_for_dark_flat = os.path.join(out_dir , "dark_flat")
        darkflat.main(data_dir, obs_date, out_dir_for_dark_flat,  gain_exp_for_obj, gain_exp_for_flat)
        
        
        ## WCS using astrometry.net
        out_folder = os.path.join(out_dir , "reduced/%s" % obs_date ) 
        obj_files = seimei_ana.get_target_files(data_dir, obs_date, cadence[i])
        failed_files = seimei_ana.main(data_dir, obs_date, obj_files, out_folder, out_dir_for_dark_flat)
        seimei_ana.save_failed_files(failed_files, out_folder )
