from preseimei import darkflat
from preseimei import install_indexfiles
from preseimei import seimei_ana
import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--date', type=str)
args = parser.parse_args()


data_dir = '/alps/south/'
obs_date = args.date
obj_files = ["/alps/center/20220524/TRCS00145662.fits"]

gain_exp_for_obj = "x8_0.996464"
gain_exp_for_flat = "x8_0.996464"

## download index files for astrometry.net
#install_indexfiles.download_index_files("/usr/local/astrometry/data")

## dark & flat
out_dir_for_dark_flat = "/alps/south/dark_flat"
darkflat.main(data_dir, obs_date, out_dir_for_dark_flat,  gain_exp_for_obj, gain_exp_for_flat)

## WCS using astrometry.net
out_folder = "/alps/south/reduced/%s" % obs_date 

failed_files = seimei_ana.main_notry(data_dir, obs_date, obj_files, out_folder, out_dir_for_dark_flat)
seimei_ana.save_failed_files(failed_files, out_folder )