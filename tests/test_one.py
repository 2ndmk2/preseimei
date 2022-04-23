from preseimei import darkflat
from preseimei import install_indexfiles
from preseimei import seimei_ana

data_dir = '/alps/south/'
obs_date = '20220111'
obj_files = ["/alps/south/20220111/TRCS00110420.fits"]

gain_exp_for_obj = "x8_0.996464"
gain_exp_for_flat = "x8_1.992928"

## download index files for astrometry.net
#install_indexfiles.download_index_files("/usr/local/astrometry/data")

## dark & flat
out_dir_for_dark_flat = "./"
darkflat.main(data_dir, obs_date, out_dir_for_dark_flat,  gain_exp_for_obj, gain_exp_for_flat)

## WCS using astrometry.net
out_folder = "/alps/south/%s_reduced" % obs_date 
seimei_ana.main(data_dir, obs_date, obj_files, out_folder, out_dir_for_dark_flat)
