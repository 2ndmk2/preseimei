import sep
import matplotlib.pyplot as plt
from preseimei import utils
from preseimei import movie_lc_ana
from preseimei import pca
import importlib
import argparse
import pandas as pd 

parser = argparse.ArgumentParser()
parser.add_argument('--date', type=str)
parser.add_argument('--band', type=str)
parser.add_argument('--dir', type=str)
args = parser.parse_args()

data_dir = '/alps/south/reduced/'
obs_date = args.date
band = args.band
wd_tpf =30

df = pd.read_csv("../dark_flat_info/dark_flat_info_all.csv")
date = df["date"].values.astype("str")
cadence = df["cadence"].values.astype("int")
dirs= df["dir"].values.astype("str")
cadence_now = cadence[date==obs_date]

obj_names = utils.get_object_names(data_dir, obs_date)
seimei_gaia_ids = utils.read_seimei_targetfile("./seimei_targets.csv")
for target in obj_names:
    out_folder  = utils.make_outdir(data_dir, band, obs_date, target)
    target_files = utils.get_target_files(data_dir, obs_date, target ,  band)
    stars, ra, dec, gaia_mag, ra_dec  = utils.get_gaia_stars_wcs_information(target_files, obs_date, target, band)
    movie_lc_ana.cut_out_images(target_files, stars, out_folder, wd_tpf, cadence_now)
    importlib.reload(movie_lc_ana)
    
    rp_sizes = [10, 15, 20, 25]
    for rp_size in rp_sizes:
        rp_mask= movie_lc_ana.make_aperture(30, rp_size)
        observed_gaia_ids = utils.get_gaiaids_from_outputfile(out_folder)
        movie_lc_ana.summarize_movie_files(observed_gaia_ids, out_folder)
        time_all, flux_all = movie_lc_ana.make_raw_lcs(observed_gaia_ids,  out_folder, rp_mask)
        
        i_target, gaia_id_target= utils.get_target_id(observed_gaia_ids, seimei_gaia_ids )
        print(i_target, gaia_id_target)
        time_target, flux_target, flux_for_comp, gaia_for_comp = movie_lc_ana.make_target_and_comp_lcs(time_all, flux_all , i_target, observed_gaia_ids)
        movie_lc_ana.save_lcs_for_all(out_folder,  time_all, flux_all, i_target, observed_gaia_ids, rp_size)
        movie_lc_ana.save_lcs_for_target(out_folder, gaia_id_target, time_target, flux_target, flux_for_comp, gaia_for_comp, rp_size)
        flux_med_arr, std_arr, least_lc = movie_lc_ana.choose_least_variable(flux_target, flux_for_comp)
        flux_target_pca = pca.remove_sys_from_pca(flux_target, flux_for_comp)
        movie_lc_ana.save_reduced_lc(out_folder, gaia_id_target, time_target,  least_lc, flux_target_pca, rp_size)
        
