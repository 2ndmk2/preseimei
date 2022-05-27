import glob
import os
import importlib 
import numpy as np
import astropy.io.fits as fits
import astropy.units as u
import sep
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from preseimei import utils
from astropy.timeseries import LombScargle
import lightkurve as lk


def back_sub_from_movies(movies ):
    bkg_func = lambda x: sep.Background(x, bw = 64, bh = 64, fw = 5, fh = 5)
    movies_bkg_sub = []
    
    for i in range(len(movies)):
        image_now = movies[i]
        image_now2 =  image_now.byteswap().newbyteorder()
        bkg_now= bkg_func(image_now2)
        movies_bkg_sub.append( image_now- bkg_now)
        
    return movies_bkg_sub
        
def choose_stars_for_analysis(x_stars, y_stars, movies_before_bkg_sub, movies_after_bkg_sub, wd_tpf = 30 ):
    wd_tpf_margin = 10 + wd_tpf
    index_for_stars = []
    movie_targets = []
    nt, ny, nx = np.shape(movies_after_bkg_sub)
    for i in range(len(x_stars)):
        
        if wd_tpf_margin<y_stars[i]<ny-wd_tpf_margin  and wd_tpf_margin<x_stars[i]<nx-wd_tpf_margin :
            
            image_now = movies_before_bkg_sub[0][int(y_stars[i])-wd_tpf:int(y_stars[i])+wd_tpf, int(x_stars[i])-wd_tpf:int(x_stars[i])+wd_tpf]
            if len(image_now[image_now==0]):
                continue
                
            else:
                movie_now = movies_after_bkg_sub[:,int(y_stars[i])-wd_tpf:int(y_stars[i])+wd_tpf, int(x_stars[i])-wd_tpf:int(x_stars[i])+wd_tpf]
                index_for_stars.append(i)
                movie_targets.append(movie_now)
    return index_for_stars, movie_targets

def make_aperture(wd = 60, aperture_size = 20):
    wd_all = 2 * wd
    x = np.linspace(0, wd_all  -1,wd_all  )
    y = np.linspace(0, wd_all -1, wd_all )
    x -= np.mean(x)
    y -= np.mean(y)
    xx, yy = np.meshgrid(x, y)
    r =( xx**2 + yy**2)**0.5
    r_ap = r<aperture_size
    return r_ap

def summarize_movie_files(obj_names, out_folder):
    for obj_name in obj_names:
        out_file_name = os.path.join(out_folder, "movie_after_bkgsub_%s" % obj_name)
        files  = utils.get_raw_movie_files_for_gaia_id(obj_name, out_folder)
        time_arr = []
        movie_arr = []        
        for file in files:
            data = np.load(file)
            time_arr.append(data["time"])
            movie_arr.append(data["movie"])
        time_sum = np.concatenate(time_arr)
        movie_sum = np.concatenate(movie_arr)
        np.savez(out_file_name, time = time_sum, movie = movie_sum)
    return None

def get_movie_files_names(obj_names, out_folder):
    file_out = []
    for obj_name in obj_names:
        out_file_name = os.path.join(out_folder, "movie_after_bkgsub_%s.npz" % obj_name)    
        file_out.append(out_file_name)
    return file_out
    
def return_lc_from_tpf(file, r_ap):

    data = np.load(file)
    time_now = data["time"]
    tpf = data["movie"]
    lc = np.sum(tpf * r_ap,axis = (1,2))
    return time_now, lc
        
def make_raw_lcs(obj_names, out_folder, r_ap=10):
    time_all = []
    flux_all = []
    files = get_movie_files_names(obj_names, out_folder)
    for file in files:
        time_now, lc_now = return_lc_from_tpf(file, r_ap)
        arg_now = np.argsort(time_now)
        time_all.append(time_now[arg_now])
        flux_all.append(lc_now[arg_now])
    return time_all, flux_all

def cut_out_images(target_files, stars, out_folder, wd_tpf =30, cadence = 60):
    
    ra_dec = utils.get_ra_dec_from_stars(stars)

    for target_file in target_files:
        print(target_file)
        file_name_header = target_file.split("/")[-1].replace(".fits", "").replace("TRCS00","")
        hdul = fits.open(target_file)
        movie_now = hdul[0].data
        fits_header = hdul[0].header
        nt, ny, nx = np.shape(movie_now)
        
        if nt==cadence:
            time_arr = float(fits_header["MJD-STR"])*24*3600 + (0.5+np.arange(fits_header["NAXIS3"])) * float(fits_header["TFRAME"] )
            w = WCS(fits_header , naxis = 2)
            center = w.pixel_to_world(nx/2, ny/2)
            xy = w.all_world2pix(ra_dec, 0)
            x = xy[:,0]
            y = xy[:,1]
            movie_now = hdul[0].data
            movie_now_bkg_sub = back_sub_from_movies(movie_now)
            movie_now_bkg_sub = np.array(movie_now_bkg_sub)
            index_for_stars, movie_targets = choose_stars_for_analysis(x,y, movie_now, movie_now_bkg_sub, wd_tpf =wd_tpf)    
            
            for (i, index_for_star) in enumerate(index_for_stars):
                gaia_id = stars["source_id"][index_for_stars[i]]
                file_name = "movie_%s_%s" % ( gaia_id, file_name_header)
                file_out_movie =os.path.join(out_folder, file_name)
                np.savez(file_out_movie,time = time_arr,  movie = movie_targets[i])
                
                
def make_target_and_comp_lcs(time_all, flux_all, target_id, gaia_ids):
    time_target = time_all[target_id]
    flux_target = flux_all[target_id]
    flux_for_comp = []
    gaia_for_comp = []
    for i in range(len(time_all)):
        if i == target_id:
            continue
        xy, x_ind, y_ind = np.intersect1d(time_target, time_all[i], return_indices=True)        
        if len(time_target) > len(xy):
            continue
        else:
            flux_for_comp.append(flux_all[i][y_ind])
            gaia_for_comp.append(gaia_ids)
    flux_for_comp = np.array(flux_for_comp)
    gaia_for_comp= np.array(gaia_for_comp)
    return time_target, flux_target, flux_for_comp, gaia_for_comp

def save_lcs_for_all(outdir, time_all, flux_all, i_target, observed_gaia_ids, rp = 10):
    
    outfile_all = os.path.join(outdir, "all_lcs%d" %  rp)
    np.savez(outfile_all, time = time_all, flux= flux_all, i_target =i_target, gaia_ids = observed_gaia_ids)
    
def save_lcs_for_target(outdir, target_id, time_target, flux_target, flux_for_comp, gaia_for_comp,rp = 10):
    
    outfile = os.path.join(outdir, "target_and_others_%s_lc%d" % (target_id,  rp))
    np.savez(outfile, time = time_target, flux_target = flux_target, flux_for_comp = flux_for_comp, gaia_for_comp  = gaia_for_comp)
    
def choose_least_variable(flux_target,flux_for_comp):
    std_arr = []
    flux_med_arr = []    
    for i in range(len(flux_for_comp)):
        flux_med = flux_target/flux_for_comp[i]
        flux_med /=np.median(flux_med)
        std_arr.append(np.std(flux_med))
        flux_med_arr.append(flux_med)
    flux_med_arr = np.array(flux_med_arr)
    least_lc = flux_med_arr[np.argmin(std_arr)]
    return flux_med_arr, std_arr, least_lc

def save_reduced_lc(outdir, target_id, time_target,  least_lc, lc_pca, rp_size):
    outfile = os.path.join(outdir, "target_processed_%s_lc%d" % (target_id,  rp_size))
    np.savez(outfile, time = time_target,  lc_least_var=  least_lc,  lc_pca =  lc_pca)
    
def make_bin_lc(time, flux, bin_dt=5):
    lc = lk.LightCurve(data=None, time=time, flux= flux, flux_err=None)
    lc_bin = lc.bin(bin_dt)
    return lc_bin
    
def show_periodogram_and_lc(time, flux, minimum_period=1, maximum_period=100, bin_dt=5):
    poly = np.polyfit(time,flux, deg=2)
    y_red = np.polyval(poly, time)
    flux_normalized = flux/y_red

    lc = lk.LightCurve(data=None, time=time, flux= flux_normalized, flux_err=None)
    pg = lc.to_periodogram(minimum_period=minimum_period, maximum_period=maximum_period)
    period = pg.period_at_max_power
    pg.plot(view='period') ;
    plt.xscale("log")
    plt.show()
    lc_bin = lc.bin(bin_dt)
    fig=plt.figure(figsize = (14, 7))
    plt.plot(lc.time.value, lc.flux.value, label = "1s")
    plt.plot(lc_bin.time.value, lc_bin.flux.value,label ="%ds" % bin_dt)
    plt.xlabel("Time [s]", fontsize  = 23)
    plt.ylabel("flux", fontsize  = 23)
    plt.legend(fontsize = 15)
    plt.show()
    
    