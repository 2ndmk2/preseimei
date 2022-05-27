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

def x_mean(x):
    ''' Calculate mean of an array
    
    Args:
        x : numpy array of shape (n,num_bamd)
        
    Returns:
        x_mean_arr : numpy array of size n
        
    '''
    n, num_band = np.shape(x)
    x_mean_arr = np.zeros(np.shape(x))
    for i in range(n):
        x_mean_arr[i] = np.mean(x[i])    
    return x_mean_arr

def flux_meddiv_and_subtract_1(lcs):
    nstar, nt = np.shape(lcs)
    flux_med_arr = x_mean(lcs)
    flux_div = np.zeros(np.shape(lcs))
    for i in range(nstar):
        flux_div[i] = (lcs[i]/flux_med_arr[i])-1
    return flux_div


def remove_sys_from_pca(flux_target, flux_for_comp, n_sub_comp=3):

    flux_med_arr = x_mean(flux_for_comp)
    flux_div = flux_meddiv_and_subtract_1(flux_for_comp)
    U, s ,V_svd = np.linalg.svd(flux_div, full_matrices = True)
   
    V_sub_from_lc = V_svd[:n_sub_comp] ## Component for subtraction 
    coeff_V = np.dot(flux_div, V_sub_from_lc.T) ## Coefficient for subtraction 
    flux_div_target = np.array([flux_target])
    flux_med_target = x_mean(flux_div_target )
    coeff_V_target = np.dot(flux_div_target, V_sub_from_lc.T) ## Coefficient for subtraction 
    sub_component_target = np.dot(coeff_V_target, V_sub_from_lc) ## Lightcurve for subtraction 
    after_subtracted_target = (flux_div_target-sub_component_target + 1) * flux_med_target ## Output after subtraction of systematics
    flux_target_pca = after_subtracted_target [0]/np.median(after_subtracted_target [0])
    return  flux_target_pca

def take_name_lc(file):
    file_sp = file.split("/")
    target_name = file_sp[-3]
    band = file_sp[-2]
    lc_name = file_sp[-1].split("_")[-1].replace(".npz", "")
    return target_name, band, lc_name

def make_lc_name(target_name, out_dir, band, lc_name, obs_date):
    out_dir_lc_target = os.path.join( out_dir, target_name)
    if not os.path.exists( out_dir_lc_target ):
        os.makedirs( out_dir_lc_target )
    file_out = os.path.join( out_dir_lc_target , "%s_%s_%s" % (band, obs_date, lc_name))
    return file_out