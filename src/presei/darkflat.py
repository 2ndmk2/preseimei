import glob
import os
import astropy.io.fits as fits
import numpy as np

def sep_files_grz(files):
    """ Sort files by 3 bands
       *0.fits -> g band
        *1.fits -> r band
        *2.fits -> z band
        Params:
            files: Arrays of files

        Returns:
            g_files: Arrays of files in g band
            r_files: Arrays of files in r band
            z_files: Arrays of files in z band
    """    
    g_files = [ file for file in files if "0.fits" in file]
    r_files = [ file for file in files if "1.fits" in file]
    z_files = [ file for file in files if "2.fits" in file]
    return g_files, r_files, z_files


def get_image_med(image_files):
    """ Make median images for arrays of movies

        Params:
            image_files: Arrays of movie files (n * n_t * n_x * n_y)

        Returns:
            image_arr: Arrays of median images  (n * n_x * n_y)
    """

    image_arr = []
    for image_file in image_files:
        hdul = fits.open(image_file)
        data = hdul[0].data
        image_med = np.median(data , axis=0)
        image_arr.append(image_med)
        hdul.close()
        
    return np.array(image_arr)

def get_dark_flat_obj_file(data_dir, obs_date, gain_exp_for_obj ="x8_0.996464", \
    gain_exp_for_flat ="x8_1.992928", dark_name = "DARK", flat_name = "FLAT"):
    """ Get dark & flat files. Note that dark for flat can be different from flat. 

        Params:
            data_dir: Directory for raw data
            obs_date: Observational date (str)
            gain_exp_for_obj: gain & exposure for target and dark for target (XX_XXXXXX), or ($gain_$exposure)
            gain_exp_for_flat: gain & exposure for flat and dark for flat (XX_XXXXXX), or ($gain_$exposure)
            dark_name name for dark in header
            flat_name: name for flat in header
        Returns:
            dark_for_flat_files: Array of files for darks for files
            dark_for_obj_files: Array of files for darks for targets
            flat_files: Array of files for flats
    """

    fns_all = glob.glob(data_dir+obs_date+'/*')
    dark_for_flat_files = []
    dark_for_obj_files = []
    flat_files = []

    for fn in fns_all:
        hdul = fits.open(fn)
        obj_type = hdul[0].header['OBJECT']
        gain = hdul[0].header['GAINCNFG']
        exposure = str(hdul[0].header['EXPTIME1'])
        obs_type = gain+'_'+exposure

        if obj_type == dark_name  and obs_type == gain_exp_for_obj:
            dark_for_obj_files.append(fn)        

        if obj_type == dark_name  and obs_type == gain_exp_for_flat :
            dark_for_flat_files.append(fn)     

        if obj_type == flat_name and obs_type == gain_exp_for_flat :
            flat_files.append(fn)        

        hdul.close()

    return dark_for_flat_files, dark_for_obj_files, flat_files
    
def get_flat_dark(flat_files, dark_for_flat_files, dark_for_obj_files):
    """ Compute median images for darks & flats. Flats are first subtracted by darks for flats

        Params:
            flat_files: Arrays of flat files
            dark_for_flat_files: Arrays of dark files for flats
            dark_for_obj_files: Arrays of dark files for targets

        Returns:
            flat: Flat image
            dark_for_obj_med: Dark image
    """
    flat_before_dark = get_image_med(flat_files)[0]
    dark_for_flat_med_arr = get_image_med(dark_for_flat_files)
    dark_for_flat_med = np.median( dark_for_flat_med_arr, axis=0)
    flat = flat_before_dark  - dark_for_flat_med
    dark_for_obs_arr = get_image_med(dark_for_obj_files)
    dark_for_obj_med = np.median( dark_for_obs_arr, axis=0)
    
    return flat, dark_for_obj_med

def three_bands_flat_dark( dark_for_flat_files, dark_for_obj_files, flat_files):
    """ Compute dark & flats in 3 bands.

        Params:
            dark_for_flat_files: Array of files for darks for files
            dark_for_obj_files: Array of files for darks for targets
            flat_files: Array of files for flats  
        Returns:
            g_flat: Flat image in g band
            r_flat: Flat image in r band
            z_flat: Flat image in z band
            g_dark_for_obj: Dark image in g band for target
            r_dark_for_obj: Dark image in r band for target
            z_dark_for_obj: Dark image in z band for target

    """
    g_dark_for_obj_files, r_dark_for_obj_files, z_dark_for_obj_files = sep_files_grz(dark_for_obj_files)
    g_dark_for_flat_files, r_dark_for_flat_files, z_dark_for_flat_files = sep_files_grz(dark_for_flat_files)
    g_flat_files, r_flat_files, z_flat_files = sep_files_grz(flat_files)
    
    g_flat, g_dark_for_obj = get_flat_dark(g_flat_files, g_dark_for_flat_files, g_dark_for_obj_files)
    r_flat, r_dark_for_obj = get_flat_dark(r_flat_files, r_dark_for_flat_files, r_dark_for_obj_files)
    z_flat, z_dark_for_obj = get_flat_dark(z_flat_files, z_dark_for_flat_files, z_dark_for_obj_files)
    
    return g_flat, r_flat, z_flat, g_dark_for_obj, r_dark_for_obj, z_dark_for_obj

def save_three_bands_flat_dark(g_flat, r_flat, z_flat, g_dark_for_obj, r_dark_for_obj, z_dark_for_obj, \
                               obs_date, out_dir="./"):
    """ Save dark & flat for 3 bands
        Params:
            g_flat: Flat image in g band
            r_flat: Flat image in r band
            z_flat: Flat image in z band
            g_dark_for_obj: Dark image in g band for target
            r_dark_for_obj: Dark image in r band for target
            z_dark_for_obj: Dark image in z band for target
            obs_date: observational date
            out_dir: output directory
        Returns:
            None
    """

    file_name = os.path.join(out_dir, "flat_dark_%s" % obs_date)
    out_dir = np.savez(file_name, g_flat = g_flat, r_flat = r_flat, z_flat = z_flat, g_dark_for_obj = g_dark_for_obj,\
                       r_dark_for_obj = r_dark_for_obj, z_dark_for_obj = z_dark_for_obj)

def main(data_dir, obs_date, out_dir, gain_exp_for_obj ="x8_0.996464", \
    gain_exp_for_flat ="x8_1.992928"):

    dark_for_flat_files, dark_for_obj_files,  flat_files = get_dark_flat_obj_file(data_dir, obs_date, gain_exp_for_obj = gain_exp_for_obj, gain_exp_for_flat = gain_exp_for_flat )
    g_flat, r_flat, z_flat, g_dark_for_obj, r_dark_for_obj, z_dark_for_obj = \
    three_bands_flat_dark(dark_for_flat_files, dark_for_obj_files,  flat_files)
    save_three_bands_flat_dark(g_flat, r_flat, z_flat, g_dark_for_obj, r_dark_for_obj, \
        z_dark_for_obj, obs_date, out_dir=out_dir)

if __name__ == "__main__":

    data_dir = '/alps/south/'
    obs_date = '20220111'
    out_dir = "./"

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    main(data_dir, obs_date, out_dir)

    