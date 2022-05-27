import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import glob
import os
import astropy.io.fits as fits
import numpy as np
from preseimei import darkflat
from astropy.wcs import WCS

def get_target_files(data_dir, obs_date, exptime=120):
    """ Get WD files from seimei data
    Params:
        data_dir: Directory for raw data
        obs_date: Observational date (str)

    Returns: 
        obj_files: Arrays for object files
    """
    fns_all = glob.glob(data_dir+obs_date+'/*.fits')
    obj_files = []
    for fn in fns_all:
        try:
            hdul = fits.open(fn)
            obj_type = hdul[0].header['OBJECT']
            n_exposure = hdul[0].header['NAXIS3']
            if obj_type=="DARK" or obj_type =="twilight":
                continue
            print(n_exposure)
            if n_exposure ==exptime:
                obj_files.append(fn)        
        except:
            pass
    return obj_files        

def dark_flat_sub_for_image(fit_name, dark_for_obj, flat_after_dark_sub, out_dir, out_header="_sub"):

    """ Dark & Flat subtraction for obj fits and output it to fit files

    Params:
        fit_name: Filename for fits (str)
        dark_for_obj: Dark image for obj, which is produced from mutiple dark frames
        flat_after_dark_sub: Flat image for obj, which is produced from mutiple flat frames
        out_dir: Directory for output
        out_header: newfile = file_name.replace("fits","") + out_header + ".fits"

    Returns: 
        images: Arrays for reduced images
        file_new_name: Filename for new file
    """

    file_name_last = fit_name.split("/")[-1]
    file_name_new_last = file_name_last.replace(".fits", "%s.fits" % out_header)
    file_name_stacked_last = file_name_last.replace(".fits", "%s_stacked.fits" % out_header)
    
    ## new
    file_new_name = os.path.join(out_dir, file_name_new_last)

    hdul = fits.open(fit_name)
    movie_now = hdul[0].data
    fits_header = hdul[0].header
    images  = (movie_now - dark_for_obj)/( flat_after_dark_sub  )

    ## For the flat values lower than threshold (flat_thr), we insert "0"
    flat_thr = np.median(flat_after_dark_sub)/4 
    images[:,flat_after_dark_sub<flat_thr] = 0    
    hdul[0].data = images

    fits_header = hdul[0].header
    hdul.writeto(file_new_name, overwrite =True)
    hdul.close()
    
    ra_cen = fits_header["RA"]
    dec_cen = fits_header["DEC"]
    
    #stacked
    """
    file_name_stacked  = os.path.join(out_dir, file_name_stacked_last)
    hdul = fits.open(fit_name)
    hdul[0].data = [np.median(images, axis=0)]
    hdul[0].header["NAXIS3"] = 1
    hdul.writeto(file_name_stacked, overwrite =True)
    hdul.close()
    """
    
    return images, file_new_name, ra_cen, dec_cen

def get_dark_flat_from_npz(obj_file, flat_dark):
    """ Get dark & Flat information for fits follwoing conventions of TriCCS filesnames:   
        *0.fits -> g band
        *1.fits -> r band
        *2.fits -> z band

        Params:
            obj_file: Filename for fits (str)
            flat_dark: np object contains flat & dark
        Returns: 
            dark_for_obj: Dark image for obj, which is produced from mutiple dark frames
            flat_after_dark_sub: Flat image for obj, which is produced from mutiple flat frames
    """

    if obj_file.endswith('0.fits'):
        dark_for_obj = flat_dark["g_dark_for_obj"]
        flat_after_dark_sub = flat_dark["g_flat"]
    if obj_file.endswith('1.fits'):
        dark_for_obj = flat_dark["r_dark_for_obj"]
        flat_after_dark_sub = flat_dark["r_flat"]        
    if obj_file.endswith('2.fits'):
        dark_for_obj = flat_dark["z_dark_for_obj"]
        flat_after_dark_sub = flat_dark["z_flat"]    
    return dark_for_obj, flat_after_dark_sub

def make_command_for_astronet(file_name, out_folder, ra_cen, dec_cen, radius = 1):

    """ Make command for solving WCS using https://astrometry.net

        Params:
            file_new_name: Name of fits file for WCS paste
            out_folder: Folder for new fits file with WCS
            ra_cen: Rough position of image center in ra: degrees or hh:mm:ss
            dec_cen: Rough position of image center in dec: degrees or hh:mm:ss
            radius: Radius for search : degrees

        Returns:
            input_command_for_at_net: command for making fits with WCS using astrometry.net, or solve-field


    """
    input_command_for_at_net = "solve-field %s -p -O -D %s --ra %s --dec %s --radius %d" % (file_name, out_folder, ra_cen, dec_cen,  radius)

    return input_command_for_at_net

def garbage_remove_for_astronet(file_name_before_fits, out_folder):

    """ Remove uncessary files from solve-filed. We only use *.new files (fits files)

        Params:
            file_name_before_fits: this part -> /"???".fits  
            out_folder: Folder for new fits file with WCS

        Returns:None
    """
    files = glob.glob("%s/%s*" % (out_folder, file_name_before_fits))
    files_new = glob.glob("%s/%s*.new" % (out_folder, file_name_before_fits))
    files_for_remove =list( set(files) - set(files_new))
    for file_for_remove in files_for_remove:
        os.remove(file_for_remove)

def put_wcs_from_stacked_image_to_target_movies(file_name_stacked, wcs_to_be_put_file):
    file_name_stacked_new = file_name_stacked.replace(".fits", ".new")
    hdul = fits.open(file_name_stacked_new)
    wcs = WCS(hdul[0].header, naxis=2)
    hdul.close()
    
    wcs_to_be_put_file_output = wcs_to_be_put_file.replace(".fits", "_wcs.fits")
    hdul_read = fits.open(wcs_to_be_put_file)
    hdul_read[0].header.update(wcs.to_header())
    hdul_read.writeto(wcs_to_be_put_file_output, overwrite=True)
    hdul.close()
def save_failed_files(failed_files, out_folder ):
    failed_file_name = os.path.join(out_folder , "failed")
    np.savez(failed_file_name, failed_files = failed_files)

def main_notry(data_dir, obs_date,  obj_files, out_folder = "/alps/south/20220111_reduced", dark_flat_folder = "./"):
    
    if not os.path.exists(out_folder ):
        os.makedirs(out_folder )
        
    flat_dark = np.load(os.path.join(dark_flat_folder, "flat_dark_%s.npz" % obs_date))
    failed_files = []

    for (i, obj_file) in enumerate(obj_files):

        ## load data and dark & flat analysis
        dark_for_obj, flat_after_dark_sub = get_dark_flat_from_npz(obj_file, flat_dark)
        images_after_sub, file_new_name, ra_cen, dec_cen = dark_flat_sub_for_image(obj_file, dark_for_obj , flat_after_dark_sub, out_folder)

        ## astrometry.net
        input_command_for_at_net = make_command_for_astronet(file_new_name , out_folder, ra_cen, dec_cen)
        os.system(input_command_for_at_net)

        #put_wcs_from_stacked_image_to_target_movies( file_name_stacked, file_new_name)

        ## remove some files
        file_name_before_fits =obj_file.split("/")[-1].replace(".fits","")
        garbage_remove_for_astronet(file_name_before_fits, out_folder)

    return failed_files

def main(data_dir, obs_date,  obj_files, out_folder = "/alps/south/20220111_reduced", dark_flat_folder = "./"):
    
    if not os.path.exists(out_folder ):
        os.makedirs(out_folder )
        
    flat_dark = np.load(os.path.join(dark_flat_folder, "flat_dark_%s.npz" % obs_date))
    failed_files = []

    for (i, obj_file) in enumerate(obj_files):
        try:
            ## load data and dark & flat analysis
            dark_for_obj, flat_after_dark_sub = get_dark_flat_from_npz(obj_file, flat_dark)
            images_after_sub, file_new_name, ra_cen, dec_cen = dark_flat_sub_for_image(obj_file, dark_for_obj , flat_after_dark_sub, out_folder)
            ## astrometry.net
            input_command_for_at_net = make_command_for_astronet(file_new_name , out_folder, ra_cen, dec_cen)
            os.system(input_command_for_at_net)

            #put_wcs_from_stacked_image_to_target_movies( file_name_stacked, file_new_name)

            ## remove some files
            file_name_before_fits =obj_file.split("/")[-1].replace(".fits","")
            garbage_remove_for_astronet(file_name_before_fits, out_folder)
        except:
            failed_files.append(obj_file)
    return failed_files

if __name__ == "__main__":
    data_dir = '/alps/south/dark_flat"'
    obs_date = '20220111'
    obj_files = get_wd_files(data_dir, obs_date)
    out_folder = "/alps/south/reduced/%s" % obs_date 
    dark_flat_folder = "./"
    
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    main(data_dir, obs_date, obj_files, out_folder, dark_flat_folder)
