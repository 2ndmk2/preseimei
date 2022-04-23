import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import glob
import os
import astropy.io.fits as fits
import numpy as np
from preseimei import darkflat



def get_wd_files(data_dir, obs_date, exptime=120):
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
        hdul = fits.open(fn)
        obj_type = hdul[0].header['OBJECT']
        gain = hdul[0].header['GAINCNFG']
        exposure = str(hdul[0].header['EXPTIME1'])
        n_exposure = hdul[0].header['NAXIS3']
        if obj_type.startswith('WD_') and n_exposure ==exptime:
            obj_files.append(fn)        
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

    return images, file_new_name

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

def make_command_for_astronet(file_name, out_folder):

    """ Make command for solving WCS using https://astrometry.net

        Params:
            file_new_name: Name of fits file for WCS paste
            out_folder: Folder for new fits file with WCS

        Returns:
            input_command_for_at_net: command for making fits with WCS using astrometry.net, or solve-field


    """
    input_command_for_at_net = "solve-field %s -p -O -D %s" % (file_name, out_folder)

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

def main(data_dir, obs_date,  obj_files, out_folder = "/alps/south/20220111_reduced", dark_flat_folder = "./"):

    flat_dark = np.load(os.path.join(dark_flat_folder, "flat_dark_%s.npz" % obs_date))

    for (i, obj_file) in enumerate(obj_files):

        ## load data and dark & flat analysis
        dark_for_obj, flat_after_dark_sub = get_dark_flat_from_npz(obj_file, flat_dark)
        images_after_sub, file_new_name = dark_flat_sub_for_image(obj_file, dark_for_obj , flat_after_dark_sub, out_folder)
        
        ## astrometry.net
        input_command_for_at_net = make_command_for_astronet(file_new_name, out_folder)
        os.system(input_command_for_at_net)

        ## remove some files
        file_name_before_fits =obj_file.split("/")[-1].replace(".fits","")
        garbage_remove_for_astronet(file_name_last, out_folder)

if __name__ == "__main__":
    data_dir = '/alps/south/'
    obs_date = '20220111'
    obj_files = get_wd_files(data_dir, obs_date)
    out_folder = "/alps/south/%s_reduced" % obs_date 
    dark_flat_folder = "./"
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    main(data_dir, obs_date, obj_files, out_folder, dark_flat_folder)
