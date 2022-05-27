from astroquery.mast import Catalogs
from astropy.wcs import WCS
import astropy.io.fits as fits
import glob
import os
import astropy.units as u
import numpy as np
import pandas as pd

def get_gaiaids_from_outputfile(out_folder):
    
    files = glob.glob(os.path.join(out_folder , "*_sub.new.npz"))
    observed_gaia_ids = set()
    for file in files:
        plate_name = file.split("_")[-3]
        observed_gaia_ids.add(int(plate_name))
    observed_gaia_ids = list(observed_gaia_ids)
   
    return observed_gaia_ids

def get_target_id(observed_gaia_ids, seimei_target_gaia_ids):
    
    for (i, observed_gaia_id) in enumerate(observed_gaia_ids):
        if observed_gaia_id in seimei_target_gaia_ids:
            return i, observed_gaia_id 
    return None, None
    
def read_seimei_targetfile(file_name="./seimei_targets.csv"):
    gaia_ids = pd.read_csv("./seimei_targets.csv")[" GAIA ID"].values
    return gaia_ids
    
def make_outdir(data_dir, band, obs_date, target):

    if band =="0":
        out_folder = data_dir  + "%s_analyzed/movies/%s/g_band/" % (obs_date, target)
    if band=="1":
        out_folder = data_dir  + "%s_analyzed/movies/%s/r_band/" % (obs_date, target)
    if band=="2":
        out_folder = data_dir  + "%s_analyzed/movies/%s/z_band/" % (obs_date, target)        
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)         
    return out_folder 
    
def get_gaia_stars(target, radius):

    """ Get nearby stars from Gaia 
        Args:
            target: Target name
            radius: Search radius (give it with unit of "u.arcsec" in Astropy)
        Returns:
            ra: Arrays of ra
            dec: Arrays of dec
            gaia_mag: Arrays of Gaia magnitude
            stars: dataframe for stars
    """  
    while(1):
        try:
            stars = Catalogs.query_object(target, radius=radius, catalog='Gaia')    
            if len(stars)>0:
                break
        except:
            pass
    ra, dec = stars['ra'], stars['dec']    
    gaia_mag = stars['phot_g_mean_mag']
    return ra, dec, gaia_mag, stars

def get_ra_dec_from_stars(stars):
    ra = stars["ra"]
    dec = stars["dec"]
    ra_dec = np.array([ra, dec]).T
    return ra_dec

def get_gaia_stars_wcs_information(target_files, obs_date, target, band_info = "0", lim_mag = 17.5):
    hdul = fits.open(target_files[0])
    fits_header = hdul[0].header
    nx = fits_header["NAXIS1"]
    ny = fits_header["NAXIS2"]
    w = WCS(fits_header , naxis = 2)
    center = w.pixel_to_world(nx/2, ny/2)
    target_cen = "%f %f" % (center.icrs.ra.deg, center.icrs.dec.deg, )
    ra, dec, gaia_mag, stars = get_gaia_stars(target_cen, 10 * 60 * u.arcsec)
    mask =gaia_mag < lim_mag 
    ra = ra[mask]
    dec = dec[mask]
    gaia_mag = gaia_mag[mask]
    stars = stars[mask]
    ra_dec = np.array([ra, dec]).T
    hdul.close()    
    return stars, ra, dec, gaia_mag, ra_dec  

def get_target_files(data_dir, obs_date, target,  band="0"):
    fns_all = glob.glob(data_dir+obs_date+'/*%s_sub.new' % band)
    fns_target= []
    for fn in fns_all:
        hdul = fits.open(fn)
        obj_type = hdul[0].header['OBJECT']
        if obj_type  == target:
            fns_target.append(fn)
    return list(set(fns_target))

def get_object_names(data_dir, obs_date):
    fns_all = glob.glob(os.path.join(os.path.join(data_dir,obs_date),'*.new'))
    obj_names = []
    for fn in fns_all:
        hdul = fits.open(fn)
        obj_type = hdul[0].header['OBJECT']
        obj_names.append(obj_type)
    return list(set(obj_names))

def get_raw_movie_files_for_gaia_id(gaia_id, out_folder):
    files = glob.glob(os.path.join(out_folder, "*%s*_sub.new.npz" % gaia_id))
    return files

