import glob
import os
import astropy.io.fits as fits
import numpy as np

def download_index_files(path_to_index_files):
    """ Download index files for astrometry.net http://astrometry.net/ accoriding to 2MASS
        Read http://astrometry.net/doc/readme.html#tricks-and-tips for the detail.
        Seimei FOV is 12 * 6 arcmin^2, so we donwload index-4201~4205. 

        "" (Cited sentences )
        Copy the files to the data directory (path_to_index_files) wherever you installed the \
        Astrometry.net code (INSTALL_DIR/data, perhaps /usr/local/astrometry/data
        ""

        Params:
            path_to_index_files: path to data_directory
        Returns:
            None
    """

    ## in total 30 GB
    #os.system("curl -o %s -O http://broiler.astrometry.net/~dstn/4200/index-4200-[00-47].fits")
    os.system("curl -o %s -O http://broiler.astrometry.net/~dstn/4200/index-4201-[00-47].fits" % path_to_index_files)
    os.system("curl -o %s -O http://broiler.astrometry.net/~dstn/4200/index-4202-[00-47].fits" % path_to_index_files)
    os.system("curl -o %s -O http://broiler.astrometry.net/~dstn/4200/index-4203-[00-47].fits" % path_to_index_files)
    os.system("curl -o %s -O http://broiler.astrometry.net/~dstn/4200/index-4204-[00-47].fits" % path_to_index_files)
    os.system("curl -o %s -O http://broiler.astrometry.net/~dstn/4200/index-4205-[00-47].fits" % path_to_index_files)


    
def download_gaia_index_files(path_to_index_files):
    """ Download index files for astrometry.net http://astrometry.net/ accoriding to GAIA DR2
        Read http://astrometry.net/doc/readme.html#tricks-and-tips for the detail.
        Seimei FOV is 12 * 6 arcmin^2, so we donwload index-4201~4205. 

        "" (Cited sentences )
        Copy the files to the data directory (path_to_index_files) wherever you installed the \
        Astrometry.net code (INSTALL_DIR/data, perhaps /usr/local/astrometry/data
        ""

        Params:
            path_to_index_files: path to data_directory
        Returns:
            None
    """

    ## in total 30 GB
    os.system("curl -O https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-5201-[00-47].fits")
    os.system("curl -O https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-5202-[00-47].fits")
    os.system("curl -O https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-5203-[00-47].fits")
    os.system("curl -O https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-5204-[00-47].fits")
    os.system("curl -O https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-5205-[00-47].fits")


    