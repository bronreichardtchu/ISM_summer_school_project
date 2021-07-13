"""
NAME:
	read_in_cubes.py

AUTHOR:
	Bronwyn Reichardt Chu
	Swinburne
	2021

EMAIL:
	<breichardtchu@swin.edu.au>

PURPOSE:
	To read in the data
	Written on MacOS Mojave 10.14.5, with Python 3.7

FUNCTIONS INCLUDED:
    read_in_data_fits
    create_wcs

MODIFICATION HISTORY:
		v.1.0 - first created August 2019

"""


from astropy.io import fits
from astropy.wcs import WCS
from astropy import units
from astropy import constants as consts



def read_in_data_fits(filename):
    """
    Reads in the data if it is contained in a fits file

    Parameters
    ----------
    filename : str
        points to the file

    Returns
    -------
    data : :obj:'~numpy.ndarray'
        the fits data as a numpy array

    header : FITS header
        the fits header
    """
    #open file and get data and header
    with fits.open(filename) as hdu:
        data = hdu[0].data
        header = hdu[0].header
    hdu.close()

    return data, header


def create_wcs(header):
    """
    Reads in the fits header and creates the wcs

    Parameters
    ----------
    header : FITS header object
        the header for the fits file to read in

    Returns
    -------
    fits_wcs : astropy WCS object
        the world coordinate system for the fits file
    """
    #create the WCS
    fits_wcs = WCS(fits_header)

    return fits_wcs
