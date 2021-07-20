"""
NAME:
	getting_coordinates.py

AUTHOR:
	Bronwyn Reichardt Chu
	Swinburne
	2021

EMAIL:
	<breichardtchu@swin.edu.au>

PURPOSE:
	To get the coordinates from the data
	Written on MacOS Mojave 10.14.5, with Python 3.7

FUNCTIONS INCLUDED:
    find_the_centre
    make_coords_array
    correct_for_rotation


"""
import numpy as np
from astropy.wcs import WCS

def find_the_centre(data, data_wcs):
    """
    Finds the centre of the galaxy

    Parameters
    ----------
    data : :obj:'~numpy.ndarray'
        the fits data as a numpy array

    data_wcs : astropy WCS object
        the world coordinate system for the fits file

    Returns
    -------

    """
    #think of a way to find the centre of the galaxy


def make_coords_array(data, data_wcs, centre_coords):
    """
    Makes a coordinate array such that it's offset by the galaxy centre

    Parameters
    ----------
    data : :obj:'~numpy.ndarray'
        the fits data as a numpy array

    data_wcs : astropy WCS object
        the world coordinate system for the fits file

    centre_coords : list of floats
        the RA and Dec coordinates of the centre of the galaxy

    Returns
    -------
    ra : :obj:'~numpy.ndarray'
        the right ascension of each pixel in degrees

    dec : :obj:'~numpy.ndarray'
        the declination of each pixel in degrees

    """
    #if the shape has too many parts get rid of the outer ones
    #this is particularly for the THINGS data
    if len(data.shape) > 2:
        print('flattening data a bit')
        data = data[0,0,:,:]
        data_wcs = data_wcs.dropaxis(2).dropaxis(2)


    #make a coordinate grid the same shape as the data
    x, y = np.indices(data.shape)

    #convert the coordinates to wcs
    ra, dec = data_wcs.all_pix2world(x, y, 0)

    #minus off the centre_coords
    ra = ra - centre_coords[0]
    dec = dec - centre_coords[1]

    return ra, dec

def correct_for_rotation(ra, dec, PA, inclination):
    """
    Corrects for the position angle and inclination of the galaxy

    Parameters
    ----------
    ra : :obj:'~numpy.ndarray'
        the right ascension of each pixel without rotation

    dec : :obj:'~numpy.ndarray'
        the declination of each pixel without rotation

    PA : float
        the position angle of the galaxy

    inclination : float
        the inclination of the galaxy

    Returns
    -------
    ra : :obj:'~numpy.ndarray'
        the right ascension of each pixel with rotation

    dec : :obj:'~numpy.ndarray'
        the declination of each pixel with rotation

    radius : :obj:'~numpy.ndarray'
        the radius array of the galaxy corrected for rotation
    """
    #rotate the ra and dec using a rotation matrix
    ra = ra*np.cos(np.radians(PA)) + dec*np.sin(np.radians(PA))
    dec = -ra*np.sin(np.radians(PA)) + dec*np.cos(np.radians(PA))

    #scale the dec for galaxy inclination
    dec = dec / np.cos(np.radians(inclination))

    #create the radius array
    radius = np.sqrt(ra**2 + dec**2)

    #find theta ======> NOTE TO SELF: WHAT IS THETA?
    theta = np.arctan2(dec, ra)

    return ra, dec, radius, theta
