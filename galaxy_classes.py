"""
NAME:
	galaxy_classes.py

AUTHOR:
	Bronwyn Reichardt Chu
	Swinburne
	2021

EMAIL:
	<breichardtchu@swin.edu.au>

PURPOSE:
	Create a class to keep all the galaxy stuff in
	Written on MacOS Mojave 10.14.5, with Python 3.7

FUNCTIONS INCLUDED:



"""
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units
from astropy import constants as consts


class Galaxy:
    """
    A class to hold all the information about the galaxy
    """
    def __init__(self, ra_centre, dec_centre, PA, inclination):
        self.ra_centre = ra_centre
        self.dec_centre = dec_centre
        self.PA = PA
        self.inclination = inclination

    def set_IR_data(self, filename):
        """
        Reads in the IR data, creates header and wcs
        """
        #read in the data
        data, header = self.read_in_data_fits(filename)

        #create the wcs
        wcs = self.create_wcs(header)

        #set these for the class
        self.IR_data = data
        self.IR_header = header
        self.IR_wcs = wcs


    def set_CO_data(self, filename_moment0, filename_moment1):
        """
        Reads in the CO data, creates header and wcs
        """
        #read in the data
        data, header = self.read_in_data_fits(filename_moment0)

        vel_data, vel_header = self.read_in_data_fits(filename_moment1)

        #create the wcs
        wcs = self.create_wcs(header)

        #set these for the class
        self.CO_data = data
        self.CO_header = header
        self.CO_vel_data = vel_data
        self.CO_vel_header = vel_header
        self.CO_wcs = wcs


    def set_HI_data(self, filename_moment0, filename_moment1):
        """
        Reads in the HI data, creates header and wcs, and corrects for the extra
        axes
        """
        #read in the data
        data, header = self.read_in_data_fits(filename_moment0)

        vel_data, vel_header = self.read_in_data_fits(filename_moment1)

        #get rid of extra axes
        data = data[0][0]
        vel_data = vel_data[0][0]

        #create the wcs
        wcs = self.create_wcs(header)
        vel_wcs = self.create_wcs(vel_header)

        #fix the wcs
        wcs = wcs.dropaxis(3).dropaxis(2)
        vel_wcs = vel_wcs.dropaxis(3).dropaxis(2)

        #fix the headers
        new_header = wcs.to_header()
        new_header['BUNIT'] = header['BUNIT']
        new_header['RESTFREQ'] = header['RESTFREQ']

        new_vel_header = vel_wcs.to_header()
        new_vel_header['BUNIT'] = "km s-1"

        #set these for the class
        self.CO_data = data
        self.CO_header = new_header
        self.CO_vel_data = vel_data/1000 #converting to km/s
        self.CO_vel_header = new_vel_header
        self.CO_wcs = wcs


    def read_in_data_fits(self, filename):
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


    def create_wcs(self, header):
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
        fits_wcs = WCS(header)

        return fits_wcs
