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
import numpy as np

from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy import constants as consts

from reproject import reproject_interp


class Galaxy:
    """
    A class to hold all the information about the galaxy
    """
    def __init__(self, galaxy_name, ra_centre, dec_centre, PA, inclination):
        self.galaxy_name = galaxy_name
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
        self.HI_data = data
        self.HI_header = new_header
        self.HI_vel_data = vel_data/1000 #converting to km/s
        self.HI_vel_header = new_vel_header
        self.HI_wcs = wcs

    def set_HI_BMIN_BMAJ(self, bmin, bmaj):
        """
        Adds these back into the header
        """
        self.HI_header['BMIN'] = bmin
        self.HI_header['BMAJ'] = bmaj
        self.HI_vel_header['BMIN'] = bmin
        self.HI_vel_header['BMAJ'] = bmaj



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


    def reproject_data(self):
        """
        Reprojects all the data to the same as IR
        """
        #reproject the CO data
        reprojected_CO, _ = reproject_interp((self.CO_data, self.CO_header), self.IR_header)
        reprojected_CO_vel, _ = reproject_interp((self.CO_vel_data, self.CO_vel_header), self.IR_header)

        #reproject the HI data
        reprojected_HI, _ = reproject_interp((self.HI_data, self.HI_header), self.IR_header)
        reprojected_HI_vel, _ = reproject_interp((self.HI_vel_data, self.HI_vel_header), self.IR_header)

        #set these for the class
        self.reprojected_CO_data = reprojected_CO
        self.reprojected_CO_vel_data = reprojected_CO_vel
        self.reprojected_HI_data = reprojected_HI
        self.reprojected_HI_vel_data = reprojected_HI_vel


    def set_systematic_velocity(self, sys_vel):
        """
        Sets the systematic velocity
        """
        self.sys_vel = sys_vel


    def correct_coordinates_for_rotation(self, data, wcs):
        """
        Creates RA, Dec and radius arrays, corrected for the rotation of the
        galaxy

        Parameters
        ----------
        data : :obj:'~numpy.ndarray'
            the fits data as a numpy array

        wcs : astropy WCS object
            the world coordinate system for the fits file

        Returns
        -------
        ra : :obj:'~numpy.ndarray'
            the right ascension of each pixel in degrees

        dec : :obj:'~numpy.ndarray'
            the declination of each pixel in degrees

        radius : :obj:'~numpy.ndarray'
            the radius array of the galaxy corrected for rotation
        """
        #make a coordinate grid the same shape as the data
        x, y = np.indices(data.shape)

        #convert the coordinates to wcs
        ra, dec = wcs.all_pix2world(x, y, 0)

        #minus off the centre_coords
        ra = ra - self.ra_centre
        dec = dec - self.dec_centre

        #degrees aren't the same in both directions
        #ra = ra * np.cos(self.dec_centre)

        #PA
        PA = 180-self.PA

        #rotate the ra and dec using a rotation matrix
        ra = ra*np.cos(np.radians(PA)) - dec*np.sin(np.radians(PA))
        dec = ra*np.sin(np.radians(PA)) + dec*np.cos(np.radians(PA))

        #scale the dec for galaxy inclination
        dec = dec / np.cos(np.radians(self.inclination))

        #create the radius array
        radius = np.sqrt(ra**2 + dec**2)

        #find theta
        theta = np.arctan2(dec, ra)

        return ra, dec, radius, theta


    def convert_HI_intensity_units(self):
        """
        Converts the HI intensity from Jy/beam to K
        """
        #first get the beam area
        beam_area = ((self.HI_header['BMAJ']*u.deg)*(self.HI_header['BMIN']*u.deg)/4) * np.pi / np.log(2)

        #get the rest frequency
        rest_freq = self.HI_header['RESTFREQ']*u.Hz

        #get the equivalency object from astropy units
        equiv = u.brightness_temperature(rest_freq)

        #convert the units
        intensity = (self.HI_data*u.Jy/beam_area).to(u.K, equivalencies=equiv).value

        self.HI_data_K = intensity

    def calculate_stellar_mass(self, mass_to_light=0.5):
        """
        Uses the IR data to calculate the stellar mass, corrected for inclination angle
        """
        stellar_mass = 330 * (mass_to_light/0.5) * np.cos(np.radians(self.inclination)) * self.IR_data

        self.stellar_mass = stellar_mass

    def calculate_atomic_gas_mass(self):
        """
        Calculates the atomic gas mass using the HI data (in K)
        """
        #convert the intensity to a number density
        n_H = 1.823e18 * self.HI_data_K * u.cm**-2 *np.cos(np.radians(self.inclination))

        #calculate the mass
        atomic_gas_mass = (n_H * consts.m_p).to(u.Msun/u.pc**2).value

        self.atomic_gas_mass = atomic_gas_mass

    def calculate_molecular_gas_mass(self):
        """
        Calculates the molecular gas mass using the CO data
        """
        molecular_gas_mass = 6.7 * self.CO_data

        self.molecular_gas_mass = molecular_gas_mass
