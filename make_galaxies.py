"""
NAME:
	make_galaxies.py

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
import galaxy_classes as galclass
import plotting_things

from astropy.coordinates import Angle
from astropy import units as u

import importlib
importlib.reload(galclass)

def initiate_galaxy_class(galaxy_name, ra_centre, dec_centre, PA, inclination, bmaj, bmin, IR_filename, CO_M0_filename, CO_M1_filename, HI_M0_filename, HI_M1_filename):
    """
    Creates the actual galaxy class and runs through all the functions so that
    everything is set up for plotting stuff.

    Parameters
    ----------
    ra_centre : float
        the RA coordinate for the centre of the galaxy in degrees

    dec_centre : float
        the Declination coordinate for the centre of the galaxy in degrees

    PA : float
        the position angle of the galaxy

    inclination : float
        the inclination angle of the galaxy
    """
    #create the galaxy class
    gal = galclass.Galaxy(galaxy_name, ra_centre, dec_centre, PA, inclination)

    #set the data arrays, headers and wcs
    gal.set_IR_data(IR_filename)
    gal.set_CO_data(CO_M0_filename, CO_M1_filename)
    gal.set_HI_data(HI_M0_filename, HI_M1_filename)

    #add the BMIN and BMAJ to the HI header
    gal.set_HI_BMIN_BMAJ(bmin, bmaj)

    #reproject the data
    gal.reproject_data()

    #convert the HI data to K
    gal.convert_HI_intensity_units()

    #calculate the masses
    gal.calculate_stellar_mass()
    gal.calculate_atomic_gas_mass()
    gal.calculate_molecular_gas_mass()

    return gal







#================================================================
if __name__ == '__main__':
    #make ngc3621 object
    ngc3621 = initiate_galaxy_class(galaxy_name='NGC3621', ra_centre=Angle('11h18m16.3').to(u.deg).value, dec_centre=-32.81192, PA=343.8, inclination=65.8, bmin=2.8437E-03*u.deg, bmaj=4.4300E-03*u.deg, IR_filename='../ismdata/ngc36213621_w1_gauss7p5_interpol.fits', CO_M0_filename='../ismdata/ngc3621_12m+7m+tp_co21_2as_broad_mom0.fits', CO_M1_filename='../ismdata/ngc3621_12m+7m+tp_co21_2as_broad_mom1.fits', HI_M0_filename='../ismdata/NGC_3621_NA_MOM0_THINGS.FITS', HI_M1_filename='../ismdata/NGC_3621_NA_MOM1_THINGS.FITS')

    plotting_things.maps_in_single_figure(ngc3621, titles=['IR data', 'CO Moment 0', 'CO Moment 1', 'HI Moment 0', 'HI Moment 1'], diverging=[False, False, True, False, True], log_data=[True, False, False, False, False])
    plt.savefig(ngc3621.galaxy_name+'_maps.png')


    
