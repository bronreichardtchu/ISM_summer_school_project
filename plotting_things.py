"""
NAME:
	plotting_things.py

AUTHOR:
	Bronwyn Reichardt Chu
	Swinburne
	2021

EMAIL:
	<breichardtchu@swin.edu.au>

PURPOSE:
	To map the data
	Written on MacOS Mojave 10.14.5, with Python 3.7

FUNCTIONS INCLUDED:
    map_cube
    maps_in_single_figure

"""
import numpy as np
from scipy import stats

from astropy.wcs import WCS
from astropy import units as u
from astropy import constants
from astropy.cosmology import WMAP9 as cosmo

from reproject import reproject_interp

import matplotlib.pyplot as plt
import cmasher as cmr

import read_in_cubes as ric
import getting_coordinates as get_coords

import importlib
importlib.reload(get_coords)

def map_cube(data, fits_wcs, title, ax=None, log_data=False, diverging=False):
    """
    Maps the fits data to an imshow on a matplotlib axis

    Parameters
    ----------
    data : :obj:'~numpy.ndarray'
        the fits data as a numpy array

    fits_wcs : astropy WCS object
        the world coordinate system for the fits file (if not inputing an axes)

    title : str
        the title for the plot

    ax : matplotlib Axes instance
        the axes you want to plot onto.  Default is None, creates its own axes.
        If you already have axes, they need to already be in the right projection.

    log_data : boolean
        whether to take the logarithm of the data (default=False)

    diverging : boolean
        whether a diverging colour map is needed (default=False)

    Returns
    -------
    map of the data
    """
    if log_data == True:
        data = np.log10(data)

    #if there's no axis provided, create one.
    if ax is None:
        ax = plt.subplot(projection=fits_wcs)

    #if radius_contours == True:


    ax.set_facecolor('white')
    #do the plotting
    if diverging == True:
        map = ax.imshow(data, origin='lower', cmap=cmr.pride)
    else:
        map = ax.imshow(data, origin='lower', cmap=cmr.ember)#, vmin=-1.75, vmax=-0.25)
    cbar = plt.colorbar(map, ax=ax)

    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')

    ax.set_title(title)

    #plt.show()

    return ax


def fix_header(data, header, wcs):
    """
    Fixes the header, particularly for HI data
    """
    print('flattening data a bit')
    data = data[0,0,:,:]
    wcs = wcs.dropaxis(2).dropaxis(2)
    new_header = wcs.to_header()
    new_header['BUNIT'] = header['BUNIT']
    new_header['RESTFREQ'] = header['RESTFREQ']

    return data, new_header, wcs


def maps_in_single_figure(IR_file, CO_M0_file, CO_M1_file, HI_M0_file, HI_M1_file, titles, diverging, log_data):
    """
    Reads in data from the filenames and maps them

    Parameters
    ----------
    filenames : list of str
        a list of the filenames of data to be included

    Returns
    -------
    All the maps in a single figure
    """
    #read in all the files
    IR_data, IR_header = ric.read_in_data_fits(IR_file)
    IR_wcs = ric.create_wcs(IR_header)

    CO_M0_data, CO_M0_header = ric.read_in_data_fits(CO_M0_file)
    CO_M0_wcs = ric.create_wcs(CO_M0_header)

    CO_M1_data, CO_M1_header = ric.read_in_data_fits(CO_M1_file)
    CO_M1_wcs = ric.create_wcs(CO_M1_header)

    HI_M0_data, HI_M0_header = ric.read_in_data_fits(HI_M0_file)
    HI_M0_wcs = ric.create_wcs(HI_M0_header)

    HI_M1_data, HI_M1_header = ric.read_in_data_fits(HI_M1_file)
    HI_M1_wcs = ric.create_wcs(HI_M1_header)

    #create a list of the data arrays
    data_arrays = [IR_data, CO_M0_data, CO_M1_data, HI_M0_data, HI_M1_data]

    #list of headers
    headers = [IR_header, CO_M0_header, CO_M1_header, HI_M0_header, HI_M1_header]

    #list of wcs
    wcs_list = [IR_wcs, CO_M0_wcs, CO_M1_wcs, HI_M0_wcs, HI_M1_wcs]

    #reproject the data to match the stellar data
    #CO_M0_data, CO_M0_footprint = reproject_interp((CO_M0_data, CO_M0_header), IR_wcs)
    #CO_M1_data, CO_M1_footprint = reproject_interp((CO_M1_data, CO_M1_header), IR_wcs)
    #HI_M0_data, HI_M0_footprint = reproject_interp((HI_M0_data, HI_M0_header), IR_wcs)
    #HI_M1_data, HI_M1_footprint = reproject_interp((HI_M1_data, HI_M1_header), IR_wcs)

    #create a figure
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(10,6), subplot_kw={'projection': IR_wcs})

    #iterate through the data
    for i, data in enumerate(data_arrays):

        #if the shape has too many parts get rid of the outer ones
        #this is particularly for the THINGS data
        if len(data.shape) > 2:
            data, headers[i], wcs_list[i] = fix_header(data, headers[i], wcs_list[i])


        #reproject the data (not for the first one)
        if i > 0:
            reprojected_data, footprint = reproject_interp((data, wcs_list[i]), IR_header)



        #make the figures
        if i < 3:
            map_cube(data, IR_wcs, title=titles[i], ax=axes[int(np.floor(i/3)), i], diverging=diverging[i], log_data=log_data[i])
        else:
            map_cube(data, IR_wcs, title=titles[i], ax=axes[int(np.floor(i/3)), i-3], diverging=diverging[i], log_data=log_data[i])

    plt.show()



def convert_intensity_mass_HI(data, CO_header):
    """
    Convert the intensity units for the HI data.
    Should go from Jy/beam m/s to K m/s
    """
    #get the major and minor axis FWHM values from the fits header
    FWHM_major = CO_header['BMAJ']
    FWHM_minor = CO_header['BMIN']

    #do the converstion to temp
    intensity = (6.07 * 10e5 * data) / (FWHM_major * FWHM_minor)

    #convert to number density
    n_H = 1.823e18 * intensity

    n_H = n_H * (u.cm*u.cm)**(-1)

    #then convert to solar masses per kpc^2
    sigma_mass = n_H * constants.m_p

    sigma_mass = sigma_mass.to(u.solMass/u.kpc**2)

    return sigma_mass



def convert_intensity_mass_IR(data, mass_to_light=0.5):
    """
    Converts from intensity to mass for the IR WISE data, using a mass-to-light
    ratio of 0.5
    """
    sigma_mass = 330*(mass_to_light/0.5) * data

    return sigma_mass



def convert_intensity_mass_CO(data, alpha_CO=4.35):
    """
    Converts from intensity to mass for the CO data
    """
    sigma_mol = 6.7 * (alpha_CO ** (1-0/4.35)) * data

    return sigma_mol



def intensity_profile(IR_file, CO_M0_file, HI_M0_file, labels, centre_coords, PA, inclination, z):
    """
    Plots the intensity vs the radius
    """
    #calculate the proper distance
    proper_dist = cosmo.kpc_proper_per_arcmin(z).to(u.kpc/u.degree)

    #read in IR data
    IR_data, IR_header = ric.read_in_data_fits(IR_file)

    #create the figure
    plt.figure()

    #iterate through the files
    for i, data_file in enumerate([IR_file, CO_M0_file, HI_M0_file]):
        #read in data
        data, header = ric.read_in_data_fits(data_file)

        #get the wcs
        fits_wcs = ric.create_wcs(header)

        if len(data.shape) > 2:
            data, header, fits_wcs = fix_header(data, header, fits_wcs)

        #reproject the data (not for the first one)
        #if i > 0:
        #    data, _ = reproject_interp((data, fits_wcs), IR_header)

        #get the coordinates
        ra, dec = get_coords.make_coords_array(data, fits_wcs, centre_coords)

        #correct for inclination and PA and get radius
        ra, dec, radius, theta = get_coords.correct_for_rotation(ra, dec, PA, inclination)

        #multiply through the radius to make it in kpc
        radius = radius * proper_dist

        #bin by radius
        binned_median_intensity, bin_edges, binnumber = stats.binned_statistic(radius.reshape(-1), data.reshape(-1), statistic='mean', bins=200)
        bin_width = (bin_edges[1] - bin_edges[0])
        bin_centers = bin_edges[1:] - bin_width/2

        #plot the figure
        plt.scatter(bin_centers, binned_median_intensity, label=labels[i])

    plt.legend()
    plt.xlabel('Radius [kpc]')
    plt.ylabel('Intensity')

    #plt.yscale('log')

    plt.show()


def mass_profile(IR_file, CO_M0_file, HI_M0_file, labels, centre_coords, PA, inclination, z):
    """
    Plots the mass vs the radius
    """
    #calculate the proper distance
    proper_dist = cosmo.kpc_proper_per_arcmin(z).to(u.kpc/u.degree)

    #create the figure
    plt.figure()

    #iterate through the files
    for i, data_file in enumerate([IR_file, CO_M0_file, HI_M0_file]):
        #read in data
        data, header = ric.read_in_data_fits(data_file)

        #get the wcs
        fits_wcs = ric.create_wcs(header)

        if len(data.shape) > 2:
            data, header, fits_wcs = fix_header(data, header, fits_wcs)

        #get the coordinates
        ra, dec = get_coords.make_coords_array(data, fits_wcs, centre_coords)

        #correct for inclination and PA and get radius
        ra, dec, radius, theta = get_coords.correct_for_rotation(ra, dec, PA, inclination)

        #multiply through the radius to make it in kpc
        radius = radius * proper_dist

        #convert the data to masses
        if i == 0:
            data = convert_intensity_mass_IR(data)

        if i == 1:
            data = convert_intensity_mass_CO(data)

        if i == 2:
            data = convert_intensity_mass_HI(data, header)

        #bin by radius
        binned_median_intensity, bin_edges, binnumber = stats.binned_statistic(radius.reshape(-1), data.reshape(-1), statistic='mean', bins=200)
        bin_width = (bin_edges[1] - bin_edges[0])
        bin_centers = bin_edges[1:] - bin_width/2

        #plot the figure
        plt.scatter(bin_centers, binned_median_intensity, label=labels[i])

    plt.legend()
    plt.xlabel('Radius [kpc]')
    plt.ylabel('Mass')

    plt.yscale('log')

    plt.show()
