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

import importlib



def map_cube(data, fits_wcs, title, cbar_label, ax=None, log_data=False, diverging=False, **formats):
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
        cbar_label = 'log ' + cbar_label

    #if there's no axis provided, create one.
    if ax is None:
        ax = plt.subplot(projection=fits_wcs)

    #if radius_contours == True:


    ax.set_facecolor('white')
    #do the plotting
    if diverging == True:
        map = ax.imshow(data, origin='lower', cmap=cmr.pride, **formats)
    else:
        map = ax.imshow(data, origin='lower', cmap=cmr.ember, **formats)
    cbar = plt.colorbar(map, label=cbar_label, ax=ax, fraction=0.05, pad=0.04, extend='max')

    #ax.invert_xaxis()

    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')

    ax.set_title(title)

    #plt.show()

    return ax



def maps_in_single_figure(galaxy, titles, diverging, log_data):
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
    #create a list of the data arrays
    data_arrays = [galaxy.IR_data, galaxy.reprojected_CO_data, galaxy.reprojected_CO_vel_data, galaxy.reprojected_HI_data, galaxy.reprojected_HI_vel_data]

    #list of wcs
    wcs_list = [galaxy.IR_wcs, galaxy.CO_wcs, galaxy.CO_wcs, galaxy.HI_wcs, galaxy.HI_wcs]

    #formating
    formats = [{}, dict(vmax=40, vmin=-10), dict(vmax=900, vmin=550), {}, {}]

    cbar_labels = [galaxy.IR_header['BUNIT'], galaxy.CO_header['BUNIT'], galaxy.CO_vel_header['BUNIT'], galaxy.HI_header['BUNIT'], galaxy.HI_vel_header['BUNIT']]


    #create a figure
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(10,6), subplot_kw={'projection': galaxy.IR_wcs})

    #iterate through the data
    for i, data in enumerate(data_arrays):
        #correct the coordinates for rotation
        ra, dec, radius, theta = galaxy.correct_coordinates_for_rotation(data, galaxy.IR_wcs)

        #make the figures
        if i < 3:
            map_cube(data, galaxy.IR_wcs, title=titles[i], cbar_label=cbar_labels[i], ax=axes[int(np.floor(i/3)), i], diverging=diverging[i], log_data=log_data[i], **formats[i])

            axes[int(np.floor(i/3)), i].contour(np.fliplr(radius), colors='gray', levels=[0.01, 0.06, 0.1, 0.2, 0.5])

        else:
            map_cube(data, galaxy.IR_wcs, title=titles[i], cbar_label=cbar_labels[i], ax=axes[int(np.floor(i/3)), i-3], diverging=diverging[i], log_data=log_data[i], **formats[i])

            axes[int(np.floor(i/3)), i-3].contour(np.fliplr(radius), colors='gray', levels=[0.01, 0.06, 0.1, 0.2, 0.5])

    axes[1,2].axis('off')

    plt.subplots_adjust(left=0.08, right=0.94, top=0.99, bottom=0.04, wspace=0.69, hspace=0.03)

    plt.show()


def quick_look_profile(galaxy):
    """
    Plots the radial profile of the IR intensity

    Parameters
    ----------
    galaxy : class
        the galaxy class
    """
    #correct the coordinates
    ra, dec, radius, theta = galaxy.correct_coordinates_for_rotation(galaxy.IR_data, galaxy.IR_wcs)

    #make the plot
    plt.figure()
    plt.scatter(radius, galaxy.IR_data)
    plt.ylabel('I '+galaxy.IR_header['BUNIT'])
    plt.xlabel('Radius')

    plt.show()



def intensity_profile(galaxy, labels, z):
    """
    Plots the intensity vs the radius
    """
    #calculate the proper distance
    proper_dist = cosmo.kpc_proper_per_arcmin(z).to(u.kpc/u.degree)

    #create a list of the data arrays
    data_arrays = [galaxy.IR_data, galaxy.reprojected_CO_data, galaxy.reprojected_HI_data]

    #list of wcs
    wcs_list = [galaxy.IR_wcs, galaxy.CO_wcs, galaxy.HI_wcs]


    #create a figure
    plt.figure()

    #iterate through the data
    for i, data in enumerate(data_arrays):
        #correct the coordinates for rotation
        ra, dec, radius, theta = galaxy.correct_coordinates_for_rotation(data, galaxy.IR_wcs)

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


def mass_profile(galaxy, labels, z):
    """
    Plots the mass vs the radius
    """
    #calculate the proper distance
    proper_dist = cosmo.kpc_proper_per_arcmin(z).to(u.kpc/u.degree)

    #create a list of the data arrays
    data_arrays = [galaxy.stellar_mass, galaxy.atomic_gas_mass, galaxy.molecular_gas_mass]

    #list of wcs
    wcs_list = [galaxy.IR_wcs, galaxy.HI_wcs, galaxy.CO_wcs]


    #create a figure
    plt.figure()

    #iterate through the data
    for i, data in enumerate(data_arrays):
        #correct the coordinates for rotation
        ra, dec, radius, theta = galaxy.correct_coordinates_for_rotation(data, galaxy.IR_wcs)

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
    plt.ylabel('Mass')

    plt.yscale('log')

    plt.show()
