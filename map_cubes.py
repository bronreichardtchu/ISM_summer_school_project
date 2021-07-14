"""
NAME:
	map_cubes.py

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

MODIFICATION HISTORY:
		v.1.0 - first created August 2019

"""
import numpy as np
from astropy.wcs import WCS
from reproject import reproject_interp

import matplotlib.pyplot as plt
import cmasher as cmr

import read_in_cubes as ric

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

    if ax is None:
        ax = plt.subplot(projection=fits_wcs)

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
            print('flattening data a bit')
            data = data[0,0,:,:]
            wcs_list[i] = wcs_list[i][0,0]

        #reproject the data (not for the first one)
        if i > 0:
            reprojected_data, footprint = reproject_interp((data, wcs_list[i]), IR_header)



        #make the figures
        if i < 3:
            map_cube(data, IR_wcs, title=titles[i], ax=axes[int(np.floor(i/3)), i], diverging=diverging[i], log_data=log_data[i])
        else:
            map_cube(data, IR_wcs, title=titles[i], ax=axes[int(np.floor(i/3)), i-3], diverging=diverging[i], log_data=log_data[i])

    plt.show()
