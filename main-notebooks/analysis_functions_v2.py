#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import astropy.units as units
from astropy import table
from astropy.coordinates import SkyCoord
from matplotlib.patches import Polygon, Ellipse
from shapely.geometry.polygon import Polygon as shapelyPolygon
from shapely.geometry import Point as shapelyPoint


#-------------------------------- General functions ---------------------------------#
# (generalized and may be useful even outside of my analysis)

def angular_radius_cut(table, center, radius):
    """Return a new table with radius cut applied. The radius cut considers angular separation.
    """
    ra0, dec0 = center
    
    alpha1 = np.radians(ra0)
    delta1 = np.radians(dec0)
    alpha2 = np.radians(table['ra'])
    delta2 = np.radians(table['dec'])
    
    # Calculating the angular separation
    angular_sep_rad = np.arccos(np.sin(delta1) * np.sin(delta2) + \
                            np.cos(delta1) * np.cos(delta2) * np.cos(alpha1 - alpha2))

    return table[np.degrees(angular_sep_rad) < radius]

def distance_to_dm(dist):
    """dist unit is pc.
    """
    return 5 * np.log10(dist / 10)

def dm_to_distance(dm):
    """Return distance in units of pc.
    """
    return 10**(dm / 5) * 10

def poly_2deg(x, a, b, c):
    """A quadratic equation for fitting some data. 
    """
    y = a * x**2 + b * x + c
    
    return y

def pm_cut_with_gradient(table_data, pmra_curve_params, pmdec_curve_params, pad=0):
    """Return an astropy table with selected stars from <table_data>, taking into
    account proper motion gradients along both RA and Dec directions. Each gradient is
    described by a polynomial curve, where <pmra_curve_params> and <pmdec_curve_params>
    contain the coefficients for the polynomial curves.
    
    The selection procedure is the following:
    - Take a star in <table_data>, use its pmra_error, pmdec_error, and pmra_pmdec_corr
    to create a covariance matrix. Then use the matrix to find semi-major and semi-minor 
    axes of the covariance matrix ellipse, and use shapelyPolygon to create the ellipse.
    - Take this star's RA coordinate, combined with coefficients given in <pmra_curve_params>
    and <pmdec_curve_params>, and input these into the quadratic function poly_2deg, from
    which we will obtain a pair of (pmra, pmdec) coordinates. This is the predicted location
    of this star in the proper motion space.
    - If the ellipse contains this predicted point, add it to the table of selected stars.
    Otherwise, move to the next star in <table_data>.
    
    Parameters
    ----------
    table_data: 
        An astropy table containing information about stars
    pad: 
        A value to increase the ellipse size. Added in quadrature with the PM errors on
        each star.
    pmra_curve_params: 
        polynomial coefficients for the function poly_2deg, 
        for the model curve in PMRA vs. RA
    pmdec_curve_params:
        polynomial coefficients for the function poly_2deg, 
        for the model curve in PMDEC vs. RA
    
    Precondition
    ------------
    table_data should include the following column names: 
    pmra, pmdec, pmra_error, pmdec_error, pmra_pmdec_corr
    
    Return
    ------
    A new astropy table with selected members.
    """    
    selected_stars = table.Table()
    selected_stars.columns = table_data.columns
    selected_stars.remove_rows(slice(0, len(selected_stars), 1))
    
    # For parametrizing the ellipse later
    t = np.arange(0, 2*np.pi, 0.01)

    for i in range(len(table_data)):
        #--- CREATE ELLIPSE ---#
        
        # adding the pad to pm error
        pmra_sigma_i = np.sqrt(table_data['pmra_error'][i]**2 + pad**2) 
        pmdec_sigma_i = np.sqrt(table_data['pmdec_error'][i]**2 + pad**2)
        
        # covariance matrix entries
        var_x = pmra_sigma_i**2
        var_y = pmdec_sigma_i**2
        cov_xy = table_data['pmra_pmdec_corr'][i] * pmra_sigma_i * pmdec_sigma_i

        # Creating the matrix for better understanding
        cov_matrix = np.array([[var_x, cov_xy], 
                               [cov_xy, var_y]])
        a = cov_matrix[0, 0]
        c = cov_matrix[1, 1]
        b = cov_matrix[0, 1]

        # The sqrt of the lambdas are the semi-major and semi-minor axes
        lambda1 = (a+c)/2 + np.sqrt(((a-c)/2)**2 + b**2)
        lambda2 = (a+c)/2 - np.sqrt(((a-c)/2)**2 + b**2)

        if b == 0 and a >= c:
            theta = 0
        elif b == 0 and a < c:
            theta = np.pi/2
        else:
            theta = np.arctan2(lambda1 - a, b)

        x_t = np.zeros(len(t) + 1)
        y_t = np.zeros(len(t) + 1)
        x_t[:-1] = lambda1**0.5 * np.cos(theta) * np.cos(t) - lambda2**0.5 * np.sin(theta) * np.sin(t) + table_data['pmra'][i]
        y_t[:-1] = lambda1**0.5 * np.sin(theta) * np.cos(t) + lambda2**0.5 * np.cos(theta) * np.sin(t) + table_data['pmdec'][i]
        # Make a closed shape to create a polygon
        x_t[-1] = x_t[0]
        y_t[-1] = y_t[0]

        ellipse = shapelyPolygon(np.array([x_t, y_t]).T)
        
        #--- EXPECTED PM CENTER FROM RA ---#
        pmra_expected = poly_2deg(table_data['ra'][i], *pmra_curve_params)
        pmdec_expected = poly_2deg(table_data['ra'][i], *pmdec_curve_params)

        pm_center = shapelyPoint(pmra_expected, pmdec_expected)
        
        #--- CHECKING FOR SELECTION ---#
        if ellipse.contains(pm_center):
            selected_stars.add_row(table_data[i])
    
    return selected_stars


#------------------------ Specific functions for my analysis ------------------------#
# (probably not able to carry over to other works)

def deredden_mag_S5(table):
    """Return the g, r, i, z dereddened magnitudes for a table from S5 data.
    
    Note: other data tables may not work because of the keywords.
    """
    g_band = table['decam_g'] - 3.185 * table['ebv']
    r_band = table['decam_r'] - 2.140 * table['ebv']
    i_band = table['decam_i'] - 1.569 * table['ebv']
    z_band = table['decam_z'] - 1.196 * table['ebv']
    
    return g_band, r_band, i_band, z_band

def crossmatch(table1, table2, sep_dist=1, get_indices=False, get_sep2d=False):
    """table1 is catalogue sky, table2 is stars to be matched to the catalogue sky.
    So the match indices are for table2.
    
    Note: sep_dist is in arcseconds.
    """
    c = SkyCoord(ra=table1['ra']*units.degree, dec=table1['dec']*units.degree)
    stars = SkyCoord(ra=table2['ra']*units.degree, dec=table2['dec']*units.degree)

    match_index, sep2d, _ = c.match_to_catalog_sky(stars)
    
    matched_stars = table2[match_index][sep2d < sep_dist*units.arcsec]

    if get_indices and not get_sep2d:
        return matched_stars, match_index
    elif not get_indices and get_sep2d:
        return matched_stars, sep2d
    elif get_indices and get_sep2d:
        return matched_stars, match_index, sep2d
    else:
        return matched_stars
    
def get_BHB_iso_color_n_mag():
    """This is an empirical BHB isochrone taken from Prof. Li's example notebook.
    Note that the distance modulus to M92 is subtracted so that this isochrone is
    shifted back to its original y-position. All we need to provide is the distance
    modulus to our target.
    """
    # Empirical BHB isochrone
    m92ebv = 0.023
    m92ag = m92ebv * 3.184
    m92ar = m92ebv * 2.130
    m92_hb_r = np.array([16.8, 15.8, 15.38, 15.1, 15.05])
    m92_hb_col = np.array([-0.36, -0.3, -0.2, -0.0, 0.1])
    m92_hb_g = m92_hb_r + m92_hb_col
    des_m92_hb_g = m92_hb_g - 0.104 * (m92_hb_g - m92_hb_r) + 0.01
    des_m92_hb_r = m92_hb_r - 0.102 * (m92_hb_g - m92_hb_r) + 0.02
    des_m92_hb_g = des_m92_hb_g - m92ag
    des_m92_hb_r = des_m92_hb_r - m92ar

    dm_m92_harris = 14.59
    
    # Here are the g and r being plotted
    g_bhb = des_m92_hb_g
    r_bhb = des_m92_hb_r
    
    # The dm to M92 is subtracted here, we have to add our BooIII dm 
    color = g_bhb - r_bhb
    mag = g_bhb - dm_m92_harris
    
    return color, mag

def apply_CMD_filter(data_table, g, r, RGB_color, RGB_g, BHB_color, BHB_g, delta_color, delta_g, dimmest_mag, get_polygon=False):
    """There are some messy things to keep track of.
    
    RGB_iso takes in the g and r isochrone arrays, and they do NOT include the dm shift
    BHB_iso takes in g and color (g - r) and this g DOES include the dm shift in the input
    It is this way because BHB is added in later and it is empirical, a bit hard to merge properly
    with the other code
    """
    g_mask = RGB_g < dimmest_mag
    g_mask2 = BHB_g < dimmest_mag
    
    RGB_color = RGB_color[g_mask]
    RGB_g = RGB_g[g_mask]
    BHB_color = BHB_color[g_mask2]
    BHB_g = BHB_g[g_mask2]
    
    # --- RGB ---
    left = RGB_color - delta_color
    right = RGB_color + delta_color
    up = RGB_g + delta_g
    down = RGB_g - delta_g
    
    # Create polygons for RGB
    x_RGB = np.append(left, np.flip(right))
    y_RGB = np.append(RGB_g, np.flip(RGB_g))
    
    RGB_df = pd.DataFrame()
    RGB_df['x_RGB'] = x_RGB
    RGB_df['y_RGB'] = y_RGB
    polygon_RGB = Polygon(RGB_df, alpha=0.2)
    
    x_RGB2 = np.append(RGB_color, np.flip(RGB_color))
    y_RGB2 = np.append(up, np.flip(down))
    
    RGB_df2 = pd.DataFrame()
    RGB_df2['x_RGB'] = x_RGB2
    RGB_df2['y_RGB'] = y_RGB2
    polygon_RGB2 = Polygon(RGB_df2, alpha=0.2)  

    # --- BHB ---
    left2 = BHB_color - delta_color
    right2 = BHB_color + delta_color
    up2 = BHB_g + delta_g
    down2 = BHB_g - delta_g
    
    # Create polygons for RGB
    x_BHB = np.append(left2, np.flip(right2))
    y_BHB = np.append(BHB_g, np.flip(BHB_g))
    
    BHB_df = pd.DataFrame()
    BHB_df['x_BHB'] = x_BHB
    BHB_df['y_BHB'] = y_BHB
    polygon_BHB = Polygon(BHB_df, alpha=0.2)
    
    x_BHB2 = np.append(BHB_color, np.flip(BHB_color))
    y_BHB2 = np.append(up2, np.flip(down2))
    
    BHB_df2 = pd.DataFrame()
    BHB_df2['x_BHB'] = x_BHB2
    BHB_df2['y_BHB'] = y_BHB2
    polygon_BHB2 = Polygon(BHB_df2, alpha=0.2)
    
    # Check for contained stars
    stars = pd.DataFrame()
    stars['color'] = g - r
    stars['mag'] = g

    contained = polygon_RGB.contains_points(stars) | polygon_RGB2.contains_points(stars) | \
                polygon_BHB.contains_points(stars) | polygon_BHB2.contains_points(stars)

    new_table = data_table[contained]

    # Return the polygon bounds if asked
    if get_polygon:
        return new_table, polygon_RGB, polygon_RGB2, polygon_BHB, polygon_BHB2
    else:
        return new_table

def cmd_cut_BHB(table0, g_column_name, r_column_name, BHB_gr, BHB_g, delta_gr, delta_g, get_polygon=False):
    selected = table.Table()
    
    #--- CREATE POLYGONS ---#
    left = BHB_gr - delta_gr
    right = BHB_gr + delta_gr
    up = BHB_g + delta_g
    down = BHB_g - delta_g
    
    gr1 = np.append(left, np.flip(right))
    g1 = np.append(BHB_g, np.flip(BHB_g))
    gr1 = np.append(gr1, gr1[0])
    g1 = np.append(g1, g1[0])
    polygon1 = shapelyPolygon(np.array([gr1, g1]).T)
    
    gr2 = np.append(BHB_gr, np.flip(BHB_gr))
    g2 = np.append(up, np.flip(down))
    gr2 = np.append(gr2, gr2[0])
    g2 = np.append(g2, g2[0])
    polygon2 = shapelyPolygon(np.array([gr2, g2]).T)

    for i in range(len(table0)):
        star = shapelyPoint(table0[g_column_name][i] - table0[r_column_name][i], table0[g_column_name][i])
        
        if polygon1.contains(star) or polygon2.contains(star):
            selected = table.vstack([selected, table0[i]])

    if get_polygon:
        return selected, polygon1, polygon2
    else:
        return selected
    
def append_BHB_dm(table, g_mag_key, r_mag_key):
    """Calculate distance modulus for a table of BHB stars and add the results to the table
    (creates new columns for the table instead of returning a new table)
    """
    gr = table[g_mag_key] - table[r_mag_key]
    
    # Mg from Belokurov & Koposov (2016)
    table['Mg'] = 0.398 - 0.392*gr + 2.729*gr**2 + 29.1128*gr**3  + 113.569*gr**4 

    # get dm from g and Mg
    table['dm'] = table[g_mag_key] - table['Mg']

def color_color_cut(data_table, stellar_locus, y_shift_top, y_shift_bottom=0, get_polygon=False):
    """Selects metal-poor stars using a stellar locus
    
    data_table : must contain headers from Legacy Survey DECaLS DR9
    stellar_locus : given file, we use this to set the bounds for the color-color selection
    """
    gr = stellar_locus['g'] - stellar_locus['r']
    rz = stellar_locus['r'] - stellar_locus['z']
    
    # Create polygon
    x = np.append(gr, np.flip(gr))
    y = np.append(rz + y_shift_top, np.flip(rz + y_shift_bottom))
    
    df = pd.DataFrame()
    df['x'] = x
    df['y'] = y
    polygon = Polygon(df, alpha=0.2)
    
    # Check for contained stars
    stars = pd.DataFrame()
    stars['g-r'] = np.array(data_table['dered_mag_g'] - data_table['dered_mag_r'])
    stars['r-z'] = np.array(data_table['dered_mag_r'] - data_table['dered_mag_z'])

    contained = polygon.contains_points(stars)

    new_table = data_table[contained]

    # Return the polygon bounds if asked
    if get_polygon:
        return new_table, polygon
    else:
        return new_table

def plot_rh(ax):
    """A convenient function to plot the half-light radius for BooIII (1, 3, and 5 times).
    
    Values taken from https://iopscience.iop.org/article/10.3847/1538-4357/ab7459
    """
    RA_BOO3, DEC_BOO3 = 209.3, 26.8
    
    # In Moskowitz & Walker, the presented radius is "circularized".
    # Its relation with the elliptical radius (semi-major axis) is:
    # R_circ = R_elliptical * sqrt(1 - e), in units of arcmin
    rh_elliptical = 33.03 / (1 - 0.33)**0.5  
    
    rh = rh_elliptical / 60  # Half-light radius [deg]
    PA = 278.91              # Position Angle [deg]
    ell = 0.33               # Ellipticity

    ells = Ellipse(xy=(RA_BOO3, DEC_BOO3), 
                   width=2*rh, 
                   height=2*rh*(1-ell),
                   angle=90-PA, fc=None, ec='k', ls='--', lw=2, fill=False, label='$r_h$, 3$r_h$, 5$r_h$') 
    ells2 = Ellipse(xy=(RA_BOO3, DEC_BOO3), 
                   width=6*rh, 
                   height=6*rh*(1-ell),
                   angle=90-PA, fc=None, ec='k', ls='--', lw=2, fill=False) 
    ells3 = Ellipse(xy=(RA_BOO3, DEC_BOO3), 
                   width=10*rh, 
                   height=10*rh*(1-ell),
                   angle=90-PA, fc=None, ec='k', ls='--', lw=2, fill=False)

    ax.add_artist(ells)
    ax.add_artist(ells2)
    ax.add_artist(ells3)
    
    
def add_delta_coordinates(table, center):
    """
    Modify the input <table> by adding the columns "ra_delta" and "dec_delta" (in degrees).
    These are the original coordinates shifted to the axes where (0, 0) is the
    center of the target galaxy.
    
    <table> should contain the (RA, Dec) coordinates in degrees.
    <center> is the (RA, Dec) of the target galaxy, also in degrees.
    """
    ra0 = np.radians(center[0])
    dec0 = np.radians(center[1])
    ra = np.radians(table['ra'])
    dec = np.radians(table['dec'])
    
    ra_delta_rad = np.cos(dec) * np.sin(ra - ra0)
    dec_delta_rad = np.sin(dec) * np.cos(dec0) - np.cos(dec) * np.sin(dec0) * np.cos(ra - ra0)
    
    table['ra_delta'] = np.degrees(ra_delta_rad)
    table['dec_delta'] = np.degrees(dec_delta_rad)
    
    
def radius_cut_with_delta_coords(table, radius):
    """Input table must have the "delta" celestial coordinates 
    obtained from add_delta_coordinates().
    """
    r = np.sqrt(table['ra_delta']**2 + table['dec_delta']**2)
    return table[r < radius]


def half_light_radius_cut_along_semimajor_axis(table, position_angle, ellipticity, radius):
    """Note radius may be multiples of the half-light radius.
    """
    q = 1 - ellipticity
    angle = np.radians(90 - position_angle)
    
    # Rotate to align with semi-major axis
    ra_delta_rotated = table['ra_delta'] * np.cos(angle) + table['dec_delta'] * np.sin(angle)
    dec_delta_rotated = - table['ra_delta'] * np.sin(angle) + table['dec_delta'] * np.cos(angle)
    
    # Elliptical radius of all stars
    r_elliptical = np.sqrt(ra_delta_rotated**2 + dec_delta_rotated**2 / q**2)
    
    return table[r_elliptical < radius]


def plot_rh_projected(ax):
    """This is an alternative version of plot_rh(), which plots the radius with respect to
    delta RA and delta Dec, with the centre at (RA, Dec) = (0, 0).
    
    A convenient function to plot the half-light radius for BooIII (1, 3, and 5 times).
    Values taken from https://iopscience.iop.org/article/10.3847/1538-4357/ab7459
    """    
    # In Moskowitz & Walker, the presented radius is "circularized".
    # Its relation with the elliptical radius (semi-major axis) is:
    # R_circ = R_elliptical * sqrt(1 - e), in units of arcmin
    rh_elliptical = 33.03 / (1 - 0.33)**0.5  
    
    rh = rh_elliptical / 60  # Half-light radius [deg]
    PA = 278.91              # Position Angle [deg]
    ell = 0.33               # Ellipticity

    ells = Ellipse(xy=(0, 0), 
                   width=2*rh, 
                   height=2*rh*(1-ell),
                   angle=90-PA, fc=None, ec='k', ls='--', lw=2, fill=False, label='$r_h$, 3$r_h$, 5$r_h$') 
    ells2 = Ellipse(xy=(0, 0), 
                   width=6*rh, 
                   height=6*rh*(1-ell),
                   angle=90-PA, fc=None, ec='k', ls='--', lw=2, fill=False) 
    ells3 = Ellipse(xy=(0, 0), 
                   width=10*rh, 
                   height=10*rh*(1-ell),
                   angle=90-PA, fc=None, ec='k', ls='--', lw=2, fill=False)

    ax.add_artist(ells)
    ax.add_artist(ells2)
    ax.add_artist(ells3)
    
    
# ------------------------------------------------------- #
# ----------This is for computing RRL distance----------- #
# ------------------------------------------------------- #

# extinction correction
import getDust # Compute the Gaia extinctions assuming relations from Babusieux (Gaia Collaboration et al. (2018))

def add_dm_to_RRLs(table, feh_mean):
    """
    This code is provided by Prof. Ting Li
    (ting.li@astro.utoronto.ca)
    
    Modify the input table (which must contain Gaia RR Lyrae table keywords)
    and add the distance modulus for those RR Lyrae stars.
    
    [Fe/H] values for those stars must be provided. The parameter feh_mean is
    the mean [Fe/H] of the sample.
    """
    # input: Gaia G, BP, RP, output: extinction in each band
    table['AG'], table['Abp'], table['Arp'] = getDust.getDust(table['phot_g_mean_mag'],
                                                              table['phot_bp_mean_mag'],
                                                              table['phot_rp_mean_mag'],
                                                              0.86*table['EBV_SFD'])

    # distance calibration from Muraveva et al. (2018)
    # (Assume the given RRLs have metallicity from our MCMC fit of the base sample)
    table['MG'] = 0.32 * feh_mean + 1.11

    # get dm from G, MG and extinction
    table['dm'] = table['phot_g_mean_mag'] - table['AG'] - table['MG']
# ------------------------------------------------------- #
    
    
   