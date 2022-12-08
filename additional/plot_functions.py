#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt

params = {'figure.titlesize': 16, #20,
          'axes.labelsize': 16,
          'axes.titlesize': 16, #20,
          'axes.linewidth': 1.5,
          'xtick.minor.visible': True,
          'xtick.major.width': 1.5,
          'xtick.major.size': 6,
          'xtick.minor.size': 4,
          'xtick.labelsize': 12,
          'xtick.direction': 'in',
          'xtick.top': True,
          'ytick.minor.visible': True,
          'ytick.major.width': 1.5,
          'ytick.major.size': 6,
          'ytick.minor.size': 4,
          'ytick.labelsize': 12,
          'ytick.direction': 'in',
          'ytick.right': True,
          'legend.fontsize': 14, #12,
          'legend.labelspacing': 0.3,
          'legend.borderaxespad': 0.8,
          'font.family': 'sans-serif'}

def plot_labels(xlabel, ylabel, title):
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

def ax_labels(ax, xlabel, ylabel, title):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

def plot_onsky(fig, ax, table, label, cbar_array=[], cbar_label='', cm_name='Reds', **kwargs):
    """Note: only works if the table headers match the ones used below.
    """
    cm = plt.cm.get_cmap(cm_name)
    
    if len(cbar_array) != 0:
        sc = ax.scatter(table['ra'], table['dec'], c=cbar_array, cmap=cm, label=label, **kwargs)
        cbar = fig.colorbar(sc, ax=ax)
        cbar.ax.set_ylabel(cbar_label)
    else:
        ax.scatter(table['ra'], table['dec'], label=label, **kwargs)

def plot_pm(table, label, **kwargs):
    plt.errorbar(table['pmra'], table['pmdec'], 
                 xerr=table['pmra_error'], yerr=table['pmdec_error'], label=label, **kwargs)

def plot_pm_ax(ax, table, label, **kwargs):
    ax.errorbar(table['pmra'], table['pmdec'], 
                 xerr=table['pmra_error'], yerr=table['pmdec_error'], label=label, **kwargs)
    
def plot_RGB_isochrone(ax, g, r, dm, label, **kwargs):
    ax.plot(g-r, g+dm, label=label, **kwargs)

def plot_BHB_isochrone(ax, dm, label='', **kwargs):
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
    
    # The dm to M92 is already subtracted here
    ax.plot(g_bhb-r_bhb, g_bhb-dm_m92_harris+dm, label=label, **kwargs)

def plot_cmd(fig, ax, g, r, label, cbar_array=[], cbar_label='', cm_name='Reds', **kwargs):
    """This function does NOT take a table, it takes the magnitudes directly.
    """
    cm = plt.cm.get_cmap(cm_name)
    
    if len(cbar_array) != 0:
        sc = ax.scatter(g-r, g, c=cbar_array, cmap=cm, label=label, **kwargs)
        cbar = fig.colorbar(sc, ax=ax)
        cbar.ax.set_ylabel(cbar_label)
    else:
        ax.scatter(g-r, g, label=label, **kwargs)    

def plot_rv_metallicity(table, label, **kwargs):
    plt.errorbar(table['feh50'], table['vel_calib'], 
                 xerr=table['feh_calib_std'], yerr=table['vel_calib_std'], label=label, **kwargs)

def plot_rv_metallicity_ax(ax, table, label, **kwargs):
    ax.errorbar(table['feh50'], table['vel_calib'], 
                 xerr=table['feh_calib_std'], yerr=table['vel_calib_std'], label=label, **kwargs)

def plot_pmra_ax(ax, table, label, **kwargs):
    ax.errorbar(table['ra'], table['pmra'], yerr=table['pmra_error'], label=label, **kwargs)

def plot_pmdec_ax(ax, table, label, **kwargs):
    ax.errorbar(table['ra'], table['pmdec'], yerr=table['pmdec_error'], label=label, **kwargs)
    
def plot_rv_s5_ax(ax, table, label, **kwargs):
    ax.errorbar(table['ra'], table['vel_calib'], yerr=table['vel_calib_std'], label=label, **kwargs)    
    
# def plot_spatial_rapid(ax, table):
#     colorlist = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'y']

#     for i in range(len(table)):
#         ax.plot(table['ra'][i], table['dec'][i], marker='o', c=colorlist[i])

# def plot_errorbar_rapid(ax, table, xname, yname, xerrname, yerrname):
#     colorlist = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'y']
    
#     for i in range(len(table)):
#         if xerrname == 'noxerrname':
#             ax.errorbar(table[xname][i], table[yname][i], 
#                         yerr=table[yerrname][i], fmt='o', c=colorlist[i], capsize=0, lw=1)
#         else:
#             ax.errorbar(table[xname][i], table[yname][i], 
#                         xerr=table[xerrname][i], yerr=table[yerrname][i], fmt='o', c=colorlist[i], capsize=0, lw=1)
        
# def plot_cmd_rapid(fig, ax, table, xname, yname):
#     colorlist = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'y']

#     for i in range(len(table)):
#         plot_cmd(fig, ax, table[xname][i], table[yname][i], '', marker='o', c=colorlist[i])

# NOTE: these ones won't work if there are more than 9 stars


def plot_spatial_rapid(ax, table, colors=True, color='none'):
    if colors:
        colorlist = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'y']

        for i in range(len(table)):
            ax.plot(table['ra'][i], table['dec'][i], marker='o', c=colorlist[i])
    else:
        for i in range(len(table)):
            ax.plot(table['ra'][i], table['dec'][i], marker='o', c=color)

def plot_errorbar_rapid(ax, table, xname, yname, xerrname, yerrname, colors=True, color='none'):
    if colors:
        colorlist = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'y']

        for i in range(len(table)):
            if xerrname == 'noxerrname':
                ax.errorbar(table[xname][i], table[yname][i], 
                            yerr=table[yerrname][i], fmt='o', c=colorlist[i], capsize=0, lw=1)
            else:
                ax.errorbar(table[xname][i], table[yname][i], 
                            xerr=table[xerrname][i], yerr=table[yerrname][i], fmt='o', c=colorlist[i], capsize=0, lw=1)
    else:
        for i in range(len(table)):
            if xerrname == 'noxerrname':
                ax.errorbar(table[xname][i], table[yname][i], 
                            yerr=table[yerrname][i], fmt='o', c=color, capsize=0, lw=1)
            else:
                ax.errorbar(table[xname][i], table[yname][i], 
                            xerr=table[xerrname][i], yerr=table[yerrname][i], fmt='o', c=color, capsize=0, lw=1)
                
def plot_cmd_rapid(fig, ax, table, xname, yname, colors=True, color='none'):
    if colors:
        colorlist = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'y']

        for i in range(len(table)):
            plot_cmd(fig, ax, table[xname][i], table[yname][i], '', marker='o', c=colorlist[i])
    else:
        for i in range(len(table)):
            plot_cmd(fig, ax, table[xname][i], table[yname][i], '', marker='o', c=color)
        