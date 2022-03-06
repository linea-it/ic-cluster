# Funções para o notebook de Object Selection

import os
import numpy as np 
import healpy as hp
from astropy.table import Table
from collections import OrderedDict
#from gavodb import DBManager
import sqlalchemy
import pylab as pl
import seaborn as sns
import pandas as pd 
from scipy.stats import gaussian_kde, linregress
import matplotlib.pyplot as plt


def get_vac(pid='6607',
            vac_schema='vac_cluster', 
            bands=['g', 'r', 'i', 'z'],
            sample_frac = 0.1):
    """
    Get data from VAC a table in database

    Parameters
    ----------
    pid: str
        Process id
    columnst: list
        List of columns to extract from pz table, besides magnitudes
    bands: list
        List of bands for magnitudes
    zmin: float
        Minimum value for z
    zmax: float
        Maximum value for z
    magmin: float
        Minimum value for mag
    magmax: float
        Maximum value for mag
    sample_frac: float
        Fraction for sampling large tables

    Returns
    -------
    Pandas DataFrame object 
    """
    
    # Connect to database    
    engine = sqlalchemy.create_engine('postgres://untrustedprod:untrusted@desdb4.linea.gov.br:5432/prod_gavo')
    conn = engine.connect()
     
    table_name = f"vac_{vac_schema}.catalog_{pid}"    
    columns=['coadd_objects_id', 'ra', 'dec', 'hpix_4096', 'z_best', 'err_z']
    for band in bands:
        columns.append(f"mag_{band}")
        columns.append(f"magerr_{band}")
    query = f"select {(',').join(columns)} from {table_name} where random() < {sample_frac}"
    stm = sqlalchemy.sql.text(query)
    result = conn.execute(stm).fetchall()  
    
    vac_dataframe = pd.DataFrame(data=result, columns=columns)
    vac_dataframe.set_index('coadd_objects_id', inplace=True)
    
    return vac_dataframe


def find_green_valley(x, y, mag_bins=np.arange(14,28,0.5)):
    green_valley = np.zeros_like(mag_bins)
    
    

def redseq_fit(x, y, z, color_cut=0.3, mag_bins=np.arange(14,28,0.5), istar_dic=None):
    counts, xbins = np.histogram(x, bins=mag_bins)
    #print(counts)
    #print(len(counts), len(mag_bins))
    #mags = np.array(mag_bins[:-1])[(counts[0]>0)]
    dmag = (mag_bins[1]-mag_bins[0])
    mags = np.arange(mag_bins[0]+(dmag/2.), mag_bins[-1], dmag)
    #print(mags)
    #print(len(mags))
    colors = np.zeros_like(mags)#_bins)
    #print(len(counts), len(mags), len(colors))
    
    #colors = np.array([    for i, mag in enumerate(mag_bins[:-1])])
    for i, mag in enumerate(mag_bins[:-1]): 
        colors_in_mag_bin = y[(x>=mag_bins[i])&(x<mag_bins[i+1])&(y>color_cut)]
        counts_in_mag_bin, ybins = np.histogram(colors_in_mag_bin)
        mode = ybins[np.argmax(counts_in_mag_bin)]
        colors[i] = mode
         
     
    #colors = np.array([np.median(y[(x>=mag_bins[i])&(x<mag_bins[i+1])&(y>color_cut)]) for i, mag in enumerate(mag_bins[:-1])])
    #mags = np.array(mag_bins[:-1])[(counts[0]>10)]
    #colors = colors[counts[0]>10]
    
    #print(len(mags), len(colors))
    
    fit_result = linregress(mags[(counts>50)&(mags<istar_dic[z]+1.5)], colors[(counts>50)&(mags<istar_dic[z]+1.5)])
    a = fit_result[0]
    b = fit_result[1]
    #print(mags)
    #print(colors)
    xfit = np.arange(10, 30, 1)
    yfit = a*xfit + b
    
    return mags, colors, a, b, xfit, yfit, counts 
    
    
def cmd_plot(x, y, bins=[200,200], plot_range=[[16,23.2],[-1, 3]], weights=None, cmin=1, cmax=None, 
             z_range=(0.0,0.1), title='', x_label='', y_label='', panel=1, istar_dic=False, 
             color_cut=0.3, dmag=1.5):
    
    p = plt.subplot(1,3,panel)
    plt.title('%s (%d gals)'%(title, len(x)), fontsize=14) 
    #plt.hexbin(x, y, C=C, gridsize=250, mincnt=1, cmap='rainbow')
    plt.hist2d(x, y, bins=bins, range=plot_range, weights=weights, cmin=cmin, cmax=cmax, cmap='rainbow')
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=12)
    if (panel==3 or panel==6):
        cbar.set_label('density', size=12)
    plt.xlim(plot_range[0])
    plt.ylim(plot_range[1])
    plt.xlabel(x_label, fontsize=16)
    plt.ylabel(y_label, fontsize=16)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    tx = (min(plot_range[0])+(max(plot_range[0])-min(plot_range[0]))*0.05)
    ty = (min(plot_range[1])+(max(plot_range[1])-min(plot_range[1]))*0.9)
    ttext='%.2f < z < %.2f'%(z_range[0], z_range[1])
    plt.text(tx, ty, ttext, fontsize=14)

    z = round((max(z_range)+min(z_range))/2., 2)
    plt.vlines(istar_dic[z], -5, 5, linestyles='dotted')    
    plt.vlines(istar_dic[z]+dmag, -5, 5, linestyles='dotted')    
    
    # Red Sequence 
    mask_red = (y>=color_cut)&(x<(istar_dic[z]+dmag))
    # first fit
    mags1, colors1, a1, b1, xfit1, yfit1, counts = redseq_fit(x[mask_red], y[mask_red], z, color_cut=color_cut, istar_dic=istar_dic)
    
    b2 = color_cut - a1*22.
    #print(b2)
    y2 = a1*xfit1 + b2
    plt.plot(xfit1, y2, 'k--', lw=0.8)
           
    #plt.hlines(color_cut, -50, 50, linestyles='dashed', lw=0.8)
    
    
    
    mask_red2 = (y>=((a1*x)+b2))#&(x<(istar_dic[z]+dmag))
    #PAREI AQUI
    
    #print(mask_red2)
    
    # second fit (after selecting red sequence cut with a slope)
    
    mags, colors, a, b, xfit, yfit, counts = redseq_fit(x[mask_red2], y[mask_red2], z, color_cut=color_cut, istar_dic=istar_dic)
    
    plt.text(tx, ty, ttext, fontsize=14)
        
     
    plt.plot(mags[(counts>50)&(mags<istar_dic[z]+dmag)], colors[(counts>50)&(mags<istar_dic[z]+dmag)], 'ro')
    
    plt.plot(xfit, yfit, 'r-', lw=1.5)
    tx = (min(plot_range[0])+(max(plot_range[0])-min(plot_range[0]))*0.05)
    ty = (min(plot_range[1])+(max(plot_range[1])-min(plot_range[1]))*0.04)
    ttext='red seq. slope: %.2f'%a 
    plt.text(tx, ty, ttext, fontsize=14)
    red_slope = a
                
                
    if weights is not None:
        #print(len(weights))
        #print(len(mask_red))
        red_frac = np.sum(weights[mask_red])/np.sum(weights)
    else:
        red_frac = float(len(x[mask_red]))/float(len(x))
    
    tx = (min(plot_range[0])+(max(plot_range[0])-min(plot_range[0]))*0.05)
    ty = (min(plot_range[1])+(max(plot_range[1])-min(plot_range[1]))*0.12)
    ttext='red fraction: %.2f'%red_frac 
    plt.text(tx, ty, ttext, fontsize=14)

        
    return p, red_frac, red_slope

def plot_loop(vacs, x, y, z_low, z_up, color_cut, x_range, y_range, titles, istar_dic=False):
    #-------------------------------------------------#
    plt.figure(dpi=300,figsize=[16,4])    
    #-------------------------------------------------#
    for j, df in enumerate(vacs):
        mask = ((df['z_best']>z_low)&(df['z_best']<=z_up)&
                (df[x]>min(x_range))&(df[x]<max(x_range))&
                (df[y]>min(y_range))&(df[y]<max(y_range)))
        try:
            p, red_frac, red_slope = cmd_plot(df[x][mask], df[y][mask], panel=j+1, istar_dic=istar_dic, bins=[100,100], #delta_x*15.,delta_y*50.], 
                             plot_range=[x_range,y_range], weights=None, cmin=1, cmax=None, z_range=(z_low,z_up), title=titles[j], x_label=x, y_label=y, color_cut=color_cut, dmag=1.5)
            if j == 0: 
                frac[np.mean([z_low,z_up])] = red_frac
                slope[np.mean([z_low,z_up])] = red_slope
        
        except:
            pass
            print('plot fail')

    #-------------------------------------------------#
    plt.tight_layout()

def CCD(mag1, mag2, mag3, mag4, name_mag1, name_mag2, name_mag3, name_mag4, maglim):
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(8, 4), dpi=150)
    bright = (mag1 <= maglim)
    mag1 = mag1[bright]
    mag2 = mag2[bright]
    mag3 = mag3[bright]
    mag4 = mag4[bright]
    x = mag1-mag2
    y = mag2-mag3
    z = mag3-mag4
    
    im1 = ax1.hexbin(x, y, None,  mincnt=1, cmap='rainbow', gridsize=[200,100])
    im2 = ax2.hexbin(z, y, None,  mincnt=1, cmap='rainbow', gridsize=[200,100])

    ax1.set_xlim([-1, 2.2])
    ax1.set_ylim([-1, 2.2])
    ax1.set_xlabel(f'{name_mag1} - {name_mag2}')
    ax1.set_ylabel(f'{name_mag2} - {name_mag3}')
    ax1.grid(True)

    ax2.set_xlim([-1, 2.2])
    ax2.set_ylim([-1, 2.2]) 
    ax2.set_xlabel(f'{name_mag3} - {name_mag4}')
    ax2.set_ylabel(f'{name_mag2} - {name_mag3}')
    ax2.grid(True)
    
    cbaxes = fig.add_axes([1.0, 0.23, 0.015, 0.7])
    cbar = fig.colorbar(im1, cax=cbaxes, cmap='rainbow', orientation='vertical')

    plt.suptitle('         Color-color diagram', fontsize=15)
    plt.tight_layout()