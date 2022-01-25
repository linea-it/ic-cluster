import os
import numpy as np
import healpy as hp
from astropy.table import Table
from collections import OrderedDict
from gavodb import DBManager
import sqlalchemy
import pylab as pl
import seaborn as sns
import pandas as pd
# seaborn plot setup
sns.set(color_codes=True, font_scale=1.5) 


###########################
##### Save/Write ##########
###########################
def compact_map_save(hpmap, outname, exclude=lambda x:x>0):
    """
    Save maps in a pixel, signal format

    Parameters
    ----------
    hpmap: healpix map
        Healpix map
    outname: str
        Output name
    exclude: function
        Function to mask pixels to be excluded, must have hpmap as argument

    Returns
    -------
    None
    """
    mask = exclude(hpmap)
    Table([np.arange(hpmap.size, dtype=int)[mask], hpmap[mask]],
             names=['pixel', 'signal']).write(outname)
    return
def save_maps(maps, DIR=None, exclude=lambda x:x>0):
    """
    Save all maps in a pixel, signal format

    Parameters
    ----------
    maps: astropy Table
        Table with maps to be saved, with colnames as the output names
    DIR: str, None
        Output directory, if None maps are saved localy
    exclude: function
        Function to mask pixels to be excluded, must have hpmap as argument

    Returns
    -------
    None
    """
    if DIR is not None:
        os.system(f'mkdir {DIR}')
    else:
        DIR = '.'
    for col in maps.colnames:
        compact_map_save(maps[col], f'{DIR}/{col}.fits',
                         exclude=exclude)
    return
def pixel_to_map(pixel, signal, nside=4096, outvalue=0):
    """
    Convert pixel, signal to map

    Parameters
    ----------
    pixel: array
        Healpy pixels
    signal: array
        Value of map
    nside: int
        Healpix NSIDE
    outvalue: float
        Value for pixels not in the map

    Returns
    -------
    out_map: numpy array
        Map
    """
    out_map = np.zeros(12*nside**2)
    out_map[:] = outvalue
    out_map[np.array(pixel, dtype=int)] = signal
    return out_map
def compact_map_read(mapname, nside=4096, outvalue=0):
    """
    Reads map from pixel, signal format

    Parameters
    ----------
    mapname: str
        Name of map file
    nside: int
        Healpix NSIDE
    outvalue: float
        Value for pixels not in the map

    Returns
    -------
    numpy array
        Map
    """
    data = Table.read(mapname)
    return pixel_to_map(data['pixel'], data['signal'],
                        nside=nside, outvalue=outvalue)
def read_maps(colnames='griz', nside=4096, DIR=None, outvalue=0):
    """
    Reads all maps in a pixel, signal format

    Parameters
    ----------
    colnames: list
        List of names of each map, will became the names in output table
    nside: int
        Healpix NSIDE
    DIR: str, None
        Maps directory, if None, maps are searched localy
    outvalue: float
        Value for pixels not in the map

    Returns
    -------
    maps: astropy Table
        Table with maps
    """
    if DIR is None:
        DIR = '.'
    maps = Table()
    for col in colnames:
        maps[col] = compact_map_read(f'{DIR}/{col}.fits',
                                     nside=nside, outvalue=outvalue)
    return maps
###########################
##### Specific Inputs #####
###########################
def get_depthmaps(pid='6005', bands = 'griz', DIR='pid',
                 save_maps=True):
    """
    Get Depth Maps

    Parameters
    ----------
    pid: str
        Process id
    bands: list, tuple, str
        Band of the depthmap
    DIR: str
        Output directory, if 'pid' assumes the value of pid
    save_maps: bool
        Saves the maps localy when read from DB
    """
    if DIR is None:
        DIR = '.'
    if DIR == 'pid':
        DIR = pid
    if np.array([os.path.isfile(f'{DIR}/{b}.fits')
                 for b in bands]).all():
        maps = read_maps(colnames=bands, nside=4096,
                         DIR=DIR, outvalue=0)
    else:
        dbm = DBManager()
        tables_dict = dbm.get_unique_band_table(pid, bands,
                                                find_format='_nest_%s_')
        print(tables_dict)
        maps = Table()
        for b in bands:
            print(f'*** {b} ***')
            maps[b] = pixel_to_map(*dbm.get_db_table(tables_dict[b], ['pixel, signal']),
                                   nside=4096, outvalue=0)
        del data
        print('ok')
        if save_maps:
            compact_save(maps,  DIR=pid)
    return maps
def get_photoz_valid(pid='6036'):
    """
    Get All validation files

    Parameters
    ----------
    pid: str
        Photo-z comput process id

    Returns
    -------
    Dict
        Dict of astropy tables with photo-z code as keys
    """
    dbm = DBManager()
    pype_input = dbm.get_pype_input(pid)

    # Get training files from DB
    pid_training = [i.split('process_id="')[1].split('"')[0]
                    for i in pype_input.split('\n') if '"Photo-z Training"' in i][0]
    training_files = dbm.get_output_files(pid_training)

    # Read files into dict
    validation_tables = {}
    for f in training_files:
        if f[-15:]=='_tsm_valid.fits':
            name = f.split('photoz_valid_1_._')[1].split('_tsm_valid.fits')[0]
            validation_tables[name] = Table.read(f'/process/production/{f}')
            for c in validation_tables[name].colnames:
                validation_tables[name].rename_column(c, c.lower())
    return validation_tables

def get_photoz_comp(pid='6036',
    columns=['z_best', 'err_z'],
    bands=['g', 'r', 'i', 'z'],
    zmin=0, zmax=2,
    magmin=0, magmax=90,
    sample_size = 100000):
    """
    Get subsample of photo-z computation

    Parameters
    ----------
    pid: str
        Process id
    columnst: list
        List of columns to extract from pz table
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
    sample_size: int
        Number of objects to be selected

    Returns
    -------
    Dict
        Dict of astropy tables with photo-z code as keys
    """
    dbm = DBManager()
    config = dbm.get_config(pid)
    pype_input = dbm.get_pype_input(pid)

    # Get coadd table name
    coadd_pid = dbm.get_pid_in_xml(pype_input, 'coadd')
    print(coadd_pid)
    coadd_tb = dbm.get_tablelist_from_pid(coadd_pid)#[0]
    print(coadd_pid)
    
    # mag type
    training_pid = dbm.get_pid_in_xml(pype_input, '"Photo-z Training"')
    training_config = dbm.get_config(training_pid)
    mag_type = dbm.get_value_in_xml(training_config, 'photo_type"').lower()
    mag_fmt = {'auto':'mag_auto_%s',
               'mof_bdf_mag_%s': 'mof_bdf_mag_%s',
               'sof_bdf_mag_%s': 'sof_bdf_mag_%s',
               'mof_bdf_mag_%s_corrected': 'mof_bdf_mag_%s_corrected',
               'sof_bdf_mag_%s_corrected': 'sof_bdf_mag_%s_corrected',
              }[mag_type]
    mags = [mag_fmt%b for b in bands]

    # Construct query
    mag_filter = 'and '.join([f'{m}<{magmax:g} and {m}>{magmin:g}' for m in mags])
    z_filter = f'z_best>={zmin:g} and z_best<={zmax:g}'

    # Get tables with mags
    colnames = columns+[f'mag_{b}' for b in bands]
    pz_tables = {}
    for tb in dbm.get_tablelist_from_pid(pid):
        name = tb.split('.')[1].split('_')[0]
        query = f'(select * from {tb} where {z_filter} limit {sample_size}) pz'
        query += f' inner join {coadd_tb} c on c.coadd_objects_id=pz.coadd_objects_id'
        query += f' where {mag_filter}'
        pz_tables[name] = Table(list(dbm.get_db_table(query, columns+mags)),
                                names=colnames)
    return pz_tables

def get_specz_sample(pid='6315', column_names=('ra', 'dec', 'z', 'source')):
    """
    Get spectroscopic sample from DES Science Portal products

    Parameters
    ----------
    pid: str
        Process id
        Default: pid='6315' # Version of 23 Nov 2020
    column_names: tuple
        Names of columns to extract from spec-z table.
        Default: columns=('ra', 'dec', 'z', 'source')
        Samples built after Nov 2020 can also include the column 'class' (str, STAR/GALAXY/QSO).      
    
    Returns
    -------
    Astropy table with selected columns.
    """
    
    
    # Connect to database    
    engine = sqlalchemy.create_engine('postgres://untrustedprod:untrusted@desdb4.linea.gov.br:5432/prod_gavo')
    conn = engine.connect()
     
    table_name = f"centralized_spec_db.spec_database_{pid}"
    query = f"select {(',').join(column_names)} from {table_name}"
    stm = sqlalchemy.sql.text(query)
    result = conn.execute(stm).fetchall()  
    specz_table = Table(rows=result, names=column_names)
   
    return specz_table


def get_train_set(pid='6322', column_names=('ra', 'dec', 'z', 'source')):
    """
    Get training set from DES Science Portal products

    Parameters
    ----------
    pid: str
        Process id
        Default: pid='6322' # Version of 23 Nov 2020 (Y3 Gold Full)
    column_names: tuple
        Names of columns to extract from training set table.
        Default: columns=('ra', 'dec', 'z', 'source')
    
    Returns
    -------
    Astropy table with selected columns.
    """
    
    
    # Connect to database    
    engine = sqlalchemy.create_engine('postgres://untrustedprod:untrusted@desdb4.linea.gov.br:5432/prod_gavo')
    conn = engine.connect()
     
    table_name = f"training_set_maker.training_set_{pid}"
    query = f"select {(',').join(column_names)} from {table_name}"
    stm = sqlalchemy.sql.text(query)
    result = conn.execute(stm).fetchall()  
    train_set_table = Table(rows=result, names=column_names)
   
    return train_set_table




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




def plot_specz_spatial_dist(ra, dec, foot_ra, foot_dec, 
                            projection='mollweide', cmap='viridis', 
                            cbar_label='log(#spec-$z$/deg$^{2}$)'): 
    """
    Plot projected spatial density of spec-zs 
    
    Parameters
    ----------
    ra: numpy array 
        Right ascention in sterradians
    dec: numpy array 
        Declination in sterradians
    
    """
    
    fig = pl.figure(figsize=[14,6], dpi=300)
    ax = fig.add_subplot(111, projection=projection)   
    pl.hexbin(ra, dec, None,  mincnt=1,
              gridsize=[180,90], bins='log', cmap=cmap)
    cbar = pl.colorbar()
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label(cbar_label, size=16)
    #pl.clim(0,3.8)
    #pl.plot(-np.radians(foot_ra), np.radians(foot_dec), 'r-')
    org=0.0
    tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tick_labels = np.remainder(tick_labels+360+org,360)
    ax.set_xticklabels(tick_labels)     # we add the scale on the x axis
    ax.set_xlabel("R.A.")
    ax.xaxis.label.set_fontsize(14)
    ax.set_ylabel("Dec.")
    ax.yaxis.label.set_fontsize(14)
    ax.grid(True)
    pl.tight_layout()
    #pl.savefig('specz_GOOD_%s_spatial_dist_mollweide.png'%version)
    
    
def plot_dist_surveys(n_obj_train, train_set_table, ra2, foot_ra, foot_dec):
    for survey in sorted(n_obj_train, key=n_obj_train.get, reverse=True):
        if n_obj_train[survey] > 1000:
            pl.figure(figsize=[14,4])
            pl.subplot(121) 
            sns.kdeplot(train_set_table['z'][train_set_table['source'] == survey], shade=True, 
                    label=f'{survey.strip()}: {n_obj_train[survey]} objects')
            pl.xlabel("spec-$z$", fontsize=13)
            pl.ylabel("counts", fontsize=14)
            pl.xlim(0,2.5)
            pl.subplot(122)
            pl.plot(ra2[train_set_table['source']  == survey], train_set_table['dec'][train_set_table['source']  == survey], '.')
            pl.plot(foot_ra, foot_dec, 'k-')
            pl.xlabel('R.A. (deg)')
            pl.ylabel('Dec. (deg)')
            pl.xlim(180,-180)
            pl.tight_layout()
            
def plot_dist_mag(train_set_table, dataset):
    pl.figure(figsize=[14,4])#, dpi=300)
    if "DEEP" in dataset:
        mag_name = "sof_bdf_mag_i_corrected"
    else:
        mag_type = "mag_auto"
    pl.subplot(121)
    mag_i = train_set_table[mag_name]
    mask = (mag_i > 12)*(mag_i < 28)
    sns.kdeplot(mag_i[mask], shade=True, label='matched')
    #pl.xlim(0,1.5)
    pl.xlabel('i-band MAG AUTO')
    pl.ylabel("counts", fontsize=14)
    pl.subplot(122)
    pl.plot(train_set_table['z'], train_set_table[mag_name], ',')
    pl.xlabel('redshift', fontsize=16)
    pl.ylabel(mag_name, fontsize=16)
    pl.xticks(fontsize=14)
    pl.yticks(fontsize=14)
    pl.xlim(0,1.8)
    pl.ylim(14,26)
    pl.tight_layout()
