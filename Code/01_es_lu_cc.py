#------------------------------------------------------------------------------
# Name:        01_es_lu_cc.py
# Purpose:     Processing for the CREW project on ES, LUC and CC.
#
# Author:      James Sample
#
# Created:     14/01/2015
# Copyright:   (c) James Sample and JHI, 2015
#------------------------------------------------------------------------------
""" Processes the Future Flows (FF) climate data and estimate climate and land 
    use change effects on Ecosystem Services (ES). Reads workshop outputs and
    performs the following steps:    
    
        1. For each ES, reads monthly rainfall and ET grids for the months
           specified for both baseline and future periods. For the seasons of
           interest, calculates the % change in rainfall and ET between 
           baseline and future.
        
        2. Combines rainfall and runoff percentage changes into a qualitative 
           grid of change in runoff.
           
        3. Estimates impacts grids for each ES for CC only, LUC only and CC &
           LUC combined.

    Inputs grids are supplied in HDF5 file format.
"""

import pandas as pd, h5py, numpy as np, matplotlib, matplotlib.pyplot as plt
import os, sys
from mpl_toolkits.axes_grid1 import ImageGrid
from osgeo import gdal, gdalconst, osr

def read_array_from_h5(h5, variable, model, year, month):
    """ Read an array from a specified location in an H5 file.
    
    Args:
        h5:       The open HDF5 file object
        variable: The variable of interest ('rainfall' or 'pet')
        model:    The code for the climate model of interest (string)
        year:     Year (integer)
        month:    Month (integer)
    
    Returns:
        array
    """
    dset_path = r'/ff_data/%s/%s/%s_%s' % (model, variable, variable, year)
    data = h5.get(dset_path)[:,:,month-1].astype(float)
    
    # Set NoData to NaN
    data[data==-99] = np.nan
    
    # Convert units
    data = data/100

    return data

def avg_rain_et(h5, st_yr, end_yr, months):
    """ Calculate average rainfall and ET grids for the specified years and 
        months.
    
    Args:
        h5:     The open HDF5 file object
        st_yr:  Start year for period of interest (integer)
        end_yr: End year for period of interest (integer)
        months: List of months of interest (integers) 
    
    Returns:
        Tuple of arrays (average rainfall, average PET)
    """
    # Empty arrays to store rainfall and ET totals
    rn_tot = np.zeros((715, 485))
    et_tot = np.zeros((715, 485))
    
    # Total number of years to average over
    years = end_yr + 1 - st_yr
    
    # Loop over rainfall and ET
    for year in range(st_yr, end_yr+1):
        for month in months:
            # Read rainfall and ET grids
            rn = read_array_from_h5(h5, 'rainfall', model, year, month)
            et = read_array_from_h5(h5, 'pet', model, year, month)
            
            # Add to totals
            rn_tot += rn
            et_tot += et
    
    # Average            
    rn_av = rn_tot/years
    et_av = et_tot/years

    return (rn_av, et_av)

def plot_avg_grids(base_rn_av, base_et_av, fut_rn_av, fut_et_av):
    """ Plot the average rainfall and ET grids. Used for testing.
    
    Args:
        base_rn_av: Average rainfall grid for baseline period.
        base_et_av: Average PET grid for baseline period.
        fut_rn_av:  Average rainfall grid for future period.
        fut_et_av:  Average PET grid for future period.
    
    Returns:
        None. Displays maps of each grid using same colour scale.
    """
    # Get min and max values from grids
    rnmin = min(np.nanmin(base_rn_av), np.nanmin(fut_rn_av))    
    rnmax = max(np.nanmax(base_rn_av), np.nanmax(fut_rn_av))
    
    etmin = min(np.nanmin(base_et_av), np.nanmin(fut_et_av))
    etmax = max(np.nanmax(base_et_av), np.nanmax(fut_et_av))
    
    # Plot
    fig = plt.figure()
    grid = ImageGrid(fig, 111,
                     nrows_ncols = (1, 4), 
                     axes_pad=0.5,
                     cbar_mode='each')
    
    im0 = grid[0].imshow(base_rn_av, vmin=rnmin, vmax=rnmax,
                         interpolation='nearest')
    grid.cbar_axes[0].colorbar(im0)
                   
    im1 = grid[1].imshow(fut_rn_av, vmin=rnmin, vmax=rnmax, 
                         interpolation='nearest')
    grid.cbar_axes[1].colorbar(im1)
                                     
    im2 = grid[2].imshow(base_et_av, vmin=etmin, vmax=etmax, 
                         interpolation='nearest')
    grid.cbar_axes[2].colorbar(im2)
                   
    im3 = grid[3].imshow(fut_et_av, vmin=etmin, vmax=etmax, 
                         interpolation='nearest')
    grid.cbar_axes[3].colorbar(im3)
    
    plt.show()

def plot_reclassified_grid(array, out_path, sup_title='Main title', 
                           title='Sub-title'):
    """ Plot and save the reclassified grid.
    
    Args:
        array:     Grid of integers in range -2 to +2
        out_path:  Output file path (PNG or PDF)
        sup_title: Main title for plot (string)
        title:     Sub-title for plot (string)
        
    Returns:
        None. Saves a plot to the specified path.
    """
    # Make a color map of fixed colors
    cmap = matplotlib.colors.ListedColormap(['Red', 'Orange', 'LimeGreen',
                                             'DeepSkyBlue', 'Blue'])

    bounds=[-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    
    # Create axes for plot (A4 size)
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(8.3,11.7))
    
    # Plot the array, using the colours specified
    img = axes.imshow(array, interpolation='nearest', origin='upper',
                      cmap=cmap, norm=norm)

    # Add labels to plot
    plt.title(title)
    plt.suptitle(sup_title, fontsize=16, y=0.95)
    plt.ylabel('Northing')
    plt.xlabel('Easting')
    plt.grid(True)

    # Reformat the axis labels (mainly change the Y values into northings)
    axes.set_yticks([35, 135, 235, 335, 435, 535, 635, 735])
    axes.set_yticklabels([1200, 1100, 1000, 900, 800, 700, 600, 500])
    axes.set_xticks([100, 200, 300, 400])

    # Add axes for the color bar
    cax = fig.add_axes([0.2, 0.785, 0.02, 0.10])

    # Add the colour bar and set labels
    cbar = fig.colorbar(img, cax=cax, cmap=cmap, norm=norm, boundaries=bounds,
                        ticks=[-2.2,-1.2,-0.2,0.8,1.8])

    cbar.set_ticklabels(['Large decrease',
                         'Small decrease',
                         'Neutral',
                         'Small increase',
                         'Large increase'], update_ticks=True)

    # Make the cbar ticks invisible
    ticks = cbar.ax.get_yticklines()
    for tick in ticks:
        plt.setp(tick, alpha=0)

    cbar_labels = plt.getp(cbar.ax.axes, 'yticklabels')
    plt.setp(cbar_labels, fontsize=10) 

    # Save fig
    plt.savefig(out_path, dpi=300)
##    plt.show()
    plt.clf()
    plt.close()

def reclass_rn_et_grid(array):
    """ Take an array of percentage changes and reclassify it according to:
    
            % change  | Class
             x<=-15   |  -2
            -15<x<=-5 |  -1
            -5<x<=5   |   0
             5<x<=15  |  +1
             15<x     |  +2 
    
    Args:
        array: Array of percentage changes to be reclassified.
        
    Returns:
        Reclassified array
    """
    # Create copy of array for reclass values         
    rc = array.copy()
    rc[array<=-15] = -2
    rc[(-15<array) & (array<=-5)] = -1
    rc[(-5<array) & (array<=5)] = 0
    rc[(5<array) & (array<=15)] = 1
    rc[15<array] = 2
    
    return rc

def reclass_ro(matrix_path, rn, et):
    """ Generate reclassification matrix for runoff based on reclassified 
        change grids for rainfall and PET and the runoff reclassification 
        matrix from the workshop.
    
    Args:
        matrix_path: Path to CSV file representing runoff matrix.
        rn:          Reclassified rainfall grid from reclass_rn_et_grid
        et:          Reclassified PET grid from reclass_rn_et_grid
    
    Returns:
        Array (grid of integers representing change in runoff)
    """
    # Read matrix
    df = pd.read_csv(matrix_path, index_col=0)
    
    # Grid of NaNs wih correct shape
    ro = rn.copy()*np.nan
    
    # Loop over inidces
    for x, y in np.ndindex(ro.shape):
        # Get values for change in rainfall and ET
        et_ch = et[x, y]
        rn_ch = rn[x, y]
        
        # If both are not nan, reclassify
        if (np.isfinite(et_ch) and np.isfinite(rn_ch)):
            rc_val = df.ix[int(et_ch), str(int(rn_ch))]
            ro[x, y] = rc_val
    
    return ro

def reclass_es_ro(es_idx, ro):
    """ Reclassify the runoff grid to estimate effects of runoff change on each
        ES.
    
    Args:
        es_idx: The ID of the ES of interest in data frame ro_df
        ro:     The runoff change grid from reclass_ro
        
    Returns:
        Array (grid of integers representing change in ES) 
    """
    # Make a copy of the ro grid to update
    es = ro.copy()
    
    # Reclassify
    for chng in [-2, -1, 0, 1, 2]:
        es[ro==chng] = ro_df.ix[es_idx, 'RO_%d' % chng]
    
    return es

def read_ascii(ascii_path,
               xmin=0,
               xmax=485000,
               ymin=520000,
               ymax=1235000,
               exptd_rows=715,
               exptd_cols=485,
               exptd_px_wd=1000,
               exptd_px_ht=-1000,
               exptd_ndv=-9999):
    """ Read an ASCII grid file, clip it to the specified bounding box and
        return a numpy array.
    
    Args:
        xmin:         Minimum Easting in OSGB1936 metres.
        xmax:         Maximum Easting in OSGB1936 metres.
        ymin:         Minimum Northing in OSGB1936 metres.
        ymax:         Maximum Northing in OSGB1936 metres.
        exptd_rows:   No. of rows expected in file.
        exptd_cols:   No. of columns expected in file.
        exptd_px_wd:  Cell width.
        exptd_px_ht:  Cell height.
        exptd_ndv:    No data value.
        
    Returns:
        Array (floats).
    """
    # Register drivers
    gdal.AllRegister()

    # Process the file with GDAL
    ds = gdal.Open(ascii_path, gdalconst.GA_ReadOnly)
    if ds is None:
    	print 'Could not open ' + ascii_path
    	sys.exit(1)

    # In order to select the first cell correctly, choose a point just within
    # the top left corner of the specified bounding box.
    x = xmin + 10
    y = ymax - 10

    # Dataset properties
    geotransform = ds.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]

    # Calculate number of rows and cols to return
    rows = abs(int((ymax-ymin)/pixelHeight))
    cols = int((xmax-xmin)/pixelWidth)

    # Select starting pixel
    xOffset = int((x - originX) / pixelWidth)
    yOffset = int((y - originY) / pixelHeight)

    band = ds.GetRasterBand(1)
    no_data_val = band.GetNoDataValue()

    # Simple checking
    assert rows == exptd_rows
    assert cols == exptd_cols
    assert pixelWidth == exptd_px_wd
    assert pixelHeight == exptd_px_ht
    assert no_data_val == exptd_ndv

    # Read the data to an array
    data = band.ReadAsArray(xOffset, yOffset, cols, rows)

    # Close the dataset
    ds = None

    return data.astype(float)

def process_land_use_change(lu_mat_path, base, fut, esid, codes_df):
    """ Estimate land use change (LUC) only effects for the specified ES.
    
    Args:
        lu_mat_path: Excel file containing land use matrices from the workshop.
        base:        Baseline land luse grid.
        fut:         Future land luse grid.
        esid:        ES ID from land use matrices Excel file
        codes_df:    Land use code look-up table (as data frame)
        
    Returns:
        Array (grid of integers representing change in ES)  
    """
    # Read matrix for this ES
    lu_mat = pd.read_excel(lu_mat_path, sheetname='Land Use')

    # Get row for start of matrix
    st_row = (lu_mat['ES_ID']==esid).nonzero()[0][0] + 2

    # Read matrix of interest
    lu_mat = pd.read_excel(lu_mat_path, sheetname='Land Use', skiprows=st_row,
                           skip_footer=(120-6-st_row), parse_cols='C:I',
                           index_col=0)

    # Perform reclassification
    # Grid of NaNs wih correct shape
    rc = base.copy()*np.nan
    
    # Loop over inidces
    for x, y in np.ndindex(base.shape):
        # Get values for baseline and future LU
        base_lu = base[x, y]
        fut_lu = fut[x, y]
        
        # If both are not nan, reclassify
        if (np.isfinite(base_lu) and np.isfinite(fut_lu)):
            # Get the base and fut LU as a string
            base_str = codes_df.ix[int(base_lu)]['LU_Class']
            fut_str = codes_df.ix[int(fut_lu)]['LU_Class']
            rc_val = lu_mat.ix[base_str, fut_str]
            rc[x, y] = rc_val
    
    return rc

def process_land_use_and_climate_change(lucc_mat_path, lugrid, ccgrid, esid):
    """ Estimate combined land use and climate change effects for the specified
        ES.
        
    Args:
        lucc_mat_path: Excel file containing matrices from the workshop.
        lugrid:        The grid of land use change effects.
        ccgrid:        The grid of climate change effects.
        esid:          ES ID from workshop matrices Excel file.
        
    Returns:
        Array (grid of integers representing change in ES) 
    """
    # Read matrix for this ES
    lucc_mat = pd.read_excel(lucc_mat_path, sheetname='CC_LU')

    # Get row for start of matrix
    st_row = (lucc_mat['ES_ID']==esid).nonzero()[0][0] + 2
    
    # Read matrix of interest
    lucc_mat = pd.read_excel(lucc_mat_path, sheetname='CC_LU', skiprows=st_row,
                             skip_footer=(108-5-st_row), parse_cols='C:I',
                             index_col=0)

    # Perform reclassification
    # Grid of NaNs wih correct shape
    rc = lugrid.copy()*np.nan
    
    # Loop over inidces
    for x, y in np.ndindex(lugrid.shape):
        # Get values for baseline and future LU
        lu = lugrid[x, y]
        cc = ccgrid[x, y]
        
        # If both are not nan, reclassify
        if (np.isfinite(lu) and np.isfinite(cc)):
            # Get the base and fut LU as a string
            rc_val = lucc_mat.ix[int(lu), int(cc)]
            rc[x, y] = rc_val
    
    return rc

def array_to_gtiff(out_path, data_array, ndv=-9999, xmin=0, ymax=1235000, 
                   cell_size=1000):
    """ Convert numpy array to 16-bit integer GeoTiff.
    
    Args:
        out_path:   The .tif file to be created.
        data_array: The (integer) data array to save.
        ndv:        No data value.
        xmin:       Minimum x (Easting) co-ordinate, in OSGB1936 metres
        ymax:       Maximim y (Northing) co-ordinate, in OSGB1936 metres
        cell_size:  Cell size (metres)
    
    Returns:
        None. Array is saved to specified path.
    """
    # Copy data_array so that it is not modified
    data = data_array.copy()
    
    # Convert NaNs to NDV
    data[np.isnan(data)] = ndv
    
    # Get array shape
    cols = data.shape[1]
    rows = data.shape[0]

    # Get driver
    driver = gdal.GetDriverByName('GTiff') # NB can't directly create ArcInfo ASCII grids in this way

    # Create a new raster data source
    out_ds = driver.Create(out_path, cols, rows, 1, gdal.GDT_Int16)

    # Get spatial ref details
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(27700) # From EPSG for OSGB36 grid

    # Write metadata
    out_ds.SetGeoTransform((xmin, cell_size, 0.0, ymax, 0.0, -1*cell_size)) #(xmin, cellsize, 0, ymax, 0, -cellsize)
    out_ds.SetProjection(srs.ExportToWkt())
    out_band = out_ds.GetRasterBand(1)
    out_band.SetNoDataValue(ndv)
    out_band.WriteArray(data)

    # Tidy up
    del out_ds, out_band
    
# #############################################################################
# User input
# Climate data
ff_h5_path = r'D:\WBM_Development_2014\WBM_2014_Monthly_Input_File.h5'

# Runoff matrices
ro_path = r'D:\Eco_Services_Impacts\Matrices_Development\03_Group_1_Matrices\Runoff_Impacts_Grp1.csv'
ro_matrix_15 = r'D:\Eco_Services_Impacts\Matrices_Development\02_Common_Matrices\Runoff_Matrix_15pct.csv'

# Land use data
base_path = r'D:\Eco_Services_Impacts\Land_Use\baseline_lu_lcm07.txt'
fut_path = r'D:\Eco_Services_Impacts\Land_Use\future_lu_2050.txt'

# Land use matrices
lu_classes_path = r'D:\Eco_Services_Impacts\Land_Use\Land_Use_Classes.csv'
lu_matrices_path = r'D:\Eco_Services_Impacts\Matrices_Development\03_Group_1_Matrices\Land_Use_Matrices_Grp1.xlsx'

# Land use and climate combined matrices
lucc_matrices_path = r'D:\Eco_Services_Impacts\Matrices_Development\03_Group_1_Matrices\Climate_And_Land_Use_Matrices_Grp1.xlsx'

# Output folders
out_pdf_fold = r'D:\Eco_Services_Impacts\Model_Output\02_Group_1_Output\PDF'
out_array_fold = r'D:\Eco_Services_Impacts\Model_Output\02_Group_1_Output\GeoTiffs'

# Time periods to compare
base_st_yr, base_end_yr = 1961, 1990
fut_st_yr, fut_end_yr = 2041, 2070

# Future Flows models of interest
models = ['afixa', 'afixc', 'afixl', 'afixm', 'afixo', 'afixh', 
          'afixi', 'afixj', 'afixk', 'afgcx', 'afixq']
# #############################################################################

# Read LU grids
base = read_ascii(base_path)
base[base==-9999] = np.nan
fut = read_ascii(fut_path)
fut[fut==-9999] = np.nan

# Read LU class codes
codes_df = pd.read_csv(lu_classes_path, index_col=0)

# Read the runoff matrices
ro_df = pd.read_csv(ro_path, index_col=0)

# Open H5 file
h5 = h5py.File(ff_h5_path, 'r')

# Iterate over each ES
for idx in ro_df.index:
    print '\nProcessing land use change impacts for %s.' % ro_df.ix[idx, 'ES']
    
    # 1. Process land use change only
    luc = process_land_use_change(lu_matrices_path, base, fut, idx, codes_df)
    
    # Prepare to save
    out_name = 'ES%02d_LUC' % idx
    
    # Save array
    out_array = os.path.join(out_array_fold, '%s.tif' % out_name)
    array_to_gtiff(out_array, luc)
    
    # Save PDF
    out_pdf = os.path.join(out_pdf_fold, '%s.pdf' % out_name)
    plot_reclassified_grid(luc, out_pdf,
                           sup_title='Change in %s' % ro_df.ix[idx, 'ES'],
                           title='(land use change only)' )

    # 2. Process climate change only    
    # Get the relevant months for this ES
    months = [int(i) for i in ro_df.ix[idx, 'Key_Months'].split(',')]

    # Loop over climate models of interest
    for model in models:
        print ('Processing climate change impacts for '
               '%s (model %s).' % (ro_df.ix[idx, 'ES'], model))
        
        # 2.1. Baseline
        base_rn_av, base_et_av = avg_rain_et(h5, base_st_yr, base_end_yr,
                                             months)

        # 2.2. Future
        fut_rn_av, fut_et_av = avg_rain_et(h5, fut_st_yr, fut_end_yr,
                                           months)

        # Plot
#        plot_avg_grids(base_rn_av, base_et_av, fut_rn_av, fut_et_av)
        
        # Calculate % change
        rn_pct = 100*(fut_rn_av - base_rn_av)/base_rn_av
        et_pct = 100*(fut_et_av - base_et_av)/base_et_av

        # Reclassify
        rn_rc = reclass_rn_et_grid(rn_pct)
        et_rc = reclass_rn_et_grid(et_pct)

#        plot_reclassified_grid(rn_rc)
#        plot_reclassified_grid(et_rc)
        
        # Generate runoff grid
        ro = reclass_ro(ro_matrix_15, rn_rc, et_rc)
        
#        # Plot runoff grid        
#        plot_reclassified_grid(ro, 
#                               sup_title='Change in runoff',
#                               title='(Model %s; %s)' % (model, months))
        
        # Reclass ro grid to estimate ES impact
        es = reclass_es_ro(idx, ro)

        # Prepare to save
        out_name = 'ES%02d_%s' % (idx, model)
        
        # Save array
        out_array = os.path.join(out_array_fold, '%s.tif' % out_name)
        array_to_gtiff(out_array, es)
        
        # Save PDF
        out_pdf = os.path.join(out_pdf_fold, '%s.pdf' % out_name)   
        plot_reclassified_grid(es, out_pdf,
                               sup_title='Change in %s' % ro_df.ix[idx, 'ES'],
                               title='(climate model %s only)' % model)
        
        # 3. Process combined land use and climate effects
        print ('Processing climate and land use change impacts for '
               '%s (model %s).' % (ro_df.ix[idx, 'ES'], model))
               
        # Reclassify to get CC and LUC effects
        cc_lu = process_land_use_and_climate_change(lucc_matrices_path, luc,
                                                    es, idx)

        # Prepare to save
        out_name = 'ES%02d_LUC_%s' % (idx, model)
        
        # Save array
        out_array = os.path.join(out_array_fold, '%s.tif' % out_name)
        array_to_gtiff(out_array, cc_lu)
        
        # Save PDF
        out_pdf = os.path.join(out_pdf_fold, '%s.pdf' % out_name)   
        plot_reclassified_grid(cc_lu, out_pdf,
                               sup_title='Change in %s' % ro_df.ix[idx, 'ES'],
                               title='(climate and land use change together)')
                               
# Close H5 file
h5.close()

print '\nFinished.'            