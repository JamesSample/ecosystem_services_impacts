#------------------------------------------------------------------------------
# Name:        02_wb_zonal_stats.py
# Purpose:     Performs zonal statistics on overlapping polygons for the 
#              CREW ES, LU and CC project.
#
# Author:      James Sample
#
# Created:     20/01/2015
# Copyright:   (c) James Sample and JHI, 2015
# License:     https://github.com/JamesSample/ecosystem_services_impacts/blob/master/LICENSE
#------------------------------------------------------------------------------
""" Takes the output from 01_es_lu_cc.py and performs zonal statistics for each
    of SEPA's nested water bodies.
    
    Modified from original code here: http://www.gdal.org/ogr/drv_filegdb.html

    (Zonal Statistics, Vector-Raster Analysis, Copyright 2013 Matthew Perry)
"""

import gdal, ogr, numpy as np, pandas as pd, os, geopandas as gpd, time
from gdalconst import *
from fiona.crs import from_epsg
from scipy.stats.mstats import mode
#gdal.UseExceptions() 

def bbox_to_pixel_offsets(gt, bbox):
    """ Get pixel offsets.
	
	Args:
		gt:   Geotransform defining the extent of a raster.
		bbox: Bounding box marking the extent of an individual feature.
	
	Returns:
		Pixel offsets for the rasterised bounding box covering the feature.
    """
    originX = gt[0]
    originY = gt[3]
    pixel_width = gt[1]
    pixel_height = gt[5]
    x1 = int((bbox[0] - originX) / pixel_width)
    x2 = int((bbox[1] - originX) / pixel_width) + 1

    y1 = int((bbox[3] - originY) / pixel_height)
    y2 = int((bbox[2] - originY) / pixel_height) + 1

    xsize = x2 - x1
    ysize = y2 - y1
    return (x1, y1, xsize, ysize)


def zonal_stats(vector_path, zones_field, raster_path, nodata_value=None):
    """ Perform zonal stats. with overlapping polygons.
    
    Args:
        vector_path:   Path to shapefile defining zones. Zones may overlap.
        zones_field:   Field name containing unique ID for each zone.
        raster_path:   Value ratser on which statistics are based.
        nodata_value:  No data value.
    
    Returns:
        List of dictionaries containing summary statistics.
    """
    # Open the value raster
    rds = gdal.Open(raster_path, GA_ReadOnly)
    assert(rds)
    rb = rds.GetRasterBand(1)
    rgt = rds.GetGeoTransform()

    if nodata_value:
        nodata_value = float(nodata_value)
        rb.SetNoDataValue(nodata_value)

    # Open the zones dataset
    vds = ogr.Open(vector_path, GA_ReadOnly)
    assert(vds)
    vlyr = vds.GetLayer(0)

    # Get in-memory drivers for creating temporary copies of each polygon and
    # its rasterised equivalent
    mem_drv = ogr.GetDriverByName('Memory')
    driver = gdal.GetDriverByName('MEM')

    # Loop through vectors
    stats = []
    feat = vlyr.GetNextFeature()
    while feat is not None:
        # Get the extent of the current feature and read the part of the raster
        # corresponding to the feature's bounding box
        src_offset = bbox_to_pixel_offsets(rgt, feat.geometry().GetEnvelope())
        src_array = rb.ReadAsArray(*src_offset)

        # Calculate new geotransform of the feature subset
        new_gt = ((rgt[0] + (src_offset[0] * rgt[1])),
                  rgt[1],
                  0.0,
                  (rgt[3] + (src_offset[1] * rgt[5])),
                  0.0,
                  rgt[5])

        # Create a temporary vector layer in memory
        mem_ds = mem_drv.CreateDataSource('out')
        mem_layer = mem_ds.CreateLayer('poly', None, ogr.wkbPolygon)
        mem_layer.CreateFeature(feat.Clone())

        # Rasterize it, setting pixels covered by the feature = 1
        rvds = driver.Create('', src_offset[2], src_offset[3], 1, gdal.GDT_Byte)
        rvds.SetGeoTransform(new_gt)
        gdal.RasterizeLayer(rvds, [1], mem_layer, burn_values=[1])
        
        # Read as array
        rv_array = rvds.ReadAsArray()

        # Mask the source data array with our current feature
        # we take the logical_not to flip 0<->1 to get the correct mask effect
        # we also mask out nodata values explictly
        masked = np.ma.MaskedArray(src_array,
                                   mask=np.logical_or(src_array==nodata_value,
                                                      np.logical_not(rv_array)))

        # Calculate the desired statistics and append to output list
        # If the number of unmasked cells for this WB is <5, the results are
        # not reliable. Set to NaN
        if np.ma.count(masked) < min_cells:
            feature_stats = {'min': np.nan,
                             'mean': np.nan,
                             'max': np.nan,
                             'med': np.nan,
                             'mode': np.nan,
                             'std': np.nan,
                             'sum': np.nan,
                             'count': 0,
                             'fid': int(feat.GetFID()),
                             zones_field: feat.GetField(zones_field)}
        else:
            feature_stats = {'min': float(masked.min()),
                             'mean': float(masked.mean()),
                             'max': float(masked.max()),
                             'med': np.ma.median(masked),
                             'mode': mode(masked, axis=None)[0][0],
                             'std': float(masked.std()),
                             'sum': float(masked.sum()),
                             'count': int(masked.count()),
                             'fid': int(feat.GetFID()),
                             zones_field: feat.GetField(zones_field)}
        stats.append(feature_stats)

        # Tidy up and move on to next feature
        rvds = None
        mem_ds = None
        feat = vlyr.GetNextFeature()

    # Tidy up
    vds = None
    rds = None
    return stats

def process_sepa_data(sepa_csv):    
    """ The SEPA ES data includes lots of duplicates. This code does some 
        cleaning of the data. Makes the following assumptions:
                  
            1. Only interested in whether a service exists or not, not how 
               big it is. This means reclassifying the 'Total Service Provided 
               Class' column.
            
            2. Assume that blanks in the 'Total Service Provided Class' column 
               are the same as 'No service'.
            
            3. Assume that 'N/A' in the 'Total Service Provided Class' column 
               is the same as 'No service'.
               
            4. Assume that if a WB has no information at all for a particular 
               ES, then that ES is not present.
    
    Args:
        sepa_csv: CSV file summarising ES data. Based on data available from:
                  http://www.sepa.org.uk/data-visualisation/benefits-of-the-water-environment/
    
    Returns:
        Data frame of ES presence/absence
    """       
    # Read data
    df = pd.read_csv(sepa_csv)
    
    # Delete unwanted columns
    del df['Service provided'], df['Water Body Service Provided Sum'], df['ES_Name']
    
    # Drop rows without a WB_ID
    df = df.dropna(subset=['WB_ID',])
    
    # Patch blanks in 'Total Service Provided Class' column with 'No service'
    df['Total Service Provided Class'].fillna('No service', inplace=True)

    # Change 'N/A' in 'Total Service Provided Class' column to 'No service'
    df['Total Service Provided Class'].replace(to_replace='N/A', value='No service', inplace=True)
    
    # Reclassify so that areas with service are 1, otherwise 0
    df['Service'] = (df['Total Service Provided Class'] != 'No service').astype(int)
    
    del df['Total Service Provided Class']
    
    # Drop duplicates
    df = df.drop_duplicates()
   
    # Pivot
    df = df.pivot(index='WB_ID', columns='ES_ID', values='Service')
    
    # Fill NaN values with 0
    df = df.fillna(value=0)
    
    return df

def write_shapefile(shp_gdf, es_df, out_path):
    """ Join the ES data frame to the shapefile geodataframe and write the 
        result to a new shapefile.
    
    Args:
        shp_gdf:  Geodataframe of SEPA waterbody catchments
        es_df:    Data frame of ES presence/absence.
        out_path: Output file path.
        
    Returns:
        None. A shapefile is written to the specified location.
    """
    # Copy the GDF so that it is not modified
    shp = shp_gdf.copy()
    
    # Delete unwanted columns
    del shp['WB_Name'], shp['Area_km2'], shp['Status2013']
    
    # Set indices
    shp.index = shp['WB_ID']
   
    # Join to "inter" WBs
    df = shp.join(es_df, how='left')
    
    # Set spatial ref
    df.crs = from_epsg(27700)

    # Fill NoData
    df.fillna(-9999, inplace=True)
    
    # Write to output
    # Convert back to GDF (even though it already is one!) to avoid gpd errors
    df = gpd.GeoDataFrame(df)
    df.to_file(out_path)

# #############################################################################
# User input
# Path to zones dataset (e.g. shapefile)
zones_path = r'D:\Eco_Services_Impacts\SEPA_Data\GIS\SEPA_Nested_Onshore_Offshore.shp'

# The name of the field used to uniquely identify each zone
zones_field = 'WB_ID'

# Path to non-nested WBs
sepa_wbs_path = r'D:\Eco_Services_Impacts\SEPA_Data\GIS\SEPA_Onshore_Offshore_WB_Inter.shp'

# Folders of GeoTiffs
gtif_fold1 = r'D:\Eco_Services_Impacts\Model_Output\01_James_Output\GeoTiffs'
gtif_fold2 = r'D:\Eco_Services_Impacts\Model_Output\02_Group_1_Output\GeoTiffs'

# Folder for CSVs
csv_fold = r'D:\Eco_Services_Impacts\Model_Output\04_Overall_Output\CSV'

# Folder for SHPs
shp_fold = r'D:\Eco_Services_Impacts\Model_Output\04_Overall_Output\Shapefiles'

# Sepa data CSV
sepa_csv = r'D:\Eco_Services_Impacts\SEPA_Data\SEPA_ES_Data_Update_Offshore_Final.csv'

# The minimum number of cells in a WB for a result to be valid
min_cells = 5
# #############################################################################

# Output from the two expert groups
gtif_folds = [gtif_fold1, gtif_fold2] 

# Future Flows models of interest
models = ['afixa', 'afixc', 'afixl', 'afixm', 'afixo', 'afixh', 
          'afixi', 'afixj', 'afixk', 'afgcx', 'afixq']

# Read SEPA data
sepa_df = process_sepa_data(sepa_csv)

# Read the SEPA WBs to a GDF
sepa_wbs = gpd.read_file(sepa_wbs_path)          

# 1. Process land use
print 'Currently processing land use.'
    
# Empty list to store DFs
lu_list = []

# Loop over ES    
for esid in range(1, 13):
    print '    ES %02d.' % esid
    
    # Loop over group outputs
    for idx, gtif_fold in enumerate(gtif_folds):
        
        # Path to GTIFF
        gtif_path = os.path.join(gtif_fold, 'ES%02d_LUC.tif' % esid)
    
        # Get stats  
        stats = zonal_stats(zones_path, zones_field, gtif_path, 
                            nodata_value=-9999)
        
        # Convert to dataframe
        df = pd.DataFrame(stats)
        
        # Get just data for mode, indexed by WB_ID
        df = df[[zones_field, 'mode']]
        
        # Rename columns and set index
        df.columns = [zones_field, 'ES%02d_%s' % (esid, idx+1)]
        df.index = df[zones_field]
        del df[zones_field]
    
        # Append to list
        lu_list.append(df)
    
# Join
lu_df = pd.concat(lu_list, axis=1, join='outer')

# Write to output
out_csv = os.path.join(csv_fold, 'LUC_Only.csv')
lu_df.to_csv(out_csv)

# 2. Process climate
print 'Currently processing climate.'
    
# Loop over ES    
for esid in range(1, 13):
    print '    ES %02d.' % esid

    # Empty list to store DFs
    cc_list = []

    # Loop over group outputs
    for idx, gtif_fold in enumerate(gtif_folds):

        # Loop over models
        for model in models:
            
            # Path to GTIFF
            gtif_path = os.path.join(gtif_fold, 'ES%02d_%s.tif' % (esid, model))
        
            # Get stats  
            stats = zonal_stats(zones_path, zones_field, gtif_path, 
                                nodata_value=-9999)
            
            # Convert to dataframe
            df = pd.DataFrame(stats)
            
            # Get just data for mode, indexed by WB_ID
            df = df[[zones_field, 'mode']]
            
            # Rename columns and set index
            df.columns = [zones_field, '%s_%s' % (model, idx+1)]
            df.index = df[zones_field]
            del df[zones_field]
        
            # Append to list
            cc_list.append(df)
    
    # Join
    cc_df = pd.concat(cc_list, axis=1, join='outer')

    # Extract just the min, median and max scores over all simulations
    cc_df = cc_df.T.describe().T[['min', '50%', 'max']] 
    
    # Rename columns for later
    cc_df.columns = ['CC_Min', 'CC_Med', 'CC_Max']
    
    # Write to output
    out_csv = os.path.join(csv_fold, 'CC_Only_ES%02d.csv' % esid)
    cc_df.to_csv(out_csv)

# 3. Process climate and land use
print 'Currently processing climate and land use combined.'
    
# Loop over ES    
for esid in range(1, 13):
    print '    ES %02d.' % esid

    # Empty list to store DFs
    cclu_list = []

    # Loop over group outputs
    for idx, gtif_fold in enumerate(gtif_folds):
        
        # Loop over models
        for model in models:
            
            # Path to GTIFF
            gtif_path = os.path.join(gtif_fold, 'ES%02d_LUC_%s.tif' % (esid, model))
        
            # Get stats  
            stats = zonal_stats(zones_path, zones_field, gtif_path, 
                                nodata_value=-9999)
            
            # Convert to dataframe
            df = pd.DataFrame(stats)
            
            # Get just data for mode, indexed by WB_ID
            df = df[[zones_field, 'mode']]
            
            # Rename columns and set index
            df.columns = [zones_field, '%s_%s' % (model, idx+1)]
            df.index = df[zones_field]
            del df[zones_field]
        
            # Append to list
            cclu_list.append(df)

    # Join
    cclu_df = pd.concat(cclu_list, axis=1, join='outer')
    
    # Extract just the min, median and max scores over all simulations
    cclu_df = cclu_df.T.describe().T[['min', '50%', 'max']] 

    # Rename columns for later
    cclu_df.columns = ['Comb_Min', 'Comb_Med', 'Comb_Max']
    
    # Write to output
    out_csv = os.path.join(csv_fold, 'CC_and_LUC_ES%02d.csv' % esid)
    cclu_df.to_csv(out_csv)

# 4. Re-arrange the data so that each ES has a single table that includes all
# the LUC, CC and LUC+CC data for that particular service.
# This whole thing could have been coded much more neatly, but I'm just
# tagging it on here as an afterthought!
print 'Re-structuring data.'

# Read the LU CSV
lu_csv = os.path.join(csv_fold, 'LUC_Only.csv')
lu_df_all = pd.read_csv(lu_csv, index_col=0)

# Loop over ES
for esid in range(1, 13):
    # Get the LUC only data for this ES
    lu_df = lu_df_all[['ES%02d_1' % esid, 'ES%02d_2' % esid]]

    # Extract just the min, median and max scores both simulations
    lu_df = lu_df.T.describe().T[['min', '50%', 'max']] 
    
    # Rename column
    lu_df.columns = ['LUC_Min', 'LUC_Med', 'LUC_Max']
    
    # Read the CC data for this ES
    cc_csv = os.path.join(csv_fold, 'CC_Only_ES%02d.csv' % esid)
    cc_df = pd.read_csv(cc_csv, index_col=0)
    
    # Read the CC+LUC data for this ES
    cclu_csv = os.path.join(csv_fold, 'CC_and_LUC_ES%02d.csv' % esid)
    cclu_df = pd.read_csv(cclu_csv, index_col=0)
    
    # Concatenate
    df = pd.concat([lu_df, cc_df, cclu_df], axis=1, join='outer')

    # Get presence/absence info for this ES
    es_pa = sepa_df[[esid]].copy()
    
    # Set zeros to NaN
    es_pa[es_pa[esid]==0]=np.nan    
    
    # Join in the SEPA service data
    df = pd.merge(df, es_pa, left_index=True, right_index=True, 
                  how='left')
       
    # Multiply through by 'Service' columns
    for col in df.columns:
        df[col] = df[col]*df[esid]
    
    # Remove esid column
    del df[esid]

    # Round values to nearest integer
    df = np.round(df, 0)
    
    # Calculate range
    df['LUC_Range'] = df['LUC_Max'] - df['LUC_Min']
    df['CC_Range'] = df['CC_Max'] - df['CC_Min']
    df['Comb_Range'] = df['Comb_Max'] - df['Comb_Min']

    # Extract columns of interest
    df = df[['LUC_Med', 'LUC_Range', 'CC_Med', 'CC_Range', 'Comb_Med',
             'Comb_Range']]    
             
    # Explicitly set NDV to -9999
    df = df.fillna(value=-9999)
    
    # Write to output CSV
    out_path = os.path.join(csv_fold, 'ES_%02d.csv' % esid)
    df.to_csv(out_path)

    # Write to output SHP
    out_path = os.path.join(shp_fold, 'ES_%02d.shp' % esid)
    df = write_shapefile(sepa_wbs, df, out_path)
    
print 'Finished.'