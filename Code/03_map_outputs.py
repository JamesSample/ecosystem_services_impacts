#------------------------------------------------------------------------------
# Name:        03_map_outputs.py
# Purpose:     Uses ArcPy to display the results from the qualitative modelling
#
# Author:      James Sample
#
# Created:     
# Copyright:   (c) James Sample and JHI, 2015
#------------------------------------------------------------------------------
""" Takes outputs from 02_wb_zonal_stats.py and uses them to produce sets of
    maps showing (i) median modelled impact on each ecosystem service and (ii)
    model confidence based on the range of modelled scores in each waterbody.
    
    Note: This code requires ArcPy and an ArcInfo licence. Must be run using
    32-bit Python to avoid DLL errors from ArcGIS.
"""
import arcpy, os, pandas as pd

# #############################################################################
# User input
template_path = r'D:\Eco_Services_Impacts\Model_Output\04_Overall_Output\MXD\Template.mxd'
shp_fold = r'D:\Eco_Services_Impacts\Model_Output\04_Overall_Output\Shapefiles'
mxd_fold = r'D:\Eco_Services_Impacts\Model_Output\04_Overall_Output\MXD'
pdf_fold = r'D:\Eco_Services_Impacts\Model_Output\04_Overall_Output\PDF'
png_fold = r'D:\Eco_Services_Impacts\Model_Output\04_Overall_Output\PNG'
es_csv = r'D:\Eco_Services_Impacts\SEPA_Data\ES_List_Simple2.csv'

# .lyr files
luc_med = r'D:\Eco_Services_Impacts\Model_Output\04_Overall_Output\LYR\LUC_Med.lyr'
luc_range = r'D:\Eco_Services_Impacts\Model_Output\04_Overall_Output\LYR\LUC_Range.lyr'
cc_med = r'D:\Eco_Services_Impacts\Model_Output\04_Overall_Output\LYR\CC_Med.lyr'
cc_range = r'D:\Eco_Services_Impacts\Model_Output\04_Overall_Output\LYR\CC_Range.lyr'
cc_luc_med = r'D:\Eco_Services_Impacts\Model_Output\04_Overall_Output\LYR\Comb_Med.lyr'
cc_luc_range = r'D:\Eco_Services_Impacts\Model_Output\04_Overall_Output\LYR\Comb_Range.lyr'
# #############################################################################

# Dict of lyr files
lyr_dict = {'LUC':[luc_med, luc_range, 'land use change only'],
            'CC':[cc_med, cc_range, 'climate change only'],
            'Comb':[cc_luc_med, cc_luc_range, 
                    'climate and land use change combined']}

# Read list of ES into df
es_df = pd.read_csv(es_csv, index_col=0)

# Loop over ES
for esid in es_df.index:
    # Loop over scenarios
    for scen in lyr_dict.keys():
        # Get the name of the current es
        es_text = es_df['ES'].ix[esid]

        print 'Currently processing: %s (%s).' % (es_text, lyr_dict[scen][2])
        
        # Open the template
        temp_mxd = arcpy.mapping.MapDocument(template_path)

        # Save a copy of the template to modify
        mxd_path = os.path.join(mxd_fold, 'ES_%02d_%s.mxd' % (esid, scen))
        temp_mxd.saveACopy(mxd_path)
        
        # open the new mxd
        mxd = arcpy.mapping.MapDocument(mxd_path)
    
        # Modify the map title
        title_text_el = arcpy.mapping.ListLayoutElements(mxd, "TEXT_ELEMENT")[-1]
        title_text_el.text = '%s: %s' % (es_text, lyr_dict[scen][2])
            
        # Get the shapefile
        shp = os.path.join(shp_fold, 'ES_%02d.shp' % esid)
        
        # Create layer
        lyr = arcpy.mapping.Layer(shp)
        
        # Get the median data frame
        df = arcpy.mapping.ListDataFrames(mxd, "Median")[0]  
        
        # Add layer to df
        arcpy.mapping.AddLayer(df, lyr, "TOP")
        
        # Set layer symbology from .lyr file
        med_lyr = arcpy.mapping.ListLayers(mxd, 'ES_%02d' % esid, df)[0]
        sym_lyr = arcpy.mapping.Layer(lyr_dict[scen][0])
        arcpy.mapping.UpdateLayer(df, med_lyr, sym_lyr, True)
        
        # Get the confidence data frame
        df = arcpy.mapping.ListDataFrames(mxd, "Confidence")[0]  
        
        # Add layer to map
        arcpy.mapping.AddLayer(df, lyr, "TOP")
        
        # Set layer symbology from .lyr file
        con_lyr = arcpy.mapping.ListLayers(mxd, 'ES_%02d' % esid, df)[0]
        sym_lyr = arcpy.mapping.Layer(lyr_dict[scen][1])
        arcpy.mapping.UpdateLayer(df, con_lyr, sym_lyr, True)
        
        # Refresh map
        arcpy.RefreshActiveView()
        arcpy.RefreshTOC()
        mxd.save()
        
        # Save to PDF
        pdf_path = os.path.join(pdf_fold, 'ES_%02d_%s.pdf' % (esid, scen))
        arcpy.mapping.ExportToPDF(mxd, pdf_path)
        
        # Save to PNG
        png_path = os.path.join(png_fold, 'ES_%02d_%s.png' % (esid, scen))
        arcpy.mapping.ExportToPNG(mxd, png_path, resolution=300)
        
        del mxd, temp_mxd

print 'Finished.'