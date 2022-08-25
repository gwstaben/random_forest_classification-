#!/usr/bin/env python

"""
Read in a raster image and point shapefile with classes and extracts values from the imput imagery and return a csv of the results for each band in the raster file.
Author: Grant Staben
Date: 14/10/2017
Modified 20/03/2019
"""
from __future__ import print_function, division
import pdb   
import fiona
import rasterio
import pandas as pd 
import argparse
from rasterstats import zonal_stats 
from rasterstats import point_query
import sys
import os
import shutil
import glob



def getCmdargs():

    p = argparse.ArgumentParser(description="""Input a single or multiband raster to extracts the values from the input point shapefile. The script currently outputs a csv file containing the unique identifyer for each point, the value of the raster is returned.""")

    p.add_argument("-i","--image", help="Input image to extract the stats from")
        
    p.add_argument("-n","--nodata",default=None, help="define the no data value for the input raster image, the default is none (default is %(default)s))")
    
    p.add_argument("-s","--shape", help="shape file contiaing the point files with classes, it aslo needs to have a field defined as uid")
    
    p.add_argument("-u","--uid", help="input the column name for the unique id field in the shapefile") 
    
    p.add_argument("-o","--csv", help="name of the output csv file containing the results")
    
    cmdargs = p.parse_args()
    
    if cmdargs.image is None:

        p.print_help()

        sys.exit()

    return cmdargs


def applyZonalstats(image,nodata, band, shape, uid):
        
    """
    function to derive stats for a single or multi band raster image
    """    
    # create an empty lists to write the results 
            
    zonestats = []
    siteID = []
    classifcation =[]
    image_Name = []
    nodata = nodata
    
    with rasterio.open(image, nodata=nodata) as srci:
        affine = srci.transform
        array = srci.read(band)
        
        with fiona.open(shape) as src:
            # added kwgs the interpolate to nearest - this ensures that only the value for the pixel under the point is returned
            # the 'bilinear' (which is the default) interploate method returns the weighted mean of the 4 closest pixels.  
            zs = point_query(src, array, interpolate='nearest', affine=affine,nodata=nodata) 
                        
            # extract the image name from the opened file from the input file read in by rasterio
            imgName1 = str(srci)[:-11]
            imgName = imgName1[-34:] 
            imgDate = imgName[16:20]      
            
            for zone in zs:
                value = zone
                           
                # put the individual results in a list and append them to the zonestats list
                result = [value]
                zonestats.append(result)
                            
            # extract out the site number for the polygon
            for i in src:
                table_attributes = i['properties'] # reads in the attribute table for each record 
                        
                site = table_attributes[uid] # reads in the id field from the attribute table and prints out the selected record 
                details = [site] #,imgDate]
                siteID.append(details)
                
                classif = table_attributes["class"] 
                Classif = [classif]
                classifcation.append(Classif)

        # join the elements in each of the lists row by row 
        finalresults =  [siteid + classi + zoneR for siteid, classi, zoneR in zip(siteID,classifcation, zonestats)] #,image_Name)]                    
        
        # close the vector and raster file 
        src.close() 
        srci.close() 

        # print out the file name of the processed image
        print (imgName1 + ' ' + 'band' + ' ' + str(band) + ' ' + 'is' + ' ' + 'complete') 
                
    return(finalresults)

def mainRoutine():
        
    # read in the command arguments
    cmdargs = getCmdargs()
    image = cmdargs.image
    nodata= int(cmdargs.nodata)
    shape = cmdargs.shape 
    uid = cmdargs.uid
    export_csv = cmdargs.csv
    
    # make a temp dir to save the individual band results 
    tempDir = './temp_individual_bands'
    
    # check if the temp dir exists and if it does remove it, otherwise create the new dir for the individual outputs
    check_if_dir_exists = os.path.isdir(tempDir)
    if check_if_dir_exists == True:
        shutil.rmtree(tempDir)
    else:    
        os.makedirs(tempDir)
    
    # create a list of headers for the final output csv file
    final_header = []
        
    with rasterio.open(image, nodata=nodata) as srci:
        
        bands = srci.indexes # this will return the number of spectral band for the input raster image as a tuple
        num_bands = len(bands)
    
    for band in bands:
        
        # creates the individual band csv file name
        bandResults = 'band_'+str(band)+'.csv'

        # run the zonal stats function 
        finalresults = applyZonalstats(image, nodata, band,shape, uid)
        
        # write out the individual band results to a list
        
        outputlist = []
        
        for i in finalresults:
        
            outputlist.append(i)    
       
        # convert the list to a pandas dataframe with a headers identifying the band number being processed
        headers = ['site', "class", "b" + str(band)] 
        
        # select out the band being assessed to make the header for the final output csv
        current_band = "b" + str(band)
        final_header.append(current_band)
                                  
        output  = pd.DataFrame.from_records(outputlist,columns=headers)
              
        output.to_csv(tempDir + '/'+ bandResults,index=False)
    
    # read in the individual band results and concatenate them to a single dataframe
    all_files = glob.glob(os.path.join(tempDir, "*.csv"))     # advisable to use os.path.join as this makes concatenation OS independent
    
    df_from_each_file = (pd.read_csv(f) for f in all_files)
    concatenated_df   = pd.concat(df_from_each_file,ignore_index=False, axis=1)
    print (list(concatenated_df))
 
    # select out the site and class rows to add to the new df  
    new_df = concatenated_df.iloc[:, 0:2]
    
    # select out the bands only from the results to add to the final csv file
    new_df2 = concatenated_df[final_header]
    
    # join the selected data frames
    join_df = new_df.join(new_df2)
    
    # export the results to a csv file
    join_df.to_csv(export_csv)  
    
    # remove the temp dir and single band csv files
    shutil.rmtree(tempDir)
    
if __name__ == "__main__":
    mainRoutine()   