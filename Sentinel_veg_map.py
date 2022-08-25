#!/usr/bin/env python

from __future__ import print_function, division
import sys
import os
import argparse
import pickle as pickle
import numpy as np
from rios import applier, fileinfo
import pdb
from sklearn.preprocessing import Imputer

def getCmdargs():
    """
    Get command line arguments
    """
    p = argparse.ArgumentParser()
    p.add_argument("--reffile_a", help="Input the annual surface reflectance file")
    
    p.add_argument("--reffile_d", help="Input the dry season surface reflectance file")
    
    p.add_argument("--DEM", help="Input the DEM file")
    
    p.add_argument("--stc", help="Input the structural formation image")
    
    p.add_argument("--composite", help="Input Tree Structure Composite")
    
    p.add_argument("--outfile", help="Name of output file (default is chm of reffile)")
    
    p.add_argument("--picklefile", help="Input pickle file")
    
    #This is the original line: #p.add_argument("--picklefile", default="U:\Working\Land_veg\AA_Vegetation\INPEX_project\PP_B\Data\Vector\Training\Zonal_Stats_Tables\Table_V4_Landsat_Result4.p",help="Input pickle file (default is %(default)s)")
    
    cmdargs = p.parse_args()
    
    if cmdargs.reffile_a is None:
        p.print_help()
        sys.exit()
        
    return cmdargs


def main():
    """
    Main routine
    
    """
    cmdargs = getCmdargs()
    controls = applier.ApplierControls()
    infiles = applier.FilenameAssociations()
    outfiles = applier.FilenameAssociations()
    otherargs = applier.OtherInputs()
    
    infiles.sfcref_a = cmdargs.reffile_a
    
    infiles.sfcref_d = cmdargs.reffile_d
    
    infiles.DEM = cmdargs.DEM
    
    infiles.stc = cmdargs.stc
    
    infiles.composite = cmdargs.composite
    
    
       
    imageExt = infiles.sfcref_d
         
    controls.setReferenceImage(infiles.sfcref_a) 
    
    outfiles.hgt = cmdargs.outfile
    
    otherargs.rf = pickle.load(open(cmdargs.picklefile, 'rb'))
    # refInfo = fileinfo.ImageInfo(infiles.sfcref)

    otherargs.refnull =  32767

    applier.apply(doModel, infiles, outfiles, otherargs,controls=controls)
    
    
    
def doModel(info, inputs, outputs, otherargs):

    nonNullmask = (inputs.sfcref_a[0] != otherargs.refnull)
    
    # get the shape of the annual image and convert it to the shape of a single band  
    a_imgshape = inputs.sfcref_a.shape
    # convert the tuple to a list to convert the 6 bands to represent 1 band and then convert it back to a tuple
    list_imgshape = list(a_imgshape)
    list_imgshape[0] = 1
    imgShape = tuple(list_imgshape)

    # input the predictor variables 
    #psB1a = inputs.sfcref_a[0][nonNullmask]
    psB2a = inputs.sfcref_a[1][nonNullmask]
    psB3a = inputs.sfcref_a[2][nonNullmask]
    psB4a = inputs.sfcref_a[3][nonNullmask]
    psB5a = inputs.sfcref_a[4][nonNullmask]
    psB6a = inputs.sfcref_a[5][nonNullmask]
    psB7a = inputs.sfcref_a[6][nonNullmask]
    psB8a = inputs.sfcref_a[7][nonNullmask]
    psB9a = inputs.sfcref_a[8][nonNullmask]
    psB10a = inputs.sfcref_a[9][nonNullmask]
    
    #psB1d = inputs.sfcref_d[0][nonNullmask]
    psB2d = inputs.sfcref_d[1][nonNullmask]
    psB3d = inputs.sfcref_d[2][nonNullmask]
    psB4d = inputs.sfcref_d[3][nonNullmask]
    psB5d = inputs.sfcref_d[4][nonNullmask]
    psB6d = inputs.sfcref_d[5][nonNullmask]
    psB7d = inputs.sfcref_d[6][nonNullmask]
    psB8d = inputs.sfcref_d[7][nonNullmask]
    psB9d = inputs.sfcref_d[8][nonNullmask]
    psB10d = inputs.sfcref_d[9][nonNullmask]
    
       
    # convert the reflectance data to floating point and calculate the band ratios and veg indicies calculations   
    
    #psB1fa = (psB1a*0.0001)+0.0 # blue
    psB2fa = (psB2a*0.0001)+0.0 # green
    psB3fa = (psB3a*0.0001)+0.0 # red
    psB4fa = (psB4a*0.0001)+0.0 # NIR
    psB5fa = (psB5a*0.0001)+0.0 # rededge1
    psB6fa = (psB6a*0.0001)+0.0 # rededge2
    psB7fa = (psB7a*0.0001)+0.0 # rededge3
    psB8fa = (psB8a*0.0001)+0.0 # narrownir
    psB9fa = (psB9a*0.0001)+0.0 # swir1
    psB10fa = (psB10a*0.0001)+0.0 # swir2
    
    
    #psB1fd = (psB1d*0.0001)+0.0 # blue
    psB2fd = (psB2d*0.0001)+0.0 # green
    psB3fd = (psB3d*0.0001)+0.0 # red
    psB4fd = (psB4d*0.0001)+0.0 # NIR
    psB5fd = (psB5d*0.0001)+0.0 # rededge1
    psB6fd = (psB6d*0.0001)+0.0 # rededge2
    psB7fd = (psB7d*0.0001)+0.0 # rededge3
    psB8fd = (psB8d*0.0001)+0.0 # narrownir
    psB9fd = (psB9d*0.0001)+0.0 # swir1
    psB10fd = (psB10d*0.0001)+0.0 # swir2
    
    #pdb.set_trace()
   
    # ratio calculations   
    
    ratio32a = np.int32(np.around((psB3a / psB2a)*10**7))
    ratio42a = np.int32(np.around((psB4a / psB2a)*10**7))
    ratio43a = np.int32(np.around((psB4a / psB3a)*10**7))
    ratio52a = np.int32(np.around((psB5a / psB2a)*10**7)) 
    ratio53a = np.int32(np.around((psB5a / psB3a)*10**7))
    ratio54a = np.int32(np.around((psB5a / psB4a)*10**7)) 
    ratio62a = np.int32(np.around((psB6a / psB2a)*10**7))
    ratio63a = np.int32(np.around((psB6a / psB3a)*10**7)) 
    ratio64a = np.int32(np.around((psB6a / psB4a)*10**7))
    ratio65a = np.int32(np.around((psB6a / psB5a)*10**7))
    ratio72a = np.int32(np.around((psB7a / psB2a)*10**7))
    ratio73a = np.int32(np.around((psB7a / psB3a)*10**7)) 
    ratio74a = np.int32(np.around((psB7a / psB4a)*10**7))
    ratio75a = np.int32(np.around((psB7a / psB5a)*10**7))
    ratio76a = np.int32(np.around((psB7a / psB6a)*10**7))
    ratio82a = np.int32(np.around((psB8a / psB2a)*10**7))
    ratio83a = np.int32(np.around((psB8a / psB3a)*10**7)) 
    ratio84a = np.int32(np.around((psB8a / psB4a)*10**7))
    ratio85a = np.int32(np.around((psB8a / psB5a)*10**7))
    ratio86a = np.int32(np.around((psB8a / psB6a)*10**7))
    ratio87a = np.int32(np.around((psB8a / psB7a)*10**7))
    ratio92a = np.int32(np.around((psB9a / psB2a)*10**7))
    ratio93a = np.int32(np.around((psB9a / psB3a)*10**7)) 
    ratio94a = np.int32(np.around((psB9a / psB4a)*10**7))
    ratio95a = np.int32(np.around((psB9a / psB5a)*10**7))
    ratio96a = np.int32(np.around((psB9a / psB6a)*10**7))
    ratio97a = np.int32(np.around((psB9a / psB7a)*10**7))
    ratio98a = np.int32(np.around((psB9a / psB8a)*10**7))
    ratio102a = np.int32(np.around((psB10a / psB2a)*10**7))
    ratio103a = np.int32(np.around((psB10a / psB3a)*10**7)) 
    ratio104a = np.int32(np.around((psB10a / psB4a)*10**7))
    ratio105a = np.int32(np.around((psB10a / psB5a)*10**7))
    ratio106a = np.int32(np.around((psB10a / psB6a)*10**7))
    ratio107a = np.int32(np.around((psB10a / psB7a)*10**7))
    ratio108a = np.int32(np.around((psB10a / psB8a)*10**7))
    ratio109a = np.int32(np.around((psB10a / psB9a)*10**7))
    
    ratio32d = np.int32(np.around((psB3d / psB2d)*10**7))
    ratio42d = np.int32(np.around((psB4d / psB2d)*10**7))
    ratio43d = np.int32(np.around((psB4d / psB3d)*10**7))
    ratio52d = np.int32(np.around((psB5d / psB2d)*10**7)) 
    ratio53d = np.int32(np.around((psB5d / psB3d)*10**7))
    ratio54d = np.int32(np.around((psB5d / psB4d)*10**7)) 
    ratio62d = np.int32(np.around((psB6d / psB2d)*10**7))
    ratio63d = np.int32(np.around((psB6d / psB3d)*10**7)) 
    ratio64d = np.int32(np.around((psB6d / psB4d)*10**7))
    ratio65d = np.int32(np.around((psB6d / psB5d)*10**7))
    ratio72d = np.int32(np.around((psB7d / psB2d)*10**7))
    ratio73d = np.int32(np.around((psB7d / psB3d)*10**7)) 
    ratio74d = np.int32(np.around((psB7d / psB4d)*10**7))
    ratio75d = np.int32(np.around((psB7d / psB5d)*10**7))
    ratio76d = np.int32(np.around((psB7d / psB6d)*10**7))
    ratio82d = np.int32(np.around((psB8d / psB2d)*10**7))
    ratio83d = np.int32(np.around((psB8d / psB3d)*10**7)) 
    ratio84d = np.int32(np.around((psB8d / psB4d)*10**7))
    ratio85d = np.int32(np.around((psB8d / psB5d)*10**7))
    ratio86d = np.int32(np.around((psB8d / psB6d)*10**7))
    ratio87d = np.int32(np.around((psB8d / psB7d)*10**7))
    ratio92d = np.int32(np.around((psB9d / psB2d)*10**7))
    ratio93d = np.int32(np.around((psB9d / psB3d)*10**7)) 
    ratio94d = np.int32(np.around((psB9d / psB4d)*10**7))
    ratio95d = np.int32(np.around((psB9d / psB5d)*10**7))
    ratio96d = np.int32(np.around((psB9d / psB6d)*10**7))
    ratio97d = np.int32(np.around((psB9d / psB7d)*10**7))
    ratio98d = np.int32(np.around((psB9d / psB8d)*10**7))
    ratio102d = np.int32(np.around((psB10d / psB2d)*10**7))
    ratio103d = np.int32(np.around((psB10d / psB3d)*10**7)) 
    ratio104d = np.int32(np.around((psB10d / psB4d)*10**7))
    ratio105d = np.int32(np.around((psB10d / psB5d)*10**7))
    ratio106d = np.int32(np.around((psB10d / psB6d)*10**7))
    ratio107d = np.int32(np.around((psB10d / psB7d)*10**7))
    ratio108d = np.int32(np.around((psB10d / psB8d)*10**7))
    ratio109d = np.int32(np.around((psB10d / psB9d)*10**7))
    
    
    GSAVIa = np.int32(np.around((((psB4fa-psB2fa)/(psB4fa+psB2fa+0.5))*(1.5))*10**7))
    GNDVIa = np.int32(np.around(((psB4fa-psB2fa)/(psB4fa+psB2fa))*10**7))
    CVIa = np.int32(np.around(((psB4fa/psB2fa)*(psB3fa/psB2fa))*10**7))
    NDGIa = np.int32(np.around(((psB2fa-psB3fa)/(psB2fa+psB3fa))*10**7))   
    RIa = np.int32(np.around(((psB3fa-psB2fa)/(psB3fa+psB2fa))*10**7))
    NBRa = np.int32(np.around(((psB4fa-psB10fa)/(psB4fa+psB10fa))*10**7))
    NDIIa =  np.int32(np.around(((psB4fa-psB9fa)/(psB4fa+psB9fa))*10**7))
    GDVIa = np.int32(np.around((psB4fa-psB2fa)*10**7))
    DVIa = np.int32(np.around((psB4fa-psB3fa)*10**7))
    SAVIa = np.int32(np.around((((psB4fa-psB3fa)/(psB4fa+psB3fa+0.5))*(1.5))*10**7))
    NDVIa = np.int32(np.around(((psB4fa-psB3fa)/(psB4fa+psB3fa))*10**7))
    MSRa = np.int32(np.around(((psB4fa/psB3fa)-1)/(((np.sqrt(psB4fa/psB3fa))+1))*10**7))
    MSAVIa = np.int32(np.around(((((2 * psB4fa) + 1) - np.sqrt((np.power(((2 * psB4fa) +1),2)) - (8 * (psB4fa - psB3fa))))/2)*10**7))
    CIgreena = np.int32(np.around(((psB7fa/psB2fa)-1)*10**7))
    NDRE1a = np.int32(np.around(((psB6fa-psB5fa)/(psB6fa+psB5fa))*10**7))
    NDRE2a = np.int32(np.around(((psB7fa-psB5fa)/(psB7fa+psB5fa))*10**7))
    CIrededgea = np.int32(np.around(((psB7fa/psB5fa)-1)*10**7))
    
    
    #NDSIa = (psB5fa-psB4fa)/(psB5fa+psB4fa).astype(np.float32)
       
    GSAVId = np.int32(np.around((((psB4fd-psB2fd)/(psB4fd+psB2fd+0.5))*(1.5))*10**7))
    GNDVId = np.int32(np.around(((psB4fd-psB2fd)/(psB4fd+psB2fd))*10**7))
    CVId = np.int32(np.around(((psB4fd/psB2fd)*(psB3fd/psB2fd))*10**7))
    NDGId = np.int32(np.around(((psB2fd-psB3fd)/(psB2fd+psB3fd))*10**7))   
    RId = np.int32(np.around(((psB3fd-psB2fd)/(psB3fd+psB2fd))*10**7))
    NBRd = np.int32(np.around(((psB4fd-psB10fd)/(psB4fd+psB10fd))*10**7))
    NDIId =  np.int32(np.around(((psB4fd-psB9fd)/(psB4fd+psB9fd))*10**7))
    GDVId = np.int32(np.around((psB4fd-psB2fd)*10**7))
    DVId = np.int32(np.around((psB4fd-psB3fd)*10**7))
    SAVId = np.int32(np.around((((psB4fd-psB3fd)/(psB4fd+psB3fd+0.5))*(1.5))*10**7))
    NDVId = np.int32(np.around(((psB4fd-psB3fd)/(psB4fd+psB3fd))*10**7))
    MSRd = np.int32(np.around(((psB4fd/psB3fd)-1)/(((np.sqrt(psB4fd/psB3fd))+1))*10**7))
    MSAVId = np.int32(np.around(((((2 * psB4fd) + 1) - np.sqrt((np.power(((2 * psB4fd) +1),2)) - (8 * (psB4fd - psB3fd))))/2)*10**7))
    #DBSId = np.int32(np.around((((psB9fd-psB2fd)/(psB9fd+psB2fd))-((psB4fd-psB3fd)/(psB4fd+psB3fd)))*10**7))
    CIgreend = np.int32(np.around(((psB7fd/psB2fd)-1)*10**7))
    NDRE1d = np.int32(np.around(((psB6fd-psB5fd)/(psB6fd+psB5fd))*10**7))
    NDRE2d = np.int32(np.around(((psB7fd-psB5fd)/(psB7fd+psB5fd))*10**7))
    CIrededged = np.int32(np.around(((psB7fd/psB5fd)-1)*10**7))
    
    #NDSId = (psB5fd-psB4fd)/(psB5fd+psB4fd).astype(np.float32)
   
    # read band composite file
    H99 = inputs.composite[0][nonNullmask]
    MCH = inputs.composite[1][nonNullmask]
    COV = inputs.composite[2][nonNullmask]
    
    
    
    DEM = inputs.DEM[0][nonNullmask]
    stc = inputs.stc[0][nonNullmask]
    
    # pass the variables into a np array and transform it to look like the pandas df 
    
    allVars= np.vstack([psB2a,psB3a,psB4a,psB5a,psB6a,psB7a,psB8a,psB9a,psB10a,ratio32a,ratio42a,ratio43a,ratio52a,ratio53a,ratio54a, ratio62a, ratio63a, ratio64a, ratio65a, ratio72a, ratio73a, ratio74a, ratio75a, ratio76a, ratio82a, ratio83a, ratio84a, ratio85a, ratio86a, ratio87a, ratio92a, ratio93a, ratio94a, ratio95a, ratio96a, ratio97a, ratio98a, ratio102a, ratio103a, ratio104a, ratio105a, ratio106a, ratio107a, ratio108a, ratio109a, GSAVIa, GNDVIa, CVIa, NDGIa, RIa,NBRa,NDIIa,GDVIa, DVIa,SAVIa, NDVIa, MSRa, MSAVIa, CIgreena, NDRE1a, NDRE2a, CIrededgea, psB2d, psB3d, psB4d, psB5d, psB6d, psB7d, psB8d, psB9d, psB10d, ratio32d, ratio42d, ratio43d, ratio52d, ratio53d, ratio54d, ratio62d, ratio63d, ratio64d, ratio65d, ratio72d, ratio73d, ratio74d, ratio75d, ratio76d, ratio82d, ratio83d, ratio84d, ratio85d, ratio86d, ratio87d, ratio92d, ratio93d, ratio94d, ratio95d, ratio96d, ratio97d, ratio98d, ratio102d, ratio103d, ratio104d, ratio105d, ratio106d, ratio107d, ratio108d, ratio109d, GSAVId, GNDVId, CVId, NDGId, RId, NBRd, NDIId, GDVId, DVId, SAVId, NDVId, MSRd, MSAVId, CIgreend, CIrededged, NDRE1d, NDRE2d, stc, DEM, H99, MCH, COV]).T
        

   
    # sets up the shape and dtype for the chm output  
    outputs.hgt = np.zeros(imgShape, dtype=np.uint8)
    
    # applies the rfr model to produce the chm layer
    
    if allVars.shape[0] > 0:
        # run check over the input data and replaces nan and infinity values
        allVars[np.isnan(allVars)] = 0.0
        allVars[np.isinf(allVars)] = 0.0
        
        hgt = otherargs.rf.predict(allVars)
        
        outputs.hgt[0][nonNullmask] = hgt
    

if __name__ == "__main__":
    main()
