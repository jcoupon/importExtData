
# coding: utf-8

# In[1]:


import pickle
import numpy as np

import lsst.daf.persistence as dafPersist
import lsst.afw.cameraGeom as camGeom
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage

from astropy.wcs import WCS
from lsst.afw import geom


# In[2]:


skyMap = pickle.load( open('/data0/desprezg/HSC/rerun/ssp_map/deepCoadd/skyMap.pickle', "rb" ) )


# In[3]:


def create_header(tract_index,patch_index1, patch_index2):
    """
    Create the header for a patch's tract
    """
    #tract = skyMap.generateTract(tract_index)
    tract = skyMap[tract_index]
    patch = tract[patch_index1, patch_index2]
        
    tract_wcs = tract.getWcs()
    #tract_fits_meta = tract.getFitsMetadata()
    
    hdr_wcs = WCS(naxis=2)
    hdr = hdr_wcs.to_header()
    
    # Get reference Sky Coordinates
    #hdr['PATCH'] = str(tract_index)+'-'+str(patch_index1)+'.'+str(patch_index2)
    #hdr['CRVAL1'] = tract_wcs.getSkyOrigin()[0].asDegrees()
    #hdr['CRVAL2'] = tract_wcs.getSkyOrigin()[1].asDegrees()
    
    # Get reference Pixel Coordinates
    #hdr['CRPIX1'] = tract_wcs.getPixelOrigin()[0] - patch.getOuterBBox().getMinX() + 1
    #hdr['CRPIX2'] = tract_wcs.getPixelOrigin()[1] - patch.getOuterBBox().getMinY() + 1
    
    # Get the WCS coordinate scale matrix
    hdr['CD1_1'] = tract_wcs.getCdMatrix()[0,0]
    hdr['CD1_2'] = tract_wcs.getCdMatrix()[0,1]
    hdr['CD2_1'] = tract_wcs.getCdMatrix()[1,0]
    hdr['CD2_2'] = tract_wcs.getCdMatrix()[1,1]
    
    
    # Get Dimension in pixels
    hdr['NAXIS1'] = patch.getOuterBBox().getDimensions()[0]
    hdr['NAXIS2'] = patch.getOuterBBox().getDimensions()[1]
    
    # Get reference Sky Coordinates
    
    xy0 = patch.getOuterBBox().getMin()
    shift = geom.Extent2D(geom.Point2I(0, 0) - xy0)
    newWcs = tract.getWcs().copyAtShiftedPixelOrigin(shift)
    wcsMetadata = newWcs.getFitsMetadata(precise=True)
    
    wcsMetadata.combine(geom.wcsUtils.createTrivialWcsMetadata("A", xy0))
    
    hdr.update(wcsMetadata.toDict())
    
    # Get patch pixel coordinates in Tract reference
    # Get LTV1 and 2 (minus of CRVAL1A and 2A )
    hdr['LTV1'] = - hdr['CRVAL1A']
    hdr['LTV2'] = - hdr['CRVAL2A']
    
    return hdr


# In[4]:


def write_header(path,hdr,tract_index,patch_index1,patch_index2):
    """
    """
    PATCH = str(tract_index)+'-'+str(patch_index1)+'.'+str(patch_index2)
    print('[INFO] Creating '+PATCH+'.head file')
    hdr.totextfile(path+PATCH+'.head',overwrite = True)
   
    return
    


# In[7]:


if __name__ == '__main__':
    path = '/data0/desprezg/SSP_header/'
    tract_list = [9569,9570,9571,9572,9812,9813,9814,10054,10055,10056] #cosmos deep tract list
    
    #tract_list = [8523,8524,8525,8765,8766,8767,8282,8283,8284] # XMM-LSS deep tract list
    
    for tract_i in tract_list:
        for p_1 in np.arange(9):
            for p_2 in np.arange(9):
                
                hdr = create_header(tract_i,p_1,p_2)
                write_header(path,hdr,tract_i,p_1,p_2)

