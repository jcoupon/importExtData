#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011, 2012 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import pyfits as fits


import lsst.pex.config as pexConfig
import lsst.afw.image as afwImage
from   lsst.pipe.tasks.makeCoaddTempExp import MakeCoaddTempExpTask

from lsst.pipe.tasks.calibrate import CalibrateTask

import lsst.daf.base as dafBase
import lsst.afw.table as afwTable
import lsst.meas.algorithms as measAlg


class EmulateHscCoaddConfig(MakeCoaddTempExpTask.ConfigClass):
    fileOutName = pexConfig.Field("Name of output file", str, "exposure.fits")
    
    calibrate = pexConfig.ConfigurableField(
        target = CalibrateTask,
        doc = "Calibration (inc. high-threshold detection and measurement)",
    )

    detection    = pexConfig.ConfigurableField(
        target = measAlg.SourceDetectionTask,
        doc = "Initial (high-threshold) detection phase for calibration",
    )
    initialMeasurement = pexConfig.ConfigurableField(
        target = measAlg.SourceMeasurementTask,
        doc = "Initial measurements used to feed PSF determination and aperture correction determination",
    )
 
    #    filterPolicy = pexPolicy.Policy()
    #    filterPolicy.add("lambdaEff", 470.0)
    #    afwImage.Filter.define(afwImage.FilterProperty("g", filterPolicy))
    # see http://hsca.ipmu.jp/doxygen/3.6.1/page_p_a_f.html
    
class EmulateHscCoaddTask(MakeCoaddTempExpTask):


    #_DefaultName = "EmulateHscCoadd"
    ConfigClass  = EmulateHscCoaddConfig

    def __init__(self, **kwargs):
        
        MakeCoaddTempExpTask.__init__(self, **kwargs)
        self.makeSubtask("calibrate")


        self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.algMetadata = dafBase.PropertyList()
        
        #self.schema = self.calibrate.schema
        #self.algMetadata = self.calibrate.algMetadata

        self.makeSubtask("detection", schema=self.schema)
        self.makeSubtask("initialMeasurement", schema=self.schema, algMetadata=self.algMetadata)


    def run(self, patchRef, selectDataList=[]):


        exposure = patchRef.get('raw')

        idFactory = afwTable.IdFactory.makeSimple()
        calib = self.calibrate.run(exposure, idFactory=idFactory, expId=1)
          



        exit(-1)











        # test
        coadd    = patchRef.get(self.config.coaddName + "Coadd")




        if False:
            coadd.getMaskedImage().getImage().writeFits(   '/Users/coupon/data/tmp/test_img.fits')
            coadd.getMaskedImage().getVariance().writeFits('/Users/coupon/data/tmp/test_var.fits')
            coadd.getMaskedImage().getMask().writeFits(    '/Users/coupon/data/tmp/test_msk.fits')
      
        skyInfo  = self.getSkyInfo(patchRef)
        exposure = afwImage.ExposureF(skyInfo.bbox, skyInfo.wcs)

        # to do: check if input wcs and reference coadd have same wcs
        maskedImage = afwImage.MaskedImageF(
            afwImage.ImageF("/Users/coupon/data/tmp/test_img.fits"), 
            afwImage.MaskU( "/Users/coupon/data/tmp/test_msk.fits"), 
            afwImage.ImageF("/Users/coupon/data/tmp/test_var.fits"))
        exposure.setMaskedImage(maskedImage)

        #exposure.writeFits(self.config.fileOutName)


        exposure.setFilter(coadd.getFilter())

       
        self.calibrate.installInitialPsf(exposure)

        idFactory = afwTable.IdFactory.makeSimple()
     
        table = afwTable.SourceTable.make(self.schema, idFactory)
        table.setMetadata(self.algMetadata)
        detRet = self.detection.makeSourceCatalog(table, exposure)
        sources = detRet.sources




     
        self.initialMeasurement.measure(exposure, sources)

        matches = None


        for s in sources:
            print s.getIxx(), s.getIyy()
     




        psfRet = self.calibrate.measurePsf.run(exposure, sources, expId=0, matches=matches)
        
        cellSet = psfRet.cellSet
        psf = psfRet.psf


        return result





       





    
























