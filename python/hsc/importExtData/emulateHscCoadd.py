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

import lsst.pex.config as pexConfig
from lsst.pipe.tasks.coaddBase import CoaddBaseTask
from lsst.pipe.tasks.calibrate import CalibrateTask

import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.afw.table as afwTable


__all__ = ["EmulateHscCoaddTask"]

class EmulateHscCoaddConfig(CoaddBaseTask.ConfigClass):
    """Config for EmulateHscCoaddTask
    """

    fileOutName = pexConfig.Field("Name of output file", str, "exposure.fits")

    calibrate = pexConfig.ConfigurableField(
        target = CalibrateTask,
        doc = "Calibration (inc. high-threshold detection and measurement)",
    )

    detection  = pexConfig.ConfigurableField(
        target = measAlg.SourceDetectionTask,
        doc = "Initial (high-threshold) detection phase for calibration",
    )

    initialMeasurement = pexConfig.ConfigurableField(
        target = measAlg.SourceMeasurementTask,
        doc = "Initial measurements used to feed PSF determination and aperture correction determination",
    )

    
    # can we add filter info here, Without changing the policy files?
    #    filterPolicy = pexPolicy.Policy()
    #    filterPolicy.add("lambdaEff", 470.0)
    #    afwImage.Filter.define(afwImage.FilterProperty("g", filterPolicy))
    # see http://hsca.ipmu.jp/doxygen/3.6.1/page_p_a_f.html

class EmulateHscCoaddTask(CoaddBaseTask):


    ConfigClass  = EmulateHscCoaddConfig
    _DefaultName = "emulateHscCoadd"

    def __init__(self, *args, **kwargs):
        
        CoaddBaseTask.__init__(self, *args, **kwargs)

        schema = afwTable.SourceTable.makeMinimalSchema()
        self.schema = schema
        self.algMetadata = dafBase.PropertyList()
      
        self.makeSubtask("detection", schema=self.schema)
        self.makeSubtask("initialMeasurement", schema=self.schema, algMetadata=self.algMetadata)
        self.makeSubtask("calibrate")
        
      
       
    def run(self, patchRef, selectDataList=[]):

        # Import existing coadd 
        coadd    = patchRef.get(self.config.coaddName + "Coadd")

        # Save image, mask and variance as if imported from external data
        if False:
            coadd.getMaskedImage().getImage().writeFits(   '/Users/coupon/data/tmp/test_img.fits')
            coadd.getMaskedImage().getVariance().writeFits('/Users/coupon/data/tmp/test_var.fits')
            coadd.getMaskedImage().getMask().writeFits(    '/Users/coupon/data/tmp/test_msk.fits')
      
        # Create a new exposure from exported data
        # to do: check if input wcs and reference coadd have same wcs
        skyInfo  = self.getSkyInfo(patchRef)
        exposure = afwImage.ExposureF(skyInfo.bbox, skyInfo.wcs)

        maskedImage = afwImage.MaskedImageF(
            afwImage.ImageF("/Users/coupon/data/tmp/test_img.fits"), 
            afwImage.MaskU( "/Users/coupon/data/tmp/test_msk.fits"), 
            afwImage.ImageF("/Users/coupon/data/tmp/test_var.fits"))
        exposure.setMaskedImage(maskedImage)

        # Until we add CFHT filter info...
        exposure.setFilter(coadd.getFilter())
       

        #print dir(self.calibrate.measurePsf.psfDeterminer)
        #return

        
        # Start PSF measurement on coadd
        self.calibrate.installInitialPsf(exposure)
        idFactory = afwTable.IdFactory.makeSimple()     
        table    = afwTable.SourceTable.make(self.calibrate.schema, idFactory)
        table.setMetadata(self.calibrate.algMetadata)
        
        # detect sources
        detRet = self.detection.makeSourceCatalog(table, exposure)
        sources = detRet.sources
        
        self.initialMeasurement.measure(exposure, sources)

        for i, s in enumerate(sources):
            print s.getIxx(), s.getIyy()
            if i == 10: break
    



        return

        psfRet  = self.calibrate.measurePsf.run(exposure, sources, expId=0, matches=None)
        
        cellSet = psfRet.cellSet
        psf = psfRet.psf

        # Write exposure
        if False:
            exposure.writeFits(self.config.fileOutName)


        #exposure = afwImage.ExposureF("/Users/coupon/data/HSC/SSP/rerun/tutorial_3.6.1/deepCoadd/HSC-I/1/5,5.fits")
        #exposure = afwImage.ExposureF("/Users/coupon/data/HSC/SSP/rerun/tutorial_3.6.1/01116/HSC-I/corr/CORR-0019666-049.fits")

        #import lsst.afw.geom as afwGeom
        #exposure.setXY0(afwGeom.Point2D(0,0))

        #calib = self.calibrate.run(exposure)

        return



    # Overload these if your task inherits from CmdLineTask
    def _getConfigName(self):
        return None
    def _getEupsVersionsName(self):
        return None
    def _getMetadataName(self):
        return None
       





    
























