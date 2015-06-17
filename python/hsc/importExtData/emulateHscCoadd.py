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


import math

import lsst.pex.config as pexConfig
from lsst.pipe.tasks.coaddBase import CoaddBaseTask
from lsst.pipe.tasks.calibrate import CalibrateTask
from lsst.pipe.tasks.measurePsf import MeasurePsfTask

import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.afw.table as afwTable


__all__ = ["EmulateHscCoaddTask"]



class InitialPsfConfig(pexConfig.Config):
    """Describes the initial PSF used for detection and measurement before we do PSF determination."""

    model = pexConfig.ChoiceField(
        dtype = str,
        doc = "PSF model type",
        default = "SingleGaussian",
        allowed = {
            "SingleGaussian": "Single Gaussian model",
            "DoubleGaussian": "Double Gaussian model",
        },
    )
    fwhm = pexConfig.Field(
        dtype = float,
        doc = "FWHM of PSF model (arcsec)",
        default = 1.0,
    )
    size = pexConfig.Field(
        dtype = int,
        doc = "Size of PSF model (pixels)",
        default = 15,
    )


class EmulateHscCoaddConfig(CoaddBaseTask.ConfigClass):
    """Config for EmulateHscCoaddTask
    """


    imgInName = pexConfig.Field("Name of input image",          str, "im.fits")
    mskInName = pexConfig.Field("Name of input mask",           str, "msk.fits")
    varInName = pexConfig.Field("Name of input variance image", str, "var.fits")

    fileOutName = pexConfig.Field("Name of output file", str, "exposure.fits")

    magLim = pexConfig.Field("Magnitude faint limit for PSF measurement", float, 23.0)

    initialPsf = pexConfig.ConfigField(dtype=InitialPsfConfig, doc=InitialPsfConfig.__doc__)

    detection  = pexConfig.ConfigurableField(
        target = measAlg.SourceDetectionTask,
        doc = "Initial (high-threshold) detection phase for calibration",
    )

    initialMeasurement = pexConfig.ConfigurableField(
        target = measAlg.SourceMeasurementTask,
        doc = "Initial measurements used to feed PSF determination and aperture correction determination",
    )

    measurePsf   = pexConfig.ConfigurableField(target = MeasurePsfTask, doc = "")


    def setDefaults(self):

        pexConfig.Config.setDefaults(self)

        self.detection.includeThresholdMultiplier = 10.0
        self.initialMeasurement.prefix = "initial."
        self.initialMeasurement.algorithms.names -= ["correctfluxes"]
        initflags = [self.initialMeasurement.prefix+x
                     for x in self.measurePsf.starSelector["catalog"].badStarPixelFlags]
        self.measurePsf.starSelector["catalog"].badStarPixelFlags.extend(initflags)
        
        # Crashes if > 0.0
        self.measurePsf.reserveFraction = 0.0

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

        self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.algMetadata = dafBase.PropertyList()
      
        self.makeSubtask("detection", schema=self.schema)
        self.makeSubtask("initialMeasurement", schema=self.schema, algMetadata=self.algMetadata)
        #self.makeSubtask("calibrate")
        self.makeSubtask("measurePsf")
        
      
    def run(self, patchRef, selectDataList=[]):

        # ---------------------------------------------- #
        # Import reference coadd and records wcs info
        # ---------------------------------------------- #

        coadd    = patchRef.get(self.config.coaddName + "Coadd")
        skyInfo  = self.getSkyInfo(patchRef)

        # ---------------------------------------------- #
        # tests
        # ---------------------------------------------- #

        if True:
    
            fluxMag0 = coadd.getCalib().getFluxMag0()[0]
    
            starSelectorConfig = self.measurePsf.starSelector.ConfigClass()
            starSelectorConfig.fluxLim = fluxMag0 * pow(10.0, -0.4*self.config.magLim)
            self.measurePsf.starSelector.config = starSelectorConfig

            imgInName, mskInName, varInName = self.writeTest(coadd, dirName="/Users/coupon/data/tmp")

        # ---------------------------------------------- #
        # Create new exposure and feed with input images
        # ---------------------------------------------- #

        exposure = afwImage.ExposureF(skyInfo.bbox, skyInfo.wcs)
        maskedImage = afwImage.MaskedImageF(
            afwImage.ImageF(imgInName), 
            afwImage.MaskU( mskInName), 
            afwImage.ImageF(varInName))
        exposure.setMaskedImage(maskedImage)

        # Until we add CFHT filter info...
        exposure.setFilter(coadd.getFilter())

        #exposure = afwImage.ExposureF("/Users/coupon/data/HSC/SSP/rerun/tutorial_3.6.1/deepCoadd/HSC-I/1/5,5.fits")
        #exposure = afwImage.ExposureF("/Users/coupon/data/HSC/SSP/rerun/tutorial_3.6.1/01116/HSC-I/corr/CORR-0019666-049.fits")
        
        # ---------------------------------------------- #
        # Start PSF measurement on coadd
        # ---------------------------------------------- #

        self.installInitialPsf(exposure)

        # prepare table        
        idFactory = afwTable.IdFactory.makeSimple()     
        table    = afwTable.SourceTable.make(self.schema, idFactory)
        table.setMetadata(self.algMetadata)
        # detect sources
        detRet = self.detection.makeSourceCatalog(table, exposure)
        sources = detRet.sources
        # measure moments
        self.initialMeasurement.measure(exposure, sources)
        # measure psf
        psfRet  = self.measurePsf.run(exposure, sources, expId=0, matches=None)
        
        cellSet = psfRet.cellSet
        psf = psfRet.psf

        exposure.setPsf(psf)


        # Write exposure
        if False:
            exposure.writeFits(self.config.fileOutName)

        return exposure


    def writeTest(self, coadd, dirName="."):
        """ This method takes the input coadd 
        and writes the images, mask and variance 
        independently in dirName
        """

        imgInName = dirName+"/"+"img.fits"
        mskInName = dirName+"/"+"msk.fits"
        varInName = dirName+"/"+"var.fits"

        # Save image, mask and variance as if imported from external data
        coadd.getMaskedImage().getImage().writeFits(   imgInName)
        coadd.getMaskedImage().getMask().writeFits(    mskInName)
        coadd.getMaskedImage().getVariance().writeFits(varInName)
     
        return imgInName, mskInName, varInName


    def installInitialPsf(self, exposure):
        """[Method taken from lsst.pipe.tasks.calibrate]

        Initialise the calibration procedure by setting the PSF to a configuration-defined guess.

        @param[in,out] exposure Exposure to process; fake PSF will be installed here.
        """
        assert exposure, "No exposure provided"
        
        wcs = exposure.getWcs()
        assert wcs, "No wcs in exposure"

        cls = getattr(measAlg, self.config.initialPsf.model + "Psf")

        fwhm = self.config.initialPsf.fwhm / wcs.pixelScale().asArcseconds()
        size = self.config.initialPsf.size
        self.log.info("installInitialPsf fwhm=%.2f pixels; size=%d pixels" % (fwhm, size))
        psf = cls(size, size, fwhm/(2*math.sqrt(2*math.log(2))))
        exposure.setPsf(psf)


    # Overload these if your task inherits from CmdLineTask
    def _getConfigName(self):
        return None
    def _getEupsVersionsName(self):
        return None
    def _getMetadataName(self):
        return None























