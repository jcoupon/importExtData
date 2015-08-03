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

    mag0   = pexConfig.Field("Magnitude zero point", float, 27.0)

    magLim = pexConfig.Field("Magnitude faint limit for PSF measurement", float, 23.0)

    filtName = pexConfig.Field("Filter name", str, None)

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


        if False:
            """
            This is HSC config, but it seems that PSFex isn't working on coadds. deblend.nchild is also
            returning an error

            """
            import os
            self.initialMeasurement.load(os.path.join(os.environ['MEAS_EXTENSIONS_SHAPEHSM_DIR'], 'config', 'enable.py'))
            #self.initialMeasurement.algorithms["shape.hsm.regauss"].deblendNChild = "deblend.nchild"
            self.initialMeasurement.slots.shape = "shape.hsm.moments"

            try:
                import lsst.meas.extensions.psfex.psfexPsfDeterminer
                self.measurePsf.psfDeterminer["psfex"].spatialOrder = 2
                self.measurePsf.psfDeterminer.name = "psfex"
            except ImportError as e:
                print "WARNING: Unable to use psfex: %s" % e
                self.measurePsf.psfDeterminer.name = "pca"


#        initflags = [self.initialMeasurement.prefix+x
#                     for x in self.measurePsf.starSelector["catalog"].badStarPixelFlags]
#        self.measurePsf.starSelector["catalog"].badStarPixelFlags.extend(initflags)


        # Crashes if > 0.0
        # TODO: implement it
        self.measurePsf.reserveFraction = 0.0

class EmulateHscCoaddTask(CoaddBaseTask):

    ConfigClass  = EmulateHscCoaddConfig
    _DefaultName = "emulateHscCoadd"

    def __init__(self, *args, **kwargs):

        CoaddBaseTask.__init__(self, *args, **kwargs)

        self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.algMetadata = dafBase.PropertyList()

        self.makeSubtask("detection", schema=self.schema)
        self.makeSubtask("initialMeasurement", schema=self.schema, algMetadata=self.algMetadata)
        self.makeSubtask("measurePsf")

    def run(self, patchRef, selectDataList=[]):
        """
        Task to import external data and transform it into LSST exposure
        object. The task also perfomrs PSF measurement
        """

        # ---------------------------------------------- #
        # Import reference coadd and records wcs info
        # ---------------------------------------------- #

        coadd    = patchRef.get(self.config.coaddName + "Coadd")
        skyInfo  = self.getSkyInfo(patchRef)


        #measAlg.utils.showPsfMosaic(coadd, coadd.getPsf(), frame=1)
        #return

        #mask = afwImage.MaskU(skyInfo.bbox)
        #mask.writeFits("mask.fits")
        #return

        # ---------------------------------------------- #
        # Create new exposure object and feed with input images
        # ---------------------------------------------- #

        fluxMag0 = pow(10.0, +0.4*self.config.mag0)

        if False:
            imgInName, mskInName, varInName = self.writeTest(coadd, dirName="/Users/coupon/data/tmp")
            fluxMag0 = coadd.getCalib().getFluxMag0()[0]
        else:
            imgInName, mskInName, varInName = self.config.imgInName, self.config.mskInName, self.config.varInName


        # set everything but the image
        imgInName, _ , _ = self.writeTest(coadd, dirName="/Users/coupon/data/tmp")


        exposure = afwImage.ExposureF(skyInfo.bbox, skyInfo.wcs)
        maskedImage = afwImage.MaskedImageF(
            afwImage.ImageF(imgInName),
            afwImage.MaskU( mskInName),
            afwImage.ImageF(varInName))
        exposure.setMaskedImage(maskedImage)

        # set dummy coadd info
        expSchema = afwTable.ExposureTable.makeMinimalSchema()
        coaddInputs = afwImage.CoaddInputs(expSchema, expSchema)
        exposure.getInfo().setCoaddInputs(coaddInputs)

        # coaddInputs = coadd.getInfo().getCoaddInputs()
        # exposure.getInfo().setCoaddInputs(coaddInputs)

        # set filter
        if self.config.filtName is None:
            filt = coadd.getFilter()
        else:
            filt = afwImage.Filter(self.config.filtName)

        exposure.setFilter(filt)

        # set calib object
        calib = afwImage.Calib()
        calib.setFluxMag0(fluxMag0)
        exposure.setCalib(calib)

        # test coadd
        # exposure = afwImage.ExposureF("/Users/coupon/data/HSC/SSP/rerun/tutorial_3.6.1/deepCoadd/HSC-I/1/5,5.fits")

        # test single exposure
        # exposure = afwImage.ExposureF("/Users/coupon/data/HSC/SSP/rerun/tutorial_3.6.1/01116/HSC-I/corr/CORR-0019666-049.fits")


        # ---------------------------------------------- #
        # Do PSF measurement on coadd
        # ---------------------------------------------- #


        #starSelectorConfig = self.measurePsf.starSelector.ConfigClass()
        #starSelectorConfig.fluxLim = fluxMag0 * pow(10.0, -0.4*self.config.magLim)
        #self.measurePsf.starSelector.config = starSelectorConfig

        self.installInitialPsf(exposure)

        starSelectorName = "objectSize"
        starSelectorClass  = measAlg.starSelectorRegistry.get(starSelectorName)
        starSelectorConfig = starSelectorClass.ConfigClass()

        starSelectorConfig.sourceFluxField = "initial.flux.psf"

        starSelectorConfig.fluxMin = fluxMag0 * pow(10.0, -0.4*self.config.magLim)

        self.measurePsf.starSelector = starSelectorClass(starSelectorConfig)

        # prepare table
        idFactory = afwTable.IdFactory.makeSimple()
        table     = afwTable.SourceTable.make(self.schema, idFactory)
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

        #fwhm = self.config.initialPsf.fwhm / exposure.getWcs().pixelScale().asArcseconds()
        #psf  = measAlg.DoubleGaussianPsf(15, 15, fwhm/(2*math.sqrt(2*math.log(2))))

        display = True
        if display:
            measAlg.utils.showPsfMosaic(exposure, psf, frame=1)
            #measAlg.utils.showPsfMosaic(coadd, coadd.getPsf(), frame=1)

        # set PSF
        exposure.setPsf(psf)

        # ---------------------------------------------- #
        # Return exposure
        # ---------------------------------------------- #

        # Write exposure
        if True:
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
