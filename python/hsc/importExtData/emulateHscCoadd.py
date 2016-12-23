#!/usr/bin/env python
#
# Jean Coupon (jean.coupon@unige.ch)
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/) and the HSC software team
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
import errno
import os

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

    imgInName = pexConfig.Field("Name of input image", str, "im.fits")
    mskInName = pexConfig.Field("Name of input mask", str, "msk.fits")
    varInName = pexConfig.Field("Name of input variance image", str, "var.fits")

    mskInRef = pexConfig.Field("Use mask image from reference image", bool, False)
    weight = pexConfig.Field("Set if variance file is weight", bool, False)

    # fileOutName = pexConfig.Field("Name of output file", str, "exposure.fits")

    fileOutName = pexConfig.Field("Name of output file", str, "")
    dirOutName  = pexConfig.Field("Name of output directory (will write output files as dirOutName/FILTER/TRACT/PATCH.fits)", str, "")

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

        self.detection.doFootprintBackground = False

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

        self.log.info("Processing %s" % (patchRef.dataId))


        if False:
            # exposure = afwImage.ExposureF("/Users/coupon/data/HSC/SSP/rerun/tutorial_3.6.1/deepCoadd/HSC-I/1/5,5.fits")
            #coadd = afwImage.ExposureF("/Volumes/dataTmp/HSC/SSP/rerun/tutorial_3.6.1/deepCoadd/HSC-I/1/5,5.fits")
            coadd = afwImage.ExposureF("/Users/coupon/Desktop/5,5.fits")
            measAlg.utils.showPsfMosaic(coadd, coadd.getPsf(), frame=1)
            return

        # ---------------------------------------------- #
        # Import reference coadd and records wcs info
        # ---------------------------------------------- #

        coadd    = patchRef.get(self.config.coaddName + "Coadd_calexp")
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
            # for tests
            imgInName, mskInName, varInName = self.writeTest(coadd, dirName="/Users/coupon/data/tmp")
            #imgInName, _ , _ = self.writeTest(coadd, dirName="/Users/coupon/data/tmp")
            fluxMag0 = coadd.getCalib().getFluxMag0()[0]
            exit(-1)
        else:
            imgInName, mskInName, varInName = self.config.imgInName, self.config.mskInName, self.config.varInName
            if self.config.mskInRef:
                mskIn = coadd.getMaskedImage().getMask()
            else:
                mskIn = afwImage.MaskU( mskInName)

        # ---------------------------------------------- #
        # record in mask where there's no data
        # ---------------------------------------------- #

        # first get bit value for NO_DATA
        mask = coadd.getMaskedImage().getMask()
        mask_labels = mask.getMaskPlaneDict()
        noDataBit = mask_labels["NO_DATA"]

        varIn = afwImage.ImageF(varInName)

        if self.config.weight:
            # noDataIn = (varIn.getArray()[:] == float('Inf')) | (varIn.getArray()[:] == -float('Inf')) | (varIn.getArray()[:] == float('NaN'))

            import numpy as np
            noDataIn = np.logical_not(np.isfinite(varIn.getArray()[:]))

            varIn.getArray()[:] = 1.0/varIn.getArray()[:]
        else:
            noDataIn = (varIn.getArray()[:] == 0) | (np.logical_not(np.isfinite(varIn.getArray()[:])))
        nopatchRefNotSet = mskIn.getArray()[:]&(1<<noDataBit) == 0 # check if not already set in ref mask

        mskIn.getArray()[noDataIn & nopatchRefNotSet ] += 2**noDataBit

        # add median variance of good pixles where it's 0
        varMedian = np.median(varIn.getArray()[mskIn.getArray() == 0])
        varIn.getArray()[noDataIn] = varMedian


        # ---------------------------------------------- #
        # create exposure
        # ---------------------------------------------- #


        ##### DEBUGGING
        # img = afwImage.ImageF("/Users/coupon/data/CLAUDS/Gwyn/DEEPv3/Mega-u_8767_8c1.fits")
        # weird = img.getArray()[:] > 8000.0
        # print len(img.getArray()[:][weird])
        # print np.max(img.getArray()), np.min(img.getArray())
        # print img.getArray()
        # return
        ##########

        exposure = afwImage.ExposureF(skyInfo.bbox, skyInfo.wcs)
        maskedImage = afwImage.MaskedImageF(
            afwImage.ImageF(imgInName),
            mskIn,
            varIn)
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

        # if len(sources) < 8:
        #    fwhm = self.config.initialPsf.fwhm / exposure.getWcs().pixelScale().asArcseconds()
        #    psf  = measAlg.DoubleGaussianPsf(15, 15, fwhm/(2*math.sqrt(2*math.log(2))))
        # else:
        #    cellSet = psfRet.cellSet
        #    psf = psfRet.psf

        display = False
        if display:
            measAlg.utils.showPsfMosaic(exposure, psf, frame=1)
            #measAlg.utils.showPsfMosaic(coadd, coadd.getPsf(), frame=1)

        # set PSF
        exposure.setPsf(psf)

        # ---------------------------------------------- #
        # Return exposure
        # ---------------------------------------------- #

        # DEBUGGING
        if False:
            image    = exposure.getMaskedImage().getImage().getArray()
            imageTmp = image.copy()
            image *= 0.0
            image[1000:2000, 1000:2000] = imageTmp[1000:2000, 1000:2000]

            image    = exposure.getMaskedImage().getMask().getArray()
            imageTmp = image.copy()
            image[:] = 1000
            image[1000:2000, 1000:2000] = imageTmp[1000:2000, 1000:2000]


        # Write exposure
        if self.config.fileOutName == "":
            if self.config.dirOutName == "" :
                dirOutName = patchRef.getButler().mapper.root+"/"+self.config.coaddName+"Coadd"
                self.log.info("WARNING: the output file will be written in {0:s}.".format(dirOutName))
            else:
                dirOutName = self.config.dirOutName

            fileOutName = "{0}/{1}/{2}/{3}.fits".format(dirOutName,self.config.filtName,patchRef.dataId["tract"],patchRef.dataId["patch"])
        else:
            fileOutName = self.config.fileOutName

        self.log.info("Writing {0:s}".format(fileOutName))

        self.mkdir_p(os.path.dirname(fileOutName))
        exposure.writeFits(fileOutName)

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

    def mkdir_p(self, path):
        try:
            os.makedirs(path)
        except OSError as exc:  # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise
