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
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask


import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.meas.algorithms as measAlg
import lsst.meas.base as measBase
import lsst.afw.table as afwTable



__all__ = ["EmulateHscCoaddTask"]


class InitialPsfConfig(pexConfig.Config):
    """Describes the initial PSF used for detection
    and measurement before we do PSF determination."""

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

    imgInName = pexConfig.Field("Name of input image", str, "")
    mskInName = pexConfig.Field("Name of input mask", str, "")
    varInName = pexConfig.Field("Name of input variance image", str, "")

    weight = pexConfig.Field("Set if variance file is weight", bool, False)

    test = pexConfig.Field("Output", bool, False)

    mag0 = pexConfig.Field("Magnitude zero point", float, 27.0)
    magLim = pexConfig.Field(
        "Magnitude faint limit for PSF measurement", float, 23.0)

    # filtName = pexConfig.Field("Filter name", str, None)

    initialPsf = pexConfig.ConfigField(
        dtype=InitialPsfConfig, doc=InitialPsfConfig.__doc__)

    detection  = pexConfig.ConfigurableField(
        target=measAlg.SourceDetectionTask,
        doc="Initial (high-threshold) detection phase for calibration",
    )

    measurement = pexConfig.ConfigurableField(
        target=measBase.SingleFrameMeasurementTask,
        doc="Initial measurements used to feed PSF determination and aperture correction determination",
    )

    measurePsf = pexConfig.ConfigurableField(target = MeasurePsfTask, doc = "")

    #
    # N.b. These configuration options only set the bitplane config.brightObjectMaskName
    # To make this useful you *must* also configure the flags.pixel algorithm, for example
    # by adding
    #   config.measurement.plugins["base_PixelFlags"].masksFpCenter.append("BRIGHT_OBJECT")
    #   config.measurement.plugins["base_PixelFlags"].masksFpAnywhere.append("BRIGHT_OBJECT")
    # to your measureCoaddSources.py and forcedPhotCoadd.py config overrides
    #
    doMaskBrightObjects = pexConfig.Field(dtype=bool, default=True,
                                          doc="Set mask and flag bits for bright objects?")
    brightObjectMaskName = pexConfig.Field(dtype=str, default="BRIGHT_OBJECT",
                                           doc="Name of mask bit used for bright objects")



    charImage = pexConfig.ConfigurableField(
        target=CharacterizeImageTask,
        doc="""Task to characterize a science exposure:
            - detect sources, usually at high S/N
            - estimate the background, which is subtracted from the image and returned as field "background"
            - estimate a PSF model, which is added to the exposure
            - interpolate over defects and cosmic rays, updating the image, variance and mask planes
            """,
    )



    def setDefaults(self):

        pexConfig.Config.setDefaults(self)

        self.detection.includeThresholdMultiplier = 10.0

        # self.detection.doFootprintBackground = False

        # self.measurePsf.reserve.fraction = 0.0
        self.measurePsf.reserve.fraction = 0.2

        fluxMag0 = pow(10.0, +0.4*self.mag0)

        self.measurePsf.starSelector['objectSize'].sourceFluxField = \
            "base_PsfFlux_flux"
        self.measurePsf.starSelector['objectSize'].fluxMin = \
            fluxMag0 * pow(10.0, -0.4*self.magLim)

        self.charImage.measurePsf.starSelector['objectSize'].fluxMin = \
            fluxMag0 * pow(10.0, -0.4*self.magLim)


class EmulateHscCoaddTask(CoaddBaseTask):

    ConfigClass  = EmulateHscCoaddConfig
    _DefaultName = "emulateHscCoadd"

    def __init__(self, *args, **kwargs):

        CoaddBaseTask.__init__(self, *args, **kwargs)

        self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.algMetadata = dafBase.PropertyList()

        self.makeSubtask("detection", schema=self.schema)
        self.makeSubtask(
            "measurement", schema=self.schema, algMetadata=self.algMetadata)
        self.makeSubtask("measurePsf", schema=self.schema)
        self.makeSubtask("charImage")


        if self.config.doMaskBrightObjects:
            mask = afwImage.Mask()
            try:
                self.brightObjectBitmask = 1 << mask.addMaskPlane(self.config.brightObjectMaskName)
            except pexExceptions.LsstCppException:
                raise RuntimeError(
                    "Unable to define mask plane for bright objects; \
                    planes used are %s" % mask.getMaskPlaneDict().keys())
            del mask




    def run(self, patchRef, selectDataList=[]):
        """ Task to import external data and transform
        it into LSST exposure object. The task also
        performs PSF measurement
        """

        # --------------------------------------------- #
        # output images for tests
        # --------------------------------------------- #

        if self.config.test:
            coadd = patchRef.get(self.config.coaddName + "Coadd_calexp")
            self.writeTest(coadd, dirName=".")
            return

        if self.config.imgInName == "":
            raise ValueError(
                'You must provide an input image (imgInName)')

        if self.config.varInName == "":
            raise ValueError(
                'You must provide an input variance image (varInName)')


        import numpy as np

        self.log.info("Processing %s" % (patchRef.dataId))

        # --------------------------------------------- #
        # wcs info
        # --------------------------------------------- #

        skyInfo = self.getSkyInfo(patchRef)

        # --------------------------------------------- #
        # load image
        # --------------------------------------------- #

        imgIn = afwImage.ImageF(self.config.imgInName)

        # --------------------------------------------- #
        # load variance
        # --------------------------------------------- #

        varIn = afwImage.ImageF(self.config.varInName)

        # if weight map provided
        # variance = inverse weight
        if self.config.weight:
            varIn.getArray()[:] = 1.0/varIn.getArray()[:]

        # get pixels with 0 or non finite variance (= no data observed)
        noDataIn = (varIn.getArray()[:] == 0) \
            | (np.logical_not(np.isfinite(varIn.getArray()[:])))

        # --------------------------------------------- #
        # load mask or build one from scratch
        # --------------------------------------------- #

        if self.config.mskInName != "":
            #mskIn = afwImage.ImageU(self.config.mskInName) # not tested
            mskIn = afwImage.Mask(self.config.mskInName) # not tested
        else:
            # create mask array
            self.log.info("Creating new mask")
            mskIn = afwImage.Mask(skyInfo.bbox)

            # set NO_DATA bit to mask
            # where noDataIn
            mask_labels = mskIn.getMaskPlaneDict()
            self.log.info("Mask labels:".format(mask_labels))
            noDataBit = mask_labels["NO_DATA"]

            # check if not already set in ref mask
            # (no longer used)
            # nopatchRefNotSet = mskIn.getArray()[:]&(1<<noDataBit) == 0
            mskIn.getArray()[noDataIn] += 2**noDataBit

            # needed for multiband
            mskIn.addMaskPlane('CLIPPED')
            mskIn.addMaskPlane('CROSSTALK')
            mskIn.addMaskPlane('INEXACT_PSF')
            mskIn.addMaskPlane('NOT_DEBLENDED')
            mskIn.addMaskPlane('REJECTED')
            mskIn.addMaskPlane('SENSOR_EDGE')
            mskIn.addMaskPlane('UNMASKEDNAN')

        # --------------------------------------------- #
        # create exposure
        # --------------------------------------------- #

        self.log.info("Creating exposure")

        # exposure object
        exposure = afwImage.ExposureF(skyInfo.bbox, skyInfo.wcs)
        maskedImage = afwImage.MaskedImageF(imgIn, mskIn, varIn)
        exposure.setMaskedImage(maskedImage)

        # set dummy coadd info
        expSchema = afwTable.ExposureTable.makeMinimalSchema()
        coaddInputs = afwImage.CoaddInputs(expSchema, expSchema)
        exposure.getInfo().setCoaddInputs(coaddInputs)

        # set filter
        filter = afwImage.Filter(patchRef.dataId['filter'])
        exposure.setFilter(filter)

        # set calib object
        fluxMag0 = pow(10.0, +0.4*self.config.mag0)
        calib = afwImage.Calib()
        calib.setFluxMag0(fluxMag0)
        exposure.setCalib(calib)

        # --------------------------------------------- #
        # add bright object mask
        # --------------------------------------------- #

        if self.config.doMaskBrightObjects:
            brightObjectMasks = self.readBrightObjectMasks(
                patchRef)
            self.setBrightObjectMasks(
                exposure, patchRef.dataId, brightObjectMasks)


        # ---------------------------------------------- #
        # Do PSF measurement on coadd
        # ---------------------------------------------- #

        results = self.charImage.characterize(exposure)
        exposure = results.exposure

        if False:
            self.log.info("Measuring PSF")

            self.installInitialPsf(exposure)

            # prepare table
            idFactory = afwTable.IdFactory.makeSimple()
            table = afwTable.SourceTable.make(self.schema, idFactory)
            table.setMetadata(self.algMetadata)

            # detect sources
            detRet = self.detection.makeSourceCatalog(table, exposure)
            sources = detRet.sources

            # measure moments
            self.measurement.measure(sources, exposure)

            # measure psf
            psfRet = self.measurePsf.run(exposure, sources, expId=0, matches=None)

            cellSet = psfRet.cellSet
            psf = psfRet.psf

            displayPsf = False
            if displayPsf:
                measAlg.utils.showPsf(psf, frame=1)
                sigma = psf.computeShape().getDeterminantRadius()
                print(sigma)

                # seems to be broken (?)
                # measAlg.utils.showPsfMosaic(exposure, psf, frame=1, showFwhm=True)

            # set PSF
            exposure.setPsf(psf)

        # ---------------------------------------------- #
        # write exposure
        # ---------------------------------------------- #

        butler = patchRef.butlerSubset.butler
        butler.put(exposure, self.config.coaddName + 'Coadd' , patchRef.dataId)

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
        coadd.getMaskedImage().getImage().writeFits(imgInName)
        coadd.getMaskedImage().getMask().writeFits(mskInName)
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

        fwhm = self.config.initialPsf.fwhm / wcs.getPixelScale().asArcseconds()
        size = self.config.initialPsf.size
        self.log.info("installInitialPsf fwhm=%.2f pixels; size=%d pixels" % (fwhm, size))
        psf = cls(size, size, fwhm/(2*math.sqrt(2*math.log(2))))
        exposure.setPsf(psf)

    def readBrightObjectMasks(self, dataRef):
        """Returns None on failure


        # ---------------------------------------------- #

        copied from AssembleCoaddTask. To be implemented so that it
        can be imported more easily

        see https://jira.lsstcorp.org/browse/DM-15030

        # ---------------------------------------------- #

        """

        try:
            return dataRef.get("brightObjectMask", immediate=True)
        except Exception as e:
            self.log.warn(
                "Unable to read brightObjectMask for %s: %s", dataRef.dataId, e)
            return None

    def setBrightObjectMasks(self, exposure, dataId, brightObjectMasks):
        """Set the bright object masks

        exposure:          Exposure under consideration
        dataId:            Data identifier dict for patch
        brightObjectMasks: afwTable of bright objects to mask

        # ---------------------------------------------- #

        copied from AssembleCoaddTask. To be implemented so that it
        can be imported more easily

        see https://jira.lsstcorp.org/browse/DM-15030

        # ---------------------------------------------- #

        """
        #
        # Check the metadata specifying the tract/patch/filter
        #
        if brightObjectMasks is None:
            self.log.warn("Unable to apply bright object mask: none supplied")
            return
        self.log.info("Applying %d bright object masks to %s", len(brightObjectMasks), dataId)
        md = brightObjectMasks.table.getMetadata()
        for k in dataId:
            if not md.exists(k):
                self.log.warn("Expected to see %s in metadata", k)
            else:
                if md.get(k) != dataId[k]:
                    self.log.warn("Expected to see %s == %s in metadata, saw %s", k, md.get(k), dataId[k])

        mask = exposure.getMaskedImage().getMask()
        wcs = exposure.getWcs()
        plateScale = wcs.getPixelScale().asArcseconds()

        for rec in brightObjectMasks:
            center = afwGeom.PointI(wcs.skyToPixel(rec.getCoord()))
            if rec["type"] == "box":
                assert rec["angle"] == 0.0, ("Angle != 0 for mask object %s" % rec["id"])
                width = rec["width"].asArcseconds()/plateScale    # convert to pixels
                height = rec["height"].asArcseconds()/plateScale  # convert to pixels

                halfSize = afwGeom.ExtentI(0.5*width, 0.5*height)
                bbox = afwGeom.Box2I(center - halfSize, center + halfSize)

                bbox = afwGeom.BoxI(afwGeom.PointI(int(center[0] - 0.5*width), int(center[1] - 0.5*height)),
                                    afwGeom.PointI(int(center[0] + 0.5*width), int(center[1] + 0.5*height)))
                spans = afwGeom.SpanSet(bbox)
            elif rec["type"] == "circle":
                radius = int(rec["radius"].asArcseconds()/plateScale)   # convert to pixels
                spans = afwGeom.SpanSet.fromShape(radius, offset=center)
            else:
                self.log.warn("Unexpected region type %s at %s" % rec["type"], center)
                continue
            spans.clippedTo(mask.getBBox()).setMask(mask, self.brightObjectBitmask)



    # Overload these if your task inherits from CmdLineTask
    def _getConfigName(self):
        return None
    def _getEupsVersionsName(self):
        return None
    def _getMetadataName(self):
        return None
