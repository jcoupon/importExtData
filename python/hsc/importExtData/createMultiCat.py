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
import numpy as np

#from lsst.pipe.base import CmdLineTask, Struct, TaskRunner, ArgumentParser, ButlerInitializedTaskRunner
from lsst.pex.config import Config, Field, ListField, ConfigurableField, RangeField, ConfigField
from lsst.pipe.tasks.multiBand import MergeSourcesTask, MergeSourcesConfig, MergeSourcesRunner, _makeGetSchemaCatalogs, _makeMakeIdFactory, getShortFilterName

from lsst.pipe.tasks.coaddBase import getSkyInfo


import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.meas.mosaic as measMosaic

from   lsst.afw.geom     import Box2D


__all__ = ["createMultiCat"]


class CreateMultiCatConfig(MergeSourcesConfig):
    """Config for createMultiCatTask
    """

    filterRef = Field("Name of reference filter for size measurement (default: HSC-I)", str, "HSC-I")

    aperId    = Field("Aperture id", int, 2)


    coaddName   = Field(dtype=str, default="deep", doc="Name of coadd")
    fileOutName = Field("Name of output file", str, "multiCat.fits")

    dustSgpFileName = Field("Name of output file", str, "/Users/coupon/data/SchlegelDust/SFD_dust_4096_sgp.fits")
    dustNgpFileName = Field("Name of output file", str, "/Users/coupon/data/SchlegelDust/SFD_dust_4096_ngp.fits")

    def setDefaults(self):
        Config.setDefaults(self)

class CreateMultiCatTask(MergeSourcesTask):

    #catName           = "forced_src"
    catName           = "src"
    ConfigClass       = CreateMultiCatConfig
    RunnerClass       = MergeSourcesRunner
    _DefaultName      = "createMultiCat"
    inputDataset      = catName
    getSchemaCatalogs = _makeGetSchemaCatalogs(catName)

    def __init__(self, butler=None, schema=None, **kwargs):

        MergeSourcesTask.__init__(self, butler=butler, schema=schema, **kwargs)
        self.schema = self.getInputSchema(butler=butler, schema=schema)

        # load dust map in South and North Galactic caps
        self.dustSgpWcs      = afwImage.makeWcs(afwImage.readMetadata(self.config.dustSgpFileName))
        self.dustSgpArray    = afwImage.ImageF(self.config.dustSgpFileName).getArray()
        self.dustNgpWcs      = afwImage.makeWcs(afwImage.readMetadata(self.config.dustNgpFileName))
        self.dustNgpArray    = afwImage.ImageF(self.config.dustNgpFileName).getArray()

    def readCoadd(self, patchRef):
        """Read input coadd
        """
        filterName = patchRef.dataId["filter"]
        coadd      = patchRef.get(self.config.coaddName + "Coadd")
        self.log.info("Read coadd for filter %s: %s" % (filterName, patchRef.dataId))
        return filterName, coadd

    def readSkyInfo(self, patchRef):
        """Read input coadd
        """
        filterName = patchRef.dataId["filter"]
        skyInfo = getSkyInfo(coaddName=self.config.coaddName, patchRef=patchRef)
        self.log.info("Read getSkyInfo for filter %s: %s" % (filterName, patchRef.dataId))
        return filterName, skyInfo

    def isOutside(self, coord, xy, skyInfo):
        """Returns True if inside patch AND tract
        """

        if (not coord == coord):
            return 0

        # inner patch
        innerFloatBBox = Box2D(skyInfo.patchInfo.getInnerBBox())
        isPatchInner   = innerFloatBBox.contains(xy)

        # inner tract
        sourceInnerTractId = skyInfo.skyMap.findTract(coord).getId()
        isTractInner       = sourceInnerTractId == skyInfo.tractInfo.getId()

        if (isPatchInner) & (isTractInner):
            return 0
        else:
            return 1


    def run(self, patchRefList=[]):
        """
        Task to create multiband catalogues
        """

        catalogs = dict(self.readCatalog(patchRef) for patchRef in patchRefList)
        coadds   = dict(self.readCoadd(patchRef)   for patchRef in patchRefList)
        skyInfo  = dict(self.readSkyInfo(patchRef) for patchRef in patchRefList)

        #print catalogs[self.config.filterRef].schema
        #print catalogs[self.config.filterRef].schema.getOrderedNames()

        filters  = coadds.keys()

        if self.config.filterRef in filters:
            wcs         = coadds[self.config.filterRef].getWcs()
            pixel_scale = wcs.pixelScale().asDegrees()*3600.0
            skyInfo     = skyInfo[self.config.filterRef]
            aperSize    = catalogs[self.config.filterRef].getMetadata().get("flux_aperture_radii")[self.config.aperId] * 2.0 * pixel_scale
            self.log.info("Diameter of flux apertures: {0:f}\"".format(aperSize))

        else:
            raise ValueError("Reference filter \"{0:s}\" not found in filter list. Exiting...".format(self.filterRef))

        fluxMag0 = {}
        for f in filters:
            fluxMag0[f] = coadds[f].getCalib().getFluxMag0()[0]
            self.log.info("Mag ZP for filter {0:s}: {1:f}".format(f, 2.5*np.log10(fluxMag0[f])))


        # create new table table
        mergedSchema = afwTable.Schema()



        fields=[]
        # define table fields
        fields.append(mergedSchema.addField("id",               type="L", doc="Unique id"))
        fields.append(mergedSchema.addField("ra",               type="F", doc="ra [deg]"))
        fields.append(mergedSchema.addField("dec",              type="F", doc="dec [deg]"))
        fields.append(mergedSchema.addField("eB_V",             type="F", doc="Galactic extinction [mag]"))
        fields.append(mergedSchema.addField("countInputs",      type="I", doc="Number of input single exposures for the reference filter"))
        fields.append(mergedSchema.addField("detRadius",        type="F", doc="Determinant radius for the object in the reference filter = sigma if gaussian [arcsec]"))
        fields.append(mergedSchema.addField("PSFDetRadius",     type="F", doc="Determinant radius for the PSF at the object position = sigma if gaussian [arcsec]"))
        fields.append(mergedSchema.addField("cmodel_fracDev",   type="F", doc="fraction of flux in de Vaucouleur component"))
        fields.append(mergedSchema.addField("isOutside",        type="I", doc="1 if outside the inner tract or patch"))
        fields.append(mergedSchema.addField("isOnTheEdge",      type="I", doc="1 if on the EDGE or NO_DATA area"))
        fields.append(mergedSchema.addField("hasBadPhotometry", type="I", doc="1 if interpolated, saturated, suspect, or has CR at center"))
        fields.append(mergedSchema.addField("isParent",         type="I", doc="1 if parent of a deblended object"))
        fields.append(mergedSchema.addField("isClean",          type="I", doc="1 if none of other flags is set"))
        fields.append(mergedSchema.addField("isExtended",       type="F", doc="1 if cmodel_flux > flux_psf"))

        photo = ["flux.aperture", "flux.kron", "flux.psf", "cmodel.flux"]
        for p in photo:
            for f in filters:
                keyName = (p+"_"+f).replace(".", "_").replace("-", "_")
                fields.append(mergedSchema.addField(keyName,        type="F", doc="{0:s} for filter {1:s}".format(p,f)))
                fields.append(mergedSchema.addField(keyName+"_err", type="F", doc="{0:s} error for filter {1:s}".format(p,f)))

#                fields.append(mergedSchema.addField(p+"_"+f, type="F", doc="Aperture flux in {0:f}\" diam. apertures for filter {1:s}".format(aperSize,f)))

        #print mergedSchema; return

        # create table object
        merged = afwTable.BaseCatalog(mergedSchema)

        N = len(catalogs[self.config.filterRef])

        for i, ref in enumerate(catalogs[self.config.filterRef]):

#       for i in range(8000, 9000):
#       i = 8506
#       if True:
#            ref = catalogs[self.config.filterRef][i]

#           if (i > 0) & (i % 1000 == 0):
#               self.log.info("Processed {0:d} objects".format(i))

            # create new record
            record = merged.addNew()
            coord = ref.get('coord')

            # record if any of the filter is flagged as bad photometry
            for f in filters:
                hasBadPhotometry = (catalogs[f][i].get('flags.pixel.interpolated.center')) | (catalogs[f][i].get('flags.pixel.saturated.center')) \
                                  | (catalogs[f][i].get('flags.pixel.suspect.center'))      | (catalogs[f][i].get('flags.pixel.cr.center'))
                if hasBadPhotometry:
                    break

            isOutside   = self.isOutside(coord, wcs.skyToPixel(coord), skyInfo)
            isParent    = ref.get('deblend.nchild') != 0
            isOnTheEdge = ref.get('flags.pixel.edge')

            isClean = (not hasBadPhotometry) & (not isOutside) & (not isParent) & (not isOnTheEdge)

            #isExtended = (ref.get('classification.extendedness')) &
            isExtended = (catalogs[f][i].get('flux.kron') > 0.8*catalogs[f][i].get('flux.psf')) | \
                          (ref.get("shape.sdss").getDeterminantRadius() > 1.1*ref.get("shape.sdss.psf").getDeterminantRadius())
            if not (isExtended == isExtended):
                isExtended = 0

            # record common info from reference filter
            record.set(mergedSchema['id'].asKey(),               ref.get('id'))
            record.set(mergedSchema['ra'].asKey(),               coord.toFk5().getRa().asDegrees())
            record.set(mergedSchema['dec'].asKey(),              coord.toFk5().getDec().asDegrees())
            record.set(mergedSchema['countInputs'].asKey(),      ref.get('countInputs'))
            record.set(mergedSchema['detRadius'].asKey(),        ref.get("shape.sdss").getDeterminantRadius()*pixel_scale)
            record.set(mergedSchema['PSFDetRadius'].asKey(),     ref.get("shape.sdss.psf").getDeterminantRadius()*pixel_scale)
            record.set(mergedSchema['cmodel_fracDev'].asKey(),   ref.get('cmodel.fracDev'))
            record.set(mergedSchema['isOutside'].asKey(),        int(isOutside))
            record.set(mergedSchema['isOnTheEdge'].asKey(),      int(isOnTheEdge))
            record.set(mergedSchema['hasBadPhotometry'].asKey(), int(hasBadPhotometry))
            record.set(mergedSchema['isParent'].asKey(),         int(isParent))
            record.set(mergedSchema['isClean'].asKey(),          int(isClean))
            record.set(mergedSchema['isExtended'].asKey(),       int(isExtended))

            # flux in micro Jansky
            for f in filters:
                for p in photo:
                    if p == "flux.aperture":
                        flux     = pow(10.0, 23.9/2.5) * catalogs[f][i].get(p)[self.config.aperId] / fluxMag0[f]
                        flux_err = pow(10.0, 23.9/2.5) * catalogs[f][i].get(p+".err")[self.config.aperId] / fluxMag0[f]
                    else:
                        flux     = pow(10.0, 23.9/2.5) * catalogs[f][i].get(p) / fluxMag0[f]
                        flux_err = pow(10.0, 23.9/2.5) * catalogs[f][i].get(p+".err") / fluxMag0[f]

                    keyName =  (p+"_"+f).replace(".", "_").replace("-", "_")
                    record.set(mergedSchema[keyName].asKey(),          flux)
                    record.set(mergedSchema[keyName+"_err"].asKey(),   flux_err)

            # Galactic extinction (currently not implemented)
            if False:
                if coord.toGalactic().getB().asDegrees() < 0.0:
                    x, y = self.dustSgpWcs.skyToPixel(coord)
                    eB_V = self.dustSgpArray[y, x]
                else:
                    x, y = self.dustNgpWcs.skyToPixel(coord)
                    eB_V = self.dustNgpArray[y, x]
            else:
                eB_V = 0.0
            record.set(mergedSchema['eB_V'].asKey(),   eB_V)



        self.log.info("Writing {0:s}".format(self.config.fileOutName))
        merged.writeFits(self.config.fileOutName)

        return


    # Overload these if your task inherits from CmdLineTask
    def _getConfigName(self):
        return None
    def _getEupsVersionsName(self):
        return None
    def _getMetadataName(self):
        return None
