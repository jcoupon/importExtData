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

import numpy as np
import errno
import os

from argparse import ArgumentError

import lsst.pex.config     as pexConfig
from lsst.pipe.tasks.coaddBase import CoaddBaseTask
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage

__all__ = ["CreateMultiCatTask"]

class CreateMultiCatConfig(CoaddBaseTask.ConfigClass):

    filters   = pexConfig.Field("Name of filters to combine [default HSC-G^HSC-R^HSC-I^HSC-Z^HSC-Y]", str, "HSC-G^HSC-R^HSC-I^HSC-Z^HSC-Y")
    dustCoefs = pexConfig.Field("Correction coefficient to compute dust correction [default 3.711^2.626^1.916^1.469^1.242]", str, "3.711^2.626^1.916^1.469^1.242")
    aperId    = pexConfig.Field("Aperture id", str, "2,3")

    fileOutName = pexConfig.Field("Name of output file", str, "")
    dirOutName  = pexConfig.Field("Name of output directory (will write output files as dirOutName/FILTER/TRACT/PATCH/multiCat-FILTER-TRACT-PATCH.fits)", str, "")

    dustSgpFileName = pexConfig.Field("Name of output file", str, "/Users/coupon/data/SchlegelDust/SFD_dust_4096_sgp.fits")
    dustNgpFileName = pexConfig.Field("Name of output file", str, "/Users/coupon/data/SchlegelDust/SFD_dust_4096_ngp.fits")

    clobber = pexConfig.Field("Replace output file", bool, False)

    def setDefaults(self):
        pexConfig.Config.setDefaults(self)

class CreateMultiCatTask(CoaddBaseTask):
    """A Task to merge catalogs
    """

    _DefaultName = 'CreateMultiCat'
    ConfigClass = CreateMultiCatConfig

    class dustMap(object):
        """Dust map info
        """
        def __init__(self):
            pass

    def __init__(self, schema=None, *args, **kwargs):
        CoaddBaseTask.__init__(self,  *args, **kwargs)

        if len(self.config.dustCoefs.split("^")) !=  len(self.config.filters.split("^")):
                raise ArgumentError(None, "filters and dustCoefs must have the same number of elements")

        # ---------------------------------------------------------- #
        # for Galactic extinction until https://hsc-jira.astro.princeton.edu/jira/browse/HSC-1350 is fixed
        # ---------------------------------------------------------- #
        from   astropy.io        import ascii,fits
        import astropy.wcs       as wcs

        sFile = fits.open(self.config.dustSgpFileName)
        nFile = fits.open(self.config.dustNgpFileName)

        self.dustMap.sMap  = sFile[0].data
        self.dustMap.nMap  = nFile[0].data

        self.dustMap.sWcs = wcs.WCS(sFile[0].header)
        self.dustMap.nWcs = wcs.WCS(nFile[0].header)
        # ---------------------------------------------------------- #
        #
        # ---------------------------------------------------------- #


    def iround(self, x):
        """iround(number) -> integer
        Round a number to the nearest integer.
        From https://www.daniweb.com/software-development/python/threads/299459/round-to-nearest-integer-
        """
        return int(round(x) - .5) + (x > 0)

    def readCatalog(self, dataRef, filterName):
        """Read input catalog

        We read the input dataset provided by the 'inputDataset'
        class variable.
        """
        dataRef.dataId["filter"] = filterName
        catalog = dataRef.get("deepCoadd_forced_src", immediate=True)
        self.log.info("Read %d sources for filter %s: %s" % (len(catalog), filterName, dataRef.dataId))
        return filterName, catalog

    def readCoadd(self, dataRef, filterName):
        """Read input coadd
        """
        dataRef.dataId["filter"] = filterName
        coadd      = dataRef.get("deepCoadd_calexp")
        self.log.info("Read coadd for filter %s: %s" % (filterName, dataRef.dataId))
        return filterName, coadd

    def getDustCorrection(self, dustMap, ra, dec):

        from astropy import units as u
        from astropy.coordinates import SkyCoord
        import astropy.wcs       as wcs

        coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='fk5')

        if coord.galactic.b.degree > 0.0:
            x, y = wcs.utils.skycoord_to_pixel(coord, dustMap.nWcs,  origin=0)
            if np.isfinite(x) and np.isfinite(y):
                return float(dustMap.nMap[self.iround(y), self.iround(x)])
            else:
                return 0.0
        else:
            x, y = wcs.utils.skycoord_to_pixel(coord, dustMap.sWcs,  origin=0)
            if np.isfinite(x) and np.isfinite(y):
                return float(dustMap.sMap[self.iround(y), self.iround(x)])
            else:
                return 0.0


    def run(self, dataRef, selectDataList=[]):

        if self.config.fileOutName == "":
            if self.config.dirOutName == "" :
                dirOutName = dataRef.getButler().mapper.root+"/"+self.config.coaddName+"Coadd-results"
                self.log.info("WARNING: the output file will be written in {0:s}.".format(dirOutName))
            else:
                dirOutName = self.config.dirOutName
            fileOutName = "{0}/{1}/{2}/{3}/multiCat-{2}-{3}.fits".format(dirOutName,"merged",dataRef.dataId["tract"],dataRef.dataId["patch"])
        else:
            fileOutName = self.config.fileOutName

        if os.path.isfile(fileOutName) and not self.config.clobber:
            self.log.info("File for  %s exists. Exiting..." % (dataRef.dataId))
            return

        self.log.info("Processing %s" % (dataRef.dataId))

        filters = self.config.filters.split("^")
        dustCoefs = [float(c) for c in self.config.dustCoefs.split("^") ]
        ref = dataRef.get("deepCoadd_ref")

        catalogs = dict(self.readCatalog(dataRef, f) for f in filters)
        coadds = dict(self.readCoadd(dataRef, f) for f in filters)

        # print ref.schema.getOrderedNames()
        # print dir(ref.schema)
        # print catalogs[filters[0]].schema.getOrderedNames()
        # return

        fluxMag0 = {}
        for f in filters:
            fluxMag0[f] = coadds[f].getCalib().getFluxMag0()[0]
            self.log.info("Mag ZP for filter {0:s}: {1:f}".format(f, 2.5*np.log10(fluxMag0[f])))

        wcs = coadds[filters[0]].getWcs()
        pixel_scale = wcs.pixelScale().asDegrees()*3600.0

        aperId = [int(a) for a in self.config.aperId.split(",")]

        """display which aperture diameter size will be used
        """
        aperSize = []
        for j, a in enumerate(aperId):
            aperSize.append(catalogs[filters[0]].getMetadata().get("flux_aperture_radii")[a] * 2.0 * pixel_scale)
            self.log.info("Diameter of flux apertures: {0:f}\"".format(aperSize[j]))

        """create new table table
        """
        mergedSchema = afwTable.Schema()

        """define table fields
        """
        fields=[]
        fields.append(mergedSchema.addField("id", type="L", doc="Unique id"))
        fields.append(mergedSchema.addField("ra", type="F", doc="ra [deg]"))
        fields.append(mergedSchema.addField("dec", type="F", doc="dec [deg]"))
        fields.append(mergedSchema.addField("tract", type="I", doc="tract number"))
        fields.append(mergedSchema.addField("patch", type="String", size=3, doc="patch number"))
        fields.append(mergedSchema.addField("countInputs", type="I", doc="Number of input single exposures for the reference filter"))
        fields.append(mergedSchema.addField("detRadius", type="F", doc="Determinant radius for the object in the reference filter = sigma if gaussian [arcsec]"))
        fields.append(mergedSchema.addField("PSFDetRadius", type="F", doc="Determinant radius for the PSF at the object position = sigma if gaussian [arcsec]"))
        fields.append(mergedSchema.addField("cmodel_fracDev", type="F", doc="fraction of flux in de Vaucouleur component"))
        fields.append(mergedSchema.addField("blendedness", type="F", doc="Ranges from 0 (unblended) to 1 (blended)"))
        fields.append(mergedSchema.addField("EB_V", type="F", doc="Milky Way dust E(B-V) [mag]"))
        fields.append(mergedSchema.addField("extendedness", type="F", doc="probability of being extended from PSF/cmodel flux difference"))
        fields.append(mergedSchema.addField("hasBadCentroid", type="I", doc="1 if has bad centroid (but not used in islean"))
        fields.append(mergedSchema.addField("isSky", type="I", doc="1 if sky object"))
        fields.append(mergedSchema.addField("isDuplicated", type="I", doc="1 if outside the inner tract or patch"))
        fields.append(mergedSchema.addField("isParent", type="I", doc="1 if parent of a deblended object"))
        fields.append(mergedSchema.addField("refFilter", type="String", size=10, doc="Name of the filter used as reference"))
        fields.append(mergedSchema.addField("isClean_refFilter", type="I", doc="1 if none of other flags is set for reference filter"))

        for f in filters:
            fields.append(mergedSchema.addField("hasBadPhotometry_{0:s}".format(f.replace(".", "_").replace("-", "_")), type="I", doc="1 if interpolated, saturated, suspect, has CR at center or near bright object for filter {0:s}".format(f)))
        for f in filters:
            fields.append(mergedSchema.addField("isEdge_{0:s}".format(f.replace(".", "_").replace("-", "_")), type="I", doc="1 if offImage or in region masked EDGE or NO_DATA for filter {0:s}".format(f)))

        """photometry estimates
        """
        photo = ["flux.aperture", "flux.kron", "flux.psf", "cmodel.flux"]
        for p in photo:
            for f in filters:
                if p == "flux.aperture":
                    for a in aperSize:
                        keyName = (p+"_"+str(int(a))+"arcsec_"+f).replace(".", "_").replace("-", "_")
                        fields.append(mergedSchema.addField(keyName,        type="F", doc="{0:s} for filter {1:s} within {2:f}\" diameter aperture".format(p,f,a)))
                        fields.append(mergedSchema.addField(keyName+"_err", type="F", doc="{0:s} error for filter {1:s}".format(p,f)))
                else:
                    keyName = (p+"_"+f).replace(".", "_").replace("-", "_")
                    fields.append(mergedSchema.addField(keyName,        type="F", doc="{0:s} for filter {1:s}".format(p,f)))
                    fields.append(mergedSchema.addField(keyName+"_err", type="F", doc="{0:s} error for filter {1:s}".format(p,f)))
            fields.append(mergedSchema.addField(p+"_flag".replace(".", "_"),   type="I", doc="Highest flag value among filters for {0:s}".format(p)))

        """dust corrections
        """
        for f in filters:
            fields.append(mergedSchema.addField(("EB_V_corr_"+f).replace(".", "_").replace("-", "_"),   type="F", doc="Milky way dust flux correction for filter {0:s}".format(f)))

        """create table object
        """
        merged = afwTable.BaseCatalog(mergedSchema)

        N = len(ref)
        # count = 0
        for i in range(N):
        # for i in range(10000,10200):
        # for i in range(1,100):

            """create new record
            """
            record = merged.addNew()
            coord = ref[i].get('coord')

            """record common info from reference filter
            """
            record.set(mergedSchema['id'].asKey(), ref[i].get('id'))
            record.set(mergedSchema['ra'].asKey(), coord.toFk5().getRa().asDegrees())
            record.set(mergedSchema['dec'].asKey(), coord.toFk5().getDec().asDegrees())
            record.set(mergedSchema['tract'].asKey(), dataRef.dataId["tract"])
            record.set(mergedSchema['patch'].asKey(), dataRef.dataId["patch"])
            record.set(mergedSchema['countInputs'].asKey(), ref[i].get('countInputs'))
            record.set(mergedSchema['detRadius'].asKey(), ref[i].get("shape.sdss").getDeterminantRadius()*pixel_scale)
            record.set(mergedSchema['PSFDetRadius'].asKey(), ref[i].get("shape.sdss.psf").getDeterminantRadius()*pixel_scale)
            record.set(mergedSchema['cmodel_fracDev'].asKey(), ref[i].get('cmodel.fracDev'))
            record.set(mergedSchema['blendedness'].asKey(), ref[i].get('blendedness.abs.flux'))
            record.set(mergedSchema['extendedness'].asKey(), ref[i].get('classification.extendedness'))
            record.set(mergedSchema['hasBadCentroid'].asKey(), int(ref[i].get('centroid.sdss.flags')))
            record.set(mergedSchema['isSky'].asKey(), int(ref[i].get('merge.footprint.sky')))
            record.set(mergedSchema['isDuplicated'].asKey(), int(not ref[i].get('detect.is-primary')))
            record.set(mergedSchema['isParent'].asKey(), int(ref[i].get('deblend.nchild') != 0))

            """record the name of the filter used as reference
            """
            for f in filters:
                name = afwImage.Filter(afwImage.Filter(f).getId()).getName()
                if ref[i].get("merge.measurement."+name):
                    refName = f.replace(".", "_").replace("-", "_")
                    record.set(mergedSchema["refFilter"].asKey(), refName)
                    break

            """photometry measurement flags
            """
            for p in photo:
                flag = 0
                for f in filters:
                    if catalogs[f][i].get(p+".flags") > flag:
                        flag = catalogs[f][i].get(p+".flags")
                record.set(mergedSchema[p+"_flag".replace(".", "_")].asKey(), flag)

            """record bad photometry flag for each filter
            """
            for f in filters:
                record.set(mergedSchema["isEdge_{0:s}".format(f.replace(".", "_").replace("-", "_"))].asKey(), int((catalogs[f][i].get('flags.pixel.offimage')) | (catalogs[f][i].get('flags.pixel.edge'))))
                record.set(mergedSchema["hasBadPhotometry_{0:s}".format(f.replace(".", "_").replace("-", "_"))].asKey(), int(
                                (catalogs[f][i].get('flags.pixel.interpolated.center')) \
                                |  (catalogs[f][i].get('flags.pixel.saturated.center')) \
                                |  (catalogs[f][i].get('flags.pixel.suspect.center'))  \
                                |  (catalogs[f][i].get('flags.pixel.cr.center')) \
                                |  (catalogs[f][i].get('flags.pixel.bad')) \
                                |  (catalogs[f][i].get('flags.pixel.bright.object.center'))))

            """record isClean flag for reference filter
            """
            hasBadPhotometry_refFilter = \
                (ref[i].get('flags.pixel.interpolated.center')) \
                |  (ref[i].get('flags.pixel.saturated.center')) \
                |  (ref[i].get('flags.pixel.suspect.center'))  \
                |  (ref[i].get('flags.pixel.cr.center')) \
                |  (ref[i].get('flags.pixel.bad')) \
                |  (ref[i].get('flags.pixel.bright.object.center'))
            isEdge_refFilter = (ref[i].get('flags.pixel.offimage')) | (ref[i].get('flags.pixel.edge'))
            isSky_refFilter = ref[i].get('merge.footprint.sky')
            isDuplicated_refFilter = not ref[i].get('detect.is-primary')
            isParent_refFilter = ref[i].get('deblend.nchild') != 0

            record.set(mergedSchema["isClean_refFilter"].asKey(), int(
                    (not hasBadPhotometry_refFilter) \
                    & (not isEdge_refFilter) \
                    & (not isSky_refFilter) \
                    & (not isDuplicated_refFilter) \
                    & (not isParent_refFilter)))

            #isExtended = (ref[i].get('classification.extendedness'))
            #isExtended = (catalogs[f][i].get('flux.kron') > 0.8*catalogs[f][i].get('flux.psf')) | \
            #              (ref.get("shape.sdss").getDeterminantRadius() > 1.1*ref.get("shape.sdss.psf").getDeterminantRadius())
            #if not (isExtended == isExtended):
            #    isExtended = 0
            #record.set(mergedSchema['isExtended'].asKey(), int(isExtended))

            """dust correction
            """
            EB_V = self.getDustCorrection(self.dustMap, record.get(mergedSchema['ra'].asKey()), record.get(mergedSchema['dec'].asKey()))
            record.set(mergedSchema['EB_V'].asKey(), EB_V)

            """flux in micro Jansky
            """
            for f, dc in zip(filters, dustCoefs):

                record.set(mergedSchema[("EB_V_corr_"+f).replace(".", "_").replace("-", "_")].asKey(), pow(10.0, +0.4* dc * EB_V))
                for p in photo:
                    if p == "flux.aperture":
                        for a in aperId:
                            flux = pow(10.0, 23.9/2.5) * catalogs[f][i].get(p)[a] / fluxMag0[f]
                            flux_err = pow(10.0, 23.9/2.5) * catalogs[f][i].get(p+".err")[a] / fluxMag0[f]
                            keyName = (p+"_"+str(int(a))+"arcsec_"+f).replace(".", "_").replace("-", "_")
                            record.set(mergedSchema[keyName].asKey(), flux)
                            record.set(mergedSchema[keyName+"_err"].asKey(), flux_err)
                    else:
                        flux = pow(10.0, 23.9/2.5)*catalogs[f][i].get(p)/fluxMag0[f]
                        flux_err = pow(10.0, 23.9/2.5)*catalogs[f][i].get(p+".err")/fluxMag0[f]
                        keyName = (p+"_"+f).replace(".", "_").replace("-", "_")
                        record.set(mergedSchema[keyName].asKey(), flux)
                        record.set(mergedSchema[keyName+"_err"].asKey(), flux_err)

            # count += 1

        """write catalog
        """
        self.log.info("Writing {0:s}".format(fileOutName))
        self.mkdir_p(os.path.dirname(fileOutName))
        merged.writeFits(fileOutName)

        return


    """ Don't forget to overload these
    """
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

if __name__ == '__main__':
    CreateMultiCatTask.parseAndRun()
