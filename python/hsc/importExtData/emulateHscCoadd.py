#!/usr/bin/env python


import random
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipebBase
import lsst.afw.table as afwTable
import lsst.meas.algorithms as measAlg


def foo():
    print "hello"
class FewerSourceDetectionConfig(measAlg.SourceDetectionConfig):
    nObjects = pexConfig.Field(doc="Number of sources to select", dtype=int, optional=False, default=200)

class FewerSourceDetectionTask(measAlg.SourceDetectionTask):
    """This task serves only to cull the source list and make measurement faster"""

    ConfigClass = FewerSourceDetectionConfig

    def makeSourceCatalog(self, table, exposure, doSmooth=True, sigma=None, clearMask=True):
        if self.negativeFlagKey is not None and self.negativeFlagKey not in table.getSchema():
            raise ValueError("Table has incorrect Schema")

        # detect the footprints as usual
        fpSets = self.detectFootprints(exposure=exposure, doSmooth=doSmooth, sigma=sigma,
                                       clearMask=clearMask)

        # shuffle the footprints to ensure they're random across the frame
        n = self.config.nObjects
        fpPos = fpSets.positive.getFootprints()
        random.shuffle(fpPos)

        # delete the excess footprints, and the negative footprints
        del fpPos[n:]
        fpSets.numPos = n
        if fpSets.negative:
            del fpSets.negative.getFootprints()[0:]
            fpSets.negative = None

        # make sources
        sources = afwTable.SourceCatalog(table)
        table.preallocate(fpSets.numPos)
        if fpSets.positive:
            fpSets.positive.makeSources(sources)

        return pipeBase.Struct(sources=sources, fpSets=fpSets)
