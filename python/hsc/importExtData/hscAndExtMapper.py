#!/usr/bin/env python

from lsst.obs.hsc import HscMapper
import lsst.afw.image as afwImage


class HscAndExtMapper(HscMapper):
    """Provides abstract-physical mapping for HSC + external data"""

    def __init__(self, **kwargs):

        HscMapper.__init__(self, **kwargs)



        # add filters
        afwImage.utils.defineFilter(name='MegaCam-u', lambdaEff=375, alias=['u',])
        for f in ['MegaCam-u']:
            self.filters[f] = afwImage.Filter(afwImage.Filter(f).getId()).getName()

        #print dir(afwImage.Filter("HSC-G"))
        #print self.filters
