#!/usr/bin/env python

from lsst.obs.hsc import HscMapper
import lsst.afw.image as afwImage
from lsst.daf.persistence import Policy, RepositoryCfg

class HscAndExtMapper(HscMapper):
    """Provides abstract-physical mapping for HSC + external data"""

    @staticmethod
    def makeNewConfig(oldConfig):
        newPolicy = Policy(Policy.defaultPolicyFile("importExtData",
                                                    "importExtData.yaml",
                                                    "policy"))
        return RepositoryCfg(root=oldConfig.root,
                             mapper=oldConfig.mapper,
                             mapperArgs=oldConfig.mapperArgs,
                             parents=oldConfig.parents,
                             policy=newPolicy)

    def __init__(self, **kwargs):

        # Inject new mappings from importExtData's policy file.
        # This is a bit of a hack; we're pretending these policy entries
        # come from the configuration inside the repository itself, since
        # those always override and extend those from the camera's definitions.
        # Luckily, no one actually uses those per-repository policy entries
        # for anything else, so this should be safe.
        kwargs["repositoryCfg"] = self.makeNewConfig(kwargs["repositoryCfg"])

        HscMapper.__init__(self, **kwargs)

        # add filters
        afwImage.utils.defineFilter(name='MegaCam-uS', lambdaEff=375, alias=['u1', 'u',])
        afwImage.utils.defineFilter(name='MegaCam-u', lambdaEff=375, alias=['u2',])
        afwImage.utils.defineFilter(name='VIRCAM-Y', lambdaEff=375, alias=['Y','y',])

        for f in ['MegaCam-uS', 'MegaCam-u', 'VIRCAM-Y']:
            self.filters[f] = afwImage.Filter(afwImage.Filter(f).getId()).getName()

        #print dir(afwImage.Filter("HSC-G"))
        #print self.filters
