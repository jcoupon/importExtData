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

from lsst.pex.config import Config, Field, ListField, ConfigurableField, RangeField, ConfigField
from lsst.pipe.tasks.multiBand import MergeSourcesTask, MergeSourcesConfig, MergeSourcesRunner, _makeGetSchemaCatalogs, _makeMakeIdFactory


__all__ = ["createMultiCat"]


class CreateMultiCatConfig(MergeSourcesConfig):
    """Config for createMultiCatTask
    """
    coaddName   = Field(dtype=str, default="deep", doc="Name of coadd")
    fileOutName = Field("Name of output file", str, "multiCat.fits")

    def setDefaults(self):
        Config.setDefaults(self)

class CreateMultiCatTask(MergeSourcesTask):

    ConfigClass = CreateMultiCatConfig
    RunnerClass = MergeSourcesRunner
    _DefaultName = "createMultiCat"
    inputDataset = "meas"
    getSchemaCatalogs = _makeGetSchemaCatalogs("meas")
    makeIdFactory = _makeMakeIdFactory("CoaddId")

    def __init__(self, butler=None, schema=None, **kwargs):

        MergeSourcesTask.__init__(self, butler=butler, schema=schema, **kwargs)
        self.schema = self.getInputSchema(butler=butler, schema=schema)
        #self.merged = afwDetect.FootprintMergeList(
        #    self.schema,
        #    [getShortFilterName(name) for name in self.config.priorityList]
        #)



    def run(self, patchRefList=[]):
        """
        Task to create multiband catalogues
        """

        catalogs = dict(self.readCatalog(patchRef) for patchRef in patchRefList)
        print catalogs

        





        return

    def readCatalog(self, patchRef):
        """Read input catalog

        We read the input dataset provided by the 'inputDataset'
        class variable.
        """
        filterName = patchRef.dataId["filter"]
        catalog = patchRef.get(self.config.coaddName + "Coadd_" + self.inputDataset, immediate=True)
        self.log.info("Read %d sources for filter %s: %s" % (len(catalog), filterName, patchRef.dataId))
        return filterName, catalog




    # Overload these if your task inherits from CmdLineTask
    def _getConfigName(self):
        return None
    def _getEupsVersionsName(self):
        return None
    def _getMetadataName(self):
        return None
