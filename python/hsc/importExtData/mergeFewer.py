import lsst.pex.config as pexConfig
import lsst.pipe.tasks as pipeTask
from lsst.pipe.tasks.coaddBase import getSkyInfo
from lsst.pipe.tasks.multiBand import getShortFilterName

__all__ = ["FewerMergeDetectionsTask"]

class FewerMergeDetectionsConfig(pipeTask.multiBand.MergeDetectionsConfig):
    nObjects = pexConfig.Field(
        doc="Number of sources to select",
        dtype=int, optional=False, default=100)

class FewerMergeDetectionsTask(pipeTask.multiBand.MergeDetectionsTask):
    """This task serves only to cull the source list and make measurement faster"""

    _DefaultName = "FewerMergeCoaddMeasurements"
    ConfigClass = FewerMergeDetectionsConfig

    def mergeCatalogs(self, catalogs, patchRef):
        """!
        \brief Merge multiple catalogs.

        After ordering the catalogs and filters in priority order,
        \ref getMergedSourceCatalog of the \ref FootprintMergeList_ "FootprintMergeList" created by
        \ref \_\_init\_\_ is used to perform the actual merging. Finally, \ref cullPeaks is used to remove
        garbage peaks detected around bright objects.

        \param[in]  catalogs
        \param[in]  patchRef
        \param[out] mergedList
        """

        # print("test")

        # Convert distance to tract coordinate
        skyInfo = getSkyInfo(coaddName=self.config.coaddName, patchRef=patchRef)
        tractWcs = skyInfo.wcs
        peakDistance = self.config.minNewPeak / tractWcs.getPixelScale().asArcseconds()
        samePeakDistance = self.config.maxSamePeak / tractWcs.getPixelScale().asArcseconds()

        # Put catalogs, filters in priority order
        orderedCatalogs = [catalogs[band] for band in self.config.priorityList if band in catalogs.keys()]
        orderedBands = [getShortFilterName(band) for band in self.config.priorityList
                        if band in catalogs.keys()]

        mergedList = self.merged.getMergedSourceCatalog(orderedCatalogs, orderedBands, peakDistance,
                                                        self.schema, self.makeIdFactory(patchRef),
                                                        samePeakDistance)

        #
        # Add extra sources that correspond to blank sky
        #
        skySeed = patchRef.get(self.config.coaddName + "MergedCoaddId")
        skySourceFootprints = self.getSkySourceFootprints(mergedList, skyInfo, skySeed)
        if skySourceFootprints:
            key = mergedList.schema.find("merge_footprint_%s" % self.config.skyFilterName).key
            for foot in skySourceFootprints:
                s = mergedList.addNew()
                s.setFootprint(foot)
                s.set(key, True)

        # pick only the first nObjects merged sources (default:100)
        self.log.info("DEBUGGING: Keep {} sources".format(self.config.nObjects))
        mergedList = mergedList[:self.config.nObjects]

        # Sort Peaks from brightest to faintest
        for record in mergedList:
            record.getFootprint().sortPeaks()
        self.log.info("Merged to %d sources" % len(mergedList))
        # Attempt to remove garbage peaks
        self.cullPeaks(mergedList)
        return mergedList
