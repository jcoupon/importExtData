from lsst.pex.config import Config, ListField, Field
from lsst.pipe.base import CmdLineTask, Struct
from lsst.ip.isr import interpolateFromMask


__all__ = ("ExternalIsrConfig", "ExternalIsrTask")


class ExternalIsrConfig(Config):
    datasetType = "external"
    maskPlanesToInterpolate = ListField(
        dtype=str,
        doc="List of mask plane names that should be interpolated.",
        default=["NO_DATA", "SAT"],
    )
    fwhm = Field(
        dtype=float,
        doc="FWHM of PSF (arcsec) for interpolation (just a placeholder, not yet used).",
        default=1.0,
    )
    doWrite = Field(
        dtype=bool,
        doc="Persist postISRCCD?",
        default=True,
    )


class ExternalIsrTask(CmdLineTask):
    """An extremely simple replacement for IsrTask for externally-detrended data.

    All this Task does is read in the 'external' dataset and interpolate mask
    planes.  It should be extended to do any additional work that is too
    expensive for ExternalImage.readFits.
    """

    def runDataRef(self, sensorRef):
        self.log.info("Performing ISR on sensor %s" % (sensorRef.dataId))
        exposure = sensorRef.get(self.config.datasetType)
        result = self.run(exposure)
        if self.config.doWrite:
            sensorRef.put(result.exposure, "postISRCCD")

    def run(self, exposure):
        for maskPlane in self.config.maskPlanesToInterpolate:
            interpolateFromMask(exposure.maskedImage, fwhm=self.config.fwhm, maskName=maskPlane)
        result = Struct(exposure=exposure)
        return result
