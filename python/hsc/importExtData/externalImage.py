import os
import numpy as np
from lsst.afw.geom import Box2I, SkyWcs
from lsst.afw.image import ExposureF, readMetadata
from lsst.meas.algorithms import SingleGaussianPsf


class ExternalImage:
    """A dummy class that's used to read in an external image and transform
    it into an LSST Exposure object.

    This should have a single static method, readFits, that returns an
    lsst.afw.image.ExposureF.

    Any significant processing (e.g. pixel interpolation) should probably
    be done in a custom ISR task instead of here so it can be done once and
    saved, instead of being done every time the image is loaded.
    """

    @staticmethod
    def readFits(path):
        directory, filename = os.path.split(path)

        # Figure out what camera/filter this is by parsing filename,
        # and customize the code below accordingly.

        bbox = Box2I(width, height)
        result = ExposureF(bbox)
        result.image.array = ...  # main image, as a [y, x] numpy.float32 array
        result.variance.array = ...  # variance image, as a [y, x] numpy.float32 array

        # This example includes masking NaN pixels as NO_DATA and pixels over
        # 1E5 counts as SAT.  External information about where bad pixels
        # should be preferred when available, and obviously that saturation
        # threshold is just an example.
        # Interpolating these bad pixels is handled by a later task.
        noDataBitMask = result.mask.getPlaneBitMask("NO_DATA")
        satBitMask = result.mask.getPlaneBitMask("SAT")
        result.mask.array |= noDataBitMask*np.isnan(result.image.array)
        result.mask.array |= satBitMask*(result.image.array > 1E5)

        # If you have a better guess at the PSF, we can find a way to use it.
        # But it'd be a good idea to at least put this in with a guess at the
        # seeing (RMS in pixels).
        result.setPsf(SingleGaussianPsf(seeingRMS))

        # Add a guess for the WCS, in this case assuming it's in the FITS
        # header of the first HDU.  Need to have something here, even if it
        # isn't very good (e.g. whatever comes from the telescope).
        metadata = readMetadata(filename)
        wcs = SkyWcs(metadata)
        result.setWcs(wcs)

        return result
