import os
import re
import numpy as np
import sqlite3
from lsst.afw.geom import Point2I, Extent2I, Box2I, SkyWcs
from lsst.afw.image import ExposureF, readMetadata
from lsst.meas.algorithms import SingleGaussianPsf
from lsst.daf.persistence import Butler


__all__ = ("VisitMetadata", "ExternalImage")


class VisitMetadata:
    """Metadata associated with a visit that is stored in a SQL Registry.

    Parameters
    ----------
    dateObs : str
        Date of the observation, formatted as YYYY-MM-DD.  Must not be None.
    filter : str
        Full name of the bandpass filter.  Must not be None.
    expId : str
        Internal ID for the full-focal-plane exposure, prefixed by the
        camera name (HSC itself uses a string like "HSCA000NNNNNN").
        May be None, but probably useful to make it meaningful.
    field : str
        Name of the field being observed.  It may be convenient to
        use the same field names as HSC SSP observations for overlapping
        fields (e.g. SSP_UDEEP_COSMOS), but this is not actually necessary.
        May not be None.
    exptime : float
        Exposure time in seconds.  Can probably be None, but probably useful
        to use the real value if available.
    pointing : int
        An integer identifier for a set of observations.  Has a special
        meaning for HSC data, but may be set to any useful number
        for external data.  May not be None, but defaults to zero and
        this default can probably be used for all external data.
    dataType : string
        A string indicating the type of observation.
        Can probably be None, but letting it default to "OBJECT" is probably
        better.
    pa : double
        Rotator angle.  Defaults to zero, and should not generally need to
        be set to anything else for external data.
    """

    def __init__(self, dateObs, filter, expId, field, expTime,
                 pointing=0, dataType="OBJECT", pa=0.0):
        self.dateObs = dateObs
        self.filter = filter
        self.expId = expId
        self.field = field
        self.expTime = expTime
        self.pointing = pointing
        self.dataType = dataType
        self.pa = pa

    @property
    def taiObs(self):
        return self.dateObs


class ExternalImage:
    """A collection of class methods for operating on external images.

    ExternalImage should not be instantiated; it has only class methods.

    This class assumes that each supported camera has an entry in the
    CAMERA_INFO dictionary, providing an ID between 1 and 9.  These IDs
    are scaled by CAMERA_ID_MULTIPLIER and added to the visit number to
    avoid conflicts with actualy HSC visits.  Note that this requires
    all external visit numbers to be less than CAMERA_ID_MULTIPLIER.

    The operations on this class include:

     - readFits: called by the Butler when reading datasets of the "external"
       DatasetType.  This should be customized to support all desired external
       data formats.  The only requirement is that this should read a single
       FITS file and return an afw.image.Exposure object that can be handled
       by ExternalIsrTask.

     - ingest: add an external image (one that can be read by readFits) to
       a data repository by symlinking the file to the right location and
       adding metadata to the SQL registry.

    """

    EXTERNAL_REGEX = re.compile(r'external/(?P<visit>\d+)-(?P<ccd>)\.fits')

    CAMERA_ID_MULTIPLIER = 1000000

    CAMERA_INFO = {
        "megacam": {
            "id": 1,   # This is multiplied by CAMERA_ID_MULTIPLIER and added
                       # to the original visit ID to avoid conflicts with HSC
                       # visit numbers.
        },
        # add more cameras here
    }

    @classmethod
    def getCameraFromVisit(cls, visitId):
        cameraId = visitId % cls.CAMERA_ID_MULTIPLIER
        for camera, info in cls.CAMERA_INFO.items():
            if info["id"] == cameraId:
                return camera
        raise LookupError("Camera with ID=%d not found" % cameraId)

    @classmethod
    def readFits(cls, path):
        """Read an external detrended image and create an LSST Exposure object
        from it.

        Any significant processing (e.g. pixel interpolation) should probably
        be done in a ExternalIsrTask instead of here so it can be done once and
        saved, instead of being done every time the image is loaded.

        THIS METHOD IS INCOMPLETE; IT MUST BE MODIFIED ACCORDING TO THE
        FORMAT OF THE DATA BEING LOADED.
        """
        directory, filename = os.path.split(path)
        match = cls.EXTERNAL_REGEX.match(filename)
        camera = cls.getCameraFromVisit(match.group("visit"))

        # Customize the code below based on the camera determined above.
        # To support more than one camera it may be useful to delegate
        # to other methods that are specific to certain cameras.

        # Read the actual image in from the given path using e.g. astropy,
        # and use it to fill in various arrays below.

        bbox = Box2I(Point2I(0, 0), Extent2I(..., ...))  # width, height
        result = ExposureF(bbox)
        # main image, as a [y, x] numpy.float32 array
        result.image.array = ...
        # variance image, as a [y, x] numpy.float32 array
        result.variance.array = ...

        # This example includes masking NaN pixels as NO_DATA and pixels above
        # 1E5 counts as SAT.  External information about where bad pixels
        # should be preferred when available, and obviously that saturation
        # threshold is just an example (saturation should actually be
        # determined before flat-fielding, of course).
        # Interpolating these bad pixels is handled by ExternalIsrTask.
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

    @classmethod
    def ingest(cls, root, camera, visit, filenames, sensors, metadata):
        """Add all images from an external visit (a full-focal-plane
        exposure) to a data repository.

        This both symlinks the external data files to the appropriate
        location in the directory structure and adds the necessary
        rows to the SQLite registry tables.

        Parameters
        ----------
        root : str
            Directory of the data repository to add data to.  Must have
            an existing "registry.sqlite3" file present directly in the
            root and a _mapper file pointing to HscAndExtMapper.
        camera : str
            Name of the camera used to produced the external observation.
            Must have an entry in ExternalImage.CAMERA_INFO.
        visit : int
            Original integer visit ID for the observation, *before* adding
            CAMERA_INFO[camera]["ID"]*CAMERA_ID_MULTIPLIER.
        filenames : list
            A list of file names containing the external data files, either
            relative to the current directory or absolute.
        sensors : list
            A list of integer sensor IDs corresponding to the filenames list.
        metadata : VisitMetadata
            An object containing additional metadata for this visit to be
            added to the registry.  See VisitMetadata for a description of
            what attributes are required.
        """
        db = sqlite3.connect(os.path.join(root, "registry.sqlite3"))
        butler = Butler(inputs=[root])
        visit += cls.CAMERA_INFO[camera]["id"]*cls.CAMERA_ID_MULTIPLIER
        ccdCols = ["filter", "dateObs", "taiObs", "field", "expId", "pointing",
                   "dataType", "pa"]
        ccdSql = "INSERT INTO raw (visit, ccd, {}) VALUES (?, ?, {})".format(
            ", ".join(ccdCols),
            ", ".join(["?"] * len(ccdCols))
        )
        ccdValues = tuple(getattr(metadata, col) for col in ccdCols)
        visitCols = ["filter", "dateObs", "taiObs", "field"]
        visitSql = "INSERT INTO raw_visit (visit, {}) VALUES (?, {})".format(
            ", ".join(visitCols),
            ", ".join(["?"] * len(visitCols))
        )
        visitValues = tuple(getattr(metadata, col) for col in visitCols)
        for filename, sensor in zip(filenames, sensors):
            outputFileName = butler.get("external_filename", visit=visit,
                                        ccd=sensor)[0]
            os.symlink(filename, outputFileName)
            db.execute(ccdSql, (visit, sensor,) + ccdValues)
        db.execute(visitSql, (visit,) + visitValues)
        db.commit()
        return visit
