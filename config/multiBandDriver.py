import os

config.measureCoaddSources.doPropagateFlags = False

config.measureCoaddSources.match.refObjLoader.load(
    os.path.join(os.environ["IMPORTEXTDATA_DIR"], "config", "MegaCam", "filterMap.py"))

config.measureCoaddSources.match.refObjLoader.load(
    os.path.join(os.environ["IMPORTEXTDATA_DIR"], "config", "VIRCAM", "filterMap.py"))

config.mergeCoaddDetections.priorityList=[
    "HSC-I2", "HSC-I", "HSC-R2", "HSC-R", "HSC-Z", "HSC-Y", "HSC-G",
    "MegaCam-u", "MegaCam-uS",
    "VIRCAM-Y", "VIRCAM-J", "VIRCAM-H", "VIRCAM-Ks",
    "NB0921", "NB0816", "NB1010", "NB0387", "NB0515",
    ]

config.mergeCoaddMeasurements.priorityList=[
    "HSC-I2", "HSC-I", "HSC-R2", "HSC-R", "HSC-Z", "HSC-Y", "HSC-G",
    "MegaCam-u", "MegaCam-uS",
    "VIRCAM-Y", "VIRCAM-J", "VIRCAM-H", "VIRCAM-Ks",
    "NB0921", "NB0816", "NB1010", "NB0387", "NB0515",
    ]

config.measureCoaddSources.deblend.maxFootprintArea=100000
config.measureCoaddSources.deblend.maskLimits["BRIGHT_OBJECT"] = 0.20
