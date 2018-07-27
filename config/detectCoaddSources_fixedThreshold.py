from lsst.meas.algorithms import SourceDetectionTask
config.detection.retarget(SourceDetectionTask)

config.detection.thresholdValue = 6.5
config.doScaleVariance = True
config.detection.thresholdType = 'pixel_stdev'

# config.detection.thresholdValue=7.0
# config.detection.thresholdType='stdev'
