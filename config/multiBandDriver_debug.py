from hsc.importExtData import mergeFewer

config.mergeCoaddDetections.retarget(mergeFewer.FewerMergeDetectionsTask)

import os
obsSubaru = os.environ["OBS_SUBARU_DIR"]
overrides = [os.path.join(obsSubaru, "config", "mergeCoaddDetections.py"),
             os.path.join(obsSubaru, "config", "hsc", "mergeCoaddDetections.py"),
             ]
for filename in overrides:
    if os.path.exists(filename):
        config.mergeCoaddDetections.load(filename)

config.load(os.path.join(os.environ["IMPORTEXTDATA_DIR"], "config", "multiBandDriver.py"))
