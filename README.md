# importExtData

This module adds non HSC data and makes it compatible with HSC pipeline.

To install it, download the latest version to path/to/importExtData
```
git clone https://github.com/jcoupon/importExtData.git
```
and run
```
cd path/to/importExtData
setup -v -j -r .
scons -Q -j 6 opt=3
cd -
```

If not present, the `data/_mapper` file must be put in `ROOT/`
and the `data/skyMap.pickle` file in `ROOT/deepCoadd/`.

## emulateHscCoadd.py

Command to emulate an image and a variance image into a LSST exposure object. For example:

```shell
emulateHscCoadd.py $ROOTDIR --rerun=$DIR_RERUN \
	--id tract=$TRACT patch=$PATCH filter=$FILTER \
	--config imgInName=$IMAGE varInName=$WEIGHT mag0=30.0 weight=True
```

Which must be followed by the detection step (if not run during the multiBand process):

```shell
detectCoaddSources.py $ROOTDIR --rerun=$DIR_RERUN:$DIR_RERUN
  	--id tract=$TRACT patch=$PATCH filter=$FILTER
```

And finally, when running `multiBandDriver.py`, the corresponding config file must be loaded:

```shell
multiBandDriver.py $ROOTDIR --rerun=emulate_test:emulate_test \
	--id tract=$TRACT patch=$PATCH filter=$FILTER  \
	-C $IMPORTEXTDATA/config/multiBandDriver[_debug].py
```

The `debug` version only keeps a fewer number of detections to speed up the tests.

# How to add a new filter?

## Mapper

Edit `python/hsc/importExtData/hscAndExtMapper.py`. For example:

```python
afwImage.utils.defineFilter(name='MegaCam-uS', lambdaEff=375, alias=['u1', 'u',])
for f in [ ..., 'MegaCam-uS']:
            self.filters[f] = afwImage.Filter(afwImage.Filter(f).getId()).getName()
```

## Filter mapping

The filter mapping helps to match the sources to the reference sources involved in the (astormetric and photometric) calibration. Note that if no re-calibration is run during the `emulateHscCoadd` phase, alternatively, `  --config measureCoaddSources.doMatchSources=False` can probably be set during the multiband process (not tested).

Create a new directory in `config` with the name of the camera and add the filter mapping to the reference objects:

```shell
mkdir config/MegaCam
touch config/MegaCam/filterMap.py 
```

And edit `config/MegaCam/filterMap.py`:

```python
for source, target in [
        ("MegaCam-u", 'g'),
        ("MegaCam-uS", 'g'),
    ]:
    config.filterMap[source] = target

```

Finally, add the filter map in the `multiBAndDriver.py` config file:

```python
config.measureCoaddSources.match.refObjLoader.load(
    os.path.join(os.environ["IMPORTEXTDATA_DIR"], "config", "MegaCam", "filterMap.py"))
```

## Filter priority list

Also in the `multiBandDriver.py` config file, update the filter priority list for `mergeCoaddDetections` and `mergeCoaddMeasurements` config. For example:

```python
config.mergeCoaddDetections.priorityList=[
    "HSC-I2", "HSC-I", "HSC-R2", "HSC-R", "HSC-Z", "HSC-Y", "HSC-G",
    "MegaCam-u", "MegaCam-uS",
    "VIRCAM-Y", "VIRCAM-J", "VIRCAM-H", "VIRCAM-Ks",
    "NB0921", "NB0816", "NB1010", "NB0387", "NB0515",
    ]
```

# Options

`imgInName`: name of input image (default: "")

`mskInName`: name of input mask, (default: "")

`varInName`: name of input variance image, (default: "")

`weight`: set if variance file is weight (default: False)

`test`: output LSST-made images (default: False)

`mag0`: magnitude AB zero point (default: 27)

`magLim`: Magnitude faint limit for PSF measurement (default: 23)

`charImage_package`: Package to load the charImage config from (default: obs_subaru)

# Coadd base options


```
usage: emulateHscCoadd.py input [options]

positional arguments:
  input                 path to input data repository, relative to
                        $PIPE_INPUT_ROOT

optional arguments:
  -h, --help            show this help message and exit
  --calib RAWCALIB      path to input calibration repository, relative to
                        $PIPE_CALIB_ROOT
  --output RAWOUTPUT    path to output data repository (need not exist),
                        relative to $PIPE_OUTPUT_ROOT
  --rerun [INPUT:]OUTPUT
                        rerun name: sets OUTPUT to ROOT/rerun/OUTPUT;
                        optionally sets ROOT to ROOT/rerun/INPUT
  -c [NAME=VALUE [NAME=VALUE ...]], --config [NAME=VALUE [NAME=VALUE ...]]
                        config override(s), e.g. -c foo=newfoo bar.baz=3
  -C [CONFIGFILE [CONFIGFILE ...]], --configfile [CONFIGFILE [CONFIGFILE ...]]
                        config override file(s)
  -L [LEVEL|COMPONENT=LEVEL [LEVEL|COMPONENT=LEVEL ...]], --loglevel [LEVEL|COMPONENT=LEVEL [LEVEL|COMPONENT=LEVEL ...]]
                        logging level; supported levels are
                        [trace|debug|info|warn|error|fatal]
  --longlog             use a more verbose format for the logging
  --debug               enable debugging output?
  --doraise             raise an exception on error (else log a message and
                        continue)?
  --noExit              Do not exit even upon failure (i.e. return a struct to
                        the calling script)
  --profile PROFILE     Dump cProfile statistics to filename
  --show SHOW [SHOW ...]
                        display the specified information to stdout and quit
                        (unless run is specified).
  -j PROCESSES, --processes PROCESSES
                        Number of processes to use
  -t TIMEOUT, --timeout TIMEOUT
                        Timeout for multiprocessing; maximum wall time (sec)
  --clobber-output      remove and re-create the output directory if it
                        already exists (safe with -j, but not all other forms
                        of parallel execution)
  --clobber-config      backup and then overwrite existing config files
                        instead of checking them (safe with -j, but not all
                        other forms of parallel execution)
  --no-backup-config    Don't copy config to file~N backup.
  --clobber-versions    backup and then overwrite existing package versions
                        instead of checkingthem (safe with -j, but not all
                        other forms of parallel execution)
  --no-versions         don't check package versions; useful for development
  --id [KEY=VALUE1[^VALUE2[^VALUE3...] [KEY=VALUE1[^VALUE2[^VALUE3...] ...]]
                        data ID, e.g. --id tract=12345 patch=1,2
  --selectId [KEY=VALUE1[^VALUE2[^VALUE3...] [KEY=VALUE1[^VALUE2[^VALUE3...] ...]]
                        data ID, e.g. --selectId visit=6789 ccd=0..9

Notes:
            * --config, --configfile, --id, --loglevel and @file may appear multiple times;
                all values are used, in order left to right
            * @file reads command-line options from the specified file:
                * data may be distributed among multiple lines (e.g. one option per line)
                * data after # is treated as a comment and ignored
                * blank lines and lines starting with # are ignored
            * To specify multiple values for an option, do not use = after the option name:
                * right: --configfile foo bar
                * wrong: --configfile=foo bar

```

## Full procedure to add external data

![alt text](https://github.com/jcoupon/importExtData/blob/master/doc/doc.001.png)
![alt text](https://github.com/jcoupon/importExtData/blob/master/doc/doc.002.png)
![alt text](https://github.com/jcoupon/importExtData/blob/master/doc/doc.003.png)
![alt text](https://github.com/jcoupon/importExtData/blob/master/doc/doc.004.png)
![alt text](https://github.com/jcoupon/importExtData/blob/master/doc/doc.005.png)
![alt text](https://github.com/jcoupon/importExtData/blob/master/doc/doc.006.png)
![alt text](https://github.com/jcoupon/importExtData/blob/master/doc/doc.007.png)
![alt text](https://github.com/jcoupon/importExtData/blob/master/doc/doc.008.png)
