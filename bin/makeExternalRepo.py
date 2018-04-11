#!/usr/bin/env python

import os
from lsst.daf.persistence import RepositoryCfg, Policy, Butler, PosixStorage


def main(parent, root):
    # Construct a Butler just to make a root directory with a
    # repositoryCfg.yaml file.
    Butler(inputs=[parent], outputs=[root])
    # Load that YAML file, and create a copy that adds the extra policy
    # items from this package.
    oldConfig = PosixStorage.getRepositoryCfg(root)
    newPolicy = Policy(Policy.defaultPolicyFile("importExtData",
                                                "importExtData.yaml",
                                                "policy"))
    newConfig = RepositoryCfg(root=oldConfig.root, mapper=oldConfig.mapper,
                              mapperArgs=oldConfig.mapperArgs,
                              parents=oldConfig.parents,
                              policy=newPolicy)
    # Remove the old YAML file.
    os.remove(os.path.join(root, "repositoryCfg.yaml"))
    # Write the new one in its place.
    PosixStorage.putRepositoryCfg(newConfig)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description=("Create an HSC data repository with extra dataset "
                     "definitions for imported external datasets.")
    )
    parser.add_argument("parent", metavar="PARENT", type=str,
                        help=("Directory containing the parent repository to "
                              "use for inputs (a regular HSC data "
                              "repository).  May be a rerun."))
    parser.add_argument("root", metavar="ROOT", type=str,
                        help=("Directory that will contain the new repository "
                              "with the new dataset definitions."))
    args = parser.parse_args()
    main(parent=args.parent, root=args.root)
