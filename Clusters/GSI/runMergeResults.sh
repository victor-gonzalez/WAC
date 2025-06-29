#!/bin/bash

if [ $# -gt 2 ]; then
  echo "usage: runMergeResults mergedfname filepattern"
  exit 1
fi

if [ $# -lt 2 ]; then
  echo "usage: runMergeResults mergedfname filepattern"
  exit 1
fi

MERGEDFNAME=$1
PATTERNNAME=$2

# no more core files
ulimit -c 0
echo ulimit `ulimit -c`

# setting the proper environment scenario
NEEDED="VO_ALICE@O2Physics::daily-20250614-0000-1"

source <( /cvmfs/alice.cern.ch/bin/alienv printenv $NEEDED)
echo $NEEDED

# Execute application code
root -l -b -q "/lustre/alice/users/$USER/CLUSTERMODELWAC/TOOLS/recursiveResultFilesMerge.C(\"${MERGEDFNAME}\",\"${PATTERNNAME}\")"
