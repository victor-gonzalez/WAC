#!/bin/bash

if [ $# -gt 2 ]; then
  echo "usage: runMergePythiaResults mergedfname filepattern"
  exit 1
fi

if [ $# -lt 2 ]; then
  echo "usage: runMergePythiaResults mergedfname filepattern"
  exit 1
fi

MERGEDFNAME=$1
PATTERNNAME=$2

# no more core files
ulimit -c 0
echo ulimit `ulimit -c`

# setting the root and pythia scenario
export ALIEN_SITE=GSI
LATEST="VO_ALICE@ROOT::v6-32-06-alice1-4"
export ALIPHYSICS_VERSION=$LATEST

source <( /cvmfs/alice.cern.ch/bin/alienv printenv $LATEST)
echo $LATEST

# Execute application code
root -l -b -q "/lustre/alice/users/vgonzale/CLUSTERMODELWAC/TOOLS/recursiveResultFilesMerge.C(\"${MERGEDFNAME}\",\"${PATTERNNAME}\")"
