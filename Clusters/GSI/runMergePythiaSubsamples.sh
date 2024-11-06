#!/bin/bash

# no more core files
ulimit -c 0
echo ulimit `ulimit -c`

if [ $# -gt 1 ]; then
  echo "usage: runMergePythiaSubsamples filename"
  exit 1
fi

if [ $# -lt 1 ]; then
  echo "usage: runMergePythiaSubsamples filename"
  exit 1
fi

FILENAME=$1

# setting the root and pythia scenario
export ALIEN_SITE=GSI
LATEST="VO_ALICE@ROOT::v6-32-06-alice1-4"
export ALIPHYSICS_VERSION=$LATEST

source <( /cvmfs/alice.cern.ch/bin/alienv printenv $LATEST)
echo $LATEST

# merge the subsamples in a single file
hadd ${FILENAME}.root BUNCH??/Output/${FILENAME}.root

