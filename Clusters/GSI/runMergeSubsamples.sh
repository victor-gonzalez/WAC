#!/bin/bash

# no more core files
ulimit -c 0
echo ulimit `ulimit -c`

if [ $# -gt 1 ]; then
  echo "usage: runMergeSubsamples filename"
  exit 1
fi

if [ $# -lt 1 ]; then
  echo "usage: runMergeSubsamples filename"
  exit 1
fi

FILENAME=$1

# setting the proper environment scenario
NEEDED="VO_ALICE@O2Physics::daily-20250614-0000-1"

source <( /cvmfs/alice.cern.ch/bin/alienv printenv $NEEDED)
echo $NEEDED

# merge the subsamples in a single file
hadd ${FILENAME}.root BUNCH??/Output/${FILENAME}.root

