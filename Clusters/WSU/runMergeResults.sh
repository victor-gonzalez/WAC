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
module load root

# Execute application code
root -l -b -q "/wsu/home/gm/gm62/gm6269/PROJECTS/CLUSTERMODELWAC/TOOLS/recursiveResultFilesMerge.C(\"${MERGEDFNAME}\",\"${PATTERNNAME}\")"
