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
module load root

# merge the subsamples in a single file
hadd ${FILENAME}.root BUNCH??/Output/${FILENAME}.root

