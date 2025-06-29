#!/bin/bash

if [ $# -lt 1 ]; then
  echo "usage: runScriptInSingularity scriptname scriptargs..."
  exit 1
fi

SCRIPTNAME=$1

shift

singularity_image=/cvmfs/alice.cern.ch/containers/fs/singularity/el9
singularity exec -B /cvmfs,/lustre ${singularity_image} ${SCRIPTNAME} "$@"
