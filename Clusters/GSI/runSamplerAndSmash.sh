#!/bin/bash

if [ $# -lt 2 ]; then
  echo "usage: runSamplerAndSmash args..."
  exit 1
fi

BESTDIRECTORY=/lustre/alice/users/vgonzale/BESTWRAPPER


# setting the proper environment
NEEDED="VO_ALICE@CMake::v3.31.6-1"

source <( /cvmfs/alice.cern.ch/bin/alienv printenv $NEEDED)
echo $NEEDED

echo "$@"
$BESTDIRECTORY/build/sampler_and_smash "$@"

