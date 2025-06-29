#!/bin/bash

if [ $# -lt 2 ]; then
  echo "usage: runSamplerAndSmash args..."
  exit 1
fi

CLUSTERMODELWACDIR=/lustre/alice/users/vgonzale/CLUSTERMODELWAC


# setting the proper environment
export ALIEN_SITE=GSI
NEEDED=VO_ALICE@O2Physics::daily-20250614-0000-1
source <( /cvmfs/alice.cern.ch/bin/alienv printenv $NEEDED)
echo $NEEDED

####################################################################################################
echo "Setting up WAC"
####################################################################################################
export WAC_ROOT=$CLUSTERMODELWACDIR
export WAC_SOURCE="$WAC_ROOT"
export WAC_BIN="$WAC_ROOT/bin"
export WAC_LIB="$WAC_ROOT/lib"
export PATH="$WAC_BIN:$PATH"
export DYLD_LIBRARY_PATH="$WAC_LIB:$DYLD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$WAC_LIB:$LD_LIBRARY_PATH"

echo "$@"
$CLUSTERMODELWACDIR/bin/RunBestSimulationTwoParticlesDiff "$@"

