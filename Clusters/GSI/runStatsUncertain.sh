#!/bin/bash

if [ $# -gt 4 ]; then
  echo "usage: runStatsUncertain productiontag ixrap ixptrange ixevtflt"
  exit 1
fi

if [ $# -lt 4 ]; then
  echo "usage: runStatsUncertain productiontag ixrap ixptrange ixevtflt"
  exit 1
fi

PRODUCTIONTAG=$1
IXRAP=$2
IXPTRANGE=$3
IXEVTFLT=$4

# no more core files
ulimit -c 0
echo ulimit `ulimit -c`

# setting the root and pythia scenario
export ALIEN_SITE=GSI
LATEST="VO_ALICE@AliGenerators::v20241101-1"
export ALIPHYSICS_VERSION=$LATEST

source <( /cvmfs/alice.cern.ch/bin/alienv printenv $LATEST)
echo $LATEST

export PYTHIA8=/cvmfs/alice.cern.ch/el7-x86_64/Packages/pythia/v8311-18

####################################################################################################
echo "Setting up WAC"
####################################################################################################
export WAC_ROOT=/lustre/alice/users/$USER/CLUSTERMODELWAC
export WAC_SOURCE="$WAC_ROOT"
export WAC_BIN="$WAC_ROOT/bin"
export WAC_LIB="$WAC_ROOT/lib"

export PATH="$WAC_BIN:$PATH"
export DYLD_LIBRARY_PATH="$WAC_LIB:$PYTHIA8/lib:$DYLD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$WAC_LIB:$PYTHIA8/lib:$LD_LIBRARY_PATH"

statUncertainPythia $PRODUCTIONTAG "" $IXRAP $IXPTRANGE $IXEVTFLT
