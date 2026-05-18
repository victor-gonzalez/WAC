#!/bin/bash

if [ $# -gt 6 ]; then
  echo "usage: runStatsUncertain pythia/best productiontag ixrap ixptrange ixevtflt [opt]"
  exit 1
fi

if [ $# -lt 5 ]; then
  echo "usage: runStatsUncertain pythia/best productiontag ixrap ixptrange ixevtflt [opt]"
  exit 1
fi

WHICHGEN=$1
PRODUCTIONTAG=$2
IXRAP=$3
IXPTRANGE=$4
IXEVTFLT=$5
# optional statUncertain option token; "det" selects the detector-effects
# (reconstructed) pass, empty selects the raw pass (unchanged behaviour)
OPT=${6:-}

# no more core files
ulimit -c 0
echo ulimit `ulimit -c`

# setting the proper environment scenario
NEEDED="VO_ALICE@O2Physics::daily-20250614-0000-1"

source <( /cvmfs/alice.cern.ch/bin/alienv printenv $NEEDED)
echo $NEEDED

####################################################################################################
echo "Setting up WAC"
####################################################################################################
export WAC_ROOT=/lustre/alice/users/$USER/CLUSTERMODELWAC
export WAC_SOURCE="$WAC_ROOT"
export WAC_BIN="$WAC_ROOT/bin"
export WAC_LIB="$WAC_ROOT/lib"

export PATH="$WAC_BIN:$PATH"

if [[ ${WHICHGEN} == "pythia" ]]
then

  export PYTHIA8=/cvmfs/alice.cern.ch/el9-x86_64/Packages/pythia/v8315-alice1-1

  export DYLD_LIBRARY_PATH="$WAC_LIB:$PYTHIA8/lib:$DYLD_LIBRARY_PATH"
  export LD_LIBRARY_PATH="$WAC_LIB:$PYTHIA8/lib:$LD_LIBRARY_PATH"

  statUncertainPythia $PRODUCTIONTAG "$OPT" $IXRAP $IXPTRANGE $IXEVTFLT
else 
  if [[ ${WHICHGEN} == "best" ]]
  then
    export DYLD_LIBRARY_PATH="$WAC_LIB:$DYLD_LIBRARY_PATH"
    export LD_LIBRARY_PATH="$WAC_LIB:$LD_LIBRARY_PATH"

    statUncertainBest $PRODUCTIONTAG "$OPT" $IXRAP $IXPTRANGE $IXEVTFLT
  else
    echo "ERROR: UNKNOWN GENERATOR ${WHICHGEN}"
    exit 1
  fi
fi
