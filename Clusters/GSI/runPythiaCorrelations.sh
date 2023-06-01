#!/bin/bash

MIXEDEVENTS=$1

# no more core files
ulimit -c 0
echo ulimit `ulimit -c`

# setting the root and pythia scenario
export ALIEN_SITE=GSI
LATEST="VO_ALICE@AliGenerators::v20221125-1"
export ALIPHYSICS_VERSION=$LATEST

source <( /cvmfs/alice.cern.ch/bin/alienv printenv VO_ALICE@Python::v3.9.12-10)
source <( /cvmfs/alice.cern.ch/bin/alienv printenv $LATEST)
echo $LATEST

export PYTHIA8=/cvmfs/alice.cern.ch/el7-x86_64/Packages/pythia/v8304-44

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

TASKIX=$SLURM_ARRAY_TASK_ID
SEED=$(( (SLURM_ARRAY_TASK_ID + SLURM_ARRAY_JOB_ID*1000) % 900000000 ))
echo "The seed is $SEED"

# Execute application code
if [ $MIXEDEVENTS == "YES" ]
then
  RunPythiaSimulationTwoParticlesDiffME $TASKIX $SEED
else
  RunPythiaSimulationTwoParticlesDiff $TASKIX $SEED
fi

