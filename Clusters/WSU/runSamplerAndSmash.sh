#!/bin/bash

if [ $# -lt 2 ]; then
  echo "usage: runSamplerAndSmash args..."
  exit 1
fi

BESTDIRECTORY=/wsu/home/gm/gm62/gm6269/PROJECTS/BEST


# setting the proper environment
module load gnu7
module load eigen
module load cmake/3.21.1
module load gsl
module load boost

echo "$@"
$BESTDIRECTORY/build/sampler_and_smash "$@"

