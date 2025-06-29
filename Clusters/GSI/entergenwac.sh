#!/bin/bash
export ALIEN_SITE=GSI
source <( /cvmfs/alice.cern.ch/bin/alienv printenv VO_ALICE@O2Physics::daily-20250614-0000-1)
source <( /cvmfs/alice.cern.ch/bin/alienv printenv VO_ALICE@CMake::v3.31.6-1)
source set-WAC-GSI

