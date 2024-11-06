#!/bin/bash
export ALIEN_SITE=GSI
source <( /cvmfs/alice.cern.ch/bin/alienv printenv VO_ALICE@AliGenerators::v20241101-1)
source <( /cvmfs/alice.cern.ch/bin/alienv printenv VO_ALICE@CMake::v3.28.1-20)
source set-WAC-GSI

