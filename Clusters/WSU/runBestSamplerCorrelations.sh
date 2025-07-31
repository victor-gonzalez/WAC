#!/bin/bash
WORKINGDIRECTORY=$1
HYPERSURFACE=$2
NUMBEROFEVENTS=$3
ME=$4

CLUSTERMODELWAC=/wsu/home/gm/gm62/gm6269/PROJECTS/CLUSTERMODELWAC
BESTDIRECTORY=/wsu/home/gm/gm62/gm6269/PROJECTS/BEST
PIPEDIRECTORY=/tmp/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
OUTPUTDIRECTORY=$WORKINGDIRECTORY
CONFIGFILE=$BESTDIRECTORY/config.yaml

# create the intermediate directory for named pipe output
mkdir -p $PIPEDIRECTORY

# Pipe name
PIPENAME="particle_lists.oscar"
TEMPIPENAME=$PIPENAME.unfinished
# We create it with the temporal name
mkfifo $PIPEDIRECTORY/$TEMPIPENAME

# Pass base filename to original_app
source $CLUSTERMODELWAC/Clusters/WSU/runSamplerAndSmash.sh -c $CONFIGFILE -o "General: {Nevents: $NUMBEROFEVENTS}" -o "Output_Directory: $PIPEDIRECTORY" -o "Output_Kind: namedpipe" -o "HyperSurface: $HYPERSURFACE/surface.dat" &

# run the data collecting engine
$CLUSTERMODELWAC/Clusters/WSU/runBestDataCollectingEngine.sh $PIPEDIRECTORY/$TEMPIPENAME $NUMBEROFEVENTS $ME ${SLURM_ARRAY_TASK_ID}

# But we clean it with the definite name
rm $PIPEDIRECTORY/$PIPENAME

# Recover a single copy of the configuration file used by the sampler
if [ ! -f "$OUTPUTDIRECTORY/config.yaml" ] 
then 
  cp -n "$PIPEDIRECTORY/config.yaml" "$OUTPUTDIRECTORY/config.yaml"
fi

rm -rf $PIPEDIRECTORY
