#!/bin/bash

if [ $# -gt 2 ]; then
  echo "usage: batchRunStatUncertain basedirectory productiondirectory"
  exit 1
fi

if [ $# -lt 2 ]; then
  echo "usage: batchRunStatUncertain basedirectory productiondirectory"
  exit 1
fi

# incorporate shell utilities
utilities_file=/lustre/alice/users/$USER/CLUSTERMODELWAC/Clusters/GSI/utilities.sh
[ ! -f "$utilities_file" ] && { echo "Error: $utilities_file not found." >&2; exit 1; }
. "$utilities_file"

BASEDIRECTORY=$1
PRODUCTIONDIRECTORY=$2

if [ ! -d $BASEDIRECTORY/$PRODUCTIONDIRECTORY ]
then
  echo Production directory $PRODUCTIONDIRECTORY does NOT exist. WE CANNOT PROCEED!!!!
  exit 1	
fi

# the configuration file
CONFIGURATIONFILE=$BASEDIRECTORY/$PRODUCTIONDIRECTORY/configuration.json

# Extract the rapidities/pseudorapidities
if ! extract_json_scaled_num_array "$CONFIGURATIONFILE" "abs_y" "CRAPIDITIES"; then
    echo "Failed to extract 'abs_y'" >&2
    exit 1
fi
# and convert them to array
read -ra CRAPIDITIES_ARRAY <<< "$CRAPIDITIES"
# Extract the pT ranges
if ! extract_json_scaled_num_array "$CONFIGURATIONFILE" "ptRangeLows" "CPTRANGELOWS"; then
    echo "Failed to extract 'ptRangeLows'" >&2
    exit 1
fi
if ! extract_json_scaled_num_array "$CONFIGURATIONFILE" "ptRangeUps" "CPTRANGEUPS"; then
    echo "Failed to extract 'ptRangeUps'" >&2
    exit 1
fi
# and convert them to arrays
read -ra CPTRANGELOWS_ARRAY <<< "$CPTRANGELOWS"
read -ra CPTRANGEUPS_ARRAY <<< "$CPTRANGEUPS"

# extract needed information
EVENTFILTERS=(`sed -n '/"teventfilter"\s*:\s*\[\(.*\)\],/p' configuration.json | sed 's/\s*"teventfilter"\s*:\s*\[\(.*\)\],/\1/' | tr ',' ' '`)
NEVENTFILTERS=${#EVENTFILTERS[@]}

# the results production tag
PRODUCTIONTIME=`date +%Y%m%d_%H%M%S`

# submit the extraction of results with statistical uncertainties
# one job per rapidity and pT range (parallel) and per multiplicity class (sequential)
for (( i=0; i<${#CRAPIDITIES_ARRAY[@]}; i++ ))
do
  for ((iPtRange=0; iPtRange<${#CPTRANGELOWS_ARRAY[@]}; iPtRange++))
  do
    PRODUCTIONTAG=`printf "%s_Rap%03d_Pt%02d%02d" ${PRODUCTIONTIME} ${CRAPIDITIES_ARRAY[i]} ${CPTRANGELOWS_ARRAY[iPtRange]} ${CPTRANGEUPS_ARRAY[iPtRange]}`
    echo Submitting for production tag ${PRODUCTIONTAG}
    JOBID=
    for (( j=0; j<NEVENTFILTERS; j++ ))
    do
      if [ ${j} -eq 0 ]
      then
        DEPENDENCY=
      else
        DEPENDENCY="-d afterany:${JOBID}"
      fi
      echo Submitting multiplicity class ${EVENTFILTERS[j]}
      JOBNAME=`printf "waitStatsUncertain_%03d_%02d_%03d" ${i} ${iPtRange} ${j}`
      cmd="sbatch -J ${JOBNAME} ${DEPENDENCY} --chdir=${BASEDIRECTORY}/${PRODUCTIONDIRECTORY} --mem-per-cpu=8000 --time=07:00:00 -o ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/merge/${JOBNAME}Job.out -e ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/merge/${JOBNAME}Job.err /lustre/alice/users/${USER}/CLUSTERMODELWAC/Clusters/GSI/runScriptInSingularity.sh /lustre/alice/users/${USER}/CLUSTERMODELWAC/Clusters/GSI/runStatsUncertain.sh ${PRODUCTIONTAG} ${i} ${iPtRange} ${j}"
      JOBID=($(eval $cmd | tee /dev/tty | awk '{print $4}'))
      echo $cmd >> ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/submit.log
      echo "" >> ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/submit.log
    done
  done
done


