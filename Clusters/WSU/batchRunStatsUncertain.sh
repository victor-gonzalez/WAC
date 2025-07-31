#!/bin/bash

if [ $# -gt 3 ]; then
  echo "usage: batchRunStatUncertain pythia/best basedirectory productiondirectory"
  exit 1
fi

if [ $# -lt 3 ]; then
  echo "usage: batchRunStatUncertain pythia/best basedirectory productiondirectory"
  exit 1
fi

# incorporate shell utilities
utilities_file=/wsu/home/gm/gm62/gm6269/PROJECTS/CLUSTERMODELWAC/Clusters/utilities.sh
[ ! -f "$utilities_file" ] && { echo "Error: $utilities_file not found." >&2; exit 1; }
. "$utilities_file"

WHICHGEN=$1
BASEDIRECTORY=$2
PRODUCTIONDIRECTORY=$3

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
      cmd="sbatch -J ${JOBNAME} ${DEPENDENCY} -q primary --chdir=${BASEDIRECTORY}/${PRODUCTIONDIRECTORY} --mem-per-cpu=8000 --time=07:00:00 -o ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/merge/${JOBNAME}Job.out -e ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/merge/${JOBNAME}Job.err /wsu/home/gm/gm62/gm6269/PROJECTS/CLUSTERMODELWAC/Clusters/GSI/runStatsUncertain.sh ${WHICHGEN} ${PRODUCTIONTAG} ${i} ${iPtRange} ${j}"
      JOBID=($(eval $cmd | tee /dev/tty | awk '{print $4}'))
      echo $cmd >> ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/submit.log
      echo "" >> ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/submit.log
    done
  done
done


