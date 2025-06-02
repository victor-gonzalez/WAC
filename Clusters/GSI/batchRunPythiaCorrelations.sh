#!/bin/bash

if [ $# -gt 4 ]; then
  echo "usage: batchRunPythiaCorrelations basedirectory nmainjobs nsubjobs {me}"
  exit 1
fi

if [ $# -lt 3 ]; then
  echo "usage: batchRunPythiaCorrelations basedirectory nmainjobs nsubjobs {me}"
  exit 1
fi

# incorporate shell utilities
utilities_file=/lustre/alice/users/$USER/CLUSTERMODELWAC/Clusters/GSI/utilities.sh
[ ! -f "$utilities_file" ] && { echo "Error: $utilities_file not found." >&2; exit 1; }
. "$utilities_file"

BASEDIRECTORY=$1
NMAINJOBS=$2
NSUBJOBS=$3
MIXEDEVENTS=${4:-NO}
if [ $MIXEDEVENTS != "me" ] && [ $MIXEDEVENTS != "NO" ]
then
  echo "usage: batchRunPythiaCorrelations basedirectory nmainjobs nsubjobs {me}"
  exit 1
fi
if [ $MIXEDEVENTS == "me" ]
then
  MIXEDEVENTS=YES
fi

PRODUCTIONDIRECTORY=OUT`date +%Y%m%d%H%M`
ALICESWREVISION="PENDING"

if [ -d $BASEDIRECTORY/$PRODUCTIONDIRECTORY ]
then
  echo Production directory $PRODUCTIONDIRECTORY does exist. WAIT AT LEAST ONE MINUTE!!!!
  exit 1	
fi

mkdir $BASEDIRECTORY/$PRODUCTIONDIRECTORY
echo Production $PRODUCTIONDIRECTORY

mkdir $BASEDIRECTORY/$PRODUCTIONDIRECTORY/log
mkdir $BASEDIRECTORY/$PRODUCTIONDIRECTORY/log/merge

# let's preserve the configuration
CONFIGURATIONFILE=/lustre/alice/users/$USER/CLUSTERMODELWAC/Clusters/GSI/configuration.json
cp $CONFIGURATIONFILE $BASEDIRECTORY/$PRODUCTIONDIRECTORY

# and extract needed information
OUTFNAME=`sed -n '/"outputfname"\s*:\s*"\(.*\)",/p' ${CONFIGURATIONFILE} | sed 's/\s*"outputfname"\s*:\s*"\(.*\)",/\1/'`
TASKNAME=`sed -n '/"taskname"\s*:\s*"\(.*\)",/p' ${CONFIGURATIONFILE} | sed 's/\s*"taskname"\s*:\s*"\(.*\)",/\1/'`

# Extract the rapidities/pseudorapidities
if ! extract_json_scaled_num_array "$CONFIGURATIONFILE" "abs_y" "CRAPIDITIES"; then
    echo "Failed to extract 'abs_y'" >&2
    exit 1
fi
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

MERGEJOBSIDS=

# submit the main jobs
for ijob in $(seq 1 $NMAINJOBS)
do
  # launching the main generation jobs
  WORKINGDIRECTORY=$BASEDIRECTORY/$PRODUCTIONDIRECTORY/$(printf "BUNCH%02d" $ijob)
  mkdir -p $WORKINGDIRECTORY/Output
  mkdir -p $WORKINGDIRECTORY/log
  mkdir -p $WORKINGDIRECTORY/log/merge

  cp $CONFIGURATIONFILE $WORKINGDIRECTORY/

  # submit the job array
  cmd="sbatch -J batch__PythiaCorr --array=1-${NSUBJOBS} --chdir=${WORKINGDIRECTORY} --time=03:00:00 -o ${WORKINGDIRECTORY}/log/Job_%A_%a.out -e ${WORKINGDIRECTORY}/log/Job_%A_%a.err /lustre/alice/users/${USER}/CLUSTERMODELWAC/Clusters/GSI/runScriptInSingularity.sh /lustre/alice/users/${USER}/CLUSTERMODELWAC/Clusters/GSI/runPythiaCorrelations.sh $MIXEDEVENTS"
  ARRAYJOBID=($(eval $cmd | tee /dev/tty | awk '{print $4}'))
  echo $cmd >> ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/submit.log
  echo "" >> ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/submit.log

  # submit the merging of the job array results per rapidity
  for crap in ${CRAPIDITIES}
  do
    # and per pT range
    for ((iPtRange=0; iPtRange<${#CPTRANGELOWS_ARRAY[@]}; iPtRange++))
    do
      sleep 2s
      FNAMEPATERN=`printf ${OUTFNAME}_??? ${crap} ${CPTRANGELOWS_ARRAY[iPtRange]} ${CPTRANGEUPS_ARRAY[iPtRange]}`
      MERGEDFNAME=`printf ${OUTFNAME} ${crap} ${CPTRANGELOWS_ARRAY[iPtRange]} ${CPTRANGEUPS_ARRAY[iPtRange]}`

      # merge the results
      mergeJOBNAME="waitMerge_$(printf "BUNCH%02d_Rap%03d_Pt%02d%02d" $ijob ${crap} ${CPTRANGELOWS_ARRAY[iPtRange]} ${CPTRANGEUPS_ARRAY[iPtRange]})"
      cmd="sbatch -J $mergeJOBNAME --chdir=${WORKINGDIRECTORY}/Output --mem-per-cpu=8000 --time=03:00:00 -d afterany:$ARRAYJOBID -o ${WORKINGDIRECTORY}/log/merge/Job_%A.out -e ${WORKINGDIRECTORY}/log/merge/Job_%A.err /lustre/alice/users/${USER}/CLUSTERMODELWAC/Clusters/GSI/runScriptInSingularity.sh /lustre/alice/users/${USER}/CLUSTERMODELWAC/Clusters/GSI/runMergePythiaResults.sh ${MERGEDFNAME} ${FNAMEPATERN}"
      MERGEJOBID=($(eval $cmd | tee /dev/tty | awk '{print $4}'))
      echo $cmd >> ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/submit.log

      echo "" >> ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/submit.log
      MERGEJOBSIDS=${MERGEJOBSIDS}:${MERGEJOBID}
    done
  done
done

# submit the final merging for both pairs and singles results
for crap in ${CRAPIDITIES}
do
  for ((iPtRange=0; iPtRange<${#CPTRANGELOWS_ARRAY[@]}; iPtRange++))
  do
    sleep 2s
    MERGEDFNAME=`printf ${OUTFNAME} ${crap} ${CPTRANGELOWS_ARRAY[iPtRange]} ${CPTRANGEUPS_ARRAY[iPtRange]}`

    cmd="sbatch -J waitFinalMerge --chdir=${BASEDIRECTORY}/${PRODUCTIONDIRECTORY} --mem-per-cpu=8000 --time=03:00:00 -d afterany${MERGEJOBSIDS} -o ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/merge/SinglesMergeJob_%A.out -e ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/merge/SinglesMergeJob_%A.err /lustre/alice/users/${USER}/CLUSTERMODELWAC/Clusters/GSI/runScriptInSingularity.sh /lustre/alice/users/${USER}/CLUSTERMODELWAC/Clusters/GSI/runMergePythiaSubsamples.sh ${MERGEDFNAME}"
    JOBID=($(eval $cmd | tee /dev/tty | awk '{print $4}'))
    echo $cmd >> ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/submit.log

    echo "" >> ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/submit.log
  done
done
sleep 2s

# submit the extraction of results with statistical uncertainties
cmd="sbatch -J waitStatsUncertain --chdir=${BASEDIRECTORY}/${PRODUCTIONDIRECTORY} --time=01:00:00 -d afterany${MERGEJOBSIDS} -o ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/merge/WaitStatsUncertainJob.out -e ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/merge/WaitStatsUncertainJob.err /lustre/alice/users/${USER}/CLUSTERMODELWAC/Clusters/GSI/batchRunStatsUncertain.sh ${BASEDIRECTORY} ${PRODUCTIONDIRECTORY}"
JOBID=($(eval $cmd | tee /dev/tty | awk '{print $4}'))
echo $cmd >> ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/submit.log
echo "" >> ${BASEDIRECTORY}/${PRODUCTIONDIRECTORY}/log/submit.log


