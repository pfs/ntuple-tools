#!/bin/bash

INFILE=${1}
OUTFILE=${2}

echo "Input for this job is ${INFILE} with output @ ${OUTFILE}"
cd /afs/cern.ch/user/p/psilva/work/HGCal/Electronics/CMSSW_9_3_2/src
eval `scram r -sh`
cd -
outDir=`dirname ${OUTFILE}`
localOut=`basename ${OUTFILE}`
mkdir -p ${outDir}
echo "cmsRun $CMSSW_BASE/src/RecoNtuples/HGCalAnalysis/test/roiAnalysisConfig.py maxEvents=-1 inputFiles=file:${INFILE} outputFile=${localOut}";
cmsRun $CMSSW_BASE/src/RecoNtuples/HGCalAnalysis/test/roiAnalysisConfig.py maxEvents=-1 inputFiles=${INFILE} outputFile=${localOut};
mv -v ${localOut} ${OUTFILE}
