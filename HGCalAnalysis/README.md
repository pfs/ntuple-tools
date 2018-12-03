# Simulations

Check test/sim directory for different generation scripts as well as DIG and RECO.
Condor files are available to submit the jobs on the batch
The integrity of the files can be checked with the following script.
Zombie or files with 0 events are removed.

```
python scripts/checkIntegrity.py input_dir
```

# Region Of Interest analyzer

The ROIanalyzer is currently set to analyze H->gg events but can be easily changed for other topologies.
To create a summary ntuple including the ROI with 3 different SR sizes in HGCAL around the generator level particles
and additional 5 control regions for noise/pileup from a rotation in phi at the same eta run:

```
cmsRun test/roiAnalysisConfig.py maxEvents=-1 inputFiles=/store/cmst3/user/psilva/HGCal/H125gg_EE/CMSSW_9_3_2/DIGI_PU0_0p0/RECO/ outputFile=ROISummary_PU0_0p0.root;
```

A second script is used to collect the energy in each ROI and in the associated noise control regions.
It produces a small ntuple with the basic inputs for calibration.

```
python scripts/summarizeROIforCalibration.py ROISummary_PU0_0p0.root ROISummary_PU0_0p0_forcalib.root 
```

The calibration can be run as follows (first argument is the no pileup file, second argument is the pileup file, last argument is a tag for the calibration).
The L0 (relative - eta), L1 (absolute - E), L2 (pileup - average noise) calibration is stored in local pickle file in a dict.

```
python scripts/runROICalibration.py ROISummary_PU0_0p0_forcalib.root ROISummary_PU140_0p0_forcalib.root 140
```

# Occupancy analyzer

To create the ntuples onw can use the following script.
The underlying directory is hardcoded to my CMST3 directory (to be changed if needed)
and tag_name is the name of the sub-directory with DIGIs to process

```
do cmsRun ../test/occAnalysisConfig.py tag=tag_name
```

Once the histogram file is available one can summarize the profile of a given distribution using the following script

```
python  scripts/createOccProfileSummaries.py --doBusy --var var_name -i occ_analysis_tag_name.root -o occ_plots;
```

The profiles can be furthermore compared with

```
python scripts/plotOccProfiles.py -o occ_plots/tag1vstag2 -i occ_plots/occ_summary_tag_name_1.root:"tag1",summary/occ_summary_tag_name_2.root:"tag2"
```
 