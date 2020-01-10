-- change atlas to one with smaller ROIs?

# README

This repo provides MATLAB code for using informational connectivity to construct whole-brain networks with functional MRI data.

Tested on macOS 10.13.6 with MATLAB R2015B and R2019A. 

Required software:
* MATLAB

Required toolboxes (provided in repo):
* CoSMoMVPA (http://github.com/CoSMoMVPA/CoSMoMVPA)
* IC toolbox (http://informationalconnectivity.org)
* Brain Connectivity toolbox (https://sites.google.com/site/bctnet/)

## Usage

The main demo is provided in **"analyze_IC_brainnetome.m"**

Details on how to set up input for IC toolbox can be found in **run_ROI_IC.m** in **toolboxes/IC_toolbox/**. **create_** scripts within directories show how inputs were created for the demo.

## Data and timing information

Data file **niftiDATA_Subject001.nii.gz** in **data/** contains functional MRI images collected while one subject viewed images of 9 different faces. Images were presented one at a time, in a pseudo-randomized order (2 secs per image, 2 secs TR), for 5 scan runs (200 TRs per run, 1000 TRs total), with blank trials intermixed. The subject's task was to press a button with each hand simultaneously when an image repeated (same face presented twice in a row). Blank and repeat trials are excluded in the IC analysis.

Text files in **timing/timing_files/** indicate the stimulus event order. There is a separate text file for each run and trial type, with binary values (0=not on; 1=on). Mat files in **timing/** have timing information in format required by the IC toolbox, **create_IC_timing_files.m** shows how the mat files were generated based on the text files.

Data was pre-processed using the CONN toolbox (https://web.conn-toolbox.org/), which included:
* motion correction
* registration to structural
* distortion correction using fieldmap B0 image
* slice-timing correction
* outlier detection
* registration to MNI space
* spatial smoothing (8mm FWHM)

Data was denoised using the CONN toolbox, which included:
* regression to remove:
  * top 5 principal components of white matter and CSF voxels (aCompCor)
  * outliers
  * linear trends
  * motion parameters
* band-pass filtering

## Atlas/ROI information












