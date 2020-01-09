-- remove need for COSMO (change nii to mat)
-- add brain connectivity toolbox
-- change atlas to one with smaller ROIs?

# README

This repo provides MATLAB code for using informational connectivity to construct whole-brain networks.

Tested on macOS 10.13.6 with MATLAB R2015B and macOS ##.##.# with MATLAB R2019A. 

Required software:
* MATLAB

Required toolboxes (provided in repo):
* CoSMoMVPA (http://github.com/CoSMoMVPA/CoSMoMVPA)
* IC toolbox (http://informationalconnectivity.org)
* Brain Connectivity toolbox (https://sites.google.com/site/bctnet/)

## Data and timing information

Data file in **data/** is functional MRI images collected while one subject viewed images of 9 different faces. Images were presented one at a time at the center of a screen, in a pseudo-randomized order (2 secs per image, 2 secs TR), for 5 runs (200 TRs per run, 1000 TRs total), with blank trials intermixed. The subject's task was to press a button with each hand simultaneously when an image repeated (same face presented twice in a row), these trials are excluded in the IC analysis.

Text files in **timing/timing_files/** indicate the event order. There is a separate text file for each run and trial type, with binary values (0=not on; 1=on). Mat files in **timing/** have timing information in format required by the IC toolbox, **create_IC_timing_files.m** shows how the mat files were generated based on the text files.

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


```bash
ds = cosmo_fmri_dataset(data_path);
data = ds.samples';
```

## ROI information



## Usage









