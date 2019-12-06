%% create ROI files for IC analysis
%   amf
%   Nov 2019
%
%   This prepares ROI data for input to the IC toolbox.
%   Toolbox requires each ROI to have a unique id number,and input to be a 
%   vector in which each element corresponds to the id number for that
%   voxel. All other voxels should be set to 0.
%
%   -- loads in Brainnetome atlas (Fan et al., 2016) data using cosmo-mvpa
%   toolbox (Oosterhof & Connolly)
%        atlas:    https://atlas.brainnetome.org/
%        toolbox:  http://www.cosmomvpa.org/
%   -- separates into left and right hemispheres based on ROI numbers
%        LH = odd ;  RH = even
%        ROI numbers are specific to the Brainnetome atlas
%   -- saves mat files
%
%%

% make sure cosmo-mvpa toolbox is on matlab path
addpath(genpath('/Users/afam/Documents/MATLAB/CoSMoMVPA-master'));

atlas_fn = 'BN_Atlas_246_2mm.nii.gz';


% load data
roi_data = cosmo_fmri_dataset(atlas_fn);
ROIs = roi_data.samples';

% remove sub-cortical ROIs
ROIs(ROIs>210) = 0;

% find indices where ROI numbers are odd/even
evenIndices = rem(ROIs, 2) == 0;
oddIndices  = rem(ROIs, 2) == 1;

% make left/right hemi vectors
LH_ROIs = ROIs.*oddIndices;
RH_ROIs = ROIs.*evenIndices;

% save to mat files
save('LH_ROIs.mat','LH_ROIs')
save('RH_ROIs.mat','RH_ROIs')

