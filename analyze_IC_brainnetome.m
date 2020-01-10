%% calculate informational connectivity between ROI pairs
%   amf
%   Nov 2019
%
%       - Loads nifti data file using cosmomvpa toolbox
%       - Calculates IC for regions in user-defined hemispheres 
%           using IC toolbox
%       - Saves result to mat file
%
%% user defined

hemispheres = {'LH','RH'}; % 'LH' and/or 'RH'

%% add current directory and sub-dir's to path
addpath(genpath(pwd))

%% load ROI data
disp('STATUS: loading atlas ROIs')
load('brainnetome_ROI/LH_ROIs.mat');
load('brainnetome_ROI/RH_ROIs.mat');
load('brainnetome_ROI/ROI_names.mat');

%% load timing info
disp('STATUS: loading timing info')
load('timing/selector.mat');
load('timing/conditions.mat');
load('timing/folds.mat');

%% run IC
IC = [];
IC.ROI_names  = ROI_names;
IC.conditions = conditions;
IC.selector   = selector;
IC.folds      = folds;

% load subject data
disp('STATUS: loading subject data')
data_file = 'niftiDATA_Subject001.nii.gz';
data_path = fullfile('data/',data_file);

ds = cosmo_fmri_dataset(data_path);
data = ds.samples';

% update IC structure
IC.data = data;

% run IC for each hemisphere separately
for h = 1:size(hemispheres,2)
	if strcmp(hemispheres{h},'LH')
        ROIs = LH_ROIs;
    elseif strcmp(hemispheres{h},'RH')
        ROIs = RH_ROIs;
	end

    % update IC structure
    IC.ROIs = [];
    IC.ROIs = ROIs;

    % run IC
    disp(['STATUS: running IC in ' hemispheres{h}])
    [~, IC_matrix] = run_ROI_IC(IC);

    % save results to mat file
    disp('STATUS: saving results')
    save(['IC_results/sub001_', hemispheres{h} '.mat'],'IC_matrix')
end

