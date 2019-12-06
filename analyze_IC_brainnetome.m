%% calculate informational connectivity between ROI pairs
%   amf
%   Nov 2019
%
%% user defined

hemispheres = {'LH','RH'}; % 'LH' and/or 'RH'

%% make sure cosmo-mvpa and IC toolboxes are on matlab path
addpath(genpath('/Users/afam/Documents/MATLAB/CoSMoMVPA-master'));
addpath(genpath('IC_toolbox'))

%% load ROI data
load('brainnetome_ROI/LH_ROIs.mat');
load('brainnetome_ROI/RH_ROIs.mat');
load('brainnetome_ROI/ROI_names.mat');

%% load timing info
load('timing/selector.mat');
load('timing/conditions.mat');
load('timing/folds.mat');

%% run IC
numSubjs = 1;

IC = [];
IC.ROI_names = ROI_names;
IC.conditions = conditions;
IC.selector = selector;
IC.folds = folds;

for s = 1:numSubjs
     disp(['subject ' int2str(s)])
     % load subject data
     if s < 10
         data_file = ['niftiDATA_Subject00' int2str(s) '_Condition000.nii'];
     elseif s < 100
         data_file = ['niftiDATA_Subject0' int2str(s) '_Condition000.nii'];
     end
     data_path = fullfile('data/',data_file);

     ds = cosmo_fmri_dataset(data_path);
     data = ds.samples';
     
     % update IC structure
     IC.data = [];
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
         [~, IC_matrix] = run_ROI_IC(IC);

         % save results
         save(['IC_results/subj' int2str(s),'_', hemispheres{h} '.mat'],'IC_matrix')
     end
     
end

