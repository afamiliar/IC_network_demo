function [IC] = import_from_mvpa_toolbox(subject,pattern,wholebrain,ROI,regressor,folds,selector,radius)

% Takes objects from the Princeton MVPA toolbox and gets them ready for
% Informational Connectivity (IC) [IC Toolbox Version 1.1]
%
% [IC] =
% IMPORT_FROM_MVPA_TOOLBOX(SUBJECT,PATTERN,WHOLEBRAIN,ROI,REGRESSOR,FOLDS,SELECTOR,RADIUS)
%
% Adds a structure containing inputs for IC. The output will later be used
% for running the Informational Connectivity analysis.
%
% This script populates the IC structure with objects from a subject
% structure of the Princeton MVPA toolbox
% (https://code.google.com/p/princeton-mvpa-toolbox/). If you already use
% the MVPA toolbox, this script will prepare you data for IC. If you don't
% already have your data and searchlight indices within Matlab, check out
% the MVPA toolbox and then run this script.
% 
% SUBJECT: The MVPA toolbox subject (e.g., subj).
%
% DATA: The MVPA toolbox pattern containing your data. This pattern should
% have been added to your MVPA toolbox subject structure through the
% WHOLEBRAIN mask. Surround the name with apostrophes.
%
% WHOLEBRAIN: An MVPA toolbox wholebrain mask. It doesn't have to be a
% wholebrain mask if you are only examining (for example) just one lobe or
% hemisphere. Surround the name with apostrophes.
%
% ROI: If you wish to perform exploratory IC based on a seed's connections
% to searchlights, you should list the MVPA toolbox mask of your seed
% region. (This will be converted to a 2-dimensional matrix and then added
% to the IC structure.) Surround the name with apostrophes. If you wish to
% measure IC between pre-defined ROIs, list the relevant mask-names in a
% cell structure: e.g., {'ATL','PFC','FFA') if you wanted to extract IC
% between these three regions.
%
% REGRESSOR: An MVPA toolbox regressor containing the conditions you want
% to analyze. The regressor should already have been shifted to account for
% the hemodynamic delay. Surround the name with apostrophes.
%
% FOLDS: A 1-by-N MVPA toolbox selector of N trials, each numbered by fold.
% This could be a selector labelled by run. This vector determines which
% trials act as testing (fold M for iteration M) and training (all-but-M
% folds for iteration M) when calculating discriminability. Surround the
% name with apostrophes.
%
% SELECTOR: A 1-by-N MVPA toolbox selector showing which of N trials
% should be included (1) or excluded (0). At a minimum, rest should be
% excluded (unless rest is one of your classified conditions). See the
% MVPA toolbox script 'create_no_rest' as one solution. Remember to base
% this selector on the same regressor given above (i.e., AFTER shifting for
% the hemodynamic delay). Surround the name with apostrophes.
%
% RADIUS (optional): If you are running an exploratory IC analysis and
% haven't already used the toolbox script 'create_adj_list', this script
% will automatically construct searchlights across your WHOLEBRAIN mask,
% using the given radius (e.g., 2).
%
%
% Example for ROI analyses:
% [IC]=import_from_mvpa_toolbox(subj,'my_data','brain_mask',{'ATL','PFC','FFA'},'conditions_convt','runs','conditions_convt_norest')
%
% Example for exploratory analysis:
% [IC]=import_from_mvpa_toolbox(subj,'my_data','brain_mask','ATL_ssed','conditions_convt','runs','conditions_convt_norest',2)
%
%
% ======================================================================
% This script is part of the freely available Informational Connectivity
% Toolbox. See http://www.informationalconnectivity.org for more
% information.
%
% Authored by Marc Coutanche, University of Pennsylvania.
%
% No responsibility is taken for any problems that may be related to the
% use of this script.
% ======================================================================

if nargin<7
    error('You are missing inputs - have you listed at least seven?')
elseif (nargin==7)&&(isfield(subject,'adj_sphere')==0)&&(ischar(ROI)==1)
    error('Since you have only included one ROI, I am guessing you wish to run an exploratory searchlight IC analysis. In this case, you either need to specify a searchlight radius or have already run ''create_adj_list'' in the Princeton MVPA toolbox')
elseif nargin>7
    if isnumeric(radius)==0
        error('The searchlight radius must be a number (with no apostrophes)...')
    elseif isfield(subject,'adj_sphere')==0
        subject.adj_sphere = create_adj_list(subject,wholebrain,'radius',radius);
    elseif isfield(subject,'adj_sphere')==1
        warning('You have listed a radius AND previously created searchlight indices for your subject. We''ll use your existing searchlight indices. Remove them and then rerun if you''d prefer not to.')
    end
end

if ischar(ROI)==1
if count(get_mat(subject,'mask',ROI))<1
    error('Your seed mask doesn''t seem to contain any voxels - each seed voxel should be labelled with 1')
end
elseif ischar(ROI)==0
if length(ROI)<2
    error('You need to list at least two ROIs for the ROI IC analysis. If you wish to run an exploratory analysis with a seed, you should be listing a string with the seed region''s mask name')
end
for a=1:length(ROI)
    if count(get_mat(subject,'mask',ROI{a}))<1
    error(['One of your ROIs (' ROI{a} ') doesn''t contain any voxels'])
    end
end
end

IC.data=get_mat(subject,'pattern',pattern);
IC.wholebrain=get_mat(subject,'mask',wholebrain);
IC.conditions=get_mat(subject,'regressors',regressor);
IC.folds=get_mat(subject,'selector',folds);
IC.selector=get_mat(subject,'selector',selector);

if ischar(ROI)==1
seed_pattern_indices=convert_mask_to_pat_idx(get_mat(subject,'mask',wholebrain),get_mat(subject,'mask',ROI)==1);
seed_pattern=zeros(count(get_mat(subject,'mask',wholebrain)),1);
seed_pattern(seed_pattern_indices)=1;
IC.seed=seed_pattern;
IC.sl_indices=subject.adj_sphere;
else
% Let's label the ROIs by number:    
IC.ROIs=zeros(count(get_mat(subject,'mask',wholebrain)),1);
IC.ROI_names=ROI;
for nROIs=1:length(ROI)
clear ROI_pattern_indices
ROI_pattern_indices=convert_mask_to_pat_idx(get_mat(subject,'mask',wholebrain),get_mat(subject,'mask',ROI{nROIs})==1);
IC.ROIs(ROI_pattern_indices)=nROIs;
end
end

if size(IC.data,2)~=size(IC.conditions,2)
    error('Your data and condition matrices need to have the same numbers of trials')
elseif size(IC.data,2)~=size(IC.selector,2)
    error('Your data and selection matrices need to have the same numbers of trials')
elseif size(IC.selector,2)~=size(IC.conditions,2)
    error('Your condition and selection matrices need to have the same numbers of trials')
elseif size(IC.folds,2)~=size(IC.conditions,2)
    error('Your folds and condition matrices need to have the same numbers of trials')
elseif size(IC.folds,2)~=size(IC.selector,2)
    error('Your folds and selection matrices need to have the same numbers of trials')
elseif size(IC.folds,2)~=size(IC.data,2)
    error('Your folds and data matrices need to have the same numbers of trials')
elseif size(IC.conditions,1)<2
    error('You need to have at least 2 conditions')
elseif length(unique(IC.folds))<2
    error('You need to have at least 2 folds')
end

if ischar(ROI)==1
if size(IC.data,1)~=size(IC.sl_indices,1)
    error('Your data and searchlight indices need to have the same number of rows (voxels)')
elseif size(IC.seed,1)~=size(IC.sl_indices,1)
    error('Your seed and searchlight indices need to have the same number of rows (voxels)')
elseif size(IC.data,1)~=size(IC.seed,1)
    warning('Your data and seed matrices have differing numbers of rows (voxels). You won''t be able to run IC unless this is fixed. If you''re having problems, try adding the seed as a new MVPA toolbox pattern (rather than mask) and adding this into the IC structure (IC.seed = get_mat...)')
end
elseif size(IC.data,1)~=size(IC.ROIs,1)
    warning('Your data and ROIs matrices have differing numbers of rows (voxels). You won''t be able to run IC unless this is fixed. If you''re having problems, try adding the ROIs as an MVPA toolbox pattern (rather than masks) with a unique number for each region, and then adding this into the IC structure (IC.ROIs = get_mat...)')   
end
end
