function [IC] = run_exploratory_IC(IC);

% Runs Informational Connectivity (IC) between a seed and searchlights [IC
% Toolbox Version 1.1]
%
% [IC] = RUN_EXPLORATORY_IC(IC)
%
% Takes in an IC structure (details below), calculates the IC between the
% seed and each searchlight, and then outputs results to the IC structure.
% IC.results is a vector with V rows (one per voxel). Each value represents
% the magnitude of IC between the seed and the searchlight centered on that
% voxel. The script outputs its progress in 10% increments (10%, 20%, 30%,
% etc.)
%
% IC: If you already use the Princeton MVPA toolbox (or don''t already have
% your data in Matlab), see 'import_from_mvpa_toolbox' to create the IC
% structure. Alternatively, you can create the structure yourself, with the
% following necessary fields:
%   - IC.data: a 2-dimensional matrix with V rows (one per voxel) and N
%       columns (one per trial).
%   - IC.seed: a column vector with V rows (one per voxel) that should
%       match the rows in IC.data (i.e., voxel 42 in the seed is also voxel
%       42 in the data). 1 represents a seed voxel, while 0 indicates a
%       non-seed voxel.
%   - IC.conditions: a matrix with N columns and one row per condition.
%       Each trial has a zero if the condition is not active and one when
%       it is. Only one condition should be active at any trial (i.e., no
%       column sums to more than one). The conditions should already have
%       been shifted for the hemodynamic delay.
%   - IC.folds: a 1-by-N vector, with each trial numbered by fold. In many
%       cases, the folds will be runs. Each trial's IC discriminability
%       compares the trial's pattern to representative patterns for each
%       condition from the training set. This vector determines which
%       trials act as testing (fold M for iteration M) and training
%       (all-but-M folds for iteration M).
%   - IC.selector: a 1-by-N vector that indicates which trials should be
%       included (1) or excluded (0). At a minimum, rest should be excluded
%       (unless rest is one of your classified conditions). Remember to
%       base this selector on trials after shifting for the hemodynamic
%       delay.
%   - IC.sl_indices: A 2-dimensional matrix with V rows and columns equal
%       to the radius of the searchlights. The numbers in a row reflect the
%       indices of the other voxels found within the searchlight. This must
%       correspond to IC.data (i.e., voxel 42 in IC.sl_indices is voxel 42
%       in IC.data).
%
%
% Example: [IC] = run_exploratory_IC(IC)
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

% First, check over the IC structure:
if isfield(IC,'data')==0
    error('Your IC structure doesn''t have a data subfield. Have you named it correctly? (.data)')
elseif isfield(IC,'seed')==0
    error('Your IC structure doesn''t have a seed subfield. Have you named it correctly? (.seed)')
elseif isfield(IC,'conditions')==0
    error('Your IC structure doesn''t have a conditions subfield. Have you named it correctly? (.conditions)')
elseif isfield(IC,'folds')==0
    error('Your IC structure doesn''t have a folds subfield. Have you named it correctly? (.folds)')
elseif isfield(IC,'selector')==0
    error('Your IC structure doesn''t have a selector subfield. Have you named it correctly? (.selector)')
elseif isfield(IC,'sl_indices')==0
    error('Your IC structure doesn''t include your searchlight indices. Did you remember to create them or run ''import_from_mvpa_toolbox''?')
elseif size(IC.data,2)~=size(IC.conditions,2)
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
elseif size(IC.data,1)~=size(IC.sl_indices,1)
    error('Your data and searchlight indices need to have the same numbers of voxels')
elseif size(IC.conditions,1)<2
    error('You need to have at least 2 conditions')
elseif length(unique(IC.folds))<2
    error('You need to have at least 2 folds')
elseif size(IC.data,1)~=size(IC.seed,1)
    error('Your data and seed matrices need to have the same number of rows (voxels)')
elseif size(IC.data,1)~=size(IC.sl_indices,1)
    error('Your data and searchlight indices need to have the same number of rows (voxels)')
elseif size(IC.seed,1)~=size(IC.sl_indices,1)
    error('Your seed and searchlight indices need to have the same number of rows (voxels)')
end

% Some basic prep:
for conds=1:size(IC.conditions,1)
conds_labelled(conds,:)=IC.conditions(conds,:).*conds;
end
conds_labelled=sum(conds_labelled);

% Calculate multi-voxel pattern discriminability (MVP-D) for each
% searchlight:
clear mvpd_searchlights
for a=1:length(IC.sl_indices)

% Displaying progress:
current_percentage = 10*(ceil(a/length(IC.sl_indices)*100/10));
next_percentage = 10*(ceil((a+1)/length(IC.sl_indices)*100/10));
if next_percentage>current_percentage
    disp(['Calculating MVP discriminability for searchlights = ' num2str(current_percentage) '% complete'])
end

clear mean_training_corrs
clear timeseries_corrs
clear mvpd
for folds=1:max(IC.folds)
nDATA=sum(IC.selector(:,IC.folds==folds));
trainlabels=conds_labelled(:,(IC.selector.*(IC.folds~=folds))==1);
testlabels=conds_labelled(:,(IC.selector.*(IC.folds==folds))==1);
train=IC.data(IC.sl_indices(a,IC.sl_indices(a,:)>0),(IC.selector.*(IC.folds~=folds))==1);
test=IC.data(IC.sl_indices(a,IC.sl_indices(a,:)>0),(IC.selector.*(IC.folds==folds))==1);

for conds=1:size(IC.conditions,1)
mean_training_corrs(:,conds)=mean(train(:,trainlabels==conds),2);
timeseries_corrs(:,conds)=corr(test,mean_training_corrs(:,conds));
end

for trials=1:size(timeseries_corrs,1)
trial_corrs=timeseries_corrs(trials,:);
% Selecting the strongest correlation of the 'incorrect' classes:
trial_corrs(testlabels(trials))=-9;
[top top_index]=max(trial_corrs);
mvpd((folds-1).*nDATA+trials,1)=fisher_z(timeseries_corrs(trials,testlabels(trials)))-fisher_z(timeseries_corrs(trials,top_index));
end
end
mvpd_searchlights(:,a)=mvpd;
end

% Calculate MVP-D for the seed:
clear mean_training_corrs
clear timeseries_corrs
clear mvpd
clear mvpd_seed
for folds=1:max(IC.folds)
nDATA=sum(IC.selector(:,IC.folds==folds));
trainlabels=conds_labelled(:,(IC.selector.*(IC.folds~=folds))==1);
testlabels=conds_labelled(:,(IC.selector.*(IC.folds==folds))==1);
train=IC.data(IC.seed==1,(IC.selector.*(IC.folds~=folds))==1);
test=IC.data(IC.seed==1,(IC.selector.*(IC.folds==folds))==1);

for conds=1:size(IC.conditions,1)
mean_training_corrs(:,conds)=mean(train(:,trainlabels==conds),2);
timeseries_corrs(:,conds)=corr(test,mean_training_corrs(:,conds));
end

for trials=1:size(timeseries_corrs,1)
trial_corrs=timeseries_corrs(trials,:);
% Selecting the strongest correlation of the 'incorrect' classes:
trial_corrs(testlabels(trials))=-9;
[top top_index]=max(trial_corrs);
mvpd((folds-1).*nDATA+trials,1)=fisher_z(timeseries_corrs(trials,testlabels(trials)))-fisher_z(timeseries_corrs(trials,top_index));
end
end
mvpd_seed(:,1)=mvpd;

% Calculate IC!
disp(['Calculating IC . . .'])
IC.results=corr(mvpd_searchlights(:,:),mvpd_seed(:,1),'type','Spearman');
disp(['Done!'])

function z=fisher_z(r)
z=0.5*log((1+r)./(1-r));
end
end
