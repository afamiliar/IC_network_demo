function [IC IC_matrix] = run_ROI_IC(IC)

% Measures Informational Connectivity (IC) between different regions (ROIs)
% [IC Toolbox Version 1.1]
%
% [IC IC_MATRIX] = RUN_ROI_IC(IC)
%
% Takes in an IC structure (details below), calculates IC between the
% chosen regions, and outputs results to the IC structure. The function
% produces IC.ROI_results: a matrix with one row and column for each ROI.
% Each value is the IC between regions: for example, IC.ROI_results(1,2) is
% the IC between the first and second ROIs. If you include a second output
% when you call the function, the results matrix is also written into your
% workspace, and a figure is opened that shows the relative IC magnitudes
% between regions (Note: This is not a confusion matrix, even if it may
% resemble one!)
%
% The script outputs its progress in 10% increments.
%
% IC: If you already use the Princeton MVPA toolbox (or don''t already have
% your data in Matlab), see 'import_from_mvpa_toolbox' to create the IC
% structure. Alternatively, you can create the structure yourself, with the
% following necessary fields (you can ignore this if you have already used
% import_from_mvpa_toolbox): 
%   - IC.data: a 2-dimensional matrix with V rows (one per voxel) and N
%       columns (one per trial).
%   - IC.ROIs: a column vector with V rows (one per voxel) that should
%       match IC.data (i.e., voxel 42 in this vector is also voxel 42 in
%       IC.data). The vector should give each region's voxels the same
%       positive integer. For example, region 1's voxels will be marked
%       with 1, region 2's voxels wil have 2, etc. Others = 0.
%   - IC.conditions: a matrix with one row per condition and N columns for
%       time-points. A time-point has a zero if that row's condition is not
%       active and 1 when it is. Only one condition should be active at any
%       time (i.e., no column sums to more than one). The conditions should
%       already be shifted for the hemodynamic delay. 
%   - IC.folds: 1-by-N vector, with each time-point numbered by fold. In many
%       cases, the folds will be runs. This vector determines which
%       time-points act as testing (fold M) and training (all-but-M folds).
%   - IC.selector: 1-by-N vector that indicates which time-points should be
%       included (1) or excluded (0). At a minimum, rest should be excluded
%       (unless rest is one of your classified conditions). Remember to
%       base this selector on time-points AFTER shifting for the hemodynamic
%       delay.
%   - IC.ROI_names: a cell structure with strings of the ROI names (in the
%       same order as the numbers in the IC.ROIs vector.
%
%
% Example: [IC IC_matrix] = run_ROI_IC(IC)
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
% if isfield(IC,'data')==0
%     error('Your IC structure doesn''t have a data subfield. Have you named it correctly? (.data)')
% elseif isfield(IC,'ROIs')==0
%     error('Your IC structure doesn''t have an ROIs subfield. Have you named it correctly? (.ROIs)')
% elseif isfield(IC,'conditions')==0
%     error('Your IC structure doesn''t have a conditions subfield. Have you named it correctly? (.conditions)')
% elseif isfield(IC,'folds')==0
%     error('Your IC structure doesn''t have a folds subfield. Have you named it correctly? (.folds)')
% elseif isfield(IC,'selector')==0
%     error('Your IC structure doesn''t have a selector subfield. Have you named it correctly? (.selector)')
% elseif size(IC.data,2)~=size(IC.conditions,2)
%     error('Your data and condition matrices need to have the same numbers of time-points')
% elseif size(IC.data,2)~=size(IC.selector,2)
%     error('Your data and selection matrices need to have the same numbers of time-points')
% elseif size(IC.selector,2)~=size(IC.conditions,2)
%     error('Your condition and selection matrices need to have the same numbers of time-points')
% elseif size(IC.folds,2)~=size(IC.conditions,2)
%     error('Your folds and condition matrices need to have the same numbers of time-points')
% elseif size(IC.folds,2)~=size(IC.selector,2)
%     error('Your folds and selection matrices need to have the same numbers of time-points')
% elseif size(IC.folds,2)~=size(IC.data,2)
%     error('Your folds and data matrices need to have the same numbers of time-points')
% elseif size(IC.data,1)~=size(IC.ROIs,1)
%     error('Your data and ROI vector need to have the same numbers of voxels')
% elseif size(IC.conditions,1)<2
%     error('You need to have at least 2 conditions')
% elseif length(unique(IC.folds))<2
%     error('You need to have at least 2 folds')
% elseif count(unique(IC.ROIs(IC.ROIs>0)))<2
%     error('You need to have at least 2 ROIs')
% end

% Some basic prep:
ROI_labels=unique(IC.ROIs(IC.ROIs>0));
nROIs=length(unique(IC.ROIs(IC.ROIs>0))); % 

for conds=1:size(IC.conditions,1)
    conds_labelled(conds,:)=IC.conditions(conds,:).*conds;
end
conds_labelled=sum(conds_labelled);

% Calculating the multi-voxel pattern discriminability (MVP-D) time-course
% for ROIs:
clear mvpd_ROIs
for a=1:nROIs
    % Displaying progress:
    current_percentage = 10*(ceil(a/nROIs*100/10));
    next_percentage = 10*(ceil((a+1)/nROIs*100/10));
    if next_percentage>current_percentage
        disp(['Calculating MVP discriminability: ' num2str(current_percentage) '% complete'])
    end

    clear mean_training_corrs
    clear timeseries_corrs
    clear mvpd
    for folds=1:max(IC.folds)
        nDATA=sum(IC.selector(:,IC.folds==folds));
        trainlabels=conds_labelled(:,(IC.selector.*(IC.folds~=folds))==1);
        testlabels=conds_labelled(:,(IC.selector.*(IC.folds==folds))==1);
        train=IC.data(IC.ROIs==ROI_labels(a),(IC.selector.*(IC.folds~=folds))==1);
        test=IC.data(IC.ROIs==ROI_labels(a),(IC.selector.*(IC.folds==folds))==1);

        timeseries_corrs=[]; %
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
    mvpd_ROIs(:,a)=mvpd;
end

IC.ROI_MVPD=mvpd_ROIs;

% Calculate IC!
disp(['Calculating IC . . .'])
IC.ROI_results=corr(mvpd_ROIs(:,:),'type','Spearman');
disp(['Done!'])
disp(['ROI results:'])
disp(IC.ROI_results)
disp(['A reminder of your ROI labels:'])
disp(IC.ROI_names)

if nargout==2
IC_matrix=IC.ROI_results;
figure
imagesc(IC.ROI_results)
colormap(hot)
set(gca,'CLim',[0 1])
colorbar
set(gca,'XTickLabel',IC.ROI_names)
set(gca,'YTickLabel',IC.ROI_names)
title('Informational Connectivity between regions');
end

function z=fisher_z(r)
z=0.5*log((1+r)./(1-r));
end
end
