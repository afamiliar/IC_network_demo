function [IC pvalues] = permute_ROI_IC(IC,block_length)

% Performs a within-subject permutation test to determine statistical
% significance of an ROI-based IC analysis [IC Toolbox Version 1.1]
%
% [IC PVALUES] = PERMUTE_ROI_IC(IC,BLOCK_LENGTH)
%
% After running an ROI-based IC analysis (through run_ROI_IC), you might
% want to know if your results are statistically significance. There are
% plenty of options for testing this at the group level, but permutation
% testing is a good choice for individuals. This script performs a
% permutation test (with 10,000 iterations) for whether ROI connectivity
% strengths are significantly above zero (a one-tailed test).
%
%
% IC: RUN_ROI_IC will ensure the IC structure is ready for this function.
%
% BLOCK_LENGTH(optional): The number of consecutive TRs within each block.
% It is VERY important that your permutation test respects temporal
% dependencies between data-points. This is most apparent in a block-design
% - we should not break-up blocks when we permute our labels. Include the
% number of TRs in each block here to ensure their trials are kept together
% in the shuffling process. Obviously, if you have had to throw out the
% occasional TR (for movement, etc), then this will not function as
% intended so make sure that EVERY block is BLOCK_LENGTH long if you are
% using script.
%
% The results are placed in the IC structure (IC.perm_pvalues).
%
% If you include a second output (PVALUE), the script will also place the
% p-values in your workspace. 
% 
%
% Example: [IC P-VALUE] = permute_ROI_IC(IC,10)
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
if isfield(IC,'ROI_MVPD')==0
    error('Your IC structure doesn''t include the results from an IC ROI analysis. You should have already used run_ROI_IC.')
elseif isfield(IC,'ROI_names')==0
    error('Your IC structure doesn''t include your ROI names (IC.ROI_names). You will need these for the p-values.')
end

if nargin<2
warning('You did not include a block_length argument when you called this function. The test will still run but assumes no temporal dependency between time-points (e.g., a sparse event-related design). If you actually have a block design, the resulting p-values will be artificially inflated without this argument!!');
end

% Some basic prep:
nROIs=count(unique(IC.ROIs(IC.ROIs>0)));
if isint(length(IC.ROI_MVPD)/block_length)==0
	error('You have given a block_length that does not fit perfectly into your set of time-points (i.e., dividing the set of time-points by the block_length leaves a fraction')
end
clear IC.perm_results

for perms=1:10000
% Displaying progress:
current_percentage = 10*(ceil(perms/10000*100/10));
next_percentage = 10*(ceil((perms+1)/10000*100/10));
if next_percentage>current_percentage
    disp(['Permutations = ' num2str(current_percentage) '% complete'])
end

clear perm_ROI_MVPD
clear new_ROI_MVPD
for a=1:nROIs

if nargin==1
perm_ROI_MVPD(:,a)=IC.ROI_MVPD(randperm(length(IC.ROI_MVPD)),a);
else
nblocks=length(IC.ROI_MVPD)/block_length;
new_order=randperm(nblocks);
new_ROI_MVPD=IC.ROI_MVPD(:,a);
for shuffling=1:nblocks
new_ROI_MVPD((((shuffling-1)*block_length)+1):(shuffling*block_length))=IC.ROI_MVPD((((new_order(shuffling)-1)*block_length)+1):(new_order(shuffling)*block_length),a);
end
perm_ROI_MVPD(:,a)=new_ROI_MVPD;
end
end

% Calculate IC for this permutation:
perm_results=corr(perm_ROI_MVPD(:,:),'type','Spearman');
IC.perm_results(perms,:)=reshape(perm_results,1,size(perm_results,1)*size(perm_results,2));
end

% Calculate p-values:
real_results=reshape(IC.ROI_results,1,size(IC.ROI_results,1)*size(IC.ROI_results,2));
IC.perm_results=sort(IC.perm_results,1,'descend');

clear pvalues
for connection=1:length(real_results)
   for perms=1:10000
    if real_results(1,connection)>IC.perm_results(perms,connection)
       pvalues(1,connection)=perms/10000;
       break
    end
   end
   if real_results(1,connection)<=IC.perm_results(10000,connection)
   pvalues(1,connection)=1;
   end
end

pvalues=reshape(pvalues,size(IC.ROI_results,1),size(IC.ROI_results,2));
pvalues(eye(size(pvalues,1))==1)=NaN;
IC.perm_pvalues=pvalues;
disp(['Permutation p-values:'])
disp(pvalues)
disp(['A reminder of your ROI labels:'])
disp(IC.ROI_names)

end
