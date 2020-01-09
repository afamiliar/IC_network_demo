function [subject] = export_to_MVPA_toolbox(IC,subject,new_pattern_name,template_pattern,template_file,file_space);

% Brings IC results into the Princeton MVPA toolbox subject structure [IC
% Toolbox Version 1.1]
%
% [SUBJECT] = EXPORT_TO_MVPA_TOOLBOX(IC,SUBJECT,NEW_PATTERN_NAME,TEMPLATE_PATTERN,TEMPLATE_FILE)
%
% Takes in an IC structure containing a results field (IC.results) and
% writes it as a new pattern in the given subject structure. This is only
% relevant if you have run an EXPLORATORY IC analysis (there is nothing to
% export from an ROI analysis...)
%
% IC: The IC structure containing a results field (IC.results).
%
% SUBJECT: The subject structure (in the Princeton MVPA toolbox) that will
% contain the new pattern.
%
% NEW_PATTERN_NAME: The name of the new pattern. Surround the name with
% apostrophes.
%
% TEMPLATE_PATTERN: A pattern to act as a template for the new results
% pattern. This should contain the same number of voxels (in the same
% order) as the analyzed data - the pattern actually containing the data 
% would do it! Surround the name with apostrophes.
%
% TEMPLATE_FILE (optional): If you want to also write out the IC results as
% AFNI BRIK and HEAD files, include the name of an existing file to use as
% a template. This should be at the same resolution as the data you
% analyzed - the actual data file would work well! Surround the name with 
% apostrophes.
%
% FILE_SPACE (optional): If you used the above argument to write out the IC
% results, and want them in Talairach space, write '+tlrc' as the last 
% argument. Without this, the file will be given a +orig suffix.
% 
% 
% Example: [subject] = export_to_MVPA_toolbox(IC,subj,'my_IC_results','data_z','data+tlrc','+tlrc')
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

if nargin<4
    error('You are missing inputs - have you listed at least four?')
elseif isfield(IC,'results')==0
    error('Your IC structure doesn''t have a results structure. You should have already run IC.')
elseif length(get_mat(subject,'pattern',template_pattern))~=length(IC.results);
    error('Your template pattern and IC results don''t seem to have the same number of voxels. These should match.')
end

subject=duplicate_object(subject,'pattern',template_pattern,new_pattern_name);
subject=set_mat(subject,'pattern',new_pattern_name,IC.results);

if nargin==5
    write_to_afni(subject,'pattern',new_pattern_name,template_file);
elseif nargin==6
    write_to_afni(subject,'pattern',new_pattern_name,template_file,'view',file_space);
end
