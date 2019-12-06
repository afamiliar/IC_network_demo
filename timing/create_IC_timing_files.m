%% create timing files for IC analysis
%   amf
%   Nov 2019
%
%   This prepares timing data for input to the IC toolbox.
%
%   -- 'conditions' is a matrix (# conditions X # time-points) with binary
%   values for each condition and time-point (1 if time-point is included
%   in that condition, 0 if it is not)
%
%   -- 'selector' is a vector with one element per time-point,
%   1 if the time-point should be included in the analysis and
%   0 if it is a time-point of no-interest.
%
%   -- 'folds' is a vector with one element per time-point, value
%   corresponds to run number
%
%%

nRuns = 5;
conditions = [];
selector = [];
folds = [];
for r = 1:nRuns
    % load original timing files, each containing a binary vector with one
    % element per time-point (1 if time-point included in that condition, 0
    % if not)
    timing_face1 = load(['timing_files/r' int2str(r) '_face1.txt'])';
    timing_face2 = load(['timing_files/r' int2str(r) '_face2.txt'])';
    timing_face3 = load(['timing_files/r' int2str(r) '_face3.txt'])';
    timing_face4 = load(['timing_files/r' int2str(r) '_face4.txt'])';
    timing_face5 = load(['timing_files/r' int2str(r) '_face5.txt'])';
    timing_face6 = load(['timing_files/r' int2str(r) '_face6.txt'])';
    timing_face7 = load(['timing_files/r' int2str(r) '_face7.txt'])';
    timing_face8 = load(['timing_files/r' int2str(r) '_face8.txt'])';
    timing_face9 = load(['timing_files/r' int2str(r) '_face9.txt'])';
    conditions = [conditions, ...
                    [timing_face1; timing_face2; timing_face3;...
                     timing_face4; timing_face5; timing_face6;...
                     timing_face7; timing_face8; timing_face9]];

    % load time-points of no-interest (blank trials & trials in which
    % the image was a repeat of the one preceeding it)
    blanks  = load(['timing_files/r' int2str(r) '_blank.txt'])';
    repeats = load(['timing_files/r' int2str(r) '_repeat.txt'])';
    exclude = blanks+repeats;
    
    selector = [selector, [ones(1,size(exclude,2)) - exclude]];

    folds = [folds, repmat(r,1,size(exclude,2))];
end

save('conditions.mat','conditions')
save('selector.mat','selector')
save('folds.mat','folds')
