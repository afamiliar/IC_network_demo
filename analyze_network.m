%% calculate network properties using graph theory
%   amf
%   Jan 2020
%
%       If you haven't already, run 'analyze_IC_brainnetome.m' to generate
%       IC-based networks.
%
%       - Loads network
%       - Normalizes correlation values
%       - Calculates a set of local and global network properties using the
%           Brain Connectivity Toolbox
%       - Saves results
%
%%
addpath(genpath(pwd))

files = struct2cell(dir('IC_results/*.mat'));

%%
for f = 1:size(files,2)
    
    % load subject data
    load(files{1,f});

    % normalize correlation values
    A = atanh(IC_matrix); % Fisher normalize w/in subj
    A(A<0)=0;
    A(A==Inf)=0;
    
    % save to group-level structure (if running group analysis)
%     group{f} = A;
    
    % calculate local and global network properties
    local_props  = analyze_net_local(A);
    global_props = analyze_net_global(A, local_props);

    % save results
    save(['network_results/' files{1,f}],'local_props','global_props')

    clear IC_matrix
    clear A
    clear local_props
    clear global_props

end

%% if aggregating at group-level, could average across networks (see below)
%   and then analyze network properties

% % average across subj's
% Y = cat(3,group{:});
% group_ave = mean(Y,3);
