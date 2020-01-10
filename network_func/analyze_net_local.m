function local_props = analyze_net_local(A)

    comm_struct = community_louvain(A);             % community modules
    cluster_coefs = clustering_coef_bu(A);          % weighted clustering coefficient (like degree centrality)
    % cluster_coefs = clustering_coef_wu(A);          % weighted clustering coefficient (like degree centrality)
    P = participation_coef(A,comm_struct);          % participation coefficient
    BC = betweenness_wei(A);                        % betweeness centrality
    core_periph_ids = [core_periphery_dir(A)+1]';   % core-periphery partition
    
    local_props.data = [comm_struct,cluster_coefs,P,BC,core_periph_ids];
    local_props.labels = {'community modules','cluster coefficient',...
                            'participation coefficient','betweeness centrality',...
                            'core-periphery ids'};
end