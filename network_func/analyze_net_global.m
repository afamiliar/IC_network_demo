function global_props = analyze_net_global(A,local_props)

    comm_struct = local_props.data(:,1);
    cluster_coefs = local_props.data(:,2);

    modularity = max(comm_struct);
    cluster_coef = mean(cluster_coefs(:));
    char_path_length = charpath(A);
    [SWP, ~, ~] = small_world_propensity(A, 'B');
    assort_val = assortativity(A,0);
    
    global_props.data = [modularity;cluster_coef;char_path_length;SWP;assort_val];
    global_props.labels = {'modularity';'cluster coefficient';...
                            'characteristic path length'; 'small worldedness';...
                            'assortativity'};
end