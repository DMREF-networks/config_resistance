% these lines are needed to run on the cluster, rather than locally
% c = parcluster; % cluster identified by default profile
% parpool(c);

% load in point clouds generated with Lloyd_create_pointclouds.m
% load paths_test_run.mat

% --- these are defined in the saved data file
% n_sims = 1;
% n_node = 200;
% crs = [0 0; 0 2000; 2000 2000; 2000 0]; % box to contain points

% Iter = [0 round(logspace(0,2,8))];
% Iter = [0 1];
% Iter = [0:100];
numIterations = max(Iter);
% ---

R_total = zeros(n_sims,length(Iter));
V_vars = cell(n_sims,length(Iter));
I_vars = cell(n_sims,length(Iter));
%path_length = zeros(n_sims,length(Iter));

for kn = 1:length(Iter)
    fprintf('iteration %d of %d\n',kn,length(Iter))
    
    R_tmp = zeros(n_sims,1);
    V_vars_tmp = cell(n_sims,1);
    I_vars_tmp = cell(n_sims,1);
    %path_length_tmp = zeros(n_sims,1);

    parfor k=1:n_sims
    % for k=1:n_sims
        x = x_loc{k,Iter(kn)+1};
        y = y_loc{k,Iter(kn)+1};

        tic
        [x,y,Amatrix] = find_corners_adjacency_B(x,y,crs);
        % [x,y,Amatrix] = find_corners_adjacency_voronoi(x,y,crs);
        % Amatrix = Adj_all{k,Iter(kn)+1};
        toc

        tic
        [R_tmp(k,1), V_vars_tmp{k,1}, I_vars_tmp{k,1}] = compute_voltage_Adj(x,y,Amatrix);
        toc
        % R_tmp(k,1) = compute_voltage_Adj(x,y,Amatrix);

        % network statistics
        % G = graph(Amatrix);
        % d = distances(G,1,n_node);  %shortest path length
        % path_length_tmp(k,1) = d;
      
    end

    R_total(:,kn) = R_tmp;
    V_vars(:,kn) = V_vars_tmp; % voltage for each node
    I_vars(:,kn) = I_vars_tmp; % current for each edge
    % path_length(:,kn) = path_length_tmp;

    % save('resistence_voltage.mat')
    % save('N200_DT_reweighted_resistance.mat')
    % save('N100_resistance_voltage.mat')
    save('N200_resistance_Ti64_test_configB_240917.mat')
end % loop over simulations


