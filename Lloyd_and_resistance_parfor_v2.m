% these lines are needed to run on the cluster, rather than locally
% c = parcluster; % cluster identified by default profile
% parpool(c);

% load in point clouds generated with Lloyd_create_pointclouds_parfor.m
load filename.mat

% --- these are defined in the saved data file
% n_sims = 1;
% n_node = 200;
% crs = [0 0; 0 2000; 2000 2000; 2000 0]; % box to contain points

% Here, which iterations to compute resistance at can be defined, if not already in filename.mat
% Iter = [0 round(logspace(0,2,8))];
% Iter = [0 1];
% Iter = [0:100];
numIterations = max(Iter);
% ---

R_total = zeros(n_sims,length(Iter));
V_vars = cell(n_sims,length(Iter));
I_vars = cell(n_sims,length(Iter));

for kn = 1:length(Iter)
    fprintf('iteration %d of %d\n',kn,length(Iter))
    
    R_tmp = zeros(n_sims,1);
    V_vars_tmp = cell(n_sims,1);
    I_vars_tmp = cell(n_sims,1);

    parfor k=1:n_sims
        x = x_loc{k,Iter(kn)+1};
        y = y_loc{k,Iter(kn)+1};

        tic
        [x,y,Amatrix] = find_corners_adjacency_B(x,y,crs);

        [R_tmp(k,1), V_vars_tmp{k,1}, I_vars_tmp{k,1}] = compute_voltage_Adj(x,y,Amatrix);
        toc
      
    end

    R_total(:,kn) = R_tmp;
    V_vars(:,kn) = V_vars_tmp; % voltage for each node
    I_vars(:,kn) = I_vars_tmp; % current for each edge

    save('resistence_voltage_current.mat')
end % loop over simulations


