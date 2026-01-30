function record = edge_flip_test(x, y, A)
% input: 
%   Adjacency matrix A
%   x: x-coordiantes of nodes
%   y: y-coordinate of nodes
% output: record is a structure containing
    % flipped_edge: [i, j] nodes of removed edge
    % new_edge: [k, l] nodes of added edge
    % delta_Reff: absolute value of change in the effective resistance
    % rel_delta_Reff: delta_Reff divided by the initial Reff
    % delta_Rtot: absolute value of change in the total resistance
    % rel_delta_Rtot: delta_Rtot divided by the initial Rtot
    % x
    % y
    % quad_nodes: [i, j, k, l] nodes of the quadrilateral 

n = size(A,1);  % number of nodes

% compute effective resistance between these two nodes
index1 = n;
index2 = 1;



% --- Get a  list of all triangles
tri_list = [];  % list of triangles in network
for i = 1:n
    for j = i+1:n
        if A(i,j)
            for k = j+1:n
                if A(i,k) && A(j,k)
                    tri_list = [tri_list; i, j, k];
                end
            end
        end
    end
end

% Enumerate Edges of the Triangles
edge2tri = containers.Map('KeyType','char','ValueType','any');

for t = 1:size(tri_list, 1)
    tri = tri_list(t, :);
    edges = [tri([1 2]); tri([2 3]); tri([3 1])];
    for e = 1:3
        i = min(edges(e,1), edges(e,2));
        j = max(edges(e,1), edges(e,2));
        key = sprintf('%d_%d', i, j);
        if isKey(edge2tri, key)
            edge2tri(key) = [edge2tri(key); t];
        else
            edge2tri(key) = t;
        end
    end
end

flipped = false;
keys_list = keys(edge2tri);

% --- Flip Edges Once at Each Time
record = struct('flipped_edge', {}, 'new_edge', {}, 'delta_Reff', {},'rel_delta_Reff', {},'delta_Rtot', {},'rel_delta_Rtot', {}, 'x', {}, 'y', {},'quad_nodes', {});

% starting values of the resistances
[Reff0,Rtot0,~,~,~] = compute_Rs(x, y, A, index1,index2);

for idx = 1:length(keys_list)
    key = keys_list{idx};
    tri_ids = edge2tri(key);
    if length(tri_ids) ~= 2
        continue  % not a shared edge
    end

    tri1 = tri_list(tri_ids(1), :);
    tri2 = tri_list(tri_ids(2), :);
    shared = intersect(tri1, tri2);  % should be 2 nodes
    if numel(shared) ~= 2
        continue
    end
    opp1 = setdiff(tri1, shared);
    opp2 = setdiff(tri2, shared);

    if isempty(opp1) || isempty(opp2)
        continue
    end
    i = shared(1); j = shared(2);
    quad = [opp1, i, j, opp2];

    if isConvex_from_coords(quad,x,y)
        % conducting flip if it's a convex quadrilateral
        Atmp = A;
        Atmp(i,j) = 0; Atmp(j,i) = 0;
        Atmp(opp1, opp2) = 1; Atmp(opp2, opp1) = 1;
        flipped = true;

        [Reff_new,Rtot_new,~,~,~] = compute_Rs(x, y, Atmp, index1,index2);

        delta = abs(Reff_new - Reff0);
        delta_Rtot = abs(Rtot_new - Rtot0);
        center = mean([x(i), x(j)]);
        center_y = mean([y(i), y(j)]); 
        % Recording the location of edge that results in Max Delta Reff
        rel_delta = 100*delta./Reff0;
        rel_delta_Rtot = 100*delta_Rtot./Rtot0;
        % Build one struct with all fields at once
        record(end+1) = struct( ...
        'flipped_edge', [i, j], ...
        'new_edge', [opp1, opp2], ...
        'delta_Reff', delta, ...
        'rel_delta_Reff',rel_delta,...
        'delta_Rtot', delta_Rtot, ...
        'rel_delta_Rtot',rel_delta_Rtot,...
        'x', center, ...
        'y', center_y, ...
        'quad_nodes', [opp1, i, j, opp2]);
    end
end

if ~flipped
    disp("No legal flip found.");
end
