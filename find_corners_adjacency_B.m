function [x,y,Amatrix] = find_corners_adjacency_B(x,y,crs)
% ensures node 1 is NW corner and node n is SE corner
% this is experimental "configuration B"
% Voronoi in bounded box to create adjacency matrix

n = numel(x);   % number of points

% --- find NW and SE corners:  this should be 'configuration B'
d2 = (x-0).^2 + (y-2000).^2;
[~,ind_v1] = min(d2);

x = [x(ind_v1); x([1:ind_v1-1 ind_v1+1:end])];
y = [y(ind_v1); y([1:ind_v1-1 ind_v1+1:end])];

d2 = (x-2000).^2 + (y-0).^2;
[~,ind_v2] = min(d2);

x = [x([1:ind_v2-1 ind_v2+1:end]); x(ind_v2)];
y = [y([1:ind_v2-1 ind_v2+1:end]); y(ind_v2)];

% ADD-IN to calculate adjacency matrix via voronoi nearest neighbors
% from https://people.sc.fsu.edu/~jburkardt/presentations/voronoi_neighbors.pdf
% re-written by Newhall to fix issue that sometimes leaves a node disconnected

rgx = max(crs(:,1))-min(crs(:,1));
rgy = max(crs(:,2))-min(crs(:,2));
rg = max(rgx,rgy);
midx = (max(crs(:,1))+min(crs(:,1)))/2;
midy = (max(crs(:,2))+min(crs(:,2)))/2;

% add 4 additional edges
xA = [x; midx + [0;0;-5*rg;+5*rg]];
yA = [y; midy + [-5*rg;+5*rg;0;0]];

[vi,ci]=voronoin([xA,yA]); % vi: verticies of the voronoi cells
                           % ci{j}: list of verticies vi for the j^th cell

% use polyshap and intersect to crop the cells
poly_box = polyshape(crs);
poly_x = cell(n,1);
poly_y = cell(n,1);
% figure(1);clf;hold on
for k = 1:n
    tmp = intersect( poly_box,polyshape(vi(ci{k},1),vi(ci{k},2)) );
    poly_x{k} = round(tmp.Vertices(:,1),6); % round due to numerical error when equating
    poly_y{k} = round(tmp.Vertices(:,2),6);
    % plot(tmp); text(x(k),y(k),sprintf('%d',k))
end

% ---
Amatrix = zeros( n, n );
for i = 1 : n
    for j = i + 1 : n
        % s = size ( intersect ( c{i}, c{j} ) ); % intersection of two sets (polygons)
        [xi,yi] = intersect( [poly_x{i},poly_y{i}], [poly_x{j},poly_y{j}] ,"rows");
        if ( 1 < length(xi) ) % two numbers ->  cells share one edge, one number -> share vertex
            Amatrix(i,j) = 1;
            Amatrix(j,i) = 1;
        end
    end
end

end % end function