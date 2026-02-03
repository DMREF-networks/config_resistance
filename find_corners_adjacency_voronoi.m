function [x,y,Amatrix] = find_corners_adjacency_voronoi(x,y,crs)
% ensures node 1 is NW corner and node n is SE corner
% Voronoi in bounded box to create adjacency matrix
% input [x,y] are the points in each Voronoi cell
% output [x,y] are the nodes of the graph (verticies of Voronoi tesselation

n = numel(x);   % number of points

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
nodes = unique( [cell2mat(poly_x) cell2mat(poly_y)] , 'rows' );
x = nodes(:,1);
y = nodes(:,2);
% find NW and SE corners:
d2 = (x-0).^2 + (y-2000).^2;
[~,ind_v1] = min(d2);

x = [x(ind_v1); x([1:ind_v1-1 ind_v1+1:end])];
y = [y(ind_v1); y([1:ind_v1-1 ind_v1+1:end])];

d2 = (x-2000).^2 + (y-0).^2;
[~,ind_v2] = min(d2);

x = [x([1:ind_v2-1 ind_v2+1:end]); x(ind_v2)];
y = [y([1:ind_v2-1 ind_v2+1:end]); y(ind_v2)];

% add edges of polygon to graph object
G = graph();

for k = 1:n
    pxy = [poly_x{k} poly_y{k}];
    for j = 1:size(pxy,1)-1 % number of rows
        indA = find( x==pxy(j,1) &  y==pxy(j,2));
        indB = find( x==pxy(j+1,1) &  y==pxy(j+1,2));
        G = addedge(G,indA,indB);
    end
    % connect back to start
    indA = find( x==pxy(j+1,1) &  y==pxy(j+1,2));
    indB = find( x==pxy(1,1) &  y==pxy(1,2));
    G = addedge(G,indA,indB);
end

G = simplify(G);    % removes redundant edges
Amatrix = adjacency(G);

end % end function