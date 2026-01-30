function [R_total, R_totaleff_full, V, I, spectralgap] = compute_Rs(x, y, A, index1, index2)
% inputs = x and y positions of nodes, A adjacency matrix
% outputs = R_total: effective resistance, R_totaleff_full: Total Effective Resistance, V voltage each node, I current
%         each edge, spectralgap: the second smallest eigenvalues of the
%         Laplacian

cur = 10; % applied current diagonally across the network

n_node = length(x); % number of nodes
%here we read in Adj and xy, so the Amatrix is not arranged as the
%find_corner_adjacency, find the two nodes closest to the poles

% It is possible to have that the node we apply current is the last node on the list, 
% eliminate such situation by changing the direction of applying current
% index1 and 2 are nodes farthest in SW(0,0) and NE poles, and their distances r1 and r2
%[r1,index1] = min(sqrt(x.^2+y.^2)); 
%[r2,index2] = max(sqrt(x.^2+y.^2));
if index1 == n_node
    index1_1 = index2;
    index2 = index1;
    index1 = index1_1;
    cur = -cur; 
end

% --- add weights to edges
% first, resistance is proportional to the length of the edge
dx = x*ones(1,n_node)-ones(n_node,1)*x';
dy = y*ones(1,n_node)-ones(n_node,1)*y';
dr = sqrt(dx.^2 + dy.^2);    % ALL pairwise distances between nodes


% second, add cross-sectional area and resistivity of material with correct units
% - these are for the steel 17-4PH prints, resistivity 80 muOhm cm
% xarea = 3*1/10*1/10;    % 3mm * 1cm/10mm * 1mm * 1cm/10mm % this is for steel
% tmp = dr*75/2000*1/10*80*1/1000 / xarea;  
% 75mm/2000 nondim * 1cm/10mm * 80\muOhm.cm * 1mOhm/1000\muOhm

% - these are for the titanium Ti-64 prints, resistivity 178 muOhm cm
xarea = 2.75*1/10*1/10;    % 3mm * 1cm/10mm * 1mm * 1cm/10mm
tmp = dr*75/2000*1/10*178*1/1000 / xarea; 
% 75mm/2000 nondim * 1cm/10mm * 178\muOhm.cm * 1mOhm/1000\muOhm

R = A./tmp; % adjacency matrix with 1/R_{ij} (dr as a factor in temp) for each edge
% can change it to 
for k=1:n_node % set the diagonal to be zero
    R(k,k)=0;
end
% ---

D = diag( R*ones(n_node,1) ); % diagonal matrix with row-sums of R

% we have one degree of freedom, set V(end)=0 wlog

Amat = D([1:index2-1 index2+1:end],[1:index2-1 index2+1:end])-R([1:index2-1 index2+1:end],[1:index2-1 index2+1:end]);%graph Lapalacian

%Debug:looked at the dimension of matrices since it didn't match: l = dim(colspace(Amat));
Ivec = zeros(n_node-1,1);
if index1 < index2
    Ivec(index1,1) = cur;
else
    Ivec(index1-1,1) = cur;
end

L = D - R; % Graph Laplacian weighted
lambda_full_tmp = eig(L); % Compute eigenvalues (unsorted)
lambda_full = sort(lambda_full_tmp); % Sort in ascending order
spectralgap = lambda_full(2); % λ_2 the second smallest eigenvalues of the Laplacian
R_totaleff_full = n_node * sum(1./lambda_full(2:end)); % Ignore λ_1 = 0


%disp(['Total Resistance (Full Laplacian): ', num2str(R_totaleff_full)]);

Vtemp = Amat\Ivec;% solve
V = [Vtemp(1:index2-1); 0 ; Vtemp(index2:end)];

R_total = (V(index1)-V(index2))/cur; % voltage difference across network is V(1) since V(n_node)=0

% get the edge list
G = graph(A);
tmp = G.Edges;
tmp = table2array(tmp);
e_list = tmp(:,1:2);    % edges without the weights

% the positive direction is defined by the edge list order, 1st->2nd
dV = V(e_list(:,2))-V(e_list(:,1)); % voltage difference across each edge
Redge = 0*e_list(:,2);
for k=1:length(e_list)
    Redge(k) = R(e_list(k,1),e_list(k,2)); % this is 1/R actually
end
I = dV.*Redge;  % current along each edge
end % end function



% 
% Extra code for plotting current along each edge using a colorbar
% figure(4);clf;hold on
% plot(G,'XData',x,'YData',y,'NodeLabel',[],...
%     'EdgeCData',I,'linewidth',4,'nodecolor','k','markersize',3)
% colorbar
% 
% title('Voltage')

