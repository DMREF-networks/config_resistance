function [R_total, V, I] = compute_voltage_Adj(x,y,A)
% inputs = x and y positions of nodes, A adjacency matrix
% outputs = R_total effective resistance, V voltage each node, I current
%         each edge

cur = 10; % applied current diagonally across the network

n_node = length(x); % number of nodes

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

R = A./tmp; % adjacency matrix with 1/R_{ij} for each edge
for k=1:n_node % set the diagonal to be zero
    R(k,k)=0;
end
% ---

D = diag( R*ones(n_node,1) ); % diagonal matrix with row-sums of R

% we have one degree of freedom, set V(end)=0 wlog
Amat = D(1:n_node-1,1:n_node-1)-R(1:n_node-1,1:n_node-1);
Ivec = [cur; zeros(n_node-2,1)];

V = [Amat\Ivec; 0];

R_total = V(1)/cur; % voltage difference across network is V(1) since V(n_node)=0

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


% Extra code for plotting voltage at each node using a colorbar
% figure(3);clf;hold on
% plot(G,'XData',x,'YData',y,'NodeLabel',[],...
%     'NodeCData',V,'markersize',10)
% colorbar
% 
% Extra code for plotting current along each edge using a colorbar
% figure(4);clf;hold on
% plot(G,'XData',x,'YData',y,'NodeLabel',[],...
%     'EdgeCData',I,'linewidth',4,'nodecolor','k','markersize',3)
% colorbar
% 
% title('Voltage')
