function draw_flip_heatmap_edges_only(x, y, record, A, valueField)
% Draw a heatmap where only the flipped edges are colored
%
% Inputs:
%   x, y        - node coordinates (column vectors)
%   record      - struct array; each entry must contain:
%                   .flipped_edge = [node_i node_j]
%                   .(valueField) = scalar value to map color
%   A           - adjacency matrix (for plotting base network in gray)
%   valueField  - string specifying which field to use for color

if nargin < 5 || isempty(valueField)
    valueField = 'deltaReff';
end

ax = gca;
hold(ax, 'on');

%% 1) Draw base network (all edges in light gray)
[i_idx, j_idx] = find(triu(A, 1));
Xbase = [x(i_idx) x(j_idx) NaN(size(i_idx))]';
Ybase = [y(i_idx) y(j_idx) NaN(size(i_idx))]';
plot(Xbase(:), Ybase(:), 'Color', [0.6 0.6 0.6], 'LineWidth', 3);
plot(Xbase(:), Ybase(:), 'Color', [0.98 0.98 0.98], 'LineWidth', 2);

%% 2) Extract flipped edges and associated values
m = numel(record);
E = nan(m, 2);
vals = nan(m, 1);
mask = false(m, 1);

for k = 1:m
    if isfield(record(k), 'flipped_edge') && numel(record(k).flipped_edge) == 2 ...
            && isfield(record(k), valueField)
        E(k, :) = record(k).flipped_edge(:).';
        vals(k) = record(k).(valueField);
        mask(k) = true;
    end
end

E = E(mask, :);
vals = vals(mask);

if isempty(E)
    warning('No flipped edges found or missing field "%s".', valueField);
    return;
end

X = [x(E(:,1)) x(E(:,2)) NaN(size(E,1),1)]';
Y = [y(E(:,1)) y(E(:,2)) NaN(size(E,1),1)]';
C = [vals vals NaN(size(vals))]';

patch('XData', X(:), ...
      'YData', Y(:), ...
      'CData', C(:), ...
      'FaceColor', 'none', ...
      'EdgeColor', 'flat', ...           
      'LineWidth', 2, ...             
      'HandleVisibility','off');

%cmap_full = flipud(hot(256));  
%N = 250;                       
%cmap_cropped = cmap_full(end-N+1:end, :);  
%colormap(gca, cmap_cropped);

cmap = purple_cmap;
colormap(flipud(cmap));

axis equal tight;
box off;
hold(ax, 'off');
end