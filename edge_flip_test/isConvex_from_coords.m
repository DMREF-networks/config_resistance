function flag = isConvex_from_coords(indices, x, y)
% ISCONVEX_FROM_COORDS Check whether four points form a convex quadrilateral.
%
% INPUTS:
%   indices : a 1×4 vector containing the indices of four nodes (order can be arbitrary)
%   x       : x-coordinates of all nodes (vector of length n)
%   y       : y-coordinates of all nodes (vector of length n)
%
% OUTPUT:
%   flag    : logical (true or false)
%             - true  → if the four points form a convex quadrilateral
%             - false → otherwise (e.g., if concave)
%
% USAGE EXAMPLE:
%   quad = [opp1, i, j, opp2];  % any 4 node indices
%   if isConvex_from_coords(quad, x, y)
%       % safe to perform edge flip
%   end

    coords = [x(indices(:)), y(indices(:))];  % Always gives 4×2
    center = mean(coords, 1);
    angles = atan2(coords(:,2) - center(2), coords(:,1) - center(1));
    [~, order] = sort(angles);
    coords_sorted = coords(order, :);

    z = zeros(4,1);
    for i = 1:4
        A = coords_sorted(i,:);
        B = coords_sorted(mod(i,4)+1,:);
        C = coords_sorted(mod(i+1,4)+1,:);
        AB = B - A;
        BC = C - B;
        z(i) = AB(1)*BC(2) - AB(2)*BC(1);
    end

    flag = all(z > 0) || all(z < 0);
end
