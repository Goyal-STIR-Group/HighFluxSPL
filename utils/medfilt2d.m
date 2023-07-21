function y = medfilt2d(x,order)
% Inputs
% x: 2D image
% order: order of median filter in each dimension. Should be positive, odd. Default = 3;
% Written by Joshua Rapp to avoid requiring Image Processing Toolbox

if nargin < 2
    order = 3;
end

ord_dist = floor((order-1)/2);

if size(order(:),1) == 1
    ord_i = ord_dist;
    ord_j = ord_dist;
elseif size(order(:),1) == 2
    ord_i = ord_dist(1);
    ord_j = ord_dist(2);
end
    

[Nr,Nc] = size(x);
y = zeros(Nr,Nc);

for ii = 1:Nr
    % rows
    min_i = max(1, ii-ord_i);
    max_i = min(Nr, ii+ord_i);
    
    for jj = 1:Nc
        % cols
        min_j = max(1, jj-ord_j);
        max_j = min(Nc, jj +ord_j);
        
        window = x(min_i:max_i,min_j:max_j);
        y(ii,jj) = median(window(:));
    
    end
end
end

