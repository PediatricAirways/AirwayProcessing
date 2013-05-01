function [median] = pointwiseMedian(data, w_j)
% data: n * p * q, n is the number of shapes, p*q is the dimension

if nargin < 2, w_j = ones(size(data, 1), 1) / size(data, 1); end

median = zeros( size(data, 2), size(data, 3) );
median_value = zeros( size(data, 2), 1 );

for iI = 1:size( data, 2 )
    for iJ = 1:size( data, 1 )
        sum_value = 0;
        for iK = 1:size( data, 1 )
            sum_value = sum_value + w_j(iK) * norm( reshape( data(iJ, iI, :) - data(iK, iI, :), size(data, 3), 1) );
        end
        if iJ == 1 || median_value(iI) > sum_value
            median(iI, :) = data(iJ, iI, :);
            median_value(iI) = sum_value;
        end
    end
end

