function [median_curve, centerInf, centerSup, minCurve, maxCurve] = wpbplot(data)
% find the weighted pointwise boxplot for data p x n

if nargin < 2, w_j = ones(size(data, 2), 1) / size(data, 2); end

factor  = 1.5;

median_curve = zeros( size(data, 1), 1 );
median_value = zeros( size(data, 1), 1 );
centerInf = zeros( size(data, 1), 1 );
centerSup = zeros( size(data, 1), 1 );
minCurve  = zeros( size(data, 1), 1 );
maxCurve = zeros( size(data, 1), 1 );
for iI = 1:size( data, 1 )
    for iJ = 1:size( data, 2 )
        sum_value = sum( w_j .* abs( data( iI, iJ ) - reshape( data( iI, : ), size(w_j) ) ) );
        if iJ == 1 || median_value( iI ) > sum_value
            median_curve( iI ) = data( iI, iJ );
            median_value( iI ) = sum_value;
        end
    end
	[data_sort, sort_id] = sort( reshape( data( iI, : ), size( data, 2 ), 1 ) );
    wj_sort = w_j( sort_id );
	wj_value_sum = 0;
    flag_quartile_25 = 0;
    flag_quartile_75 = 0;
    flag_percentile_993 = 0;
    for iJ = 1:length( wj_sort )
            wj_value_sum = wj_value_sum + wj_sort( iJ );
            if flag_quartile_25 == 0 && wj_value_sum >= 0.25
                centerInf( iI ) = data_sort( iJ );
                flag_quartile_25 = 1;
            end
            if flag_quartile_75 == 0 && wj_value_sum >= 0.75
                centerSup( iI ) = data_sort( iJ );
                flag_quartile_75 = 1;
            end
            if flag_percentile_993 == 0 && wj_value_sum >= 0.993
                flag_percentile_993 = iJ;
            end
    end
    dist = factor * ( centerSup(iI) - centerInf(iI) );
    minValue = centerInf(iI) - dist;
    maxValue = centerSup(iI) + dist;
    outly = or( reshape( data(iI, :), size(data, 2), 1 ) < minValue, reshape( data(iI, :), size(data, 2), 1 ) > maxValue );
    outly(sort_id(flag_percentile_993+1:end)) = 1;
    minCurve(iI) = min( data( iI, sort_id( find(1-outly) ) ) );
    maxCurve(iI) = max( data( iI, sort_id( find(1-outly) ) ) );
end