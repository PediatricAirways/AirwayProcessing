function [mean0, left1, right1, left2, right2] = pointwisePCA(data)
% data: n * p * q
% n: number of cases, p: number of points in each case, q: dimension of each point

mean0 = zeros( size(data, 2), size(data, 3) );
left1 = zeros( size(data, 2), size(data, 3) );
left2 = zeros( size(data, 2), size(data, 3) );
right1 = zeros( size(data, 2), size(data, 3) );
right2 = zeros( size(data, 2), size(data, 3) );

dataAlign = zeros( size(data, 1), size(data, 3) );
meanData = zeros( size(data, 3), 1 );
for iJ = 1:size(data, 2)
    for iK = 1:size(data, 3)
        meanData(iK) = mean( data(:, iJ, iK) );
        dataAlign(:, iK) = data(:, iJ, iK) - meanData(iK);
    end
    [coeff, score, latent] = princomp( dataAlign );
    mean0(iJ, :) = meanData;
    left1(iJ, :) = meanData - sqrt( latent(1) ) * coeff(:, 1);
    right1(iJ, :) = meanData + sqrt( latent(1) ) * coeff(:, 1);
    left2(iJ, :) = meanData - 2*sqrt( latent(1) ) * coeff(:, 1);
    right2(iJ, :) = meanData + 2*sqrt( latent(1) ) * coeff(:, 1); 
    plot( data(:, iJ, 1), data(:, iJ, 2), 'b*' );
    hold on
    plot( meanData(1), meanData(2), 'gs' );
    plot( [left2(iJ, 1) right2(iJ, 1)], [left2(iJ, 2) right2(iJ, 2)], 'r-' );
    hold off
end

