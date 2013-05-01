function [mean0, left1, right1, left2, right2] = vectorPCA(data)
% data: n * p

mean0 = (mean(data, 1))';
dataAlign = zeros( size(data) );
for iJ = 1:size( data, 2 )
    dataAlign(:, iJ) = data(:, iJ) - mean0(iJ);
end
[coeff, score, latent] = princomp( dataAlign, 'econ' );
left1 = mean0 - sqrt( latent(1) ) * coeff(:, 1);
right1 = mean0 + sqrt( latent(1) ) * coeff(:, 1);
left2 = mean0 - 2*sqrt( latent(1) ) * coeff(:, 1);
right2 = mean0 + 2*sqrt( latent(1) ) * coeff(:, 1);

cumsum(latent) ./ sum(latent)