function [pc] = compPrincipalDirection(x)

avg = mean(x, 1);
for iJ = 1:size(x, 2)
    x_center(:, iJ) = x(:, iJ) - avg(iJ);
end
[U, S, pc] = svd( x_center ./ sqrt(size(x, 1)-1), 0 );
