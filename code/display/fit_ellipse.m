% x, y are vectors of coordinates
function a = fit_ellipse(x,y)
% Build design matrix
D = [x.*x x.*y y.*y x y ones(size(x))];
% Build scatter matrix
S = D'*D;
% Build 6X6 constraint matrix
C(6,6) = 0; C(1,3) = -2; C(2,2) = 1; C(3,1) = -2;
% Solve generalized eigensystem
[gevec, geval] = eig(S, C);
% Find the only negative eigenvalue
[NegR, NegC] = find(geval<0 & ~isinf(geval));
% Get fitted parameters
a = gevec(:, NegC);