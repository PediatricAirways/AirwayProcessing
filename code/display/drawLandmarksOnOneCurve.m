function drawLandmarksOnOneCurve( meanNormFile, areaFile, perimeterFile, ellipseAreaFile, landmarksIdOnCenterlineFile, curveAndLandmarksFile, outputPrefix )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Pediatric Airway Atlas Processing Code  ---
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%
%      http://www.apache.org/licenses/LICENSE-2.0 
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.
%   
%   Author: Yi Hong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read the mean points and compute dist in 1D
numf = 1;
numP = zeros(numf, 1);
for n = 1:numf
    fid = fopen( meanNormFile, 'r');
    tline = fgets(fid);
    count = 0;
    while ~feof(fid)
        tline = fgets(fid);
        count = count + 1;
        tmp = sscanf(tline, '%f');
        A(:, count, n) = tmp(1:3);
    end
    fclose(fid);
    numP(n) = count;
    dist(1,n) = 0;
    for i=2:count
        p1_x = A(1,i,n); p1_y = A(2,i,n); p1_z = A(3,i,n);
        p2_x = A(1,i-1,n); p2_y = A(2,i-1,n); p2_z = A(3,i-1,n);
        dist(i,n) = sqrt( (p1_x-p2_x)^2 + (p1_y-p2_y)^2 + (p1_z-p2_z)^2 ) + dist(i-1,n);
    end
    dist(1:count,n) = dist(1:count,n) ./ dist(count,n);
end

% read landmarks 
%for n = 1:numf
%    [key, val] = textread( landmarksFile, '%s:%[^\n]' );
%    for iI = 1:size( key, 1 )
%        landmarks_pnt( 1:3, iI, n ) = str2num( val{iI} );
%        minDistPnt = 1;
%        p1_x = landmarks_pnt(1,iI,n); p1_y = landmarks_pnt(2,iI,n); p1_z = landmarks_pnt(3,iI,n);
%        p2_x = A(1,1,n); p2_y = A(2,1,n); p2_z = A(3,1,n);
%        minDist = sqrt( (p1_x - p2_x)^2 + (p1_y - p2_y)^2 + (p1_z - p2_z)^2 );
%        for k = 2:numP(n)
%            p2_x = A(1,k,n); p2_y = A(2,k,n); p2_z = A(3,k,n);
%            tmpDist = sqrt( (p1_x - p2_x)^2 + (p1_y - p2_y)^2 + (p1_z - p2_z)^2 );
%            if tmpDist < minDist
%                minDist = tmpDist;
%                minDistPnt = k;
%            end
%        end
%        Landmarks(n, iI) = dist( minDistPnt, n);
%    end
%end

% read landmarks id on centerline
for n = 1:numf
    fidLandmarksIdFile = fopen( landmarksIdOnCenterlineFile, 'r' );
    tline = fgets( fidLandmarksIdFile );
    nLandmarksId = sscanf( tline, '%d' );
    for iI = 1:nLandmarksId
        tline = fgets( fidLandmarksIdFile );
		landmarksIdOnCenterline = sscanf( tline, '%d' );
        Landmarks(n, iI) = dist( landmarksIdOnCenterline, n);
    end
end

% read area 
for n = 1:numf
    fid = fopen(areaFile, 'r');
    tline = fgets(fid);
    nNumArea = sscanf(tline, '%d');
    for i = 1:nNumArea
        tline = fgets(fid);
        if(tline == -1) break; end
        area(i,n) = sscanf(tline, '%f');
    end
    fclose(fid);
end

% read area of ellipse
for n = 1:numf
    fid = fopen(ellipseAreaFile, 'r');
    if fid == -1
	continue;
    end
    tline = fgets(fid);
    nNumArea = sscanf(tline, '%d');
    for i = 1:nNumArea
        tline = fgets(fid);
        if(tline == -1) break; end
        areaEllipse(i,n) = sscanf(tline, '%f');
    end
    fclose(fid);
end

% read perimeter 
for n = 1:numf
    fid = fopen(perimeterFile, 'r');
    tline = fgets(fid);
    nNumPerimeter = sscanf(tline, '%d');
    for i = 1:nNumPerimeter
        tline = fgets(fid);
        if(tline == -1) break; end
        perimeter(i,n) = sscanf(tline, '%f');
    end
    fclose(fid);
end

addpath ('./Matlabfunctions/fdaM')
%  set up a basis for the functions W(t) that define the warping
%  functions
agefine = linspace(0,1,501)';
for n = 1:numf
    rng      = [0,1];
    knots    = dist(1:numP(n),n)';
    norder   = 6;
    nbasis   = length(knots) + norder - 2;
    areabasis = create_bspline_basis(rng, nbasis, norder, knots);
    Lfdobj   = int2Lfd(2);
    lambda   = 0.5*1e-6;
    areafdPar = fdPar(areabasis, Lfdobj, lambda);

    xValue = dist(1:numP(n), n);
    xValue(1) = xValue(1) + 100 * 1e-7;
    xValue(end) = xValue(end) - 100 * 1e-7;
    yValue = area(1:numP(n), n);
    for iTmp = 1:100
        xValue = [ xValue(1) - 1e-7; xValue ; xValue(end) + 1e-7 ];
        yValue = [ yValue(1); yValue; yValue(end) ];
    end
    areafd = smooth_basis( xValue, yValue, areafdPar );
    
    clear yValue
    yValue = areaEllipse(1:numP(n), n);
    for iTmp = 1:100
        yValue = [ yValue(1); yValue; yValue(end) ];
    end 
    areaEllipsefd = smooth_basis( xValue, yValue, areafdPar );
	
	clear yValue
	yValue = perimeter(1:numP(n), n);
    for iTmp = 1:100
        yValue = [ yValue(1); yValue; yValue(end) ];
    end 
	perimeterfd = smooth_basis(xValue, yValue, areafdPar );

	clear yValue 
	yValue = 4 .* area(1:numP(n), n) ./ perimeter(1:numP(n), n);
    for iTmp = 1:100
        yValue = [ yValue(1); yValue; yValue(end) ];
    end 
	hydraulicDiameterfd = smooth_basis(xValue, yValue, areafdPar );
    
    areafine(:,n) = eval_fd(agefine, areafd);
    areaEllipsefine(:,n) = eval_fd(agefine, areaEllipsefd);
	perimeterfine(:,n) = eval_fd(agefine, perimeterfd);
	hydraulicDiameter(:,n) = eval_fd(agefine, hydraulicDiameterfd);
end

rng_new = [0,1];
knots_new = agefine;
norder_new = 6;
nbasis_new = length(knots_new) + norder_new - 2;
basis_new = create_bspline_basis(rng_new, nbasis_new, norder_new, knots_new);
Lfdobj   = int2Lfd(2);
lambda   = 0.5*1e-6;
fdPar_new = fdPar(basis_new, Lfdobj, lambda);
areafd_new = smooth_basis(knots_new, areafine, fdPar_new);
areaEllipsefd_new = smooth_basis(knots_new, areaEllipsefine, fdPar_new);
perimeterfd_new = smooth_basis(knots_new, perimeterfine, fdPar_new);
hydraulicDiameterfd_new = smooth_basis(knots_new, hydraulicDiameter, fdPar_new);

for icase = 1:numf
    valueLM(icase, :) = eval_fd( Landmarks(icase, :), areafd_new(icase));
    valueLM_ellipse(icase, :) = eval_fd( Landmarks(icase, :), areaEllipsefd_new(icase) );
	valueLM_perimeter(icase, :) = eval_fd( Landmarks(icase, :), perimeterfd_new(icase) );
	valueLM_hydraulicDiameter(icase, :) = eval_fd( Landmarks(icase, :), hydraulicDiameterfd_new(icase) );
end

% draw curve and landmarks
style = {'bx', 'gx', 'rx', 'co'};
areafine_new = eval_fd(agefine, areafd_new);
figure, h = plot(agefine, areafine_new, 'LineWidth', 2);
hold on
% plot(agefine, mean(areafine_new, 2), 'm--', 'LineWidth', 3);
for icase = 1:numf
    plot(Landmarks(icase, :), valueLM(icase, :), char(style(icase)), 'LineWidth', 10);
    plot( dist(1:numP(icase), icase), area(1:numP(icase), icase), '.' );
end
hold off
xlabel( 'Depth along the centerline', 'FontSize', 24 );
ylabel( 'Cross-sectional area(mm^2)', 'FontSize', 24 );
ylim( [0 1000] );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
saveas( h, curveAndLandmarksFile );

filename = sprintf( '%s/curveWithLandmarks.png', outputPrefix );
saveas( h, filename );

% draw curve of perimeter
perimeterfile_new = eval_fd(agefine, perimeterfd_new);
figure, plot(agefine, perimeterfile_new, 'LineWidth', 2);
hold on
for icase = 1:numf
    plot(Landmarks(icase, :), valueLM_perimeter(icase, :), char(style(icase)), 'LineWidth', 10);
    plot( dist(1:numP(icase), icase), perimeter(1:numP(icase), icase), '.' );
end
hold off
xlabel( 'Depth along the centerline', 'FontSize', 24 );
ylabel( 'Perimeter of the cross section(mm)', 'FontSize', 24 );
filename = sprintf( '%s/curveOfPerimeter.png', outputPrefix );
saveas( gca, filename );

% draw curve of hydraulic diameter
hydraulicDiameter_new = eval_fd(agefine, hydraulicDiameterfd_new);
figure, plot(agefine, hydraulicDiameter_new, 'LineWidth', 2);
hold on
for icase = 1:numf
    plot(Landmarks(icase, :), valueLM_hydraulicDiameter(icase, :), char(style(icase)), 'LineWidth', 10);
    plot( dist(1:numP(icase), icase), 4.*area(1:numP(icase), icase)./perimeter(1:numP(icase), icase), '.' );
end
hold off
xlabel( 'Depth along the centerline', 'FontSize', 24 );
ylabel( 'Hydraulic diameter (mm)', 'FontSize', 24 );
filename = sprintf( '%s/curveOfHydraulicDiameter.png', outputPrefix );
saveas( gca, filename );

% draw curve of ellipse area
areaEllipsefine_new = eval_fd(agefine, areaEllipsefd_new);
figure, h = plot(agefine, areaEllipsefine_new, 'LineWidth', 2);
hold on
% plot(agefine, mean(areafine_new, 2), 'm--', 'LineWidth', 3);
for icase = 1:numf
    plot(Landmarks(icase, :), valueLM_ellipse(icase, :), char(style(icase)), 'LineWidth', 10);
    plot( dist(1:numP(icase), icase), areaEllipse(1:numP(icase), icase), '.' );
end
hold off
xlabel( 'Position' );
ylabel( 'Area' );
ylim( [0 1000] );
filename = sprintf( '%s/areaCurveWithLandmarksEllipse.png', outputPrefix );
saveas( h, filename );
