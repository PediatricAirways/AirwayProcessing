function drawLandmarksOnMultiCurvesAndSGS( meanNormFile, areaFile, ellipseAreaFile, landmarksIdOnCenterlineFile, curveAndLandmarksFile, outputPrefix, cases, ages, weights, gaussian_sigma, uniform_winSize, pos_SGS )

% meanNormFile: the centerline for each subject
% areaFile: cross-sectional areas for each subject
% ellipseAreaFile: the area of ellipse fitting cross sections for each subject
% landmarksIdOnCenterlineFile: the landmarks position on the centerline
% curveAndLandmarksFile: output file of plotting area curve and landmarks
% outputPrefix: the prefix for outputting files
% cases: subject id 
% ages: age for each subject
% weights: weight for each subject
% gaussian_sigma: the sigma for gaussian kernel
% uniform_winsize: the half width for the uniform window
% pos_SGS: 1-pos_SGS-1 are the control subject, pos_SGS-end are the SGS subjects

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: read files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read the mean points and compute dist in 1D, from 0 to 1
numf = max( size( meanNormFile ) );
numP = zeros(numf, 1);
for n = 1:numf
    fid = fopen( meanNormFile{n}, 'r');
    meanNormFile{n}
    tline = fgets(fid);
    count = 0;
    while ~feof(fid)
        tline = fgets(fid);
        count = count + 1;
        tmp = sscanf(tline, '%f');
        A(:, count, n) = tmp(1:3);
	normCases(:, count, n) = tmp(4:6);
    end
    fclose(fid);
    numP(n) = count;
    dist(1,n) = 0;
    for i=2:count
        p1_x = A(1,i,n); p1_y = A(2,i,n); p1_z = A(3,i,n);
        p2_x = A(1,i-1,n); p2_y = A(2,i-1,n); p2_z = A(3,i-1,n);
        dist(i,n) = sqrt( (p1_x-p2_x)^2 + (p1_y-p2_y)^2 + (p1_z-p2_z)^2 ) + dist(i-1,n);
    end
    dist_copy( 1:count, n ) = dist( 1:count, n );
    dist(1:count,n) = dist(1:count,n) ./ dist(count,n);
end

% read landmarks file
%for n = 1:numf
%    [key, val] = textread( landmarksFile{n}, '%s:%[^\n]' );
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
    %landmarksIdOnCenterlineFile{n}
    fidLandmarksIdFile = fopen( landmarksIdOnCenterlineFile{n}, 'r' );
    tline = fgets( fidLandmarksIdFile );
    nLandmarksId = sscanf( tline, '%d' );
    nInc = 0;
    for iI = 1:nLandmarksId
        tline = fgets( fidLandmarksIdFile );
        landmarksIdOnCenterline = sscanf( tline, '%d' );
        % remove the subglottic landmark
	if nLandmarksId > 5 && iI == 5 
		nInc = nInc + 1;
		continue;
	end
 	LandmarksId( n, iI-nInc ) = landmarksIdOnCenterline;
        Landmarks(n, iI-nInc) = dist( landmarksIdOnCenterline, n);
    end
end

% make the landmarks increase with its corresponding value
%LandmarksId = sort( LandmarksId, 2 );
%Landmarks = sort( Landmarks, 2 );


% output the length from TVC to trachea carina
trachea_length = zeros(numf, 1);
for n = 1:numf
	trachea_length( n ) = dist_copy( LandmarksId( n, end ), n ) - dist_copy( LandmarksId( n, end-1 ), n );
end
filename_tmp = sprintf( '%s/trachea_length.mat', outputPrefix );
save( filename_tmp, 'trachea_length' );

% read area file
for n = 1:numf
    fid = fopen(areaFile{n}, 'r');
    tline = fgets(fid);
    nNumArea = sscanf(tline, '%d');
    for i = 1:nNumArea
        tline = fgets(fid);
        if(tline == -1) break; end
        area(i,n) = sscanf(tline, '%f');
    end
    fclose(fid);
end

% write out the area and normal for landmarks of each subject
%{
fid_landmarks_area = fopen( 'landmarks_area.txt', 'wt' );
fid_landmarks_normal = fopen( 'landmarks_normal.txt', 'wt' );
for n = 1:numf
    fprintf( fid_landmarks_area, '%s : ', cases{n} );
    fprintf( fid_landmarks_normal, '%s : ', cases{n} );
    for nLM = 1:size( LandmarksId, 2 )
        fprintf( fid_landmarks_area, '%f ', area( LandmarksId( n, nLM ), n ) );
        fprintf( fid_landmarks_normal, '(%f, %f, %f) ', normCases( 1, LandmarksId(n, nLM), n ), normCases( 2, LandmarksId(n, nLM), n ), normCases( 3, LandmarksId(n, nLM),  n ) );
    end
    fprintf( fid_landmarks_area, '\n' );
    fprintf( fid_landmarks_normal, '\n' );
end
fclose( fid_landmarks_area );
fclose( fid_landmarks_normal );
%}

% read area of ellipse if needed
for n = 1:numf
    %ellipseAreaFile{n}
    fid = fopen(ellipseAreaFile{n}, 'r');
    tline = fgets(fid);
    nNumArea = sscanf(tline, '%d');
    for i = 1:nNumArea
        tline = fgets(fid);
        if(tline == -1) break; end
        areaEllipse(i,n) = sscanf(tline, '%f');
    end
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curve representation and registration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ('./Matlabfunctions/fdaM')
addpath ('./FunctionalBoxplot')
%  set up a basis for the functions W(t) that define the warping functions
%  create curves to fit the area of contour and of ellipse
for n = 1:numf
    rng      = [0,1];
    knots    = dist(1:numP(n),n)';
    norder   = 6;
    nbasis   = length(knots) + norder - 2;
    areabasis = create_bspline_basis(rng, nbasis, norder, knots);
    Lfdobj   = int2Lfd(2);
    lambda   = 1e-5;  % adjust this number to make the curve fit the scattered point
    areafdPar = fdPar(areabasis, Lfdobj, lambda);

    % curve fitting function for the area
    % try to make the first and last points on the curve.
    xValue = dist(1:numP(n), n);
    xValue(1) = xValue(1) + 100 * 1e-7;
    xValue(end) = xValue(end) - 100 * 1e-7;
    yValue = area(1:numP(n), n);
    for iTmp = 1:100
        xValue = [ xValue(1) - 1e-7; xValue ; xValue(end) + 1e-7 ];
        yValue = [ yValue(1); yValue; yValue(end) ];
    end
    %areafd = smooth_basis(dist(1:numP(n),n), area(1:numP(n),n), areafdPar);
    areafd_tmp = smooth_basis( xValue, yValue, areafdPar );
    % to make the nasal spine to TVC smooth, and TVC to trachea carina fit
    % points well
    areaSmooth(1:LandmarksId(n, 4)-1, n) = eval_fd(dist(1:LandmarksId(n, 4)-1, n), areafd_tmp );
    areaSmooth(LandmarksId(n, 4):numP(n), n) = area(LandmarksId(n, 4):numP(n), n);
end

agefine = linspace(0,1,501)';
for n = 1:numf
    rng      = [0,1];
    knots    = dist(1:numP(n),n)';
    norder   = 6;
    nbasis   = length(knots) + norder - 2;
    areabasis = create_bspline_basis(rng, nbasis, norder, knots);
    Lfdobj   = int2Lfd(2);
    lambda   = 0.5 * 1e-6;  % adjust this number to make the curve fit the scattered point
    areafdPar = fdPar(areabasis, Lfdobj, lambda);

    % curve fitting function for the area
    % try to make the first and last points on the curve.
    xValue = dist(1:numP(n), n);
    xValue(1) = xValue(1) + 100 * 1e-7;
    xValue(end) = xValue(end) - 100 * 1e-7;
    yValue = areaSmooth(1:numP(n), n);
    for iTmp = 1:100
        xValue = [ xValue(1) - 1e-7; xValue ; xValue(end) + 1e-7 ];
        yValue = [ yValue(1); yValue; yValue(end) ];
    end
    %areafd = smooth_basis(dist(1:numP(n),n), area(1:numP(n),n), areafdPar);
    areafd = smooth_basis( xValue, yValue, areafdPar );
    
    clear yValue
    yValue = areaEllipse(1:numP(n), n);
    for iTmp = 1:100
        yValue = [ yValue(1); yValue; yValue(end) ];
    end    
    %areaEllipsefd = smooth_basis(dist(1:numP(n),n), areaEllipse(1:numP(n),n), areafdPar);
    areaEllipsefd = smooth_basis( xValue, yValue, areafdPar );
    
%     figure, plotfit_fd(area(1:numP(n),n), dist(1:numP(n),n), areafd);
%     hold on
%     plotfit_fd(areaEllipse(1:numP(n),n), dist(1:numP(n),n), areaEllipsefd);
%     hold off
%     legend( 'AreaDots', 'AreaCurve', 'EllipseDots', 'EllipseCurve' );
%     filename_area = sprintf('%s/areaMultiPlot%02d.png', outputPrefix, n);
%     saveas( gcf, filename_area );
    areafine(:,n) = eval_fd(agefine, areafd);
    areaEllipsefine(:,n) = eval_fd(agefine, areaEllipsefd);
end

% use the same number of points to represent the curves for alll subjects
rng_new = [0,1];
knots_new = agefine;
norder_new = 6;
nbasis_new = length(knots_new) + norder_new - 2;
basis_new = create_bspline_basis(rng_new, nbasis_new, norder_new, knots_new);
Lfdobj   = int2Lfd(2);
lambda   = 0.5 * 1e-6;
fdPar_new = fdPar(basis_new, Lfdobj, lambda);
areafd_new = smooth_basis(knots_new, areafine, fdPar_new);
areaEllipsefd_new = smooth_basis(knots_new, areaEllipsefine, fdPar_new);

% find the value for landmarks on curves
for icase = 1:numf
    valueLM(icase, :) = eval_fd( Landmarks(icase, :), areafd_new(icase));
    valueLM_ellipse(icase, :) = eval_fd( Landmarks(icase, :), areaEllipsefd_new(icase) );
end

% write the area of Landmarks on curve
%{
fid_landmarks_area_onCurve = fopen( 'landmarks_area_on_curve.txt', 'wt' );
for icase = 1:numf
    fprintf( fid_landmarks_area_onCurve, '%s : ', cases{icase} );
    for nLM = 1:size( valueLM, 2 )
        fprintf( fid_landmarks_area_onCurve, '%f ', valueLM( icase, nLM ) );
    end
    fprintf( fid_landmarks_area_onCurve, '\n' );
end
fclose( fid_landmarks_area_onCurve );
%}

% draw curve and landmarks
style = {'yx', 'rx', 'gx', 'bx', 'cx'};
areafine_new = eval_fd(agefine, areafd_new);
%figure, plot(agefine, areafine_new, 'LineWidth', 2);
%hold on
% draw the mean of the curve if needed
%plot(agefine, mean(areafine_new, 2), 'm--', 'LineWidth', 3);
%maxAge = max( ages );
%minAge = min( ages );
%cmap = hot(64);
%cmap = autumn(64);

% the color to display different ages
nBins = ceil( max(ages)/100.0 ) * 100;
cmap = zeros( nBins, 3 );
stepColor = 1.0/(nBins/6.0);
for iI = 1:nBins
    if iI < (nBins/12.0)
        cmap( iI, 1 ) = 1;
        cmap( iI, 2 ) = 0.5 - iI * stepColor;
        cmap( iI, 3 ) = 1;
    elseif iI < (nBins/12.0*3)
		cmap( iI, 1 ) = 1 - (iI-nBins/12.0) * stepColor;
		cmap( iI, 2 ) = 0;
		cmap( iI, 3 ) = 1;
	elseif iI < (nBins/12.0*6)
		cmap( iI, 1 ) = 0;
		cmap( iI, 2 ) = (iI - nBins/12.0*3) * stepColor * 2/3;
		cmap( iI, 3 ) = 1;
	elseif iI < (nBins/12.0*7)
		cmap( iI, 1 ) = 0;
		cmap( iI, 2 ) = 1;
		cmap( iI, 3 ) = 1 - (iI - nBins/12.0*6.0) * stepColor * 2;
    elseif iI < (nBins/12.0*8)
        cmap( iI, 1 ) = (iI - nBins/12.0*7.0) * stepColor * 2;
        cmap( iI, 2 ) = 1;
        cmap( iI, 3 ) = 0;
    elseif iI < (nBins/12.0*11)
		cmap( iI, 1 ) = 1;
		cmap( iI, 2 ) = 1 - (iI - nBins/12.0*8.0) * stepColor * 2 / 3;
		cmap( iI, 3 ) = 0;
    else
        cmap( iI, 1 ) = 1 - (iI - nBins/12.0*11.0) * stepColor;
        cmap( iI, 2 ) = 0;
        cmap( iI, 3 ) = 0;
	end
end
%cmap = jet(nBins);
colormap(cmap);

% normalized weight for coloring
minWeight = 0;
maxWeight = ceil(max( weights ));
colorPosW = round( ( weights - minWeight ) / ( maxWeight - minWeight + 1e-10 ) * nBins + 1 );
colorPosW = min( colorPosW, nBins );
colorPosW = max( colorPosW, 1 );

% display all the unregistered curves
figure('Colormap', cmap), hold on
for icase = 1:numf
    if icase < pos_SGS
    	plot( agefine, areafine_new( :, icase ), 'Color', cmap(round(ages(icase)), :), 'LineWidth', 1 );
    else
	plot( agefine, areafine_new( :, icase ), 'Color', cmap(round(ages(icase)), :), 'LineWidth', 2, 'LineStyle', '--' );
    end
end
% display landmarks on curves
for icase = 1:size(Landmarks, 2)
     plot( Landmarks(:, icase), valueLM(:, icase), char(style(mod(icase, 5)+1)), 'LineWidth', 3);
end
hold off
caxis( [0 nBins] );
h = colorbar( 'peer', gca );
set( get( h, 'ylabel' ), 'String', 'Age: month(s)', 'FontSize', 20, 'FontWeight', 'Bold' );
set( h, 'ylim', [0, nBins] );
set( gca, 'XTick', 0:0.1:1 );
xlabel( 'Depth along the centerline', 'FontSize', 20, 'FontWeight', 'Bold' );
ylabel( 'Cross-sectional area (mm^2)', 'FontSize', 20, 'FontWeight', 'Bold' );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
title( 'Unregistered Curves', 'FontSize', 20, 'FontWeight', 'Bold' );
%set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
%{
set( gcf, 'PaperUnits', 'inches' );
set( gcf, 'PaperSize', [11 8.5] );
set( gcf, 'PaperPositionMode', 'manual' );
set( gcf, 'PaperPosition', [0 0 11 8.5] );
set( gcf, 'renderer', 'painters' );
print( gcf, '-dpng', curveAndLandmarksFile );
%}
saveas( gca, curveAndLandmarksFile );

% display all the unregistered curves for ellipse's area
areaEllipsefine_new = eval_fd(agefine, areaEllipsefd_new);
figure('Colormap', cmap), hold on
for icase = 1:numf
    plot( agefine, areaEllipsefine_new( :, icase ), 'Color', cmap(round(ages(icase)), :), 'LineWidth', 1 );
end
for icase = 1:size(Landmarks, 2)
    plot( Landmarks(:, icase), valueLM_ellipse(:, icase), char(style(mod(icase, 5)+1)), 'LineWidth', 3);
end
hold off
caxis( [0 nBins] );
h = colorbar;
set( h, 'ylim', [0, nBins] );
xlabel( 'Depth along the centerline', 'FontSize', 20, 'FontWeight', 'Bold' );
ylabel( 'Ellipse area', 'FontSize', 20, 'FontWeight', 'Bold' );
set( gca, 'XTick', 0:0.1:1 );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
title( 'Unregistered Curves  (Ellipse)', 'FontSize', 20, 'FontWeight', 'Bold' );
filename = sprintf( '%s/areaMultiCurveWithLandmarksEllipse.png', outputPrefix );
%set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
set( gcf, 'PaperUnits', 'inches' );
set( gcf, 'PaperSize', [11 8.5] );
set( gcf, 'PaperPositionMode', 'manual' );
set( gcf, 'PaperPosition', [0 0 11 8.5] );
set( gcf, 'renderer', 'painters' );
print( gcf, '-dpng', filename );
%saveas( gca, filename );



% do the rigistration with landmarks
Landmarks( :, 1 ) = Landmarks( :, 1 ) + 1e-7;
Landmarks( :, end ) = Landmarks( :, end ) - 1e-7;
nLandmarks = size( Landmarks, 2 );
 
wbasisLM = create_bspline_basis([0,1], max(nLandmarks+3-2,4), 3);
WfdLM    = fd(zeros(max(nLandmarks+3-2,4),1),wbasisLM);
WfdParLM = fdPar(WfdLM,1,1e-12);
 
%  carry out the landmark registration
LandmarksMean = mean(Landmarks);
[areafdLM, areawarpfdLM, WfdLM] = ...
        landmarkreg(areafd_new, Landmarks, LandmarksMean, WfdParLM, 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot registeration results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% %  plot registered accelerations along with warping functions
areamatUR = eval_fd(agefine, areafd_new);
areamatLM = eval_fd(agefine, areafdLM);
areawarpmatLM  = eval_fd(agefine, areawarpfdLM);
areawarpmatLM(1,:) = 0; areawarpmatLM(length(agefine),:) = 1;

figure('Colormap', cmap), hold on
for icase = 1:numf
    plot( agefine, areawarpmatLM( :, icase ), 'Color', cmap(round(ages(icase)), :), 'LineWidth', 1 );
end
for icase = 1:numf
    for iLM = 1:size( Landmarks, 2 )
	plot(mean(Landmarks(:,iLM)), Landmarks(icase, iLM), char(style(mod(iLM, 5)+1)), 'LineWidth', 3);
    end
end
ylabel( 'Physical position of the landmark', 'FontSize', 20, 'FontWeight', 'Bold' );
xlabel( 'Mean position of the landmark', 'FontSize', 20, 'FontWeight', 'Bold' );
title( 'Warping function', 'FontSize', 20, 'FontWeight', 'Bold' );
set( gca, 'XTick', 0:0.1:1 );
set( gca, 'YTick', 0:0.1:1 );
caxis( [0 nBins] );
h = colorbar;
set( get( h, 'ylabel' ), 'String', 'Age: month(s)', 'FontSize', 20, 'FontWeight', 'Bold' );
set( h, 'ylim', [0, nBins] );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
filename = sprintf( '%s/warpFunctionLandmarks.png', outputPrefix );
%set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
%{
set( gcf, 'PaperUnits', 'inches' );
set( gcf, 'PaperSize', [11 8.5] );
set( gcf, 'PaperPositionMode', 'manual' );
set( gcf, 'PaperPosition', [0 0 11 8.5] );
set( gcf, 'renderer', 'painters' );
print( gcf, '-dpng', filename );
%}
saveas( gca, filename );
%filename = sprintf( '%s/warpFunctionLandmarks.fig', outputPrefix );
%saveas( gca, filename );



% display all the registered curves, colored by age
figure('Colormap', cmap); hold on
for icase = pos_SGS:size( areamatLM, 2 )
	plot( agefine, areamatLM( :, icase ), 'Color', cmap(round(ages(icase)), :), 'LineWidth', 2, 'LineStyle', '-.' );
	hold on
end
for icase = 1:numf
    if icase < pos_SGS
	plot( agefine, areamatLM( :, icase ), 'Color', cmap(round(ages(icase)), :), 'LineWidth', 1 );
    else
	plot( agefine, areamatLM( :, icase ), 'Color', cmap(round(ages(icase)), :), 'LineWidth', 2, 'LineStyle', '-.' );
    end
end
for icase = 1:numf
    warpInvfd = smooth_basis(areawarpmatLM(:,icase), agefine, fdPar_new);
    warpedLM(icase,:) = eval_fd(Landmarks(icase,:), warpInvfd);
    for iLM = 1:size(warpedLM, 2)
	warpedLM(icase, iLM) = max( warpedLM(icase, iLM), 0 );
	warpedLM(icase, iLM) = min( warpedLM(icase, iLM), 1 );
    end
    valueLM(icase, :) = eval_fd(warpedLM(icase,:), areafdLM(icase));
    for iLM = 1:size(warpedLM, 2)
    	plot(warpedLM(icase, iLM), valueLM(icase, iLM), char(style(mod(iLM, 5)+1)), 'LineWidth', 3);
    end
end
hold off
caxis( [0 nBins] );
h = colorbar;
set( get( h, 'ylabel' ), 'String', 'Age: month(s)', 'FontSize', 20, 'FontWeight', 'Bold' );
set( h, 'ylim', [0 nBins] );
xlabel( 'Depth along the centerline', 'FontSize', 20, 'FontWeight', 'Bold' );
ylabel( 'Cross-sectional area (mm^2)', 'FontSize', 20, 'FontWeight', 'Bold' );
if pos_SGS <= size( areamatLM, 2 )
	legend( cases(pos_SGS:end) );
end
set( gca, 'XTick', 0:0.1:1 );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
title( 'Registered Curves', 'FontSize', 20, 'FontWeight', 'Bold' );
filename = sprintf( '%s/landmark_registration_age.png', outputPrefix );
%set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
%{
set( gcf, 'PaperUnits', 'inches' );
set( gcf, 'PaperSize', [11 8.5] );
set( gcf, 'PaperPositionMode', 'manual' );
set( gcf, 'PaperPosition', [0 0 11 8.5] );
set( gcf, 'renderer', 'painters' );
print( gcf, '-dpng', filename );
%}
saveas( gca, filename );
filename = sprintf( '%s/landmark_registration_age.fig', outputPrefix );
saveas( gca, filename );

% display all the registered curves, colored by weight
minWeight = min( weights );
maxWeight = max( weights );
figure('Colormap', cmap), hold on
for icase = 1:numf
    if icase < pos_SGS
    	plot( agefine, areamatLM( :, icase ), 'Color', cmap(colorPosW(icase), :), 'LineWidth', 1 );
    else
	plot( agefine, areamatLM( :, icase ), 'Color', cmap(colorPosW(icase), :), 'LineWidth', 2, 'LineStyle', '--' );
    end
end
%plot(agefine, mean(areamatLM, 2), 'm--', 'LineWidth', 3);
for icase = 1:numf
    for iLM = 1:size(warpedLM, 2)
        plot(warpedLM(icase, iLM), valueLM(icase, iLM), char(style(mod(iLM, 5)+1)), 'LineWidth', 3);
    end
end
hold off
caxis( [minWeight maxWeight] );
h = colorbar
set( h, 'ylim', [minWeight, maxWeight] );
xlabel( 'Depth along the centerline', 'FontSize', 20, 'FontWeight', 'Bold' );
ylabel( 'Cross-sectional area', 'FontSize', 20, 'FontWeight', 'Bold' );
set( gca, 'XTick', 0:0.1:1 );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
title( 'Registered Curves (purple-red : light-heavy, kg)', 'FontSize', 20, 'FontWeight', 'Bold' );
filename = sprintf( '%s/landmark_registration_weight.png', outputPrefix );
%set(gca, 'FontSize', 20 );
%set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
set( gcf, 'PaperUnits', 'inches' );
set( gcf, 'PaperSize', [11 8.5] );
set( gcf, 'PaperPositionMode', 'manual' );
set( gcf, 'PaperPosition', [0 0 11 8.5] );
set( gcf, 'renderer', 'painters' );
print( gcf, '-dpng', filename );
%saveas( gca, filename );

filename = sprintf( '%s/landmark_registration_weight.fig', outputPrefix );
saveas( gca, filename );


% display the mean curve and variations
figure, plot( agefine, mean(areamatLM(:,1:pos_SGS-1), 2), 'm--', 'LineWidth', 3 );
hold on
plot( agefine, max(areamatLM(:,1:pos_SGS-1), [], 2), 'b-', 'LineWidth', 2 );
plot( agefine, min(areamatLM(:,1:pos_SGS-1), [], 2), 'b-', 'LineWidth', 2 );
for iLM = 1:size(warpedLM, 2)
    plot( mean( warpedLM(1:pos_SGS-1, iLM) ), mean( valueLM(1:pos_SGS-1, iLM) ), char(style(mod(iLM, 5)+1)), 'LineWidth', 10);
end
hold off
xlabel( 'Depth along the centerline', 'FontSize', 20, 'FontWeight', 'Bold' );
ylabel( 'Cross-sectional area', 'FontSize', 20, 'FontWeight', 'Bold' );
set( h, 'ylim', [0, 1000] );
set( gca, 'XTick', 0:0.1:1 );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
title( 'Mean Curve with Variations', 'FontSize', 20, 'FontWeight', 'Bold' );
%{
[valueMin, idMin] = min(min(areamatLM, [], 2));
tmpstr = sprintf( '\\leftarrow minimal:(%f, %f)', agefine(idMin), valueMin );
text( agefine(idMin), valueMin, tmpstr, 'HorizontalAlignment', 'left' );
[valueMax, idMax] = max(max(areamatLM, [], 2));
tmpstr = sprintf( '\\leftarrow maximal:(%f, %f)', agefine(idMax), valueMax );
text( agefine(idMax), valueMax, tmpstr, 'HorizontalAlignment', 'left' );
%}
filename = sprintf( '%s/landmark_registration_mean_variation.png', outputPrefix );
set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
set( gcf, 'PaperUnits', 'inches' );
set( gcf, 'PaperSize', [11 8.5] );
set( gcf, 'PaperPositionMode', 'manual' );
set( gcf, 'PaperPosition', [0 0 11 8.5] );
set( gcf, 'renderer', 'painters' );
print( gcf, '-dpng', filename );
%saveas( gca, filename );
%filename = sprintf( '%s/landmark_registration_mean_variation.fig', outputPrefix );
%saveas( gca, filename );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistical analysis: functional boxplot, pca, point-wise boxplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default settings
factor = 1.5;
fullout = false;
barcol = 'b';
outliercol = 'r';
%color = [ [0.8784, 0.6902, 1.0]; [1, 0, 1]; [0.5451, 0, 0.8] ];
%prob = [0.75, 0.5, 0.25];
color = 'm';
prob = 0.5;
show = true;
method = 'MBD';
depth = [];

% display the functional boxplot for all the curves
figure;
for icase = pos_SGS : size( areamatLM, 2 )
	plot( agefine, areamatLM(:, icase), 'LineWidth', 2, 'Color', cmap(round(ages(icase)), :), 'LineStyle', '-.' );
 	hold on
end
if pos_SGS <= size( areamatLM, 2 )
	legend( cases(pos_SGS:end), 'FontSize', 20, 'FontWeight', 'Bold' );
end
fbplot( areamatLM(:, 1:pos_SGS-1), agefine, depth, method, show, prob, color, outliercol, barcol, fullout, factor );
hold on
for icase = pos_SGS : size( areamatLM, 2 )
	plot( agefine, areamatLM(:, icase), 'LineWidth', 2, 'Color', cmap(round(ages(icase)), :), 'LineStyle', '-.' );
end
hold off
filename = sprintf( '%s/boxplot/functional_boxplot_allCurves_%s.png', outputPrefix, method );
set( gca, 'ylim', [0, 1000] );
title( 'Functional boxplot for all curves' );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
%set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
set( gcf, 'PaperUnits', 'inches' );
set( gcf, 'PaperSize', [11 8.5] );
set( gcf, 'PaperPositionMode', 'manual' );
set( gcf, 'PaperPosition', [0 0 11 8.5] );
set( gcf, 'renderer', 'painters' );
print( gcf, '-dpng', filename );
%saveas( gca, filename );

% display the point-wise boxplot for all the curves
%tmpLinePoint = 1:length(agefine);
%figure, boxplot( areamatLM( tmpLinePoint( logical( mod( tmpLinePoint, 30 ) == 0 ) ), :)' );
figure; 
for icase = pos_SGS : size( areamatLM, 2 )
	plot( 1:size(areamatLM, 1), areamatLM(:, icase), 'LineWidth', 2, 'Color', cmap(round(ages(icase)), :), 'LineStyle', '-.' );
	hold on
end
if pos_SGS <= size( areamatLM, 2 )
	legend( cases(pos_SGS:end), 'FontSize', 20, 'FontWeight', 'Bold');
end
boxplot( (areamatLM(:, 1:pos_SGS-1))' );
hold on
for icase = pos_SGS : size( areamatLM, 2 )
	plot( 1:size(areamatLM, 1), areamatLM(:, icase), 'LineWidth', 2, 'Color', cmap(round(ages(icase)), :), 'LineStyle', '-.' );
end
hold off
set( gca, 'ylim', [0, 1000] );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
title( 'Pointwise boxplot for all curves' );
filename = sprintf( '%s/boxplot/pointwise_boxplot_allCurves_allPoints.png', outputPrefix );
%set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
set( gcf, 'PaperUnits', 'inches' );
set( gcf, 'PaperSize', [11 8.5] );
set( gcf, 'PaperPositionMode', 'manual' );
set( gcf, 'PaperPosition', [0 0 11 8.5] );
set( gcf, 'renderer', 'painters' );
print( gcf, '-dpng', filename );
%saveas( gca, filename );


% PCA and boxplot performed on parts of the curves

% find the id for each landmark after registration
landmarksNewId = zeros(size(warpedLM,2), 1);
for iSample = 1:size(warpedLM,2)
    [nValue, nId] = min( abs( agefine - mean( warpedLM( :, iSample ) ) ) );
    landmarksNewId( iSample ) = nId;
end

% subsection of curves: from the current landmark to the end
part_sample = [1, size( warpedLM, 2 ) - 1 ];
for iSample = 1:length( part_sample )
	nStartLM = part_sample(iSample);
	[nStartValue, nStartId] = min( abs( agefine - mean( warpedLM( :, nStartLM ) ) ) );
	[nEndValue, nEndId] = min( abs( agefine - mean( warpedLM( :, nStartLM+1 ) ) ) );
	area_part = areamatLM( nStartId:end, : );
	agefine_part = agefine( nStartId:end ); 
	yLim = ceil( max( max(area_part) ) / 100.0 ) * 100;
	rng = [ min(agefine_part), max(agefine_part) ];
	knots = agefine_part;
	norder = 6;
	nbasis = length(knots) + norder - 2;
	areabasis_part = create_bspline_basis( rng, nbasis, norder, knots );
	Lfdobj = int2Lfd(2);
	lambda = 0.5 * 1e-6;       % adjust for fitting
	if part_sample == 1
		lambda = 0.5 * 1e-5;
	end
	areafdPar_part = fdPar( areabasis_part, Lfdobj, lambda );
	areafd_part = smooth_basis( agefine_part, area_part, areafdPar_part );
	areamatLM_part = eval_fd( agefine_part, areafd_part );

	% display the subsection of the curves
	figure, hold on
	for icase = pos_SGS : size( areamatLM_part, 2 )
		plot( agefine_part, areamatLM_part(:, icase), 'Color', cmap( round(ages(icase)), : ), 'LineWidth', 2, 'LineStyle', '-.' );
	end
	for icase = 1:numf
	    if icase < pos_SGS
	    	plot( agefine_part, areamatLM_part(:, icase), 'Color', cmap( round(ages(icase)), : ), 'LineWidth', 1 );
	    else
		plot( agefine_part, areamatLM_part(:, icase), 'Color', cmap( round(ages(icase)), : ), 'LineWidth', 1, 'LineStyle', '-.' );
	    end
	end
	hold off
	set( gca, 'xlim', [min(agefine_part), max(agefine_part)] );
	set( gca, 'ylim', [0, yLim] );
	if pos_SGS <= size( areamatLM_part, 2 )
		legend( cases{pos_SGS:end} );
	end
	set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
	filename = sprintf( '%s/PCA/Curve_Part_StartFromLandmark%d.png', outputPrefix, part_sample(iSample) );
	set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	set( gcf, 'PaperUnits', 'inches' );
	set( gcf, 'PaperSize', [11 8.5] );
	set( gcf, 'PaperPositionMode', 'manual' );
	set( gcf, 'PaperPosition', [0 0 11 8.5] );
	set( gcf, 'renderer', 'painters' );
	print( gcf, '-dpng', filename );
	%saveas( gca, filename );

	% principle component analysis
	nharm = 5;
	areafd_part_CRL = areafd_part(1:pos_SGS-1);
	areafdLM_PCA_part = pca_fd( areafd_part_CRL-mean(areafd_part_CRL), nharm );
	%areafdLM_PCA_part = pca_fd( areafd_part - mean( areafd_part ), nharm );
	areamatLM_PCA_part = eval_fd( agefine_part, areafdLM_PCA_part.harmfd );
	for iharm = 1:nharm
	    newmat1_part = areamatLM_PCA_part( :, iharm ) * sqrt( areafdLM_PCA_part.values( iharm ) ) + mean( areamatLM_part(:, 1:pos_SGS-1), 2 );
	    newmat2_part = areamatLM_PCA_part( :, iharm ) * -sqrt( areafdLM_PCA_part.values( iharm ) ) + mean( areamatLM_part(:, 1:pos_SGS-1), 2 );
	    figure, plot( agefine_part, mean( areamatLM_part(:, 1:pos_SGS-1), 2 ), 'r-', 'LineWidth', 2 );
	    hold on
	    plot( agefine_part, newmat1_part, 'g--', 'LineWidth', 2 );
	    plot( agefine_part, newmat2_part, 'b-.', 'LineWidth', 2 );
	    hold off
	    title(['Cross-sectional Area: Part PC', int2str(iharm), ' (', num2str(areafdLM_PCA_part.varprop(iharm)*100, '%10.2f'), '%)'], ...
		'FontSize', 20, 'FontWeight', 'Bold');
	    filename = sprintf( '%s/PCA/PC%02d_part.png', outputPrefix, iharm );
	    set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
	    set( gca, 'xlim', [min(agefine_part), max(agefine_part)] );
	    set( gca, 'ylim', [0, yLim] );
	    legend( 'Mean', 'Mean+\lambda*v', 'Mean-\lambda*v' );
	    set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	    set( gcf, 'PaperUnits', 'inches' );
	    set( gcf, 'PaperSize', [11 8.5] );
	    set( gcf, 'PaperPositionMode', 'manual' );
	    set( gcf, 'PaperPosition', [0 0 11 8.5] );
	    set( gcf, 'renderer', 'painters' );
	    print( gcf, '-dpng', filename );
	    %saveas( gca, filename );

	    % display all the curves with just one component
	    figure,
	    for icase = 1:size(areafdLM_PCA_part.harmscr, 1)
		curveComponent_part = areamatLM_PCA_part( :, iharm ) * areafdLM_PCA_part.harmscr( icase, iharm ) + mean( areamatLM_part(:, 1:pos_SGS-1), 2 );
		plot( agefine_part, curveComponent_part, 'Color', cmap(round(ages( icase )), :), 'LineWidth', 1 );
		hold on
	    end
	    hold off
	    title( ['Curve component (Part) ', int2str(iharm)], 'FontSize', 20, 'FontWeight', 'Bold' );
	    caxis( [0 nBins] );
	    h = colorbar
	    set( h, 'ylim', [0, nBins] );
	    set( gca, 'xlim', [min(agefine_part), max(agefine_part)] );
	    set( gca, 'ylim', [0, yLim] );
	    set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
	    filename = sprintf( '%s/PCA/CurveComponentPart%02d_StartFromLandmark%d.png', outputPrefix, iharm, part_sample(iSample) );
	    set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	    set( gcf, 'PaperUnits', 'inches' );
	    set( gcf, 'PaperSize', [11 8.5] );
	    set( gcf, 'PaperPositionMode', 'manual' );
	    set( gcf, 'PaperPosition', [0 0 11 8.5] );
	    set( gcf, 'renderer', 'painters' );
	    print( gcf, '-dpng', filename );
	    %saveas( gca, filename );

	    % display scores
	    figure, hold on
	    for icase = 1:size(areafdLM_PCA_part.harmscr, 1 )
		plot( areafdLM_PCA_part.harmscr( icase, iharm ), rand(1)/100, 'Color', cmap( round(ages( icase )), : ), 'Marker', 'o', 'MarkerSize', 6, 'LineWidth', 2 );
	    end
	    hold off
	    set( gca, 'ylim', [-0.05, 0.05] );
	    set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
	    title( 'Scores of part of the curves' );
	    filename = sprintf( '%s/PCA/Score_part%02d_StartFromLandmark%d.png', outputPrefix, iharm, part_sample(iSample) );
	    set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	    set( gcf, 'PaperUnits', 'inches' );
	    set( gcf, 'PaperSize', [11 8.5] );
	    set( gcf, 'PaperPositionMode', 'manual' );
	    set( gcf, 'PaperPosition', [0 0 11 8.5] );
	    set( gcf, 'renderer', 'painters' );
	    print( gcf, '-dpng', filename );
	    %saveas( gca, filename );
	end
	figure, hold on
	for icase = 1:size(areafdLM_PCA_part.harmscr, 1 )
	    plot( areafdLM_PCA_part.harmscr( icase, 1 ), areafdLM_PCA_part.harmscr( icase, 2 ), 'Color', cmap( round(ages( icase )), : ), 'Marker', 'o', 'MarkerSize', 6, 'LineWidth', 2 );
	end
	hold off
	filename = sprintf( '%s/PCA/Score1_Score2_part_StartFromLandmark%d.png', outputPrefix, part_sample(iSample) );
 	set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
	set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	set( gcf, 'PaperUnits', 'inches' );
	set( gcf, 'PaperSize', [11 8.5] );
	set( gcf, 'PaperPositionMode', 'manual' );
	set( gcf, 'PaperPosition', [0 0 11 8.5] );
	set( gcf, 'renderer', 'painters' );
	print( gcf, '-dpng', filename );
	%saveas( gca, filename );

        % display PCA results in one figure
	displayPCA(agefine_part, areafdLM_PCA_part, areamatLM_PCA_part, areamatLM_part, areafd_part, cmap, ages, outputPrefix, part_sample(iSample), pos_SGS);

  	% display boxplot results 
	outputPrefix_part = sprintf( '%s/Age/boxplotPart%d', outputPrefix, part_sample(iSample) );
	% boxplot for control subjects
	displayBoxplot( agefine_part, areamatLM_part(:, 1:pos_SGS-1), ages(1:pos_SGS-1), outputPrefix_part, yLim, gaussian_sigma, uniform_winSize, cases );
	% boxplot for SGS subjects
        displayBoxplotForSGS( agefine_part, areamatLM(landmarksNewId(part_sample(iSample)):end, :), ages, cases, outputPrefix_part, yLim, gaussian_sigma, uniform_winSize, cmap, pos_SGS, areamatLM, landmarksNewId, nStartLM );
	
	%outputPrefix_part = sprintf( '%s/Weight/boxplotPart%d', outputPrefix, part_sample(iSample) );
	%displayBoxplotForSGS( agefine_part, areamatLM(landmarksNewId(part_sample(iSample)):end, :), weights, cases, outputPrefix_part, yLim, gaussian_sigma/3, uniform_winSize/3, cmap, pos_SGS, areamatLM, landmarksNewId, nStartLM );
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction: plot results of pca
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displayPCA(agefine_part, areafdLM_PCA_part, areamatLM_PCA_part, areamatLM_part, areafd_part, cmap, ages, outputPrefix, partId, pos_SGS)
	close all
	if pos_SGS <= size(areamatLM_part, 2)
		areaCoef = getcoef( (areafd_part - mean(areafd_part(1:pos_SGS-1))) );
		sgsScores = (areaCoef(:,pos_SGS:end))' * (areafdLM_PCA_part.Lmat)' * (areafdLM_PCA_part.eigvecs);
		%sgsScores
		%areafdLM_PCA_part.harmscr( pos_SGS:end, : )
	end
	fprintf( 'size part ....' );
	size( areafdLM_PCA_part.harmscr, 1 )
	figure;
        for iSub = 1:4
                ix = floor( (iSub-1)/2 ) + 1;
                iy = mod( iSub-1, 2 ) + 1;
                subplot( 2, 2, iSub ),
                if ix == iy
                        for icase = 1:size(areafdLM_PCA_part.harmscr, 1)
			        curveComponent_part = areamatLM_PCA_part( :, ix ) * areafdLM_PCA_part.harmscr( icase, ix ) + mean( areamatLM_part(:, 1:pos_SGS-1), 2 );
				%if icase < pos_SGS
                                	plot( agefine_part, curveComponent_part, 'Color', cmap(round(ages( icase )), :), 'LineWidth', 1 );
				%else
				%	plot( agefine_part, curveComponent_part, 'Color', cmap(ages( icase ), :), 'LineWidth', 1, 'LineStyle', '--' );
				%end
                                hold on
                        end
			for icase = pos_SGS:size(areamatLM_part, 2)
				%score_ix = areamatLM_PCA_part( :, ix )' * (areamatLM_part(:, icase) - mean( areamatLM_part(:, 1:pos_SGS-1), 2 ));
				%score_ix = (areaCoef(ix,icase))' * (areafdLM_PCA_part.Lmat)' * (areafdLM_PCA_part.eigvecs);
				curveComponent_part = areamatLM_PCA_part( :, ix ) * sgsScores(icase-pos_SGS+1, ix) + mean( areamatLM_part( :, 1:pos_SGS-1 ), 2 );
				plot( agefine_part, curveComponent_part, 'Color', cmap( round(ages( icase )), : ), 'LineWidth', 2, 'LineStyle', '-.' );
			end
                        hold off
			axis tight
			xlim( [min(agefine_part) max(agefine_part)] );
			ylim( [min(min(areamatLM_part)) max(max(areamatLM_part))] );
			set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
			xlabel( 'Depth along the centerline', 'FontSize', 20, 'FontWeight', 'Bold' );
			ylabel( 'X-sectional area (mm^2)', 'FontSize', 20, 'FontWeight', 'Bold'  );
			title_name = sprintf( 'PCA Component %d, %3.2f%%', ix, areafdLM_PCA_part.varprop(ix)*100 );
			title( title_name, 'FontSize', 20, 'FontWeight', 'Bold' );
                else
                        for icase = 1:size(areafdLM_PCA_part.harmscr,1)
				%if icase < pos_SGS
			      		plot( areafdLM_PCA_part.harmscr( icase, ix ), areafdLM_PCA_part.harmscr( icase, iy ), 'Color', cmap( round(ages( icase )), : ), 'Marker', 'o', 'MarkerSize', 8, 'LineWidth', 2 );
				%else
				%	plot( areafdLM_PCA_part.harmscr( icase, ix ), areafdLM_PCA_part.harmscr( icase, iy ), 'Color', cmap( ages( icase ), : ), 'Marker', '+', 'MarkerSize', 5 );
				%end
                               hold on
                        end
			for icase = pos_SGS:size(areamatLM_part, 2)
					%score_ix = areamatLM_PCA_part( :, ix )' * (areamatLM_part(:, icase) - mean( areamatLM_part(:, 1:pos_SGS-1), 2 ));
					%score_iy = areamatLM_PCA_part( :, iy )' * (areamatLM_part(:, icase) - mean( areamatLM_part(:, 1:pos_SGS-1), 2 ));
					%score_ix = (areaCoef(ix,icase))' * (areafdLM_PCA_part.Lmat)' * (areafdLM_PCA_part.eigvecs);
					%score_iy = (areaCoef(iy,icase))' * (areafdLM_PCA_part.Lmat)' * (areafdLM_PCA_part.eigvecs);
					plot( sgsScores(icase-pos_SGS+1, ix), sgsScores(icase-pos_SGS+1, iy), 'Color', cmap( round(ages(icase)), : ), 'Marker', '*', 'MarkerSize', 9, 'LineWidth', 2 ); 
			end
                        hold off
			axis tight
			xlabel_name = sprintf( 'PCA %d', ix );
			ylabel_name = sprintf( 'PCA %d', iy );
			xlabel( xlabel_name, 'FontSize', 20, 'FontWeight', 'Bold'  );
			ylabel( ylabel_name, 'FontSize', 20, 'FontWeight', 'Bold'  );
			set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
			title( 'Scores', 'FontSize', 20, 'FontWeight', 'Bold'  );
                end
                hold off
        end
	%{
        set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
        set( gcf, 'PaperUnits', 'inches' );
        set( gcf, 'PaperSize', [45 38] );
        set( gcf, 'PaperPositionMode', 'manual' );
        set( gcf, 'PaperPosition', [0 0 45 38] );
        set( gcf, 'renderer', 'painters' );
	%}
        %filename = sprintf( '%s/PCA/Part_StartFromLandmark%d.png', outputPrefix, partId );
        %saveas( gca, filename );
	%print( gcf, '-dpng', filename );
    
    filename = sprintf( '%s/PCA/Part_StartFromLandmark%d.fig', outputPrefix, partId );
    saveas( gca, filename );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction: plot results of boxplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displayBoxplot( agefine, areamatLM, ages, outputPrefix, yLim, gaussian_sigma, uniform_winSize, cases )
close all
bflag = true;
% default settings
factor = 1.5;
fullout = false;
barcol = 'b';
outliercol = 'r';
%color = [ [0.8784, 0.6902, 1.0]; [1, 0, 1]; [0.5451, 0, 0.8] ];
%prob = [0.75, 0.5, 0.25];
color = 'm';
prob = 0.5;
show = false;
method = 'MBD';
depth = [];

nBins = ceil(max(ages)/10.0) * 10;

filename = sprintf('%s/areamatLM.mat', outputPrefix);
save( filename, 'areamatLM' );
filename = sprintf('%s/agefine.mat', outputPrefix);
save( filename, 'agefine' );

% fbplot of all the curves
size( areamatLM' )
size( agefine )
figure, fbplot( areamatLM, agefine, depth, method, show, prob, color, outliercol, barcol, fullout, factor );
filename = sprintf( '%s/functional_boxplot_allCurves_%s.png', outputPrefix, method );
set( gca, 'xlim', [min(agefine), max(agefine)] );
set( gca, 'ylim', [0, yLim] );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
title( 'Functional boxplot for all curves' );
if bflag 
	%set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	set( gcf, 'PaperUnits', 'inches' );
	set( gcf, 'PaperSize', [11 8.5] );
	set( gcf, 'PaperPositionMode', 'manual' );
	set( gcf, 'PaperPosition', [0 0 11 8.5] );
	set( gcf, 'renderer', 'painters' );
	print( gcf, '-dpng', filename );
else
	saveas( gca, filename );
end
%tmpLinePoint = 1:length(agefine);
%figure, boxplot( areamatLM( tmpLinePoint( logical( mod( tmpLinePoint, 30 ) == 0 ) ), :)' );
figure, boxplot( areamatLM', 'boxstyle', 'filled' );
set( gca, 'ylim', [0, yLim] );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
title( 'Pointwise boxplot for all curves' );
filename = sprintf( '%s/pointwise_boxplot_allCurves_allPoints.png', outputPrefix );
	if bflag 
	%set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	set( gcf, 'PaperUnits', 'inches' );
	set( gcf, 'PaperSize', [11 8.5] );
	set( gcf, 'PaperPositionMode', 'manual' );
	set( gcf, 'PaperPosition', [0 0 11 8.5] );
	set( gcf, 'renderer', 'painters' );
	print( gcf, '-dpng', filename );
	else
	saveas( gca, filename );
	end

%num_frames_per_second = 10;
%filename = sprintf( '%s/functional_pointwise.avi', outputPrefix );
%aviobj1 = avifile( filename, 'fps', num_frames_per_second );

%filename = sprintf( '%s/weighted_pointwise.avi', outputPrefix );
%aviobj2 = avifile( filename, 'fps', num_frames_per_second );

%filename = sprintf( '%s/weighted_functional.avi', outputPrefix );
%aviobj3 = avifile( filename, 'fps', num_frames_per_second );

median_age_wbp = [];
median_case_wbp = {};

% window sliding
estimated_ages = min(ages):3:max(ages)
for estimatedId = 1:length( estimated_ages )
  ageId = estimated_ages( estimatedId );
  pickedCurves = (ages >= ageId - uniform_winSize) .* (ages <= ageId + uniform_winSize);
  
  
  close all
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% figure 1: square window, functional boxplot, top left %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  handle1 = figure, subplot(2,2,1), [depthTmp, outliers, medianCurveId] = fbplot( areamatLM( :, logical( pickedCurves ) ), ...
                  agefine, depth, method, show, prob, color, outliercol, barcol, fullout, factor );
  set(handle1, 'Position', [0 0 1000 800] );
  iCountId = 0;
  originalMedianCurveId = -1;
  for iCurveId = 1:length(pickedCurves)
      if pickedCurves( iCurveId ) > 0
          iCountId = iCountId + 1;
      	  if iCountId == medianCurveId
	      originalMedianCurveId = iCurveId;
              break;
          end
      end
  end
  %{
  if originalMedianCurveId > 0
     titlename = sprintf( 'Functional boxplot: %d months, %d curves included, age of median curve: %d' , ageId, sum( pickedCurves ), ages(originalMedianCurveId) );
  else
     titlename = sprintf( 'Functional boxplot: %d months, %d curves included, no median curve', ageId, sum( pickedCurves ) );
  end

  %}
  title( 'Functional boxplot', 'FontSize', 20, 'FontWeight', 'Bold' );
  set( gca, 'ylim', [0, yLim] );
  set( gca, 'xlim', [min(agefine), max(agefine)] );
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% figure 1: square window, histogram, bottom left %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,3), 
  %hist(ages, max(ages));
  hist(ages, nBins);
  hold on
  %bar( ages(originalMedianCurveId), length(find(ages==ages(originalMedianCurveId))), 'm' );
  plot( [ages(originalMedianCurveId), ages(originalMedianCurveId)], [0, length(find(ages==ages(originalMedianCurveId)))], 'm', 'LineWidth', 2 );
  textString = sprintf( '\\downarrow %d (median)', ages(originalMedianCurveId) );
  text( ages(originalMedianCurveId)-2, 1.35, textString, 'Color', 'm', 'FontSize', 20, 'FontWeight', 'Bold');
  plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
  textString = sprintf( '\\uparrow %d months', ageId );
  text( ageId-2, -1.0, textString, 'Color', 'r', 'FontSize', 20, 'FontWeight', 'Bold');
  rectangle('Position', [max(ageId-uniform_winSize,0.1), 0, min(ageId+uniform_winSize, nBins)-max(ageId-uniform_winSize, 0), 1.1], 'Curvature', [0, 0], 'LineWidth', 2, 'LineStyle', '-.', 'EdgeColor', 'r');
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
  hold off

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% figure 1: square window, point-wise, top right %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % square window, point-wise boxplot
  %figure, boxplot( areamatLM( tmpLinePoint( logical( mod( tmpLinePoint, 30 ) == 0 ) ), logical( (ages >= ageId - sizeWindow) .* (ages <= ageId + sizeWindow) ) )' );
  subplot(2,2,2), boxplot( areamatLM( 1:size(areamatLM, 1), logical( pickedCurves ) )', 'boxstyle', 'filled' );
  %{
  hold on
  for iNum = 1:20:size(areamatLM, 1)
      iNumArray = ones( 1, sum( pickedCurves ) );
      plot( iNumArray * (floor(iNum/20)+1), areamatLM( iNum, logical( pickedCurves ) ), 'g.' );
  end
  hold off
  %}
  set( gca, 'ylim', [0, yLim] );
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
  %titlename = sprintf( 'Pointwise boxplot: %d months, Curves in the window: %d', ageId, sum((ages >= ageId - sizeWindow) .* (ages <= ageId + sizeWindow)) );
  title( 'Pointwise boxplot', 'FontSize', 20, 'FontWeight', 'Bold' );

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% figure 1: square window, histogram, bottom right %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,4), 
  %hist(ages, max(ages));
  hist(ages, nBins);
  hold on
  plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
  textString = sprintf( '\\uparrow %d months', ageId );
  text( ageId-2, -0.35, textString, 'Color', 'r', 'FontSize', 20, 'FontWeight', 'Bold');
  rectangle('Position', [max(ageId-uniform_winSize,0.1), 0, min(ageId+uniform_winSize, nBins)-max(ageId-uniform_winSize, 0), 1.1], 'Curvature', [0, 0], 'LineWidth', 2, 'LineStyle', '-.', 'EdgeColor', 'r');
  hold off
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
  filename = sprintf( '%s/functional_pointwise_boxplot_age%d.png', outputPrefix, ageId );
	if bflag
	%set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	set( gcf, 'PaperUnits', 'inches' );
	set( gcf, 'PaperSize', [11 8.5] );
	set( gcf, 'PaperPositionMode', 'manual' );
	set( gcf, 'PaperPosition', [0 0 11 8.5] );
	set( gcf, 'renderer', 'painters' );
	print( gcf, '-dpng', filename );
  	else
	saveas( gca, filename );
	end

  %aviobj1 = addframe( aviobj1, handle1 );

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% figure 2: gaussian window, weighted functional boxplot, top left %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % weighted gaussian window
  w_j_unNormalized = gaussmf( ages, [ gaussian_sigma, ageId ] );
  w_j = w_j_unNormalized / sum( w_j_unNormalized ) ;
  handle2 = figure, subplot(2,2,1), [depthTmp, outliers, medianCurveId, posPercentile] = wfbplot( areamatLM, agefine, w_j, depth, method, show, prob, color, outliercol, barcol, fullout, factor, 0 );
  set(handle2, 'Position', [0 0 1000 800] );
  titlename = sprintf( 'Weighted functional boxplot' );
  title( titlename, 'FontSize', 20, 'FontWeight', 'Bold' );
  set( gca, 'xlim', [min(agefine), max(agefine)] );
  set( gca, 'ylim', [0, yLim] );
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% figure 2: gaussian window, histogram, bottom left %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,3), 
  %hist(ages, max(ages));
  hist(ages, nBins);
  hold on
  %bar( ages(medianCurveId), length(find(ages==ages(medianCurveId))), 'm' );
  plot( [ages(medianCurveId), ages(medianCurveId)], [0, length(find(ages==ages(medianCurveId)))], 'm', 'LineWidth', 2);
  textString = sprintf( '\\downarrow %d (median)', ages(medianCurveId) );
  text( ages(medianCurveId)-2, 1.35, textString, 'Color', 'm', 'FontSize', 20, 'FontWeight', 'Bold');
  plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
  textString = sprintf( '\\uparrow %d months', ageId );
  text( ageId-2, -1.0, textString, 'Color', 'r', 'FontSize', 20, 'FontWeight', 'Bold');
  [sort_age, sort_Id] = sort(ages);
  plot(sort_age, w_j_unNormalized(sort_Id), 'r-.*', 'LineWidth', 2);
  hold off
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );

  median_age_wbp( estimatedId, 1 ) = ageId;
  median_age_wbp( estimatedId, 2 ) = ages(medianCurveId);
  median_case_wbp{ estimatedId } = cases{medianCurveId};

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% figure 2: gaussian window, wighted point-wise plot, top right %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  median_curve = zeros( size(areamatLM, 1), 1 );
  median_value = zeros( size(areamatLM, 1), 1 );
  quartile_25_curve = zeros( size(areamatLM, 1), 1 );
  quartile_75_curve = zeros( size(areamatLM, 1), 1 );
  quartile_5_curve  = zeros( size(areamatLM, 1), 1 );
  quartile_95_curve = zeros( size(areamatLM, 1), 1 );
  for iI = 1:size( areamatLM, 1 )
        for iJ = 1:size( areamatLM, 2 )
            sum_value = sum( w_j .* abs( areamatLM( iI, iJ ) - reshape( areamatLM( iI, : ), size(w_j) ) ) );
            if iJ == 1 || median_value( iI ) > sum_value
                median_curve( iI ) = areamatLM( iI, iJ );
                median_value( iI ) = sum_value;
            end
        end
        [area_sort, sort_id] = sort( reshape( areamatLM( iI, : ), size( areamatLM, 2 ), 1 ) );
        wj_sort = w_j( sort_id );
	wj_value_sum = 0;
        flag_quartile_25 = 0;
        flag_quartile_75 = 0;
        flag_percentile_5 = 0;

        flag_percentile_95 = 0;
        for iJ = 1:length( wj_sort )
            wj_value_sum = wj_value_sum + wj_sort( iJ );
            if flag_quartile_25 == 0 && wj_value_sum >= 0.25
                quartile_25_curve( iI ) = area_sort( iJ );
                flag_quartile_25 = 1;
            end
            if flag_quartile_75 == 0 && wj_value_sum >= 0.75
                quartile_75_curve( iI ) = area_sort( iJ );
                flag_quartile_75 = 1;
            end
            if flag_percentile_5 == 0 && wj_value_sum >= 0.05
                quartile_5_curve( iI ) = area_sort( iJ );
                flag_percentile_5 = 1;
            end
            if flag_percentile_95 == 0 && wj_value_sum >= 0.95
                quartile_95_curve( iI ) = area_sort( iJ );
                flag_percentile_95 = 1;
            end

            if flag_quartile_25 == 1 && flag_quartile_75 == 1 && flag_percentile_95 == 1 && flag_percentile_5 == 1
                break;
            end
	end
  end
  subplot(2,2,2), 
  plot(agefine, quartile_95_curve, 'r-', 'LineWidth', 2 );
  hold on
  plot(agefine, quartile_75_curve, 'g-.', 'LineWidth', 2 );
  plot(agefine, median_curve, 'm-', 'LineWidth', 2);
  plot(agefine, quartile_25_curve, 'c-.', 'LineWidth', 2 );
  plot(agefine, quartile_5_curve, 'b-', 'LineWidth', 2);
  hold off
  legend( '95%', '75%', 'Median', '25%', '5%' );
  set( gca, 'xlim', [min(agefine), max(agefine)] );
  set( gca, 'ylim', [0, yLim] );
  set( gca, 'FontSize', 24, 'FontWeight', 'Normal' );
  title( 'Weighted percentile', 'FontSize', 24, 'FontWeight', 'Normal' );

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% figure 2: gaussian window, histogram, bottom right %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,4), 
  %hist(ages, max(ages));
  hist(ages, nBins);
  hold on
  plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
  textString = sprintf( '\\uparrow %d months', ageId );
  text( ageId-0.7, -0.5, textString, 'Color', 'r', 'FontSize', 20, 'FontWeight', 'Bold');
  plot(sort_age, w_j_unNormalized(sort_Id), 'r-.*', 'LineWidth', 2);
  hold off
  set( gca, 'FontSize', 24, 'FontWeight', 'Normal' );
  filename = sprintf( '%s/weighted_functional_pointwise_boxplot_age%d.fig', outputPrefix, ageId );
	%{
	if bflag
	%set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	set( gcf, 'PaperUnits', 'inches' );
	set( gcf, 'PaperSize', [11 8.5] );
	set( gcf, 'PaperPositionMode', 'manual' );
	set( gcf, 'PaperPosition', [0 0 11 8.5] );
	set( gcf, 'renderer', 'painters' );
	print( gcf, '-dpng', filename );
 	else
	%}
	saveas( gca, filename );
	%end

  %aviobj2 = addframe( aviobj2, handle2 );

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% figure 3: square window, functional boxplot, top left %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  h=figure, subplot(2,2,1), [depthTmp, outliers, medianCurveId] = fbplot( areamatLM( :, logical( pickedCurves ) ), ...
                  agefine, depth, method, show, prob, color, outliercol, barcol, fullout, factor );
  set(h, 'Position', [0 0 1000 800] );
  iCountId = 0;
  originalMedianCurveId = -1;
  for iCurveId = 1:length(pickedCurves)
      if pickedCurves( iCurveId ) > 0
          iCountId = iCountId + 1;
          if iCountId == medianCurveId
              originalMedianCurveId = iCurveId;
              break;
          end
      end
  end
  title( 'Functional boxplot', 'FontSize', 20, 'FontWeight', 'Bold' );
  set( gca, 'xlim', [min(agefine), max(agefine)] );
  set( gca, 'ylim', [0, yLim] );
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% figure 3: square window, histogram, bottom left %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,3), 
  %hist(ages, max(ages));
  hist(ages, nBins);
  hold on
  %bar( ages(originalMedianCurveId), length(find(ages==ages(originalMedianCurveId))), 'm' );
  plot( [ages(originalMedianCurveId), ages(originalMedianCurveId)], [0, length(find(ages==ages(originalMedianCurveId)))], 'm', 'LineWidth', 2 );
  textString = sprintf( '\\downarrow %d (median)', ages(originalMedianCurveId) );
  text( ages(originalMedianCurveId)-2, 1.35, textString, 'Color', 'm', 'FontSize', 20, 'FontWeight', 'Bold');
  plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
  textString = sprintf( '\\uparrow %d months', ageId );
  text( ageId-2, -1.0, textString, 'Color', 'r', 'FontSize', 20, 'FontWeight', 'Bold');
  rectangle('Position', [max(ageId-uniform_winSize,0.1), 0, min(ageId+uniform_winSize, nBins)-max(ageId-uniform_winSize, 0), 1.1], 'Curvature', [0, 0], 'LineWidth', 2, 'LineStyle', '-.', 'EdgeColor', 'r');
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
  hold off

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% figure 3: gaussian window, weighted functional boxplot, top right %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,2), [depthTmp, outliers, medianCurveId, posPercentile] = wfbplot( areamatLM, agefine, w_j, depth, method, show, prob, color, outliercol, barcol, fullout, factor, 0 );
  titlename = sprintf( 'Weighted functional boxplot' );
  title( titlename, 'FontSize', 20, 'FontWeight', 'Bold' );
  set( gca, 'xlim', [min(agefine), max(agefine)] );
  set( gca, 'ylim', [0, yLim] );
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% figure 3: gaussian window, histogram, bottom right %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,4), 
  %hist(ages, max(ages));
  hist(ages, nBins);
  hold on
  %bar( ages(medianCurveId), length(find(ages==ages(medianCurveId))), 'm' );
  plot( [ages(medianCurveId), ages(medianCurveId)], [0, length(find(ages==ages(medianCurveId)))], 'm', 'LineWidth', 2);
  textString = sprintf( '\\downarrow %d (median)', ages(medianCurveId) );
  text( ages(medianCurveId)-2, 1.35, textString, 'Color', 'm', 'FontSize', 20, 'FontWeight', 'Bold');
  plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
  textString = sprintf( '\\uparrow %d months', ageId );
  text( ageId-2, -1.0, textString, 'Color', 'r', 'FontSize', 20, 'FontWeight', 'Bold');
  plot(sort_age, w_j_unNormalized(sort_Id), 'r-.*', 'LineWidth', 2);
  hold off 
  set( h, 'Position', [100, 100, 1024, 768] );
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
  filename = sprintf( '%s/weighted_unweighted_functional_boxplot_age%d.png', outputPrefix, ageId );
	%if bflag
	%set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	%set( gcf, 'PaperUnits', 'inches' );
	%set( gcf, 'PaperSize', [11 8.5] );
	%set( gcf, 'PaperPositionMode', 'manual' );
	%set( gcf, 'PaperPosition', [0 0 11 8.5] );
	%set( gcf, 'renderer', 'painters' );
	%print( gcf, '-dpng', filename );
  	%else
	saveas( h, filename );
	%end

   %aviobj3 = addframe( aviobj3, h );

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % store boxplot's data
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %dataBoxplotMin( :, estimatedId ) = dataBoxplot(1, :);
   %dataBoxplotMax( :, estimatedId ) = dataBoxplot(5, :);
   %dataBoxplotInf( :, estimatedId ) = dataBoxplot(2, :);
   %dataBoxplotSup( :, estimatedId ) = dataBoxplot(4, :);
   %dataBoxplotMedian( :, estimatedId ) = dataBoxplot(3, :);
   %age_median( 1, estimatedId ) = ages( medianCurveId );
end
%aviobj1 = close( aviobj1 );
%aviobj2 = close( aviobj2 );
%aviobj3 = close( aviobj3 );

filename = sprintf( '%s/median_wbp.mat', outputPrefix );
save(filename, 'median_age_wbp', 'median_case_wbp');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prediction of weighted functional boxplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
new_ages = [5; 210; 50; 100; 120; 180; 190; 200];
nDegree = 1;
new_dataBoxplot = curve_prediction( [dataBoxplotMin; dataBoxplotInf; dataBoxplotMedian; dataBoxplotSup; dataBoxplotMax], age_median, new_ages, nDegree, outputPrefix, 0 );
new_dataBoxplotMin = new_dataBoxplot( 1:length(agefine), : );
new_dataBoxplotInf = new_dataBoxplot( 1+length(agefine):2*length(agefine), : );
new_dataBoxplotMedian = new_dataBoxplot( 1+2*length(agefine):3*length(agefine), : );
new_dataBoxplotSup = new_dataBoxplot( 1+3*length(agefine):4*length(agefine), : );
new_dataBoxplotMax = new_dataBoxplot( 1+4*length(agefine):end, : );
new_dataBoxplotMin = min( new_dataBoxplotMin, new_dataBoxplotInf );
new_dataBoxplotMax = max( new_dataBoxplotMax, new_dataBoxplotSup );
new_dataBoxplotMedian = min( new_dataBoxplotMedian, new_dataBoxplotSup );
new_dataBoxplotMedian = max( new_dataBoxplotMedian, new_dataBoxplotInf );
%new_dataBoxplotMin = curve_prediction( dataBoxplotMin, age_median, new_ages, nDegree, outputPrefix, 1 );
%new_dataBoxplotMax = curve_prediction( dataBoxplotMax, age_median, new_ages, nDegree, outputPrefix, 2 );
%new_dataBoxplotInf = curve_prediction( dataBoxplotInf, age_median, new_ages, nDegree, outputPrefix, 3 );
%new_dataBoxplotSup = curve_prediction( dataBoxplotSup, age_median, new_ages, nDegree, outputPrefix, 4 );
%new_dataBoxplotMedian = curve_prediction( dataBoxplotMedian, age_median, new_ages, nDegree, outputPrefix, 5 );

for iTmp = 1:size( new_dataBoxplotMedian, 2 )
    figure, hold on
    [xinv,xindex]=sort(agefine,'descend');
    xx=[agefine; xinv];
    supinv=new_dataBoxplotSup(xindex, iTmp);
    yy=[new_dataBoxplotInf(:, iTmp)', supinv'];
    fill(xx,yy,'m');
    barval=(agefine(1)+agefine(end))/2;
    loc=find(sort([agefine;barval])==barval);
    bar=loc(1);
    line([agefine(bar) agefine(bar)],[new_dataBoxplotMax(bar, iTmp) new_dataBoxplotSup(bar, iTmp)], 'Color', 'b', 'LineWidth', 2);
    line([agefine(bar) agefine(bar)],[new_dataBoxplotMin(bar, iTmp) new_dataBoxplotInf(bar, iTmp)], 'Color', 'b', 'LineWidth', 2);

    plot( agefine, new_dataBoxplotMedian( :, iTmp ), 'k-', 'LineWidth', 2 );
    plot( agefine, new_dataBoxplotMin( :, iTmp ), 'b-', 'LineWidth', 2 );
    plot( agefine, new_dataBoxplotMax( :, iTmp ), 'b-', 'LineWidth', 2 );
    plot( agefine, new_dataBoxplotInf( :, iTmp ), 'b-', 'LineWidth', 2 );
    plot( agefine, new_dataBoxplotSup( :, iTmp ), 'b-', 'LineWidth', 2 );

    set( gca, 'xlim', [min(agefine), max(agefine)] );
    set( gca, 'ylim', [0, yLim] );
    filename = sprintf( '%s/prediction_age%02d.png', outputPrefix, new_ages(iTmp) );
    saveas( gca, filename );
end

for iI = 1:length(new_ages)
	ageId = new_ages(iI);
	w_j_unNormalized = gaussmf( ages, [ gaussian_sigma, ageId ] );
	w_j = w_j_unNormalized / sum( w_j_unNormalized ) ;

	h=figure('Position', [100, 100, 1024, 786] );
	subplot(2,2,1), [depthTmp, outliers, medianCurveId, posPercentile] = wfbplot( areamatLM, agefine, w_j, depth, method, show, prob, color, outliercol, barcol, fullout, factor, 0 );
	titlename = sprintf( 'Weighted functional boxplot' );
	title( titlename, 'FontSize', 20, 'FontWeight', 'Bold' );
	set( gca, 'xlim', [min(agefine), max(agefine)] );
	set( gca, 'ylim', [0, yLim] );

	subplot(2,2,3), hist(ages, nBins);
	hold on
	%bar( ages(medianCurveId), length(find(ages==ages(medianCurveId))), 'm' );
	plot( [ages(medianCurveId), ages(medianCurveId)], [0, length(find(ages==ages(medianCurveId)))], 'm', 'LineWidth', 2);
	textString = sprintf( '\\downarrow %d (median)', ages(medianCurveId) );
	text( ages(medianCurveId)-2, 1.25, textString, 'Color', 'm', 'FontSize', 15, 'FontWeight', 'Bold');
	plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
	textString = sprintf( '\\uparrow %d (months)', ageId );
	text( ageId-2, -0.35, textString, 'Color', 'r', 'FontSize', 15, 'FontWeight', 'Bold');
	plot(sort_age, w_j_unNormalized(sort_Id), 'r-.*', 'LineWidth', 2);
	hold off 
	set( gca, 'ylim', [0, 2] );

	subplot(2,2,2), hold on
	[xinv,xindex]=sort(agefine,'descend');
	xx=[agefine; xinv];
	supinv=new_dataBoxplotSup(xindex, iI);
	yy=[new_dataBoxplotInf(:, iI)', supinv'];
	fill(xx,yy,'m');
	barval=(agefine(1)+agefine(end))/2;
	loc=find(sort([agefine;barval])==barval);
	bar=loc(1);
	line([agefine(bar) agefine(bar)],[new_dataBoxplotMax(bar, iI) new_dataBoxplotSup(bar, iI)],'Color', 'b', 'LineWidth', 2);
	line([agefine(bar) agefine(bar)],[new_dataBoxplotMin(bar, iI) new_dataBoxplotInf(bar, iI)],'Color', 'b', 'LineWidth', 2);
	plot( agefine, new_dataBoxplotMedian( :, iI ), 'k-', 'LineWidth', 2 );
	plot( agefine, new_dataBoxplotMin( :, iI ), 'b-', 'LineWidth', 2 );
	plot( agefine, new_dataBoxplotMax( :, iI ), 'b-', 'LineWidth', 2 );
	plot( agefine, new_dataBoxplotInf( :, iI ), 'b-', 'LineWidth', 2 );
	plot( agefine, new_dataBoxplotSup( :, iI ), 'b-', 'LineWidth', 2 );
	titlename = sprintf( 'Prediction at age %03d', ageId );
	title( titlename, 'FontSize', 20, 'FontWeight', 'Bold' );
	set( gca, 'xlim', [min(agefine), max(agefine)] );
	set( gca, 'ylim', [0, yLim] );
	hold off

	filename = sprintf( '%s/comparison_predict_age%d.png', outputPrefix, ageId );
	saveas( h, filename );
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction: plot boxplot for SGS subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function displayBoxplotForSGS( agefine, areamatLM, ages, cases, outputPrefix, yLim, gaussian_sigma, uniform_winSize, cmap, pos_SGS, areamatLM_whole, landmarksId, iSample )
close all
bflag = true;
% default settings
factor = 1.5;
fullout = false;
barcol = 'b';
outliercol = 'r';
%color = [ [0.8784, 0.6902, 1.0]; [1, 0, 1]; [0.5451, 0, 0.8] ];
%prob = [0.75, 0.5, 0.25];
color = 'm';
prob = 0.5;
show = true;
method = 'MBD';
depth = [];

nBins = ceil(max(ages)/10.0) * 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% functional boxplot for all curves %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fbplot of all the curves
figure;
for icase = pos_SGS : size( areamatLM, 2 )
	plot( agefine, areamatLM(:, icase), 'Color', cmap(round(ages(icase)), :), 'LineWidth', 2, 'LineStyle', '-.' );
	hold on
end
fbplot( areamatLM(:, 1:pos_SGS-1), agefine, depth, method, show, prob, color, outliercol, barcol, fullout, factor );
hold on
for icase = pos_SGS : size( areamatLM, 2 )
	plot( agefine, areamatLM(:, icase), 'Color', cmap(round(ages(icase)), :), 'LineWidth', 2, 'LineStyle', '-.' );
end
if pos_SGS <= size( areamatLM, 2 )
	legend( cases( pos_SGS:end ), 'FontSize', 20, 'FontWeight', 'Bold' );
end
hold off
filename = sprintf( '%s/functional_boxplot_allCurvesWithSGS_%s.png', outputPrefix, method );
set( gca, 'xlim', [min(agefine), max(agefine)] );
set( gca, 'ylim', [0, yLim] );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
title( 'Functional boxplot for all curves' );
if bflag 
	%set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	set( gcf, 'PaperUnits', 'inches' );
	set( gcf, 'PaperSize', [11 8.5] );
	set( gcf, 'PaperPositionMode', 'manual' );
	set( gcf, 'PaperPosition', [0 0 11 8.5] );
	set( gcf, 'renderer', 'painters' );
	print( gcf, '-dpng', filename );
else
	saveas( gca, filename );
end
%tmpLinePoint = 1:length(agefine);
%figure, boxplot( areamatLM( tmpLinePoint( logical( mod( tmpLinePoint, 30 ) == 0 ) ), :)' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% point-wise boxplot for all curves %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
for icase = pos_SGS : size( areamatLM, 2 )
	plot( 1:size(areamatLM, 1), areamatLM(:, icase), 'Color', cmap(round(ages(icase)), :), 'LineWidth', 2, 'LineStyle', '-.' );
	hold on
end
boxplot( areamatLM(:, 1:pos_SGS-1)' );
hold on
for icase = pos_SGS : size( areamatLM, 2 )
	plot( 1:size(areamatLM, 1), areamatLM(:, icase), 'Color', cmap(round(ages(icase)), :), 'LineWidth', 2, 'LineStyle', '-.' );
end
if pos_SGS <= size( areamatLM, 2 )
	legend( cases( pos_SGS:end ), 'FontSize', 20, 'FontWeight', 'Bold' );
end
hold off
set( gca, 'ylim', [0, yLim] );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
title( 'Pointwise boxplot for all curves' );
filename = sprintf( '%s/pointwise_boxplot_allCurves_allPointsWithSGS.png', outputPrefix );
	if bflag 
	%set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	set( gcf, 'PaperUnits', 'inches' );

	set( gcf, 'PaperSize', [11 8.5] );
	set( gcf, 'PaperPositionMode', 'manual' );
	set( gcf, 'PaperPosition', [0 0 11 8.5] );
	set( gcf, 'renderer', 'painters' );
	print( gcf, '-dpng', filename );
	else
	saveas( gca, filename );
	end

%if pos_SGS <= size( areamatLM, 2 )
%	below_portion = zeros( size( areamatLM, 2 )-pos_SGS+1, 5 );
%	max_block = zeros( size( areamatLM, 2 )-pos_SGS+1, 6 );
%end

below_portion = zeros( size( areamatLM, 2 ), 5 );
max_block = zeros( size( areamatLM, 2 ), 7 );

% 1: percent < 5%
% 2: percent reduction in area w.r.t. 5%
% 3: percent reduction in area w.r.t. median
% 4: percent reduction in area w.r.t. 95%
% 5: percent reduction in area w.r.t. itself, b/a 
% 6: percent reduction in area w.r.t. itself's median
% 7: percent reduction in area w.r.t. itself, b/epiglottis
% 8: percent reduction in area w.r.t. itself, b/c
% 9-10: store the starting and ending points for stenosis
% 11: position of stenosis
% 12: percent reduction in area w.r.t. mean(5%)
% 13: percent reduction in area w.r.t. mean(median)
% 14: percent reduction in area w.r.t. mean(95%)
% 
diagnosis = zeros( size( areamatLM, 2 ), 14 );
posInAtlas = zeros( size( areamatLM, 2 ), 2 );
% window sliding
start_id_tmp = 1;
for icase = start_id_tmp : size( areamatLM, 2 )
  ageId = ages( icase );
  pickedCurves = (ages(1:pos_SGS-1) >= ageId - uniform_winSize) .* (ages(1:pos_SGS-1) <= ageId + uniform_winSize);
 
  close all;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%  figure 1, square window, functional boxplot, left %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % square window, functional boxplot
  figure('Position', [100, 100, 1024, 786] ), subplot(2,2,1), 
  plot( agefine, areamatLM(:, icase ), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap(round(ageId), :) );
  hold on
  [depthTmp, outliers, medianCurveId] = fbplot( areamatLM( :, logical( pickedCurves ) ), ...
                  agefine, depth, method, show, prob, color, outliercol, barcol, fullout, factor );
  hold on
  plot( agefine, areamatLM(:, icase ), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap(round(ageId), :) );
  hold off
  legend( cases{icase});
  iCountId = 0;
  originalMedianCurveId = -1;
  for iCurveId = 1:length(pickedCurves)
      if pickedCurves( iCurveId ) > 0
          iCountId = iCountId + 1;
      	  if iCountId == medianCurveId
	      originalMedianCurveId = iCurveId;
              break;
          end
      end
  end
  title( 'Functional boxplot', 'FontSize', 20, 'FontWeight', 'Bold' );
  set( gca, 'ylim', [0, yLim] );
  set( gca, 'xlim', [min(agefine), max(agefine)] );
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
   
  subplot(2,2,3), 
  hist(ages(1:pos_SGS-1), nBins);
  hold on
  %bar( ages(originalMedianCurveId), length(find(ages==ages(originalMedianCurveId))), 'm' );
  plot( [ages(originalMedianCurveId), ages(originalMedianCurveId)], [0, length(find(ages(1:pos_SGS-1)==ages(originalMedianCurveId)))], 'm', 'LineWidth', 2 );
  textString = sprintf( '\\downarrow %d (median)', ages(originalMedianCurveId) );
  text( ages(originalMedianCurveId)-2, 1.35, textString, 'Color', 'm', 'FontSize', 20, 'FontWeight', 'Bold');
  plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
  textString = sprintf( '\\uparrow %d months', ageId );
  text( ageId-2, -1.0, textString, 'Color', 'r', 'FontSize', 20, 'FontWeight', 'Bold');
  rectangle('Position', [max(ageId-uniform_winSize,0.1), 0, min(ageId+uniform_winSize, nBins)-max(ageId-uniform_winSize, 0), 1.1], 'Curvature', [0, 0], 'LineWidth', 2, 'LineStyle', '-.', 'EdgeColor', 'r');
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
  hold off

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%  figure 1, square window, point-wise boxplot, right %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % square window, point-wise boxplot
  %figure, boxplot( areamatLM( tmpLinePoint( logical( mod( tmpLinePoint, 30 ) == 0 ) ), logical( (ages >= ageId - sizeWindow) .* (ages <= ageId + sizeWindow) ) )' );
  subplot(2,2,2),
  plot( 1:size(areamatLM, 1), areamatLM( :, icase ), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap( round(ageId), : ) ); 
  hold on
  boxplot( areamatLM( 1:size(areamatLM, 1), logical( pickedCurves ) )', 'boxstyle', 'filled' );
  hold on
  plot( 1:size(areamatLM, 1), areamatLM( :, icase ), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap( round(ageId), : ) ); 
  hold off
  legend( cases{icase} );
  
  %hold on
  %for iNum = 1:20:size(areamatLM, 1)
  %    iNumArray = ones( 1, sum( pickedCurves ) );
  %    plot( iNumArray * (floor(iNum/20)+1), areamatLM( iNum, logical( pickedCurves ) ), 'g.' );
  %
  %end
  %hold off
  
  set( gca, 'ylim', [0, yLim] );
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
  %titlename = sprintf( 'Pointwise boxplot: %d months, Curves in the window: %d', ageId, sum((ages >= ageId - sizeWindow) .* (ages <= ageId + sizeWindow)) );
  title( 'Pointwise boxplot', 'FontSize', 20, 'FontWeight', 'Bold' );
  
  subplot(2,2,4), 
  hist(ages(1:pos_SGS-1), nBins);
  hold on
  plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
  textString = sprintf( '\\uparrow %d months', ageId );
  text( ageId-2, -1.0, textString, 'Color', 'r', 'FontSize', 20, 'FontWeight', 'Bold');
  rectangle('Position', [max(ageId-uniform_winSize,0.1), 0, min(ageId+uniform_winSize, nBins)-max(ageId-uniform_winSize, 0), 1.1], 'Curvature', [0, 0], 'LineWidth', 2, 'LineStyle', '-.', 'EdgeColor', 'r');
  hold off
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
  filename = sprintf( '%s/functional_pointwise_boxplot_%s_age%d.png', outputPrefix, cases{icase}, ageId );
	%if bflag
	%set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	%set( gcf, 'PaperUnits', 'inches' );
	%set( gcf, 'PaperSize', [11 8.5] );
	%set( gcf, 'PaperPositionMode', 'manual' );
	%set( gcf, 'PaperPosition', [0 0 11 8.5] );
	%set( gcf, 'renderer', 'painters' );
	%print( gcf, '-dpng', filename );
  	%else
	saveas( gca, filename );
	%end
  filename = sprintf( '%s/functional_pointwise_boxplot_%s_age%d.fig', outputPrefix, cases{icase}, ageId );
  saveas( gca, filename );

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%  figure 2, gaussian window, weighted functional boxplot, left %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % weighted gaussian window
  w_j_unNormalized = gaussmf( ages(1:pos_SGS-1), [ gaussian_sigma, ageId ] );
  w_j = w_j_unNormalized / sum( w_j_unNormalized ) ;
  figure('Position', [100, 100, 1024, 786]), subplot(2,2,1), 
  plot( agefine, areamatLM( :, icase ), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap( round(ageId), : ) );
  hold on
  [depthTmp, outliers, medianCurveId, posPercentile, centerId] = wfbplot( areamatLM(:, 1:pos_SGS-1), agefine, w_j, depth, method, show, prob, color, outliercol, barcol, fullout, factor, areamatLM( :, icase ) );
  posInAtlas( icase, 1 ) = posPercentile;
  posInAtlas( icase, 2 ) = ageId;
  textString = sprintf( 'Percentile: %f, Age: %d', posPercentile, round(ageId) );
  text( mean(agefine), yLim*0.8, textString );
  hold on
  plot( agefine, areamatLM( :, icase ), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap( round(ageId), : ) );
  hold off
  legend( cases{icase} );
  titlename = sprintf( 'Weighted functional boxplot' );
  title( titlename, 'FontSize', 20, 'FontWeight', 'Bold' );
  set( gca, 'xlim', [min(agefine), max(agefine)] );
  set( gca, 'ylim', [0, yLim] );
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );

  % save the ids in the 50% central region
  filenameTmp = sprintf( '%s/idsInCentralRegion_%s.txt', outputPrefix, cases{icase} );
  fid = fopen( filenameTmp, 'w' );
  fprintf(fid, '%d\n', length(centerId) );
  for iCenterId = 1:length(centerId)
  	fprintf(fid, '%d\n', centerId(iCenterId) );
  end
  fclose(fid);

  subplot(2,2,3), 
  hist(ages(1:pos_SGS-1), nBins);
  hold on
  %bar( ages(medianCurveId), length(find(ages==ages(medianCurveId))), 'm' );
  plot( [ages(medianCurveId), ages(medianCurveId)], [0, length(find(ages(1:pos_SGS-1)==ages(medianCurveId)))], 'm', 'LineWidth', 2);
  textString = sprintf( '\\downarrow %d (median)', ages(medianCurveId) );
  text( ages(medianCurveId)-2, 1.35, textString, 'Color', 'm', 'FontSize', 20, 'FontWeight', 'Bold');
  plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
  textString = sprintf( '\\uparrow %d months', ageId );
  text( ageId-2, -1.0, textString, 'Color', 'r', 'FontSize', 20, 'FontWeight', 'Bold');
  [sort_age, sort_Id] = sort(ages(1:pos_SGS-1));
  plot(sort_age, w_j_unNormalized(sort_Id), 'r-.*', 'LineWidth', 2);
  hold off
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%  figure 2, gaussian window, wighted point-wise boxplot, right %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  median_curve = zeros( size(areamatLM, 1), 1 );
  median_value = zeros( size(areamatLM, 1), 1 );
  quartile_25_curve = zeros( size(areamatLM, 1), 1 );
  quartile_75_curve = zeros( size(areamatLM, 1), 1 );
  quartile_5_curve  = zeros( size(areamatLM, 1), 1 );
  quartile_95_curve = zeros( size(areamatLM, 1), 1 );
  for iI = 1:size( areamatLM, 1 )
        for iJ = 1:size( areamatLM(:, 1:pos_SGS-1), 2 )
            sum_value = sum( w_j .* abs( areamatLM( iI, iJ ) - reshape( areamatLM( iI, 1:pos_SGS-1 ), size(w_j) ) ) );
            if iJ == 1 || median_value( iI ) > sum_value
                median_curve( iI ) = areamatLM( iI, iJ );
                median_value( iI ) = sum_value;
            end
        end
        [area_sort, sort_id] = sort( reshape( areamatLM( iI, 1:pos_SGS-1 ), size( areamatLM(:, 1:pos_SGS-1), 2 ), 1 ) );
        wj_sort = w_j( sort_id );
	wj_value_sum = 0;
        flag_quartile_25 = 0;
        flag_quartile_75 = 0;
        flag_percentile_5 = 0;

        flag_percentile_95 = 0;
        for iJ = 1:length( wj_sort )
            wj_value_sum = wj_value_sum + wj_sort( iJ );
            if flag_quartile_25 == 0 && wj_value_sum >= 0.25
                quartile_25_curve( iI ) = area_sort( iJ );
                flag_quartile_25 = 1;
            end
            if flag_quartile_75 == 0 && wj_value_sum >= 0.75
                quartile_75_curve( iI ) = area_sort( iJ );
                flag_quartile_75 = 1;
            end
            if flag_percentile_5 == 0 && wj_value_sum >= 0.05
                quartile_5_curve( iI ) = area_sort( iJ );
                flag_percentile_5 = 1;
            end
            if flag_percentile_95 == 0 && wj_value_sum >= 0.95
                quartile_95_curve( iI ) = area_sort( iJ );
                flag_percentile_95 = 1;
            end

            if flag_quartile_25 == 1 && flag_quartile_75 == 1 && flag_percentile_95 == 1 && flag_percentile_5 == 1
                break;
            end
	end
  end

  below_portion( icase, : ) = 0;
  max_block( icase, : ) = 0;
  max_block_index = 0;
  %for iI = 1:length( median_curve )
  %	if iI == 1 || areamatLM( iI, icase ) < max_block( icase, 1 )
  %		max_block( icase, 1 ) = areamatLM( iI, icase );
  %		max_block_index = iI;
  %	end	
  %end

  for iI = 2:length( median_curve )-1
        if areamatLM( iI-1, icase ) >= areamatLM( iI, icase ) && areamatLM( iI+1, icase ) >= areamatLM( iI, icase )
		if max_block_index == 0 || max_block( icase, 1 ) > areamatLM( iI, icase )
			max_block_index = iI;
			max_block( icase, 1 ) = areamatLM( iI, icase );
		end
	end
  end

  if max_block_index == 0
	for iI = 1:length( median_curve )
  		if iI == 1 || areamatLM( iI, icase ) < max_block( icase, 1 )
  			max_block( icase, 1 ) = areamatLM( iI, icase );
  	             	max_block_index = iI;
  		end
  	end
  end

  max_block( icase, 2 ) = agefine( max_block_index );
  max_block( icase, 3 ) = (quartile_5_curve(max_block_index) - max_block(icase, 1 )) / quartile_5_curve(max_block_index);
  max_block( icase, 4 ) = (quartile_25_curve(max_block_index) - max_block(icase, 1)) / quartile_25_curve(max_block_index);
  max_block( icase, 5 ) = (median_curve(max_block_index) - max_block(icase, 1)) / median_curve(max_block_index);
  max_block( icase, 6 ) = (quartile_75_curve(max_block_index) - max_block(icase, 1)) / quartile_75_curve(max_block_index);
  max_block( icase, 7 ) = (quartile_95_curve(max_block_index) - max_block(icase, 1)) / quartile_95_curve(max_block_index);

  iStartPoint = max_block_index;
  while iStartPoint >= 1 && areamatLM( iStartPoint, icase ) < quartile_5_curve( iStartPoint ) 
	iStartPoint = iStartPoint - 1;
	if iStartPoint >=1 && areamatLM( iStartPoint, icase ) < areamatLM( iStartPoint+1, icase )
		break;
	end
  end
  iStartPoint = iStartPoint + 1;
  iEndPoint = max_block_index;
  while iEndPoint <= length( median_curve ) && areamatLM( iEndPoint, icase ) < quartile_5_curve( iEndPoint )
	iEndPoint = iEndPoint + 1;
	if iEndPoint <= length( median_curve ) && areamatLM( iEndPoint, icase ) < areamatLM( iEndPoint-1, icase )
		break;
	end
  end
  iEndPoint = iEndPoint - 1;
  if iEndPoint >= iStartPoint
	below_portion( icase, 1 ) = agefine( iEndPoint ) - agefine(iStartPoint);
  end

  iStartPoint = max_block_index;
  while iStartPoint >= 1 && areamatLM( iStartPoint, icase ) < quartile_25_curve( iStartPoint ) 
	iStartPoint = iStartPoint - 1;
	if iStartPoint >=1 && areamatLM( iStartPoint, icase ) < areamatLM( iStartPoint+1, icase )
                break;
        end
  end
  iStartPoint = iStartPoint + 1;
  iEndPoint = max_block_index;
  while iEndPoint <= length( median_curve ) && areamatLM( iEndPoint, icase ) < quartile_25_curve( iEndPoint )
	iEndPoint = iEndPoint + 1;
	if iEndPoint <= length( median_curve ) && areamatLM( iEndPoint, icase ) < areamatLM( iEndPoint-1, icase )
                break;
        end
  end
  iEndPoint = iEndPoint - 1;
  if iEndPoint >= iStartPoint
	below_portion( icase, 2 ) = agefine( iEndPoint ) - agefine(iStartPoint);
  end

  iStartPoint = max_block_index;
  while iStartPoint >= 1 && areamatLM( iStartPoint, icase ) < median_curve( iStartPoint ) 
	iStartPoint = iStartPoint - 1;
	if iStartPoint >=1 && areamatLM( iStartPoint, icase ) < areamatLM( iStartPoint+1, icase )
                break;
        end
  end
  iStartPoint = iStartPoint + 1;
  iEndPoint = max_block_index;
  while iEndPoint <= length( median_curve ) && areamatLM( iEndPoint, icase ) < median_curve( iEndPoint )
	iEndPoint = iEndPoint + 1;
	if iEndPoint <= length( median_curve ) && areamatLM( iEndPoint, icase ) < areamatLM( iEndPoint-1, icase )
                break;
        end
  end
  iEndPoint = iEndPoint - 1;
  if iEndPoint >= iStartPoint
	below_portion( icase, 3 ) = agefine( iEndPoint ) - agefine(iStartPoint);
  end

  iStartPoint = max_block_index;
  while iStartPoint >= 1 && areamatLM( iStartPoint, icase ) < quartile_75_curve( iStartPoint ) 
	iStartPoint = iStartPoint - 1;
	if iStartPoint >=1 && areamatLM( iStartPoint, icase ) < areamatLM( iStartPoint+1, icase )
                break;
        end
  end
  iStartPoint = iStartPoint + 1;
  iEndPoint = max_block_index;
  while iEndPoint <= length( median_curve ) && areamatLM( iEndPoint, icase ) < quartile_75_curve( iEndPoint )
	iEndPoint = iEndPoint + 1;
	if iEndPoint <= length( median_curve ) && areamatLM( iEndPoint, icase ) < areamatLM( iEndPoint-1, icase )
                break;
        end
  end
  iEndPoint = iEndPoint - 1;
  if iEndPoint >= iStartPoint
	below_portion( icase, 4 ) = agefine( iEndPoint ) - agefine(iStartPoint);
  end

  iStartPoint = max_block_index;
  while iStartPoint >= 1 && areamatLM( iStartPoint, icase ) < quartile_95_curve( iStartPoint ) 
	iStartPoint = iStartPoint - 1;
	if iStartPoint >=1 && areamatLM( iStartPoint, icase ) < areamatLM( iStartPoint+1, icase )
                break;
        end
  end
  iStartPoint = iStartPoint + 1;
  iEndPoint = max_block_index;
  while iEndPoint <= length( median_curve ) && areamatLM( iEndPoint, icase ) < quartile_95_curve( iEndPoint )
	iEndPoint = iEndPoint + 1;
	if iEndPoint <= length( median_curve ) && areamatLM( iEndPoint, icase ) < areamatLM( iEndPoint-1, icase )
                break;
        end
  end
  iEndPoint = iEndPoint - 1;
  if iEndPoint >= iStartPoint
	below_portion( icase, 5 ) = agefine( iEndPoint ) - agefine(iStartPoint);
  end

  below_portion( icase, : ) = below_portion( icase, : ) / ( agefine(end) - agefine(1) );

  diagnosis( icase, 1 ) = below_portion( icase, 1 );
  diagnosis( icase, 2 ) = max_block( icase, 3 );
  diagnosis( icase, 3 ) = max_block( icase, 5 );
  diagnosis( icase, 4 ) = max_block( icase, 7 );
  diagnosis( icase, 6 ) = ( median(areamatLM(:, icase)) - max_block( icase, 1 ) ) / median(areamatLM( :, icase ) ); 
  diagnosis( icase, 7 ) = ( areamatLM_whole(landmarksId(3), icase) - max_block( icase, 1 ) ) / areamatLM_whole( landmarksId(3), icase );

  diagnosis( icase, 12 ) = ( mean( quartile_5_curve ) - max_block(icase, 1) ) / mean( quartile_5_curve);
  diagnosis( icase, 13 ) = ( mean( median_curve ) - max_block(icase, 1 ) ) / mean( median_curve );
  diagnosis( icase, 14 ) = ( mean( quartile_95_curve ) - max_block(icase, 1) ) / mean( quartile_95_curve );

  iStartPoint = max_block_index + landmarksId(iSample) - 1;
  while iStartPoint >= landmarksId(4)
	iStartPoint = iStartPoint - 1;
	if iStartPoint >=1 && areamatLM_whole( iStartPoint, icase ) < areamatLM_whole( iStartPoint+1, icase )
                break;
        end
  end
  iStartPoint = iStartPoint + 1;
  flagStart = true;
  if iStartPoint == landmarksId(4)
	if areamatLM_whole( iStartPoint-1, icase ) > areamatLM_whole( iStartPoint, icase )
		flagStart = false;
	end
  end

  iEndPoint = max_block_index + landmarksId(iSample) - 1;
  while iEndPoint <= size(areamatLM_whole, 1)
        iEndPoint = iEndPoint + 1;
        if iEndPoint <= size(areamatLM_whole, 1) && areamatLM_whole( iEndPoint, icase ) < areamatLM_whole( iEndPoint-1, icase )
                break;
        end
  end
  iEndPoint = iEndPoint - 1;  

  %if iStartPoint == max_block_index + landmarksId(iSample) - 1
  % 	baseTmp = areamatLM_whole( iEndPoint, icase );
  %else
	%baseTmp = 0.5 * ( areamatLM_whole( iEndPoint, icase ) + areamatLM_whole( iStartPoint, icase ) );
  	%baseTmp = areamatLM_whole( iStartPoint, icase );
  %end

  %baseTmp = mean( areamatLM_whole(landmarksId(3):landmarksId(4), icase) );
  %if iStartPoint < landmarksId(4)
  %	iStartPointNew = max( iStartPoint, 2 * (max_block_index + landmarksId(iSample) - 1) - iEndPoint );
  %end
  if flagStart == true
	%baseTmp = max( areamatLM_whole(iEndPoint, icase ), areamatLM_whole(iStartPoint, icase ) );
  	baseTmp = 0.5 * ( areamatLM_whole(iEndPoint, icase ) + areamatLM_whole( iStartPoint, icase ) );
  else
	baseTmp = areamatLM_whole( iEndPoint, icase );
  end

  diagnosis( icase, 9 ) = iStartPoint;
  diagnosis( icase, 10 ) = max_block_index + landmarksId(iSample) - 1;
  diagnosis( icase, 11 ) = areamatLM(max_block_index, icase);
  %maxBaseTmp = max( areamatLM( iStartPoint, icase ), areamatLM( iEndPoint, icase ) );
  diagnosis( icase, 8 ) = ( baseTmp - max_block( icase, 1 ) ) / baseTmp; 
  
  %baseTmp_for5 = areamatLM_whole( iStartPoint, icase );
  %if areamatLM_whole( iStartPoint, icase ) < areamatLM_whole( landmarksId(4), icase )
  %	baseTmp_for5 = areamatLM_whole( landmarksId(4), icase ); 
  %end
  %baseTmp_for5 = max( areamatLM_whole( landmarksId(4):end, icase ) );
  %baseTmp_for5 = mean( areamatLM_whole( landmarksId(3):landmarksId(4), icase ) );
  baseTmp_for5 = 0.5 * ( areamatLM_whole( landmarksId(3), icase ) + areamatLM_whole( landmarksId(4), icase ) );
  diagnosis( icase, 5 ) = ( baseTmp_for5 - max_block( icase, 1 ) ) / baseTmp_for5;
  
  baseTmp_for6 = median( areamatLM_whole( landmarksId(3):end, icase ) );
  diagnosis( icase, 6 ) = ( baseTmp_for6 - max_block( icase, 1 ) ) / baseTmp_for6; 
  %diagnosis( icase, 5 ) = ( areamatLM_whole( iStartPoint, icase ) - max_block(icase, 1) ) / areamatLM_whole( iStartPoint, icase );
  %diagnosis( icase, 8 ) = ( areamatLM_whole( iEndPoint, icase ) - max_block(icase, 1) ) / areamatLM_whole( iEndPoint, icase );  

  %iStartPoint = max_block_index + landmarksId(iSample) - 1;
  %while iStartPoint >= landmarksId(3)
  %      iStartPoint = iStartPoint - 1;
  %      if iStartPoint >=1 && areamatLM_whole( iStartPoint, icase ) < areamatLM_whole( iStartPoint+1, icase )
  %              break;
  %      end
  %end
  %iStartPoint = iStartPoint + 1; 
  %diagnosis( icase, 6 ) = ( median( areamatLM_whole( iStartPoint:end, icase ) ) - max_block(icase, 1 ) ) / median( areamatLM_whole( iStartPoint:end, icase ) );
 
%{
  for iI = max_block_index:-1:1  	
	if areamatLM( iI, icase ) < quartile_5_curve( iI )
		reductionAreaCur = quartile_5_curve(iI) - areamatLM(iI, icase);
		if reductionAreaCur < 
			below_portion( icase, 1 ) = below_portion( icase, 1 ) + 1;
		end	
	end
	if areamatLM( iI, icase ) < quartile_25_curve( iI )
		below_portion( icase, 2 ) = below_portion( icase, 2 ) + 1;
	end
	if areamatLM( iI, icase ) < median_curve( iI ) 
		below_portion( icase, 3 ) = below_portion( icase, 3 ) + 1;
	end
	if areamatLM( iI, icase ) < quartile_75_curve( iI )
		below_portion( icase, 4 ) = below_portion( icase, 4 ) + 1;
	end
	if areamatLM( iI, icase ) < quartile_95_curve( iI )
		below_portion( icase, 5 ) = below_portion( icase, 5 ) + 1;
	end 
  end  
  below_portion( icase, : ) = below_portion( icase, : ) ./ length( median_curve );
%}

  subplot(2,2,2), 
  plot(agefine, quartile_95_curve, 'r-', 'LineWidth', 2 );
  hold on
  plot(agefine, quartile_75_curve, 'g-.', 'LineWidth', 2 );
  plot(agefine, median_curve, 'm-', 'LineWidth', 2);
  plot(agefine, quartile_25_curve, 'c-.', 'LineWidth', 2 );
  plot(agefine, quartile_5_curve, 'b-', 'LineWidth', 2);
  legend( '95%', '75%', 'Median', '25%', '5%' );
  plot( agefine, areamatLM(:, icase), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap(round(ageId), :) );
  hold off
  %{
  hold on
  plot(agefine, median_curve, 'k-', 'LineWidth', 2);
  plot(agefine, quartile_5_curve, 'g-', 'LineWidth', 2);
  plot(agefine, quartile_95_curve, 'g-', 'LineWidth', 2);
  plot(agefine, quartile_25_curve, 'b-', 'LineWidth', 2);
  plot(agefine, quartile_75_curve', 'b-', 'LineWidth', 2);
  hold off
  %}
  legend( '95%', '75%', 'Median', '25%', '5%', cases{icase} );
  set( gca, 'xlim', [min(agefine), max(agefine)] );
  set( gca, 'ylim', [0, yLim] );
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
  title( 'Weighted percentile', 'FontSize', 20, 'FontWeight', 'Bold' );


  subplot(2,2,4), 
  hist(ages(1:pos_SGS-1), nBins);
  hold on
  plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
  textString = sprintf( '\\uparrow %d months', ageId );
  text( ageId-0.7, -1.0, textString, 'Color', 'r', 'FontSize', 20, 'FontWeight', 'Bold');
  plot(sort_age, w_j_unNormalized(sort_Id), 'r-.*', 'LineWidth', 2);
  hold off
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
  %if icase < pos_SGS
  %	filename = sprintf( '%s/weighted_functional_pointwise_boxplot_%s_age%0.2f.png', outputPrefix, cases{icase}, ageId );
  %else
  %	filename = sprintf( '%s/weighted_functional_pointwise_boxplot_%s_age%0.2f.png', outputPrefix, cases{icase}, ageId );
  %end	
        %if bflag
	%set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	%set( gcf, 'PaperUnits', 'inches' );
	%set( gcf, 'PaperSize', [11 8.5] );
	%set( gcf, 'PaperPositionMode', 'manual' );
	%set( gcf, 'PaperPosition', [0 0 11 8.5] );
	%set( gcf, 'renderer', 'painters' );
	%print( gcf, '-dpng', filename );
 	%else
 	%
	%saveas( gca, filename );
	%end

  if icase < pos_SGS
        filename = sprintf( '%s/weighted_functional_pointwise_boxplot_%s_age%0.2f.fig', outputPrefix, cases{icase}, ageId );
  else
        filename = sprintf( '%s/weighted_functional_pointwise_boxplot_%s_age%0.2f.fig', outputPrefix, cases{icase}, ageId );
  end
  saveas( gca, filename );

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%  figure 3, square window, functional boxplot, left %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  h=figure('Position', [100, 100, 1024, 786] ), subplot(2,2,1), 
  plot( agefine, areamatLM(:, icase), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap(round(ageId), :) );
  hold on
  [depthTmp, outliers, medianCurveId] = fbplot( areamatLM( :, logical( pickedCurves ) ), ...
                  agefine, depth, method, show, prob, color, outliercol, barcol, fullout, factor );
  hold on
  plot( agefine, areamatLM(:, icase), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap(round(ageId), :) );
  hold off
  iCountId = 0;
  originalMedianCurveId = -1;
  for iCurveId = 1:length(pickedCurves)
      if pickedCurves( iCurveId ) > 0
          iCountId = iCountId + 1;
          if iCountId == medianCurveId
              originalMedianCurveId = iCurveId;
              break;
          end
      end
  end
  legend( cases{icase} );
  title( 'Functional boxplot', 'FontSize', 20, 'FontWeight', 'Bold' );
  set( gca, 'xlim', [min(agefine), max(agefine)] );
  set( gca, 'ylim', [0, yLim] );
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );

  subplot(2,2,3), 
  %hist(ages, max(ages));
  hist(ages(1:pos_SGS-1), nBins);
  hold on
  %bar( ages(originalMedianCurveId), length(find(ages==ages(originalMedianCurveId))), 'm' );
  plot( [ages(originalMedianCurveId), ages(originalMedianCurveId)], [0, length(find(ages(1:pos_SGS-1)==ages(originalMedianCurveId)))], 'm', 'LineWidth', 2 );
  textString = sprintf( '\\downarrow %d (median)', ages(originalMedianCurveId) );
  text( ages(originalMedianCurveId)-2, 1.35, textString, 'Color', 'm', 'FontSize', 20, 'FontWeight', 'Bold');
  plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
  textString = sprintf( '\\uparrow %d months', ageId );
  text( ageId-2, -1.0, textString, 'Color', 'r', 'FontSize', 20, 'FontWeight', 'Bold');
  rectangle('Position', [max(ageId-uniform_winSize,0.1), 0, min(ageId+uniform_winSize, nBins)-max(ageId-uniform_winSize, 0), 1.1], 'Curvature', [0, 0], 'LineWidth', 2, 'LineStyle', '-.', 'EdgeColor', 'r');
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
  hold off

  % compute the distance between the median curve and the estimated curve
  title_ssd = sprintf( 'SSD: %f', norm( areamatLM( :, originalMedianCurveId ) - areamatLM( :, icase ), 2 ) );	
  title( title_ssd );

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%  figure 3, gaussian window, weighted functional boxplot, right %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,2), 
  plot( agefine, areamatLM(:, icase), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap(round(ageId), :) );
  hold on
  [depthTmp, outliers, medianCurveId, posPercentile] = wfbplot( areamatLM(:, 1:pos_SGS-1), agefine, w_j, depth, method, show, prob, color, outliercol, barcol, fullout, factor, areamatLM(:, icase) );
  textString = sprintf( 'Percentile: %f, Age: %d', posPercentile, round(ageId) );
  text( mean(agefine), mean(areamatLM(:, icase)), textString );
  hold on
  plot( agefine, areamatLM(:, icase), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap(round(ageId), :) );
  hold off
  legend( cases{icase} );
  titlename = sprintf( 'Weighted functional boxplot' );
  title( titlename, 'FontSize', 20, 'FontWeight', 'Bold' );
  set( gca, 'xlim', [min(agefine), max(agefine)] );
  set( gca, 'ylim', [0, yLim] );
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );

  subplot(2,2,4), 
  %hist(ages, max(ages));
  hist(ages(1:pos_SGS-1), nBins);
  hold on
  %bar( ages(medianCurveId), length(find(ages==ages(medianCurveId))), 'm' );
  plot( [ages(medianCurveId), ages(medianCurveId)], [0, length(find(ages(1:pos_SGS-1)==ages(medianCurveId)))], 'm', 'LineWidth', 2);
  textString = sprintf( '\\downarrow %d (median)', ages(medianCurveId) );
  text( ages(medianCurveId)-2, 1.35, textString, 'Color', 'm', 'FontSize', 20, 'FontWeight', 'Bold');
  plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
  textString = sprintf( '\\uparrow %d months', ageId );
  text( ageId-2, -1.0, textString, 'Color', 'r', 'FontSize', 20, 'FontWeight', 'Bold');
  plot(sort_age, w_j_unNormalized(sort_Id), 'r-.*', 'LineWidth', 2);
  hold off 
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );

  % compute the distance between the median curve and the estimated curve
  title_ssd = sprintf( 'SSD: %f', norm( areamatLM( :, medianCurveId ) - areamatLM( :, icase ), 2 ) );	
  title( title_ssd );

  %set( h, 'Position', [100, 100, 1024, 768] );
  filename = sprintf( '%s/weighted_unweighted_functional_boxplot_%s_age%d.png', outputPrefix, cases{icase}, ageId );
	%if bflag
	%set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	%set( gcf, 'PaperUnits', 'inches' );
	%set( gcf, 'PaperSize', [11 8.5] );
	%set( gcf, 'PaperPositionMode', 'manual' );
	%set( gcf, 'PaperPosition', [0 0 11 8.5] );
	%set( gcf, 'renderer', 'painters' );
	%print( gcf, '-dpng', filename );

  	%else
	saveas( gca, filename );
	%end
  filename = sprintf( '%s/weighted_unweighted_functional_boxplot_%s_age%d.fig', outputPrefix, cases{icase}, ageId );
  saveas( gca, filename );
end

figure, hold on
for icase = 1:size(posInAtlas, 1)
	if icase < pos_SGS
		plot( posInAtlas( icase, 1), posInAtlas( icase, 2), 'Marker', 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', cmap(round(ages(icase)), :), 'MarkerSize', 8 ); 
	elseif icase < pos_SGS+5   % the number of the pre-surgery subjects
		plot( posInAtlas( icase, 1), posInAtlas( icase, 2), 'Marker', 'p', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', cmap(round(ages(icase)), :), 'MarkerSize', 12 );
	else
		plot( posInAtlas( icase, 1), posInAtlas( icase, 2), 'Marker', 's', 'MarkerEdgeColor', 'm', 'MarkerFaceColor', cmap(round(ages(icase)), :), 'MarkerSize', 10 );
	end
end
hold off
filename = sprintf( '%s/posInAtlas.png', outputPrefix );
saveas( gca, filename );

filename = sprintf( '%s/posInAtlas.mat', outputPrefix );
save( filename, 'posInAtlas' );

quartile_x = zeros( 1, 5 );
quartile_x( 1 ) = 0.05; quartile_x( 2 ) = 0.25; quartile_x( 3 ) = 0.5; quartile_x( 4 ) = 0.75; quartile_x( 5 ) = 0.95;
figure;
for icase = size( areamatLM, 2 ):-1:start_id_tmp
	if icase >= pos_SGS
		plot( quartile_x, below_portion( icase, : )*100, 'LineStyle', '--', 'LineWidth', 3, 'Marker', '*', 'MarkerSize', 10, 'Color', cmap( round( ages(icase) ), : ) );
	else
		plot( quartile_x, below_portion( icase, : )*100, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', 10, 'Color', cmap( round( ages( icase ) ), : ) );
	end
	hold on
end
hold off
xlabel( 'Weighted quartile' );
ylabel( 'Percent below quartile (%)' );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
filename = sprintf( '%s/percent_below_quartile.png', outputPrefix );
saveas( gca, filename );


figure;
for icase = 1:pos_SGS-1
	plot( quartile_x, below_portion( icase, : )*100, 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10, 'Color', cmap( round( ages( icase ) ), : ) );
        hold on
end
hold off
xlabel( 'Weighted quartile' );
ylabel( 'Percent below quartile (%)' );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
filename = sprintf( '%s/percent_below_quartile_crl.png', outputPrefix );
saveas( gca, filename );

figure;
for icase = pos_SGS:size( areamatLM, 2 )
	plot( quartile_x, below_portion( icase, : )*100, 'LineStyle', '--', 'LineWidth', 3, 'Marker', '*', 'MarkerSize', 10, 'Color', cmap( round( ages(icase) ), : ) );
        hold on
end
hold off
xlabel( 'Weighted quartile' );
ylabel( 'Percent below quartile (%)' );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
filename = sprintf( '%s/percent_below_quartile_sgs.png', outputPrefix );
saveas( gca, filename );

filename = sprintf( '%s/percent_below_quartile.mat', outputPrefix );
save( filename, 'below_portion' );

figure;
for icase = size( areamatLM, 2 ): -1 : start_id_tmp
	for iTmp = 3:7
		if max_block( icase, iTmp ) < 0 
			max_block( icase, iTmp ) = 0;
		end
	end
	if icase >= pos_SGS
		plot( quartile_x, max_block( icase, 3:7 )*100, 'LineStyle', '--', 'LineWidth', 3, 'Marker', '*', 'MarkerSize', 10, 'Color', cmap( round(ages(icase) ), : ) );
	else
		plot( quartile_x, max_block( icase, 3:7 )*100, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', 10, 'Color', cmap( round( ages( icase ) ), : ) );
	end
	hold on
end
hold off
xlabel( 'Weighted quartile' );
ylabel( 'Percent relative reduction in area (%)' );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
filename = sprintf( '%s/percent_reduciton_area.png', outputPrefix );
saveas( gca, filename );

figure;
for icase = 1:pos_SGS-1
        plot( quartile_x, max_block( icase, 3:7 )*100, 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10, 'Color', cmap( round( ages( icase ) ), : ) );
	hold on
end
hold off
xlabel( 'Weighted quartile' );
ylabel( 'Percent relative reduction in area (%)' );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
filename = sprintf( '%s/percent_reduciton_area_crl.png', outputPrefix );
saveas( gca, filename );

figure;
for icase = pos_SGS : size( areamatLM, 2 )
	plot( quartile_x, max_block( icase, 3:7 )*100, 'LineStyle', '--', 'LineWidth', 3, 'Marker', '*', 'MarkerSize', 10, 'Color', cmap( round(ages(icase) ), : ) );
        hold on
end
hold off
xlabel( 'Weighted quartile' );
ylabel( 'Percent relative reduction in area (%)' );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
filename = sprintf( '%s/percent_reduciton_area_sgs.png', outputPrefix );
saveas( gca, filename );

filename = sprintf( '%s/percent_reduction_area.mat', outputPrefix );
save( filename, 'max_block' );

filename = sprintf( '%s/area.mat', outputPrefix );
save( filename, 'areamatLM_whole' );

filename = sprintf( '%s/diagnosis.mat', outputPrefix );
save( filename, 'diagnosis' );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end ----
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
filename = sprintf( '%s/mean.mat', outputPrefix);
meanAreaLM = mean(areamatLM, 2);
maxAreaLM = max(areamatLM, [], 2 );
minAreaLM = min(areamatLM, [], 2 );
for iLM = 1:size(warpedLM, 2)
    meanXLM( iLM ) = mean( warpedLM( :, iLM ) );
    meanYLM( iLM ) = mean( valueLM( :, iLM ) );
end
save( filename, 'meanAreaLM', 'maxAreaLM', 'minAreaLM', 'meanXLM', 'meanYLM' ); 
%}

%{
% continuous registration after landmark-based registration
areaLMmeanfd = mean(areafdLM);
[areaLMCregfd, areaLMCwarpfd] = registerfd(areaLMmeanfd, areafdLM, WfdParLM);
areaLMCfine = eval_fd(agefine, areaLMCregfd);
areaLMCwarpmat = eval_fd(agefine, areaLMCwarpfd);
%areaLMCwarpmat(1,:) = 0; areaLMCwarpmat(length(agefine), :) = 1;
figure, plot(agefine, areaLMCfine, 'LineWidth', 2);
hold on
plot(agefine, mean(areaLMCfine, 2), 'm--', 'LineWidth', 3);
for icase = 1:numf
     for iLM = 1:size(areaLMCwarpmat, 1)
	 areaLMCwarpmat(iLM, icase) = max( areaLMCwarpmat(iLM, icase), 0 );
	 areaLMCwarpmat(iLM, icase) = min( areaLMCwarpmat(iLM, icase), 1 );
     end 
     warpInvfd = smooth_basis(areaLMCwarpmat(:,icase), agefine, fdPar_new);
     warpedLMC(icase,:) = eval_fd(warpedLM(icase,:), warpInvfd);
     for iLM = 1:size(warpedLMC, 2)
	 warpedLMC(icase, iLM) = max( warpedLMC(icase, iLM), 0 );
	 warpedLMC(icase, iLM) = min( warpedLMC(icase, iLM), 1 );
     end
     valueLMC(icase,:) = eval_fd(warpedLMC(icase,:), areaLMCregfd(icase));
     for iLM = 1:size(warpedLMC, 2)    
     	 plot(warpedLMC(icase, iLM), valueLMC(icase, iLM), char(style(mod(iLM, 5)+1)), 'LineWidth', 10);
     end
end
hold off
xlabel( 'Position' );
ylabel( 'Area' );
legend( cases );
title( 'Continuous Registered curves After Landmark-based Registration' );
filename = sprintf( '%s/continuous_landmark_registration.png', outputPrefix ); 
saveas( gca, filename );

figure, plot( agefine, mean(areaLMCfine, 2), 'm--', 'LineWidth', 3 );
hold on
plot( agefine, max(areaLMCfine, [], 2), 'b-', 'LineWidth', 2 );
plot( agefine, min(areaLMCfine, [], 2), 'b-', 'LineWidth', 2 );
for iLM = 1:size(warpedLMC, 2)
    plot( mean( warpedLMC(:, iLM) ), mean( valueLMC(:, iLM) ), char(style(mod(iLM, 5)+1)), 'LineWidth', 10);
end
hold off
xlabel( 'Position' );
ylabel( 'Area' );
title( 'Mean Curve After Continuous and Landmark-based Registration' );
[valueMin, idMin] = min(min(areaLMCfine, [], 2));
tmpstr = sprintf( '\\leftarrow minimal:(%f, %f)', agefine(idMin), valueMin );
text( agefine(idMin), valueMin, tmpstr, 'HorizontalAlignment', 'left' );
[valueMax, idMax] = max(max(areaLMCfine, [], 2));
tmpstr = sprintf( '\\leftarrow maximal:(%f, %f)', agefine(idMax), valueMax );
text( agefine(idMax), valueMax, tmpstr, 'HorizontalAlignment', 'left' );
filename = sprintf( '%s/continuous_landmark_registration_mean_variation.png', outputPrefix );
saveas( gca, filename );
%}

% index = 1:numf;
% nbasisw = 15;
% norder  =  5;
% basisw  = create_bspline_basis([0,1], nbasisw, norder);
% coef0 = zeros(nbasisw, length(index));
% Wfd0 = fd(coef0, basisw);
% Lfdobj = int2Lfd(2);
% lambda = 1;
% WfdPar = fdPar(Wfd0, Lfdobj, lambda);
% % register landmark-registered area
% areaLMmeanfd = mean(areafdLM);
% [areaLMCregfd, areaLMCwarpfd] = registerfd(areaLMmeanfd, areafdLM, WfdPar);
% areaLMCfine = eval_fd(agefine, areaLMCregfd);
% areaLMCwarpmat = eval_fd(agefine, areaLMCwarpfd);
% areaLMCwarpmat(1,:) = 0; areaLMCwarpmat(length(agefine), :) = 1;
% figure, plot(agefine, areaLMCfine, 'LineWidth', 2);
% hold on
% for icase = 1:numf
%     warpInvfd = smooth_basis(areaLMCwarpmat(:,icase), agefine, fdPar_new);
%     warpedLMC(icase,:) = eval_fd(warpedLM(icase,:), warpInvfd);
%     valueLM(icase,:) = eval_fd(warpedLMC(icase,:), areaLMCregfd(icase));
%     plot(warpedLMC(icase, :), valueLM(icase, :), char(style(icase)));
% end
% for nLM = 1:nLandmarks 
%     linex(:) = LandmarksMean(nLM);
%     plot(linex, liney, 'm-');
% end
% hold off
% 
% 
% % landmarks registration for radius
% nLandmarks = 4;
% Landmarks = zeros(numf, nLandmarks);
% valueLM = zeros(numf, nLandmarks);
% figure, subplot(1,1,1),
% for icase = 1:numf
%     radius = eval_fd(agefine, radiusfd_new(icase));
%     plot(agefine, radius);
%     for nLM = 1:nLandmarks
%         value = ginput(1);
%         Landmarks(icase, nLM) = value(1);
%         valueLM(icase, nLM) = eval_fd(value(1), radiusfd_new(icase));
%     end
%     pause;
% end
% 
% style = {'bo', 'go', 'ro', 'co'};
% radiusfine_new = eval_fd(agefine, radiusfd_new);
% figure, plot(agefine, radiusfine_new, 'LineWidth', 2);
% hold on
% for icase = 1:numf
%     plot(Landmarks(icase, :), valueLM(icase, :), char(style(icase)));
% end
% hold off
% 
% wbasisLM = create_bspline_basis([0,1], max(nLandmarks+3-2,4), 3);
% WfdLM    = fd(zeros(max(nLandmarks+3-2,4),1),wbasisLM);
% WfdParLM = fdPar(WfdLM,1,1e-12);
% 
% %  carry out the landmark registration
% LandmarksMean = mean(Landmarks);
% [radiusfdLM, radiuswarpfdLM, WfdLM] = ...
%        landmarkreg(radiusfd_new, Landmarks, LandmarksMean, WfdParLM, 1);
%    
% %  plot registered accelerations along with warping functions
% radiusmatUR = eval_fd(agefine, radiusfd_new);
% radiusmatLM = eval_fd(agefine, radiusfdLM);
% radiuswarpmatLM  = eval_fd(agefine, radiuswarpfdLM);
% radiuswarpmatLM(1,:) = 0; radiuswarpmatLM(length(agefine),:) = 1;
% % figure, plot(agefine, radiusmatUR, 'LineWidth', 2);
% figure, plot(agefine, radiusmatLM, 'LineWidth', 2);
% rangeY = get(gca,'ylim');
% liney = linspace(rangeY(1), rangeY(2), 10);
% linex = liney;
% hold on
% for icase = 1:numf
%     warpInvfd = smooth_basis(radiuswarpmatLM(:,icase), agefine, fdPar_new);
%     warpedLM(icase,:) = eval_fd(Landmarks(icase,:), warpInvfd);
%     valueLM(icase, :) = eval_fd(warpedLM(icase,:), radiusfdLM(icase));
%     plot(warpedLM(icase, :), valueLM(icase, :), char(style(icase)));
% end
% for nLM = 1:nLandmarks 
%     linex(:) = LandmarksMean(nLM);
%     plot(linex, liney, 'm-');
% end
% hold off
% % figure, plot(agefine, radiuswarpmatLM);
% 
% index = 1:numf;
% nbasisw = 15;
% norder  =  5;
% basisw  = create_bspline_basis([0,1], nbasisw, norder);
% coef0 = zeros(nbasisw, length(index));
% Wfd0 = fd(coef0, basisw);
% Lfdobj = int2Lfd(2);
% lambda = 1;
% WfdPar = fdPar(Wfd0, Lfdobj, lambda);
% % register landmark-registered radius
% radiusLMmeanfd = mean(radiusfdLM);
% [radiusLMCregfd, radiusLMCwarpfd] = registerfd(radiusLMmeanfd, radiusfdLM, WfdPar);
% radiusLMCfine = eval_fd(agefine, radiusLMCregfd);
% radiusLMCwarpmat = eval_fd(agefine, radiusLMCwarpfd);
% radiusLMCwarpmat(1,:) = 0; radiusLMCwarpmat(length(agefine), :) = 1;
% figure, plot(agefine, radiusLMCfine, 'LineWidth', 2);
% hold on
% for icase = 1:numf
%     warpInvfd = smooth_basis(radiusLMCwarpmat(:,icase), agefine, fdPar_new);
%     warpedLMC(icase,:) = eval_fd(warpedLM(icase,:), warpInvfd);
%     valueLM(icase,:) = eval_fd(warpedLMC(icase,:), radiusLMCregfd(icase));
%     plot(warpedLMC(icase, :), valueLM(icase, :), char(style(icase)));
% end
% for nLM = 1:nLandmarks 
%     linex(:) = LandmarksMean(nLM);
%     plot(linex, liney, 'm-');
% end
% hold off
