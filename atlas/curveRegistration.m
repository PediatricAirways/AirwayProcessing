function curveRegistration( meanNormFile, areaFile, landmarksIdOnCenterlineFile, ...
    outputPrefix, cases, ages, weights )

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


% read area file
for n = 1:numf
    fid = fopen(areaFile{n}, 'r');
    tline = fgets(fid);
    nNumArea = sscanf(tline, '%d');
    for i = 1:nNumArea
        tline = fgets(fid);
        if(tline == -1) break; end
        area(i,n) = sscanf(tline, '%f');
        % deal with the missed part in the virtual surgery result
        %if(area(i,n) <= 0 && n > 1) 
        %    area(i,n) = area(i, n-1);
        %end
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

% smooth the area before TVC
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

    areafd_tmp = smooth_basis( xValue, yValue, areafdPar );
    areaSmooth(1:LandmarksId(n, 4)-1, n) = eval_fd(dist(1:LandmarksId(n, 4)-1, n), areafd_tmp);
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
    lambda   = 0.5* 1e-6;  % adjust this number to make the curve fit the scattered point
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

    areafd = smooth_basis( xValue, yValue, areafdPar );
    areafine(:,n) = eval_fd(agefine, areafd);
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

% find the value for landmarks on curves
for icase = 1:numf
    valueLM(icase, :) = eval_fd( Landmarks(icase, :), areafd_new(icase));
end

% draw curve and landmarks
style = {'yx', 'rx', 'gx', 'bx', 'cx'};
areafine_new = eval_fd(agefine, areafd_new);

% the color to display different ages
nBins = ceil( max(ages)/10.0 ) * 10;
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
    if mod(icase, 3) == 1
        plot( agefine, areafine_new( :, icase ), 'Color', cmap(round(ages(icase)), :), 'LineStyle', '-', 'LineWidth', 2 );
    else
        plot( agefine, areafine_new( :, icase ), 'Color', cmap(round(ages(icase)), :), 'LineStyle', '--', 'LineWidth', 2 );
    end
end
for icase = 1:numf
    plot( dist(1:numP(icase), icase), area(1:numP(icase), icase), '.', 'Color', cmap(round(ages(icase)), :) );
end
% display landmarks on curves
for icase = 1:size(Landmarks, 2)
     plot( Landmarks(:, icase), valueLM(:, icase), char(style(mod(icase, 5)+1)), 'LineWidth', 3);
end
hold off
caxis( [0 nBins] );
h = colorbar( 'peer', gca );
set( get( h, 'ylabel' ), 'String', 'Age: month(s)', 'FontSize', 24, 'FontWeight', 'Normal' );
set( h, 'ylim', [0, nBins] );
set( gca, 'XTick', 0:0.1:1 );
xlabel( 'Depth along the centerline', 'FontSize', 24, 'FontWeight', 'Normal' );
ylabel( 'Cross-sectional area (mm^2)', 'FontSize', 24, 'FontWeight', 'Normal' );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
legend( cases );
title( 'Unregistered Curves', 'FontSize', 24, 'FontWeight', 'Normal' );
filename = sprintf( '%s/unregistered_curves.png', outputPrefix );
saveas( gca, filename );
filename = sprintf( '%s/unregistered_curves.fig', outputPrefix );
saveas( gca, filename );


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
    
% %  plot registered area along with warping functions
areamatUR = eval_fd(agefine, areafd_new);
areamatLM = eval_fd(agefine, areafdLM);
areawarpmatLM  = eval_fd(agefine, areawarpfdLM);
areawarpmatLM(1,:) = 0; areawarpmatLM(length(agefine),:) = 1;

figure('Colormap', cmap), hold on
for icase = 1:numf
    plot( agefine, areawarpmatLM( :, icase ), 'Color', cmap(round(ages(icase)), :), 'LineWidth', 2);
end
for icase = 1:numf
    for iLM = 1:size( Landmarks, 2 )
	plot(mean(Landmarks(:,iLM)), Landmarks(icase, iLM), char(style(mod(iLM, 5)+1)), 'LineWidth', 3);
    end
end
ylabel( 'Physical position of the landmark', 'FontSize', 24, 'FontWeight', 'Normal' );
xlabel( 'Mean position of the landmark', 'FontSize', 24, 'FontWeight', 'Normal' );
title( 'Warping function', 'FontSize', 24, 'FontWeight', 'Normal' );
set( gca, 'XTick', 0:0.1:1 );
set( gca, 'YTick', 0:0.1:1 );
caxis( [0 nBins] );
h = colorbar;
set( get( h, 'ylabel' ), 'String', 'Age: month(s)', 'FontSize', 24, 'FontWeight', 'Normal' );
set( h, 'ylim', [0, nBins] );
set( gca, 'FontSize', 24, 'FontWeight', 'Normal' );
legend( cases );
filename = sprintf( '%s/warpFunctionLandmarks.png', outputPrefix );
saveas( gca, filename );


% display all the registered curves, colored by age
figure('Colormap', cmap); hold on
for icase = 1:numf
    if mod( icase, 3 ) == 1
        plot( agefine, areamatLM( :, icase ), 'Color', cmap(round(ages(icase)), :), 'LineStyle', '-', 'LineWidth', 2 );
    else
        plot( agefine, areamatLM( :, icase ), 'Color', cmap(round(ages(icase)), :), 'LineStyle', '--', 'LineWidth', 2 );
    end
end
legend(cases);
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
set( get( h, 'ylabel' ), 'String', 'Age: month(s)', 'FontSize', 24, 'FontWeight', 'Normal' );
set( h, 'ylim', [0 nBins] );
xlabel( 'Depth along the centerline', 'FontSize', 24, 'FontWeight', 'Normal' );
ylabel( 'Cross-sectional area (mm^2)', 'FontSize', 24, 'FontWeight', 'Normal' );
set( gca, 'XTick', 0:0.1:1 );
set( gca, 'FontSize', 24, 'FontWeight', 'Normal' );
title( 'Registered Curves', 'FontSize', 24, 'FontWeight', 'Normal' );
filename = sprintf( '%s/landmark_registration_age.png', outputPrefix );
saveas( gca, filename );
filename = sprintf( '%s/landmark_registration_age.fig', outputPrefix );
saveas( gca, filename );

% generate aligned the centerlines
samples = linspace(0, 1, 100)';
samples_org = eval_fd( samples, areawarpfdLM );
centerlines = zeros( numf, length(samples), 6 );
for icase = 1:numf 
    startId = 1;
    endId = numP(icase);
    for iI = 1:length(samples_org(:, icase))
        while startId < endId && samples_org(iI, icase) > dist(startId, icase)
            startId = startId + 1;
            continue;
        end
        if startId <= 1
            centerlines( icase, iI, 1:3 ) = A(:, startId, icase); 
            centerlines( icase, iI, 4:6 ) = normCases(:, startId, icase);
        else
            percent = ( dist(startId, icase) - samples_org(iI, icase) ) / ...
                (dist(startId, icase) - dist(startId-1, icase));
            centerlines( icase, iI, 1:3) = percent * A(:, startId-1, icase) + ...
                (1-percent) * A(:, startId, icase);
            
            angle = acos( dot(normCases(:, startId-1, icase), normCases(:, startId, icase) ) ...
                / (norm(normCases(:, startId-1, icase)) * norm(normCases(:, startId, icase) ) ) );
            if angle > pi/2 
                normCases(:, startId, icase) = -normCases(:, startId, icase);
            end
            centerlines( icase, iI, 4:6) = percent * normCases(:, startId-1, icase) + ...
                (1-percent) * normCases(:, startId, icase);
        end
    end
    figure, plot3(A(1, 1:numP(icase), icase), A(2, 1:numP(icase), icase), A(3, 1:numP(icase), icase), 'bo', 'LineWidth', 2 );
    hold on
    plot3(centerlines(icase, :, 1), centerlines(icase, :, 2), centerlines(icase, :, 3), 'r*', 'LineWidth', 2);
    for iPnt = 1:size(centerlines, 2)
        plot3( [centerlines(icase, iPnt, 1), centerlines(icase, iPnt, 1) + 5*centerlines(icase, iPnt, 4) ], ...
               [centerlines(icase, iPnt, 2), centerlines(icase, iPnt, 2) + 5*centerlines(icase, iPnt, 5) ], ...
               [centerlines(icase, iPnt, 3), centerlines(icase, iPnt, 3) + 5*centerlines(icase, iPnt, 6) ], 'm-' );
     
    end
    hold off
    axis equal
    fname = sprintf( '%s/centerlineAligned_%s.txt', outputPrefix, cases{icase} );
    fid = fopen( fname, 'wt' );
    fprintf( fid, '%d\n', size(centerlines, 2) );
    for iPnt = 1:size(centerlines, 2) 
        fprintf( fid, '%f %f %f %f %f %f\n', centerlines(icase, iPnt, 1), centerlines(icase, iPnt, 2), ...
            centerlines(icase, iPnt, 3), centerlines(icase, iPnt, 4), centerlines(icase, iPnt, 5), centerlines(icase, iPnt, 6) );
    end
    fclose(fid);
end
