function drawLandmarksOnOneCurve( meanNormFile, areaFile, ellipseAreaFile, landmarksIdOnCenterlineFile, curveAndLandmarksFile, outputPrefix )

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
    
    % deal with some computing errors
%     for i = 2:nNumArea-1
%         if area(i-1, n) / area(i, n) > 2.5 & area(i+1, n) / area(i,n) > 2.5
%             area(i, n) = ( area(i-1,n) + area(i+1,n) ) / 2.0;
%         end
%     end
%     if area(nNumArea-1, n) / area(nNumArea, n) > 2.5
%         area(nNumArea, n) = area(nNumArea-1, n);
%     end
%     for i = 1:nNumArea
%         if area(i,n) > 600 
%             area(i,n) = area(i,n) / 2;
%         end
%     end
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
    %areafd = smooth_basis(dist(1:numP(n),n), area(1:numP(n),n), areafdPar);
    areafd = smooth_basis( xValue, yValue, areafdPar );
    
    clear yValue
    yValue = areaEllipse(1:numP(n), n);
    for iTmp = 1:100
        yValue = [ yValue(1); yValue; yValue(end) ];
    end 
    areaEllipsefd = smooth_basis( xValue, yValue, areafdPar );
    
    %areafd = smooth_basis(dist(1:numP(n),n), area(1:numP(n),n), areafdPar);
    %areaEllipsefd = smooth_basis(dist(1:numP(n),n), areaEllipse(1:numP(n),n), areafdPar);
    %figure, plotfit_fd(area(1:numP(n),n), dist(1:numP(n),n), areafd);
    %hold on
    %plotfit_fd(areaEllipse(1:numP(n),n), dist(1:numP(n),n), areaEllipsefd);
    %hold off
    %legend( 'AreaDots', 'AreaCurve', 'EllipseDots', 'EllipseCurve' );
    %filename_area = sprintf('%s/areaPlot%02d.png', outputPrefix, n);
    %saveas( gcf, filename_area );
    areafine(:,n) = eval_fd(agefine, areafd);
    areaEllipsefine(:,n) = eval_fd(agefine, areaEllipsefd);
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

for icase = 1:numf
    valueLM(icase, :) = eval_fd( Landmarks(icase, :), areafd_new(icase));
    valueLM_ellipse(icase, :) = eval_fd( Landmarks(icase, :), areaEllipsefd_new(icase) );
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
xlabel( 'Depth Along Centerline', 'FontSize', 20, 'FontWeight', 'Bold' );
ylabel( 'Cross Sectional Area', 'FontSize', 20, 'FontWeight', 'Bold' );
ylim( [0 1000] );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
saveas( h, curveAndLandmarksFile );

filename = sprintf( '%s/curveWithLandmarks.png', outputPrefix );
saveas( h, filename );

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

% nLandmarks = size( Landmarks, 2 );
% 
% wbasisLM = create_bspline_basis([0,1], max(nLandmarks+3-2,4), 3);
% WfdLM    = fd(zeros(max(nLandmarks+3-2,4),1),wbasisLM);
% WfdParLM = fdPar(WfdLM,1,1e-12);
% 
% %  carry out the landmark registration
% LandmarksMean = mean(Landmarks);
% [areafdLM, areawarpfdLM, WfdLM] = ...
%        landmarkreg(areafd_new, Landmarks, LandmarksMean, WfdParLM, 1);
%    
% %  plot registered accelerations along with warping functions
% areamatUR = eval_fd(agefine, areafd_new);
% areamatLM = eval_fd(agefine, areafdLM);
% areawarpmatLM  = eval_fd(agefine, areawarpfdLM);
% areawarpmatLM(1,:) = 0; areawarpmatLM(length(agefine),:) = 1;
% % figure, plot(agefine, radiusmatUR, 'LineWidth', 2);
% figure, plot(agefine, areamatLM, 'LineWidth', 2);
% rangeY = get(gca,'ylim');
% liney = linspace(rangeY(1), rangeY(2), 10);
% linex = liney;
% hold on
% plot(agefine, mean(areamatLM, 2), 'm--', 'LineWidth', 3);
% for icase = 1:numf
%     warpInvfd = smooth_basis(areawarpmatLM(:,icase), agefine, fdPar_new);
%     warpedLM(icase,:) = eval_fd(Landmarks(icase,:), warpInvfd);
%     valueLM(icase, :) = eval_fd(warpedLM(icase,:), areafdLM(icase));
%     plot(warpedLM(icase, :), valueLM(icase, :), char(style(icase)), 'LineWidth', 10);
% end
% warpedLMSGS08 = eval_fd( LandmarksSGS08, warpInvfd );
% valueLMSGS08 = eval_fd( warpedLMSGS08, areafdLM(3) );
% plot( warpedLMSGS08, valueLMSGS08, 'ro', 'LineWidth', 5 );
% % for nLM = 1:nLandmarks 
% %     linex(:) = LandmarksMean(nLM);
% %     plot(linex, liney, 'm-');
% % end
% hold off
% xlabel( 'Position' );
% ylabel( 'Area' );
% legend( legends );

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
