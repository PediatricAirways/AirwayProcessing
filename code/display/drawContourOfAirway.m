function drawContourOfAirway( meanPointsFile, normCheckFile, contourPrefix, outputPrefix, ellipseAreaFile, landmarksFile )

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(meanPointsFile, 'r');
tline = fgets(fid);
count = sscanf(tline, '%d');
for i=1:count
    tline = fgets(fid);
    newMeanPoints(:, i) = sscanf(tline, '%f');
end
fclose(fid);

fid = fopen( normCheckFile, 'r' );
tline = fgets(fid);
nNormCheck = sscanf( tline, '%d' );
for i=1:nNormCheck 
    tline = fgets(fid);
    normCheck(i) = sscanf(tline, '%d' );
end

if count ~= nNormCheck
    error( 'The number of points on the centerline is not kept consistent.\n' );
end

handle = figure(1), hold on
plot3( newMeanPoints(1, :), newMeanPoints(2, :), newMeanPoints(3, :), 'ro', 'LineWidth', 2 );


%{
for i=1:count
    if normCheck(i) == 1
        plot3(newMeanPoints(1, i), newMeanPoints(2, i), newMeanPoints(3, i), 'ro', 'LineWidth', 2);    
    else
        plot3(newMeanPoints(1, i), newMeanPoints(2, i), newMeanPoints(3, i), 'co', 'LineWidth', 2); 
    end
end
%}

%fidArea = fopen( areaFile, 'wt');
%fprintf( fidArea, '%d\n', count );

fidAreaEllipse = fopen( ellipseAreaFile, 'wt' );
fprintf( fidAreaEllipse, '%d\n', count );

fileNameEllipseInfo = sprintf( '%sEllipseInfo.txt', outputPrefix );
fileNameEllipseInfo
fidEllipseInfo = fopen( fileNameEllipseInfo, 'wt' )
fprintf( fidEllipseInfo, '%d\n', count );

h_ellipse = figure(2), 
hold on 
for n = 1:count
    % read contours
    fname =sprintf('%s%03d.txt', contourPrefix, n);
    fid = fopen(char(fname), 'r');
    if fid == -1 
	fprintf( fidAreaEllipse, '%f\n', 0 );
        % output ellipse info
        fprintf( fidEllipseInfo, '%f %f %f %f %f\n', 0, 0, 0, 0, 0 );
	continue;
    end
    tline = fgets(fid);
    num = sscanf(tline, '%d');
    nbCluster = zeros( num(2), 1 );
    clear A
    curPos = 0;
    clear posTmp;
    for iI = 1:num(2)
        if feof(fid)
            break;
        end
        tline = fgets(fid);
        nbCluster(iI) = sscanf(tline, '%d');
        for iJ = 1:nbCluster(iI)
            tline = fgets(fid);
            curPos = curPos + 1;
            A(curPos, :) = sscanf(tline, '%lf %lf %lf');
        end
        nStart = curPos - nbCluster(iI) + 1;
        figure(1), plot3( [A( nStart:curPos, 1 );A(nStart,1)], [A( nStart:curPos, 2);A(nStart,2)] , [A( nStart:curPos, 3 );A(nStart,3)] );
        hold on
        
%         hTmp = figure(100), 
%         xmeanTmp = mean( [A( nStart:curPos, 1 );A(nStart,1)] );
%         ymeanTmp = mean( [A( nStart:curPos, 2 );A(nStart,2)] );
%         zmeanTmp = mean( [A( nStart:curPos, 3 );A(nStart,3)] );
%         [coeffTmp, scoreTmp] = princomp( [ [A( nStart:curPos, 1 );A(nStart,1)] - xmeanTmp, ...
%                            [A( nStart:curPos, 2 );A(nStart,2)] - ymeanTmp, [A( nStart:curPos, 1 );A(nStart,1)] - zmeanTmp ] );
%         plot( scoreTmp(:, 1), scoreTmp(:, 2) );
%         hold on
      
        posTmp(iI*2 -1 ) = nStart;
        posTmp(iI*2) = curPos;
    end
    fclose(fid);
    
    if curPos > 0
    	% compute the ellipse
    	xmean = mean( A( :, 1 ) );
    	ymean = mean( A( :, 2 ) );
    	zmean = mean( A( :, 3 ) );
    	[coeff, score] = princomp( [ A( :, 1 ) - xmean, A( :, 2 ) - ymean, A( :, 3 ) - zmean ] );
    
    	hTmp = figure(100);
    	hold on
    	for iTmp = 1:round(length(posTmp)/2)
        	plot( score( posTmp(iTmp*2-1) : posTmp(iTmp*2), 1 ), score( posTmp(iTmp*2-1) : posTmp(iTmp*2), 2 ), 'LineWidth', 2 );
    	end
	set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
    	filenameTmp = sprintf( '%s_contours/fig_contour%03d.png', outputPrefix, n );
    	saveas( hTmp, filenameTmp );
    	close( hTmp );
    
    
    	[a_long, b_short, pts_center, pts_t] = getEllipse2DTo3D( score(:, 1:2), coeff, xmean, ymean, zmean );
    	figure(2), plot3( pts_t( :, 1 ), pts_t( :, 2 ), pts_t( :, 3 ) );
    	hold on
	figure(2), plot3( pts_center(1), pts_center(2), pts_center(3), 'ro', 'LineWidth', 2 );
	hold on
	fprintf( fidAreaEllipse, '%f\n', pi * a_long * b_short );

	% output ellipse info
	fprintf( fidEllipseInfo, '%f %f %f %f %f\n', pts_center(1), pts_center(2), pts_center(3), a_long, b_short );
	clear score coeff
    end
end

fclose( fidAreaEllipse );
fclose( fidEllipseInfo );

% read landmarks
%[key, val] = textread( landmarksFile, '%s:%[^\n]' );
%for iI = 1:size( key, 1 )
%    landmarks_pnt( iI, 1:3 ) = str2num( val{iI} );
%end
landmarks = zeros( 6, 3 );
flagLM = zeros(6, 1);
[key, val] = textread( landmarksFile, '%s:%[^\n]' );
for iI = 1:size( key, 1 )
    landmarksTmp = str2num( val{iI} );
    landmarksIdTmp = 0;
    if strcmp( key{iI}, 'NasalSpine' )
	landmarksIdTmp = 1;
    elseif strcmp( key{iI}, 'PosteriorInferiorVomerCorner' )
        landmarksIdTmp = 2;
    elseif strcmp( key{iI}, 'TracheaCarina' )
        landmarksIdTmp = 6;
    elseif strcmp( key{iI}, 'EpiglottisTip' )
        landmarksIdTmp = 3;
    elseif strcmp( key{iI}, 'TVC' )
        landmarksIdTmp = 4;
    elseif strcmp( key{iI}, 'Subglottic' )
		landmarksIdTmp = 5;
    end
    if landmarksIdTmp >= 1 && landmarksIdTmp <= size(key, 1)
        landmarks_pnt( landmarksIdTmp, : ) = landmarksTmp;
		flagLM( landmarksIdTmp ) = 1;
    end
end
for iI = 1:length(flagLM)
	if flagLM( iI ) == 0 
		landmarks_pnt( iI, : ) = [];
	end
end
figure(1), plot3( landmarks_pnt(:,1), landmarks_pnt(:,2), landmarks_pnt(:,3), 'gx', 'LineWidth', 5 );
%axis( [ xmin xmax ymin ymax zmin zmax ] );
view( 90, 0 );
axis equal
axis tight

figure(2), plot3( landmarks_pnt(:,1), landmarks_pnt(:,2), landmarks_pnt(:,3), 'gx', 'LineWidth', 5 );
%axis( [ xmin xmax ymin ymax zmin zmax] );
view( 90, 0 );
axis equal
axis tight

hold off
hold off

filenamePng = sprintf('%s.png', outputPrefix);
saveas(handle, filenamePng);
filenameFig = sprintf('%s.fig', outputPrefix);
saveas(handle, filenameFig);
filenamePng = sprintf( '%s_ellipse.png', outputPrefix );
saveas( h_ellipse, filenamePng );
filenameFig = sprintf( '%s_ellipse.fig', outputPrefix );
saveas( h_ellipse, filenameFig );
