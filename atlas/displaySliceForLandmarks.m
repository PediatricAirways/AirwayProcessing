function score = displaySliceForLandmarks( landmarksIdFile, contourPrefix, clr, flag)

fid = fopen( landmarksIdFile, 'r' );
tline = fgets(fid);
count = sscanf( tline, '%d' );
landmarksId = zeros( count, 1 );
for iI = 1:count
	tline = fgets(fid);
	landmarksId(iI) = sscanf( tline, '%d' );
end
fclose(fid);

%handle = figure;
%set(handle, 'Position', [0 0 1200 240] );
selectedPoints = [landmarksId(length(landmarksId)-1), landmarksId(length(landmarksId)-1), landmarksId(length(landmarksId)-1)];
                %round( (landmarksId(length(landmarksId)-1) + landmarksId(length(landmarksId)))/2 )];
for iK = 1:length(selectedPoints)
    if iK == 1 && flag == 1 || iK == 2 && flag == 0
	continue;
    end
	n = selectedPoints(iK);
    % read contours
    fname =sprintf('%s/contour%03d.txt', contourPrefix, n);
    fid = fopen(char(fname), 'r');
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
        posTmp(iI*2 -1 ) = nStart;
        posTmp(iI*2) = curPos;
    end
    if curPos > 0
    	xmean = mean( A( :, 1 ) );
    	ymean = mean( A( :, 2 ) );
    	zmean = mean( A( :, 3 ) );
    	[coeff, score] = princomp( [ A( :, 1 ) - xmean, A( :, 2 ) - ymean, A( :, 3 ) - zmean ] );
    
        if flag == 0
            lineStyle = {'-'};
	    lineWidth = 2;
        else
            lineStyle = {'--'};
	    lineWidth = 2;
        end
    	for iTmp = 1:round(length(posTmp)/2)
        	subplot(1, length(selectedPoints), iK), 
            if max( score(:, 1) ) - min( score(:, 1) ) >= max( score(:, 2) ) - min( score(:, 2) )
                plot( [score( posTmp(iTmp*2-1) : posTmp(iTmp*2), 1 ); score(posTmp(iTmp*2-1), 1)], ...
                      [score( posTmp(iTmp*2-1) : posTmp(iTmp*2), 2 ); score(posTmp(iTmp*2-1), 2)], ...
                      'Color', clr, 'LineWidth', lineWidth, 'LineStyle', lineStyle{1} );
                %{  
                fill( [score( posTmp(iTmp*2-1) : posTmp(iTmp*2), 1 ); score(posTmp(iTmp*2-1), 1)], ...
                      [score( posTmp(iTmp*2-1) : posTmp(iTmp*2), 2 ); score(posTmp(iTmp*2-1), 2)], ...
                      [1.0 192/255 203/255] );
                %}
            else
                plot( [score( posTmp(iTmp*2-1) : posTmp(iTmp*2), 2 ); score(posTmp(iTmp*2-1), 2)], ...
                      [score( posTmp(iTmp*2-1) : posTmp(iTmp*2), 1 ); score(posTmp(iTmp*2-1), 1)], ...
                  'Color', clr, 'LineWidth', lineWidth, 'LineStyle', lineStyle{1} );
                
                 %{
                 fill( [score( posTmp(iTmp*2-1) : posTmp(iTmp*2), 2 ); score(posTmp(iTmp*2-1), 2)], ...
                       [score( posTmp(iTmp*2-1) : posTmp(iTmp*2), 1 ); score(posTmp(iTmp*2-1), 1)], ...
                       [1.0 192/255 203/255] );
                %}
            end
			axis equal
			hold on			
    	end
		%if iK == length(landmarksId)
			xlim([-10 10]);
			ylim([-10 10]);
		%end
    end
end

%{
filenamePng = sprintf('%s/subglottic.png', outputPrefix);
saveas(gca, filenamePng);
filenameFig = sprintf('%s/subglottic.fig', outputPrefix);
saveas(gca, filenameFig);
%}
