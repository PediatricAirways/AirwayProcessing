
clear all
close all

% parameter settings
subject_info = readtext( '../../data/selected_subjects_CRL_SGS03.csv' );
data_prefix = '../../results';
results_prefix = '../../results/Atlas30', 
gaussian_sigma = 30;

% read the subject information from the csv file of data flow
[nRow, nCol] = size( subject_info );
for iI = 1:nCol
    if strcmp( subject_info{1, iI}, 'PatientId' ) == 1
        idPatient = iI;
    end
    if strcmp( subject_info{1, iI}, 'Age' ) == 1
        idAge = iI;
    end
    if strcmp( subject_info{1, iI}, 'Weight' ) == 1
        idWeight = iI;
    end
end

cases = {};
ages = {};
weights = {};
for icase = 1:nRow-1
    cases{icase} = subject_info{icase+1, idPatient};
    ages{icase} = subject_info{icase+1, idAge};
    weights{icase} = subject_info{icase+1, idWeight};
end

ages = cell2mat( ages );
weights = cell2mat( weights );

% display ages distribution
figure, hold on
tmp = gaussmf( ages, [gaussian_sigma, floor(max(ages)/2)] );
[tmpAges, tmpId] = sort( ages );
tmp = tmp( tmpId );
%plot( tmpAges, tmp, '--r*', 'LineWidth', 2);
hist( ages, ceil(max(ages)/10.0)*10 );
title_text = sprintf( 'Histogram of %d subjects'' ages', length(ages) );
title( title_text, 'FontSize', 20, 'FontWeight', 'Bold' );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
%legend_text = sprintf( 'Gaussian = [%d, %d]', gaussian_sigma, floor(max(ages)/2) );
%legend( legend_text );
xlabel( 'Age: month(s)' );
ylabel( 'Subject number' );
filename = sprintf( '%s/ageHist.fig', results_prefix );
saveas( gca, filename ); 

for icase = 1:length(cases)
    tmp = cases{icase};
    tmpStr = sprintf( '%s/IsoSurface%s/%s_MeanAndNormal.txt', data_prefix, cases{icase}, cases{icase} );
    meanNormFile{icase} = tmpStr;
    tmpStr = sprintf( '%s/Contour%s/%s_Area.txt', data_prefix, cases{icase}, cases{icase} );
    areaFile{icase} = tmpStr;
    tmpStr = sprintf( '%s/Curves%s/%s_AreaEllipse.txt', data_prefix, cases{icase}, cases{icase} );
    ellipseAreaFile{icase} = tmpStr;
    tmpStr = sprintf( '%s/IsoSurface%s/%s_LandmarksIdOnCenterline.txt', data_prefix, cases{icase}, cases{icase} );
    landmarksFile{icase} = tmpStr;
end

pos_SGS = 53;
ages_CRL = ages( 1:pos_SGS-1 );
%for iFactor = 1:17
	% calculate the compariable uniform window 
	x = 1:max(ages_CRL);
	fx_gaussian = zeros( 1, length(x));
	%gaussian_sigma = 18*(sqrt(2).^(iFactor-1));
	%results_prefix_sub = sprintf( '%s/test%02d', results_prefix, iFactor );
	for iI = 1:length(ages_CRL)
	    tmp = gaussmf( x, [gaussian_sigma, ages_CRL(iI)] );
	    fx_gaussian = fx_gaussian + tmp;
	end
	fx_gaussian = fx_gaussian ./ length(ages_CRL);
	figure, hist( ages_CRL, ceil(max(ages_CRL)/10.0)*10 );
	hold on
	plot( x, fx_gaussian, 'r--', 'LineWidth', 2);

	winSizeMatched = 0;
	minValueMatched = -1;
	for winSize = gaussian_sigma:2*gaussian_sigma
		fx_uniform = zeros( 1, length(x) );
		for iI = 1:length(ages_CRL)
		    pickedCurves = ( (x >= ages_CRL(iI) - winSize) .* (x <= ages_CRL(iI) + winSize) );
		    fx_uniform = fx_uniform + pickedCurves;
		end
		fx_uniform = fx_uniform ./ length(ages_CRL);
		tmpValue = norm( fx_gaussian-fx_uniform, 2 );
		if minValueMatched < 0 || tmpValue < minValueMatched
			winSizeMatched = winSize;
			minValueMatched = tmpValue;
		end
	end
	fx_uniform = zeros( 1, length(x) );
	for iI = 1:length(ages_CRL)
		pickedCurves = ( (x >= ages_CRL(iI) - winSizeMatched) .* (x <= ages_CRL(iI) + winSizeMatched) );
		fx_uniform = fx_uniform + pickedCurves;
	end
	fx_uniform = fx_uniform ./ length(ages_CRL);
	plot( x, fx_uniform, 'b-.', 'LineWidth', 2 );
	hold off
	title_name = sprintf( 'Gaussian window size: %f, uniform window size: %d', gaussian_sigma, winSizeMatched );
	title( title_name, 'FontSize', 15, 'FontWeight', 'Bold' );
	%filename = sprintf( '%s/windowSize.png', results_prefix_sub );
	filename = sprintf( '%s/windowSize.png', results_prefix );
	saveas( gca, filename );

	curveAndLandmarksFile = sprintf('%s/multCurvesAndLandmarks.fig', results_prefix);

	drawLandmarksOnMultiCurvesAndSGS( meanNormFile, areaFile, ellipseAreaFile, landmarksFile, curveAndLandmarksFile, results_prefix, cases, ages, weights, gaussian_sigma, winSizeMatched, pos_SGS );
%end

