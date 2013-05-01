%tmp = load( '../../results/Atlas30/Age/boxplotPart4/median_wbp.mat' );
%cases = tmp.median_case_wbp;
%ages = tmp.median_age_wbp;
%ages = ages(:, 2);

% add two cases: SGS03 pre and post surgery
%{
numf = length(cases);
numf_crl = numf;
numf = numf + 1;
cases{numf} = 'SGS03ManPre';
ages(numf) = 9;
numf = numf + 1;
cases{numf} = 'SGS03VirPost';
ages(numf) = 9;
numf = numf + 1;
cases{numf} = 'SGS03ManPost';
ages(numf) = 20;
%}

subject_info = readtext( '../../data/selected_subjects_CRL_SGS03.csv' );

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

% read the central region id
filename = '../../results/Atlas30/Age/boxplotPart4/idsInCentralRegion_SGS03ManPre.txt';
fid = fopen( filename, 'r' );
tline = fgets(fid);
numIds = sscanf(tline, '%d');
for iI = 1:numIds
    tline = fgets(fid);
    idTmp = sscanf( tline, '%d' );
    casesNew{iI} = cases{idTmp};
    agesNew(iI) = ages(idTmp);
    weightsNew(iI) = weights(idTmp);
end
fclose(fid);

clear cases ages weights
cases = casesNew;
ages = agesNew;
weights = weightsNew;
clear casesNew agesNew weightsNew

numf = length(cases);
numf_crl = numf;
numf = numf + 1;
cases{numf} = 'SGS03ManPre';
ages(numf) = 9;
numf = numf + 1;
cases{numf} = 'SGS03VirPost';
ages(numf) = 9;
numf = numf + 1;
cases{numf} = 'SGS03ManPost';
ages(numf) = 20;


% the color to display different ages
%nBins = ceil( max(ages)/10.0 ) * 10;
nBins = 200;
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
%figure('Colormap', cmap), hold on

inputPrefix = '../../results';
outputPrefix = '../../results/Atlas30/Slice';

for icase = 1:length(cases)
	landmarksIdFile = sprintf('%s/IsoSurface%s/%s_LandmarksIdOnCenterline.txt', inputPrefix, cases{icase}, cases{icase} );
	contourPrefix = sprintf('%s/Contour%s', inputPrefix, cases{icase} );
	%filenamePrefix = sprintf('%s/slice%02d', outputPrefix, icase );
    if icase <= numf_crl
        displaySliceForLandmarks( landmarksIdFile, contourPrefix, 'b', 0 );
    else
        displaySliceForLandmarks( landmarksIdFile, contourPrefix, cmap( ages(icase), :), 1 );
    end
end

%filenamePng = sprintf('%s/subglottic.png', outputPrefix);
%saveas(gca, filenamePng);
filenameFig = sprintf('%s/subglottic.fig', outputPrefix);
saveas(gca, filenameFig);
