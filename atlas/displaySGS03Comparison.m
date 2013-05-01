clear all
close all

% parameter settings
subject_info = readtext( '../../data/SGS03_Virtual.csv' );
data_prefix = '../../results';
results_prefix = '../../results/SGS03', 

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


for icase = 1:length(cases)
    tmp = cases{icase};
    tmpStr = sprintf( '%s/IsoSurface%s/%s_MeanAndNormal.txt', data_prefix, cases{icase}, cases{icase} );
    meanNormFile{icase} = tmpStr;
    tmpStr = sprintf( '%s/Contour%s/%s_Area.txt', data_prefix, cases{icase}, cases{icase} );
    areaFile{icase} = tmpStr;
    tmpStr = sprintf( '%s/IsoSurface%s/%s_LandmarksIdOnCenterline.txt', data_prefix, cases{icase}, cases{icase} );
    landmarksFile{icase} = tmpStr;
end

curveRegistration( meanNormFile, areaFile, landmarksFile, results_prefix, cases, ages, weights );