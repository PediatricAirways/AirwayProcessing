%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Pediatric Airway Measurement Code  ---
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

clear all
close all

%% Initialisation of POI Libs
% Add Java POI Libs to matlab javapath
javaaddpath('poi_library/poi-3.8-20120326.jar');
javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
javaaddpath('poi_library/xmlbeans-2.3.0.jar');
javaaddpath('poi_library/dom4j-1.6.1.jar');
javaaddpath('poi_library/stax-api-1.0.1.jar');

% the form with the information for all cases, including the id, age, gender and weight
subject_info = readtext( '../data/selected_subjects_CRL_SGS_Carina.csv' );
% the location of the results after airway processing
data_prefix = '../results';
% the location of the measurement output
results_prefix = '.';
% the id where SGS starts, the cases before it are all controls
SGS_startId = 69;

% this is the location of the atlas output
atlas_score = '../results_perimeter/Atlas24_68_ScorewrtMin/Age/boxplotPart4/diagnosis_data';

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
    tmpStr = sprintf( '%s/IsoSurface%s/%s_MeanAndNormal.txt', data_prefix, cases{icase}, cases{icase} );
	meanNormFile{icase} = tmpStr;
	tmpStr = sprintf( '%s/Contour%s/%s_Area.txt', data_prefix, cases{icase}, cases{icase} );
	areaFile{icase} = tmpStr;
	tmpStr = sprintf( '%s/IsoSurface%s/%s_LandmarksIdOnCenterline.txt', data_prefix, cases{icase}, cases{icase} );
	landmarksFile{icase} = tmpStr;
end

[areaLandmarks, midPos] = measureArea( meanNormFile,  areaFile, landmarksFile, SGS_startId, cases, atlas_score);
save( 'areaLandmarks.mat', 'areaLandmarks' );

% plot based on ages or weights
outputResult = sprintf('%s/age', results_prefix);
plotRatio( areaLandmarks, midPos, cases, ages, SGS_startId, outputResult );
outputResult = sprintf('%s/weight', results_prefix);
plotRatio( areaLandmarks, midPos, cases, weights, SGS_startId, outputResult );

