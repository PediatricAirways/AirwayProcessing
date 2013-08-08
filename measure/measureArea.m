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

function [areaLandmarks, midPos] = measureArea( meanNormFile,  areaFile, landmarksFile, SGS_StartId, cases, atlas_score_dir)
% meanNormFile: the file names of centerline for all cases
% areaFile: the file names of cross-sectional area for all cases
% landmarksFile: the file names of landmark id for all cases
% SGS_StartId: the id where SGS starts
% cases: the patient ids
% atlas_score_dir: the location of the atlas output

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
    fidLandmarksIdFile = fopen( landmarksFile{n}, 'r' );
    tline = fgets( fidLandmarksIdFile );
    nLandmarksId = sscanf( tline, '%d' );
    for iI = 1:nLandmarksId
        tline = fgets( fidLandmarksIdFile );
        landmarksIdOnCenterline = sscanf( tline, '%d' );
        LandmarksId( n, iI ) = landmarksIdOnCenterline;
        Landmarks(n, iI) = dist( landmarksIdOnCenterline, n);
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
        area(i, n) = sscanf(tline, '%f');
    end
    fclose(fid);
end


addpath ('../code/display/Matlabfunctions/fdaM')
%  set up a basis for the functions W(t) that define the warping
%  functions
%  create curves to fit the area of contour and of ellipse
areaLandmarks = zeros( numf, 7 ); % subglottic, mid of trachea carina, ratio

filenameXLS = ('geometry_measurments_XA.xlsx');
headXLS = {'PatientId', 'XA_TVC', 'XA_Subglottis', 'XA_MidTrachea', 'XA_Ratio', 'XA_Atlas_Score', 'XA_Ratio_Score', 'XA_Stenosis'};
dataXLS = cell(numf, length(headXLS));
for n = 1:numf
    rng      = [0,1];
    knots    = dist(1:numP(n),n)';
    norder   = 6;
    nbasis   = length(knots) + norder - 2;
    areabasis = create_bspline_basis(rng, nbasis, norder, knots);
    Lfdobj   = int2Lfd(2);
    lambda   = 1e-6;
    areafdPar = fdPar(areabasis, Lfdobj, lambda);

    % curve fitting function for the area
    % try to make the first and last points on the curve.
    xValue = dist(1:numP(n), n);
    xValue(1) = xValue(1) + 100 * 1e-7;
    xValue(end) = xValue(end) - 100 * 1e-7;
    clear yValue
    yValue = area(1:numP(n), n);
    for iTmp = 1:100
        xValue = [ xValue(1) - 1e-7; xValue ; xValue(end) + 1e-7 ];
        yValue = [ yValue(1); yValue; yValue(end) ];
    end
    areafd = smooth_basis( xValue, yValue, areafdPar );
    
    if n < SGS_StartId   % deal with controls
        % measure 1.5mm around suglottics
        length_margin_suglottics = 1.5 / dist_copy(numP(n), n);
        areaLandmarks(n, 1) = sum( eval_fd( linspace(dist(LandmarksId(n, 5), n) - length_margin_suglottics/2, ...
                    dist(LandmarksId(n, 5), n) + length_margin_suglottics/2, 10 ), areafd ) ) / 10;

        % measure 15mm around middle of subglottis and carina
        length_margin = 15 / dist_copy(numP(n), n);
        midTrachea = dist(LandmarksId(n, 5), n) + ( dist(LandmarksId(n, 6), n) - dist(LandmarksId(n, 5), n) )/2;
        [midValue, midId] = min( abs( dist(:, n) - midTrachea ) );
        midPos(n, 1:3) = A(1:3, midId, n)';
        areaLandmarks(n, 2) = sum( eval_fd( linspace( max(midTrachea - length_margin/2, 0), min(midTrachea + length_margin/2, 1),  10 ), areafd ) ) / 10;

	% ratio of subglottics and mid-trachea
        areaLandmarks(n, 3) = areaLandmarks(n, 1) / areaLandmarks(n, 2);
        
	% keep the same for controls, will be different for SGSs
        areaLandmarks(n, 4) = areaLandmarks(n, 2);    
        areaLandmarks(n, 5) = areaLandmarks(n, 3);
        
    else
        % find the stenosis location, no subglottic landmark, remove it if it exists, should have only 5 landmarks
        if LandmarksId(n, 6) > 0
            LandmarksId(n, 5) = LandmarksId(n, 6);
        end

	% 1. take the dip region with minimal area as a potential stenosis
        stenosisId = -1;
        areaValueTmp = eval_fd( dist(:, n), areafd );
        
        for iI = LandmarksId(n, 4)+1 : LandmarksId(n, 5)-1
            if areaValueTmp(iI) <= areaValueTmp(iI-1) && areaValueTmp(iI) <= areaValueTmp(iI+1)
                if stenosisId < 0 || areaValueTmp(iI) < stenosisValue
                    stenosisId = iI;
                    stenosisValue = areaValueTmp(iI);
                end
            end
        end
        cases{n}
        stenosisId

	% 2. take the global minimal as another potential stenosis
        [valueTmp, idTmp] = min( areaValueTmp(LandmarksId(n, 4):LandmarksId(n, 5)) );
        stenosisId2 = idTmp + LandmarksId(n, 4) - 1;
        
        if stenosisId == stenosisId2   % we found the correct stenosis location
            % do nothing
        else
            if stenosisId > 0.5 * ( LandmarksId(n, 4) + LandmarksId(n, 5) )    % outside of the subglottic region
                stenosisId = stenosisId2;
                stenosisValue = valueTmp;
            end
        end
	% write the 3D points into stenosis
	posStenosis = sprintf('(%.2f, %.2f, %.2f)', A(1, stenosisId, n), A(2, stenosisId, n), A(3, stenosisId, n) );
	dataXLS(n, end) = {posStenosis};

        %stenosisId
        %stenosisId2
        
        areaLandmarks(n, 1) = stenosisValue;
        
        % compute the middle region 
        
        length_margin = 15 / dist_copy(numP(n), n);
        midTrachea = dist(stenosisId, n) + ( dist(LandmarksId(n, 5), n) - dist(stenosisId, n) )/2;
        [midValue, midId] = min( abs( dist(:, n) - midTrachea ) );
        midPos(n, 1:3) = A(1:3, midId, n)';
        areaLandmarks(n, 2) = sum( eval_fd( linspace( max(midTrachea - length_margin/2, 0), min(midTrachea + length_margin/2, 1),  10 ), areafd ) ) / 10;
        areaLandmarks(n, 3) = areaLandmarks(n, 1) / areaLandmarks(n, 2);

        
        % compute the ramp after stenosis
        for iI = stenosisId+1:LandmarksId(n, 5)
            if areaValueTmp(iI) > areaValueTmp(iI-1)
                continue;
            else
                break;
            end
        end
        areaLandmarks(n, 4) = areaValueTmp(iI-1);
 
        if iI-1 > (stenosisId + LandmarksId(n, 5))/2
            length_margin = 15 / dist_copy(numP(n), n);
            midTrachea = dist(stenosisId, n) + ( dist(LandmarksId(n, 5), n) - dist(stenosisId, n) )/2;
            [midValue, midId] = min( abs( dist(:, n) - midTrachea ) );
            midPos(n, 1:3) = A(1:3, midId, n)';
            areaLandmarks(n, 4) = sum( eval_fd( linspace( max(midTrachea - length_margin/2, 0), min(midTrachea + length_margin/2, 1),  10 ), areafd ) ) / 10;
        end
       
	% the ratio of stenosis and ramp
        areaLandmarks(n, 5) = areaLandmarks(n, 1) / areaLandmarks(n, 4);
        
    end

    % store cross-sectional area to xlsx file
    dataXLS(n, 1) = {cases{n}};   % patient Id
    dataXLS(n, 2) = num2cell(round(eval_fd(dist(LandmarksId(n, 4), n), areafd)*100)/100.0);   % TVC
    dataXLS(n, 3) = num2cell(round(areaLandmarks(n, 1)*100)/100.0);   % Subglottic
    dataXLS(n, 4) = num2cell(round(areaLandmarks(n, 2)*100)/100.0);   % Mid-Trachea
    dataXLS(n, 5) = num2cell(round(areaLandmarks(n, 3)*100)/100.0);   % Ratio
end

% load cases and scores based on atlas, if you have the results, uncomment it
%{
filename = sprintf('%s/cases.mat', atlas_score_dir);
tmp = load(filename);
cases_atlas = tmp.cases;

filename = sprintf('%s/posInAtlas.mat', atlas_score_dir);
tmp = load(filename);
posInAtlas = tmp.posInAtlas;
score_atlas = posInAtlas(:, 5);

for n = 1:numf
    for iI = 1:length(cases_atlas)
        if strcmp( cases{n}, cases_atlas{iI} )
            dataXLS(n, 6) = num2cell(round(score_atlas(iI)*100)/100.0);   % score based on atlas
            break;
        end
    end
end
%}

% convert the ratio to score 
areaLandmarks(:, 6) = max(1 - areaLandmarks(:, 3) ./ median(areaLandmarks(1:SGS_StartId-1, 3) ), 0);
%areaLandmarks(:, 7) = max(1 - areaLandmarks(:, 5 ), 0);
dataXLS(:, 7) = num2cell(round(areaLandmarks(:, 6) * 100 ) / 100.0 );   % score based on ratio

xlwrite(filenameXLS, [headXLS; dataXLS]);

