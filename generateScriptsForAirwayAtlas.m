function generateScriptsForAirwayAtlas( scriptsFile, dataFoldName, resultsDirName, sampleNum )

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
%   scriptsFile: the name of the file, with chmod +x , you can run it at command line
%   dataFoldName: the folder name of the prepared data, like ./data
%   resultsDirName: the folder name of the generated results, like ./results
%   sampleNum: the number of the points on the centerline, usually 100
% 
%   the name of the following prepared data can be changed:
%   binary segmentation: *_Segmented_NoSinuses.nrrd
%   Landmarks: *_Points.txt
%   The sphere used to remove the mouth: *_CuttingPlanes.txt
%   airway geometry: *_Airway_NoSinuses_NoMouth_Smooth.vtk
% 
%   Author: Yi Hong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%{
cases = {'CRL02', 'CRL04', 'CRL06', 'CRL16', 'CRL21', 'CRL23', 'CRL27', 'CRL31', 'CRL32', 'CRL34', 'CRL36', 'CRL39', 'CRL40', 'CRL42', 'CRL44', 'CRL46', 'CRL48', 'CRL49', 'CRL50', 'CRL51', ...
		 'CRL53', 'CRL54', 'CRL59', 'CRL60', 'CRL61', 'CRL62', 'CRL63', 'CRL64', 'CRL65', 'CRL66', 'CRL67', 'CRL68', 'CRL71', 'CRL72', 'CRL73', 'CRL74', 'CRL76', 'CRL77', 'CRL78', 'CRL80', ...
		 'CRL81', 'CRL82', 'CRL83', 'CRL84', 'CRL85', 'CRL86', 'CRL88', 'CRL89', 'CRL90', 'CRL91', 'CRL92', 'CRL93', 'CRL94', 'CRL95', 'CRL96', 'CRL97', 'CRL98', 'CRL99', 'CRL102', 'CRL104', ...
		 'CRL105', 'CRL107', 'CRL108', 'CRL109', 'CRL111', 'CRL113', 'CRL114', 'CRL115', 'CRL116', 'CRL117', 'CRL118', 'CRL119', 'CRL120', 'CRL121', 'CRL122', 'CRL123', 'CRL124', 'CRL127', ...
		 'SGS03ManPre', 'SGS07', 'SGS11', 'SGS12', 'SGS13', 'SGS18', 'SGS03ManPost', 'SGS05', 'SGS06', 'SGS08', 'SGS09', 'SGS10', 'SGS14', 'SGS17', 'SGS07_V3', 'SGS04_V3', 'SGS01_V3'};
%}
%cases = {'2024'};
cases = { '1002', '1003', '1004', '1005' }

fid = fopen( scriptsFile, 'wt' );
for iI = 1:max(size( cases))
    fprintf( fid, '######  %s  ######\n', cases{iI} );
    
    % scripts for computing centerline and its tangent normals    
    fprintf( fid, 'rm -r ./%s/IsoSurface%s\n', resultsDirName, cases{iI} );
    fprintf( fid, 'mkdir ./%s/IsoSurface%s\n', resultsDirName, cases{iI} );
    fprintf( fid, 'matlab -nodesktop -nosplash -nodisplay -r "cd ./code/centerline; generateCenterlineOfAirway(' );
    fprintf( fid, '''%s'', ''../../%s/%s_OUTPUT.nrrd'', ''../../%s/IsoSurface%s'', %d, ''../../%s/%s_LANDMARKS.txt'', ''../../%s/%s_CLIPPINGS.txt'' ); ', ...
                   cases{iI}, dataFoldName, cases{iI}, resultsDirName, cases{iI}, sampleNum, dataFoldName, cases{iI}, dataFoldName, cases{iI});
    fprintf( fid, 'cd ../..; quit;"\n\n' );
    
    % scripts for computing the cross-sectional area based on above centerline
    fprintf( fid, 'rm -r ./%s/Contour%s\n', resultsDirName, cases{iI} );
    fprintf( fid, 'mkdir ./%s/Contour%s\n', resultsDirName, cases{iI} );
    fprintf( fid, './code/crossSections/bin/computeAreaAndContourWithMeanNorm ./%s/%s_OUTPUT.vtk ./%s/IsoSurface%s/%s_MeanAndNormal.txt ./%s/%s_CuttingPlanes.txt ', ...
                   dataFoldName, cases{iI}, resultsDirName, cases{iI}, cases{iI}, dataFoldName, cases{iI} );
    fprintf( fid, './%s/Contour%s/%s_Area.txt ./%s/Contour%s/%s_Perimeter.txt ./%s/Contour%s/contour\n\n', resultsDirName, cases{iI}, cases{iI}, resultsDirName, cases{iI}, cases{iI}, resultsDirName, cases{iI} );
  
    % scripts for displaying the cross-sections and the cross-sectional curve
    fprintf( fid, 'rm -r ./%s/Curves%s\n', resultsDirName, cases{iI} );
    fprintf( fid, 'mkdir ./%s/Curves%s\n', resultsDirName, cases{iI} );
    
    fprintf( fid, 'mkdir ./%s/Curves%s/contour%s_contours\n', resultsDirName, cases{iI}, cases{iI} );

    fprintf( fid, 'matlab -nodesktop -nosplash -nodisplay -r "cd ./code/display; drawContourOfAirway(' );
    fprintf( fid, '''../../%s/IsoSurface%s/%s_MeanAndNormal.txt'', ''../../%s/IsoSurface%s/%s_NormReliableCheck.txt'', ', ...
                   resultsDirName, cases{iI}, cases{iI}, resultsDirName, cases{iI}, cases{iI});
    fprintf( fid, '''../../%s/Contour%s/contour'', ''../../%s/Curves%s/contour%s'', ''../../%s/Curves%s/%s_AreaEllipse.txt'', ''../../%s/%s_LANDMARKS.txt''); ', ...
                   resultsDirName, cases{iI}, resultsDirName, cases{iI}, cases{iI}, resultsDirName, cases{iI}, cases{iI}, dataFoldName, cases{iI} );
    fprintf( fid, 'cd ../..; quit;"\n\n' );
    
    fprintf( fid, 'matlab -nodesktop -nosplash -nodisplay -r "cd ./code/display; drawLandmarksOnOneCurve(' );
    fprintf( fid, '''../../%s/IsoSurface%s/%s_MeanAndNormal.txt'', ''../../%s/Contour%s/%s_Area.txt'', ''../../%s/Contour%s/%s_Perimeter.txt'', ''../../%s/Curves%s/%s_AreaEllipse.txt'', ', ...
                   resultsDirName, cases{iI}, cases{iI}, resultsDirName, cases{iI}, cases{iI}, resultsDirName, cases{iI}, cases{iI}, resultsDirName, cases{iI}, cases{iI} );
    fprintf( fid, '''../../%s/IsoSurface%s/%s_LandmarksIdOnCenterline.txt'', ''../../%s/Curves%s/%s_CurveAndLandmarks.png'', ''../../%s/Curves%s''); ', ...
                   resultsDirName, cases{iI}, cases{iI}, resultsDirName, cases{iI}, cases{iI}, resultsDirName, cases{iI} );
    fprintf( fid, 'cd ../..; quit;"\n\n' );
    
end
fclose( fid );
