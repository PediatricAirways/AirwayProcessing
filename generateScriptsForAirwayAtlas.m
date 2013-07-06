function generateScriptsForAirwayAtlas( scriptsFile, dataFoldName, resultsDirName, sampleNum )

%
% scriptsFile: the name of the file, with chmod +x , you can run it at command line
% dataFoldName: the folder name of the prepared data, like ./data
% resultsDirName: the folder name of the generated results, like ./results
% sampleNum: the number of the points on the centerline, usually 100
% 
% the name of the following prepared data can be changed:
% binary segmentation: *_Segmented_NoSinuses.nrrd
% Landmarks: *_Points.txt
% The sphere used to remove the mouth: *_CuttingPlanes.txt
% airway geometry: *_Airway_NoSinuses_NoMouth_Smooth.vtk
% 
 

cases = {'CRL02', 'CRL16'};
fid = fopen( scriptsFile, 'wt' );
for iI = 1:max(size( cases))
    fprintf( fid, '######  %s  ######\n', cases{iI} );
    
    % scripts for computing centerline and its tangent normals    
    fprintf( fid, 'rm -r ./%s/IsoSurface%s\n', resultsDirName, cases{iI} );
    fprintf( fid, 'mkdir ./%s/IsoSurface%s\n', resultsDirName, cases{iI} );
    fprintf( fid, 'matlab -nodesktop -nosplash -nodisplay -r "cd ./code/centerline; generateCenterlineOfAirway(' );
    fprintf( fid, '''%s'', ''../../%s/%s_Segmented_NoSinuses.nrrd'', ''../../%s/IsoSurface%s'', %d, ''../../%s/%s_Points.txt'', ''../../%s/%s_CuttingPlanes.txt'' ); ', ...
                   cases{iI}, dataFoldName, cases{iI}, resultsDirName, cases{iI}, sampleNum, dataFoldName, cases{iI}, dataFoldName, cases{iI});
    fprintf( fid, 'cd ../..; quit;"\n\n' );
    
    % scripts for computing the cross-sectional area based on above centerline
    fprintf( fid, 'rm -r ./%s/Contour%s\n', resultsDirName, cases{iI} );
    fprintf( fid, 'mkdir ./%s/Contour%s\n', resultsDirName, cases{iI} );
    fprintf( fid, './code/crossSections/bin/computeAreaAndContourWithMeanNorm ./%s/%s_Airway_NoSinuses_NoMouth_Smooth.vtk ./%s/IsoSurface%s/%s_MeanAndNormal.txt ./%s/%s_CuttingPlanes.txt ', ...
                   dataFoldName, cases{iI}, resultsDirName, cases{iI}, cases{iI}, dataFoldName, cases{iI} );
    fprintf( fid, './%s/Contour%s/%s_Area.txt ./%s/Contour%s/contour\n\n', resultsDirName, cases{iI}, cases{iI}, resultsDirName, cases{iI} );
  
    % scripts for displaying the cross-sections and the cross-sectional curve
    fprintf( fid, 'rm -r ./%s/Curves%s\n', resultsDirName, cases{iI} );
    fprintf( fid, 'mkdir ./%s/Curves%s\n', resultsDirName, cases{iI} );
    
    fprintf( fid, 'mkdir ./%s/Curves%s/contour%s_contours\n', resultsDirName, cases{iI}, cases{iI} );

    fprintf( fid, 'matlab -nodesktop -nosplash -nodisplay -r "cd ./code/display; drawContourOfAirway(' );
    fprintf( fid, '''../../%s/IsoSurface%s/%s_MeanAndNormal.txt'', ''../../%s/IsoSurface%s/%s_NormReliableCheck.txt'', ', ...
                   resultsDirName, cases{iI}, cases{iI}, resultsDirName, cases{iI}, cases{iI});
    fprintf( fid, '''../../%s/Contour%s/contour'', ''../../%s/Curves%s/contour%s'', ''../../%s/Curves%s/%s_AreaEllipse.txt'', ''../../%s/%s_Points.txt''); ', ...
                   resultsDirName, cases{iI}, resultsDirName, cases{iI}, cases{iI}, resultsDirName, cases{iI}, cases{iI}, dataFoldName, cases{iI} );
    fprintf( fid, 'cd ../..; quit;"\n\n' );
    
    fprintf( fid, 'matlab -nodesktop -nosplash -nodisplay -r "cd ./code/display; drawLandmarksOnOneCurve(' );
    fprintf( fid, '''../../%s/IsoSurface%s/%s_MeanAndNormal.txt'', ''../../%s/Contour%s/%s_Area.txt'', ''../../%s/Curves%s/%s_AreaEllipse.txt'', ', ...
                   resultsDirName, cases{iI}, cases{iI}, resultsDirName, cases{iI}, cases{iI}, resultsDirName, cases{iI}, cases{iI} );
    fprintf( fid, '''../../%s/IsoSurface%s/%s_LandmarksIdOnCenterline.txt'', ''../../%s/Curves%s/%s_CurveAndLandmarks.png'', ''../../%s/Curves%s''); ', ...
                   resultsDirName, cases{iI}, cases{iI}, resultsDirName, cases{iI}, cases{iI}, resultsDirName, cases{iI} );
    fprintf( fid, 'cd ../..; quit;"\n\n' );
    
end
fclose( fid );
