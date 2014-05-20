function generateCenterlineOfAirway( subjectName, inputAirwaySegmentation, outputIsosurfacePrefix, nSurface, landmarksFile, clipFile )

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
%   given the airway segmentation, automatically generate the centerline
%   
%   subjectName: the name of the subject
%   inputAirwaySegmentation: the segmentation volume of airway
%   outputIsosurfacePrefix: the output directory
%   nSurface: the number of points on the centerline
%   landmarksFile: the input file with landmarks
%   clipFile: the file with sphere to remove the mouth and other unwanted parts
%
%   Author: Yi Hong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read the airway segmentation results
[im_data, im_info] = nrrdRead(inputAirwaySegmentation);
sizes = sscanf(im_info.sizes, '%d');

dims = length(sizes);
if dims == 2
    space_direction = sscanf(im_info.spacedirections, '(%f,%f) (%f,%f)');
    space_origin = sscanf(im_info.spaceorigin, '(%f,%f)');
elseif dims == 3
    space_direction = sscanf(im_info.spacedirections, '(%f,%f,%f) (%f,%f,%f) (%f,%f,%f)');
    space_origin = sscanf(im_info.spaceorigin, '(%f,%f,%f)');
elseif dims == 4
    space_direction = sscanf(im_info.spacedirections, '(%f,%f,%f,%f) (%f,%f,%f,%f) (%f,%f,%f,%f) (%f,%f,%f,%f)');
    space_origin = sscanf(im_info.spaceorigin, '(%f,%f,%f,%f)');
end
space_direction = reshape(space_direction, dims, dims);

% compute the boundray box for the airway
bounding_box = zeros( dims, 1 );
for i = 1:2
    space_origin(i) = -space_origin(i);
    space_direction(i, i) = -space_direction(i, i);
end
bounding_box(1:dims) = space_origin;
bounding_box(dims+1:2*dims) = space_direction * (sizes - ones(dims,1)) + bounding_box(1:dims);

spacing = zeros(dims,1);
for i = 1:dims
    spacing(i) = abs(space_direction(i,i));
end

% landmarks are from slicer3d, so it's in RAS coordinate system
% fidLandmarks = fopen( landmarksFile );
% tline = fgets( fidLandmarks );
% nLandmarks = sscanf( tline, '%d' );
% for i = 1:nLandmarks
%     tline = fgets( fidLandmarks );
%     landmarks( i, : ) = sscanf( tline, '%f' );
%     % change RAS to world system
%     tmp = landmarks( i, 1 );
%     landmarks( i, 1 ) = landmarks( i, 2 );
%     landmarks( i, 2 ) = tmp;
% end
% fclose( fidLandmarks );

nLandmarks = 6;
[key, val] = textread( landmarksFile, '%s:%[^\n]' );
landmarks = zeros( nLandmarks, 3 );
flagLM = zeros( nLandmarks, 1 );
nosePlane = zeros( 4, 3 );
for iI = 1:length(key)
    landmarksTmp = str2num( val{iI} );
    landmarksIdTmp = 0;
    disp( strcat('"', key{iI}, '"' ));
    if strcmp( key{iI}, 'NasalSpine' )
	landmarksIdTmp = 1;
	pntNose = landmarksTmp';
	idNoseRecord = 1;
    elseif strcmp( key{iI}, 'PosteriorInferiorVomerCorner' )
	pntVomerCorner = landmarksTmp';
	landmarksIdTmp = 2;
    elseif strcmp( key{iI}, 'EpiglottisTip' )
	landmarksIdTmp = 3;
    elseif strcmp( key{iI}, 'TVC' )
	landmarksIdTmp = 4;
    elseif strcmp( key{iI}, 'Subglottic' )
	landmarksIdTmp = 5;
    elseif strcmp( key{iI}, 'TracheaCarina' )
	landmarksIdTmp = 6;
	pntBranch = landmarksTmp';
    elseif strcmp( key{iI}, 'NoseTip' )
	nosePlane( 1, : ) = landmarksTmp;
    elseif strcmp( key{iI}, 'Columella' )
    	nosePlane( 2, : ) = landmarksTmp;
    elseif strcmp( key{iI}, 'LeftAlaRim' )
	nosePlane( 3, : ) = landmarksTmp;
    elseif strcmp( key{iI}, 'RightAlaRim' )
	nosePlane( 4, : ) = landmarksTmp;
    end
    if landmarksIdTmp >= 1 && landmarksIdTmp <= nLandmarks
    	landmarks( landmarksIdTmp, : ) = landmarksTmp;
    	% change RAS to world system
    	tmp = landmarks( landmarksIdTmp, 1 );
    	landmarks( landmarksIdTmp, 1 ) = landmarks( landmarksIdTmp, 2 );
    	landmarks( landmarksIdTmp, 2 ) = tmp;
	flagLM( landmarksIdTmp ) = 1;
    end
end

for iI = 1:length( flagLM )
	if flagLM( iI ) == 0 
		landmarks( iI, : ) = [];
	end
end
nLandmarks = size( landmarks, 1 );

% output the landmarks we used for registration
landmarks

% output the landmarks on the nose's plane
nosePlane

% try to calculate the start/end planes 
[coeff] = compPrincipalDirection( landmarks );

% take the axis with the largest eigenvalue as the normal of the
% cutting plane for the branches, and the second for the noses
[value_1, axisIdBranch] = max( abs( coeff( :, 1 ) ) );
[value_2, axisIdNose] = max( abs( coeff( :, 2 ) ) );
%axisIdBranch
%axisIdNose
if abs( value_1 - 0.7 ) < 0.5
    if idNoseRecord == 1
        [coeff_1] = compPrincipalDirection( landmarks( 2:end, : ) );
    elseif idNoseRecord == size( key, 1 )
        [coeff_1] = compPrincipalDirection( landmarks( 1:end-1, : ) );
    else
        [coeff_1] = compPrincipalDirection( [ landmarks( 1:idNoseRecord-1, :); landmarks( idNoseRecord+1:end, : ) ] );
    end 
    [value_1, axisIdBranch] = max( abs( coeff_1( :, 1 ) ) );
    coeff( axisIdBranch, 2 ) = 0;
    [value_2, axisIdNose] = max( abs( coeff( :, 2 ) ) ); 
end

axisId_1 = min( axisIdNose, axisIdBranch );
axisId_2 = max( axisIdNose, axisIdBranch );
%im_slice = zeros( sizes( axisId_1 ), sizes( axisId_2 ) );
im_slice = zeros( size( im_data, axisId_1 ), size( im_data, axisId_2 ) );
for iI = 1:size(im_data, axisId_1)
    for iJ = 1:size(im_data, axisId_2)
        if 1 ~= axisId_1 && 1 ~= axisId_2
            im_slice( iI, iJ ) = sum( double( im_data( :, iI, iJ ) ) );
        elseif 2 ~= axisId_1 && 2 ~= axisId_2
            im_slice( iI, iJ ) = sum( double( im_data( iI, :, iJ ) ) );
        elseif 3 ~= axisId_1 && 3 ~= axisId_2
            im_slice( iI, iJ ) = sum( double( im_data( iI, iJ, : ) ) );
        end
    end
end
im_slice = im_slice ./ max( max( im_slice ) );
figure(1), h = imagesc( (im_slice), [0 1] );
filename = sprintf('%s/%s_slice.png', outputIsosurfacePrefix, subjectName);
h
filename
saveas(h, filename);

ijkNose = round( inv( space_direction ) * ( pntNose - space_origin ) );
ijkBranch = round( inv( space_direction ) * ( pntBranch - space_origin ) );
ijkNose = max( ijkNose, [1; 1; 1] );
ijkNose = min( ijkNose, sizes );
ijkBranch = max( ijkBranch, [1; 1; 1] );
ijkBranch = min( ijkBranch, sizes );

ijkNose
ijkBranch

tmp = ijkNose(1); ijkNose(1) = ijkNose(2); ijkNose(2) = tmp;
tmp = ijkBranch(1); ijkBranch(1) = ijkBranch(2); ijkBranch(2) = tmp;

bLow = sign( ijkNose( axisIdNose ) - ijkBranch( axisIdNose) ) * axisIdNose;
vLow = ijkNose( axisIdNose );
bHigh = sign( ijkBranch( axisIdBranch ) - ijkNose( axisIdBranch ) ) * axisIdBranch;
vHigh = ijkBranch( axisIdBranch );

bLow 
vLow
bHigh
vHigh

% clip the airway's mouth
%clipFile
[key, val] = textread( clipFile, '%s:%[^\n]' );
nNumSphere = 0;
mouthCenter = [];
mouthRadius = [];
for iI = 1:length(key)
    if strcmp( key{iI}, 'ClipSphereCenter' )
	nNumSphere = nNumSphere + 1;
        mouthCenter( nNumSphere, : ) = ( str2num( val{iI} ) )';
    elseif strcmp( key{iI}, 'ClipSphereRadius' )
        mouthRadius( nNumSphere ) = str2num( val{iI} );
    end
end

% the airway's starting plane
%startPlaneFile
%[key, val] = textread( startPlaneFile, '%s:%[^\n]' );
%for iI = 1:size( key, 1 )
%    if strcmp( key{iI}, 'NasalPlaneCenter' )
%        nasalCenter = ( str2num( val{iI} ) )';
%    elseif strcmp( key{iI}, 'NasalPlaneNormal' )
%        nasalNormal = str2num( val{iI} );
%    end
%end
%nasalCenter
%nasalNormal

% calcuate the start plane using the landmarks
nasalCenter = pntNose;
nasalNormal = cross(nosePlane(1,:) - nosePlane(2,:), nosePlane(4,:) - nosePlane(3,:));
nasalNormal = nasalNormal ./ (norm(nasalNormal));

% write out the start/end planes
start_end_file = sprintf( '%s/%s_Planes.txt', outputIsosurfacePrefix, subjectName ); 
start_end_fileId = fopen( start_end_file, 'wt' );
fprintf( start_end_fileId, 'NasalPlaneCenter : %f %f %f\n', nasalCenter(1), nasalCenter(2), nasalCenter(3) ); 
fprintf( start_end_fileId, 'NasalPlaneNormal : %f %f %f\n', nasalNormal(1), nasalNormal(2), nasalNormal(3) );
fprintf( start_end_fileId, 'CarinaPlaneCenter : %f %f %f\n', pntBranch(1), pntBranch(2), pntBranch(3) ); 
fprintf( start_end_fileId, 'CarinaPlaneNormal : %f %f %f\n', sign(bHigh)*(abs(bHigh) == 1), sign(bHigh)*(abs(bHigh) == 2), sign(bHigh)*(abs(bHigh) == 3) );
fclose(start_end_fileId);

% remove the mouth and label nose
im_data = clipAirwayVolume( im_data, space_direction, bounding_box, mouthCenter, mouthRadius, nasalCenter, nasalNormal ); 

%imLabel = labelMap(im_data, bLow, vLow, bHigh, vHigh);
imLabel = labelMapForClipVolume( im_data, bHigh, vHigh );
%writeMETA( imLabel, 'testLabel.mhd', 'MET_FLOAT', space_origin, spacing );
airwayCrossSections( subjectName, imLabel, spacing, bounding_box, outputIsosurfacePrefix, nSurface, landmarks );


% function of removing the mouth and label the nose
function im_clip = clipAirwayVolume( im_data, space_direction, bounding_box, mouthCenter, mouthRadius, nasalCenter, nasalNormal )

% get rid of the mouth using a sphere
% a cutting plane is started from noses

dx = abs( space_direction( 1, 1 ) );
dy = abs( space_direction( 2, 2 ) );
dz = abs( space_direction( 3, 3 ) );

if bounding_box(4) < bounding_box(1)
   dx = -dx;
end
if bounding_box(5) < bounding_box(2)
   dy = -dy;
end
if bounding_box(6) < bounding_box(3)
   dz = -dz;
end

sz = size( im_data );
im_clip = im_data;
d_plane = -dot( nasalNormal, nasalCenter );
pos = zeros( 3, 1 );

for iI = 1:sz(1)
    for iJ = 1:sz(2)
        for iK = 1:sz(3)
            if im_clip( iI, iJ, iK ) == 0
            	continue;
            end
            pos(1) = bounding_box(1) + dy * iJ;
            pos(2) = bounding_box(2) + dx * iI;
            pos(3) = bounding_box(3) + dz * iK;
            if dot( pos, nasalNormal ) + d_plane  >= 0 && norm( pos - nasalCenter ) <= 40
            	im_clip( iI, iJ, iK ) = 4;   % label the nose with 4
            end
	    % remove the mouth
	    for iP = 1:length(mouthRadius)
            	if norm( pos - mouthCenter(iP,:)' ) <= mouthRadius(iP)
            	    im_clip( iI, iJ, iK ) = 0;
		    break;
            	end
	    end
        end
    end
end


function im_labelmap = labelMap(im_data, bLow, vLow, bHigh, vHigh)
% bLow: show which plane for setting dirichlet low among x, y, z plane
%       x left plane: -1, x right plane: 1, 
%       y left plane: -2, y right plane: 2,
%       z left plane: -3, z right plane: 3
% vLow: the value for the plane shown by bLow
% bHigh: show which plane for setting dirichlet high among x, y, z plane
% vHigh: the value for the plane shown by bHigh
% Solution domain: 11
% Dirichlet low: 4
% Dirichlet high: 5
% Neumann: 6

dimsize = size(im_data);
im_labelmap = zeros(dimsize(1), dimsize(2), dimsize(3));
for i=1:dimsize(1)
    for j=1:dimsize(2)
        for k=1:dimsize(3) 
            if im_data(i,j,k) > 0
                if bLow == -1 && i <= vLow || bLow == 1 && i >= vLow || ...
                   bLow == -2 && j <= vLow || bLow == 2 && j >= vLow || ...
                   bLow == -3 && k <= vLow || bLow == 3 && k >= vLow
                    im_labelmap(i,j,k) = 4;
                elseif bHigh == -1 && i <= vHigh || bHigh == 1 && i >= vHigh || ...
                   bHigh == -2 && j <= vHigh || bHigh == 2 && j >= vHigh || ...
                   bHigh == -3 && k <= vHigh || bHigh == 3 && k >= vHigh
                    im_labelmap(i,j,k) = 5;
                else
                    im_labelmap(i,j,k) = 11;
                end
            elseif i-1 >= 1 && im_data(i-1,j,k) > 0 || ...
                   i+1 <= dimsize(1) && im_data(i+1,j,k) > 0 || ...
                   j-1 >= 1 && im_data(i,j-1,k) > 0 || ...
                   j+1 <= dimsize(2) && im_data(i,j+1,k) > 0 || ...
                   k-1 >= 1 && im_data(i,j,k-1) > 0 || ...
                   k+1 <= dimsize(3) && im_data(i,j,k+1) > 0 
               im_labelmap(i,j,k) = 6;
           end
        end
    end
end


function im_labelmap = labelMapForClipVolume(im_data, bHigh, vHigh)
% bLow: show which plane for setting dirichlet low among x, y, z plane
%       x left plane: -1, x right plane: 1,
%       y left plane: -2, y right plane: 2,
%       z left plane: -3, z right plane: 3
% vLow: the value for the plane shown by bLow
% bHigh: show which plane for setting dirichlet high among x, y, z plane
% vHigh: the value for the plane shown by bHigh
% Solution domain: 11
% Dirichlet low: 4
% Dirichlet high: 5
% Neumann: 6

dimsize = size(im_data);
im_labelmap = zeros(dimsize(1), dimsize(2), dimsize(3));
for i=1:dimsize(1)
    for j=1:dimsize(2)
        for k=1:dimsize(3)
   	    if im_data(i,j,k) == 4
		im_labelmap(i,j,k) = 4;
		continue;
	    end
            if im_data(i,j,k) > 0
                if bHigh == -1 && i <= vHigh || bHigh == 1 && i >= vHigh || ...
                   bHigh == -2 && j <= vHigh || bHigh == 2 && j >= vHigh || ...
                   bHigh == -3 && k <= vHigh || bHigh == 3 && k >= vHigh
                    im_labelmap(i,j,k) = 5;
                else
                    im_labelmap(i,j,k) = 11;
                end
            elseif i-1 >= 1 && im_data(i-1,j,k) > 0 || ...
                   i+1 <= dimsize(1) && im_data(i+1,j,k) > 0 || ...
                   j-1 >= 1 && im_data(i,j-1,k) > 0 || ...
                   j+1 <= dimsize(2) && im_data(i,j+1,k) > 0 || ...
                   k-1 >= 1 && im_data(i,j,k-1) > 0 || ...
                   k+1 <= dimsize(3) && im_data(i,j,k+1) > 0
               im_labelmap(i,j,k) = 6;
           end
        end
    end   
end


function airwayCrossSections( subjectName, imLabel, spacing, boundbox, outputIsosurfacePrefix, nSurface, landmarks )

% set spacings in x,y, and z direction
dx = spacing(1);
dy = spacing(2);
dz = spacing(3);

solIm = solveLaplaceLinearSystem( imLabel, dx, dy, dz );


filename = sprintf( '%s/testLabel.mhd', outputIsosurfacePrefix );
writeMETA( imLabel, filename, 'MET_FLOAT', [0, 0, 0], spacing );
filename = sprintf( '%s/testSolIm.mhd', outputIsosurfacePrefix );
writeMETA( solIm, filename, 'MET_FLOAT', [0, 0, 0], spacing );


binVals = linspace( 0, 1, 8*nSurface );
N = hist( solIm(:), binVals );
%N = sqrt(N ./ min(N));
N = sqrt(N);
N = N / sum(N);
cdfN = cumsum( N );
cdfN(length(N)) = 1;

isovals = linspace(0,1,nSurface);
isovalsNormalized = zeros( size( isovals ) );

for iI=1:length( isovals )
  indx = find( cdfN >= isovals(iI) );
  isovalsNormalized( iI ) = binVals( indx(1) );
  if isovalsNormalized( iI ) == 1
      isovalsNormalized( iI ) = isovalsNormalized( iI ) - 0.000000001;
  end
end

iBegin = 1;
while( iBegin <= length(isovals) )
  if( isovalsNormalized( iBegin ) ~= 0 ) break; end
  iBegin = iBegin + 1;
end
iBegin = max( iBegin-1, 1 );

iEnd = length(isovals);
while( iEnd >= 1 )
  if( 1 - isovalsNormalized( iEnd ) > 0.00000001 ) break; end
  iEnd = iEnd - 1;
end
iEnd = min( iEnd+1, length(isovals) );

isovalsNormalized( iBegin:iEnd )

visualizeCrossSections( subjectName, solIm, isovalsNormalized(iBegin:iEnd), dx, dy, dz, boundbox, outputIsosurfacePrefix, landmarks );
%set(gcf,'color',[1 1 1] );
%figName = sprintf( '%s/surface.fig', outputIsosurfacePrefix );
%saveas(h, figName);



function solIm = solveLaplaceLinearSystem( labelIm, dx, dy, dz )
% function solIm = solveLaplaceLinearSystem( labelIm, dx, dy, dz )
% Returns the solution of Laplace's equation given the labelIm and
% the spacings in x, y, and z direction

% input image is in labelIm

sz = size( labelIm );

% Do Laplace solution by iterative solver using 6 neighborhood

% Solution domain: 11
% Dirichlet low: 4
% Dirichlet high: 5
% Neumann: 6

solId = 11;
neumannId = 6;
dirLowId = 4;
dirHighId = 5;

% first get all the solution domain indices

indxS = find( labelIm == solId );

nrOfUnknowns = length( indxS );

% create a matrix which holds the indices to the unknowns
solInd = zeros( size( labelIm ) );
solInd( indxS ) = 1:nrOfUnknowns;

[I1,I2,I3] = ind2sub( sz, indxS );

indxP1 = sub2ind( sz, min( sz(1), I1+1 ), I2, I3 );
indxM1 = sub2ind( sz, max( 1, I1-1 ), I2, I3 );

indxP2 = sub2ind( sz, I1, min( sz(2), I2+1 ), I3 );
indxM2 = sub2ind( sz, I1, max( 1, I2-1 ), I3 );

indxP3 = sub2ind( sz, I1, I2, min( sz(3), I3+1 ) );
indxM3 = sub2ind( sz, I1, I2, max( 1, I3-1 ) );

% set up the linear system by going through all the indices of the solution domain and creating the sparse matrix A and the right hand side b

A = -2*(1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz ) )*speye(nrOfUnknowns, nrOfUnknowns);
b = sparse( nrOfUnknowns, 1 );

[A,b] = updateAB( A, b, indxP1, labelIm, solInd, dx );
[A,b] = updateAB( A, b, indxM1, labelIm, solInd, dx );
[A,b] = updateAB( A, b, indxP2, labelIm, solInd, dy );
[A,b] = updateAB( A, b, indxM2, labelIm, solInd, dy );
[A,b] = updateAB( A, b, indxP3, labelIm, solInd, dz );
[A,b] = updateAB( A, b, indxM3, labelIm, solInd, dz );

% now solve it

sol = A\b;

% create the solution image

solIm = NaN*zeros( size( labelIm ) );
solIm( indxS ) = sol;
solIm( labelIm==dirLowId ) = 0;
solIm( labelIm==dirHighId ) = 1;


function [A,b] = updateAB( Ain, bin, indx, labelIm, solInd, dx )

nrOfUnknowns = length( bin );

fprintf('Executing updateAB ...');

b = bin;

solId = 11;
neumannId = 6;
dirLowId = 4;
dirHighId = 5;

lowVal = 0;
highVal = 1;

o_dxs = 1/(dx*dx);

iA = zeros( nrOfUnknowns, 1 );
jA = zeros( nrOfUnknowns, 1 );
vA = zeros( nrOfUnknowns, 1 );
setInd = logical( zeros( nrOfUnknowns, 1 ) );

for iI=1:nrOfUnknowns
  
  %A( iI, iI ) = A( iI, iI ) - 2*o_dxs;

  curLabel =  labelIm( indx( iI ) );
  if ( curLabel == solId )
    curId = solInd( indx( iI ) );
    iA( iI ) = iI;
    jA( iI ) = curId;
    vA( iI ) = o_dxs;
    setInd( iI ) = true;
    %A( iI, curId ) = A( iI, curId ) + o_dxs;
  elseif ( curLabel == neumannId )
    iA( iI ) = iI;
    jA( iI ) = iI;
    vA( iI ) = o_dxs;
    setInd( iI ) = true;
    %A( iI, iI ) = A( iI, iI ) + o_dxs;
  elseif ( curLabel == dirLowId )
    b( iI ) = b( iI ) - lowVal*o_dxs;
  elseif ( curLabel == dirHighId )
    b( iI ) = b( iI ) - highVal*o_dxs;
  else
    disp( 'Unknown id');
  end

  if ( mod( iI, round( nrOfUnknowns/100 ) ) == 0 )
    fprintf('#');
  end
  
end

% only extract the ones that were actually set

iA = iA( setInd );
jA = jA( setInd );
vA = vA( setInd );

A = Ain + sparse( iA, jA, vA, nrOfUnknowns, nrOfUnknowns );

fprintf('done.\n');


function visualizeCrossSections( subjectName, solIm, isovals, dx, dy, dz, boundbox, outputIsosurfacePrefix, landmarks)

sz = size( solIm );

if boundbox(4) < boundbox(1)
   dx = -dx;
end
if boundbox(5) < boundbox(2)
   dy = -dy; 
end
if boundbox(6) < boundbox(3)
   dz = -dz;
end

solImTmp = solIm;
clear solIm;
for iK = 1:sz(3)
    solIm( :, :, iK ) = solImTmp( :, :, iK )';
end
clear solImTmp;
sz = size( solIm );
tmp = dx; dx = dy; dy = tmp;

[x,y,z] = meshgrid( dy*(0:sz(2)-1)+boundbox(2), dx*(0:sz(1)-1)+boundbox(1), ...
                    dz*(0:sz(3)-1)+boundbox(3) );

landmarksIsoval = zeros( size( landmarks, 1 ), 1 );
% find isovalues for landmarks
for ix = 1:size( landmarks, 1 )
    grid_I = round( ( landmarks(ix,2) - boundbox(1) ) / dx );
    grid_J = round( ( landmarks(ix,1) - boundbox(2) ) / dy );
    grid_K = round( ( landmarks(ix,3) - boundbox(3) ) / dz );
    grid_I = min( max( grid_I, 1 ), sz(1) );
    grid_J = min( max( grid_J, 1 ), sz(2) );
    grid_K = min( max( grid_K, 1 ), sz(3) );
    landmarksIsoval( ix ) = solIm( grid_I, grid_J, grid_K ); 
end

startIsoval = landmarksIsoval(3);
endIsoval = landmarksIsoval(4);
if isnan( startIsoval ) || isnan( endIsoval )
    startIsoval = 0.5;
    endIsoval = 0.5;
end
  
referDir_section = zeros( 3, 1 );
referDir_section(1) = landmarks(3, 2) - landmarks(4, 2);
referDir_section(2) = landmarks(3, 1) - landmarks(4, 1);
referDir_section(3) = landmarks(3, 3) - landmarks(4, 3);
referDir_section = referDir_section ./ ( norm(referDir_section) + 1e-15 );

% store the closest points and isoval  on centerline for landmarks
closestPointsLM = zeros( size( landmarks, 1 ), 2 );

meanPnt = zeros( length(isovals), 3 );
normPnt = zeros( length(isovals), 3 );
normCheck = zeros( length(isovals), 1 );
meanIsovals = zeros( length( isovals ), 1 );

cmap = jet(length(isovals)+1);

handle_surface = figure(2);
handle_contour = figure(3);


for iI=1:length( isovals )

  fprintf('Rendering patch %d/%d\n', iI, length( isovals ) );
  
  handle_surface = figure(2); [f, v] = isosurface( x, y, z, solIm, isovals( iI ) );
  hold on
  if sum(size(v)) == 0
    continue
  end
  tmp = v(:,1); v(:,1) = v(:,2); v(:,2) = tmp;
  figure(2), p = patch( 'Faces', f, 'Vertices', v );
  %set(p, 'FaceColor', 'blue', 'EdgeColor', 'none');  
  set(p, 'FaceColor', cmap(round(isovals(iI) * length(isovals))+1, :), 'EdgeColor', 'none'); 
  
  filename = sprintf('%s/fv%03d.txt', outputIsosurfacePrefix, iI); 
  fid = fopen(char(filename), 'wt');
  szv = size(v);
  szf = size(f);
  fprintf(fid, '%d %d\n', szv(1), szf(1));
  for ix = 1:szv(1)
      fprintf(fid, '%f %f %f\n', v(ix,1), v(ix,2), v(ix,3));
  end
  for ix = 1:szf(1)
      fprintf(fid, '%d %d %d\n', f(ix,1)-1, f(ix,2)-1, f(ix,3)-1);
  end
  fclose(fid);

  % find the edge of isosurface
  vertexOnEdge = findEdgeOfIsosurface( v, szf(1), f );
  handle_contour = figure(3); plot3( vertexOnEdge(:,1), vertexOnEdge(:,2), vertexOnEdge(:,3), 'b.', 'MarkerSize', 0.5 );
  hold on

  % compute the mean of the contour
  planeMatrix = vertexOnEdge;
  for ix = 1:3
      meanPnt(iI, ix) = mean( vertexOnEdge(:, ix ) );
      planeMatrix(:, ix) = planeMatrix(:,ix) - meanPnt(iI, ix);
  end	
  figure(3), plot3( meanPnt(iI, 1), meanPnt(iI, 2), meanPnt(iI, 3), 'ro' );
  
  % compute a plane to fit the contour
  [U_svd, S_svd, V_svd] = svd( planeMatrix, 0 );
  normPnt(iI, :) = V_svd(:,end);
  normCheck(iI) = 1;

  % compare with referDir_section to decide change the normal or not
  if  iI > 1 && iI < length(isovals) && isovals(iI) > startIsoval && isovals(iI) < endIsoval
      angle = acos( dot( normPnt(iI, :), referDir_section ) / ( norm( normPnt(iI, :) ) * norm( referDir_section ) ) );
      if angle > pi/2
	  angle = pi - angle;
      end
      if angle > pi/3
	  normCheck(iI) = 0;
      end
  end
	
  figure(4), plot3( meanPnt(iI, 1), meanPnt(iI, 2), meanPnt(iI, 3), 'ro' );
  hold on
  tmpPnt = meanPnt(iI, :) + 10 .* normPnt(iI, :);
  figure(4), plot3( [ meanPnt(iI, 1) tmpPnt(1) ], [ meanPnt(iI, 2) tmpPnt(2) ], [ meanPnt(iI, 3) tmpPnt(3) ], 'm-');

  % plot the points of the contour onto the fitting plane without changes
  avg_a = normPnt(iI, 1); avg_b = normPnt(iI, 2); avg_c = normPnt(iI, 3);
  avg_abc = sum( normPnt(iI, :) .^2 );
  avg_d = -dot( normPnt(iI, :), meanPnt(iI, :) );
  planeMatrix = vertexOnEdge;
  for ix = 1:size(planeMatrix, 1)
      avg_t0 = ( avg_a * planeMatrix(ix, 1) + avg_b * planeMatrix(ix, 2) + avg_c * planeMatrix(ix, 3) + avg_d ) / avg_abc;
      planeMatrix(ix, 1) = planeMatrix(ix, 1) - avg_a * avg_t0;
      planeMatrix(ix, 2) = planeMatrix(ix, 2) - avg_b * avg_t0;
      planeMatrix(ix, 3) = planeMatrix(ix, 3) - avg_c * avg_t0;
  end
  handle_plane = figure(4); plot3( planeMatrix(:, 1), planeMatrix(:, 2), planeMatrix(:, 3), 'b.', 'MarkerSize', 0.5 );
  
  % compute the closest point and isoval for landmarks
  grid_I = round( ( meanPnt(iI, 1) - boundbox(1) ) / dx );
  grid_J = round( ( meanPnt(iI, 2) - boundbox(2) ) / dy );
  grid_K = round( ( meanPnt(iI, 3) - boundbox(3) ) / dz );
  grid_I = min( max( grid_I, 1 ), sz(1) );
  grid_J = min( max( grid_J, 1 ), sz(2) );
  grid_K = min( max( grid_K, 1 ), sz(3) );
  distTest = norm( [ meanPnt(iI, 1) - y(grid_I, grid_J, grid_K), meanPnt(iI, 2) - x(grid_I, grid_J, grid_K), meanPnt(iI, 3) - z(grid_I, grid_J, grid_K) ] );
  meanIsovals( iI ) = solIm( grid_I, grid_J, grid_K );
  if isnan( meanIsovals( iI ) )
      meanIsovals( iI ) = isovals( iI );
  end
  for ix = 1:size(landmarks, 1)
      if iI == 1
          closestPointsLM( ix, 1 ) = norm( [ landmarks(ix,2)-meanPnt(iI,1), landmarks(ix,1)-meanPnt(iI,2), landmarks(ix,3)-meanPnt(iI,3) ] );
          closestPointsLM( ix, 2 ) = meanIsovals( iI );
      else
          tmpDist = norm( [ landmarks(ix,2)-meanPnt(iI,1), landmarks(ix,1)-meanPnt(iI,2), landmarks(ix,3)-meanPnt(iI,3) ] );
          if tmpDist < closestPointsLM( ix, 1 )
              closestPointsLM( ix, 1 ) = tmpDist;
              closestPointsLM( ix, 2 ) = meanIsovals( iI );
          end
      end
  end
end

figure(2), plot3( landmarks(:,2), landmarks(:,1), landmarks(:,3), 'gx', 'LineWidth', 5 );
view( 90, 0 );
material metal
lighting phong
light
axis equal
axis tight

figure(3), plot3( landmarks(:,2), landmarks(:,1), landmarks(:,3), 'gx', 'LineWidth', 5 );
view( 90, 0 );
axis equal
axis tight

figure(4), plot3( landmarks(:,2), landmarks(:,1), landmarks(:,3), 'gx', 'LineWidth', 5 );
hold on
view( 90, 0 );
axis equal
axis tight

%hold off
%hold off
%hold off

% figName = sprintf( '%s/surface.png', outputIsosurfacePrefix );
% saveas(handle_surface, figName);
figName = sprintf( '%s/%s_surface.fig', outputIsosurfacePrefix, subjectName );
saveas(handle_surface, figName);

%figName = sprintf( '%s/%s_contour.png', outputIsosurfacePrefix, subjectName );
%saveas(handle_contour, figName);
figName = sprintf( '%s/%s_contour.fig', outputIsosurfacePrefix, subjectName );
saveas(handle_contour, figName);

%figName = sprintf( '%s/plane.png', outputIsosurfacePrefix );
%saveas(handle_plane, figName);
%figName = sprintf( '%s/plane.png', outputIsosurfacePrefix );
%saveas(handle_plane, figName);

% using smoothed centerline to check the normal's reliability 
if length(isovals) >= 5
    meanPnt_tmp = zeros( length(isovals)-3, 3 );
    meanPnt_smooth = meanPnt;
    meanPnt_smoothId = zeros( length(isovals) );
    meanPnt_smoothId(1) = 1;
    meanPnt_smoothId(end) = length(isovals);
    for iter = 1:200
    	for ix = 2:length(isovals)-2
            meanPnt_tmp(ix-1, :) = ( meanPnt_smooth(ix, :) + meanPnt_smooth(ix+1, :) ) * 0.5;
        end
	    for ix = 2:length(isovals)-3
            meanPnt_smooth(ix+1, :) = ( meanPnt_tmp(ix-1, :) + meanPnt_tmp(ix, :) ) * 0.5;
        end
    end
    for ix = 2:length(isovals)-1
	% find the closest point on the smoothed centerline for each point on the original centerline
	minDistValue = 0;
	minDistId = 0;
	for iy = 1:length(isovals)
	    distTmp = norm( meanPnt(ix, :) - meanPnt_smooth(iy, :) );
	    if iy == 1 || distTmp < minDistValue 
		minDistValue = distTmp;
		minDistId = iy;
	    end
	end
	meanPnt_smoothId( ix ) = minDistId;
        if minDistId == 1
            referDir_point = meanPnt_smooth(minDistId, :) - meanPnt_smooth(minDistId+1, :);
    	elseif minDistId == length(isovals)
            referDir_point = meanPnt_smooth(minDistId-1, :) - meanPnt_smooth(minDistId, :);
    	else
            referDir_point = ( meanPnt_smooth(minDistId-1, :) - meanPnt_smooth(minDistId+1, :) ) * 0.5;
        end
        referDir_point = referDir_point ./ ( norm( referDir_point) + 1e-15 );
        angle = acos( dot( normPnt(ix, :), referDir_point ) / ( norm( normPnt(ix,:) ) * norm( referDir_point) + 1e-15 ) );
        if angle > pi/2
            angle = pi - angle;
        end 
        if angle > pi/3
            normCheck(ix) = 0;
            %fprintf( 'refer point, angle: %f\n', angle );
        end
    end 
end

normCheckRecorded = normCheck;

ix = 1;
while ix <= length(isovals)
    while ix <= length(isovals) && normCheck(ix) == 1
        ix = ix + 1;
    end
    if ix > length(isovals)
        break;
    end
    startId = ix;
    endId = startId;
    while endId <= length(isovals) && normCheck(endId) == 0
        endId = endId + 1; 
    end	
    endId = endId - 1;
    % use the quaternion to interplate the normal
    norm1 = normPnt( max( 1, startId-1 ), : );
    norm2 = normPnt( min( length(isovals), endId+1 ), : );
    angle = acos( dot( norm1, norm2 ) / ( norm( norm1 ) * norm( norm2 ) + 1e-15 ) );
    if angle > pi/2
        norm2 = -norm2;
        angle = pi - angle;
    end
    
    distTmp = zeros( min( length(isovals), endId+1 ) - max( 1, startId-1 ) + 1, 1 );
    for iy = 1:length(distTmp)
        if iy == 1
            distTmp(iy) = 0;
            continue;
        end
        distTmp(iy) = norm( meanPnt( max( 1, startId-1 )+iy-2, : ) - meanPnt( max( 1, startId-1)+iy-1, : ) ) + distTmp(iy-1); 
    end
    distTmp = distTmp ./ distTmp(end);
    
    axis_v = cross( norm1, norm2 );
    axis_v = axis_v / ( norm( axis_v ) + 1e-15 );
    for iy = startId:endId
        rotateAngle = angle * distTmp( iy - max(1, startId-1) + 1 );
        quaternion = [ cos(rotateAngle/2), sin(rotateAngle/2) .* axis_v ];
        quaternion_inv = [ quaternion(1), -quaternion(2), -quaternion(3), -quaternion(4) ] / ( norm( quaternion ) ^2); 
        vector = [ 0, norm1 ];
        tmpVector(1) = quaternion(1) * vector(1) - dot( quaternion(2:4), vector(2:4) );
        tmpVector(2:4) = quaternion(1) * vector(2:4) + vector(1) * quaternion(2:4) ...
                       + cross( quaternion(2:4), vector(2:4) );
        newVector(1) = tmpVector(1) * quaternion_inv(1) - dot( tmpVector(2:4), quaternion_inv(2:4) );
        newVector(2:4) = tmpVector(1) * quaternion_inv(2:4) + quaternion_inv(1) * tmpVector(2:4) ...
                       + cross( tmpVector(2:4), quaternion_inv(2:4) );
	referId = meanPnt_smoothId( iy );
	if referId == 1
            referDir_point = meanPnt_smooth(referId, :) - meanPnt_smooth(referId+1, :);
        elseif minDistId == length(isovals)
            referDir_point = meanPnt_smooth(referId-1, :) - meanPnt_smooth(referId, :);
        else
            referDir_point = ( meanPnt_smooth(referId-1, :) - meanPnt_smooth(referId+1, :) ) * 0.5;
        end

	newAngle = acos( dot( newVector(2:4), referDir_point ) / ( norm( newVector(2:4) ) * norm( referDir_point ) + 1e-15 ) );
	if newAngle > pi/2
	     newAngle = pi - newAngle;
	end
	if newAngle > pi/3 && ( startId-1 > 1 || endId+1 < length(isovals) )
	     for iz = startId:endId
		 normCheck( iz ) = 0;
	     end
	     if startId-1 > 1
	         normCheck( startId-1 ) = 0;
	     end
	     if endId+1 < length(isovals)
		 normCheck( endId+1 ) = 0;
	     end
	     if startId-1 > 1
		 endId = startId-2;
	     else
		 endId = startId-1;
	     end
	     break;
	end
	
	normCheck(iy) = 1;
        normPnt(iy, :) = newVector(2:4);
        normPnt(iy, :) = normPnt(iy, :) / ( norm( normPnt(iy, :) ) + 1e-15 );
        normPnt(iy, :)
        figure(4), plot3( meanPnt(iy, 1), meanPnt(iy, 2), meanPnt(iy, 3), 'co', 'LineWidth', 2 );
        hold on
        tmpPnt = meanPnt(iy, :) + 10 .* normPnt(iy, :);
        figure(4), plot3( [ meanPnt(iy, 1) tmpPnt(1) ], [ meanPnt(iy, 2) tmpPnt(2) ], [ meanPnt(iy, 3) tmpPnt(3) ], 'c-', 'LineWidth', 2);
    end
    ix = endId + 1;
%     % using spline to interpolate the normal
%     if startId-2 >= 1 && normCheck(startId-2) == 1
%         splineStartId = startId-2;
%     else
%         splineStartId = max( 1, startId-1 );
%     end
%     if endId+2 <= length(isovals) && normCheck(endId+2) == 1
%         splineEndId = endId+2;
%     else
%         splineEndId = min( length(isovals), endId+1 );
%     end
%     %sf = spline( 1:splineEndId-splineStartId+1, meanPnt(splineStartId:splineEndId, :)' );
%     sf = spline( 1:splineEndId-splineStartId-endId+startId, [ meanPnt(splineStartId:startId-1, :); meanPnt(endId+1:splineEndId, :) ]' );
%     startId
%     endId
%     stepLength = 1.0 / ( endId - startId + 2 );
%     for iy = startId:endId
%         t_para = startId - splineStartId + (iy-startId+1) * stepLength;
%         meanPnt(iy, :) = ppval( sf, t_para ); 
%         dt_pre  = max( startId-splineStartId, t_para - 0.1 );
%         dt_post = min( startId-splineStartId+1, t_para + 0.1 );
%         pt_pre = ppval( sf, dt_pre );
%         pt_post = ppval( sf, dt_post );
%         normPnt(iy, :) = pt_pre - pt_post;
%         normPnt(iy, :) = normPnt(iy, :) / ( norm( normPnt(iy,:) ) + 1e-15 );
%         normPnt(iy, :)
% 	
%         figure(4), plot3( meanPnt(iy, 1), meanPnt(iy, 2), meanPnt(iy, 3), 'ko' );
%         tmpPnt = meanPnt(iy, :) + 10 .* normPnt(iy, :);
%         figure(4), plot3( [ meanPnt(iy, 1) tmpPnt(1) ], [ meanPnt(iy, 2) tmpPnt(2) ], [ meanPnt(iy, 3) tmpPnt(3) ], 'k-', 'LineWidth', 2);
%     end
%     ix = endId+1;
end

%figName = sprintf( '%s/%s_plane.png', outputIsosurfacePrefix, subjectName );
%saveas(handle_plane, figName);
figName = sprintf( '%s/%s_plane.fig', outputIsosurfacePrefix, subjectName );
saveas(handle_plane, figName);

meanNormFile_plane = sprintf( '%s/%s_MeanAndNormal.txt', outputIsosurfacePrefix, subjectName );
fidMeanNormPlane = fopen( meanNormFile_plane, 'wt' );
fprintf( fidMeanNormPlane, '%d\n', length( isovals ) );
for iI = 1:length(isovals)
    fprintf( fidMeanNormPlane, '%f %f %f %f %f %f\n', meanPnt(iI,1), meanPnt(iI,2), meanPnt(iI,3), normPnt(iI,1), normPnt(iI,2), normPnt(iI,3) );
end
fclose( fidMeanNormPlane);

normFile_reliableCheck = sprintf( '%s/%s_NormReliableCheck.txt', outputIsosurfacePrefix, subjectName );
fidNormCheck = fopen( normFile_reliableCheck, 'wt' );
fprintf( fidNormCheck, '%d\n', length( isovals ) );
for iI = 1:length(isovals)
    fprintf( fidNormCheck, '%d\n', normCheckRecorded(iI) );
end
fclose( fidNormCheck );

% find the id for landmarks' isoval
txtName = sprintf( '%s/%s_LandmarksIdOnCenterline.txt', outputIsosurfacePrefix, subjectName );
fidTxtFile = fopen( txtName, 'wt' );
fprintf( fidTxtFile, '%d\n', size(landmarks, 1 ) );
for ix = 1:size(landmarks,1)
    if isnan( landmarksIsoval( ix ) )
	landmarksIsoval( ix ) = closestPointsLM( ix, 2 );
    end
    [minValue, minId] = min( abs( meanIsovals - landmarksIsoval( ix ) ) );
    if ix == 1
        minId = 1;
    end
    if ix == size( landmarks, 1 )
        minId = length( isovals );
    end
    fprintf( fidTxtFile, '%d\n', minId );
end
fclose( fidTxtFile );

hold off
hold off
hold off

