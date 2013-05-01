// ComputeArea.cpp : Defines the entry point for the console application.
//

#include "vtkCutter.h"
#include "vtkPlane.h"
#include "vtkPolyDataReader.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkSmartPointer.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCellArray.h"
#include "vtkPolygon.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "computeAreaAndContourWithMeanNormCLP.h"
#include "vtkCellLocator.h"
#include <math.h>
#include <algorithm>
#include <string.h> 
using namespace std;

#define VTK_CREATE(type, var) vtkSmartPointer<type> var = vtkSmartPointer<type>::New()
typedef struct {
	double x;
	double y;
	double z;
}Points;

typedef struct {
	int id1;
	int id2;
	bool bFlag;
}Lines;

int main(int argc, char **argv) 
{
	PARSE_ARGS;

	/*
	string szFilePathName;
	//CString szFileName;
	string szExtName;
	szFilePathName = inputModelVTK.c_str();
	//szFileName = szFilePathName.Mid( szFilePathName.ReverseFind('\\') + 1 );
	szExtName = szFilePathName.Mid( szFilePathName.ReversEFind('.') + 1 );
	*/

	//string strFileName = inputModelVTK.c_str();
	int nPosExt = inputModelVTK.rfind( "." );
	string strExtName = inputModelVTK.substr( nPosExt + 1 );

	std::cout << nPosExt << " " << strExtName << std::endl;
	
	//vtkPolyData* airway = NULL;
	VTK_CREATE( vtkPolyData, airway );
	if( !strcmp( strExtName.c_str(), "vtp" ) )
	{
		// read .vtp model
		std::cout << " read .vtp model: " << inputModelVTK << std::endl;
		VTK_CREATE(vtkXMLPolyDataReader, vtkReader);
		//vtkXMLPolyDataReader * vtkReader = vtkXMLPolyDataReader::New();
		vtkReader->SetFileName(inputModelVTK.c_str());
		vtkReader->Update();
		airway->CopyStructure( vtkReader->GetOutput() );
		if( airway == NULL ) std::cout << "Airway is null." << std::endl;
		//vtkReader->Delete();
	}
	else if( !strcmp( strExtName.c_str(), "vtk" ) )
	{
		std::cout << " read .vtk model ... " << std::endl;
		VTK_CREATE( vtkPolyDataReader, vtkReader );
		//vtkPolyDataReader * vtkReader = vtkPolyDataReader::New();
		vtkReader->SetFileName( inputModelVTK.c_str() );
		vtkReader->Update();
		airway->CopyStructure( vtkReader->GetOutput() );
		//vtkReader->Delete();
	}
	else
	{
		std::cout << "Unknown file format!" << std::endl;
		return 0;
	}
		
	// cell locator and intersection parameters
	/*VTK_CREATE( vtkCellLocator, locator );
	locator->SetDataSet(airway);
	locator->CacheCellBoundsOn();
	locator->BuildLocator();
	double tolerance = 0.001, lineParameter, intersect[3], paraCoord[3];
	int sub_id;*/

	double bounds[6];
	airway->GetBounds( bounds );
	std::cout << bounds[0] << " " << bounds[1] << " " << bounds[2] << " " 
		   << bounds[3] << " " << bounds[4] << " " << bounds[5] << std::endl; 

	// read points from center and norm file
	std::cout << "Read points from centerline file ... " << std::endl;
	const char * meanPointsFile = inputMeanNorm.c_str();
	FILE * pFile = fopen(meanPointsFile, "r");
	char strLine[256];
	fgets(strLine, 256, pFile);
	int nNumPoints;
	sscanf(strLine, "%d", &nNumPoints);
	printf("The number of points on centerline: %d\n", nNumPoints);
	std::vector<Points> vecCenter;
	std::vector<Points> vecNorm;
	for(int i=0; i<nNumPoints; i++)
	{
		if(feof(pFile)) break;
		fgets(strLine, 256, pFile);
		Points tmpCenter;
		Points tmpNorm;
		sscanf(strLine, "%lf %lf %lf %lf %lf %lf", &tmpCenter.x, &tmpCenter.y, &tmpCenter.z, &tmpNorm.x, &tmpNorm.y, &tmpNorm.z);
		vecCenter.push_back(tmpCenter);
		vecNorm.push_back(tmpNorm);
	}
	fclose(pFile);

 	// read mouth's position from the file
	std::cout << "Read mouth's position from the file ... " << std::endl;
	double mouthCenter[3];
	double mouthRadius;
	bool bMouth = false;
	const char * mouthPosFile = inputMouthPosition.c_str();
	pFile = fopen( mouthPosFile, "r" );
	if( pFile )
	{
		bMouth = true;
		while( fgets( strLine, 256, pFile ) )
		{
			char strKey[256];
			sscanf( strLine, "%s :", strKey);
			if( !strcmp( strKey, "ClipSphereCenter" ) )
			{
				sscanf( strLine, "%s : %lf %lf %lf", strKey, &mouthCenter[0], &mouthCenter[1], &mouthCenter[2] );
			}
			else if( !strcmp( strKey, "ClipSphereRadius" ) )
			{
				sscanf( strLine, "%s : %lf", strKey, &mouthRadius );
			}
		}
		fclose( pFile );
	}
	if( bMouth )
	{
		std::cout<< mouthCenter[0] << " " << mouthCenter[1] << " " << mouthCenter[2] << std::endl;
		std::cout<< mouthRadius << std::endl;
	}

	// read landmarks for starting searching the shortest parth
/*	std::cout << "Read landmarks ... " << std::endl;
	double sourcePntLeft[3];
	double sourcePntRight[3];
	int bSourcePnts = 0;
	const char * landmarksFile = inputLandmarkPosition.c_str();
	pFile = fopen( landmarksFile, "r" );
	if( pFile )
	{
		while( fgets( strLine, 256, pFile ) )
		{
			char strKey[256];
			sscanf( strLine, "%s :", strKey );
			if( !strcmp( strKey, "LeftAlaRim" ) )
			{	
				sscanf( strLine, "%s : %lf %lf %lf", strKey, &sourcePntLeft[0], &sourcePntLeft[1], &sourcePntLeft[2] );
				bSourcePnts++;
			}
			else if( !strcmp( strKey, "RightAlaRim" ) )
			{
				sscanf( strLine, "%s : %lf %lf %lf", strKey, &sourcePntRight[0], &sourcePntRight[1], &sourcePntRight[2] );
				bSourcePnts++;
			}
		}
		if( bSourcePnts < 2 )
		{
			std::cout << "Left or/and right alaRims are missing." << std::endl;
		}
		else
		{
			std::cout << "Left: "  << sourcePntLeft[0] << ", " << sourcePntLeft[1] << ", " << sourcePntLeft[2] << ", "
				  << "Right: " << sourcePntRight[0] << ", " << sourcePntRight[1] << ", " << sourcePntRight[2] << std::endl; 
		}
		fclose( pFile );
	}*/
	
	// vtkCutter 
	VTK_CREATE(vtkCutter, cutter);
	cutter->SetInput(airway);
	VTK_CREATE(vtkPlane, plane);
	
	const char * areaFile = outputArea.c_str();
	FILE *pFileArea = fopen( areaFile, "wt" );
	fprintf( pFileArea, "%d\n", nNumPoints );
	// for every point on centerline, get a cut
	std::cout << " Begin to cut ... " << std::endl;
	for(int i=0; i<nNumPoints; i++)
	{
		printf("Begin to cut %d\n", i+1);
		plane->SetOrigin(vecCenter[i].x, vecCenter[i].y, vecCenter[i].z);
		plane->SetNormal(vecNorm[i].x, vecNorm[i].y, vecNorm[i].z);
		cutter->SetCutFunction(plane);
		cutter->Update();

                // compute the connected regions
		VTK_CREATE( vtkPolyDataConnectivityFilter, connectFilter );
		connectFilter->SetExtractionModeToAllRegions();
		connectFilter->SetInputConnection( cutter->GetOutputPort() );
		connectFilter->Update();

		int nbContours = connectFilter->GetNumberOfExtractedRegions();
		// nbContours, connectFilter->GetOutput()->GetNumberOfPoints(), connectFilter->GetOutput()->GetNumberOfVerts(), connectFilter->GetOutput()->GetNumberOfLines(), connectFilter->GetOutput()->GetNumberOfPolys(), connectFilter->GetOutput()->GetNumberOfStrips() );
		

		char fName[256];
                sprintf(fName, "%s%03d.txt", outputContourPrefix.c_str(), i+1);
                FILE * pFileTmp = fopen(fName, "wt");
                fprintf(pFileTmp, "%d %d\n", connectFilter->GetOutput()->GetNumberOfPoints(), nbContours);
	        
 		sprintf(fName, "%sOriginal%03d.txt", outputContourPrefix.c_str(), i+1);
                FILE * pFileOriginal = fopen(fName, "wt");
                fprintf(pFileOriginal, "%d %d\n", connectFilter->GetOutput()->GetNumberOfPoints(), nbContours);
	
		connectFilter->SetExtractionModeToSpecifiedRegions();
		double dCrossSectionalArea = 0;
		Lines *pContourLines = new Lines[ connectFilter->GetOutput()->GetNumberOfLines() ];
		int nTotalOfPoints = connectFilter->GetOutput()->GetNumberOfPoints();
                int *pContourIdPoints = new int[ nTotalOfPoints ];
		Points * pContourPoints = new Points[ nTotalOfPoints ];
		int * pContourPointsFlag = new int[ nTotalOfPoints ];
		memset( pContourPointsFlag, 0, sizeof(int)*nTotalOfPoints );

	     	std::cout << "nTotalOfPoints: " << nTotalOfPoints << " , nLines: " << connectFilter->GetOutput()->GetNumberOfLines() << std::endl;
		double * pXSectionalArea = new double[ nbContours ];
		bool * pFlagShortestPath = new bool[ nbContours ];
		memset( pFlagShortestPath, 0, sizeof(bool)*nbContours );
		Points * pMinBoundary = new Points[ nbContours ];
		Points * pMaxBoundary = new Points[ nbContours ];
		//double * pDistLeft = new double[ nbContours ];
		//Points * pNewSourceLeft = new Points[ nbContours ];
		//Points * pNewSourceRight = new Points[ nbContours ];
		//double * pDistRight = new double[ nbContours ];
		int * pPosPoints = new int [ nbContours ];
		//int bFlagLeft = 0;
		//int bFlagRight = 0;
		for( int idx = 0; idx < nbContours; idx++ )
		{
			connectFilter->InitializeSpecifiedRegionList();
			connectFilter->AddSpecifiedRegion( idx );
			connectFilter->Update();
		 	
			VTK_CREATE( vtkCellArray, lines );
                	lines = connectFilter->GetOutput()->GetLines();
			int nbLines = lines->GetNumberOfCells();
                        int nbPoints = connectFilter->GetOutput()->GetNumberOfPoints();
			
			std::cout << "nbLines: " << nbLines << ", nbPoints: " << nbPoints << std::endl;
			
			if( nbLines <= 1 ) continue;
			vtkIdType npts, *pts;
			int nLineCount = 0;
			fprintf(pFileOriginal, "%d\n", nbLines );
                	for( lines->InitTraversal(); lines->GetNextCell( npts, pts); nLineCount++ )
                	{		
				pContourLines[nLineCount].id1 = pts[0];
				pContourLines[nLineCount].id2 = pts[1];
				pContourLines[nLineCount].bFlag = 0;
				pContourPointsFlag[ pts[0] ]++;
				pContourPointsFlag[ pts[1] ]++;
				fprintf(pFileOriginal, "%d %d\n", pts[0], pts[1] );
                	}

			 // deal with the same edges
                        for( int iM = 0; iM < nLineCount; iM++ )
                        {
                                for( int iN = 0; iN < nLineCount; iN++ )
                                {
                                        if( iN == iM ) continue;
                                        if( pContourLines[iM].id1 == pContourLines[iN].id1 && pContourLines[iM].id2 == pContourLines[iN].id2 || pContourLines[iM].id1 == pContourLines[iN].id2 && pContourLines[iM].id2 == pContourLines[iN].id1 )
                                        {
                                                pContourLines[iM].bFlag = 1;
                                                pContourLines[iN].bFlag = 1;
                                        }
                                }
                        }

			std::cout << "nLineCounts: " << nLineCount << std::endl;
						
			int curPos = 0;
			pContourIdPoints[curPos] = pContourLines[0].id1;
			curPos++;
			pContourIdPoints[curPos] = pContourLines[0].id2;
			pContourLines[0].bFlag = 1;
		
			int nCountLines = 1;
			while( nCountLines < nbLines )
			{
				bool bFound = 0;
				for( int iTmp = 0; iTmp < nbLines; iTmp++ )
				{
					if( pContourLines[iTmp].bFlag == 1  ) continue;
					if( pContourLines[iTmp].id1 == pContourIdPoints[curPos] )
					{
						nCountLines++;
						pContourLines[iTmp].bFlag = 1;
						if( pContourLines[iTmp].id2 != pContourIdPoints[0] )
						{
							curPos++;
							pContourIdPoints[curPos] = pContourLines[iTmp].id2;
						}
						bFound = 1;
						break;
						
					}
					else if( pContourLines[iTmp].id2 == pContourIdPoints[curPos] )
					{
						nCountLines++;
						pContourLines[iTmp].bFlag = 1;
						if( pContourLines[iTmp].id1 != pContourIdPoints[0] )
						{
							curPos++;
							pContourIdPoints[curPos] = pContourLines[iTmp].id1;
						}
						bFound = 1;
						break;
					}
				}
				if( bFound == 0 ) 
				{
					if( pContourPointsFlag[ pContourIdPoints[curPos] ] == 1 )
					{
						for( int iTmp = 0; iTmp < nTotalOfPoints; iTmp++ )
						{
							if( pContourIdPoints[curPos] != iTmp && pContourPointsFlag[ iTmp ] == 1 ) 
							{
								std::cout << iTmp << " " << pContourIdPoints[curPos] << std::endl;
								if( iTmp != pContourIdPoints[0] )
								{
									curPos++;
									pContourIdPoints[curPos] = iTmp;
								}
								pContourPointsFlag[ pContourIdPoints[curPos] ] = 2;
								pContourPointsFlag[ iTmp ] = 2;
								bFound = 1;
								break;
							}
						}
					} 
				}
				if( bFound == 0 ) 
				{
					std::cout << "Error: Couldn't find a loop." <<std::endl;
					break;
				}
			}
			
			for( int iTmp = 0; iTmp < nTotalOfPoints; iTmp++ )
			{
				if( pContourPointsFlag[iTmp] == 1 ) pContourPointsFlag[iTmp] = 2;
			}
			
			std::cout << "Begin to compute area... " << curPos+1 << std::endl;
			int nPointsOnContour = curPos+1;

			// record the position of the last point in current contour
			if( idx == 0 ) pPosPoints[idx] = curPos;
			else pPosPoints[idx] = pPosPoints[idx-1] + nPointsOnContour;
			// start position of the current contour
			int nPosTmp;
			if( idx == 0 ) nPosTmp = 0;
			else nPosTmp = pPosPoints[idx-1] + 1;

			bool bInMouth = true;
			double dXmean = 0;
			double dYmean = 0;
			double dZmean = 0;
			for( int tmp = 0; tmp < nPointsOnContour; tmp++ )
                        {
                                double pointPos[3];
                                connectFilter->GetOutput()->GetPoint( pContourIdPoints[tmp], pointPos );
                                pContourPoints[tmp+nPosTmp].x = pointPos[0];
				pContourPoints[tmp+nPosTmp].y = pointPos[1];
				pContourPoints[tmp+nPosTmp].z = pointPos[2];

				if( sqrt( ( pointPos[0] - mouthCenter[0] ) * ( pointPos[0] - mouthCenter[0] ) + 
					  ( pointPos[1] - mouthCenter[1] ) * ( pointPos[1] - mouthCenter[1] ) + 
					  ( pointPos[2] - mouthCenter[2] ) * ( pointPos[2] - mouthCenter[2] ) ) > mouthRadius )
				{
					bInMouth = false;
				}

				dXmean += pointPos[0];
				dYmean += pointPos[1];
				dZmean += pointPos[2];
                        }
			if( nPointsOnContour != 0 )
			{
				dXmean /= nPointsOnContour;
				dYmean /= nPointsOnContour;
				dZmean /= nPointsOnContour;
			}
			double distTmp = sqrt( ( vecCenter[i].x - dXmean ) * ( vecCenter[i].x - dXmean ) + 
					       ( vecCenter[i].y - dYmean ) * ( vecCenter[i].y - dYmean ) + 
					       ( vecCenter[i].z - dZmean ) * ( vecCenter[i].z - dZmean ) ); 

			std::cout<< "dist " << idx << ": " << distTmp << std::endl;
			std::cout << "nbCounter: " << nbContours << " idx: " << idx << std::endl;
			pXSectionalArea[idx] = 0;
			pFlagShortestPath[idx] = 0;
			pMinBoundary[idx].x = 0; pMinBoundary[idx].y = 0; pMinBoundary[idx].z = 0;
			pMaxBoundary[idx].x = 0; pMaxBoundary[idx].y = 0; pMaxBoundary[idx].z = 0;
			if( !bInMouth && distTmp < 40 )
			{ 
				std::cout << "points on Contour: " << nPointsOnContour << std::endl;		
				VTK_CREATE( vtkPolygon, cutPolygon );
				//vtkPolygon * cutPolygon = vtkPolygon::New();
				cutPolygon->GetPointIds()->SetNumberOfIds( nPointsOnContour );
				cutPolygon->GetPoints()->SetNumberOfPoints( nPointsOnContour );
				//fprintf( pFileTmp, "%d\n", nPointsOnContour );
				for( int tmp = 0; tmp < nPointsOnContour; tmp++ )
				{
					//fprintf( pFileTmp, "%lf %lf %lf\n", pContourPoints[tmp].x, pContourPoints[tmp].y, pContourPoints[tmp].z ); 
					// record the boundary box for current contour
					if( tmp == 0 || pContourPoints[tmp+nPosTmp].x < pMinBoundary[idx].x ) pMinBoundary[idx].x = pContourPoints[tmp+nPosTmp].x;
					if( tmp == 0 || pContourPoints[tmp+nPosTmp].y < pMinBoundary[idx].y ) pMinBoundary[idx].y = pContourPoints[tmp+nPosTmp].y;
					if( tmp == 0 || pContourPoints[tmp+nPosTmp].z < pMinBoundary[idx].z ) pMinBoundary[idx].z = pContourPoints[tmp+nPosTmp].z;
					if( tmp == 0 || pContourPoints[tmp+nPosTmp].x > pMaxBoundary[idx].x ) pMaxBoundary[idx].x = pContourPoints[tmp+nPosTmp].x;
					if( tmp == 0 || pContourPoints[tmp+nPosTmp].y > pMaxBoundary[idx].y ) pMaxBoundary[idx].y = pContourPoints[tmp+nPosTmp].y;
					if( tmp == 0 || pContourPoints[tmp+nPosTmp].z > pMaxBoundary[idx].z ) pMaxBoundary[idx].z = pContourPoints[tmp+nPosTmp].z;
					cutPolygon->GetPointIds()->SetId( tmp, tmp );
					cutPolygon->GetPoints()->SetPoint( tmp, pContourPoints[tmp+nPosTmp].x, pContourPoints[tmp+nPosTmp].y, pContourPoints[tmp+nPosTmp].z );
				}
				pFlagShortestPath[idx] = 1;
			}
		}
		
		/*double dMaxAreaTmp = 0;
		int nMaxAreaIdTmp = -1;
		double dSecAreaTmp = 0;
		int nSecAreaIdTmp = -1;
		for( int idx = 0; idx < nbContours; idx++ )
		{
			if( pXSectionalArea[idx] > dMaxAreaTmp ) 
			{
				nSecAreaIdTmp = nMaxAreaIdTmp;
				dSecAreaTmp = dMaxAreaTmp;
				nMaxAreaIdTmp = idx;
				dMaxAreaTmp = pXSectionalArea[idx];
			}
			else if( pXSectionalArea[idx] > dSecAreaTmp )
			{
				nSecAreaIdTmp = idx;
				dSecAreaTmp = pXSectionalArea[idx];
			}
			if( pXSectionalArea[idx] > 100 ) pFlagShortestPath[idx] = 1;
		}	
		pFlagShortestPath[nMaxAreaIdTmp] = 1;
		if( dSecAreaTmp / dMaxAreaTmp >= 0.2 )	
		{
			pFlagShortestPath[nSecAreaIdTmp] = 1;
			for( int idx = 0; idx < nbContours; idx++ )
			{
				if( idx == nMaxAreaIdTmp || idx == nSecAreaIdTmp ) continue;
				if( pXSectionalArea[idx] / ( dSecAreaTmp + 1e-10 ) >= 0.5 ) pFlagShortestPath[idx] = 1;
			}
		}*/
		
		double firstMaxArea = 0; 
		double secondMaxArea = 0;
		int idFirstMaxArea = -1;
		int idSecondMaxArea = -1;
		for( int idx = 0; idx < nbContours; idx++ )
		{
			if( pXSectionalArea[idx] * pFlagShortestPath[idx] > firstMaxArea )
			{
				secondMaxArea = firstMaxArea;
				idSecondMaxArea = idFirstMaxArea;
				firstMaxArea = pXSectionalArea[idx] * pFlagShortestPath[idx];
				idFirstMaxArea = idx;
			}
			else if( pXSectionalArea[idx] * pFlagShortestPath[idx] > secondMaxArea )
                        {
                                secondMaxArea = pXSectionalArea[idx] * pFlagShortestPath[idx];
                                idSecondMaxArea = idx;
                        }		
		}
			
		std::cout << "test1 ..." << std::endl;

		if( idFirstMaxArea >= 0 )
		{
			int nPosTmp, nPntNumberTmp;
                        if( idFirstMaxArea == 0 ) { nPosTmp = 0; nPntNumberTmp = pPosPoints[idFirstMaxArea] + 1; }
                        else { nPosTmp = pPosPoints[idFirstMaxArea-1]+1; nPntNumberTmp = pPosPoints[idFirstMaxArea] - pPosPoints[idFirstMaxArea-1]; }
                        VTK_CREATE( vtkPolygon, polygonTmp );
                        polygonTmp->GetPointIds()->SetNumberOfIds( nPntNumberTmp );
                        polygonTmp->GetPoints()->SetNumberOfPoints( nPntNumberTmp );
			for( int tmp = 0; tmp < nPntNumberTmp; tmp++ )
                        {
                        	polygonTmp->GetPointIds()->SetId( tmp, tmp );
                               	polygonTmp->GetPoints()->SetPoint( tmp, pContourPoints[tmp+nPosTmp].x, pContourPoints[tmp+nPosTmp].y, pContourPoints[tmp+nPosTmp].z );
                        }
			std::cout << "id first: " << idFirstMaxArea << " " << nPosTmp << " " << nPntNumberTmp << std::endl;
			for( int idx = 0; idx < nbContours; idx++ )
			{
				if( idx == idFirstMaxArea ) continue;
				if( pMinBoundary[idx].x >= pMinBoundary[idFirstMaxArea].x && pMinBoundary[idx].y >= pMinBoundary[idFirstMaxArea].y && 
				    pMinBoundary[idx].z >= pMinBoundary[idFirstMaxArea].z && pMaxBoundary[idx].x <= pMaxBoundary[idFirstMaxArea].x && 
				    pMaxBoundary[idx].y <= pMaxBoundary[idFirstMaxArea].y && pMaxBoundary[idx].z <= pMaxBoundary[idFirstMaxArea].z )
				{
					double pntTmpTest[3];
					pntTmpTest[0] = pContourPoints[ pPosPoints[idx] ].x;
					pntTmpTest[1] = pContourPoints[ pPosPoints[idx] ].y;
					pntTmpTest[2] = pContourPoints[ pPosPoints[idx] ].z;
					double normTmp[3];
					normTmp[0] = vecNorm[i].x;
					normTmp[1] = vecNorm[i].y;
					normTmp[2] = vecNorm[i].z;
					if( polygonTmp->PointInPolygon( pntTmpTest, (int)(polygonTmp->GetPoints()->GetNumberOfPoints()), (double *)(polygonTmp->GetPoints()), polygonTmp->GetPoints()->GetBounds(), normTmp ) )
					{
						pFlagShortestPath[idx] = 1;
						pXSectionalArea[idx] *= (-1);		
					}
				}
			}
		}

		std::cout << "test2 ..." << std::endl;

		if( idSecondMaxArea >= 0 )
		{
			int nPosTmp, nPntNumberTmp;
                       	if( idSecondMaxArea == 0 ) { nPosTmp = 0; nPntNumberTmp = pPosPoints[idSecondMaxArea] + 1; }
                        else { nPosTmp = pPosPoints[idSecondMaxArea-1]+1; nPntNumberTmp = pPosPoints[idSecondMaxArea] - pPosPoints[idSecondMaxArea-1]; }
                        VTK_CREATE( vtkPolygon, polygonTmp );
                        polygonTmp->GetPointIds()->SetNumberOfIds( nPntNumberTmp );
                        polygonTmp->GetPoints()->SetNumberOfPoints( nPntNumberTmp );
			for( int tmp = 0; tmp < nPntNumberTmp; tmp++ )
                        {
                        	polygonTmp->GetPointIds()->SetId( tmp, tmp );
                                polygonTmp->GetPoints()->SetPoint( tmp, pContourPoints[tmp+nPosTmp].x, pContourPoints[tmp+nPosTmp].y, pContourPoints[tmp+nPosTmp].z );
                        }

			std::cout << "id first: " << idSecondMaxArea << " " << nPosTmp << " " << nPntNumberTmp << std::endl;

			for( int idx = 0; idx < nbContours; idx++ )
			{
				if( idx == idSecondMaxArea ) continue;
				if( pMinBoundary[idx].x >= pMinBoundary[idSecondMaxArea].x && pMinBoundary[idx].y >= pMinBoundary[idSecondMaxArea].y &&
				    pMinBoundary[idx].z >= pMinBoundary[idSecondMaxArea].z && pMaxBoundary[idx].x <= pMaxBoundary[idSecondMaxArea].x &&
				    pMaxBoundary[idx].y <= pMaxBoundary[idSecondMaxArea].y && pMaxBoundary[idx].z <= pMaxBoundary[idSecondMaxArea].z )
				{
                                        double pntTmpTest[3];
                                        pntTmpTest[0] = pContourPoints[ pPosPoints[idx] ].x;
                                        pntTmpTest[1] = pContourPoints[ pPosPoints[idx] ].y;
                                        pntTmpTest[2] = pContourPoints[ pPosPoints[idx] ].z;
                                        double normTmp[3];
                                        normTmp[0] = vecNorm[i].x;
                                        normTmp[1] = vecNorm[i].y;
                                        normTmp[2] = vecNorm[i].z;
                                        if( polygonTmp->PointInPolygon( pntTmpTest, (int)(polygonTmp->GetPoints()->GetNumberOfPoints()), (double *)(polygonTmp->GetPoints()), polygonTmp->GetPoints()->GetBounds(), normTmp ) )
                                        {
						pFlagShortestPath[idx] = 1;
						pXSectionalArea[idx] *= (-1);
					}
				}
			}
		}

		std::cout << "test3 ..." << std::endl;

		if( idFirstMaxArea < 0  && idSecondMaxArea < 0 )
		{
			for( int idx = 0; idx < nbContours; idx++ )
			{
				pFlagShortestPath[idx] = 1;
			}
		}

		std::cout << "test 4 ..." << std::endl;
		
		dCrossSectionalArea = 0;
		for( int idx = 0; idx < nbContours; idx++ )
		{
			if( pXSectionalArea[idx] != 0 && pFlagShortestPath[idx] == 1 )
			{
				dCrossSectionalArea += pXSectionalArea[idx];
				if( idx == 0 )
				{
					fprintf( pFileTmp, "%d\n", pPosPoints[idx] + 1 );
					for( int tmp = 0; tmp <= pPosPoints[idx]; tmp++ )
					{
						fprintf( pFileTmp, "%lf %lf %lf\n", pContourPoints[tmp].x, pContourPoints[tmp].y, pContourPoints[tmp].z );
					}
				}
				else 
				{
					fprintf( pFileTmp, "%d\n", pPosPoints[idx] - pPosPoints[idx-1] );
                                	for( int tmp = pPosPoints[idx-1]+1; tmp <= pPosPoints[idx]; tmp++ )
                                	{
                                        	fprintf( pFileTmp, "%lf %lf %lf\n", pContourPoints[tmp].x, pContourPoints[tmp].y, pContourPoints[tmp].z ); 
					}
				}
			}
		}

		std::cout << "test5 ..." << std::endl;

		if( pContourIdPoints ) { delete [] pContourIdPoints; pContourIdPoints = NULL; }
                if( pContourLines ) { delete [] pContourLines; pContourLines = NULL; }
                if( pContourPoints ) { delete [] pContourPoints; pContourPoints = NULL; }
		if( pContourPointsFlag ) { delete [] pContourPointsFlag; pContourPointsFlag = NULL; }
		if( pXSectionalArea ) { delete [] pXSectionalArea; pXSectionalArea = NULL; }
		if( pFlagShortestPath ) { delete [] pFlagShortestPath; pFlagShortestPath = NULL; }
		if( pMinBoundary ) { delete [] pMinBoundary; pMinBoundary = NULL; }
		if( pMaxBoundary ) { delete [] pMaxBoundary; pMaxBoundary = NULL; }
		if( pPosPoints ) { delete [] pPosPoints; pPosPoints = NULL; }
			
                fprintf( pFileArea, "%lf\n", dCrossSectionalArea);
		fclose( pFileTmp );
		fclose( pFileOriginal );
		printf("point %d finished\n", i+1);
	}
	fclose( pFileArea );	
	return 0;
}

