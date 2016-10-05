////////////////////////////////////////////////////////////////////////////
// Display the Noise level obtained after multi scales analysis 
//////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include <vector>
#include <algorithm>
#include "ImaGene/base/Proxy.h"
#include "ImaGene/dgeometry2d/C4CGeometry.h"
#include "ImaGene/dgeometry2d/C4CSegmentPencil.h"
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/mathutils/Mathutils.h"
#include "ImaGene/mathutils/Statistics.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/dgeometry2d/C4CTangentialCover.h"
#include "ImaGene/helper/C4CTangentialCoverGeometry.h"
#include "ImaGene/dgeometry2d/FreemanChain.h"
#include "ImaGene/dgeometry2d/FreemanChainTransform.h"
#include "ImaGene/dgeometry2d/C4CIteratorOnFreemanChain.h"
#include "ImaGene/digitalnD/C4CIteratorOnFreemanChainSurface.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/C4CIteratorOnSurface.h"
#include "ImaGene/helper/ShapeHelper.h"
#include "ImaGene/helper/ScaleProfile.h"
#include "ImaGene/helper/MultiscaleProfile.h"
#include "ImaGene/helper/DrawingXFIG.h"
#include "ImaGene/helper/CurveVariationsHelper.h"
#include "ImaGene/helper/ContourHelper.h"

#include "ImaGene/mathutils/SimpleLinearRegression.h"
#include "ImaGene/dgeometry2d/BlurredSegmentTgtCover.h"






using namespace std;
using namespace ImaGene;

static Arguments args;
static const int RESOLUTION=1200;





///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
// M A I N
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char** argv ) 
{
  // -------------------------------------------------------------------------
  // Prepare arguments.
  //  StandardArguments::addIOArgs( args, true, false );

  args.addOption("-contourSet", "-contourSet <filename> ", " " );
  args.addOption("-sourceImage", "-sourceImage <filename> <LevelStep> (default: levelStep=10) ", " ", "10" );
  args.addBooleanOption( "-enteteXFIG", "-enteteXFIG: add xfig entete" );
  args.addOption( "-drawContourSet", "-drawContourSet <color> <linewidth> <depth>: draw the original source contour.",
		  "10", "1", "40" );
  args.addOption( "-meaningfulScales", "-meaningfulScales <min_size> <max_slope>: specifies parameters for defining meaningful scales: minimum size of the interval of scales and maximum slopes between consecutivesamples within.", "2", "0.0" );  
  args.addOption("-setSampling", "-setSampling <maxScale> <samplingStep>: set the maximal scale and sampling step used for contour analysis (default: maxScale=15.0, samplingStep=1.0)","15.0", "1.0");
  args.addOption("-drawMeaningfulContours", 
		 "-drawMeaningfulContours <MeaningfulTreshold> <color> <linewidth> <depth> draw the set of meaningful contours", "1.0", "4", "2", "0"); 
  
  args.addOption("-filterMeaningfulSlopes", "-filterMeaningfulSlopes  <flatTresholdMin>  <flatTresholdMax> filter flat meaningful contour according the slope of the linear regression defined on the meaningful intervall of the multiscale profile", "-10000", "-0.52" ); 
  
  
  
  args.addOption("-samplingDistance", "-samplingDistance <size> used for comparisons with [Cao2004] (default set as Nyquist distance=2 [Cao2004])","1");
  args.addOption("-saveResult", "-saveResult <file> export the epsilon MC in a single file with a segment represented at each line", "resultMS.dat");
  args.addBooleanOption("-invertVerticalAxis", "-invertVerticalAxis used to transform the contour representation (need for DGtal) (only used with option -saveResult)");

  args.addOption("-minSegmentSize", "-minSegmentSize <threshold> used for comparisons with [Cao2004] (the min size 10 seems good) ", "10" );
  
  args.addOption("-drawImage","-drawImage <imageName> <width> <height> draw image as background in XFIG format ", "", "", "");
  args.addBooleanOption("-displaySourcePoints", "-displaySourcePoints used when a large sampling distance is used in order to display the associated contour points " );
  
  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) || (!args.check("-contourSet") && !args.check("-sourceImage"))   )
    {

	cerr << args.usage( "displayMeaningfulContours", 
			    "displayMeaningfulContours"
			    ,"" ) << endl;
	return 1;
      }
  
     
  static int    agrandissementEPS = 100;
  
  if(args.check("-enteteXFIG")){
    cout << "#FIG 3.2 \nLandscape \nCenter \nInches \nLetter  \n" <<agrandissementEPS<< "\nSingle \n1" 
	 << " \n"<<RESOLUTION<<" 1" << endl;
  }
  
  uint mscales_min_size = args.getOption( "-meaningfulScales" )->getIntValue( 0 );
  double mscales_max_slope = args.getOption( "-meaningfulScales" )->getDoubleValue( 1 );
  
  

  vector<vector<Vector2D> > contourSet;
  BlurredSegmentTgtCover tgc;


  
  if(args.check("-contourSet")){
    string fileName = args.getOption("-contourSet")->getValue(0);
    ifstream ifs ( fileName.c_str() , ifstream::in );
    contourSet = ContourHelper::getPolygonsFromStream(ifs);
  }
  

  if(args.check("-sourceImage")){
    string fileName = args.getOption("-sourceImage")->getValue(0);
    uint stepLevelSet = args.getOption("-sourceImage")->getIntValue(1); 
    cerr << "Extracting iso contours with step: "<< stepLevelSet << " [ ";
    
    for(int i=0; i<255; i+=stepLevelSet){
      cerr << ".";
      ifstream ifs ( fileName.c_str() , ifstream::in );
      vector< vector<Vector2D> >  vSet = ContourHelper::getImageIsoContours(ifs, i,0, args.check("-invertVerticalAxis"));
      ifs.close();
      for(int k=0; k< vSet.size(); k++){
	contourSet.push_back(vSet.at(k));
      }
    }
    cerr << " ]" << endl;
  }
  cerr << "Number of contours:" <<contourSet.size() << endl;
  
  
  if(args.check("-drawContourSet")){
    uint color = args.getOption("-drawContourSet")->getIntValue(0);
    uint linewidth = args.getOption("-drawContourSet")->getIntValue(1);
    uint depth = args.getOption("-drawContourSet")->getIntValue(2);
    uint minSize= args.getOption("-minSegmentSize")->getIntValue(0);
    for(int i=0; i< contourSet.size(); i++){      
      if(contourSet.at(i).size()>minSize){
	DrawingXFIG::drawContour(cout, contourSet.at(i), color, color , linewidth, false, false, depth, 1.0, 0.0, 0.0);    
	
      }    
    }    
  }
  

  
  if(args.check("-drawImage")){
    string name = args.getOption("-drawImage")->getValue(0);
    uint width = args.getOption("-drawImage")->getIntValue(1);
    uint height = args.getOption("-drawImage")->getIntValue(2);
    DrawingXFIG::drawImage(cout, name, width, height, 50);
  }
  
  

  if(args.check("-drawMeaningfulContours")){
    double samplingSizeMax = args.getOption("-setSampling")->getFloatValue(0);    
    double samplingStep = args.getOption("-setSampling")->getFloatValue(1);      
    
    double meaningfulThreshold = args.getOption("-drawMeaningfulContours")->getFloatValue(0);
    uint color = args.getOption("-drawMeaningfulContours")->getIntValue(1);
    uint linewidth = args.getOption("-drawMeaningfulContours")->getIntValue(2);
    uint depth = args.getOption("-drawMeaningfulContours")->getIntValue(3);

    bool exportData = args.check("-saveResult");
    ofstream fileExport;
    if(exportData){
      string name= args.getOption("-saveResult")->getValue(0);
      fileExport.open(name.c_str());
      fileExport << "# Meaningful Segments (each line represent the segment (Xi Yi)" <<endl;
    }
    
    vector<double> vectScales;
    vector<double> vectNoiseReal;
    for(double i=samplingStep; i<samplingSizeMax; i=i+samplingStep){
      vectScales.push_back((double)i );
    }
    Statistics statMeaningfulParts(1, false);
    //Meaningfulscales detection 
    for(int i=0;i<contourSet.size(); i++){
      // Sampling contour according the author depends of Nyquist distance =2 pixel
      uint samplingDist= args.getOption("-samplingDistance")->getIntValue(0);
      vector<Vector2D> newContour;      
      for(int k =0; k< contourSet.at(i).size(); k++){
	if(k%samplingDist==0){
	  newContour.push_back(contourSet.at(i).at(k));
	}
      }
      
      Vector2D ptMin, ptMax;
      ContourHelper::computeContourBoundingBox (newContour, ptMin,ptMax);
      if((ptMax.x()-ptMin.x()<samplingSizeMax*3.0) &&
	 (ptMax.y()-ptMin.y()<samplingSizeMax*3.0))
	continue;
      BlurredSegmentTgtCover tgc;
      tgc.init(newContour,true);
      vector<double> noiseLevels;
      vector<double> slopeValues;
      double slopeThresholdMin = args.getOption("-filterMeaningfulSlopes")->getFloatValue(0);
      double slopeThresholdMax = args.getOption("-filterMeaningfulSlopes")->getFloatValue(1);
      

      if(args.check("-filterMeaningfulSlopes")){
	tgc.computeNoiseAndSlope(noiseLevels, slopeValues, vectScales, mscales_min_size, mscales_max_slope);
      }else{
	noiseLevels = tgc.getNoiseLevels(vectScales, mscales_min_size, mscales_max_slope);	
      }
      
      vector<double> vectNoiseReal;
      vector<uint> constraintIndex;
      vector<ContourHelper::DiskConstraint> vectConstraint;
      for(int j=0; j< newContour.size(); j++){
	Vector2D ptA = newContour.at(j);
	ContourHelper::DiskConstraint cc;
	double noiseLevelA;
	if(noiseLevels.at(j)==0){
	  noiseLevelA=vectScales.at(vectScales.size()-1)/2.0;
	}else{
	  noiseLevelA= (vectScales.at(noiseLevels.at(j)-1)/2.0);
	}
	cc.center=ptA;
	cc.radius=noiseLevelA;
	constraintIndex.push_back(vectConstraint.size());	
	vectConstraint.push_back(cc);       	
      }
      
      vector<vector<Vector2D> > contourSetFiltered = ContourHelper::filterMeaningfulParts(newContour,vectConstraint, 
											  constraintIndex, meaningfulThreshold,
											  slopeValues, slopeThresholdMin,  
											  slopeThresholdMax ); 
      
      
      
      uint minSize= args.getOption("-minSegmentSize")->getIntValue(0);
     
      
      for(int k=0; k< contourSetFiltered.size(); k++){
	ContourHelper::computeContourBoundingBox (contourSetFiltered.at(k), ptMin,ptMax);
	if(args.check("-minSegmentSize")&& contourSetFiltered.at(k).size()<minSize){
	  continue;
	}
	if((ptMax.x()-ptMin.x()<samplingSizeMax*3.0) &&
	  (ptMax.y()-ptMin.y()<samplingSizeMax*3.0))
	  continue;
	//if(contourSetFiltered.at(k).size()> minSize){
	statMeaningfulParts.addValue(0, contourSetFiltered.at(k).size());
	
	vector<Vector2D> setPtToDraw;
	if(args.check("-displaySourcePoints")){
	  setPtToDraw = ContourHelper::selectPolygonPart(contourSet.at(i), 
							 contourSetFiltered.at(k).at(0),
							 contourSetFiltered.at(k).at(contourSetFiltered.at(k).size()-1) );
	}else{
	  setPtToDraw = contourSetFiltered.at(k);
	}
	
	if(setPtToDraw.size()!=0){
	DrawingXFIG::drawContour(cout, setPtToDraw, 4,4, 1, false, false, 10, 1.0, 0.0, 0.0);    

	}
	if(exportData){
	  for(int l=0; l< setPtToDraw.size(); l++){
	    fileExport << setPtToDraw.at(l).x() << " " << setPtToDraw.at(l).y() << " ";
	  }	
	
	fileExport << endl;

	}
	
      }    

      
    }
    
    statMeaningfulParts.terminate();
  
    cerr << "Num Meaningful contours = " << statMeaningfulParts.samples(0)  << endl;
    cerr << "Mean size = " << statMeaningfulParts.mean(0)  << endl;
    cerr << "Max size = " << statMeaningfulParts.max(0)  << endl;
    cerr << "Min size = " << statMeaningfulParts.min(0)  << endl;
    cerr << "Variance size = " << statMeaningfulParts.variance(0)  << endl;
    if(exportData){
      fileExport.close();
    }

  }
  


}






