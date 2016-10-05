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
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/base/Arguments.h"
#include "ImaGene/mathutils/Mathutils.h"
#include "ImaGene/mathutils/Statistics.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/helper/ShapeHelper.h"
#include "ImaGene/helper/DrawingXFIG.h"
#include "ImaGene/helper/ContourHelper.h"

#include "ImaGene/helper/MeaningfulContinuation.h"



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
  
  args.addBooleanOption( "-enteteXFIG", "-enteteXFIG: add xfig entete" );
  args.addOption("-contourSet", "-contourSet <filename> init data from a file containing all the contour set ", " " );
  args.addOption("-sourceImage", "-sourceImage <filename> <LevelStep> (default: levelStep=10) ", " ", "10" );
  
  args.addOption( "-drawContourSet", "-drawContourSet <color> <linewidth> <depth>: draw the original source contour.",
		  "10", "1", "40" );
  args.addOption("-drawEpsilonMC", 
		 "-drawEpsilonMC <epsilon> <color> <linewidth>  draw the set of meaningful contours", "1.0","4", "2"); 
  args.addOption("-samplingDistance", "-samplingDistance <size> (default set as Nyquist distance=2 [Cao2004])","2");
  args.addOption("-minSegmentSize", "-minSegmentSize <threshold> (the min size 10 seems good) ", "10" );
  args.addOption("-saveResult", "-saveResult <file> export the epsilon MC in a single file with a segment represented at each line", "resultMS.dat");
  args.addBooleanOption("-displaySourcePoints", "-displaySourcePoints used when a large sampling distance is used in order to display the associated contour points " );

  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) || (!args.check("-contourSet") && !args.check("-sourceImage"))   )
    {
      cerr << args.usage( "displayEpsilonMC", 
			  "displayEpsilonMC"
			  ,"" ) << endl;
      return 1;
    }
  
     
  static int    agrandissementEPS = 100;
  
  if(args.check("-enteteXFIG")){
    cout << "#FIG 3.2 \nLandscape \nCenter \nInches \nLetter  \n" <<agrandissementEPS<< "\nSingle \n1" 
	 << " \n"<<RESOLUTION<<" 1" << endl;
  }
  
  
  vector<vector<Vector2D> > contourSet;

  
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
      vector< vector<Vector2D> >  vSet = ContourHelper::getImageIsoContours(ifs, i);
      ifs.close();
      for(int k=0; k< vSet.size(); k++){
	contourSet.push_back(vSet.at(k));
      }
    }
    cerr << " ]" << endl;
  }

  
  if(args.check("-drawContourSet")){
    uint color = args.getOption("-drawContourSet")->getIntValue(0);
    uint linewidth = args.getOption("-drawContourSet")->getIntValue(1);
    uint depth = args.getOption("-drawContourSet")->getIntValue(2);
    for(int i=0; i< contourSet.size(); i++){      
      if(i%2==0){
	DrawingXFIG::drawContour(cout, contourSet.at(i), color, color , linewidth, false, false, depth, 1.0, 0.0, 0.0);    
      }
    }
  }
  
  
  cerr << "Number of iso contours = " << contourSet.size() << endl;
  

  if(args.check("-drawEpsilonMC")){
    double epsilon = args.getOption("-drawEpsilonMC")->getFloatValue(0);
    
    uint color = args.getOption("-drawEpsilonMC")->getIntValue(2);
    uint linewidth = args.getOption("-drawEpsilonMC")->getIntValue(3);
    uint depth = 0;
    Statistics statGoodContinuation(1, false);
    
    uint minSegmentSize= args.getOption("-minSegmentSize")->getIntValue(0);
    cerr << "Using min segment size threshold=" << minSegmentSize << endl;
    
    bool exportData = args.check("-saveResult");
    
    ofstream fileExport;
    if(exportData){
      string name= args.getOption("-saveResult")->getValue(0);
      fileExport.open(name.c_str());
      fileExport << "# Meaningful Good continuation segments (each line represent the segment (Xi Yi)" <<endl;
    }
    
    
    //Meaningfulscales detection 
    for(int k=0; k<contourSet.size(); k++){
      
      // Sampling contour according the author depends of Nyquist distance =2 pixel
      uint samplingDist= args.getOption("-samplingDistance")->getIntValue(0);
      vector<Vector2D> newContour;      
      for(int i =0; i< contourSet.at(k).size(); i++){
	if((i%samplingDist)==0){
	  newContour.push_back(contourSet.at(k).at(i));
	}
      }
      MeaningfulContinuation mc;
      uint Nc = contourSet.size();
      mc.init(newContour, epsilon, Nc, Mathutils::two_pi_d/4.0, false );
      vector<MeaningfulContinuation::EpsilonMeaningfulSegment> vectSeg= mc.getMeaningfulMaxSegments();
      for(int i=0; i< vectSeg.size(); i++){
	if(vectSeg.at(i).size()>minSegmentSize){
	  uint begin = vectSeg.at(i).back();
	  uint end = vectSeg.at(i).front();
	  statGoodContinuation.addValue(0, end-begin);
	  vector<Vector2D> setPtToDraw;
	  if(args.check("-displaySourcePoints")){
	    Vector2D ptA = newContour.at(begin);
	    Vector2D ptB = newContour.at(end);
	    setPtToDraw = ContourHelper::selectPolygonPart(contourSet.at(k), ptA, ptB);
	    for(int m=0; m< setPtToDraw.size(); m++){
	      fileExport << setPtToDraw.at(m).x() << " " << setPtToDraw.at(m).y() << " ";
	    }
	  }else{
	    for(uint j=begin; j!= end; j++){
	      setPtToDraw.push_back(newContour.at(j));
	      if(exportData){
		fileExport << newContour.at(j).x() << " " << newContour.at(j).y() << " ";
	      }
	    } 
	  }
	  if(exportData){
	    fileExport << endl;
	  }
	  DrawingXFIG::drawContour(cout, setPtToDraw, 4,4, 3, false, false, 10, 1.0, 0.0, 0.0);    
	  
	}
      }
    }      
     
    
    statGoodContinuation.terminate();
  
    cerr << "Num Meaningful good continuation = " << statGoodContinuation.samples(0)  << endl;
    cerr << "Mean size = " << statGoodContinuation.mean(0)  << endl;
    cerr << "Max size = " << statGoodContinuation.max(0)  << endl;
    cerr << "Min size = " << statGoodContinuation.min(0)  << endl;
    cerr << "Variance size = " << statGoodContinuation.variance(0)  << endl;
    
    if(exportData){
      fileExport.close();
    }
    
  }
  

  


}
