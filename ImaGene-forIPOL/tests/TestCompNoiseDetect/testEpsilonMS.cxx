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
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/mathutils/Mathutils.h"
#include "ImaGene/mathutils/Statistics.h"
#include "ImaGene/timetools/Clock.h"

#include "ImaGene/helper/ShapeHelper.h"
#include "ImaGene/helper/DrawingXFIG.h"
#include "ImaGene/helper/ContourHelper.h"
#include "ImaGene/mathutils/Mathutils.h"

#include "ImaGene/helper/MeaningfulContinuation.h"

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
  StandardArguments::addIOArgs( args, true, false );


  args.addOption("-srcPolygon", "-srcPolygon <contour.sdp> <posX> <posY> <CLOSED|OPEN> (default CLOSED)", " ", "0", "1", "CLOSED" );  
  args.addBooleanOption( "-enteteXFIG", "-enteteXFIG: add xfig entete" );
  args.addOption( "-drawContourSRC", "-drawContourSRC <color> <linewidth> <depth>: draw the original source contour.",
		  "10", "30", "40" );
  args.addOption("-drawContourPoints", "-drawContourPoints <color> <pointSize> <depth>", "2", "1.0", "20"  );
  args.addOption("-drawEpsilonMS", "-drawEpsilonMS <epsilon> <Nc> ", "1", "100", "3");
  args.addOption("-samplingDistance", "-samplingDistance <size> (default set as Nyquist distance=2 [Cao2004])","2");
  args.addOption("-minSegmentSize", "-minSegmentSize <threshold> (default value=3 min size to a segment for being Meaningful continuation) ", "3" );


  
    if ( ( argc <= 0 ) 
	 || ! args.readArguments( argc, argv ) ||  !args.check("-srcPolygon"))
    {
      cerr << args.usage( "displayNoiseBS", 
			  "displayNoiseBS   -srcPolygon."
			  ,"" ) << endl;
      return 1;
    }
  
     
  static int    agrandissementEPS = 100;
  
  if(args.check("-enteteXFIG")){
    cout << "#FIG 3.2 \nLandscape \nCenter \nInches \nLetter  \n" <<agrandissementEPS<< "\nSingle \n1" 
	 << " \n"<<RESOLUTION<<" 1" << endl;
  }
  
  vector<Vector2D> polygon;
  
  uint minSegmentSize= args.getOption("-minSegmentSize")->getIntValue(0);
  cerr << "Using min segment size threshold=" << minSegmentSize << endl;
  
  
  if(args.check("-srcPolygon")){
    string fileName = args.getOption("-srcPolygon")->getValue(0);
    uint posX = args.getOption("-srcPolygon")->getIntValue(1);
    uint posY = args.getOption("-srcPolygon")->getIntValue(2);
    string closedOption =args.getOption("-srcPolygon")->getValue(3);
    ifstream ifs ( fileName.c_str() , ifstream::in );
    polygon =  ContourHelper::getPolygonFromStream(ifs, posX, posY);
  }

  
  
  if(args.check("-drawContourSRC")){
    uint color = args.getOption("-drawContourSRC")->getIntValue(0);
    uint linewidth = args.getOption("-drawContourSRC")->getIntValue(1);
    uint depth = args.getOption("-drawContourSRC")->getIntValue(2);
    DrawingXFIG::setFillIntensity(100);
    DrawingXFIG::drawContour(cout, polygon, color,color, linewidth, true, false, depth, 1.0, 0.0, 0.0);    
  }



 

  if(args.check("-drawContourPoints")){
    uint color = args.getOption("-drawContourPoints")->getIntValue(0);
    double  pointSize = args.getOption("-drawContourPoints")->getFloatValue(1);
    uint depth = args.getOption("-drawContourPoints")->getIntValue(2);
    for(int i=0; i< polygon.size(); i++){
      DrawingXFIG::drawCircle (cout, polygon.at(i), color, pointSize, depth, true, 0);
    }    
  }
  
  if(args.check("-samplingDistance")){
    uint samplingDist= args.getOption("-samplingDistance")->getIntValue(0);
     vector<Vector2D> newContour;      
     for(int i =0; i< polygon.size(); i++){
       if((i%samplingDist)==0){
	 newContour.push_back(polygon.at(i));
       }
     }
     polygon=newContour;
  }
  
  DrawingXFIG::drawCircle (cout, polygon.at(0), 2, 1, 0, true, 0);
  
  
  if(args.check("-drawEpsilonMS")){
    MeaningfulContinuation mc;
    double epsilon = args.getOption("-drawEpsilonMS")->getFloatValue(0);
    uint Nc = args.getOption("-drawEpsilonMS")->getIntValue(1);
    double lineWidth = args.getOption("-drawEpsilonMS")->getFloatValue(2);
    mc.init(polygon, epsilon, Nc, Mathutils::two_pi_d/4.0, false );

    vector<MeaningfulContinuation::EpsilonMeaningfulSegment> vectSeg= mc.getMeaningfulMaxSegments();
    
    for(int i=0; i< vectSeg.size(); i++){
      uint begin = vectSeg.at(i).back();
      uint end = vectSeg.at(i).front();
      vector<Vector2D> setPtToDraw;
      for(uint j=begin; j<= end; j++){
	setPtToDraw.push_back(polygon.at(j));
      } 
      if(setPtToDraw.size()>=minSegmentSize){
	DrawingXFIG::drawContour(cout, setPtToDraw, 4,4, lineWidth, false, false, 10, 1.0, 0.0, 0.0);    
      }
    }

    
  }
    

}















