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
  StandardArguments::addIOArgs( args, true, false );

  args.addOption("-srcFC", "-srcFC <contour.fc>", " " );
  args.addOption("-srcPolygon", "-srcPolygon <contour.sdp> <posX> <posY> <CLOSED|OPEN> (default CLOSED)", " ", "0", "1", "CLOSED" );  
  args.addBooleanOption( "-enteteXFIG", "-enteteXFIG: add xfig entete" );
  args.addOption( "-drawContourSRC", "-drawContourSRC <color> <linewidth> <depth>: draw the original source contour.",
		  "10", "30", "40" );
  args.addOption( "-drawCntPointSRC", "-drawCntPointSRC <color> <linewidth> <depth>: draw the original source contour.",
		  "10", "30", "40" );
  args.addBooleanOption("-estimClosedContour", "-estimClosedContour estimate is the contour is closed or not" );
  
  args.addOption("-drawNoiseLevelBoxes", "-drawNoiseLevelBoxes <color>", "1");
  
  args.addOption("-drawContourPoints", "-drawContourPoints <color> <pointSize> <depth>", "2", "1.0", "20"  );
  args.addOption("-setSampling", "-setSampling <maxScale> <samplingStep>: set the maximal scale and sampling step used for contour analysis (default: maxScale=15.0, samplingStep=1.0)",
		 "15.0", "1.0");
  args.addOption("-displayTgtCover", "-displayTgtCover <width> <color> <linewidth>: display the set of tangential cover of width <width>  ", "2.0", "1","1"  );
  
  args.addOption("-selectPointTgtCover","-selectPointTgtCover <position>: (used with option displayTgtCover) draw the set of maximal segements which contains the point of index <position> ", "0");


  // def de bruit:
  args.addOption( "-meaningfulScales", "-meaningfulScales <min_size> <max_slope>: specifies parameters for defining meaningful scales: minimum size of the interval of scales and maximum slopes between consecutivesamples within.", "1", "-0.2" );  
  
  args.addOption("-setCstNoiseConstraint", "-setCstNoiseConstraint <diskSize>", "1" );
  args.addBooleanOption("-drawNoiseConstraints","-drawNoiseConstraints" );
  args.addOption("-drawPointSets", "-drawPointSets  <filename> <circleSize>", "", "0.5" );
  args.addBooleanOption("-auto", "-auto");
  args.addBooleanOption("-autoScale", "-autoScale set automaticly the scale interval by estimating the distance betwwen each vertex");

  args.addOption("-drawPolygonsSets", "-drawPolygonsSets <filename> draw a set of polygon (each line should represent the set of the coordinates", ""); 
  
  args.addOption("-extractScaleProfile", "extractScaleProfile <index> <file> extract Scale Profile of index point p ", "0", "profile.dat");

  args.addOption("-displaySmoothContour", "-displaySmoothContour <color> display the smoothed contour with the noise constraints", "0" );
  args.addOption("-exportNoiseLevel", "-exportNoiseLevel <filename> export noise level", "noise.dat" );



  // pour afficher l'image source.
  args.addOption("-setPosImage", "-setPosImage <x> <y>", "0", "0");
  args.addOption("-afficheImage", "-afficheImage <strinf nomImage> <largeur> <hauteur>","0", "0", "0");


    if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) || (!args.check("-srcFC") && !args.check("-srcPolygon"))  ) 
    {
      cerr << args.usage( "displayNoiseBS", 
			  "displayNoiseBS  -srcFC or -srcPolygon."
			  ,"" ) << endl;
      return 1;
    }
  
     
  static int    agrandissementEPS = 30;
  
  if(args.check("-enteteXFIG")){
    cout << "#FIG 3.2 \nLandscape \nCenter \nInches \nLetter  \n" <<agrandissementEPS<< "\nSingle \n1" 
	 << " \n"<<RESOLUTION<<" 1" << endl;
  }
  
  uint mscales_min_size = args.getOption( "-meaningfulScales" )->getIntValue( 0 );
  double mscales_max_slope = args.getOption( "-meaningfulScales" )->getDoubleValue( 1 );
  
  

  
  BlurredSegmentTgtCover tgc;

  
  // -------------------------------------------------------------------------
  // Read Freeman chain and creates space/contour.  

  bool exportNoise = args.check("-exportNoiseLevel");
  string nameExport = args.getOption("-exportNoiseLevel")->getValue(0);
  ofstream outNoiseLevel;
  if(exportNoise){
    outNoiseLevel.open( nameExport.c_str() , ofstream::out );
  }

  if(args.check("-srcFC")){
    FreemanChain fc; 
    string fileName = args.getOption("-srcFC")->getValue(0);
    ifstream ifs ( fileName.c_str() , ifstream::in );    
    FreemanChain::read( ifs, fc );
    if ( ! ifs.good() )
      {
	cerr << "Error reading Freeman chain code." << endl;
	return 2;
      }
    tgc.init(fc);
  }
  
  bool isOpenPolygon;
  if(args.check("-srcPolygon")){
    string fileName = args.getOption("-srcPolygon")->getValue(0);
    uint posX = args.getOption("-srcPolygon")->getIntValue(1);
    uint posY = args.getOption("-srcPolygon")->getIntValue(2);
    string closedOption =args.getOption("-srcPolygon")->getValue(3);
    ifstream ifs ( fileName.c_str() , ifstream::in );
    vector<Vector2D> polygon =  ContourHelper::getPolygonFromStream(ifs, posX, posY);
    bool loop = ContourHelper::containsLoop(polygon);
    if(loop) {
      cerr << "contains loop aborting..." << endl;
      return 0;
    }else{
      cerr << "no loop ok" << endl;
    }   
    if(args.check("-estimClosedContour")){
      isOpenPolygon = ContourHelper::isOpenPolygon(polygon,3.0);
    }else{
      isOpenPolygon=!(closedOption=="CLOSED");
    }
    cerr << "contour considered as :" << (isOpenPolygon ? "Open" :"Closed") << 
      "polygon size " << polygon.size() << endl;
    tgc.init(polygon, !isOpenPolygon);    
  }


  
  if(args.check("-drawPointSets")){
    string filename = args.getOption("-drawPointSets")->getValue(0);
    double pointSize = args.getOption("-drawPointSets")->getFloatValue(1);
    ifstream ifs ( filename.c_str() , ifstream::in );
    vector<Vector2D> pointSet =  ContourHelper::getPolygonFromStream(ifs, 0, 1);
    
    for(int i=0; i< pointSet.size(); i++){
      DrawingXFIG::drawCircle (cout, pointSet.at(i), 2, pointSize, 30, true, 0);
    }  
  }

  
  
  if(args.check("-drawContourSRC")){
    uint color = args.getOption("-drawContourSRC")->getIntValue(0);
    uint linewidth = args.getOption("-drawContourSRC")->getIntValue(1);
    uint depth = args.getOption("-drawContourSRC")->getIntValue(2);
    DrawingXFIG::setFillIntensity(100);
    
    DrawingXFIG::drawContour(cout, tgc.getPointsContour(), color,color, linewidth, 
    			     !isOpenPolygon, false, depth, 1.0, 0.0, 0.0);    


    cout << endl;
  }
  
  if(args.check("-drawCntPointSRC")){
    uint color = args.getOption("-drawCntPointSRC")->getIntValue(0);
    uint linewidth = args.getOption("-drawCntPointSRC")->getIntValue(1);
    uint depth = args.getOption("-drawCntPointSRC")->getIntValue(2);
    DrawingXFIG::setFillIntensity(100);
    
    vector<Vector2D> contour= tgc.getPointsContour(); 
    for(uint i=0; i< contour.size(); i++){
      DrawingXFIG::drawCircle (cout, contour.at(i), color, 0.2,depth);
    }
  }



  if(args.check("-afficheImage")){    
    string nomImage = args.getOption( "-afficheImage" )->getValue( 0 );
    int largeur = args.getOption( "-afficheImage" )->getIntValue( 1 );
    int hauteur = args.getOption( "-afficheImage" )->getIntValue( 2 );
    if(args.check("-setPosImage")){
      int posX = args.getOption("-setPosImage")->getIntValue(0);
      int posY = args.getOption("-setPosImage")->getIntValue(1);
      DrawingXFIG::drawImage(cout, nomImage, posX, posY, largeur, hauteur, 100);
    }else{
      DrawingXFIG::drawImage(cout, nomImage, 0, 0, largeur, hauteur, 100);
    }        
  }



  if(args.check("-drawPolygonsSets")){
    string fileName = args.getOption("-drawPolygonsSets")->getValue(0);
    cerr << "ok" << endl;
    ifstream ifs ( fileName.c_str() , ifstream::in );
    vector<vector<Vector2D> > contourSet = ContourHelper::getPolygonsFromStream(ifs);
    for(int i=0; i< contourSet.size(); i++){
      // if(contourSet.at(i).at(0).x()>0  && contourSet.at(i).at(0).y()>0){
      //DrawingXFIG::drawCircle (cout, contourSet.at(i).at(0), 2, 1, 30, true, 0);
	
	DrawingXFIG::drawContour(cout, contourSet.at(i), 0,0 , 1, false, false, 0, 1.0, 0.0, 0.0);    
	//}
    }    
  }


  if(args.check("-drawContourPoints")){
    uint color = args.getOption("-drawContourPoints")->getIntValue(0);
    double  pointSize = args.getOption("-drawContourPoints")->getFloatValue(1);
    uint depth = args.getOption("-drawContourPoints")->getIntValue(2);
    vector<Vector2D> contour = tgc.getPointsContour(); 
    for(int i=0; i< contour.size(); i++){
      DrawingXFIG::drawCircle (cout, contour.at(i), color, pointSize, depth, true, 0);
    }    
  }
  
  

  

  if(args.check("-displayTgtCover")){
    double width=  args.getOption("-displayTgtCover")->getFloatValue(0);
    uint lineColor = args.getOption("-displayTgtCover")->getIntValue(1);
    uint lineWidth = args.getOption("-displayTgtCover")->getIntValue(2);
    uint posSelected = args.getOption("-selectPointTgtCover")->getIntValue(0);
    if(args.check("-selectPointTgtCover")){
      Vector2D pt = tgc.getPointsContour().at(posSelected); 
      DrawingXFIG::drawCircle (cout, Vector2i(pt.x(), pt.y()), 4, 0.5, 0, true, 1);
    }

    tgc.ExtractSegment( width, 0, 0 );
    vector<ImaGene::BlurredSegmentTgtCover::Segment> segmentVect = tgc.getSegmentContour();  
    for(int i=0; i< segmentVect.size(); i++){
      ImaGene::BlurredSegmentTgtCover::Segment seg = segmentVect.at(i); 
      tgc.computeRealBounds(seg, width);
      cerr << "segment num " << i << " P= (" << seg.ptP.x()<< ", "<< seg.ptP.y() << " )" 
	   << " Q= (" << seg.ptQ.x()<< ", "<< seg.ptQ.y() << " )"
	   << " S= (" << seg.ptS.x()<< ", "<< seg.ptS.y() << " )" << "Nb point :"<<  
	seg.nbElementsAfterInitialisation<< " " << seg.nbElements<<   endl;
      if(!args.check("-selectPointTgtCover") || (seg.first()<=posSelected && seg.last() >= posSelected)
	 || (seg.last()<seg.first() && (posSelected> seg.first()||posSelected< seg.last()))){
	DrawingXFIG::drawBlurredSegment(cout, seg, lineColor, lineWidth, 10);
      }
    }       
  }
  
  


  
  double samplingSizeMax = args.getOption("-setSampling")->getFloatValue(0);    
  double samplingStep = args.getOption("-setSampling")->getFloatValue(1);    
  
  if(args.check("-auto") || args.check("-autoScale")){
    double meanVertexDistance= ContourHelper::getMeanVertexDistances(tgc.getPointsContour());
    cerr << "Mean distance between vertex: "  <<  meanVertexDistance << endl;
    samplingStep=meanVertexDistance;
  }
  
  
  
  
  
  
  cerr << "Computing Meaningful Thickness" << endl;
  vector<ContourHelper::DiskConstraint> vectConstraint;


  
  vector<double> vectScales;

  for(double i=1.0; i<samplingSizeMax; i=i+samplingStep){
    vectScales.push_back((double)i );
  }

  Statistics statNoise(1, false); 
  
  vector<uint> constraintIndex;
  Clock::startClock();
  vector<double> noiseLevels= tgc.getNoiseLevels(vectScales, 1, -0.0);
  long time = Clock::stopClock();
  cerr << "Meaningful Thickness computed in " << time << " ms."<<  endl;
  
  double nbSampleFill=0.1;
  double nbSampleFillSingle=0.1;
  double nbSampleFillFull=3;
  vector<Vector2D> contour = tgc.getPointsContour(); 
  DrawingXFIG::setFillIntensity(100);

  
  if(args.check("-extractScaleProfile")){
    uint pos = args.getOption("-extractScaleProfile")->getIntValue(0);
    string filename = args.getOption("-extractScaleProfile")->getValue(1);
    ScaleProfile sc;
    tgc.getScaleProfile(vectScales, pos, sc);
    ofstream ofs ( filename.c_str() , ifstream::out );
    vector<double> vectX;
    vector<double> vectY;
    sc.getProfile (vectX, vectY);
    ofs << "# Scale profile generated by displayNoiseBS on index point: " << pos << endl;    
    for(int i=0; i< vectX.size(); i++){
      ofs << exp(vectX.at(i)) << " " << exp(vectY.at(i)) << endl;
    }

    DrawingXFIG::drawCircle (cout, contour.at(pos), 4, 1, 0, true, 0);
  }



  double maxNoiseValue;
  if(args.check("-auto")){    
    for(int i=0; i< contour.size(); i++){
      Vector2D ptA = contour.at(i);
      int noiseAindex = (noiseLevels.at(i)==0)? (vectScales.size()-1): noiseLevels.at(i) ;
      double noiseLevelA= (vectScales.at(noiseAindex-1)/2.0);
      statNoise.addValue(0, noiseLevelA);
    }
    statNoise.terminate();
    maxNoiseValue=statNoise.max(0);
    cerr << "Mean noise=" << maxNoiseValue << endl; 
  }
  
  if(exportNoise){
    outNoiseLevel << "# Noise levels exported from displayNoiseBS (source ImaGene): " << std::endl;
    outNoiseLevel << "# Format X Y noiseLevel " << std::endl;
  }
     
  
  for(int i=0; i< contour.size(); i++){
     Vector2D ptA = contour.at(i);
     if(isOpenPolygon && (i+1==contour.size()) )
       break;
     
     Vector2D ptB = contour.at((i+1)%contour.size());

     
     // (index=0 si aucune échelle => bruit MAX)
     int noiseAindex = (noiseLevels.at(i)==0)? (vectScales.size()-1): noiseLevels.at(i) ;
     int noiseBindex = (noiseLevels.at((i+1)%contour.size())==0)?(vectScales.size()-1): 
       noiseLevels.at((i+1)%contour.size())  ;

     if(!(noiseBindex==1 && noiseAindex ==1)){
       nbSampleFill=nbSampleFillFull;
     }else{
       nbSampleFill=nbSampleFillSingle;
     }

     double noiseLevelA= (vectScales.at(noiseAindex-1)/2.0);
     double noiseLevelB= (vectScales.at(noiseBindex-1)/2.0);

     if(exportNoise){
       outNoiseLevel << ptA.x() << " " << ptA.y() << " "<<  noiseLevelA << std::endl;
     }
   
     if(args.check("-drawNoiseLevelBoxes")){
       //uint color = args.getOption("-drawNoiseLevelBoxes")->getIntValue(0);
       if(noiseLevelA>=0){
	 // DrawingXFIG::drawCircle (cout,ptA, 1, noiseLevelA,50,true, 0);
	 DrawingXFIG::setFillIntensity(38);
	 DrawingXFIG::drawPixel (cout, ptA, 1,1 ,  noiseLevelA*2.0 , 50);
       }
     }

  
     if(args.check("-setCstNoiseConstraint") ){
       double diskSize= args.getOption("-setCstNoiseConstraint")->getFloatValue(0);
       noiseLevelA=diskSize;
       noiseLevelB=diskSize;
     }else if(args.check("-auto")){
       noiseLevelA=maxNoiseValue;
       noiseLevelB=maxNoiseValue;
     }  
     
     
     

     double distance = VectorUtils::norm(Vector2D(ptB.x()- ptA.x(), ptB.y()- ptA.y())); 
     uint nbAdd = distance*nbSampleFill;
     if(args.check("-drawNoiseConstraints")){
       DrawingXFIG::drawCircle (cout,ptA, 1, noiseLevelA,50,true, 0);
     }


     ContourHelper::DiskConstraint cc;
     cc.center=ptA;
     cc.radius=noiseLevelA;
     vectConstraint.push_back(cc);
     constraintIndex.push_back(vectConstraint.size());

     for (int k=0; k<nbAdd; k++){
       Vector2D ptNew((ptA.x()/(nbAdd+1.0))*(k+1)+(ptB.x()/(nbAdd+1.0))*(nbAdd-k),
		      (ptA.y()/(nbAdd+1.0))*(k+1)+(ptB.y()/(nbAdd+1.0))*(nbAdd-k) );
       if(args.check("-drawNoiseConstraints")){
	 DrawingXFIG::drawCircle (cout,ptNew, 1,(noiseLevelA/(nbAdd+1.0))*(k+1)+(noiseLevelB/(nbAdd+1.0))*(nbAdd-k),50,true, 0);
       }
       ContourHelper::DiskConstraint circC;
       circC.center=ptNew;
       circC.radius=(noiseLevelA/(nbAdd+1.0))*(k+1)+(noiseLevelB/(nbAdd+1.0))*(nbAdd-k);
       vectConstraint.push_back(circC);
       
     }
     if(args.check("-drawNoiseConstraints")){
       DrawingXFIG::drawCircle (cout,ptB, 1, noiseLevelB,50,true, 0);
     }



     
  }
  if(exportNoise){
    outNoiseLevel.close();
  }
       
  
  if(args.check("-displaySmoothContour")){
      Clock::startClock();
      uint color = args.getOption("-displaySmoothContour")->getIntValue(0);
      vector<Vector2D> contourSmooth = ContourHelper::smoothContourFromNoise (contour,vectConstraint, 
									      constraintIndex,0.0001,0.25,
									      1.0, isOpenPolygon );
      long time = Clock::stopClock();
      cerr << "Noise removed  in " << time << " ms."<<  endl;
      
      
      DrawingXFIG::drawContour(cout, contourSmooth, color,color, 1, 
			       !isOpenPolygon, false, 0, 1.0, 0.0, 0.0);    
      
  }
  
  
  statNoise.terminate ();
     
  
    

}















