/////////////////////////////////////////////////////////////////////////////
// Display the Noise level obtained after multi scales analysis 
//////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
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
#include "ImaGene/helper/MultiscaleProfile.h"
#include "ImaGene/helper/DrawingXFIG.h"
#include "ImaGene/helper/CurveVariationsHelper.h"
#include "ImaGene/mathutils/SimpleLinearRegression.h"



using namespace std;
using namespace ImaGene;


static Arguments args;
void  extractField(string fileName, int indice, vector <double> &vectField );

int  getNumberOfLines( istream &inStream);







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
  
  
  args.addOption("-file", "-file <name> select the inmput file", " ");
  args.addOption("-extractAlignTgt", "-extractAlignTgt <posIndice> <posTgt> <threshold init Tgt> align the tangent values extract from <posIndice> <posTgt>  in order to avoid passage from 0 to 2PI (initial value put to I-2PI if less than <threshold init Tgt>  ", "1", "2", "6.0");
  args.addOption("-afficheXFIGTgt", "-afficheXFIGTgt <colX> <colY> <colTgt> <freq> draw the contour with tangent (of colums <col>) with a frequency <freq> " , "1","2","3", "5" );
  
  
  
  if ( ( argc <= 1 )|| ! args.readArguments( argc, argv )) 
    {
      cerr << args.usage( "test_Multiscale", 
			  "Analyse some multiscale properties of a shape with respect to its expected asymptotic properties"
			  ,"" ) << endl;
      return 1;
    }
  
  

  string fileName = args.getOption("-file")->getValue(0); 
  
  
  vector<double> vectIndices;
  vector<double> vectTangentes;
  
  
  if(args.check("-extractAlignTgt")){
    uint posIndice = args.getOption("-extractAlignTgt")->getIntValue(0);
    uint posTgt = args.getOption("-extractAlignTgt")->getIntValue(1);    
    double threshold = args.getOption("-extractAlignTgt")->getDoubleValue(2);
    
    extractField(fileName , posIndice,  vectIndices);
    extractField(fileName, posTgt,  vectTangentes);
    
    double currentTgt=vectTangentes.at(0);
    if(currentTgt>threshold){
      currentTgt-=2*M_PI;
    }
    
    vector<double> vectTgtAjusted;
    vectTgtAjusted.push_back(currentTgt);
    for(uint i =1; i < vectIndices.size(); i++){
      currentTgt+=Mathutils::AngleComputer::deviation(vectTangentes.at(i), currentTgt );
      vectTgtAjusted.push_back(currentTgt);
    }
    
    for(uint i =0; i < vectTgtAjusted.size(); i++){
      cout << i << " " << vectTgtAjusted.at(i) << endl;
    }    
  }
  

  


  if(args.check("-afficheXFIGTgt")){
    DrawingXFIG::includeXFIGHeader(cout, 1024, 2);
    
    uint posX = args.getOption("-afficheXFIGTgt")->getIntValue(0);
    uint posY = args.getOption("-afficheXFIGTgt")->getIntValue(1);    
    uint posTgt = args.getOption("-afficheXFIGTgt")->getIntValue(2);    
    int freq = args.getOption("-afficheXFIGTgt")->getIntValue(3);
    
    vector<double> contourX; 
    vector<double> contourY; 
    vector<double> contourTgt;
    vector<Vector2i> contour;
    
    extractField(fileName , posX,  contourX);
    extractField(fileName , posY,  contourY);
    extractField(fileName, posTgt, contourTgt);
    
    for (int i = 0; i< contourX.size(); i++){
      contour.push_back(Vector2i(contourX.at(i),contourY.at(i)));      
    }
      
    
    DrawingXFIG::drawContour (cout, contour, 0, 0, 20, true, false, 49, 1.0, 0.0,0.0);
    double tailleVector = 10.0;
    for(int i=0; i < contourX.size(); i++){
      if(i%freq==0){	
	double tangente= contourTgt.at(i);
	Vector2i pt1 = contour.at(i);
	Vector2i pt2 (pt1.x()+cos(tangente+M_PI/2.0)*tailleVector, pt1.y()+sin(tangente+M_PI/2.0)*tailleVector);
	Vector2i pt3 (pt1.x()+cos(tangente-M_PI/2.0)*tailleVector, pt1.y()+sin(tangente-M_PI/2.0)*tailleVector);
	DrawingXFIG::drawLine(cout, pt2, pt3, 1, 10, 0);
      }
    }
    
    
    
    

  }


}








  /**
   * Extract the field indice for all the lines of the file.
   * @param indice The initial inidice beginnig from 0.  
   * @return the tabular containing the extracted field   
   */


void 
extractField(string fileName, int indice, vector <double> &vectField ){
  fstream flux;
  flux.open (fileName.c_str(), ios::in);
  string str;
  getline(flux, str );
  int position=0;
  while ( flux.good() )
    {
      if ( ( str != "" ) 
	   && ( str[ 0 ] != '#' ) )
	{
	  istringstream in_str( str );
	  int idx = 0;
	  double val;
	  while ( in_str.good() )
	    {
	      in_str >> val;
	      if ( indice == idx ){
		vectField.push_back(val);
		position++;
		break;
		
	      }
	      ++idx;
	    }
	}
      getline(flux, str );      
    }
  
}



