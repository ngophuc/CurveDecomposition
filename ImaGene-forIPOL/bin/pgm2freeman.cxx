///////////////////////////////////////////////////////////////////////////////
// Generates contours from pgm image
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream> 

#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/Proxy.h"
#include "ImaGene/base/Vector.h"
#include "ImaGene/base/Shapes.h"
#include "ImaGene/digitalnD/GridEmbedder.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/KnRCellSet.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/C4CIteratorOnBdry.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/ObjectBoundary.h"
#include "ImaGene/digitalnD/Frame2D.h"
#include "ImaGene/helper/ContourHelper.h"
#include "ImaGene/helper/ShapeHelper.h"

using namespace std;
using namespace ImaGene;


static Arguments args;


/*
 * Return the histogram associated to the 2D image 
 */
std::vector<unsigned int>
importWithHistoFromPGM
( std::istream & in,
  std::iostream & inImage,
  KnSpace* & ks,
  uint bordXY, unsigned int &imageSize )
{

  bool quiet=true;
  ks = 0;
  string str;
  getline( in, str );
  if ( ! in.good() ) return std::vector<unsigned int>();
  if ( str != "P5" ) return std::vector<unsigned int>();
  inImage << "P5" << endl;
  do
    {
      getline( in, str );
      if ( ! in.good() ) return std::vector<unsigned int>();
    }
  while ( str[ 0 ] == '#' || str=="");
  istringstream str_in( str );
  Kn_size sizes[ 2 ];
  str_in >> sizes[ 0 ] >> sizes[ 1 ];
 
  uint xMaxObj=sizes[0];
  uint yMaxObj=sizes[1];

  imageSize = sizes[0]*sizes[1];
  uint sizeInit0 = sizes[0]; 
  uint sizeInit1 = sizes[1]; 

  sizes[0]+=2*bordXY;
  sizes[1]+=2*bordXY;
    
  getline( in, str );
  istringstream str2_in( str );
  int max_value;
  str2_in >> max_value;
  inImage << sizeInit0 << " " << sizeInit1 << endl << max_value << endl;


  std::vector<unsigned int> histo(max_value+1, 0);
  if ( ! in.good() ) return std::vector<unsigned int>();
  if(!quiet){
    cerr << "# PGM " << sizes[ 0 ] << " " << sizes[ 1 ] 
	 << " " << max_value << " from <" << str << ">" << endl;
  }
  ks = new KnSpace( 2, sizes );
  if ( ks == 0 ) return std::vector<unsigned int>();
  	
  KnSpaceScanner scan2( *ks, 
 			ks->ufirstCell( ks->dim() ),
			ks->ulastCell( ks->dim() ) );
  uint nb_read = 0;
  Kn_uid last_y, last_x;
  Kn_uid p = scan2.begin();
  uint x,y;
  
  in >> noskipws;
  for ( y=0, last_y = scan2.last( p, 1 );
	p <= last_y; y++,
	  p += scan2.gotonext( 1 ) )
    {
      for ( x=0, last_x = scan2.last( p, 0 ); 
	    p <= last_x; x++,
	      p++ ) // NB: 'scan.gotonext( 0 )' == 1;
	{ //... whatever
	  if(!((x<bordXY) ||( y<bordXY)  || 
	       (x>=(xMaxObj+bordXY))|| (y>=(yMaxObj+bordXY)) )){
	    unsigned char c; 
	    in >> c;
            inImage <<(unsigned char) c;
            histo[c]++;
	    if ( in.good() ) ++nb_read;	    
	  }
	}
    }
  if ( in.fail() || in.bad() )
    {
      cerr << "# nbread=" << nb_read << endl;
      delete ks;
      ks = 0;
      return std::vector<unsigned int>();
    }
  in >> skipws;
  inImage <<endl;
  inImage.flush();
  return histo;
}

unsigned int 
getThreshold(std::vector<unsigned int> histo, unsigned int imageSize){
  unsigned int sumA = 0;
  unsigned int sumB = imageSize;
  unsigned int muA=0;
  unsigned int muB=0;
  unsigned int sumMuAll= 0;
  for( unsigned int t=0; t< histo.size();t++){
    sumMuAll+=histo[t]*t;
  }
  
  unsigned int thresholdRes=0;
  double valMax=0.0;
  for( unsigned int t=0; t< histo.size(); t++){
    sumA+=histo[t];
    if(sumA==0)
      continue; 
    sumB=imageSize-sumA;
    if(sumB==0){
      break;
    }
    
    muA+=histo[t]*t;
    muB=sumMuAll-muA;
    double muAr=muA/(double)sumA;
    double muBr=muB/(double)sumB;
    double sigma=  (double)sumA*(double)sumB*(muAr-muBr)*(muAr-muBr);
    if(valMax<=sigma){
      valMax=sigma;
      thresholdRes=t;
    }
  }
  return thresholdRes;
}


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

  args.addOption( "-threshold", "-threshold <val>: threshold value for binarizing PGM gray values (def. is 128).", "128" );
  args.addOption( "-min_size", "-min_size <m>: minimum digital length of contours for output (def. is 4).", "4" );
  args.addOption( "-badj", "-badj <0/1>: 0 is interior bel adjacency, 1 is exterior (def. is 0).", "0" );
  args.addOption("-selectContour", "-selectContour <x0> <y0> <distanceMax>: select the contours for which the first point is near (x0, y0) with a distance less than <distanceMax>","0", "0", "0" );
  args.addBooleanOption("-invertVerticalAxis", "-invertVerticalAxis used to transform the contour representation (need for DGtal), used only for the contour displayed, not for the contour selection (-selectContour). ");

  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "pgm2freeman", 
			  "Extracts all 2D contours from a PGM image given on the standard input and writes them on the standard output as FreemanChain's. By default the threshold is computed with the Otsu algorithm but you can also set it manually (using -threshold)",
			  "" )
	   << endl;
      return 1;
    }
  
  bool yInverted = args.check("-invertVerticalAxis");
  KnSpace* ks;
  KnCharSet* voxset;
  uint threshold = (uint) args.getOption( "-threshold" )->getIntValue( 0 );
  uint imageSize = 0;
  

  if(!args.check("-threshold")){
    std::fstream fs;
    fs.open ("tmpPgm2Freeman.pgm", std::ofstream::out | std::ofstream::app);
    std::vector<unsigned int> histo= importWithHistoFromPGM( cin, fs,  ks, 1, imageSize ); 
    threshold = getThreshold(histo, imageSize);
    cerr << "Automatic threshold with otsu=" <<  threshold << endl;
    if ( histo.size()==0 )
      {
	cerr << "Error reading PGM file histo." << endl;
	return 2;
      }
    fs.close();
    fs.open ("tmpPgm2Freeman.pgm", std::ofstream::in);
    if ( ! ShapeHelper::importFromPGM( fs, ks, voxset, threshold, 1 , true) )
      {
	cerr << "Error reading PGM file." << endl;
	return 2;
      } 
    fs.close();
    remove( "tmpPgm2Freeman.pgm" );
  }else{
    if ( ! ShapeHelper::importFromPGM( cin, ks, voxset, threshold, 1, true) )
      {
	cerr << "Error reading PGM file." << endl;
	return 2;
      } 
    
  }

 
  Vector2i ptReference;
  double distanceMax=0.0;
  //Rajout (BK) 
  if(args.check("-selectContour")){
    ptReference.x()= args.getOption("-selectContour")->getIntValue(0);
    ptReference.y()= args.getOption("-selectContour")->getIntValue(1);
    distanceMax= args.getOption("-selectContour")->getIntValue(2);
  }

  bool interior = args.getOption( "-badj" )->getIntValue( 0 ) == 0;
  uint min_size = args.getOption( "-min_size" )->getIntValue( 0 );
  BelAdjacency badj( *ks, interior );
  KnRCellSet bdry = KnShapes::smakeBoundary( *ks, *voxset );
  KnRCellSet not_visited( bdry );
  uint num_contour = 0;
  for ( KnRCellSet::cell_iterator cell_it = bdry.begin();
	cell_it != bdry.end();
	++cell_it )
    {
      Kn_sid bel = *cell_it;
      uint k = *( ks->sbegin_dirs( bel ) );
      C4CIteratorOnBdry c4c_it( badj, bel, k, *voxset );
      bool is_open;
      uint nb_surfels = C4CIterator::size( c4c_it, is_open );
      if ( nb_surfels >= min_size )
	{
	  Proxy<C4CIteratorOnSurface> cp
	    ( (C4CIteratorOnSurface*) c4c_it.clone() );
	  if(!args.check("-selectContour")){
	    ContourHelper::displayFreemanChain( cout, ks, cp, 0, 1, yInverted );
	  }else{
	    //Rajout option (BK 29/07/09)
	    Frame2D frame;
	    frame.init( ks, 0, 1 );
	    Kn_sid sbel = cp->current();
	    frame.setSurfelFrame( sbel, cp->trackDir() );
	    Vector2i p1( frame.transformPoint( Vector2i( 0, 0 ) ) );
	    double distance = sqrt((p1.x() - ptReference.x())*(p1.x() - ptReference.x())+
				   (p1.y() - ptReference.y())*(p1.y() - ptReference.y()));
	    if(distance< distanceMax){
	      ContourHelper::displayFreemanChain( cout, ks, cp, 0, 1, yInverted );
	    }
	  }
	}
      // Clear contour from set of bels.
      bel = c4c_it.current();
      Kn_sid sbel = bel;
      do
	{
	  bdry[ bel ] = false;
	  if ( c4c_it.next() == 0 ) break;
	  bel = c4c_it.current();
	}
      while ( bel != sbel );

      num_contour++;
    }



  return 0;
}
  
  
