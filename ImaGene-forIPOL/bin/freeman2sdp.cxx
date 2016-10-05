///////////////////////////////////////////////////////////////////////////////
// Generates contours from pgm image
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/dgeometry2d/FreemanChain.h"

using namespace std;
using namespace ImaGene;


static Arguments args;

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
  args.addBooleanOption( "-info", "-info: adds some info as comments at the beginning of the file." );
  args.addBooleanOption("-oneLine","-oneLine: output the digital contour in one line like: X0 Y0 X1 Y1 ... XN YN" );

  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "freeman2sdp", 
			  "Converts a contour defined as a Freeman chain code to a contour defined by a sequence of digital points. It writes them on the standard output.",
			  "" )
	   << endl;
      return 1;
    }
  FreemanChain c;
  FreemanChain::read( cin, c );
  bool oneLine = args.check("-oneLine");
  if ( cin.good() )
    {
      int min_x;
      int min_y;
      int max_x;
      int max_y;
      c.computeBoundingBox( min_x, min_y, max_x, max_y );
      if ( args.check( "-info" ) )
	cout << "# SDP contour" << endl
	     << "# size=" << c.chain.size() << endl
	     << "# bbox=" << min_x << " " << min_y
	     << " " << max_x << " " << max_y << endl;
      for ( FreemanChain::const_iterator it = c.begin();
	    it != c.end();
	    ++it )
	{
	  cout << (*it).x() << " " << (*it).y() ;
	  if(oneLine)
	    cout << " " ; 
	  else 
	    cout << endl;
	}
    }
  cout << endl;
  return 0;
}
