///////////////////////////////////////////////////////////////////////////////
// Generates contours from pgm image
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/mathutils/Statistics.h"

using namespace std;
using namespace ImaGene;


static Arguments args;


/**
 * Reads the input stream [in] line by line, excluding line beginning
 * with '#' and selects some columns according to [indices] to fill
 * the array of samples [samples]. For instance, if indices == { 3, 4,
 * 1 }, then samples( 0 ) will contain the third column, samples( 1 ),
 * the fourth column and samples( 2 ) the first column.
 *
 * @param in (modified) the input stream.
 * @param samples (updates) stores the data.
 * @param indices specifies in which columns data are read.
 */
void readShapeGeometry( istream & in, Statistics & samples,
			const vector<uint> & indices )
{
  string str;
  getline( in, str );
  while ( in.good() )
    {
      if ( ( str != "" ) 
	   && ( str[ 0 ] != '#' ) )
	{
	  istringstream in_str( str );
	  uint idx = 1;
	  double val;
	  while ( in_str.good() )
	    {
	      in_str >> val;
	      for ( uint j = 0; j < indices.size(); ++j )
		if ( indices[ j ] == idx )
		  {
		    samples.addValue( j, val );
		    // cout << "Adding " << val << " to " << j << endl;
		  }
	      ++idx;
	    }
	}
      getline( in, str );
    }
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
  args.addBooleanOption( "-centering", "-centering: centers the shape so that its center of gravity is (0,0)." );
  args.addBooleanOption( "-close", "-close: close the shape by spreading the last difference vector between all the points." );
  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "sdp2pts", 
			  "Reads a sequence of digital points on the standrad inputand converts it into a sequence of points with some optional transformations.",
			  "" )
	   << endl;
      return 1;
    }

  // Force storing of samples.
  Statistics stats( 3, true );
  vector<uint> indices( 2 );
  indices[ 0 ] = 1;
  indices[ 1 ] = 2;
  readShapeGeometry( cin, stats, indices );

  uint n = stats.samples( 0 );
  vector<double> x( n );
  vector<double> y( n );
  for ( uint i = 0; i < n; ++i )
    {
      x[ i ] = stats.value( 0, i );
      y[ i ] = stats.value( 1, i );
    }

  if ( args.check( "-close" ) )
    {
      double x0 = ( x[ 0 ] - x[ n - 1 ] );
      double y0 = ( y[ 0 ] - y[ n - 1 ] );
      for ( uint i = 0; i < n; ++i )
	{
	  x[ i ] += ( x0 * (double) i ) / n;
	  y[ i ] += ( y0 * (double) i ) / n;
	}
    }

  if ( args.check( "-centering" ) )
    {
      double x0 = 0.0;
      double y0 = 0.0;
      for ( uint i = 0; i < n; ++i )
	{
	  x0 += x[ i ];
	  y0 += y[ i ];
	}
      x0 /= n;  
      y0 /= n;
      for ( uint i = 0; i < n; ++i )
	{
	  x[ i ] -= x0;
	  y[ i ] -= y0;
	}
    }
  for ( uint i = 0; i < n; ++i )
    cout << x[ i ] << " " << y[ i ] << endl;
  return 0;
}
