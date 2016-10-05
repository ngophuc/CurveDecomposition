///////////////////////////////////////////////////////////////////////////////
// Generates a shape (set of points) from curvature information
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/Proxy.h"
#include "ImaGene/base/Vector.h"
#include "ImaGene/base/Shapes.h"
#include "ImaGene/mathutils/Statistics.h"

// #include "ImaGene/dgeometry2d/FreemanChain.h"
// #include "ImaGene/dgeometry2d/C4CIteratorOnFreemanChain.h"
// #include "ImaGene/digitalnD/GridEmbedder.h"
// #include "ImaGene/digitalnD/KnSpace.h"
// #include "ImaGene/digitalnD/KnSpaceScanner.h"
// #include "ImaGene/digitalnD/KnCharSet.h"
// #include "ImaGene/digitalnD/KnRCellSet.h"
// #include "ImaGene/digitalnD/KnShapes.h"
// #include "ImaGene/digitalnD/C4CIteratorOnBdry.h"
// #include "ImaGene/digitalnD/C4CIteratorOnFreemanChainSurface.h"
// #include "ImaGene/digitalnD/BelAdjacency.h"
// #include "ImaGene/digitalnD/ObjectBoundary.h"
// #include "ImaGene/digitalnD/Frame2D.h"
// #include "ImaGene/helper/ContourHelper.h"
// #include "ImaGene/helper/ShapeHelper.h"

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

void integrateShapeFromCurvatureAndCAbs
( vector<double> & x,
  vector<double> & y,
  vector<double> & theta,
  const Statistics & samples,
  uint curv_idx,
  uint length_idx )
{
  double x0 = x[ 0 ];
  double y0 = y[ 0 ];
  double theta0 = theta[ 0 ];
  double plength = 0.0;
  for ( uint i = 0; i < samples.samples( 0 ); )
    {
      double curv = samples.value( curv_idx, i );
      double nlength = samples.value( length_idx, i );
      double ds = nlength - plength;
      // cout << "k=" << curv << " ds=" << ds << endl;
      double theta1 = theta0 + ds * curv;
      if ( curv != 0.0 )
	{
	  x0 += ( sin( theta1 ) - sin( theta0 ) ) / curv;
	  y0 += ( cos( theta0 ) - cos( theta1 ) ) / curv;
	}
      else
	{
	  x0 += ds * cos( theta0 );
	  y0 += ds * sin( theta0 );
	}
      ++i;
      x[ i ] = x0;
      y[ i ] = y0;
      theta[ i ] = theta1;
      plength = nlength;
      theta0 = theta1;
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
  args.addOption( "-curv_idx", "-curv_idx: specifies in which column is the curvature value.", "1" );
  args.addOption( "-cabs_idx", "-cabs_idx: specifies in which column is the curvilinear abscissa value.", "2" );
  args.addOption( "-p0", "-p0 <x0> <y0>: specifies starting point.", "0.0", "0.0" );
  args.addOption( "-t0", "-t0 <theta0>: specifies starting angle to x-axis in radian.", "0.0" );
  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "integrate_shape", 
			  "Create a shape approximation as a set of points from curvature and length/tangent information.",
			  "" )
	   << endl;
      return 1;
    }
  
  uint curv_idx = args.getOption( "-curv_idx" )->getIntValue( 0 );
  uint cabs_idx = args.getOption( "-cabs_idx" )->getIntValue( 0 );
  // Force storing of samples.
  Statistics stats( 3, true );
  vector<uint> indices( 2 );
  indices[ 0 ] = curv_idx;
  indices[ 1 ] = cabs_idx;
  readShapeGeometry( cin, stats, indices );
  uint nb = stats.samples( 0 );
  vector<double> x( nb + 1 );
  vector<double> y( nb + 1 );
  vector<double> theta( nb + 1 );
  x[ 0 ] = args.getOption( "-p0" )->getDoubleValue( 0 );
  y[ 0 ] = args.getOption( "-p0" )->getDoubleValue( 1 );
  theta[ 0 ] = args.getOption( "-t0" )->getDoubleValue( 0 );
  integrateShapeFromCurvatureAndCAbs( x, y, theta,
				      stats, 0, 1 );
  for ( uint i = 0; i <= nb; ++i )
    cout << x[ i ] << " " << y[ i ] << " " << theta[ i ] << endl;

  return 0;
}
