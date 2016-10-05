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
#include "ImaGene/mathutils/Mathutils.h"
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
		    // cerr << "Adding " << val << " to " << j << endl;
		  }
	      ++idx;
	    }
	}
      getline( in, str );
    }
}

struct DifferenceFunction 
{
  virtual double operator()( double x, double y ) const
  {
    return x - y;
  }
};

struct AngleDifferenceFunction : public DifferenceFunction 
{
  virtual double operator()( double x, double y ) const
  {
    return cos( x ) * sin( y ) - sin( x ) * cos( y );
  }
};

/**
 * @param values stores the three samples (expected values, estimated
 * values, measure of each value )
 */
void computeErrors( ostream & out, Statistics & values, uint ch, 
		    bool is_angle = false )
{
  Statistics stats( 8 );
  uint nb = values.samples( 0 );
  DifferenceFunction f1;
  AngleDifferenceFunction f2;
  const DifferenceFunction & diff = is_angle ? f2 : f1;

  for ( uint i = 0; i < nb; ++i )
    {
//       cerr << values.value( 0, i ) << " "
// 	   << values.value( 1, i ) << " "
// 	   << values.value( 2, i ) << endl;
      // Y-X
      stats.addValue( 0, diff( values.value( 1, i ), 
			       values.value( 0, i ) ) );

      // |Y-X|
      stats.addValue( 1, fabs( diff( values.value( 1, i ),
				     values.value( 0, i ) ) ) );
      
      // mu(Y-X)^2
      stats.addValue( 2,
		      Mathutils::sqr( diff( values.value( 1, i ),
					    values.value( 0, i ) ) )
		      * values.value( 2, i ) );

      // (Y-X)^2
      stats.addValue( 3,
		      Mathutils::sqr( diff( values.value( 1, i ),
					    values.value( 0, i ) ) )
		      );

      // (Y-X)/X
      stats.addValue( 4, ( diff( values.value( 1, i ),
				 values.value( 0, i ) ) )
		      / values.value( 0, i ) );

      // |Y-X|/X
      stats.addValue( 5, fabs( ( diff( values.value( 1, i ),
				       values.value( 0, i ) ) )
			       / values.value( 0, i ) ) );
      // X^2
      stats.addValue( 6, Mathutils::sqr( values.value( 0, i ) ) );
      // |X|
      stats.addValue( 7, fabs( values.value( 0, i ) ) );
    }
  stats.terminate();

  out << "# X=expected Y=estimated M=measure" << endl;
  if ( ( ch == 1 ) || ( ch == 0 ) )
    out << "# V=X-Y mean, std deviation, variance, unbiased variance, max, min" << endl
	<< stats.mean( 0 ) << " "
	<< sqrt( stats.variance( 0 ) ) << " "
	<< stats.variance( 0 ) << " "
	<< stats.unbiasedVariance( 0 ) << " "
	<< stats.max( 0 ) << " "
	<< stats.min( 0 ) << " "
	<< endl;
  if ( ( ch == 2 ) || ( ch == 0 ) )
  out << "# V=|X-Y| mean, std deviation, variance, unbiased variance, max, min" << endl
      << stats.mean( 1 ) << " "
      << sqrt( stats.variance( 1 ) ) << " "
      << stats.variance( 1 ) << " "
      << stats.unbiasedVariance( 1 ) << " "
      << stats.max( 1 ) << " "
      << stats.min( 1 ) << " "
      << endl;
  if ( ( ch == 3 ) || ( ch == 0 ) )
    out << "# V=(X-Y)/Y mean, std deviation, variance, unbiased variance, max, min" << endl
	<< stats.mean( 4 ) << " "
	<< sqrt( stats.variance( 4 ) ) << " "
	<< stats.variance( 4 ) << " "
	<< stats.unbiasedVariance( 4 ) << " "
	<< stats.max( 4 ) << " "
	<< stats.min( 4 ) << " "
	<< endl;
  if ( ( ch == 4 ) || ( ch == 0 ) )
    out << "# V=|(X-Y)/Y| mean, std deviation, variance, unbiased variance, max, min" << endl
	<< stats.mean( 5 ) << " "
	<< sqrt( stats.variance( 5 ) ) << " "
	<< stats.variance( 5 ) << " "
	<< stats.unbiasedVariance( 5 ) << " "
	<< stats.max( 5 ) << " "
	<< stats.min( 5 ) << " "
	<< endl;
  if ( ( ch == 5 ) || ( ch == 0 ) )
    out << "# V=M*(X-Y)^2 integral (L2-norm of difference)" << endl
	<< stats.mean( 2 ) * nb << endl;
  if ( ( ch == 6 ) || ( ch == 0 ) )
    out << "# V=(X-Y)^2 integral (sample squared error)" << endl
	<< stats.mean( 3 ) * nb << endl;
  if ( ( ch == 7 ) || ( ch == 0 ) )
    out << "# V=(X-Y)^2 / X^2 integral (relative sample squared error)" << endl
	<< ( stats.mean( 3 ) * nb ) / ( stats.mean( 6 ) * nb ) << endl;
  if ( ( ch == 8 ) || ( ch == 0 ) )
    out << "# V=|X-Y|_oo / |X|_oo integral (relative sample L-oo error)" << endl
	<< ( stats.max( 1 ) ) / ( stats.max( 7 ) ) << endl;
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
  args.addOption( "-expected_idx", "-expected_idx: specifies in which column is the expected value.", "1" );
  args.addOption( "-estimated_idx", "-estimated_idx: specifies in which column is the estimated value.", "2" );
  args.addOption( "-measure_idx", "-measure_idx: specifies in which column is the measure used for computing the error in L2-sense. 0 means that all values have measure 1.", "0" );
  args.addOption( "-choice", "-choice <n>: specifies which measure to perform, assuming X=expected Y=estimated M=measure: 0 all, 1 X-Y, 2: |X-Y|, 3: (X-Y)/Y, 4 |(X-Y)/Y|, 5: M*(X-Y)^2 integral (L2-norm of difference), 6 (X-Y)^2 integral (sample squared error), 7 (X-Y)^2 / X^2 integral (relative sample squared error), 8 |X-Y|_oo / |X|_oo integral (relative sample L-oo error).", "0" );
  args.addBooleanOption( "-angle_values", "-angle_values: specifies that the values are angles, errors are thus computed by transforming angles into unit vectors, and computing the determinant of the vectors." );
  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "error_analysis", 
			  "Computes several statistical errors given a set of samples wrt expected value.",
			  "" )
	   << endl;
      return 1;
    }
  
  uint exp_idx = args.getOption( "-expected_idx" )->getIntValue( 0 );
  uint est_idx = args.getOption( "-estimated_idx" )->getIntValue( 0 );
  uint m_idx = args.getOption( "-measure_idx" )->getIntValue( 0 );
  bool is_angle = args.check( "-angle_values" );
  // Force storing of samples.
  Statistics stats( 3, true );
  vector<uint> indices( 2 );
  indices[ 0 ] = exp_idx;
  indices[ 1 ] = est_idx;
  if ( m_idx != 0 )
    indices.push_back( m_idx );
  readShapeGeometry( cin, stats, indices );
  uint nb = stats.samples( 0 );
  cout << "# Read " << nb << " samples." << endl;
  if ( m_idx == 0 )
    for ( uint i = 0; i < nb; ++i )
      stats.addValue( 2, 1.0 );

  stats.terminate();

  computeErrors( cout, stats, args.getOption( "-choice" )->getIntValue( 0 ),
		 is_angle );

  return 0;
}
