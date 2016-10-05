///////////////////////////////////////////////////////////////////////////////
// Generates polygons from digital contours
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/Proxy.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/base/Vector.h"
#include "ImaGene/mathutils/Line2D.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/dgeometry2d/C4CSegment.h"
#include "ImaGene/dgeometry2d/C4CGeometry.h"
#include "ImaGene/dgeometry2d/C4CSegmentPencil.h"
#include "ImaGene/dgeometry2d/FreemanChain.h"
#include "ImaGene/dgeometry2d/C4CTangentialCover.h"
#include "ImaGene/dgeometry2d/C4CIteratorOnFreemanChain.h"
#include "ImaGene/digitalnD/C4CIteratorOnFreemanChainSurface.h"
#include "ImaGene/helper/C4CTangentialCoverGeometry.h"
#include "ImaGene/helper/ShapeHelper.h"


using namespace std;
using namespace ImaGene;

static Arguments args;

class LengthEstimator 
{
public:
  LengthEstimator() {}
  virtual void init( const C4CIterator & it,
		     const C4CIterator & it_end ) {}
  virtual ~LengthEstimator() {}
  virtual double length( double dh ) 
  {
    return 0.0;
  }
};

class LengthEstimatorBySymmetricTangent : public LengthEstimator
{
private:
  double m_length;
public:
  LengthEstimatorBySymmetricTangent() {}
  virtual ~LengthEstimatorBySymmetricTangent() {}

  virtual void init( const C4CIterator & it,
		     const C4CIterator & it_end )
  {
    m_length = 0.0;
    LengthEstimator::init( it, it_end );
    Proxy<C4CIterator> ptr_it_cur( it.clone() );
    do {
      Proxy<C4CIterator> ptr_it_front( ptr_it_cur->clone() );
      Proxy<C4CIterator> ptr_it_back( ptr_it_cur->clone() );
      C4CSegment segment = C4CGeometry::symmetricTangent( *ptr_it_front, 
							  *ptr_it_back, 
							  0 );
      double a = (double) segment.a();
      double b = (double) segment.b();
      double l = b / sqrt( a*a + b*b );
      //      cerr << "a=" << a << " b=" << b << " l=" << l << endl;
      m_length += l;
      // C4CGeometry::nextSymmetricSegment( segment, 
      // 					 *ptr_it_cur, 
      // 					 *ptr_it_front, 
      // 					 *ptr_it_back );
      ptr_it_cur->next();
    } while ( ! ptr_it_cur->equals( it_end ) );

  }
  virtual double length( double dh ) 
  {
    return dh * m_length;
  }

};

class LengthEstimatorByLambdaMSTangent : public LengthEstimator
{
private:
  double m_length;
  KnSpace & myKS;
public:
  LengthEstimatorByLambdaMSTangent( KnSpace & ks ) : myKS( ks ) {}
  virtual ~LengthEstimatorByLambdaMSTangent() {}

  virtual void init( const C4CIterator & it,
		     const C4CIterator & it_end )
  {
    m_length = 0.0;
    // const C4CIteratorOnSurface & input_it =
    //   dynamic_cast<const C4CIteratorOnSurface &>( it );
    const C4CIterator & input_it = it;

    // Functions used for combining all the angle of the maximal
    // segments inside the pencil.
    TriangleFunction l;
    DTriangleFunction lp;
    
    // Memorizes starting point.
    Proxy<C4CIterator> cur_it
      ( (C4CIterator*) input_it.clone() );
    // The following loop will display the estimated tangent direction
    // at each surfel. It needs a Frame to embed the digital contour.
    //  Frame2D frame;
    //frame.init( &myKS, 0, 1 );
    // A pencil has experimentally no more than 7 segments.
    const uint m = 20;
    C4CSegment segments[ m ];
    uint idx = 0;
    do
      {
	// Builds the pencil of maximal segments.
	uint j = 0;
	uint k;
	if ( C4CGeometry::maximalSegments( *cur_it, segments, j, k, m ) )
	  {
	    // All geometric computations are made in the local frame of
	    // the current boundary element.
	    C4CSegmentPencil pencil( segments, j, k, m, l, lp );
	    double theta = pencil.angleToX( Vector2D( 0.5, 0.0 ) );
	    m_length += cos( theta );
	    // // Cast angle in the global frame.
	    // Kn_sid surfel = cur_it->current();
	    // frame.setSurfelFrame( surfel, ks.stanDir( surfel ) );
	    // cout << idx << " " << frame.angleToX( theta ) << endl;
	  }
	// Go to next element.
	++idx;
	if ( cur_it->next() == 0) break;
      }
    while ( ! cur_it->equals( input_it ) );
  }

  virtual double length( double dh ) 
  {
    return dh * m_length;
  }
};


class LengthEstimatorByNbSteps : public LengthEstimator
{
private:
  double m_length;
public:
  LengthEstimatorByNbSteps() {}
  virtual ~LengthEstimatorByNbSteps() {}

  virtual void init( const C4CIterator & it,
		     const C4CIterator & it_end )
  {
    m_length = 0.0;
    LengthEstimator::init( it, it_end );
    Proxy<C4CIterator> ptr_it_cur( it.clone() );
    do {
      m_length += 1.0;
      ptr_it_cur->next();
    } while ( ! ptr_it_cur->equals( it_end ) );

  }
  virtual double length( double dh ) 
  {
    return dh * m_length;
  }
};

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
  StandardArguments::addDebugArgs( args, true, true, true );
  //  StandardArguments::addDigitalArgs( args, 2, false, false );
  StandardArguments::addIOArgs( args, true, false );
  args.addOption( "-step", "-step <h>: the discretization step or resolution, the closer to 0, the finer.", "1.0" );
  args.addBooleanOption( "-lenMLP", "-lenMLP: the length estimation is the length of the minimal length polygon of the digital contour." );
  args.addBooleanOption( "-lenTanSym", "-lenTanSym: the length estimation is the integration of the tangent orientations estimated by the symmetric tangent estimator." );
  args.addBooleanOption( "-lenLambdaMST", "-lenLambdaMST: the length estimation is the integration of the tangent orientations estimated by the lambda-Maximal Segment Tangent estimator." );
  args.addBooleanOption( "-lenNbSteps", "-lenNbSteps: the length estimation is the number of steps." );

  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "length_tan", 
			  "Estimates the length of a digital contour defined as a Freeman chain code. Several estimators are possible.",
			  "" )
	   << endl;
      return 1;
    }

  uint trace = args.getOption( "-trace" )->getIntValue( 0 );
  bool debug = args.check( "-debug" );
  uint timing = args.getOption( "-timing" )->getIntValue( 0 );
  double dh = (double) args.getOption( "-step" )->getDoubleValue( 0 );
  // -------------------------------------------------------------------------
  // Read Freeman chain and creates space/contour.

  FreemanChain c;
  istream & in_str = StandardArguments::openInput( args );
  FreemanChain::read( in_str, c );
  if ( ! in_str.good() )
    {
      if ( debug )
	cerr << "Error reading Freeman chain code." << endl;
      return 2;
    }

  KnSpace* ks = ShapeHelper::makeSpaceFromFreemanChain( c );

  C4CIteratorOnFreemanChain itfc;
  itfc.init( c.begin(), true );
  C4CIteratorOnFreemanChainSurface itfcs;
  itfcs.init( ks, itfc );

  cout << setprecision(15);

  if ( args.check( "-lenMLP" ) )
    {
      double length = dh * FreemanChain::lengthMLP( c );
      // cout << "# length by lenMLP" << endl;
      cout << length << " ";
    }
  if ( args.check( "-lenTanSym" ) )
    {
      LengthEstimatorBySymmetricTangent LE;
      LE.init( itfc, itfc );
      // cout << "# length by lenTanSym" << endl;
      cout << LE.length( dh ) << " ";
    }
  if ( args.check( "-lenLambdaMST" ) )
    {
      LengthEstimatorByLambdaMSTangent LE( *ks );
      LE.init( itfc, itfc );
      // cout << "# length by lenTanSym" << endl;
      cout << LE.length( dh ) << " ";
    }
  if ( args.check( "-lenNbSteps" ) )
    {
      LengthEstimatorByNbSteps LE;
      LE.init( itfc, itfc );
      // cout << "# length by lenTanSym" << endl;
      cout << LE.length( dh ) << " ";
    }
  cout << endl;
  return 0;
}
