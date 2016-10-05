#include <iostream>
#include <cmath>
#include <cstdlib>

//#include "LinAlg/LinAlg2D/Vector2D.hpp"
//#include "ImageLib/Gauss/G.hpp"
#include "ImaGene/base/Vector.h"
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/base/Shapes.h"
#include "ImaGene/mathutils/G.h"
#include "ImaGene/mathutils/Functions.h"
#include "ImaGene/mathutils/Polynomial.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/dgeometry2d/C4CSegment.h"
#include "ImaGene/dgeometry2d/C4CIterator.h"
#include "ImaGene/dgeometry2d/C4CGeometry.h"
#include "ImaGene/dgeometry2d/C4CSegmentPencil.h"
#include "ImaGene/dgeometry2d/EuclideanGeometry.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/KnRCellSet.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/ObjectBoundary.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/Embedder.h"
#include "ImaGene/digitalnD/GridEmbedder.h"
#include "ImaGene/digitalnD/KnSpaceScanner.h"

using namespace std;
using namespace ImaGene;

static Arguments args;


class Statistics
{
  uint m_nb;
  uint* m_samples;
  float* m_exp;
  float* m_exp2;
  float* m_var;
  float* m_unbiased_var;

public:
  Statistics( uint size )
  {
    m_samples = 0;
    m_exp = 0;
    m_exp2 = 0;
    m_var = 0;
    m_unbiased_var = 0;
    init( size );
  }

  ~Statistics()
  {
    erase();
  }

  uint nb() const
  {
    return m_nb;
  }

  uint samples( uint k ) const
  {
    return m_samples[ k ];
  }

  float mean( uint k ) const
  {
    return m_exp[ k ];
  }

  float variance( uint k ) const
  {
    return m_var[ k ];
  }

  float unbiasedVariance( uint k ) const
  {
    return m_unbiased_var[ k ];
  }
  

  void addValue( uint k, float v )
  {
    m_samples[ k ] += 1;
    m_exp[ k ] += v;
    m_exp2[ k ] += v*v;
  }
  
  void terminate()
  {
    for ( uint k = 0; k < m_nb; ++k )
      {
	m_exp[ k ] /= m_samples[ k ];
	m_exp2[ k ] /= m_samples[ k ];
	m_var[ k ] = m_exp2[ k ] - m_exp[ k ] * m_exp[ k ];
	m_unbiased_var[ k ] = m_samples[ k ] * m_var[ k ] 
	  / ( m_samples[ k ] - 1 );
      }
  }

  void init( uint size )
  {
    erase();
    m_nb = size;
    m_samples = new uint[ size ];
    m_exp = new float[ size ];
    m_exp2 = new float[ size ];
    m_var = new float[ size ];
    m_unbiased_var = new float[ size ];
    clear();
  }

  void clear()
  {
    if ( m_nb == 0 ) return;
    for ( uint i = 0; i < m_nb; ++ i )
      {
	m_samples[ i ] = 0;
	m_exp[ i ] = 0.0;
	m_exp2[ i ] = 0.0;
	m_var[ i ] = 0.0;
	m_unbiased_var[ i ] = 0.0;
      }
  }
  
  void erase() 
  {
    if ( m_samples != 0 ) delete[] m_samples;
    if ( m_exp != 0 ) delete[] m_exp;
    if ( m_exp2 != 0 ) delete[] m_exp2;
    if ( m_var != 0 ) delete[] m_var;
    if ( m_unbiased_var != 0 ) delete[] m_unbiased_var;
    m_samples = 0;
    m_exp = 0;
    m_exp2 = 0;
    m_var = 0;
    m_unbiased_var = 0;
    m_nb = 0;
  }
  
};


/**
 * @return a uniform random value between 0.0 (included) and (1.0) excluded.
 */
double random1()
{
  return ((double) rand())/(RAND_MAX+1.0);
}

float 
norm( const Vector2D & v )
{
  return sqrt( v.ro( 0 ) * v.ro( 0 ) + v.ro( 1 ) * v.ro( 1 ) );
}

void
normalize( Vector2D & v )
{
  float n = norm( v );
  v /= n;
}

float
distance( const Vector2D & v1, const Vector2D & v2 )
{ 
  Vector2D v3( v1 );
  v3 -= v2;
  return norm( v3 );
}

float
dotProduct( const Vector2D & v1, const Vector2D & v2 )
{ 
  return v1.ro( 0 ) * v2.ro( 0 ) + v1.ro( 1 ) * v2.ro( 1 );
}

float
det( const Vector2D & v1, const Vector2D & v2 )
{ 
  return v1.ro( 0 ) * v2.ro( 1 ) - v1.ro( 1 ) * v2.ro( 0 );
}



uint 
getNbSurfels( KnSpace & ks, KnCharSet voxset, Kn_sid starting_bel )
{
  // Get tangent plane.
  BelAdjacency badj( ks, true );
  ObjectBoundary bdry( badj, voxset );
  uint track_dir = *( ks.sbegin_dirs( starting_bel ) );
  C4CIteratorOnSurface* cp = bdry.newC4CIterator( starting_bel, track_dir );
  Kn_sid bel = starting_bel;
  uint nb = 0;
  do 
    {
      // Go to next one.
      cp->next();
      bel = cp->current();
      ++nb;
    }
  while ( bel != starting_bel );
  return nb;
}



struct ShapeParameters 
{
  KnSpace* p_ks;
  const Embedder* p_embedder;
  KnCharSet voxset;
  Kn_sid bel;
  uint nb_bels;
  StarShaped* sshape;
  ShapeParameters( KnSpace & ks )
    : p_ks ( &ks ), p_embedder( 0 ), voxset( KnCharSet::ucreate( ks, 0, 0 ) )
  {}
  
    
};


Kn_sid
surfelOnStarShapeBoundary( KnSpace & ks,
			   const Embedder & embedder,
			   const StarShaped & starshape,
			   const Vector2D & dir )
{
  Vector2D ll; 
  Vector2D ur;
  embedder.uembed( ks.uspel( ks.ufirst() ), ll );
  embedder.uembed( ks.uspel( ks.ulast() ), ur );
  Vector2D lr( ur.ro( 0 ), ll.ro( 1 ) );
  Vector2D ul( ll.ro( 0 ), ur.ro( 1 ) );
  
  Vector2D shapec = starshape.center();
  
  Vector2D rel_ur( ur ); rel_ur -= shapec;
  Vector2D rel_ul( ul ); rel_ul -= shapec;
  Vector2D rel_ll( ll ); rel_ll -= shapec;
  Vector2D rel_lr( lr ); rel_lr -= shapec;
  
  float ur_ccw_to_dir = det( dir, rel_ur );
  float ul_ccw_to_dir = det( dir, rel_ul );
  float ll_ccw_to_dir = det( dir, rel_ll );
  float lr_ccw_to_dir = det( dir, rel_lr );

  float x, y;
  int orientation = -1;
  if ( ( ur_ccw_to_dir >= 0.0 ) && ( lr_ccw_to_dir <= 0.0 ) )
    // Test with Right.
    {
      x = ur.ro( 0 );
      y = shapec.ro( 1 ) + ( x - shapec.ro( 0 ) ) * dir.ro( 1 ) / dir.ro( 0 ) ;
      orientation = 0;
    }
  else if ( ( ul_ccw_to_dir >= 0.0 ) && ( ur_ccw_to_dir <= 0.0 ) )
    // Test with Up.
    {
      y = ur.ro( 1 );
      x = shapec.ro( 0 ) + ( y - shapec.ro( 1 ) ) * dir.ro( 0 ) / dir.ro( 1 ) ;
      orientation = 1;
    }
  else if ( ( ll_ccw_to_dir >= 0.0 ) && ( ul_ccw_to_dir <= 0.0 ) )
    // Test with Left.
    {
      x = ll.ro( 0 );
      y = shapec.ro( 1 ) + ( x - shapec.ro( 0 ) ) * dir.ro( 1 ) / dir.ro( 0 ) ;
      orientation = 2;
    }
  else if ( ( lr_ccw_to_dir >= 0.0 ) && ( ll_ccw_to_dir <= 0.0 ) )
    // Test with Low.
    {
      y = ll.ro( 1 );
      x = shapec.ro( 0 ) + ( y - shapec.ro( 1 ) ) * dir.ro( 0 ) / dir.ro( 1 ) ;
      orientation = 3;
    }
  else
    // ??
    {
      cerr << "Geometric ERROR." << endl;
      x = y = 0.0;
    }

  Vector2D limit( x, y );
  // Starting dichotomy.
  Kn_uid lintv = embedder.uinverseSpel( shapec );
  Kn_uid lextv = embedder.uinverseSpel( limit );
  embedder.uembed( lintv, shapec );
  embedder.uembed( lextv, limit );
  if ( ( ! starshape.isInside( shapec ) ) 
       || ( starshape.isInside( limit ) ) )
      cerr << "Inside/Outside ERROR." << endl;
  while ( true )
    {
      Vector2D middle( shapec );
      middle += limit;
      middle *= 0.5;
      Kn_uid middlev = embedder.uinverseSpel( middle );
      embedder.uembed( middlev, middle );
      if ( ( middlev == lintv ) || ( middlev == lextv ) )
	break;
      
//       cerr << "M( " << middle.ro( 0 ) << " , " << middle.ro( 1 ) << " )" 
// 	   << middle 
//  	   << endl;
      if ( starshape.isInside( middle ) )
	{
	  shapec = middle;
	  lintv = middlev;
	}
      else
	{
	  limit = middle;
	  lextv = middlev;
	}
    }

   cerr << "VoxC=";
   ks.displayKn_uid( lintv, cerr );
   cerr << endl << "VoxL=";
   ks.displayKn_uid( lextv, cerr );
   cerr << endl;

  // 2nd loop.
  Kn_sid bsurfel = 0;
  uint lvl;
  if ( ! ks.uareAdjacent( lintv, lextv, lvl ) )
    cerr << "Adjacency error" << endl;
  Kn_uid interm = lintv;
  Vector xinterm( 2 );
  uint i = 0;
  int dx = 0;
  for ( ; i < 2; ++i )
    {
      dx = ( (int) ks.udecodeCoord( lextv, i ) )
	- ( (int) ks.udecodeCoord( lintv, i ) );
      cerr << "dx[" << i << "]=" << dx << endl;
      if ( dx == 0 ) continue;
      if ( dx > 0 ) interm = ks.ugetIncr( interm, i );
      else          interm = ks.ugetDecr( interm, i );
      embedder.uembed( interm, xinterm );
      cerr << "I( " << xinterm.ro( 0 ) << " , " << xinterm.ro( 1 ) << "|"
	   << starshape.isInside( xinterm ) << ")" << endl;
      if ( ! starshape.isInside( xinterm ) )
	  break;
    }
  bsurfel = ks.sincident( ks.signsNeg( interm ), i, dx < 0 );
  
//   Kn_sid bdry_v = ks.signsPos( lintv );
//   Kn_sid bsurfel;
//   switch ( orientation )
//     {
//     case 0: bsurfel = ks.sincident( bdry_v, 0, true ); break;
//     case 1: bsurfel = ks.sincident( bdry_v, 1, true ); break;
//     case 2: bsurfel = ks.sincident( bdry_v, 0, false ); break;
//     case 3: bsurfel = ks.sincident( bdry_v, 1, false ); break;
//     default: bsurfel = 0;
//     }
  uint j = ks.sorthDir( bsurfel );
  Kn_uid check_int = ks.unsigns( ks.sincident( bsurfel, j,
					       ks.sdirect( bsurfel, j ) ) );
  Kn_uid check_ext = ks.unsigns( ks.sincident( bsurfel, j,
					       ! ks.sdirect( bsurfel, j ) ) );
  Vector2D check_x;
  embedder.uembed( check_int, check_x );
  if ( ! starshape.isInside( check_x ) )
    cerr << "INSIDE ERROR." << endl;
  embedder.uembed( check_ext, check_x );
  if ( starshape.isInside( check_x ) )
    cerr << "OUTSIDE ERROR." << endl;

  return bsurfel;
}



// StarShaped* shapeFromArgs( KnSpace & ks, const Embedder & embedder,
// 			   KnCharSet & voxset,
// 			   Kn_sid & bel, uint & nb_bels )
void shapeFromArgs( KnSpace & ks, 
		    const Embedder & embedder,
		    ShapeParameters & shape_params )
{
  bool circle = args.check( "-circle" );
  bool flower = args.check( "-flower" );
  bool accflower = args.check( "-accflower" );
  bool rsquare = args.check( "-rsquare" );
  bool fuzzy_uniform = args.check( "-fuzzy_uniform" );
  float dh = args.getOption( "-step" )->getFloatValue( 0 );
  
  float xc = 0.0; 
  float yc = 0.0;
  if ( fuzzy_uniform )
    { 
      xc += ( random1() - 0.5 ) * dh;
      yc += ( random1() - 0.5 ) * dh;
    }
  Vector2D dir( 1.0, 0.0 );

  cerr << "Shapes " << xc << " " << yc << endl;
  StarShaped* shape = 0;
  if ( circle )
    {
      shape =
	new Circle( xc, yc, args.getOption( "-circle" )->getFloatValue( 0 ) );
    }
  else if ( flower )
    {
      shape =
	new Flower( xc, yc, 
		    args.getOption( "-flower" )->getFloatValue( 0 ),
		    args.getOption( "-flower" )->getFloatValue( 1 ),
		    args.getOption( "-flower" )->getIntValue( 2 ),
		    args.getOption( "-flower" )->getFloatValue( 3 ) );
    }
  else if ( accflower )
    {
      shape =
	new AccFlower( xc, yc, 
		       args.getOption( "-accflower" )->getFloatValue( 0 ),
		       args.getOption( "-accflower" )->getFloatValue( 1 ),
		       args.getOption( "-accflower" )->getIntValue( 2 ) );
    }
  else if ( rsquare )
    {
      shape =
	new RoundedSquare( xc, yc, 
			   args.getOption( "-rsquare" )->getFloatValue( 0 ),
			   args.getOption( "-rsquare" )->getFloatValue( 1 ) );
    }
  shape_params.p_embedder = &embedder;

  cerr << "makeShape" << endl;
  
  // Filling it by scanning the space.
  Kn_uid low = ks.uspel( ks.ufirst() );
  Kn_uid up = ks.uspel( ks.ulast() );
  KnCharSet s = KnCharSet::ucreate( ks, low, 0 );
  KnSpaceScanner scanner( ks, low, up );
  Kn_uid p = scanner.lower_bound;
  Vector vp( ks.dim() );
  do 
    {
      embedder.uembed( p, vp );
      if ( shape->isInside( vp ) )
	s += p;
    }
  while ( scanner.increment( p ) );

  shape_params.voxset = s;
  shape_params.bel = surfelOnStarShapeBoundary( ks, embedder, *shape, dir );
  shape_params.nb_bels = getNbSurfels( ks, s, shape_params.bel );



  //cerr << "Set=" << voxset << endl;
  //cerr << "Bel=";
  //ks.displayKn_sid( bel, cerr );
  //cerr << endl;

  shape_params.sshape = shape;
}





///////////////////////////////////////////////////////////////////////////////
// TANGENT COMPUTATION
///////////////////////////////////////////////////////////////////////////////

class AnalogComputer
{
public:
  AnalogComputer() {}
  virtual ~AnalogComputer() {}
  virtual Vector2D computeAnalog( ShapeParameters & shape,
				  C4CIteratorOnSurface & iter )
  {
    return shape.p_ks->scentroid( iter.current() );
  }
};


class AnalogComputerWithSegmentPencil : public AnalogComputer
{
  C4CSegment m_segments[ 100 ];
  R2RFunction* m_lambda;
  R2RFunction* m_lambda_p;

public:
  enum InterpolationFct {
    TRIANGLE, BELLSHAPE2 
  };

  AnalogComputerWithSegmentPencil( InterpolationFct type_fct )
  {
    if ( type_fct == TRIANGLE )
      {
	m_lambda = new TriangleFunction;
	m_lambda_p = new DTriangleFunction;
      }
    else
      {
	m_lambda = new Polynomial( Polynomial::bellShape2( 1.0f ) );
	m_lambda_p = new Polynomial( Polynomial::bellShape2( 1.0f )
				     .derivative() );
      }
  }
  
  ~AnalogComputerWithSegmentPencil()
  {
    if ( m_lambda != 0 ) delete m_lambda;
    if ( m_lambda_p != 0 ) delete m_lambda_p;
  }
  virtual Vector2D computeAnalog( ShapeParameters & shape,
				  C4CIteratorOnSurface & iter )
  {
    C4CIterator* it = iter.clone();
    uint j = 0;
    uint k;
    uint m = 100;
    if ( ! C4CGeometry::maximalSegments( *it, m_segments, j, k, m ) )
      cerr << "[AnalogComputerWithSegmentPencil::computeAnalog]"
	   << " Not enough segments." << endl;
    // Polynomial lambda = Polynomial::bellShape2( 1.0f );
    // Polynomial lambda_p = lambda.derivative();
    C4CSegmentPencil pencil( m_segments, j, k, m, *m_lambda, *m_lambda_p );
    Vector2D r = pencil.continuousAnalog( Vector2D( 0.5f, 0.0f ) );
    Kn_sid b = iter.current();
    KnSpace* space = shape.p_ks;
    uint jdir = space->sorthDir( b );
    bool orth_pos = space->sdirect( b, jdir );
    uint idir = *( space->sbegin_dirs( b ) );
    bool track_pos = space->sdirect( b, idir );
    // Shift of 0.5 to get surfel centroid at (0...0).
    r.rw( 0 ) = track_pos ? r.ro( 0 ) - 0.5f  : 0.5f - r.ro( 0 );
    r.rw( 1 ) = orth_pos ? - r.ro( 1 ) : r.ro( 1 );
    Vector2D xc;
    space->scentroid( b, xc );
    xc.rw( idir ) += r.ro( 0 );
    xc.rw( jdir ) += r.ro( 1 );
    shape.p_embedder->embedVector( xc, r );
    // cout << "(" << r.x() << " " << r.y() << ")";
    delete it;
    return r;
  }
  
  
};


void experiment( ShapeParameters & shapeinfo, AnalogComputer* computer )
{
  // Get tangent plane.
  KnSpace & ks = *shapeinfo.p_ks;
  KnCharSet & voxset = shapeinfo.voxset;
  Kn_sid starting_bel = shapeinfo.bel;
  
  BelAdjacency badj( ks, true );
  ObjectBoundary bdry( badj, voxset );
  uint track_dir = *( ks.sbegin_dirs( starting_bel ) );
  C4CIteratorOnSurface* cp = bdry.newC4CIterator( starting_bel, track_dir );
  Kn_sid bel = starting_bel;
  uint nb = 0;
  cout << endl;
  Vector2D xbel;
  Vector2D p;
  Vector2D m;
  
  do 
    {
      // Make experiment.
      shapeinfo.p_embedder->sembed( bel, xbel );
      p = computer->computeAnalog( shapeinfo, *cp );
      float t = shapeinfo.sshape->parameter( p );
      m = shapeinfo.sshape->x( t );
      
      cout << nb << " " << xbel.x() << " " << xbel.y()
	   << " " << p.x() << " " << p.y()
	   << " " << m.x() << " " << m.y() << endl;

      // Go to next one.
      cp->next();
      bel = cp->current();
      ++nb;
    }
  while ( bel != starting_bel );
}



int
main( int argc, char** argv ) 
{
  StandardArguments::addDigitalArgs( args, 2, false, false );

  args.addOption( "-step", "-step <h>: the discretization step or resolution, the closer to 0, the finer.", "1.0" );
  args.addOption( "-circle", "-circle <R>: the test shape is a circle of radius R", "10.0" );
  args.addOption( "-flower", "-flower <R> <r> <k> <phi>: the test shape is a flower with k extremeties with mean radius R and variability of radius r, phi is the phase of the flower", "3", "10.0", "5.0", "0.0" );
  args.addOption( "-accflower", "-flower <R> <r> <k>: the test shape is a phase accelerating flower with k extremeties with mean radius R and variability of radius r", "4.0", "2.0", "4" );
  args.addOption( "-rsquare", "-rsquare <R> <r>: the test shape is rounded square of big radius R and small corner radius r", "4.0", "1.0" );
  
  args.addBooleanOption( "-fuzzy_uniform", "-fuzzy_uniform: possible uniform random shift of the shape center." );

  args.addOption( "-nbtrials", "-nbtrials <N>: the number of experiments over which the result is averaged", "1" );
  
  args.addOption( "-interp", "-interp <B|T>: choose interpolation function between bell shape (B) and triangle shape (T).", "B" );
  

  if ( ( argc <= 1 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "test_Analog", 
			  "Tests several discrete tangent estimators and output some statistics.",
			  "-h -x -y -fuzzy_uniform -nbtrials -circle -flower -accflower -rsquare -step -thetas -v_angle_dev -v_window -v_angle_to_x -v_curvature -tgt_computer -blur -v_estcurv -v_r -curv_computer -integrate_curve -integrate_curve2 -integrate_truecurve2" ) 
	   << endl;
      return 1;
    }

  // Display command line.
  cout << "#";
  for ( int i = 0; i < argc; ++i )
    cout << " " << argv[ i ];
  cout << endl;

  // -------------------------------------------------------------------------
  // Build space.
  Kn_size sizes[ 2 ];
  StandardArguments::fillSizes( args, sizes );
  KnSpace ks( 2, sizes );
  Kn_uid vcenter = ks.uspel( ks.ukcode( sizes ) ); // center
  Vector xcenter = ks.ucentroid( vcenter );
  
  cerr << "--- Space: " << ks << endl;

  // -------------------------------------------------------------------------
  // Get some parameters.
  float dh = args.getOption( "-step" )->getFloatValue( 0 );
  uint nb_trials = args.getOption( "-nbtrials" )->getIntValue( 0 );

  cerr << "Embedder" << endl;
  GridEmbedder embedder;
  embedder.init( &ks );
  embedder.setCenter( xcenter );
  embedder.setScale( dh );

  // Analog computer
  AnalogComputer* analog_computer;
  if ( args.getOption( "-interp" )->getValue( 0 ) == "T" )
    analog_computer = new AnalogComputerWithSegmentPencil
      ( AnalogComputerWithSegmentPencil::BELLSHAPE2 );
  else // Default is "B"
    analog_computer = new AnalogComputerWithSegmentPencil
      ( AnalogComputerWithSegmentPencil::TRIANGLE );

  for ( uint n = 0; n < nb_trials; ++n )
    {
      cerr << " ---------------- Trial " << n << " ---------------" << endl;

      ShapeParameters shape_params( ks );
      
      shapeFromArgs( ks, embedder, shape_params );

      experiment( shape_params, analog_computer );
    }
}




  
			  
