///////////////////////////////////////////////////////////////////////////////
// Generates classical 2D contours from standard arguments.
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/base/Vector.h"
#include "ImaGene/base/Shapes.h"
#include "ImaGene/digitalnD/GridEmbedder.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/C4CIteratorOnSurface.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/ObjectBoundary.h"
#include "ImaGene/digitalnD/Frame2D.h"
#include "ImaGene/helper/ContourHelper.h"
#include "ImaGene/helper/ShapeHelper.h"

using namespace std;
using namespace ImaGene;


static Arguments args;

StarShaped* 
shapeFromArgs( KnSpace & ks, const Embedder & embedder,
	       KnCharSet & voxset,
	       Kn_sid & bel, uint & nb_bels )
{
  StarShaped* shape = ShapeHelper::makeStarShapedFromArgs( args );
  Vector2D dir( 1.0, 0.0 );
  cerr << "# --- digitizing shape " << endl;
  voxset = ShapeHelper::makeShape( ks, embedder, *shape );
  cerr << "# --- find start bel " << endl;
  bel = ShapeHelper::getBelOnStarShapedBoundary( ks, embedder, *shape, dir );
  cerr << "# --- compute contour size " << endl;
  nb_bels = KnShapes::sgetContourSize( ks, voxset, bel,
				       *( ks.sbegin_dirs( bel ) ) );
  
  return shape;
}

void
displayContour( KnSpace* ks, C4CIteratorOnSurface* cp, bool loop )
{
  cout << "#######################################################" << endl;
  cout << "# displayContour" << endl
       << "# displays the points on the boundary as couples." << endl
       << "# x y" << endl;
  ContourHelper::displayContour( cout, ks, cp, loop, 0, 1 );
}

void
displayContourWord( KnSpace* ks, C4CIteratorOnSurface* cp )
{
  cout << "#######################################################" << endl;
  cout << "# displayContourWord" << endl
       << "# displays the contour as a word with freeman codes." << endl
       << "# x y" << endl;
  ContourHelper::displayContourWord( cout, ks, cp, 0, 1 );
}

void
displayFreemanChain( KnSpace* ks, C4CIteratorOnSurface* cp )
{
  ContourHelper::displayFreemanChain( cout, ks, cp, 0, 1 );
}


void displayShapeGeometry( ostream & out, 
			   KnSpace* ks, 
			   C4CIteratorOnSurface* cp, 
			   StarShaped* shape,
			   GridEmbedder & embedder, 
			   double dh )
{
  Kn_sid sbel = cp->current();
  Frame2D frame;
  frame.init( ks, 0, 1 );
  frame.setSurfelFrame( sbel, cp->trackDir() );
  bool is_open;
  uint nb_surfels = C4CIterator::size( *cp, is_open );

  out << setprecision(15);
  out << "#######################################################" << endl;
  out << "# displayShapeGeometry" << endl
      << "# displays the geometry of the real shape." << endl 
      << "# Nb_surfels=" << nb_surfels << endl
      << "# dh=" << dh << endl
      << "# num xs ys curv tgt_angle t dabs cprev cabs surfel_el x(t) y(t)" << endl;
  
  uint i = 0;
  Kn_sid bel = sbel;
  double dabs = 0.5;
  float oldtheta = -1.0f;
  float len = 0.0f;
  do
    {
      Vector xbel( 2 );
      embedder.sembed( bel, xbel );
      // surfel's centroid angle to origin
      float theta = shape->parameter( xbel );  
      // theoretical tangent computed at surfel's centroid angle
      Vector2D t = shape->tangent( theta );     
      Vector2D coords = shape->x( theta );     
      // curvature
      double curv = (double) shape->curvature( theta );

      out << i 
	  << " " << xbel.ro( 0 )
	  << " " << xbel.ro( 1 )
	  << " " << curv 
	  << " " << atan2( t.y(), t.x() )
	  << " " << theta
	  << " " << dabs;
      if ( oldtheta < 0.0f )
	out << " 0.0";
      else
	{
	  float elen;
	  if ( theta >= oldtheta ) 
	    elen = shape->arclength( oldtheta, theta, 50 );
	  else
	    {
	      elen = 0.0f; //shape->arclength( oldtheta, theta, 40 );
	      theta = oldtheta;
	    }
	  out << " " << elen;
	  len += elen;
	}
      out << " " << len;
      dabs += 1.0;
      if ( dabs >= (double) nb_surfels )
	dabs -= (double) nb_surfels;

      // Compute surfel elementary length
      uint j = *( ks->sbegin_dirs( bel ) );
      bool jdirect = ks->sdirect( bel, j ); 
      Kn_sid base = ks->sincident( bel, j, ! jdirect );
      Kn_sid tip = ks->sincident( bel, j, jdirect );
      embedder.sembed( base, xbel );
      float theta1 = shape->parameter( xbel );  
      embedder.sembed( tip, xbel );
      float theta2 = shape->parameter( xbel );  
      out << " " << shape->arclength( theta1, theta2, 50 );
      out << " " << coords.x();
      out << " " << coords.y();
      out << endl;
      oldtheta = theta;
      if ( cp->next() == 0 ) break;
      
      ++i;
      bel = cp->current();
    }
  while ( bel != sbel );

  out << "# last_arclength total_length" << endl;
  Vector xbel( 2 );
  embedder.sembed( bel, xbel );
  // surfel's centroid angle to origin
  float theta = shape->parameter( xbel );  
  float elen = shape->arclength( oldtheta, theta, 50 );
  out << "# " << elen << " " << ( len + elen ) << endl;
  // TOFINISH
}


double perimeter( StarShaped* shape,
		uint nb_points )
{
  return (double)
    ( shape->arclength( 0, M_PI, nb_points / 2 )
      + shape->arclength( M_PI, M_PI+M_PI, nb_points / 2 ) );
}

void displayLength( ostream & out, 
		    KnSpace* ks, 
		    C4CIteratorOnSurface* cp, 
		    StarShaped* shape,
		    GridEmbedder & embedder, 
		    double dh,
		    uint nb_points )
{
  out << setprecision(15);
  Kn_sid sbel = cp->current();
  Frame2D frame;
  frame.init( ks, 0, 1 );
  frame.setSurfelFrame( sbel, cp->trackDir() );
  bool is_open;
  uint nb_surfels = C4CIterator::size( *cp, is_open );
  Kn_sid bel = sbel;
  double oldtheta = -1.0f;
  double len = 0.0f;
  do
    {
      Vector xbel( 2 );
      embedder.sembed( bel, xbel );
      // surfel's centroid angle to origin
      double theta = shape->parameter( xbel );  
      // theoretical tangent computed at surfel's centroid angle
      Vector2D t = shape->tangent( theta );     
      if ( oldtheta >= 0.0f )
	{
	  double elen;
	  if ( theta >= oldtheta ) 
	    elen = shape->arclength( oldtheta, theta, nb_points );
	  else
	    {
	      elen = 0.0f; //shape->arclength( oldtheta, theta, 40 );
	      theta = oldtheta;
	    }
	  len += elen;
	}
      oldtheta = theta;
      if ( cp->next() == 0 ) break;
      bel = cp->current();
    }
  while ( bel != sbel );

  Vector xbel( 2 );
  embedder.sembed( bel, xbel );
  // surfel's centroid angle to origin
  double theta = shape->parameter( xbel );  
  double elen = shape->arclength( oldtheta, theta, nb_points );
  out << ( len + elen ) << endl;
}


void drawShape( ostream & out,
		StarShaped* shape, uint nb )
{
  out << "#######################################################" << endl;
  out << "# drawShape" << endl
      << "# nbpoints=" << nb << endl;
  out << setprecision( 15 );
  for ( uint i = 0; i <= nb; ++i )
    {
      float t = i*2.0*M_PI / (float) nb;
      Vector2D pos( shape->x( t ) );
      out << pos.ro( 0 ) << " " << pos.ro( 1 ) << endl;
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
  StandardArguments::addDigitalArgs( args, 2, false, false );
  ShapeHelper::addStarShapedArgs( args );
  ShapeHelper::addNoiseArgs( args );
  args.addBooleanOption( "-v", "-v: verbose mode." );
  args.addOption( "-dB", "-dB <loop>: display boundary. loop=1 : rewrite the first point.", "1" );
  args.addBooleanOption( "-dBW", "-dBW: display boundary word with freeman codes." );
  args.addBooleanOption( "-dFC", "-dFC: display contour as a Freeman chain." );
  args.addBooleanOption( "-dPGM", "-dPGM: output shape in pgm format on the standard output stream." );
  args.addBooleanOption( "-dSG", "-dSG: output shape geometry (curvature, tangent, length, etc)." );
  args.addOption( "-dS", "-dS <nb>: output shape as a list of <nb> points.", "100" );
  args.addOption( "-seed", "-seed <nb>: change the seed of the random number generator to <nb>.", "0" );
  args.addOption( "-dL", "-dL <N>: output shape length estimated with a sampling of <N> per pixel.", "50" );

  if ( ( argc <= 1 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "gencontour", 
			  "Generates a 2D contour. The contour is the digitization of some classical shape as specified by the user. The contour may then be exported in different formats (list of points, contour word, Freeman chain, PGM image)",
			  "" )
	   << endl;
      return 1;
    }

  // -------------------------------------------------------------------------
  bool verbose = args.check( "-v" );
  if ( args.check( "-seed" ) )
    {
      srand( args.getOption( "-seed" )->getIntValue( 0 ) );
    }

  // -------------------------------------------------------------------------
  // Build space.
  uint d = StandardArguments::dim( args );
  if ( d != 2 )
    {
      cerr << "Dimension is 2." << endl;
      return 2;
    }
  Kn_size sizes[ d ];
  StandardArguments::fillSizes( args, sizes );
  KnSpace ks( 2, sizes );
  Kn_uid vcenter = ks.uspel( ks.ukcode( sizes ) ); // center
  Vector xcenter = ks.ucentroid( vcenter );
  
  // -------------------------------------------------------------------------
  // Get some parameters.
  double dh = (double) args.getOption( "-step" )->getDoubleValue( 0 );
  GridEmbedder embedder;
  embedder.init( &ks );
  embedder.setCenter( xcenter );
  embedder.setScale( dh );

  // -------------------------------------------------------------------------
  // Make shape.
  uint nb_bels = 0;
  Kn_sid bel = 0;
  KnCharSet voxset = KnCharSet::create( ks, ks.dim(), false, 0 );
  KnCharSet aux_voxset( voxset );
  
  StarShaped* shape = shapeFromArgs( ks, embedder, voxset, bel, nb_bels );

  // -------------------------------------------------------------------------
  // Change shape.
  Kn_uid in, out;
  bool interior = true;
  if ( args.check( "-noisify" ) )
    {
      ShapeHelper::addNoiseFromArgs( args, &ks, voxset, bel, in, out );
      bel = ShapeHelper::findInnerObject( &ks, voxset, in, aux_voxset );

      if ( args.check( "-inner_object" ) )
	{
	  // NB: change contour and object.
	  voxset = aux_voxset;
	  interior = true;
	}
      if ( args.check( "-outer_object" ) )
	{
	  // NB: change contour and object.
	  bel = ShapeHelper::findOuterObject( &ks, voxset, out, 
					      aux_voxset, in );
	  voxset = aux_voxset;
	  interior = false;
	}
    }
  else if ( args.check( "-outer_object" ) )
    interior = false;

  if ( verbose )
    {
      cerr << "# Bel:";
      ks.displayKn_sid( bel, cerr );
      cerr << endl;
    }

  // -------------------------------------------------------------------------
  // Create boundary.
  BelAdjacency badj( ks, interior );
  ObjectBoundary bdry( badj, voxset );
  uint dir = *( ks.sbegin_dirs( bel ) );
  C4CIteratorOnSurface* cp = bdry.newC4CIterator( bel, dir );

  bool open;
  uint nb = C4CIterator::size( *cp, open );
  if ( verbose )
    {
      cerr << "# Size boundary Nb=" << nb << " open=" << open << endl;
    }
  if ( args.check( "-dB" ) ) 
    displayContour( &ks, cp, 
		     (bool) args.getOption( "-dB" )->getIntValue( 0 ) );

  if ( args.check( "-dBW" ) ) 
    displayContourWord( &ks, cp );
  if ( args.check( "-dFC" ) ) 
    displayFreemanChain( &ks, cp );

  if ( args.check( "-dPGM" ) )
    ShapeHelper::exportToPGM( cout, &ks, voxset );

  if ( args.check( "-dSG" ) )
    displayShapeGeometry( cout, &ks, cp, shape, embedder, dh );

  if ( args.check( "-dS" ) )
    drawShape( cout, shape, args.getOption( "-dS" )->getIntValue( 0 ) );

  if ( args.check( "-dL" ) )
    {
      cout << setprecision( 15 ) 
	   << perimeter( shape, 
			 args.getOption( "-dL" )->getIntValue( 0 ) ) << endl;
    }
    // displayLength( cout, &ks, cp, shape, embedder, dh,
    // 		   args.getOption( "-dL" )->getIntValue( 0 ) );

  // Free some stuff.
  delete cp;
  delete shape;
   
  return 0;
}
