///////////////////////////////////////////////////////////////////////////////
// Test module for digitalnD and Kn_uid, Kn_sid
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <cstring>

#include "ImaGene/base/RnMaps.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/KnSpaceScanner.h"
#include "ImaGene/digitalnD/KnSpaceCoordScanner.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/KnRCellSet.h"


using namespace ImaGene;
using namespace std;

/**
 * Creates a set of unsigned volumic cells whose elements represents a
 * sphere of center [center] and radius [inclusive_radius].
 * @param ks any space.
 * @param center a cell that is the center of the sphere.
 * @param inclusive_radius the radius of the sphere.
 * @return a set of spel-like cells.
 */
KnRCellSet umakeVolumicSphere2( const KnSpace & ks,
			     Kn_uid center,
			     float inclusive_radius )
{
  uint k;

  // Computing sphere bounding cube.
  inclusive_radius = inclusive_radius >= 0 ? inclusive_radius : 0.0f;
  uint r = (int) ceil( inclusive_radius );
  Kn_uid low_uid = ks.ufirst( center );
  Kn_uid up_uid = ks.ulast( center );
//   Kn_uid low_uid = KnShapes::ugetCubeLowerBound( ks, center, r );
//   Kn_uid up_uid = KnShapes::ugetCubeUpperBound( ks, center, r );

  // Instantiating set.
  KnRCellSet s = KnRCellSet::create( ks, ks.dim(), false, 0 );
  Vector vcenter = ks.ucentroid( center );

  RnMap f = RnMap::makeImplicitSphere( inclusive_radius );
  Vector scale( vcenter.size() );
  for ( k = 0; k < scale.size(); ++k )
    scale.rw( k ) = 1.0;
  scale.rw( 0 ) = 0.5;
  RnMap g = RnMap::makeScalingMap( scale );
  Vector neg_vcenter( vcenter );
  neg_vcenter *= -1.0;
  RnMap h = RnMap::makeTranslationMap( neg_vcenter );
  RnMap fgh = f( g( h ) );
  
  //float r2 = inclusive_radius * inclusive_radius;

  // Filling it by scanning the bounding cube.
  KnSpaceCoordScanner scanner( ks, low_uid, up_uid );
  Kn_uid p = scanner.lower_bound;
  Vector vp = ks.ucentroid( p );
  do 
    {
      //cout << ( f.eval( vp ) ).ro( 0 ) << " ";
      
      if ( ( fgh.eval( vp ) ).ro( 0 ) <= 0.0 )
	s += p;
    }
  while ( scanner.increment( p, vp ) );
  return s;
}


/**
 * Creates a set of unsigned volumic cells whose elements represents a
 * sphere of center [center] and radius [inclusive_radius].
 * @param ks any space.
 * @param center a cell that is the center of the sphere.
 * @param inclusive_radius the radius of the sphere.
 * @return a set of spel-like cells.
 */
KnRCellSet umakeVolumicSphere( const KnSpace & ks,
			    Kn_uid center,
			    float inclusive_radius )
{
  // Computing sphere bounding cube.
  inclusive_radius = inclusive_radius >= 0 ? inclusive_radius : 0.0f;
  uint r = (int) ceil( inclusive_radius );
  Kn_uid low_uid = KnShapes::ugetCubeLowerBound( ks, center, r );
  Kn_uid up_uid = KnShapes::ugetCubeUpperBound( ks, center, r );

  // Instantiating set.
  KnRCellSet s = KnRCellSet::create( ks, ks.dim(), false, 0 );
  Vector vcenter = ks.ucentroid( center );
  float r2 = inclusive_radius * inclusive_radius;

  // Filling it by scanning the bounding cube.
  KnSpaceCoordScanner scanner( ks, low_uid, up_uid );
  uint k;
  Kn_uid p = scanner.lower_bound;
  Vector vp = ks.ucentroid( p );
  do 
    {
      float dist2 = 0.0f;
      cout << "(" << vp.ro( 0 ) << " " << vp.ro( 1 ) << ") ";
      for ( k = 0; k < ks.dim(); ++k )
	dist2 += KnUtils::sqr( vcenter.ro( k ) - vp.ro( k ) );
      if ( dist2 <= r2 )
	s += p;
    }
  while ( scanner.increment( p, vp ) );
  return s;
}


void testRSet( int argc, char** argv )
{
  Kn_size sizes[ 100 ];
  sizes[ 0 ] = 128;
  sizes[ 1 ] = 128;
  sizes[ 2 ] = 128;
  sizes[ 3 ] = 8;

  int i;
  uint d = 3;
  uint r = 0;
  // dim "-d"
  for ( i = 1; ( i < argc ) && ( strcmp( argv[ i ], "-d" ) != 0 ); ++i ) ;
  if ( ( i + 1 ) < argc ) d = atoi( argv[ i + 1 ] );
  // radius "-r"
  for ( i = 1; ( i < argc ) && ( strcmp( argv[ i ], "-r" ) != 0 ); ++i ) ;
  if ( ( i + 1 ) < argc ) r = atoi( argv[ i + 1 ] );
  // dim "-x"
  for ( i = 1; ( i < argc ) && ( strcmp( argv[ i ], "-x" ) != 0 ); ++i ) ;
  if ( ( i + 1 ) < argc ) sizes[ 0 ] = atoi( argv[ i + 1 ] );
  // dim "-y"
  for ( i = 1; ( i < argc ) && ( strcmp( argv[ i ], "-y" ) != 0 ); ++i ) ;
  if ( ( i + 1 ) < argc ) sizes[ 1 ] = atoi( argv[ i + 1 ] );
  // dim "-z"
  for ( i = 1; ( i < argc ) && ( strcmp( argv[ i ], "-z" ) != 0 ); ++i ) ;
  if ( ( i + 1 ) < argc ) sizes[ 2 ] = atoi( argv[ i + 1 ] );
  // dim "-t"
  for ( i = 1; ( i < argc ) && ( strcmp( argv[ i ], "-t" ) != 0 ); ++i ) ;
  if ( ( i + 1 ) < argc ) sizes[ 3 ] = atoi( argv[ i + 1 ] );

  KnSpace ks( d, sizes );
    
  // center
  Kn_uid center = ks.uspel( ks.ukcode( sizes ) );
  cout << "--- Space: " << ks << endl;
  cout << "--- creating sphere (r=" << r << ") -----" << endl;
  
  KnRCellSet s = umakeVolumicSphere2( ks, center, (float) r );
  cout << s << endl;
  KnRCellSet::cell_iterator ip = s.begin();
  KnRCellSet::cell_iterator iq = s.end();
  cout << "first=";
  ks.displayKn_uid( *ip, cout );
  cout << endl;
  uint nb_s = 0;
  while ( ip != iq ) 
    {
      ++nb_s;
      ++ip;
    }
  
  KnRCellSet invs = ~s;
  cout << invs << endl;

  ip = invs.begin();
  iq = invs.end();
  uint nb_invs = 0;
  cout << "first inv=";
  ks.displayKn_uid( *ip, cout );
  cout << endl;
  while ( ip != iq ) 
    {
      ks.displayKn_uid( *ip, cout );
      cout << endl;
      ++nb_invs;
      ++ip;
    }

  cout << nb_s << " " << nb_invs << " " << nb_s + nb_invs << " "
       << ks.size() << endl;
  

  // this is the surfel set that we intend to create.
  KnRCellSet bdry = KnRCellSet::create( ks, d - 1, false, 0 );

  // Filling it by scanning the bounding cube.
  bool in_here, in_further;
  uint k;
  for ( k = 0; k < ks.dim(); ++k )
    {
      Kn_uid dir_low_uid = ks.uspel( ks.ufirst() );
      Kn_uid dir_up_uid = ks.ugetDecr( ks.uspel( ks.ulast() ), k );
      
      Kn_uid p = dir_low_uid;
      do 
 	{
 	  in_here = s[ p ];
 	  in_further = s[ ks.ugetIncr( p, k ) ];
 	  if ( in_here != in_further ) // boundary element
 	    { // add it to the set.
	      bdry += ks.uincident( p, k, true );
 	    }
 	}
      while ( ks.uincrInBounds( p, dir_low_uid, dir_up_uid ) );
    }
  cout << bdry << endl;

  KnCharSet cbdry = KnRCellSet::makeKnCharSet( bdry );
  ip = bdry.begin();
  iq = bdry.end();
  KnCharSet::cell_iterator cp = cbdry.begin();
  KnCharSet::cell_iterator cq = cbdry.end();
  uint nb_ok = 0;
  while ( ip != iq )
    {
      if ( *ip != *cp )
	cout << "ERROR ! Differs: " << *ip << " != " << *cp << endl;
      else 
	nb_ok++;
      ++ip;
      ++cp;
    }
  if ( cp != cq ) cout << "ERROR at end of scan" << endl;
  cout << "Comparison: " << nb_ok << " / rSet=" << bdry.nbElements() 
       << " CharSet=" << cbdry.nbElements() << endl;
  //cout << cbdry << endl;

}

int main( int argc, char** argv )
{
  testRSet( argc, argv );
  return 0;
}


