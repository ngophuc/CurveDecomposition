/**
 * @file   test_Freeman.cxx
 * @author Jacques-Olivier Lachaud
 * @date   Fri Apr 22 13:58:37 2005
 * 
 * @brief  Test module for geometry on 4-connected Freeman chains.
 * 
 * 
 */
///////////////////////////////////////////////////////////////////////////////
// Module de test : chemin 4-connexe par code de Freeman
///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <cassert>
#include <string>
#include "ImaGene/base/HashTable.h"
#include "ImaGene/dgeometry2d/C4CIterator.h"
#include "ImaGene/dgeometry2d/C4CGeometry.h"
#include "ImaGene/dgeometry2d/C4CSegment.h"
#include "ImaGene/dgeometry2d/C4CSegmentPencil.h"
#include "ImaGene/dgeometry2d/FreemanFrame2D.h"
#include "ImaGene/digitalnD/KnTypes.h"
#include "ImaGene/digitalnD/KnSpace.h"

using namespace std;
using namespace ImaGene;

// Freeman move : 0=east 1=north, 2=west, 3=south
const string circle = "1121212223232333030300010101";
const string rectangle = "1112222222233333300000000111";
const string lozenge = "1212121232323230303030101010";


// Freeman chain turns around object counterclockwise.
struct FreemanChain 
{
  string m_codes;
  inline void setCodes( const string & s )
  {
    m_codes = s;
  }
  inline uint code( uint pos ) const
  {
    assert( pos < m_codes.size() );
    return m_codes[ pos ] - '0';
  }
  inline uint next( uint pos ) const
  {
    ++pos;
    if ( pos >= m_codes.size() )
      pos = 0;
    return pos;
  }
  inline uint previous( uint pos ) const
  {
    if ( pos == 0 ) pos = m_codes.size() - 1;
    else --pos;
    return pos;
  }
  
};

// Freeman forward movement table (source*4 + dest). From two
// consecutive directions (source, dest), returns 0 if the movement
// was impossible, 1 if it was a turn toward the inside, 2 if it was a
// go straight, 3 if it was a turn toward the outside, 4 if it was a
// backtrack.
static const uint FFMTable[ 16 ] =
  { 2, 1, 0, 3,
    3, 2, 1, 0,
    0, 3, 2, 1,
    1, 0, 3, 2 
  };
// Ie ( 6 + source - dest ) mod 4

// Freeman backward movement table (source*4 + dest). From two
// consecutive directions (source, dest), returns 0 if the movement
// was impossible, 1 if it was a turn toward the inside, 2 if it was a
// go straight, 3 if it was a turn toward the outside, 4 if it was a
// backtrack.
static const uint FBMTable[ 16 ] =
  { 2, 3, 0, 1,
    1, 2, 3, 0,
    0, 1, 2, 3,
    3, 0, 1, 2
  };
// Ie ( 6 - source + dest ) mod 4



namespace ImaGene 
{

  /**
   * An exemple of iterator on a 4-connected path : an iterator on a
   * Freeman chain.
   */
  class C4CIteratorOnFreemanChain : public C4CIterator
  {
  private:

    const FreemanChain* m_chain;
    uint m_index;

  public:
    /**
     * Default Constructor.
     * The object is not valid.
     */
    inline C4CIteratorOnFreemanChain()
      : m_chain( 0 ), m_index( 0 )
    {
    }
  
    /**
     * Initializes the iterator as following the sequence of Freeman
     * moves [s] with initial index position [pos].
     *
     * @param s any non-null string composed of '0', '1', '2', '3',
     * possibly '4' for open digital curves.
     * 
     * @param pos an index smaller than the size of [s].
     */
    inline void init( const FreemanChain & s, uint pos )
    {
      m_chain = &s;
      m_index = pos;
      assert( m_chain->code( m_index ) <= 3 && "Non valid index in chain.");
    }


    /**
     * @return the current Freeman move.
     */
    inline uint currentCode() const
    {
      return m_chain->code( m_index );
    }

    /** 
     * @return the current index in the Freeman chain.
     */
    inline uint currentIndex() const
    {
      return m_index;
    }
    
    /**
     * @return a clone of the current iterator.
     */
    inline virtual C4CIterator* clone() const
    {
      C4CIteratorOnFreemanChain* it = new C4CIteratorOnFreemanChain;
      it->init( *m_chain, m_index );
      return it;
    }
  

    /**
     * @param other any other iterator.
     * @return 'true' if [other] points to the same location as 'this'.
     */
    inline virtual bool equals( const C4CIterator & other ) const
    {
      // NB: iterators must lie on the same chain.
      const C4CIteratorOnFreemanChain & cother = 
	dynamic_cast<const C4CIteratorOnFreemanChain &>( other );
      return ( m_index == cother.m_index ) && ( m_chain == cother.m_chain );
    }
  

    /**
     * Moves the iterator on the 4-connected contour to the next position.
     * @return 0 if the move was impossible, 1 if it was a move toward the interior, 2 if it was a straight movement, 3 if it was a move toward the exterior.
     * NB: If the C4CIterator is moving over a digital surface, then 'next'
     * means moving along a direct tracking direction.
     */
    inline virtual uint next()
    {
      uint ccur = m_chain->code( m_index );
      uint cnext = m_chain->code( m_chain->next( m_index ) );
      if ( cnext > 3 ) return 0;
      m_index = m_chain->next( m_index );
      return ( 6 + ccur - cnext ) % 4; // cf FFMTable
    }
  

    /**
     * Moves the iterator on the 4-connected contour to the previous position.
     * @return 0 if the move was impossible, 1 if it was a move toward the interior, 2 if it was a straight movement, 3 if it was a move toward the exterior.
     * NB: If the C4CIterator is moving over a digital surface, then 
     * 'previous' means moving along an indirect tracking direction.
     */
    inline virtual uint previous()
    {
      uint ccur = m_chain->code( m_index );
      uint cprev = m_chain->code( m_chain->previous( m_index ) );
      if ( cprev > 3 ) return 0;
      m_index = m_chain->previous( m_index );
      return ( 6 - ccur + cprev ) % 4; // cf FBMTable
    }



  };
  
}


/** 
 * Computes an estimated tangent orientation to a shape boundary (lMST
 * estimator).
 * 
 * @param input_it the iterator on the shape boundary.
 * @param input_starting_point the coordinates of the starting iterator.
 */
void lMSTangentEstimation( const C4CIteratorOnFreemanChain & input_it,
			   const Vector2i & input_starting_point )
{
  cout << "# angle-direction ds_corrected ds_averaged x_surfel"
       << " y_surfel polar_angle_surfel idx_surfel"
       << endl;

  // Functions used for combining all the angle of the maximal
  // segments inside the pencil.
  TriangleFunction l;
  DTriangleFunction lp;
  
  FreemanFrame2D frame; // local frame
  Vector2i p( input_starting_point );
  C4CIteratorOnFreemanChain cur_it( input_it );
  // A pencil has experimentally no more than 7 segments.
  const uint m = 10;
  C4CSegment segments[ m ];
  uint idx = 0;
  float length_corrected = 0.0;
  float length_averaged = 0.0;
  do
    {
      // Builds the pencil of maximal segments.
      uint j = 0;
      uint k;
      if ( C4CGeometry::maximalSegments( cur_it, segments, j, k, m ) )
	{
	  // All geometric computations are made in the local frame of
	  // the current boundary element.
	  C4CSegmentPencil pencil( segments, j, k, m, l, lp );
	  float theta = pencil.angleToX( Vector2D( 0.5, 0.0 ) );
	  // Cast angle in the global frame.
	  frame.init( cur_it.currentCode() );
	  cout << " " << frame.angleToX( theta );
	  // Compute elementary length
	  float ds_corrected = cos( theta );
	  float ds_averaged = 1.0
	    / ( fabs( cos( theta ) ) + fabs( sin( theta ) ) );
	  cout << " " << ds_corrected << " " << ds_averaged;
	  length_corrected += ds_corrected;
	  length_averaged += ds_averaged;
	}

      // Computes the next base point on the contour.
      Vector2i q( p );
      q.move4( cur_it.currentCode() );

      // Displays some coordinates.
      Vector2D mid( ( p.x() + q.x() ) / 2.0, ( p.y() + q.y() ) / 2.0 );
      // t : polar angle of point of interest.
      float t = atan2( mid.ro( 1 ), mid.ro( 0 ) );
      if ( t < 0.0 ) t += 2*M_PI;
      cout << " " << mid.ro( 0 )
	   << " " << mid.ro( 1 )
	   << " " << t
	   << " " << idx;

      // Go to next element.
      cout << endl;
      p = q;
      ++idx;
      if ( cur_it.next() == 0) break;
    }
  while ( ! cur_it.equals( input_it ) );

  cout << "# length_corrected=" << length_corrected << endl;
  cout << "# length_averaged=" << length_averaged << endl;
  cout << "# length_discrete=" << idx << endl;
}


/** 
 * Computes an estimated tangent orientation to a shape boundary
 * (symmetric tangent estimator.
 * 
 * @param input_it the iterator on the shape boundary.
 * @param input_starting_point the coordinates of the starting iterator.
 */
void STangentEstimation( const C4CIteratorOnFreemanChain & input_it,
			 const Vector2i & input_starting_point )
{
  cout << "# angle-direction ds_corrected ds_averaged x_surfel"
       << " y_surfel polar_angle_surfel idx_surfel idx_front idx_back"
       << endl;
  
  FreemanFrame2D frame; // local frame
  Vector2i p( input_starting_point );
  C4CIteratorOnFreemanChain cur_it( input_it );
  uint idx = 0;
  float length_corrected = 0.0;
  float length_averaged = 0.0;
  do
    {
      // Builds symmetric tangent.
      C4CIteratorOnFreemanChain front_it( cur_it );
      C4CIteratorOnFreemanChain back_it( cur_it );
      C4CSegment symtgt = 
	C4CGeometry::symmetricTangent( front_it, // frontmost element
				       back_it,  // backmost element
				       0 // max_size 0: no
				       );
      Vector2i tgtv( symtgt.getTangent() );
      float theta = atan( (float) tgtv.y() / (float) tgtv.x() );
      // Cast angle in the global frame.
      frame.init( cur_it.currentCode() );
      cout << " " << frame.angleToX( theta );
      float ds_corrected = symtgt.correctedLength();
      float ds_averaged = symtgt.averagedLength();
      cout << " " << ds_corrected << " " << ds_averaged;
      length_corrected += ds_corrected;
      length_averaged += ds_averaged;

      // Computes the next base point on the contour.
      Vector2i q( p );
      q.move4( cur_it.currentCode() );

      // Displays some coordinates.
      Vector2D mid( ( p.x() + q.x() ) / 2.0, ( p.y() + q.y() ) / 2.0 );
      // t : polar angle of point of interest.
      float t = atan2( mid.ro( 1 ), mid.ro( 0 ) );
      if ( t < 0.0 ) t += 2*M_PI;
      cout << " " << mid.ro( 0 )
	   << " " << mid.ro( 1 )
	   << " " << t;

      // Displays indices of the symmetric tangent extremeties.
      cout << " " << cur_it.currentIndex() 
	   << " " << front_it.currentIndex() 
	   << " " << back_it.currentIndex();
      
      // Go to next element.
      cout << endl;
      p = q;
      ++idx;
      if ( cur_it.next() == 0) break;
    }
  while ( ! cur_it.equals( input_it ) );

  cout << "# length_corrected=" << length_corrected << endl;
  cout << "# length_averaged=" << length_averaged << endl;
  cout << "# length_discrete=" << idx << endl;
}


/** 
 * Code a surfel (ie linel in 2d) from an iterator on a Freeman chain
 * and the location of the pointel that is the base of the Freeman code.
 * 
 * @param ks the digital space
 * @param it the iterator on a Freeman chain
 * @param base_pt the position of the pointel that is the base of the Freeman code.
 * 
 * @return an oriented surfel.
 *
 * @todo Should not be done like that. A new class of 'C4CIterator'
 * should be defined to represent iterators on Freeman chains that are
 * embedded in a digital space, ie an Iterator that derives from
 * 'C4CIteratorOnSurface'.
 */
Kn_sid
surfelFromFreeman( const KnSpace & ks,
		   const C4CIteratorOnFreemanChain & it,
		   const Vector2i & base_pt )
{
  Kn_size xy[ 2 ];
  // Khalimsky coordinates of pixel.
  xy[ 0 ] = 2*base_pt.x() + 1;
  xy[ 1 ] = 2*base_pt.y() + 1;
  Kn_sid pixel;
  Kn_sid bel;
  switch ( it.currentCode() )
    {
    case 0: 
      pixel = ks.skcode( xy, KnTypes::POS );
      bel = ks.sincident( pixel, 1, false );
      break;
    case 1:
      xy[ 0 ] -= 2;
      pixel = ks.skcode( xy, KnTypes::POS );
      bel = ks.sincident( pixel, 0, true );
      break;
    case 2:
      xy[ 0 ] -= 2;
      xy[ 1 ] -= 2;
      pixel = ks.skcode( xy, KnTypes::POS );
      bel = ks.sincident( pixel, 1, true );
      break;
    case 3:
      xy[ 1 ] -= 2;
      pixel = ks.skcode( xy, KnTypes::POS );
      bel = ks.sincident( pixel, 0, false );
      break;
    default:
      bel = 0;
    }
  return bel;
}



/**
 * Energy: surfel -> float number
 */
typedef HashTable<float> EnergyFct;


float
computeStretchingEnergy( const KnSpace & ks,
			 const C4CIteratorOnFreemanChain & input_it,
			 const Vector2i & input_starting_point,
			 EnergyFct & energy )
{
  float engtotal = 0.0;
  FreemanFrame2D frame; // local frame
  C4CIteratorOnFreemanChain cur_it( input_it );
  Vector2i p( input_starting_point );
  uint idx = 0;
  do
    {
      // Builds symmetric tangent.
      C4CIteratorOnFreemanChain front_it( cur_it );
      C4CIteratorOnFreemanChain back_it( cur_it );
      C4CSegment symtgt = 
	C4CGeometry::symmetricTangent( front_it, // frontmost element
				       back_it,  // backmost element
				       0 // max_size 0: no
				       );

      // Choose an elementary length estimator.
      float ds = symtgt.correctedLength();
      // float ds = symtgt.averagedLength();

      // Attribute energy to surfel.
      Kn_sid bel = surfelFromFreeman( ks, cur_it, p );
      energy[ bel ] = ds;
      engtotal += ds;
      cout << " " << ds;
      // Computes the next base point on the contour.
      p.move4( cur_it.currentCode() );

      // Go to next element.
      ++idx;
      if ( cur_it.next() == 0) break;
    }
  while ( ! cur_it.equals( input_it ) );
  return engtotal;
}

class DummyImage
{
public:
  inline float value( uint x, uint y ) const
  {
    float d = sqrt( (double) (x-32)*(x-32) + (y-32)*(y-32) );
    if ( d < 3.5 ) return 1.0;
    else return 0.0;
  }
  
};


typedef DummyImage ImageFct;

float
computeImageEnergy( const KnSpace & ks,
		    const C4CIteratorOnFreemanChain & input_it,
		    const Vector2i & input_starting_point,
		    EnergyFct & energy,
		    ImageFct & image )
{
  float engtotal = 0.0;
  FreemanFrame2D frame; // local frame
  C4CIteratorOnFreemanChain cur_it( input_it );
  Vector2i p( input_starting_point );
  uint idx = 0;
  do
    {
      // Compute surfel and inner and outer voxels.
      Kn_sid bel = surfelFromFreeman( ks, cur_it, p );
      uint j = ks.sorthDir( bel );
      Kn_sid inner_voxel = ks.sincident( bel, j, ks.sdirect( bel, j ) );
      Kn_sid outer_voxel = ks.sincident( bel, j, ! ks.sdirect( bel, j ) );
      Kn_size xy[ 2 ];
      ks.sdecodeCoords( inner_voxel, xy );
      float img_inside = image.value( xy[ 0 ], xy[ 1 ] );
      ks.sdecodeCoords( outer_voxel, xy );
      float img_outside = image.value( xy[ 0 ], xy[ 1 ] );
      
      // Attribute energy to surfel.
      float e = - ( img_inside - img_outside ) * (img_inside - img_outside );
      energy[ bel ] = e;
      engtotal += e;
      cout << " " << e;
      // Computes the next base point on the contour.
      p.move4( cur_it.currentCode() );

      // Go to next element.
      ++idx;
      if ( cur_it.next() == 0) break;
    }
  while ( ! cur_it.equals( input_it ) );
  return engtotal;
}


void
testEnergy( const KnSpace & ks, 
	    C4CIteratorOnFreemanChain & input_it,
	    const Vector2i & input_starting_point )
{
  bool is_open;
  uint nb_surfel = C4CIterator::size( input_it, is_open );

  EnergyFct eng_stretch;
  eng_stretch.init( ks.slast(), // maximal code
		    nb_surfel   // expected number of values
		    );
  
  float total_eng_stretch =
    computeStretchingEnergy( ks,
			     input_it,
			     input_starting_point,
			     eng_stretch );
  cout << endl
       << "Stretching energy = " << total_eng_stretch
       << " for " << nb_surfel << " elements, avg = "
       << ( total_eng_stretch / nb_surfel ) << endl;

  EnergyFct eng_image;
  eng_image.init( ks.slast(), // maximal code
		    nb_surfel   // expected number of values
		    );
  DummyImage img;
  float total_eng_image =
    computeImageEnergy( ks,
			input_it,
			input_starting_point,
			eng_image,
			img);
  cout << endl
       << "Image energy = " << total_eng_image
       << " for " << nb_surfel << " elements, avg = "
       << ( total_eng_image / nb_surfel ) << endl;
}



int main( int argc, char** argv )
{
  FreemanChain chain;
  C4CIteratorOnFreemanChain start_it;

  // chain.setCodes( lozenge );
  chain.setCodes( circle );
  start_it.init( chain, 0 );

  // Test lambda-Maximal Segment Tangent estimator.
  lMSTangentEstimation( start_it, Vector2i( 4, 0 ) );

	exit(0);

  // Test Symmetric Tangent estimator.
  // STangentEstimation( start_it, Vector2i( 4, 0 ) );
  uint sizes[ 2 ] = { 64, 64 };
  KnSpace ks( 2, sizes );
  
  cout << " ----- circle ----- " << endl;
  chain.setCodes( circle );
  start_it.init( chain, 0 );
  testEnergy( ks, start_it, Vector2i( 36, 32 ) );

  cout << " ----- lozenge ----- " << endl;
  chain.setCodes( lozenge );
  start_it.init( chain, 0 );
  testEnergy( ks, start_it, Vector2i( 36, 32 ) );

  cout << " ----- rectangle ----- " << endl;
  chain.setCodes( rectangle );
  start_it.init( chain, 0 );
  testEnergy( ks, start_it, Vector2i( 36, 32 ) );
}


