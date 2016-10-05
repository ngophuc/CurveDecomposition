///////////////////////////////////////////////////////////////////////////////
// Test module for digitalnD and Kn_uid, Kn_sid
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <cstdlib>

#include "ImaGene/dgeometry2d/C4CSegment.h"
#include "ImaGene/timetools/Clock.h"

#include "ImaGene/dgeometry2d/C4CIterator.h"
#include "ImaGene/dgeometry2d/C4CGeometry.h"


using namespace std;
using namespace ImaGene;
using namespace ImaGene;



#define CIRCLE
//#define SQUARE

#if defined(CIRCLE)
#define NAME "circle"
#define RADIUS 9
#define NB 72
static int cercle[ NB ] = 
  {
    1, 3, 2, 2, 1, 3, 2, 1, 3, 1, 3, 2, 3, 1, 2, 2, 3, 1,
    1, 3, 2, 2, 1, 3, 2, 1, 3, 1, 3, 2, 3, 1, 2, 2, 3, 1,
    1, 3, 2, 2, 1, 3, 2, 1, 3, 1, 3, 2, 3, 1, 2, 2, 3, 1,
    1, 3, 2, 2, 1, 3, 2, 1, 3, 1, 3, 2, 3, 1, 2, 2, 3, 1 
  };

#define PREVIOUS( x ) x = ( x + NB - 1 ) % NB
#define NEXT( x ) x = ( x + 1 ) % NB

#define FWD_VALUE( x ) cercle[ ( x + 1 ) % NB ]
#define BWD_VALUE( x ) cercle[ x ] 
#endif

#if defined(SQUARE)
#define NAME "square"
#define NB 1000
#define RADIUS NB

#define PREVIOUS( x ) x = ( x + NB - 1 ) % NB
#define NEXT( x ) x = ( x + 1 ) % NB

#define FWD_VALUE( x ) ( ( x + 1 ) % NB == 0 ) ? 1 : 2
#define BWD_VALUE( x ) ( x == 0 ) ? 1 : 2
#endif

typedef int Iterator;

void test_tangent_sym( int nb )
{

  for ( int i = 0; i < nb; ++i )
    {
      
      C4CSegment segment;
      Iterator begin = 0;
      Iterator current = begin;
  
      do 
	{
	  // Computes symmetric tangent.
	  Iterator fwd = current;
	  Iterator bwd = current;
	  segment.init();
      
	  while ( segment.extends( BWD_VALUE( bwd ), FWD_VALUE( fwd ) ) )
	    {
	      PREVIOUS( bwd );
	      NEXT( fwd );
	    }

	  //cout << "[ " << current << " ] : " << segment << endl;

	  NEXT( current );
	}
      while ( current != begin );
    }
}

void test_tangent_sym_opt( int nb )
{
  C4CSegment segment;
  Iterator begin = 0;
  Iterator current = begin;
  segment.init();
  Iterator fwd = current;
  Iterator bwd = current;
  
  for ( int i = 0; i < nb; ++ i )
    {
      
      do 
	{
	  // Computes symmetric tangent.
      
	  while ( segment.extends( BWD_VALUE( bwd ), FWD_VALUE( fwd ) ) )
	    {
	      PREVIOUS( bwd );
	      NEXT( fwd );
	    }

	  //cout << "[ " << current << " ] : " << segment << endl;

	  // Slides to the next one.
	  if ( segment.isBackRetractable() )
	    {
	      segment.retractsBack( FWD_VALUE( bwd ) );
	      //cout << "      R1: " << segment << endl;
	      NEXT( bwd );
	      segment.slidesForward( FWD_VALUE( current ) );
	      //cout << "     Slid: " << segment << endl;
	      NEXT( current );
	      if ( ! segment.isBackRetractable() )
		cout << "      ERROR: " << segment << endl;
	      segment.retractsBack( FWD_VALUE( bwd ) );
	      //cout << "      R2: " << segment << endl;
	      NEXT( bwd );
	    }
	  else
	    {
	      NEXT( bwd );
	      NEXT( current );
	      NEXT( fwd );
	    }
	}
      while ( current != begin );
    }
}

void test_tangent_sym_opt2( int nb )
{
  C4CSegment segment;
  Iterator begin = 0;
  Iterator current = begin;
  segment.init();
  Iterator fwd = current;
  Iterator bwd = current;
  
  for ( int i = 0; i < nb; ++ i )
    {
      
      do 
	{
	  // Computes symmetric tangent.
      
	  while ( segment.extends( BWD_VALUE( bwd ), FWD_VALUE( fwd ) ) )
	    {
	      PREVIOUS( bwd );
	      NEXT( fwd );
	    }

	  //cout << "[ " << current << " ] : " << segment << endl;

	  // Slides to the next one.
	  if ( segment.isBackRetractable() )
	    {
	      segment.retractsBack( FWD_VALUE( bwd ) );
	      //cout << "      R1: " << segment << endl;
	      NEXT( bwd );
	      segment.slidesForward( FWD_VALUE( current ) );
	      //cout << "     Slid: " << segment << endl;
	      NEXT( current );

	      if ( segment.extendsFront( FWD_VALUE( fwd ) ) != 0 )
		{
		  if ( ! segment.isBackRetractable() )
		    cout << "      ERROR: " << segment << endl;
		  segment.retractsBack( FWD_VALUE( bwd ) );
		  //cout << "      R2: " << segment << endl;
		  NEXT( bwd );
		}
	      else NEXT( fwd );
	    }
	  else
	    {
	      NEXT( bwd );
	      NEXT( current );
	      NEXT( fwd );
	    }
	}
      while ( current != begin );
    }
}

class TC4Iterator : public C4CIterator
{
public:
  TC4Iterator( int idx ) : m_idx( idx ) {}
  TC4Iterator( const TC4Iterator & other ) : m_idx( other.m_idx ) {}
  TC4Iterator& operator=( const TC4Iterator & other )
  {
    m_idx = other.m_idx;
    return *this;
  }
  
  virtual ~TC4Iterator() {}
  virtual C4CIterator* clone() const 
  {
    return new TC4Iterator( m_idx );
  }
  virtual bool equals( const C4CIterator & other ) const 
  {
    const TC4Iterator* it = dynamic_cast<const TC4Iterator*>( &other );
    return ( it != 0 ) && ( m_idx == it->m_idx );
  }

  uint current() const
  {
    return m_idx;
  }
  virtual uint next()
  {
    uint code = FWD_VALUE( m_idx );
    NEXT( m_idx );
    return code;
  }
  
  virtual uint previous()
  {
    uint code = BWD_VALUE( m_idx );
    PREVIOUS( m_idx );
    return code;
  }
  
private:
  int m_idx;
};

void test_curvature( int nb )
{
  for ( int i = 0; i < nb; ++i )
    {
      TC4Iterator begin( 0 );
      TC4Iterator curf = begin;
      TC4Iterator curb = begin;
      do 
	{
	  TC4Iterator copyf( curf );
	  TC4Iterator copyb( curb );
	  float c = C4CGeometry::curvatureByAngleVariation( copyf, copyb );
	  cout << "[" << curf.current() << "] K=" << c << endl;
	  curf.next();
	  curb.next();
	}
      while ( curf.current() != begin.current() );
    }
}


int main( int argc, char** argv ) 
{
  if ( argc <= 1 )
    {
      cerr << "Usage: " << argv[ 0 ] << " <nb> " << endl
	   << "       test the class DSegment with symmetric"
	   << " tangent computations." << endl
	   << "       nb: number of times the computation is made"
	   << " (for timings)." << endl;
      return 0;
    }

  int nb = atoi( argv[ 1 ] );

  cout << NAME << " " << RADIUS << " with " << NB << " elements." << endl;
  
  Clock::startClock();
  test_tangent_sym( nb );
  long t1 = Clock:: stopClock();
  cout << nb * NB << " sym. tangent computation in " << t1 << " ms." << endl;
  Clock::startClock();
  test_tangent_sym_opt( nb );
  long t2 = Clock:: stopClock();
  cout << nb * NB << " incr. sym. tangent computation in " << t2 << " ms." << endl;

  Clock::startClock();
  test_tangent_sym_opt2( nb );
  long t3 = Clock:: stopClock();
  cout << nb * NB << " incr2. sym. tangent computation in " << t3 << " ms." << endl;

  test_curvature( 1 );

  return 0;
}

 


