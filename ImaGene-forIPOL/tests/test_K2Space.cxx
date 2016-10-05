///////////////////////////////////////////////////////////////////////////////
// Test class K2Space.
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/dgeometry2d/FreemanChain.h"
#include "ImaGene/dgeometry2d/C4CIteratorOnFreemanChain.h"
#include "ImaGene/digitalnD/GridEmbedder.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/K2Space.h"
#include "ImaGene/helper/ShapeHelper.h"

using namespace std;
using namespace ImaGene;


static Arguments args;

int test_K2Space_slinel( K2Space & K2, 
			 KnRCellSet & contour, 
			 FreemanChain & c )
{
  FreemanChain::const_iterator it;
  uint nb_ok = 0;
  uint nb_ko = 0;
  uint timing = args.getOption( "-timing" )->getIntValue( 0 );
  uint trace = args.getOption( "-trace" )->getIntValue( 0 );
  long t1;
  if ( trace >= 3 )
    {
      cerr << "---- contour from makeContourFromFreemanChain ----" << endl;
      for ( KnRCellSet::cell_iterator p = contour.begin(); 
	    p != contour.end(); 
	    ++p )
	{
	  K2.displayKn_sid( *p, cerr );
	  cerr << endl;
	}
    }
  if ( trace >= 3 )
    cerr << "---- contour from K2Space::slinel ----" << endl;
  if ( timing > 0 )
    Clock::startClock();
  for ( it = c.begin(); it != c.end(); it.next() )
    {
      uint code = it.getCode();
      Vector2i base( *it );
      Kn_sid s = K2.slinel( base.x(), base.y(), code );
      if ( trace >= 3 )
	{
	  K2.displayKn_sid( s, cerr );
	  cerr << endl;
	}
      if ( contour[ s ] )
	++nb_ok;
      else 
	++nb_ko;
    }
  if ( timing > 0 )
    {
      long t1 = Clock::stopClock();
      cerr << "[test_K2Space_slinel] in " << t1 << " ms."
	   << endl;
    }
  if ( trace >= 1 )
    {
      cerr << "[test_K2Space::test<slinel>] passed " << nb_ok << "/" 
	   << ( nb_ok + nb_ko ) << endl;
    }
  if ( trace >= 2 )
    cerr << "\tThis test computes the linels of a Freeman chain code with optimized K2Space methods and compared the result with ShapeHelper::makeContourFromFreemanChain." << endl;
  return nb_ko == 0 ? 0 : 1;
}

int main( int argc, char** argv )
{
  // -------------------------------------------------------------------------
  // Prepare arguments.
  StandardArguments::addDebugArgs( args, true, true, true );
  StandardArguments::addIOArgs( args, true, false );
  args.addBooleanOption( "-test_slinel", "-test_slinel: test K2Space::slinel." );
  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "test_K2Space", 
			  "Test class K2Space.",
			  "" )
	   << endl;
      return 1;
    }
  uint trace = args.getOption( "-trace" )->getIntValue( 0 );
  bool debug = args.check( "-debug" );
  uint timing = args.getOption( "-timing" )->getIntValue( 0 );
  // -------------------------------------------------------------------------
  // Read Freeman chain and creates space/contour.

  FreemanChain c;
  //ifstream in_str( args.getOption( "-input" )->getValue( 0 ).c_str() );
  istream & in_str = StandardArguments::openInput( args );
  FreemanChain::read( in_str, c );
  if ( ! in_str.good() )
    {
      if ( debug )
	cerr << "Error reading Freeman chain code." << endl;
      return 2;
    }

  KnSpace* ks = ShapeHelper::makeSpaceFromFreemanChain( c );
  if ( timing > 0 )
    Clock::startClock();
  KnRCellSet contour = ShapeHelper::makeContourFromFreemanChain( ks, c, true );
  if ( timing > 0 )
    {
      long t1 = Clock::stopClock();
      cerr << "[ShapeHelper::makeContourFromFreemanChain] in " << t1 << " ms."
	   << endl;
    }
  K2Space K2( ks->size( 0 ), ks->size( 1 ) );

  int err = 0;
  if ( args.check( "-test_slinel" ) )
    {
      err += test_K2Space_slinel( K2, contour, c );
    }
  return err;
}
