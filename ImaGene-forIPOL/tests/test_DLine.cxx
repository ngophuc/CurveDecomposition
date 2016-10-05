///////////////////////////////////////////////////////////////////////////////
// Test module for digitalnD and Kn_uid, Kn_sid
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <cstdlib>

#include "ImaGene/dgeometry2d/DLine.h"
//#include "ImaGene/timetools/Clock.h"

using namespace std;
using namespace ImaGene;



int main( int argc, char** argv ) 
{
  DLine l1( 10, 1, 0 );
  DLine l2( -10, 1, -10 );
  cout << "INPUT: " << l1 << l2 << endl;
  DLine m1 = DLine::medianLine( l1, l2 );
  DLine m2 = DLine::medianLine( l2, l1 );
  cout << "MEDIAN: " << m1 << m2 << endl;
  
  return 0;
  
}
