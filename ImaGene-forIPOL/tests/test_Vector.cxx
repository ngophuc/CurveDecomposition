///////////////////////////////////////////////////////////////////////////////
// Test module for Vector
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <cstdlib>

//#include "LinAlg/LinAlg2D/Vector2D.hpp"
#include "ImaGene/base/Vector.h"

using namespace std;
using namespace ImaGene;

int 
main( int argc, char** argv ) 
{
  Vector2D x( 1.0, -3.0 );
  cout << "x = " << x << " <= 1 -3" << endl;
  Vector2D y( x );
  cout << "y = " << y << " <= 1 -3" << endl;
  y += x;
  cout << "y = " << y << " <= 2 -6" << endl;
  y *= 3.0;
  cout << "y = " << y << " <= 6 -18" << endl;
  x = y;
  cout << "x = " << x << " <= 6 -18" << endl;
  return 0;
}
