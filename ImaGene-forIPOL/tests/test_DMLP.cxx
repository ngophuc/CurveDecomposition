///////////////////////////////////////////////////////////////////////////////
// Test dynamic minimum length polygon
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/base/UndoableList.h"
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/mathutils/CFraction.h"
#include "ImaGene/dgeometry2d/FreemanChain.h"
#include "ImaGene/dgeometry2d/OctantChange.h"
#include "ImaGene/dgeometry2d/DMLPContour.h"

using namespace std;
using namespace ImaGene;

class DVector2i;
int detp( const DVector2i & v1, const DVector2i & v2 );

class DVector2i
{
protected:
  DVector2i()
  {
  }

public:
  INLINE DVector2i( int i, int j ) 
  {
    x = i; 
    y = j;
    regularize();
  }

  INLINE uint absx() const
  {
    return (uint) ( x >= 0 ? x : -x );
  }
  INLINE uint absy() const
  {
    return (uint) ( y >= 0 ? y : -y );
  }

  INLINE bool isTrivial() const
  {
    uint ax = absx();
    uint ay = absy();
    return ( ( ax == 1 ) && ( ay == 0 ) )
      || ( ( ax == 0 ) && ( ay == 1 ) );
  }

  INLINE int depth( ) const
  {
    return ((int) z.size()) - 1;
  }
  INLINE uint uk( int k ) const
  {
    return z[ k ];
  }

  INLINE uint pk( int k ) const
  {
    if ( k >= 0 ) return p[ k ];
    else if ( k == -1 ) return 1;
    else return 0;
  }

  INLINE uint qk( int k ) const
  {
    if ( k >= 0 ) return q[ k ];
    else if ( k == -1 ) return 0;
    else return 1;
  }

  DVector2i partial( int k ) const
  {
    DVector2i v;
    // TODO
  }

  void regularize()
  {
    uint a = absy();
    uint b = absx();
    p.clear();
    q.clear();
    z.clear();
    int k = 0;
    while ( b != 0 )
      {
	uint u = a / b;
	uint r = a - b * u;
	z.push_back( u );
	p.push_back( u * pk( k - 1 ) + pk( k - 2 ) );
	q.push_back( u * qk( k - 1 ) + qk( k - 2 ) );
	a = b;
	b = r;
	++k;
      }
    d = a;
    if ( d != 0 )
      {
	x /= d;
	y /= d;
      }
  }

  void selfDisplay( ostream & that_stream, int verbose = 0 )
  {
    that_stream << "(" << x << "," << y << ")^" << d;
    if ( ( verbose > 0 ) && ( depth() >= 0 ) )
      {
	that_stream << " [" << z[ 0 ];
	for ( uint k = 1; k < z.size(); ++k )
	  that_stream << "," << z[ k ];
	that_stream << "]";
      }
    if ( ( verbose > 1 ) && ( depth() >= 0 ) )
      {
	that_stream << " {" << p[ 0 ] << "/" << q[ 0 ];
	for ( uint k = 1; k < p.size(); ++k )
	  that_stream << " " << p[ k ] << "/" << q[ k ];
	that_stream << "}";
      }
  }

public:
  int x;
  int y;
  int d;
  std::vector<uint> p;
  std::vector<uint> q;
  std::vector<uint> z;
};

// determinant of primitive vectors.
int detp( const DVector2i & v1, const DVector2i & v2 )
{
  return v1.x * v2.y - v2.x * v1.y;
}


void affVector( const string & s, int x, int y )
{
  DVector2i v( x, y );
  cout << "- " << s << "=" << "(" << x << "," << y << ")"
       << " " << y << "/" << x << " ";
  v.selfDisplay( cout, 2 ); 
  cout << endl;
}

void testFractions() 
{
  affVector( "v1", 35, 55 );
  affVector( "v2", -28, 42 );
  affVector( "v3", 68, 21 );
  affVector( "v4", 10, 0 );
  affVector( "v5", 0, 0 );
  affVector( "v6", 0, 8 );

  CFraction c1( 5, 8 );
  CFraction c2( 3, 5 );
  CFraction c3( 2, 3 );
  CFraction c4 = c2.mediant( c3 );
  CFraction c5 = c1.mediant( c3 );
  CFraction c6 = c5.mediant( c3 );
  cout << c1 << endl
       << c2 << endl
       << c3 << endl
       << c4 << endl
       << c5 << endl
       << c6 << endl;
}

// class OctantChange {
//   bool sym_x_eq_0;
//   bool sym_y_eq_0;
//   bool sym_x_eq_y;
// public:
//   OctantChange( bool sx0 = false, bool sy0 = false, bool sxy = false )
//     : sym_x_eq_0( sx0 ), sym_y_eq_0( sy0 ), sym_x_eq_y( sxy )
//   {}
//   void selfDisplay( ostream & out ) const
//   {
//     out << "{OC" << sym_x_eq_0 << sym_y_eq_0 << sym_x_eq_y << "}";
//   }

//   void initByCastIntoFirstOctant( Vector2i v )
//   {
//     if ( v.x() < 0 )
//       {
// 	v.x() = -v.x();
// 	sym_x_eq_0 = true;
//       }
//     else
//       sym_x_eq_0 = false;

//     if ( v.y() < 0 )
//       {
// 	v.y() = -v.y();
// 	sym_y_eq_0 = true;
//       }
//     else
//       sym_y_eq_0 = false;

//     if ( v.x() < v.y() )
//       sym_x_eq_y = true;
//     else
//       sym_x_eq_y = false;
//   }

//   Vector2i cast( Vector2i u ) const
//   {
//     if ( sym_x_eq_0 ) u.x() = -u.x();
//     if ( sym_y_eq_0 ) u.y() = -u.y();
//     if ( sym_x_eq_y ) std:swap( u.x(), u.y() );
//     return u;
//   }

//   Vector2i castBack( Vector2i u ) const
//   {
//     if ( sym_x_eq_y ) std:swap( u.x(), u.y() );
//     if ( sym_y_eq_0 ) u.y() = -u.y();
//     if ( sym_x_eq_0 ) u.x() = -u.x();
//     return u;
//   }

//   bool isSymXY() const
//   {
//     return sym_x_eq_y;
//   }

//   bool isDirect() const
//   {
//     if ( ( ( sym_x_eq_0 ) && ( sym_y_eq_0 ) )
// 	 || ( ( ! sym_x_eq_0 ) && ( ! sym_y_eq_0 ) )
// 	 )
//       return ! sym_x_eq_y;
//     else 
//       return sym_x_eq_y;
//   }
//   bool isIndirect() const
//   {
//     return ! isDirect();
//   }
// };

// ostream &
// operator<<( ostream & out, const OctantChange & n )
// {
//   n.selfDisplay( out );
//   return out;
// }


/**
 * Represents a digital contour that support some local modification
 * methods, like flip or bumps. It is able to represent the digital
 * vector as a set of edges whose length is minimal (with the
 * constraint that it must separate inside pixels from outside
 * pixels).
 */
class DynamicContour {

public:
  /**
   * Determines if the vertex lies on the pixels of the shape (INSIDE)
   * or outside.
   */
  enum TouchType { INSIDE, OUTSIDE };

  /**
   * The type of connection between two consecutive nodes. For an
   * INSIDE vertex in a clockwise contour, the determinant \c d between the
   * two vectors is:
   * - d < 0 : WEDGE
   * - d = 0 : FLAT
   * - d = 1 : AFLAT
   * - d > 1 : VEE
   */
  enum ConnectionType { WEDGE, FLAT, AFLAT, VEE };

  /**
   * RIGHT means to the right of the vector (a seen from above), LEFT
   * to the left.
   */
  enum SideType { RIGHT, LEFT };

  /**
   * A node is an elementary part of a dynamic contour. It is defined
   * by a displacement vector [v], whose base type [base] indicates
   * whether the source of the vector lies inside or outside the
   * shape. The vector is always primitive, but may be repeated [n]
   * times. The transformation to put [v] in the first octant is given
   * by [oc]. In this first octant, the continued fraction of the
   * slope of the vector is stored in [z].
   */
  struct Node {
    /**
     * The displacement vector.
     */
    Vector2i v;

    /**
     * The continued fraction of the slope of the vector [v] cast in the
     * first octant.
     */
    CFraction z;

    /**
     * The transformation to bring [v] in the first octant.
     */
    OctantChange oc;
    
    /**
     * The type of pixel at the source of the Node vector.
     */
    TouchType base;

    /**
     * The number of times this vector is repeated.
     */
    uint n;

    /**
     * Default constructor. The object is invalid.
     */
    Node() {}

    /**
     * Constructor. Other members are computed automatically.
     *
     * @param v0 an irreducible vector.
     * @param base0 the type of pixel at the source of the Node vector.
     * @param n0 the number of times this vector is repeated.
     */
    Node( Vector2i v0, TouchType base0, uint n0 )
      : v( v0 ), base( base0 ), n( n0 )
    {
      oc.initByCastIntoFirstOctant( v0 );
      Vector2i vo = oc.cast( v );
      z.init( vo.y(), vo.x() );
    }

    /**
     * Constructor. Other members are computed automatically.
     *
     * @param v0 an irreducible vector.
     * @param z0 the continued fraction of [v0] cast in the first octant.
     * @param oc0 the transformation to cast [v0] in the first octant.
     * @param base0 the type of pixel at the source of the Node vector.
     * @param n0 the number of times this vector is repeated.
     */
    Node( Vector2i v0, 
	  const CFraction & z0, 
	  const OctantChange & oc0, 
	  TouchType base0, 
	  uint n0 )
      : v( v0 ), z( z0 ), oc( oc0 ), base( base0 ), n( n0 )
    {
    }

    /**
     * Flip the source of the vector inside/out. Useful for dynamic
     * contour where the shape of interest is deformed by flipping
     * some border pixels inside/out.
     */
    void flip() 
    {
      base = base == INSIDE ? OUTSIDE : INSIDE;
    }

    /**
     * Displays the node on the given output stream.
     *
     * @param out the output stream.
     */
    void selfDisplay( ostream & out ) const
    {
      out << ( base == INSIDE ? 'x' : 'o' ) 
	  << "(" << v.x() << "," << v.y() << ")";
      if ( n > 1 )
	out << "^" << n;
    }

    /**
     * Displays the node on the given output stream, with upper case letters.
     *
     * @param out the output stream.
     */
    void selfDisplayUpperCase( ostream & out ) const
    {
      out << ( base == INSIDE ? 'X' : 'O' ) 
	  << "(" << v.x() << "," << v.y() << ")";
      if ( n > 1 )
	out << "^" << n;
    }
  };

public:
  typedef list<Node> list_of_nodes;
  typedef list<Node>::iterator iterator;
  typedef list<Node>::const_iterator const_iterator;


private:
  /**
   * The dynamic contour is essentially a list of nodes.
   */
  list_of_nodes contour;

  /** 
   * The dynamic contour has a starting position, specified by the
   * iterator [start_node], and with a position in the plane
   * [start_pos].
   */
  iterator start_node;

  /** 
   * The dynamic contour has a starting position, specified by the
   * iterator [start_node], and with a position in the plane
   * [start_pos].
   */
  Vector2i start_pos;

public:

  /**
   * 'true' if contour turns around inside clockwise (INSIDE to
   * RIGHT), 'false' otherwise (INSIDE TO LEFT).
   */
  bool cw;


public:
  /**
   * Default constructor. The contour is invalid
   */
  DynamicContour() 
    : cw( true )
  {}
  
  /**
   * @return the number of nodes of the contour.
   */
  uint size() const
  {
    return contour.size();
  }

  /**
   * @return 'true' is the contour has zero nodes.
   */
  bool null() const
  {
    return contour.begin() == contour.end();
  }


  /**
   * @return the position in the plane of the start Node.
   * @see start_node
   */
  Vector2i getStartPosition() const
  {
    return start_pos;
  }

  /**
   * Displays the contour on the output stream for debugging purposes.
   *
   * @param out the output stream.
   */
  void selfDisplay( ostream & out ) const
  {
    for ( const_iterator it = contour.begin();
	  it != contour.end();
	  ++it )
      {
	if ( it == start_node )
	  out << "{S" << start_pos.x() << "," << start_pos.y() << "}";
	(*it).selfDisplay( out );
      }
    out << endl;
  }

  /**
   * Displays the contour on the output stream for debugging
   * purposes. Emphasizes the given iterator [it] in upper case
   * letters.
   *
   * @param out the output stream.
   * @param it an iterator in the list of nodes of the contour.
   */
  void selfDisplayWithIterator( ostream & out, const iterator & it ) 
  {
    for ( iterator itc = contour.begin();
	  itc != contour.end();
	  ++itc )
      {
	if ( itc == start_node )
	  out << "{S" << start_pos.x() << "," << start_pos.y() << "}";
	if ( it == itc) (*itc).selfDisplayUpperCase( out );
	else            (*itc).selfDisplay( out );
      }
    out << endl;
  }

  /**
   * Given a circulator [it] on the contour, returns the one on the
   * next element.
   * @param it an iterator in the list of nodes [contour].
   * @return the next one (assuming the contour is a loop).
   */
  iterator next( const iterator & it )
  {
    iterator it2 = it;
    ++it2;
    if ( it2 == contour.end() )
      it2 = contour.begin();
    return it2;
  }

  /**
   * Given a circulator [it] on the contour, returns the one on the
   * previous element.
   * @param it an iterator in the list of nodes [contour].
   * @return the previous one (assuming the contour is a loop).
   */
  iterator prev( const iterator & it )
  {
    iterator it2 = it;
    if ( it2 == contour.begin() )
      return --contour.end();
    --it2;
    return it2;
  }

  /**
   * Given a circulator [it] on the contour, returns the one on the
   * next element.
   * @param it an iterator in the list of nodes [contour].
   * @return the next one (assuming the contour is a loop).
   */
  const_iterator next( const const_iterator & it ) const
  {
    const_iterator it2 = it;
    ++it2;
    if ( it2 == contour.end() )
      it2 = contour.begin();
    return it2;
  }

  /**
   * Given a circulator [it] on the contour, returns the one on the
   * previous element.
   * @param it an iterator in the list of nodes [contour].
   * @return the previous one (assuming the contour is a loop).
   */
  const_iterator prev( const const_iterator & it ) const
  {
    const_iterator it2 = it;
    if ( it2 == contour.begin() )
      return contour.end();
    --it2;
    return it2;
  }

  /**
   * Inserts the node [n] before the element pointed by [it]. The user
   * should take care beforehands that the starting position is not
   * modified.
   *
   * @param it an iterator in the contour.
   * @param n a node to insert (which is copied).
   * @return an iterator pointing on the created node within the contour.
   */
  iterator insertNode( const iterator & it, const Node & n )
  {
    return contour.insert( it, n );
  }

  /**
   * The node pointed by [it] is assigned the node [n].  The user
   * should take care beforehands that the starting position is not
   * modified.
   *
   * @param it an iterator in the contour.
   * @param n the node to copy.
   */
  void modifyNode( const iterator & it, const Node & n )
  {
    *it = n;
  }

  /**
   * Remove the node pointed by [it] from the dynamic contour.  The user
   * should take care beforehands that the starting position is not
   * modified.
   *
   * @param it an iterator in the contour.
   * @return an iterator pointing on the node which follows the suppressed one.
   */
  iterator removeNode( const iterator & it )
  {
    return contour.erase( it );
  }

  /**
   * Moves the start node of the contour one node forward. Useful
   * before modifying a node which happens to be the start node of the
   * contour.
   */
  void moveStartNodeForward()
  {
    start_pos += start_node->n * start_node->v;
    start_node = next( start_node );
  }

  /**
   * Moves the start node of the contour one node backward. Useful
   * before modifying a node which happens to be the start node of the
   * contour.
   */
  void moveStartNodeBackward()
  {
    start_node = prev( start_node );
    start_pos -= start_node->n * start_node->v;
  }

  /**
   * Moves the start node of the contour one node backward only if it
   * is equal to [it]. Useful before modifying a node which happens to
   * be the start node of the contour.
   *
   * @param it an iterator in the dynamic contour.
   */
  void moveStartNodeBackwardIfEqual( const iterator & it )
  {
    if ( start_node == it )
      {
	start_node = prev( start_node );
	start_pos -= start_node->n * start_node->v;
      }
  }

  /**
   * Moves the start Node to the beginning of the contour,
   * i.e. 'contour.begin()'. Useful so that the the logical start
   * reflects the physical start.
   */
  void moveStartNodeToBegin()
  {
    while ( start_node != contour.begin() )
      moveStartNodeBackward();
  }

  /**
   * The start node should not be either on [it] or on 'prev(it)'.
   * Models an elementary swap at the base vertex of [it]. It means
   * that the vectors of 'prev(it)' of [it] are exchanged. If [flip]
   * is true, then the vertex touch type is flipped INSIDE/OUTSIDE,
   * otherwise it it left unchanged.
   *
   * @param it an iterator on a node, whose base vertex is the one
   * being swapped.
   *
   * @param flip when 'true', the touch type of [it] is reversed (no
   * modification of the underlying digital shape), otherwise it is
   * left unchanged (the digital shape has one more or one less
   * pixel).
   */
  void swapNode( const iterator & it, bool flip = true )
  {
    iterator itp = prev( it );
    std::swap( itp->v, it->v );
    std::swap( itp->oc, it->oc );
    std::swap( itp->z, it->z );
    if ( flip ) it->flip();
  }

  /**
   * The node pointed by [it] has a direction vector v and a base type
   * t (an INSIDE pixel or an OUTSIDE pixel).
   *
   * @param it an iterator on a node of the contour.
   *
   * @return RIGHT if it is to the right of v that the pixels are of
   * type t, otherwise left. Depends of course if the contour is
   * considered clockwise or counterclockwise.
   *
   * @see cw
   */
  SideType getSameSide( const iterator & it ) const
  {
    return ( ( it->base == INSIDE && cw ) ||  ( it->base == OUTSIDE && !cw ) )
      ? RIGHT : LEFT;
  }

  /**
   * The node pointed by [it] has a direction vector v and a base type
   * t (an INSIDE pixel or an OUTSIDE pixel).
   *
   * @param it an iterator on a node of the contour.
   *
   * @return RIGHT if it is to the right of v that the pixels are \b
   * not of type t, otherwise left. Depends of course if the contour
   * is considered clockwise or counterclockwise.
   *
   * @see cw
   */
  SideType getOtherSide( const iterator & it ) const
  {
    return ( ( it->base == INSIDE && cw ) ||  ( it->base == OUTSIDE && !cw ) )
      ? LEFT : RIGHT;
  }

  /**
   * Given an iterator on a node, returns what is the shape of the
   * vertex (base of the node) according to the directions of its
   * surrounding nodes, but also to the fact it is an INSIDE or
   * OUTSIDE vertex, and according to the fact that the contour is
   * clockwise or counterclockwise.
   *
   * @param it any iterator on a node of the contour, whose base
   * vertex is the vertex of interest.
   *
   * @return the connection type between 'prev(it)' and [it], in WEDGE
   * (i.e. convex-like), FLAT, AFLAT (i.e. almost flat), VEE
   * (concave).
   */
  ConnectionType getConnectionType( const iterator & it )
  {
    // Assumes local righthandside rule.
    int d = prev( it )->v.det( it->v );
    // inverts it according to touch type of [it] and contour
    // orientation.
    if ( getSameSide( it ) == LEFT )
      d = -d;
    if ( d < 0 ) return WEDGE;
    else if ( d == 0 ) return FLAT;
    else if ( d == 1 ) return AFLAT;
    else return VEE;
  }


  /**
   * The node should be primitive non trivial.  Splits the node E(zk)
   * into E(zk-1)E(zk-2) according to the parity of the slope. Returns
   * an iterator on the second part of the fraction.
   */
  iterator splitBerstelLowerChristoffel( iterator & it, TouchType split_tt )
  {
    if ( it->n != 1 ) 
      cerr << "[DynamicContour::splitBerstelLowerChristoffel] Node should be primitive."
	   << endl;
    iterator itn = next( it );
    //    TouchType split_tt = getTouchTypeForSplit( true );
    Node nr;
    Node & nl = *it;
    CFraction z( it->z );
    if ( ! z.splitBerstel( nl.z, nl.n, nr.z, nr.n ) ) 
      {
	// cerr << "[splitLC] " << z << " no split" << endl;
	return it;
      }
    nr.base = split_tt;
    nr.oc = nl.oc;
    nl.v = nl.oc.castBack( Vector2i( nl.z.q(), nl.z.p() ) );
    nr.v = nr.oc.castBack( Vector2i( nr.z.q(), nr.z.p() ) );
    // cerr << "[splitLC] " << z << " = (" << nl.z << ")^" << nl.n
    // 	 << " + " << " (" << nr.z << ")^" << nr.n << endl;
    if ( ( split_tt != it->base )
	 && ( it->n > 1 ) )
      {
	// Should subdivide first node when the aligned inner vertices
	// have not the same touch type as the base of the vector.
	iterator it_inter = subdivide( it, it->n - 1 );
	it_inter->base = split_tt;
      }
    return insertNode( itn, nr );
  }

  /**
   * The node should be primitive non trivial.  Splits the node E(zk)
   * into E(zk-1)E(zk-2) according to the parity of the slope. Returns
   * an iterator on the second part of the fraction.
   */
  iterator splitBerstelUpperChristoffel( iterator & it, TouchType split_tt )
  {
    if ( it->n != 1 ) 
      cerr << "[DynamicContour::splitBerstelUpperChristoffel] Node should be primitive."
	   << endl;
    iterator itn = next( it );
    //    TouchType split_tt = getTouchTypeForSplit( false );
    Node nr;
    Node & nl = *it;
    CFraction z( it->z );
    if ( ! z.splitBerstel( nr.z, nr.n, nl.z, nl.n ) )
      {
	// cerr << "[splitUC] " << z << " no split" << endl;
	return it;
      }
    nr.base = split_tt;
    nr.oc = nl.oc;
    nl.v = nl.oc.castBack( Vector2i( nl.z.q(), nl.z.p() ) );
    nr.v = nr.oc.castBack( Vector2i( nr.z.q(), nr.z.p() ) );
    // cerr << "[splitUC] " << z << " = (" << nl.z << ")^" << nl.n
    // 	 << " + " << " (" << nr.z << ")^" << nr.n << endl;
    if ( ( split_tt != it->base )
	 && ( it->n > 1 ) )
      {
	// Should subdivide first node when the aligned inner vertices
	// have not the same touch type as the base of the vector.
	iterator it_inter = subdivide( it, it->n - 1 );
	it_inter->base = split_tt;
	// subdivide( it, it->n - 1 );
      }
    return insertNode( itn, nr );
  }

  /**
   * Given some iterator [it] on a node, split it at its target
   * (front) so that there is a trivial node at its front (in fact,
   * depth <= 1). The node vector is split on its lower side (standard
   * Berstel/splitting decomposition) when [lower] is 'true',
   * otherwise it is split on its upper side (the reverse). [t] gives
   * the base type assigned to each created base vertex.
   *
   * @param itn an iterator on a node of the contour.
   * @param lower the type of split.
   * @param t the base type for the created base vertices.
   * @return an iterator on the created trivial node.
   */
  iterator splitAtTargetTilQuasiTrivial( iterator itn, bool lower,
					 TouchType t )
  {
    // iterator itn = prev( it );
    while ( true )
      {
	if ( itn->n > 1 ) 
	  itn = subdivide( itn, 1 );
	iterator itr = 
	  ( itn->oc.isIndirect() ? ! lower : lower )
	  ? splitBerstelLowerChristoffel( itn, t )
	  : splitBerstelUpperChristoffel( itn, t );
	if ( itr == itn ) break;
	itn = itr;
      }
    return itn;
  }

  /**
   * Given some iterator [it] on a node, split it at its source/base (back)
   * so that [it] points to a trivial node (in fact, depth <= 1). The
   * node vector is split on its lower side (standard
   * Berstel/splitting decomposition) when [lower] is 'true',
   * otherwise it is split on its upper side (the reverse). [t] gives
   * the base type assigned to each created base vertex.
   *
   * @param itn an iterator on a node of the contour.
   * @param lower the type of split.
   * @param t the base type for the created base vertices.
   * @return an iterator on the created trivial node (in fact [it]).
   */
  iterator splitAtBaseTilQuasiTrivial( iterator itn, bool lower,
				       TouchType t )
  {
    // iterator itn = it;
    while ( true )
      {
	if ( itn->n > 1 ) 
	  subdivide( itn, itn->n - 1 );
	iterator itr = 
	  ( itn->oc.isIndirect() ? ! lower : lower )
	  ? splitBerstelLowerChristoffel( itn, t )
	  : splitBerstelUpperChristoffel( itn, t );
	if ( itr == itn ) break;
      }
    return itn;
  }

  /**
   * Given some iterator [it] on a node, split it at its target
   * (front) so that there is a trivial node at its front (in fact,
   * vector (1,0) repeated once). The node vector is split on its
   * lower side (standard Berstel/splitting decomposition) when
   * [lower] is 'true', otherwise it is split on its upper side (the
   * reverse). [t] gives the base type assigned to each created base
   * vertex.
   *
   * @param itn an iterator on a node of the contour.
   * @param lower the type of split.
   * @param t the base type for the created base vertices.
   * @return an iterator on the created trivial node.
   */
  iterator splitAtTargetTilTrivial( iterator itn, bool lower,
				    TouchType t )
  {
    iterator itt = splitAtTargetTilQuasiTrivial( itn, lower, t );
    iterator itt2 = splitQuasiTrivial( itt, lower, t );
    return subdivide( itt2, 1 );
  }

  /**
   * Given some iterator [it] on a node, split it at its source/base
   * (back) so that [it] points to a trivial node (in fact, vector
   * (1,0) repeated once). The node vector is split on its lower side
   * (standard Berstel/splitting decomposition) when [lower] is
   * 'true', otherwise it is split on its upper side (the
   * reverse). [t] gives the base type assigned to each created base
   * vertex.
   *
   * @param itn an iterator on a node of the contour.
   * @param lower the type of split.
   * @param t the base type for the created base vertices.
   * @return an iterator on the created trivial node (in fact [it]).
   */
  iterator splitAtBaseTilTrivial( iterator itn, bool lower,
				  TouchType t )
  {
    iterator itb = splitAtBaseTilQuasiTrivial( itn, lower, t );
    splitQuasiTrivial( itb, lower, t );
    subdivide( itb, itb->n - 1 );
    return itb;
  }
 

  /**
   * Checks if vertex is neutral, meaning that the two edges around
   * have opposite vectors and base types. 
   *
   * It may also remove an inconsistent vertex: prev( it ) and it are
   * opposite, of same base type, but next( it ) has a base type
   * different from prev( it ), although it is in the same place. For
   * now, only inconsistent vertex with trivial surrounding nodes are
   * removed. They may appear in a VEE splitting with two nodes in
   * different quadrants.
   *
   * @param it (modified) an iterator in this contour, which is moved
   * to the next valid nodes.
   *
   * @return 'true' if the contour was changed, 'false' otherwise.
   */
  bool removeNeutral( iterator & it )
  {
    iterator itp = prev( it );
    if ( ( it->v.x() == -itp->v.x() )
	 && ( it->v.y() == -itp->v.y() )
	 )
      if ( it->base != itp->base )
	{	
	  moveStartNodeBackwardIfEqual( it );
	  moveStartNodeBackwardIfEqual( itp );
	  iterator itn = next( it );
	  if ( itp->n == it->n )
	    {
	      removeNode( itp );
	      removeNode( it );
	    }
	  else if ( itp->n > it->n )
	    {
	      itp->n -= it->n;
	      removeNode( it );
	    }
	  else
	    {
	      // it->base = itp->base;
	      it->n -= itp->n;
	      removeNode( itp );
	    }
	  it = itn;
	  return true;
	}
      else 
	{
	  iterator itn = next( it );
	  if ( ( itn->base != itp->base ) 
	       && ( it->n == 1 ) && ( itp->n == 1 ) )
	    {
	      // cerr << "[removeNeutral]: remove inconsistency." << endl;
	      moveStartNodeBackwardIfEqual( it );
	      moveStartNodeBackwardIfEqual( itp );
	      TouchType t = it->base;
	      subdivide( itn, itn->n - 1 );
	      removeNode( itp );
	      removeNode( it );
	      itn->base = t;
	      it = itn;
	      return true;
	    }
	}
    return false;
  }

  /**
   * Checks if vertex is back and forth, meaning that the two edges
   * around have opposite vectors. The new base type depends on the
   * surrounding base types. If prev( it )->base == next( it )-> base
   * then the new base type is this one, otherwise, the new base type
   * is the same as it->base.
   *       
   @verbatim
          !      would give    !
          v                    v
   it x<==x                itn x==>o==>x
       ==>(o==>)^2x             
   @endverbatim
   *
   * Contrary to removeNeutral, this method will always remove the two
   * nodes surrounding the vertex (if they are opposite) and assign a
   * some base type for next(it). It this base type is changed, next(
   * it ) is subdivided. 
   * 
   * Note also that the user should have moved away the start node
   * from this zone.
   *
   * @param it (modified) an iterator in this contour, which is moved
   * to the next valid node.
   *
   * @return 'true' if the contour was changed, 'false' otherwise.
   */
  bool removeOpposite( iterator & it )
  {
    iterator itp = prev( it );
    if ( ( it->v.x() == -itp->v.x() )
	 && ( it->v.y() == -itp->v.y() ) )
      {	
	iterator itn = next( it );
	TouchType t = ( itn->base == itp->base )
	  ? itp->base : it->base;
	removeNode( itp );
	removeNode( it );
	it = itn;
	subdivide( it, it->n-1 );
	it->base = t;
	return true;
      }
    return false;
  }

  /**
   * This method decomposes a quasi trivial vector ( (1,x) or (x,1) )
   * into trivial vectors (1,0), (0,1) in different quadrant.
   *
   * [it] points to a primitive node with some vector (1,x) or (x,1)
   * or (-1,x) or (x,-1).  Decomposes it into two nodes in different
   * quadrant.
   *
   * The node vector is split on its lower side (standard
   * Berstel/splitting decomposition) when [lower] is 'true',
   * otherwise it is split on its upper side (the reverse). [t] gives
   * the base type assigned to each created base vertex.
   *
   * Note also that the user should have moved away the start node
   * from this zone.
   *
   * @param it an iterator in this contour which points to a quasi trivial node.
   * @param lower indicates if the split must be lower or upper.
   * @param t the base type for the created base vertices.
   * @return an iterator on the created trivial node (which follows [it]).
   */
  iterator
  splitQuasiTrivial( const iterator & it, bool lower, TouchType t ) 
  {
    if ( it->n > 1 )
      cerr << "[DynamicContour::splitQuasiTrivial] node should be primitive."
	   << endl;
    if ( it->z.depth() > 1 )
      cerr << "[DynamicContour::splitQuasiTrivial] vector is too complex."
	   << it->v << " z=" << it->z << endl;
    if ( it->z.p() == 0 ) return it;
    if ( it->z.p() != 1 ) 
      cerr << "[DynamicContour::splitQuasiTrivial] vector should be some 1/k."
	   << it->v << " z=" << it->z << endl;
    // TouchType t = getTouchTypeForSplit( lower );
    lower = ( it->oc.isIndirect() ? ! lower : lower );
    // cerr << "[DynamicContour::splitQuasiTrivial] " 
    // 	 << it->oc << " "
    // 	 << ( lower ? "LOWER" : "UPPER" )
    // 	 << " v=" << it->v << " z=" << it->z 
    // 	 << " t_base=" << ( it->base == INSIDE ? "I" : "O" )
    // 	 << " t=" << ( t == INSIDE ? "I" : "O" )
    // 	 << endl;

    if ( lower )
      {
	uint repetition = it->z.q();
	Node nv( it->oc.castBack( Vector2i( 0, 1 ) ), t, 1 );
	Node nh( it->oc.castBack( Vector2i( 1, 0 ) ), it->base, repetition );
	// cerr << "(nh " << nh.v << " nv " << nv.v << ")";
	modifyNode( it, nh );
	iterator itn = insertNode( next( it ), nv ); 
	// cerr << "[DynamicContour::splitQuasiTrivial] " << it->base << " " << t << repetition << endl;
	if ( ( it->base != t ) && ( repetition > 1 ) )
	  {
	    // cerr << "[DynamicContour::splitQuasiTrivial] subdivision for correct touch type." << endl;
	    // takes care that new nodes have the correct TouchType. 
	    iterator itn2 = subdivide( it, repetition - 1 );
	    itn2->base = t;
	    return itn2;
	  }
	return itn;
      }
    else
      {
	Node nv( it->oc.castBack( Vector2i( 0, 1 ) ), it->base, 1 );
	Node nh( it->oc.castBack( Vector2i( 1, 0 ) ), t, it->z.q() );
	// cerr << "(nv " << nv.v << " nh " << nh.v << ")";
	modifyNode( it, nv );
	return insertNode( next( it ), nh ); 
      }
  }

  /**
   * This is the main simplification method when computing dynamically
   * a minimum length polygon. Its objective is to have only WEDGE or
   * FLAT vertices. If some vertex (pointed by [it]) is some VEE, then
   * it changes the vertex (and the surrounding nodes) to make the
   * dynamic contour pass through the symetric vertex (diagonally
   * opposite to it in the one-wide band of the digital contour),
   * which is of course of opposite base type. After that, this vertex
   * is a WEDGE (since the base type has changed). Furthermore, this
   * method will remove neutral vertices (FLAT, but opposite and with
   * opposite base type) and AFLAT vertices.
   *
   * This is not generally a O(1) operation. Split takes O(log(n))
   * (and for now O(log^2(n))).
   *
   * @param it (modifed) an iterator on a node of the contour
   * (specifies the base vertex). The iterator may be moved to point
   * on a valid vertex.
   *
   * @return 'true' if the contour was changed, 'false' otherwise.
   *
   * @see ConnectionType.
   */
  bool fusion( iterator & it );


  /**
   * Given a node pointed by [it] with n repetitions, decompose it
   * into two nodes, one with n-[k] repetitions to the left, one with
   * [k] repetitions to the right. Return the iterator in the middle
   * if the node has been split, otherwise [it] itself.
   *
   * O(1) operation.
   *
   * Note also that the user should have moved away the start node
   * from this zone.
   *
   * @param it any iterator on the contour (points on the node that is
   * subdivided).
   * 
   * @param k indicates how it is subdivided into two aligned nodes.
   *
   * @return an iterator on the created node (after [it]), or [it]
   * itself if nothing was changed.
   */
  iterator subdivide( const iterator & it, uint k )
  {
    iterator it2 = it;
    if ( ( k != 0 ) && ( k < it->n ) )
      {
	it2 = next( it );
	it2 = insertNode( it2, *it );
	it->n -= k;
	it2->n = k;
      }
    return it2;
  }

  /**
   * This method merges an AFLAT vertex designated by the two
   * consecutive nodes [it] and [itn].
   *
   * Merge the two nodes pointed by [it] and [itn], provided they are
   * mergeable. Not very efficient for now, since we do not exploit
   * the repetitions of the vectors. Assume that the vectors of the
   * two nodes have determinant 1 in a clockwise contour with an
   * INSIDE join.
   *
   * Note also that the user should have moved away the start node
   * from this zone.
   *
   * @param it an iterator in the contour (the node is subdivided if
   * its repetition number is greater than 1).
   *
   * @param itn an iterator in the contour equal to 'next( it )' (the
   * node is subdivided if its repetition number is greater than 1).
   *
   * @return an iterator on the merged node.
   *
   * @TODO do it fast !
   */
  iterator lowerFusion( const iterator & it, const iterator & itn )
  {
    // we already know that these two consecutive nodes (with an
    // INSIDE between) are mergeable (det=1).
    iterator itn_left = subdivide( itn, itn->n - 1 );
    iterator it_right = subdivide( it, 1 );
    Node n( it_right->v + itn_left->v, 
	    it_right->base, 
	    1 );
    removeNode( itn );
    modifyNode( it_right, n );
    return it_right;
  }

  /**
   * This method merges an AFLAT vertex designated by the two
   * consecutive nodes [it] and [itn].
   *
   * Merge the two nodes pointed by [it] and [itn], provided they are
   * mergeable. Not very efficient for now, since we do not exploit
   * the repetitions of the vectors. Assume that the vectors of the
   * two nodes have determinant 1 in a clockwise contour with an
   * INSIDE join. For now, calls lowerFusion.
   *
   * Note also that the user should have moved away the start node
   * from this zone.
   *
   * @param it an iterator in the contour (the node is subdivided if
   * its repetition number is greater than 1).
   *
   * @param itn an iterator in the contour equal to 'next( it )' (the
   * node is subdivided if its repetition number is greater than 1).
   *
   * @return an iterator on the merged node.
   *
   * @TODO do it fast !
   */
  iterator upperFusion( const iterator & it, const iterator & itn )
  {
    return lowerFusion( it, itn );
  }

  /**
   * This method initializes consistently a dynamic contour from a
   * Freeman chain. The chain is supposed to model either the inner
   * 4-connected 8-border of an object ([inside] is 'true') or its
   * outer 4-connected 8-border ([inside] is 'false'). If you wish to
   * initialize the dynamic contour with a digital interpixel contour,
   * you should use the method FreemanChain::innerContour beforehands.
   *
   * The iterators should have a getCode() method returning between 0-3
   * according to the Freeman chain code and a nextInLoop() method to
   * visit the whole Freeman chain code.
   *
   * @param itb an iterator on the beginning of Freeman chain.
   * @param ite an iterator on the end of Freeman chain.
   * @param inside indicates if it is an inner or an outer border.
   */
  template <class FreemanChainIterator>
  void initFromFreemanChain( FreemanChainIterator itb, 
			     FreemanChainIterator ite, 
			     bool inside )
  {
    contour.clear();
    start_pos = *itb;
    list<Node>::iterator itnode;
    Vector2i v;
    TouchType base = DynamicContour::getTouchType( inside );
    do
      {
	itnode = contour.end();
	switch ( itb.getCode() ) 
	  {
	  case 0: v = Vector2i( 1, 0 ); break;
	  case 1: v = Vector2i( 0, 1 ); break;
	  case 2: v = Vector2i( -1, 0 ); break;
	  case 3: v = Vector2i( 0, -1 ); break;
	  }
	insertNode( itnode, Node( v, base, 1 ) );
	itb.nextInLoop();
      }
    while ( itb != ite );
    start_node = contour.begin();
  }

  uint simplifyAt( iterator & it )
  {
    uint nb = 0;
    bool mr = true;
    bool ml = true;
    while ( mr && ! null() ) 
      {
	mr = fusion( it );
	nb += mr ? 1 : 0;
	// mr = mergeRight( it );
	// nb += mr ? 1 : 0;
	// ml = mergeLeft( it );
	// nb += ml ? 1 : 0;
      }
    // it->selfDisplay( cerr );
    return nb;
  }

  uint simplify()
  {
    // cerr << "[DynamicContour::simplify]";
    uint nb = 0;
    uint nb0;
    do 
      {
	nb0 = 0;
	iterator itb = start_node; // contour.begin();
	if ( ! null() )
	  nb0 += simplifyAt( itb );
	iterator it = next( start_node );
	while ( ( ! null() ) && ( it != start_node ) )
	  {
	    nb0 += simplifyAt( it );
	    // cerr << "+" << nb0 << flush;
	    // selfDisplay( cerr );
	    it = next( it );
	  }
	nb += nb0;
      }
    while ( nb0 != 0 );

    return nb;
  }


  void getPoints( vector<Vector2i> & pts,
		  vector<bool> & inside ) const
  {
    if ( ! null() ) 
      {
	Vector2i start_pt( start_pos );
	const_iterator it = start_node;
	do
	  { 
	    pts.push_back( start_pt );
	    inside.push_back( it->base );
	    start_pt += it->n * it->v;
	    it = next( it );
	  }
	while ( it != start_node );
      }
  } 
		  

  static TouchType getTouchType( bool inside )
  {
    return inside ? INSIDE : OUTSIDE;
  }

};

ostream &
operator<<( ostream & out, const DynamicContour::Node & n )
{
  n.selfDisplay( out );
  return out;
}

ostream &
operator<<( ostream & out, const DynamicContour & n )
{
  n.selfDisplay( out );
  return out;
}


/**
 * Given a node [it], tries to merge it with the previous node. The
 * fusion is related to the connection type of the two nodes.
 */
bool
DynamicContour::fusion( iterator & it )
{
  iterator itp = prev( it ); 
  ConnectionType vtype = getConnectionType( it );
  bool change = false;
  // cerr << "[fusion] before:";
  // selfDisplayWithIterator( cerr, it );

  switch ( vtype ) 
    {
    case FLAT:
      // vectors are aligned. Merge is O(1).
      // cerr << "[fusion] FLAT " << (*itp) << " - " << (*it) << endl; 
      if ( ( itp->v.x() == it->v.x() ) 
	   && ( itp->v.y() == it->v.y() )
	   && ( itp->base == it->base ) )
	{
	  // TODO ! changes things
	  moveStartNodeBackwardIfEqual( it );
	  moveStartNodeBackwardIfEqual( itp );
	  it->n += itp->n;
	  //itp->n = 0;
	  removeNode( itp ); //contour.erase( itp ); //removeNode( itp );
	  change = true;
	}
      else 
	{ // vectors are opposite. Check if they can be removed.
	  // cerr << "[fusion] NON STD FLAT : ISNEUTRAL " << (*itp) << " - " << (*it) << endl; 
	  if ( removeNeutral( it ) )
	    change = true;
	}
      break;
    case AFLAT:
      // cerr << "[fusion] AFLAT " << (*itp) << " - " << (*it) << endl; 
      moveStartNodeBackwardIfEqual( it );
      moveStartNodeBackwardIfEqual( itp );
      if ( it->base == INSIDE )
	it = lowerFusion( itp, it );
      else
	it = upperFusion( itp, it );
      change = true;
      break;
    case WEDGE:
      // nothing to do.
      // cerr << "[fusion] WEDGE " << (*itp) << " - " << (*it) << endl; 
      break;
    case VEE:
      // In this case, we have to flip INSIDE/OUTSIDE at this position.
      // cerr << "[fusion] VEE " << (*itp) << " - " << (*it) << endl; 
      bool lower = getOtherSide( it ) == RIGHT;
      // The created vertices are on the other side of the shape, so
      // they should lie on OUTSIDE pixels if [it] was lying on an
      // INSIDE pixel, and inversely.
      TouchType t = it->base == INSIDE ? OUTSIDE : INSIDE;
      moveStartNodeBackwardIfEqual( it );
      moveStartNodeBackwardIfEqual( itp );

      iterator itl = splitAtTargetTilTrivial( prev( it ), lower, t ); 
      iterator itr = splitAtBaseTilTrivial( it, lower, t );
      // in case the nodes where not in the same quadrant.
      if ( removeOpposite( itr ) )
	it = itr;
      else 
	{ // standard case.
	  swapNode( itr, true );
	  it = itl;
	}
      change = true;
      break;

    }
  // if ( change )
  //   {
  //     // cerr << "[fusion] post:  ";
  //     // selfDisplayWithIterator( cerr, it );
  //   }
  // else
  //   cerr << "[fusion] post: no change." << endl;
  return change;
}




static Arguments args;



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
  // StandardArguments::addDigitalArgs( args, 2, false, false );
  StandardArguments::addIOArgs( args, true, false );
  args.addOption( "-orientation", "-orientation <CW,CCW,COMP>: gives the orientation of the input Freeman chain. COMP means that the orientation is given by the number of cw or ccw loops made by the contour. Default is COMP.", "COMP" );
  args.addOption( "-contour_type", "-contour_type <INNER,OUTER,INTERPIXEL>: indicates whether the input freeman chain codes the inner border of a digital object made of inside pixels (INNER), its outer border made of outside pixels (OUTER) or its interpixel digital contour (INTERPIXEL). Default is INNER.", "INNER" );
  args.addBooleanOption( "-mlp", "-mlp: simplifies the input contour to get the dynamic contour that is its minimum length polygon." );
  args.addBooleanOption( "-dDC", "-dDC: displays the dynamic contour as its vector representation." );
  args.addBooleanOption( "-dPTS", "-dPTS: displays the dynamic contour as a list of points, with 0/1 for inside/outside." );
  args.addBooleanOption( "-dPTShalf", "-dPTShalf: displays the dynamic contour as a list of points in the half-integer plane, with 0/1 for inside/outside." );
  args.addBooleanOption( "-dPointels", "-dPointels: displays the dynamic contour as its list of pointels. Should give back the input *interpixel* contour." );
  args.addBooleanOption( "-dRPointels", "-dRPointels: displays the dynamic contour as its list of pointels, read in reverse order. " );
  args.addBooleanOption( "-msplit", "-msplit: split every edge at middle position." );
  args.addBooleanOption( "-close", "-close: close output contours by repeating the first point after the last." );
  args.addOption( "-splitAt", "-splitAt <pos>: split at position <pos>.", "0" );
  args.addOption( "-flipAt", "-flipAt <pos>: flip at position <pos>.", "0" );
  StandardArguments::addDebugArgs( args, true, true, true );

  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "test_DMLP", 
			  "Reads a freeman contour. Tests dynamic minimum length polygon.",
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
  istream & in_str = StandardArguments::openInput( args );
  FreemanChain::read( in_str, c );
  if ( ! in_str.good() )
    {
      if ( debug )
	cerr << "Error reading Freeman chain code." << endl;
      return 2;
    }

  if ( debug )
    {
      bool ok = UndoableList<int>::selfTest( cerr );
      if ( ok )
	cerr << "---OK--- UndoableList<int>::selfTest()" << endl;
      else
	cerr << "--ERROR- UndoableList<int>::selfTest()" << endl;
    }

  // DynamicContour dc;
  string str_orientation = args.getOption( "-orientation" )->getValue( 0 );
  bool cw = true;
  int nb_ccw_loops = c.isClosed();
  cerr << "Number of CCW loops = " << nb_ccw_loops << endl;
  if ( str_orientation == "CW" )       cw = true;
  else if ( str_orientation == "CCW" ) cw = false;
  else                                 cw = nb_ccw_loops < 0;
  // DMLPContour dc( args.getOption( "-orientation" )->getValue( 0 ) == "CW" );
  DMLPContour dc( cw );
  // OrderedAlphabet A( '0', 4 );
  // FreemanChain::const_iterator it = c.findQuadrantChange4( A );
  // Vector2i start_pt( *it );
  FreemanChain::const_iterator it = c.begin();
  string contour_type = args.getOption( "-contour_type" )->getValue( 0 );
  if ( contour_type == "INNER" )
    dc.initFromFreemanChain( it, it, true );
  else if ( contour_type == "OUTER" )
    dc.initFromFreemanChain( it, it, false );
  else if ( contour_type == "INTERPIXEL" )
    {
      // FreemanChain inner_chain;
      // vector<uint> outer2inner;
      // vector<uint> inner2outer;
      // FreemanChain::innerContour( inner_chain, 
      // 				  outer2inner, inner2outer,
      // 				  c, ! dc.cw() );
      // it = inner_chain.begin();
      // dc.initFromFreemanChain( it, it, true );
      dc.initFromInterpixelContour( c );
    }
  if ( trace > 0 )
    {
      cerr << "[INPUT] size=" << dc.size() << endl;
      if ( trace > 1 )
	{
	  dc.selfDisplay( cerr );
	}
    }

  if ( args.check( "-mlp" ) )
    {
      cerr << setprecision( 10 );
      cerr << "<before mlp> length=" << dc.getLength() 
	   << " nbup=" << dc.nbOfUpdates() << endl;
      uint nb = dc.simplify();
      if ( trace > 0 ) 
	cerr << "[DynamicContour::simplify] " << nb << " changes." << endl;
      // nb += dc.simplify();
      // if ( trace > 0 ) 
      // 	cerr << "[DynamicContour::simplify] " << nb << " changes." << endl;
      dc.moveStartEdgeToBegin();
      cerr << "<after mlp> length=" << dc.getLength()
	   << " nbup=" << dc.nbOfUpdates() << endl;
      dc.updateLength();
      cerr << "<after mlp> updated_length=" << dc.getLength() 
	   << " nbup=" << dc.nbOfUpdates() << endl;
    }

  if ( trace > 0 )
    {
      cerr << "[OUTPUT] size=" << dc.size() << endl;
      if ( trace > 1 )
	{
	  dc.selfDisplay( cerr );
	}
    }

  if ( args.check( "-dDC" ) )
    {
      dc.selfDisplay( cout );
    }


  if ( args.check( "-msplit" ) )
    {
      cerr << "[Splitting at middle of each edge]" << endl;
      for ( DMLPContour::edge_iterator it = dc.beginEdge();
	    it != dc.endEdge();
	    )
	{
	  cerr << "  * Edge: " << it->v << endl;
	  // memorizes next iterator to split the "true" next edge.
	  DMLPContour::edge_iterator next_it = it;
	  ++next_it;
	  dc.moveStartEdgeBackwardIfEqual( it );
	  // Measure length of edge.
	  uint l = it->digitalLength();
	  bool lower = dc.getOtherSide( it ) == DMLPContour::RIGHT;
	  // The created vertices are on the other side of the shape, so
	  // they should lie on OUTSIDE pixels if [it] was lying on an
	  // INSIDE pixel, and inversely.
	  DMLPContour::TouchType t = 
	    it->base == DMLPContour::INSIDE 
	    ? DMLPContour::OUTSIDE 
	    : DMLPContour::INSIDE;
 	  dc.splitEdgeAt( it, l/2, lower, t );
	  it = next_it;
	}
    }

  if ( args.check( "-splitAt" ) )
    {
      uint pos = args.getOption( "-splitAt" )->getIntValue( 0 );
      cerr << "[Splitting at position " << pos << "]" << endl;
      DMLPContour::iterator itbegin ( dc.first( true ) );
      // itbegin.previous();
      for ( ; pos != 0; --pos ) ++itbegin;
      cerr << "  * Edge: " << itbegin.edge()->v << endl;
      Vector2i pt = itbegin.pointel();
      cerr << "  * Pointel " << pt.x() << "," << pt.y() << endl;
      // memorizes next iterator to split the "true" next edge.
      dc.moveStartEdgeBackwardIfEqual( itbegin.edge() );
      bool lower = dc.getOtherSide( itbegin.edge() ) == DMLPContour::RIGHT;
      // The created vertices are on the other side of the shape, so
      // they should lie on OUTSIDE pixels if [it] was lying on an
      // INSIDE pixel, and inversely.
      DMLPContour::TouchType t = 
	itbegin.edge()->base == DMLPContour::INSIDE 
	? DMLPContour::OUTSIDE 
	: DMLPContour::INSIDE;
      DMLPContour::iterator itnew = dc.splitAround( itbegin, lower, t );
      cerr << "  * Pointel " << itnew.pointel() << endl;
    }

  if ( args.check( "-flipAt" ) )
    {
      uint pos = args.getOption( "-flipAt" )->getIntValue( 0 );
      cerr << "[Flipping at position " << pos << "]" << endl;
      DMLPContour::iterator itbegin ( dc.first( true ) );
      // itbegin.previous();
      for ( ; pos != 0; --pos ) ++itbegin;
      cerr << "  * Edge: " << itbegin.edge()->v << endl;
      Vector2i pt = itbegin.pointel();
      cerr << "  * Pointel " << pt.x() << "," << pt.y() << endl;
      // memorizes next iterator to split the "true" next edge.
      bool ok = dc.flip( itbegin, true );
      cerr << "  * [" << ( ok ? "OK" : "ERR" ) << "]" << endl;
    }

  if ( args.check( "-dPTS" ) )
    {
      vector<Vector2i> pts;
      vector<bool> ins;
      dc.getPoints( pts, ins );
      cout << "# -dPTS" << endl
	   << "# " << pts.size() << " points." << endl;
      uint l = pts.size();
      for ( uint i = 0; i < l; ++i )
	cout << pts[ i ].x() << " " << pts[ i ].y() 
	     << " " << ins[ i ] 
	     << " " << pts[ ( i + 1 ) % pts.size() ].x()
	     << " " << pts[ ( i + 1 ) % pts.size() ].y()
	     << endl;
      if ( args.check( "-close" ) )
	cout << pts[ 0 ].x() << " " << pts[ 0 ].y() 
	     << " " << ins[ 0 ] 
	     << " " << pts[ ( 0 + 1 ) % pts.size() ].x()
	     << " " << pts[ ( 0 + 1 ) % pts.size() ].y()
	     << endl;
 
    }
  if ( args.check( "-dPTShalf" ) )
    {
      vector<Vector2i> pts;
      vector<bool> ins;
      dc.getPoints( pts, ins );
      cout << "# -dPTShalf" << endl
	   << "# " << pts.size() << " points." << endl;
      uint l = pts.size();
      Vector2i tdv = dc.twiceDV();
      double dx = ((double) tdv.x()) / 2.0;
      double dy = ((double) tdv.y()) / 2.0;
      for ( uint i = 0; i < l; ++i )
	cout << dx + pts[ i ].x() 
	     << " " << dy + pts[ i ].y() 
	     << " " << ins[ i ] 
	     << " " << dx + pts[ ( i + 1 ) % pts.size() ].x()
	     << " " << dy + pts[ ( i + 1 ) % pts.size() ].y()
	     << endl;
      if ( args.check( "-close" ) )
	cout << dx + pts[ 0 ].x() 
	     << " " << dy + pts[ 0 ].y() 
	     << " " << ins[ 0 ] 
	     << " " << dx + pts[ ( 0 + 1 ) % pts.size() ].x()
	     << " " << dy + pts[ ( 0 + 1 ) % pts.size() ].y()
	     << endl;

    }



  if ( args.check( "-dPointels" ) )
    {
      // DMLPContour::iterator itbegin( dc, dc.firstEdge(),
      // 				     0, dc.getStartPosition(),
      // 				     true );
      DMLPContour::iterator itbegin = dc.first( true );
      DMLPContour::iterator it = itbegin;
      do {
	Vector2i pt = it.pointel();
	cout << pt.x() << " " << pt.y() << endl;
	++it;
      } while ( it != itbegin );
      if ( args.check( "-close" ) )
	{
	  Vector2i pt = it.pointel();
	  cout << pt.x() << " " << pt.y() << endl;
	}

    }
  if ( args.check( "-dRPointels" ) )
    {
      // DMLPContour::iterator itbegin( dc, dc.firstEdge(),
      // 				     0, dc.getStartPosition(),
      // 				     true );
      DMLPContour::reverse_iterator itbegin ( dc.first( true ) );
      DMLPContour::reverse_iterator it = itbegin;
      do {
	Vector2i pt = it.pointel();
	cout << pt.x() << " " << pt.y() << endl;
	++it;
      } while ( it != itbegin );
      if ( args.check( "-close" ) )
	{
	  Vector2i pt = it.pointel();
	  cout << pt.x() << " " << pt.y() << endl;
	}
//       cerr << "--- test position in iterator." << endl;
//       for (uint pos = 0; pos < 100; ++pos )
// 	{
// 	  DMLPContour::iterator it( dc, true, 
// 				    dc.firstEdge(),
// 				    dc.getStartPosition(),
// 				    pos );
// 	  Vector2i pt = it.pointel();
// 	  cerr << pt.x() << " " << pt.y() << endl;
// 	}
      
    }

  return 0;
}
