/** @file FreemanChain.h */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// File name : FreemanChain.h
//
// Creation : 2008/05/20
//
// Version : 2008/05/20
//
// Author : JOL
//
// Summary : Header file for files FreemanChain.ih and FreemanChain.cxx
//
// History :
//	2008/05/20 : ?Name? : ?What?
//
// Rcs Id : "@(#)class FreemanChain declaration."
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if defined(FreemanChain_RECURSES)
#error Recursive header files inclusion detected in FreemanChain.h
#else // defined(FreemanChain_RECURSES)
#define FreemanChain_RECURSES

#if !defined FreemanChain_h
#define FreemanChain_h

//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <string>
#include <vector>
#include "ImaGene/base/BasicTypes.h"
#include "ImaGene/base/OrderedAlphabet.h"
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/dgeometry2d/C4CIterator.h"
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Forces inline if nothing is provided by the compiler.
#ifndef INLINE
#define INLINE inline
#endif

namespace ImaGene 
{
  
  /////////////////////////////////////////////////////////////////////////////
  // class FreemanChain
  /////////////////////////////////////////////////////////////////////////////
  /** 
   * Description of class 'FreemanChain' <p> Aim: Describes a digital
   * 4-connected contour as a string of '0', '1', '2', and '3' and the
   * coordinate of the first point. When it is a loop, it is the
   * counterclockwise boundary of the shape.
   */
  class FreemanChain
  {
    // ------------------------- iterator ------------------------------
  public:
    /**
     * This class represents an iterator on the freeman chain, storing
     * the current coordinate.
     */
    class const_iterator
    {
      // ------------------------- data -----------------------
    private:
      /**
       * The Freeman chain visited by the iterator.
       */
      const FreemanChain* m_fc;

      /**
       * The current position in the word.
       */
      uint m_pos;

      /**
       * The current coordinates of the iterator.
       */
      Vector2i m_xy;
      

      // ------------------------- Standard services -----------------------
    public:
      /**
       * Default Constructor.
       * The object is not valid.
       */
      INLINE const_iterator();

      /**
       * Constructor.
       * Nb: complexity in O(n).
       *
       * @param chain a Freeman chain,
       * @param n the position in [chain] (within 0 and chain.size()-1).
       */
      INLINE const_iterator( const FreemanChain & chain, uint n = 0 );

      /**
       * Copy constructor.
       * @param other the iterator to clone.
       */
      INLINE const_iterator( const const_iterator & other );
    
      /**
       * Assignment.
       * @param other the iterator to copy.
       * @return a reference on 'this'.
       */
      INLINE const_iterator& operator=( const const_iterator & other );
    
      /**
       * Destructor. Does nothing.
       */
      INLINE ~const_iterator();
    
      // ------------------------- iteration services -------------------------
    public:
      
      /**
       * @return the current coordinates.
       */
      INLINE Vector2i operator*() const;

      /**
       * @return the current coordinates.
       */
      INLINE Vector2i get() const;

      /**
       * Pre-increment.
       * Goes to the next point on the chain.
       */
      INLINE const_iterator& operator++();
      
      /**
       * Goes to the next point on the chain.
       */
      INLINE void next();

      /**
       * Goes to the next point on the chain as if on a loop.
       */
      INLINE void nextInLoop();

      /**
       * @return the current position (as an index in the Freeman chain).
       */
      INLINE uint getPosition() const;

      /**
       * @return the associated Freeman chain.
       */
      INLINE const FreemanChain* getChain() const;

      /**
       * @return the current Freeman code (specifies the movement to
       * the next point).
       */
      INLINE uint getCode() const;

      /**
       * Pre-decrement.
       * Goes to the previous point on the chain.
       */
      INLINE const_iterator& operator--();
      
      /**
       * Goes to the previous point on the chain if possible.
       */
      INLINE void previous();

      /**
       * Goes to the previous point on the chain as if on a loop.
       */
      INLINE void previousInLoop();

      /**
       * Equality operator.
       *
       * @param other the iterator to compare with (must be defined on
       * the same chain).
       *
       * @return 'true' if their current positions coincide.
       */
      INLINE bool operator==( const const_iterator & other ) const;

      /**
       * Inequality operator.
       *
       * @param other the iterator to compare with (must be defined on
       * the same chain).
       *
       * @return 'true' if their current positions differs.
       */
      INLINE bool operator!=( const const_iterator & other ) const;

      /**
       * Inferior operator.
       *
       * @param other the iterator to compare with (must be defined on
       * the same chain).
       *
       * @return 'true' if the current position of 'this' is before
       * the current position of [other].
       */
      INLINE bool operator<( const const_iterator & other ) const;
      
    };
      


    // ------------------------- static services ------------------------------
  public:

    /**
     * Outputs the chain [c] to the stream [out].
     * @param out any output stream,
     * @param c a Freeman chain.
     */
    static void write( std::ostream & out, const FreemanChain & c );

    /**
     * Reads a chain from the stream [in] and updates [c].
     * @param in any input stream,
     * @param c (returns) the Freeman chain.
     */
    static void read( std::istream & in, FreemanChain & c );

    /**
     * Creates a Freeman chaincode [chain] and a chain coding the
     * quadrant of each step [qchain], given an iterator on a
     * 4-connected path [it], a number of step [nb], the first step
     * [freeman_code] and the first quadrant [quadrant].
     *
     * @param chain (returns) a string of '0', '1', '2', or '3'
     * (Freeman chaincode)
     *
     * @param qchain (returns) a string of '0', '1', '2', or '3'
     * giving the quadrants.
     *
     * @param it any iterator 
     *
     * @param nb the number of 'next' performed on a copy of [it].
     *
     * @param freeman the first code or step.
     *
     * @param quadrant the first quadrant (equal to freeman_code or one below).
     *
     * @param start_index the starting index in [chain] and [qchain],
     * default is 0.
     */
    static void create( std::string & chain,
			std::string & qchain,
			const C4CIterator & it, 
			uint nb, uint freeman, uint quadrant,
			uint start_index = 0 );

    /**
     * @param zero (returns) the '0' or 'x' letter for quadrant [quadrant].
     * @param one (returns) the '1' or 'y' letter for quadrant [quadrant].
     * @param quadrant the quadrant as any of '0', '1', '2', or '3'.
     */
    static void alphabet( char & zero, char & one, char quadrant );

    /**
     * Given two consecutive moves on a Freeman chain code, this
     * method returns the type of movement: 0: return move, 1: turning
     * toward the interior, 2: going straight, 3: turning toward
     * exterior. Interior/exterior is specified by [ccw].
     *
     * @param code1 the code of the first step as an integer in 0..3.
     * @param code2 the code of the second step as an integer in 0..3.
     * @param ccw 'true' if the contour is seen counterclockwise with
     * its inside to the left.
     */
    static uint movement( uint code1, uint code2, bool ccw = true );

    /**
     * Returns the displacement vector of a Freeman code.
     *
     * @param dx (returns) the x-displacement.
     * @param dy (returns) the y-displacement.
     * @param c the code.
     */
    static void displacement( int & dx, int & dy, uint code );

    /**
     * @param c a Freeman code (between 0-3).
     * Returns the displacement vector of the Freeman code.
     */
    static Vector2i displacement( uint code );

    /**
       Given a 2D vector, return its basis.
       The basis is given by the returned two letters a<b.
       @param a (modified) the small letter in 0-3 (horizontal axis)
       @param b (modified) the big letter in 0-3 (vertical axis)
       @param v any non-null 2D vector.
       @param inner 'true' if the basis is related to the inner contour. 
       @param cw 'true' if the basis is related to a clockwise contour (inside to the right). 
    */
    static void basis( uint & a, uint & b, 
		       const Vector2i & v, bool inner, bool cw  );

    /**
       Given two letters forming a basis, returns the quadrant
       (ie. the smallest letter). This function returns the same value
       as Vector2i::quadrant for vectors that are not aligned with
       axes. Otherwise, this function better handle those cases.

       @code
       Vector2i v(4,0);
       uint a, b;
       bool ccw;
       FreemanChain::basis( a, b, true, true );
       uint q = quadrant( ccw, a, b );

       @param[out] ccw returns 'true' if a < b, false otherwise.
       @param[in] a the first (small) letter of the basis (in 0,1,2,3).
       @param[in] b the second (big) letter of the basis (in 0,1,2,3).
       @return the quadrant number in 0,1,2,3.
    */
    static uint quadrant( bool & ccw, uint a, uint b );

    /**
     * @param code any Freeman code.
     *
     * @param ccw when 'true' turns counterclockwise (or left),
     * otherwise turns clockwise (right).
     *
     * @return the turned code.
     */
    static uint turnedCode( uint code, bool ccw = true );

    /**
     * @param code any Freeman code.
     *
     * @return the opposite code.
     */
    static uint oppositeCode( uint code );
    
    /**
     * From the Freeman chain [pl_chain] representing a pointel
     * 4-connected contour, constructs the Freeman chain [pix_chain]
     * that represents its inner 4-connected border of pixels. The
     * Freeman chain [pl_chain] has its inside to the left (ie. ccw).
     * 
     * Note that chain codes going back and forth are @b never considered
     * useless: it means that the chain is always supposed to have its
     * interior to the left (ccw) or right (cw) even at configurations
     * "02", "13", "20", "31".
     * 
     * @param pix_chain (output) the code of the 4-connected inner border. 
     *
     * @param pl2pix (output) the mapping associating pointels to
     * pixels as indices in their respective Freeman chain.
     *
     * @param pix2pl (output) the inverse mapping associating pixels to
     * pointels as indices in their respective Freeman chain.
     *
     * @param pl_chain the input code of the 4-connected pointel contour.
     */
    static void pointel2pixel( FreemanChain & pix_chain,
			       std::vector<uint> & pl2pix,
			       std::vector<uint> & pix2pl,
			       const FreemanChain & pl_chain );

    /**
     * From the Freeman chain [outer_chain] representing a 4-connected
     * contour, constructs the Freeman chain [inner_chain] that
     * represents its inner 4-connected contour (which lies in its
     * interpixel space). The boolean [ccw] specifies if the inside is
     * to the left (ccw) or to the right (cw).
     * 
     * Note that chain codes going back and forth are @b never considered
     * useless: it means that the chain is always supposed to have its
     * interior to the left (ccw) or right (cw) even at configurations
     * "02", "13", "20", "31".
     * 
     * @param inner_chain (output) the code of the 4-connected inner
     * border, with starting coordinates that are floored to the closest
     * integer.
     *
     * @param outer2inner (output) the mapping associating outer to
     * inner elements as indices in their respective Freeman chain.
     *
     * @param inner2outer (output) the mapping associating inner to
     * outer elements as indices in their respective Freeman chain.
     *
     * @param outer_chain the input code of the 4-connected contour.
     *
     * @param ccw 'true' if the contour is seen counterclockwise with
     * its inside to the left.
     */
    static void innerContour( FreemanChain & inner_chain,
			      std::vector<uint> & outer2inner,
			      std::vector<uint> & inner2outer,
			      const FreemanChain & outer_chain,
			      bool ccw = true );

    /**
     * Reads the 4-connected contour [c] so that meaningless back and
     * forth steps are removed. These operations may create one or
     * several 4-connected contours (stored in [clean_cs]), whether
     * these removals cuts the contour in several loops. Because of
     * that, the mappings are more complex.

     * @param clean_cs (output) the array of cleaned 4-connected contours.
     *
     * @param c2clean (output) the mapping associating an element to
     * its clean element as a pair (n,i) where n is the index of the
     * cleaned contour and i the indice of the element in this Freeman
     * chain.
     *
     * @param clean2c (output) the array of mapping associating a
     * clean element to its non-clean element. clean2c[n][j] gives the
     * index of the non-clean element on c corresponding to the clean
     * element of index j in the n-th contour.
     *
     * @param c the input code of the 4-connected contour.
     *
     * @param ccw 'true' if the contour is seen counterclockwise with
     * its inside to the left.
     *
     * @todo This method is not implemented.
     */
    static void cleanContour( std::vector<FreemanChain> & clean_cs,
			      std::vector< std::pair<uint,uint> > & c2clean,
			      std::vector< std::vector<uint> > & clean2c,
			      const FreemanChain & c,
			      bool ccw = true );
    /**
     * Removes outer spikes along a 4-connected contour, meaning steps
     * "02", "13", "20" or "31", which point outside the shape. The
     * inside is given by parameter [ccw]. Note that 4-connected
     * pointel contours should not have any outer spikes, while
     * 4-connected pixel contours should not have any inner spikes.
     *
     * @param clean_c (output) the cleaned 4-connected contour.
     *
     * @param c2clean (output) the mapping associating an element to
     * its clean element.
     *
     * @param clean2c (output) the inverse mapping associating a
     * clean element to its non-clean element. 
     *
     * @param c the input code of the 4-connected contour (should be a loop !).
     *
     * @param ccw 'true' if the contour is seen counterclockwise with
     * its inside to the left.
     *
     * @return 'true' if the contour add an interior, 'false' otherwise.
     */
    static bool cleanOuterSpikes( FreemanChain & clean_c,
				  std::vector<uint> & c2clean,
				  std::vector<uint> & clean2c,
				  const FreemanChain & c,
				  bool ccw = true );

    /**
     * Given a Freeman chain [c] coding a 4-connected pixel loop, computes
     * its subsampling by the transformation:
     * X = ( x - x0 ) div h, 
     * Y = ( y - y0 ) div v.
     *
     * @param subc (output) the subsampled Freeman chain code (may
     * contain spikes)
     * 
     * @param c2subc (output) the mapping associating an element to
     * its subsampled element.
     *
     * @param subc2c (output) the inverse mapping associating a
     * subsampled element to its element. More precisely, subc2c[ j ]
     * is the last pointel to be in j.
     *
     * @param c the input chain code.
     *
     * @param h the subsampling along x
     * @param v the subsampling along y
     * @param x0 the x-origin of the frame (X,Y) in (x,y)
     * @param y0 the y-origin of the frame (X,Y) in (x,y)
     *
     * @return 'false' if initial contour was empty or if [subc] is empty,
     * 'true' otherwise.
     */
    static bool subsample( FreemanChain & subc,
			   std::vector<uint> & c2subc,
			   std::vector<uint> & subc2c,
			   const FreemanChain & c,
			   uint h, uint v,
			   int x0, int y0 );

    /**
     * Computes the Minimum Length Polygon (MLP) of the Freeman chain [fc].
     *
     * @param vx (returns) the x-coordinates of the MLP vertices.
     * @param vy (returns) the y-coordinates of the MLP vertices.
     * @param vi (returns) the indices of the MLP vertices in [fc].
     * @param fc the input Freeman chain.
     *
     * @return 'true' if the MLP has made a correct loop, 'false'
     * otherwise (error ?).
     */
    static double computeMLP( std::vector<int> & vx,
			      std::vector<int> & vy,
			      std::vector<uint> & vi,
			      const FreemanChain & fc );

    /**
     * Computes the Minimum Length Polygon (MLP) of the Freeman chain
     * [fc]. Tells also if vertices are inside or outside.
     *
     * @param vx (returns) the x-coordinates of the MLP vertices.
     * @param vy (returns) the y-coordinates of the MLP vertices.
     * @param vi (returns) the indices of the MLP vertices in [fc].
     *
     * @param vt (returns) the type (inside=true, outside=false) of the
     * MLP vertices in [fc].
     *
     * @param fc the input Freeman chain.
     *
     * @param cw when 'true', the inside is considered to be to the right
     * of the contour (the contour turns clockwise around the inside),
     * when 'false' the inside is considered to be to the left of the
     * contour (the contour turns counterclockwise around the inside).
     *
     * @return 'true' if the MLP has made a correct loop, 'false'
     * otherwise (error ?).
     */
    static bool computeMLP( std::vector<int> & vx,
			    std::vector<int> & vy,
			    std::vector<uint> & vi,
			    std::vector<bool> & vt,
			    const FreemanChain & fc,
			    bool cw );

    /**
     * Computes the Minimum Length Polygon (MLP) of the Freeman chain
     * [fc]. Tells also if vertices are inside or outside. The MLP lies in
     * the half-integer plane compared with the Freeman Chain. Vertex
     * coordinates are given as integers. This function gives the
     * displacement vector to go from the digital contour to the MLP.
     *
     * @param vx (returns) the x-coordinates of the MLP vertices.
     * @param vy (returns) the y-coordinates of the MLP vertices.
     * @param vi (returns) the indices of the MLP vertices in [fc].
     *
     * @param vt (returns) the type (inside=true, outside=false) of the
     * MLP vertices in [fc].
     *
     * @param twice_dv (returns) twice the displacement vector to go from
     * the integer plane of the Freeman chain to the half-integer plane of
     * the MLP.
     *
     * @param fc the input Freeman chain.
     *
     * @param cw when 'true', the inside is considered to be to the right
     * of the contour (the contour turns clockwise around the inside),
     * when 'false' the inside is considered to be to the left of the
     * contour (the contour turns counterclockwise around the inside).
     *
     * @return 'true' if the MLP has made a correct loop, 'false'
     * otherwise (error ?).
     */
    static const_iterator computeMLP( std::vector<int> & vx,
				      std::vector<int> & vy,
				      std::vector<uint> & vi,
				      std::vector<bool> & vt,
				      Vector2i & twice_dv,
				      const FreemanChain & fc,
				      bool cw );

			    
    /**
     * Computes the Minimum Length Polygon (MLP) of the Freeman chain [fc]
     * and return its length.
     *
     * @param fc the input Freeman chain.
     *
     * @return the Euclidean length of the MLP.
     */
    static double lengthMLP( const FreemanChain & fc );
    static double testCMLP( const FreemanChain & fc );





    /**
     * Return a vector containing all the interger points of the freemanchain.
     *
     * @param fc the FreemanChain
     * @param vContour (returns) the vector containing all the integer contour points.
     */
    static void getContourPoints(const FreemanChain & fc, std::vector<Vector2i> & vContour); 
    
    /**
     * Equality operator.
     *
     * @param other the FreemanChain to compare with.
     *
     * @return 'true' if the two FreemanChain describe the same contour
     *         (whatever the initial position).
     */
    bool operator==( const FreemanChain & other ) const;


    // ------------------------- public datas ---------------------------------
  public:
    /**
     * The chain code.
     */
    std::string chain;

    /**
     * the x-coordinate of the first point.
     */
    int x0;

    /**
     * the y-coordinate of the first point.
     */
    int y0;
    
    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor. 
     */
    ~FreemanChain();

    /**
     * Constructor.
     * @param s the chain code.
     * @param x the x-coordinate of the first point.
     * @param y the y-coordinate of the first point.
     */
    INLINE FreemanChain( const std::string & s = "", int x = 0, int y = 0 );

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    INLINE FreemanChain( const FreemanChain & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    INLINE FreemanChain & operator=( const FreemanChain & other );


    /**
     * Iterator service.
     * @return an iterator pointing on the first point of the chain.
     */
    INLINE FreemanChain::const_iterator begin() const;

    /**
     * Iterator service.
     * @return an iterator pointing after the last point of the chain.
     */
    INLINE FreemanChain::const_iterator end() const;

    /**
     * @param pos a position in the chain code.
     * @return the code at position [pos].
     */ 
    INLINE uint code( uint pos ) const;

    /**
     * @param pos a position in the chain code.
     * @return the next position.
     */ 
    INLINE uint next( uint pos ) const;

    /**
     * @param pos a position in the chain code.
     * @return the previous position.
     */ 
    INLINE uint previous( uint pos ) const;

    /**
     * @return the length of the Freeman chain code.
     */
    INLINE uint size() const;

    /**
     * Computes a bounding box for the Freeman chain code.
     *
     * @param min_x (returns) the minimal x-coordinate.
     * @param min_y (returns) the minimal y-coordinate.
     * @param max_x (returns) the maximal x-coordinate.
     * @param max_y (returns) the maximal y-coordinate.
     */
    void computeBoundingBox( int32 & min_x, int32 & min_y, 
			     int32 & max_x, int32 & max_y ) const;
    /**
     * Computes a bounding box for the Freeman chain code.
     *
     * @param min_x (returns) the minimal x-coordinate.
     * @param min_y (returns) the minimal y-coordinate.
     * @param max_x (returns) the maximal x-coordinate.
     * @param max_y (returns) the maximal y-coordinate.
     */
    void computeBoundingBox( int64 & min_x, int64 & min_y, 
			     int64 & max_x, int64 & max_y ) const;

    /**
     * Finds a quadrant change in 'this' Freeman chain and returns the
     * position as an iterator. A quadrant change is some
     <code>
     abb..bc
      |
     iterator
     <endcode>
     *
     * The alphabet is possibly re-ordered so that a > b > c.
     *
     * @param A (possibly updated) a Freeman chain alphabet, possibly
     * re-ordered so that a > b > c.
     *
     * @return an iterator on 'this' that points on the first letter b.
     */
    FreemanChain::const_iterator
    findQuadrantChange( OrderedAlphabet & A ) const;

    /**
     * Finds a quadrant change in 'this' Freeman chain and returns the
     * position as an iterator. A quadrant change is some
     @code
     [abc]*bc...cd
            |
     iterator
     @endcode
     *
     * This quadrant change also guarantees that is not a place where a
     * convexity change occurs in the combinatorial MLP algorithm.
     *
     * The alphabet is possibly re-ordered so that b > c > d > a.
     *
     * @param A (possibly updated) a Freeman chain alphabet, possibly
     * re-ordered so that b > c > d > a.
     *
     * @return an iterator on 'this' that points on the first letter c.
     */
    FreemanChain::const_iterator
    findQuadrantChange4( OrderedAlphabet & A ) const;

    /**
     * This method takes O(n) operations and works only for Freeman
     * chains whose successive codes are between +1/-1. It determines
     * if the FreemanChain corresponds to a closed contour, and if
     * this is the case, determines how many counterclockwise loops the
     * contour has done. Of course, it the contour has done
     * clockwise loops, then the given number is accordingly
     * negative.
     *
     * @return the number of counterclockwise loops, or '0' is the contour
     * is open or invalid.
     */
    int isClosed() const;

    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param that_stream the output stream where the object is written.
     */
    void selfDisplay( std::ostream & that_stream ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool OK() const;
  

    // ------------------------- Datas ----------------------------------------
  private:


    // ------------------------- Hidden services ------------------------------
  protected:


  private:

  
    // ------------------------- Internals ------------------------------------
  private:
  
  };

  /**
   * Overloads 'operator<<' for displaying objects of class 'FreemanChain'.
   * @param that_stream the output stream where the object is written.
   * @param that_object_to_display the object of class 'FreemanChain' to write.
   * @return the output stream after the writing.
   */
  INLINE std::ostream&
  operator<<( std::ostream & that_stream, 
	      const FreemanChain & that_object_to_display );

  
} // namespace ImaGene


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods if necessary.
#if defined(INLINE)
#include "ImaGene/dgeometry2d/FreemanChain.ih"
#endif

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined FreemanChain_h

#undef FreemanChain_RECURSES
#endif // else defined(FreemanChain_RECURSES)
