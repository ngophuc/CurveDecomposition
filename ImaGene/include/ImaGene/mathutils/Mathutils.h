/** @file Mathutils.h */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// File name : Mathutils.h
//
// Creation : 2004/11/30
//
// Version : 2004/11/30
//
// Author : JOL
//
// Summary : Header file for files Mathutils.ih and Mathutils.cxx
//
// History :
//	2004/11/30 : ?Name? : ?What?
//
// Rcs Id : "@(#)class Mathutils declaration."
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if defined(Mathutils_RECURSES)
#error Recursive header files inclusion detected in Mathutils.h
#else // defined(Mathutils_RECURSES)
#define Mathutils_RECURSES

#if !defined Mathutils_h
#define Mathutils_h

//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include "ImaGene/base/BasicTypes.h"
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Forces inline if nothing is provided by the compiler.
#ifndef INLINE
#define INLINE inline
#endif


namespace ImaGene {
  

///////////////////////////////////////////////////////////////////////////////
// class Mathutils
///////////////////////////////////////////////////////////////////////////////
/** 
 * Description of class 'Mathutils' <p>
 * Aim: A utility class to provide simple math functions.
 */
class Mathutils
{
  // ------------------------- helper classes ------------------------------
public:

  /**
   * A simple class to perform integer computation modulo a given value.
   */
  struct ModuloComputer 
  {
    /**
     * Modulo of all computations.
     */
    uint k;
    
    /**
     * Initializes the modulo computer with the value [m].
     * @param m any non-zero integer.
     */
    INLINE ModuloComputer( uint m );
    
    /**
     * Increment the value [i] modulo.
     * @param i any value between 0 and [k] (excluded).
     * @see k
     */
    INLINE void increment( uint & i ) const;

    /**
     * Decrement the value [i] modulo.
     * @param i any value between 0 and [k] (excluded).
     * @see k
     */
    INLINE void decrement( uint & i ) const;

    /**
     * @param i any value between 0 and [k] (excluded).
     * @return the incremented value of [i] modulo [k].
     * @see k
     */
    INLINE uint next( uint i ) const;

    /**
     * @param i any value between 0 and [k] (excluded).
     * @return the decremented value of [i] modulo [k].
     * @see k
     */
    INLINE uint previous( uint i ) const;

    /**
     * @param i any integer value.
     * @return the value of [i] modulo [k].
     * @see k
     */
    INLINE uint cast( int i ) const;
    
    /**
     * Less comparator modulo. Be careful, modulo comparisons have no
     * sense when the absolute difference of the values are around k / 2.
     *
     * @param i any value between 0 and [k] (excluded).
     * @param j any value between 0 and [k] (excluded).
     * @return 'true' if [i] strictly precedes [j] in a window 'floor([k]/2)'.
     * @see k
     */
    INLINE bool less( uint i, uint j ) const;

    /**
     * Performs j - i modulo, assuming less(i,j) is true.
     *
     * @param j any value between 0 and [k] (excluded).
     * @param i any value between 0 and [k] (excluded).
     * @return the value j - i, always positive. 
     * @see k
     */
    INLINE uint posDiff( uint j, uint i ) const;
    
    
  };



  /**
   * A simple class to perform angle computations. All angles are in [0:2pi[
   */
  struct AngleComputer 
  {
    /**
     * @param i any angle.
     * @return the corresponding angle in [0:2pi[
     */
    INLINE static float cast( float i );
    
    /**
     * Less comparator modulo. Be careful, modulo comparisons have no
     * sense when the absolute difference of the values are around pi.
     *
     * @param i any angle in [0:2pi[
     * @param j any angle in [0:2pi[
     * @return 'true' if [i] strictly precedes [j] in a window 'pi'.
     */
    INLINE static bool less( float i, float j );

    /**
     * Performs j - i modulo 2pi, assuming less(i,j) is true.
     *
     * @param j any angle in [0:2pi[
     * @param i any angle in [0:2pi[
     * @return the value j - i, always positive. 
     */
    INLINE static float posDiff( float j, float i );

    /**
     * Performs j - i, assuming th result is in [-pi:pi[
     *
     * @param j any angle in [0:2pi[
     * @param i any angle in [0:2pi[
     * @return the value j - i, always positive. 
     */
    INLINE static float deviation( float j, float i );

    /**
     * Equivalent to 'less( i, j ) ? i : j'.
     *
     * @param i any angle in [0:2pi[
     * @param j any angle in [0:2pi[
     * @return the smallest angle of [i] and [j] in a window 'pi'.
     */
    INLINE static float min( float i, float j );

    /**
     * Equivalent to 'less( i, j ) ? j : i'.
     *
     * @param i any angle in [0:2pi[
     * @param j any angle in [0:2pi[
     * @return the greatest angle of [i] and [j] in a window 'pi'.
     */
    INLINE static float max( float i, float j );


    /**
     * @param i any angle.
     * @return the corresponding angle in [0:2pi[
     */
    INLINE static double cast( double i );
    
    /**
     * Less comparator modulo. Be careful, modulo comparisons have no
     * sense when the absolute difference of the values are around pi.
     *
     * @param i any angle in [0:2pi[
     * @param j any angle in [0:2pi[
     * @return 'true' if [i] strictly precedes [j] in a window 'pi'.
     */
    INLINE static bool less( double i, double j );

    /**
     * Performs j - i modulo 2pi, assuming less(i,j) is true.
     *
     * @param j any angle in [0:2pi[
     * @param i any angle in [0:2pi[
     * @return the value j - i, always positive. 
     */
    INLINE static double posDiff( double j, double i );

    /**
     * Performs j - i, assuming th result is in [-pi:pi[
     *
     * @param j any angle in [0:2pi[
     * @param i any angle in [0:2pi[
     * @return the value j - i, always positive. 
     */
    INLINE static double deviation( double j, double i );

    /**
     * Equivalent to 'less( i, j ) ? i : j'.
     *
     * @param i any angle in [0:2pi[
     * @param j any angle in [0:2pi[
     * @return the smallest angle of [i] and [j] in a window 'pi'.
     */
    INLINE static double min( double i, double j );

    /**
     * Equivalent to 'less( i, j ) ? j : i'.
     *
     * @param i any angle in [0:2pi[
     * @param j any angle in [0:2pi[
     * @return the greatest angle of [i] and [j] in a window 'pi'.
     */
    INLINE static double max( double i, double j );

    
  };
  

  // ------------------------- Static services ------------------------------
public:

  /**
   * 2*pi.
   */
  static const float two_pi_f;

  /**
   * 2*pi.
   */
  static const double two_pi_d;

  /**
   * @param angle an angle in radian (arbitrary).
   * @return the equivalent angle in [0;2PI[.
   */
  static float angleMod( float angle );

  /**
   * @param angle an angle in radian (arbitrary).
   * @return the equivalent angle in [-PI;PI[.
   */
  static float angleAroundZero( float angle );

  /**
   * @param angle_ref an angle in radian (arbitrary).
   * @param other_angle an other angle in radian (arbitrary).
   * @return the angle equivalent to [other_angle] but expressed in a window [-PI;PI[ wrt [angle_ref].
   */
  static float tune( float angle_ref, float other_angle );

  /**
   * @param angle1 an angle in radian (arbitrary).
   * @param angle2 an angle in radian (arbitrary).
   * @return 'true' if angle1 strictly precedes angle2 when both are brought back in the same window of size PI.
   */
  static bool less( float angle1, float angle2 );


  /**
   * Average length of an angle t: dl(t)=1/(|cos(t)|+|sin(t)|).
   * @param t an angle in radian (arbitrary).
   * @return the corresponding "averaged length" (discrete geometry).
   */
  static float averagedLength( float t );


  // ------------------------- arithmetic services ---------------------------
public:

  /**
   * Returns the simple continued fraction of a/b.
   * 
   * @param z (returns) the partial coefficients.
   * @param a the numerator
   * @param b the denominator
   *
   * @return the gcd of a and b.
   */
  static uint cfrac( std::vector<uint> & z, uint a, uint b );

  /**
   * Returns the greatest common divisor of a and b.
   * 
   * @param a the numerator (possibly negative).
   * @param b the denominator (possibly negative).
   *
   * @return the gcd of a and b.
   */
  static int gcd( int a, int b );


  // ------------------------- random services -----------------------------
public:

  /**
   * @param p any integer below 2^31.
   * @return a random number between 0 and p-1.
   */
  static unsigned int random( unsigned int p );
  
  
  /**
   * @param p any integer below 2^63.
   * @return a random number between 0 and p-1.
   */
  static unsigned long long random( unsigned long long p );

  /**
   * @return a uniform random value between 0.0 (included) and (1.0) excluded.
   */
  static double random1();
    

  // ------------------------- bit services ------------------------------
public:

  /**
   * NB: Relatively constant time and generally faster than the others.
   * @param v any value.
   * @return the most-significant bit set to 1 in [v], (0) if 0
   */
  static int getMSBbyLog( unsigned long long v );

  /**
   * Calls 'getMSBbyLog'.
   * @param v any value.
   * @return the most-significant bit set to 1 in [v], (0) if 0
   */
  static int getMSB( unsigned long long v );


  // ------------------------- prime services ------------------------------
public:

  /**
   * @param n any positive integer.
   * @return a prime number greater than [n].
   */
  static unsigned long long greaterPrime( unsigned long long n );


  // ------------------------- elementary services ----------------------------
public:

  /**
   * @param x any value.
   * @return the square of [x].
   */
  INLINE static float sqr( float x );

  /**
   * @param x any value.
   * @return the square of [x].
   */
  INLINE static double sqr( double x );

  /**
   * Compute x^n with fast exponentiation.
   *
   * @param x the number
   * @param n the exponent
   * @return the value x^n.
   */
  template <typename Number>
  INLINE static
  Number power( Number x, uint n );
  
  // ------------------------- Standard services ------------------------------
public:
  /**
   * Destructor. 
   */
  ~Mathutils();


 

  // ------------------------- Interface --------------------------------------
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

  
  // ------------------------- Datas ------------------------------------------
private:

  /**
   * Size of the array [m_primes].
   * @see m_primes
   */
  static const unsigned int m_size_primes = 33;
  
  /**
   * An array of prime number such that m_primes[ i ] >= 2^i
   */
  static const unsigned long long m_primes[ m_size_primes ];
  

  // ------------------------- Hidden services --------------------------------
protected:
  /**
   * Constructor.
   * Forbidden by default (protected to avoid g++ warnings).
   */
  INLINE Mathutils();
private:
  /**
   * Copy constructor.
   * @param other the object to clone.
   * Forbidden by default.
   */
  INLINE Mathutils( const Mathutils & other );
  /**
   * Assignment.
   * @param other the object to copy.
   * @return a reference on 'this'.
   * Forbidden by default.
   */
  INLINE Mathutils & operator=( const Mathutils & other );
  
  // ------------------------- Internals --------------------------------------
private:
  
};

/**
 * Overloads 'operator<<' for displaying objects of class 'Mathutils'.
 * @param that_stream the output stream where the object is written.
 * @param that_object_to_display the object of class 'Mathutils' to write.
 * @return the output stream after the writing.
 */
INLINE std::ostream&
operator<<( std::ostream & that_stream, const Mathutils & that_object_to_display );


}  // namespace ImaGene




///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods if necessary.
#if defined(INLINE)
#include "ImaGene/mathutils/Mathutils.ih"
#endif

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined Mathutils_h

#undef Mathutils_RECURSES
#endif // else defined(Mathutils_RECURSES)
