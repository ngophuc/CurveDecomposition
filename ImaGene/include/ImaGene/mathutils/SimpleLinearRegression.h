/** @file SimpleLinearRegression.h */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// File name : SimpleLinearRegression.h
//
// Creation : 2009/09/01
//
// Version : 2009/09/01
//
// Author : JOL
//
// Summary : Header file for files SimpleLinearRegression.ih and SimpleLinearRegression.cxx
//
// History :
//	2009/09/01 : ?Name? : ?What?
//
// Rcs Id : "@(#)class SimpleLinearRegression declaration."
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if defined(SimpleLinearRegression_RECURSES)
#error Recursive header files inclusion detected in SimpleLinearRegression.h
#else // defined(SimpleLinearRegression_RECURSES)
#define SimpleLinearRegression_RECURSES

#if !defined SimpleLinearRegression_h
#define SimpleLinearRegression_h

//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Forces inline if nothing is provided by the compiler.
#ifndef INLINE
#define INLINE inline
#endif

namespace ImaGene 
{
  
  /////////////////////////////////////////////////////////////////////////////
  // class SimpleLinearRegression
  /////////////////////////////////////////////////////////////////////////////
  /** 
   * Description of class 'SimpleLinearRegression' <p> Aim: Represents
   * a simple linear regression model with 1 regressor constant and 1
   * variable regressor and n data: Y = X B + U, where U follows
   * Gaussian law N(0, sigma^2 I_n).  Y, U are n-vectors, B is a
   * 2-vector, X is the nx2-matrix [(1 x1) ... (1 xn) ] with rank 2.
   *
   * This class can compute the linear regression coefficients and
   * also performs some tests to check if the data corresponds to a
   * linear model.
   */
  class SimpleLinearRegression
  {

    // ----------------------- Static services ------------------------------
  public:

    /**
     * Several tests for the class.
     */
    static bool test();

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor. 
     */
    ~SimpleLinearRegression();

    /**
     * Constructor.
     * The object is empty (and invalid for regression).
     *
     * @param eps_zero the value below which the absolute value of the
     * determinant is considered null.
     */
    SimpleLinearRegression( double eps_zero = 1e-8 );

    /**
     * Clear all datas.
     */
    void clear();

    /**
     * Adds the samples (y,x). Does not compute immediately the
     * regression. See 'computeRegression' for computing the
     * regression with the current samples.
     *
     * @param begin_x an iterator on the first x-data
     * @param end_x an iterator after the last x-data
     * @param begin_y an iterator on the first y-data
     *
     * @see computeRegression
     */
    template <class XIterator, class YIterator>
    void addSamples( XIterator begin_x, XIterator end_x, YIterator begin_y );

    /**
     * Adds the sample (y,x). Does not compute immediately the
     * regression. See 'computeRegression' for computing the
     * regression with the current samples.
     *
     * @param the x data.
     * @param the y data.
     *
     * @see computeRegression
     */
    INLINE void addSample( double x, double y );

    /**
     * Computes the regression of the current parameters.
     *
     * @return 'true' if the regression was valid (non null number of
     * samples, rank of X is 2), 'false' otherwise.
     */
    bool computeRegression();

    /**
     * @return the slope of the linear regression (B1 in Y=B0+B1*X).
     */
    INLINE double slope() const;

    /**
     * @return the intercept of the linear regression (B0 in Y=B0+B1*X).
     */
    INLINE double intercept() const;

    /**
     * Given a new x, predicts its y (hat{y}) according to the linear
     * regression model.
     * 
     * @param x any value.
     * @return the estimated y value, ie hat{y} = B0 + B1*x.
     */
    INLINE double estimateY( double x ) const;

    /**
     * @return the current estimation of the variance of the Gaussian
     * perturbation (i.e. variance of U).
     */
    INLINE double estimateVariance() const;

    /**
     * Given a test confidence value (1-[a]), return the expected interval
     * of value for Y, given a new [x], so that the model is still
     * linear. One may thus check if a new pair (y,x) is still in the
     * current linear model or not.
     *
     * @param x any x value.
     *
     * @param a the expected confidence value for the test (a=0.05
     * means 95% of confidence).
     *
     * @return the expected interval [min_y, max_y] such that any
     * value y within confirms the current linear model.
     */
    std::pair<double,double> trustIntervalForY( double x, double a ) const;

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
  public:
    const double epsilon_zero;

  public:

    unsigned int m_n;
    std::vector<double> m_Y;
    double m_B[ 2 ];
    std::vector<double> m_X;
    std::vector<double> m_U;
    double m_sum_x;
    double m_sum_x2;
    double m_sum_y;
    double m_sum_xy;
    double m_d;
    double m_norm_U2;

    // ------------------------- Hidden services ------------------------------
  protected:

  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    INLINE SimpleLinearRegression( const SimpleLinearRegression & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    INLINE SimpleLinearRegression & operator=( const SimpleLinearRegression & other );
  
    // ------------------------- Internals ------------------------------------
  private:

  
  };

  /**
   * Overloads 'operator<<' for displaying objects of class 'SimpleLinearRegression'.
   * @param that_stream the output stream where the object is written.
   * @param that_object_to_display the object of class 'SimpleLinearRegression' to write.
   * @return the output stream after the writing.
   */
  INLINE std::ostream&
  operator<<( std::ostream & that_stream, 
	      const SimpleLinearRegression & that_object_to_display );

  
} // namespace ImaGene


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods if necessary.
#if defined(INLINE)
#include "ImaGene/mathutils/SimpleLinearRegression.ih"
#endif

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined SimpleLinearRegression_h

#undef SimpleLinearRegression_RECURSES
#endif // else defined(SimpleLinearRegression_RECURSES)
