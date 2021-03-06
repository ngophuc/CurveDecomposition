  /**
  * @file MeasureOfStraightLines.cpp
  * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
  * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
  *
  * @date 2010/03/04
  *
  * Implementation of methods defined in MeasureOfStraightLines.h
  *
  * This file is part of the ImaGene library.
  */

  ///////////////////////////////////////////////////////////////////////////////
  #include <cmath>
  #include <assert.h>
  #include "ImaGene/mathutils/MeasureOfStraightLines.h"

  ///////////////////////////////////////////////////////////////////////////////

  using namespace std;

  ///////////////////////////////////////////////////////////////////////////////
  // class MeasureOfStraightLines
  ///////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////
  // Standard services - public :




  ///////////////////////////////////////////////////////////////////////////////
  // Interface - public :


  /**
    * Constructor.
    */
  ImaGene::MeasureOfStraightLines::MeasureOfStraightLines()
  {
    ///Default value
    myEpsilon = 0.005;
  }


  /**
    * Destructor.
    */
  ImaGene::MeasureOfStraightLines:: ~MeasureOfStraightLines()
  {
  }

  /**
  * Writes/Displays the object on an output stream.
  * @param out the output stream where the object is written.
  */
  void
  ImaGene::MeasureOfStraightLines::selfDisplay( std::ostream & out ) const
  {
    out << "[MeasureOfStraightLines]";
  }

  /**
  * Checks the validity/consistency of the object.
  * @return 'true' if the object is valid, 'false' otherwise.
  */
  bool
  ImaGene::MeasureOfStraightLines::isValid() const
  {
    return true;
  }


  ///
  /// Compute the measure associated to the Edge (a0,b0) -- (a1,b1) in the (a,b)-space
  ///
  double
  ImaGene::MeasureOfStraightLines::computeMeasureEdge ( double a0,double b0, double a1, double b1 )
  {
    double measure=0;

    //Aligned with the 'b' axis
    if ( ( b0 == 0 ) and ( b1 == 0 ) )
      return 0;
    //Aligned with the 'b' axis
    if ( ( a0 == 0 ) and ( a1 == 0 ) )
      return 0;

    double delta = ( a0*b1 - a1*b0 );

    if ( a0 == 0 )
      return delta * ( 1.0/ ( 1.0 + sqrt ( 1 + pow ( a1,2 ) ) ) );
    if ( a1 == 0 )
      return delta * ( ( -1.0 + sqrt ( 1.0 + pow ( a0,2 ) ) ) /pow ( a0,2 ) );

    if ( a0 != a1 )
    {
      return  delta * ( ( a0 - a1 + sqrt ( 1 + pow ( a0,2 ) ) *a1 - a0*sqrt ( 1 + pow ( a1,2 ) ) ) / ( pow ( a0,2 ) *a1 - a0*pow ( a1,2 ) ) );
    }
    else
      return  delta * a1/ ( a0* ( 1 + pow ( a1,2 ) + sqrt ( 1 + pow ( a1,2 ) ) ) );
  }


  /**
  * Set the epsilon threshold in the numerical approximation.
  *
  * @param aValue the new epsilon value
  */
  void
  ImaGene::MeasureOfStraightLines::setEpsilon(const double aValue)
  {
    myEpsilon = aValue;
  }

  ///////////////////////////////////////////////////////////////////////////////
  // Implementation of inline methods                                          //

  /**
  * Compute the measure of the polygon {(a_i,b_i)} in the (a,b)-parameter space
  * @param a the a-value of polygon vertices
  * @param b the b-value of polygon vertices
  *
  * @return the measure value (positive value)
  */
  double
  ImaGene::MeasureOfStraightLines::computeMeasure(const std::vector<double> &a,
  const std::vector<double> &b)
  {
    double measure = 0;

    assert(a.size() == b.size());

    for (unsigned int i=0 ; i < a.size() ; i++)
    {
      measure += computeMeasureEdge ( a[i], b[i],
	      a[ ( i+1 ) % a.size()],b[ ( i+1 ) %a.size()] );
    }

    return measure;
  }

  /**
  * Compute the abscissa of the centroid of the polygon {(a_i,b_i)} in the (a,b)-parameter space
  * @param a the a-value of polygon vertices
  * @param b the b-value of polygon vertices
  *
  * @return the measure value (positive value)
  */
  double
  ImaGene::MeasureOfStraightLines::computeCentroidA(const std::vector<double> &a,
    const std::vector<double> &b)
    {
      double measure = 0;
      double C_a = 0;

      assert(a.size() == b.size());

      for (unsigned int i=0 ; i < a.size() ; i++)
      {
	measure += computeMeasureEdge( a[i], b[i],
	      a[ ( i+1 ) % a.size()],b[ ( i+1 ) %a.size()] );
	      C_a += computeCentroidEdge_a ( a[i], b[i],
		      a[ ( i+1 ) % a.size()],b[ ( i+1 ) %a.size()] );
      }

      return C_a/measure;
    }

    /**
    * Compute the abscissa of the centroid of the polygon {(a_i,b_i)} in the (a,b)-parameter space
    * @param a the a-value of polygon vertices
    * @param b the b-value of polygon vertices
    *
    * @return the measure value (positive value)
    */
    double
    ImaGene::MeasureOfStraightLines::computeCentroidB(const std::vector<double> &a,
      const std::vector<double> &b)
      {
	double measure = 0;
	double C_b = 0;

	assert(a.size() == b.size());

	for (unsigned int i=0 ; i < a.size() ; i++)
	{
	  measure += computeMeasureEdge( a[i], b[i],
		a[ ( i+1 ) % a.size()],b[ ( i+1 ) %a.size()] );
		C_b += computeCentroidEdge_b ( a[i], b[i],
			a[ ( i+1 ) % a.size()],b[ ( i+1 ) %a.size()] );
	}

	return C_b/measure;
      }


      ///
      /// Compute the centroid on 'b' on the rectangular domain with vertices (x1,,y1) - (x2,y2)
      /// PRECONDITION: y1<y2
      ///
      ///

      double
      ImaGene::MeasureOfStraightLines::__computeCentroidSquare_b ( double x1, double y1, double x2,double y2 )
      {
	double val;
	val = ( ( -1.0* ( sqrt ( 1.0 + pow ( x1,2 ) ) *x2 ) + x1*sqrt ( 1 + pow ( x2,2 ) ) ) * ( pow ( y1,2 ) - pow ( y2,2 ) ) ) / ( 2.0*sqrt ( ( 1.0 + pow ( x1,2 ) ) * ( 1.0 + pow ( x2,2 ) ) ) );

	assert ( val>=0 );
	return val;
      }



      int
      ImaGene::MeasureOfStraightLines::sign ( double a )
      {
	if ( a>=0 )
	  return 1;
	else
	  return -1;
      }


      ///
      /// Approximate  the centroid on 'b' on the trapezioid  (a0,0)-(a0,b0)-(a1,b1)-(a1,0)
      /// (internal function)
      ///

      double
      ImaGene::MeasureOfStraightLines::__computeCentroidEdgeApprox_b ( double a0, double b0,double a1,double b1 )
      {
	double measure,y1,y2;
	int nb_step;
	/// A0 -> A1
	if ( a1 < a0 )
	{
	  double tmp = a0;
	  a0 = a1;
	  a1 = tmp;
	  tmp = b0;
	  b0 = b1;
	  b1 = tmp;
	}

	measure = 0;
	nb_step = ( int ) floor ( fabs ( a1-a0 ) /myEpsilon );
	double slope = ( b1-b0 ) / ( a1-a0 );
	double decal = b1 - ( b1-b0 ) / ( a1-a0 ) *a1;

	if ( slope == 0 )
	  return __computeCentroidSquare_b ( a0, 0, a1, b0 );

	for ( unsigned int i=0; i < nb_step ; i++ )
	{
	  y2 = ( b1-b0 ) / ( a1-a0 ) * ( a0+ ( i+1 ) *myEpsilon ) + decal;
	  if ( y2>0 )
	    measure += __computeCentroidSquare_b ( a0+i*myEpsilon, 0, a0+ ( i+1 ) *myEpsilon, y2 );
	  else
	    measure += __computeCentroidSquare_b ( a0+i*myEpsilon, y2, a0+ ( i+1 ) *myEpsilon, 0 );
	}
	return measure;
      }

      ///
      /// Approximate  the centroid on 'b' on the triangle (0,0)-(a0,b0)-(a1,b1)
      /// (internal function)
      ///

      double
      ImaGene::MeasureOfStraightLines::__computeCentroidTriApprox_b ( double a0, double b0,double a1,double b1 )
      {
	int signe = sign ( a0*b1 - a1*b0 );
	double measure = 0;

	double m_A0A1=0, m_0A0=0, m_0A1=0;
	if ( ( b0 != 0 ) and ( a0 != 0 ) )
	  m_0A0 = __computeCentroidEdgeApprox_b ( 0,0,a0,b0 );
	if ( ( b1 != 0 )  and ( a1 != 0 ) )
	  m_0A1 = __computeCentroidEdgeApprox_b ( 0,0,a1,b1 );
	if ( ( a1 != a0 ) )
	  m_A0A1 = __computeCentroidEdgeApprox_b ( a0,b0,a1,b1 );

	assert ( m_0A0>=0 );
	assert ( m_0A1>=0 );
	assert ( m_A0A1>=0 );

	if ( a0 < a1 )
	{
	  double det = ( a1 * ( b0 - b1 ) - b1* ( a0 - a1 ) );
	  if ( det >0 )
	    measure = m_0A0 + m_A0A1 - m_0A1;
	  else
	    measure = m_0A1 - m_0A0 - m_A0A1;
	}
	else
	{
	  double det = ( a0 * ( b1 - b0 ) - b0* ( a1 - a0 ) );
	  if ( det >0 )
	    measure = m_0A1 + m_A0A1 - m_0A0;
	  else
	    measure = m_0A0 - m_0A1 - m_A0A1;
	}
	assert ( measure>=0 );
	return signe*measure;
      }

      ///
      /// Compute the centroid on 'b' on the triangle (0,0)-(a0,b0)-(a1,b1)
      ///

      double
      ImaGene::MeasureOfStraightLines::computeCentroidEdge_b ( double a0,double b0, double a1, double b1 )
      {
	double measure=0;
	double delta= ( a0*b1 - a1*b0 );

	if ( ( b0 == 0 ) and ( b1 == 0 ) )
	  return 0;
	if ( ( a0 == 0 ) and ( a1 == 0 ) )
	  return 0;

	if ( a0 == 0 )
	  return delta* ( a1* ( ( -2 + sqrt ( 1 + a1*a1 ) ) *b0 + 2*b1 ) + ( b0 - 2*b1 ) *asinh ( a1 ) ) / ( 2*a1*a1*a1 );
	if ( a1 == 0 )
	  return delta* ( a0* ( 2*b0 + ( -2 + sqrt ( 1 + a0*a0 ) ) *b1 ) + ( -2*b0 + b1 ) *asinh ( a0 ) ) / ( 2*a0*a0*a0 );

	if ( a0 == a1 )
	  return delta* ( ( - ( ( a1* ( b0 + a1* ( -a1 + sqrt ( 1 + pow ( a1,-2 ) ) *sqrt ( 1 + pow ( a1,2 ) ) ) *
	  b1 ) ) /sqrt ( 1.0 + pow ( a1,2 ) ) ) + ( b0 + b1 ) *asinh ( a1 ) ) /
	  ( 2.*pow ( a1,3 ) ) );

	return __computeCentroidTriApprox_b ( a0,b0,a1,b1 );

      }

      ///
      /// Compute the centroid on 'a' on the triangle (0,0)-(a0,b0)-(a1,b1)
      ///

      double
      ImaGene::MeasureOfStraightLines::computeCentroidEdge_a ( double a0,double b0, double a1, double b1 )
      {

	double measure=0;
	double delta= ( a0*b1 - a1*b0 );

	if ( ( b0 == 0 ) and ( b1 == 0 ) )
	  return 0;
	if ( ( a0 == 0 ) and ( a1 == 0 ) )
	  return 0;

	if ( a0 == 0 )
	  return delta* ( a1 - asinh ( a1 ) ) / ( a1*a1 );
	if ( a1 == 0 )
	  return delta* ( a0 - asinh ( a0 ) ) / ( a0*a0 );

	if ( a0 == a1 )
	  return delta* ( ( -sqrt ( 1 + pow ( a1,-2 ) ) - 1.0/ ( a1*sqrt ( 1 + pow ( a1,2 ) ) ) + a1/sqrt ( 1 + pow ( a1,2 ) ) +
	  ( 2*asinh ( a1 ) )  /pow ( a1,2 ) ) /2. );

	return delta* ( ( - ( a1*asinh ( a0 ) ) + a0*asinh ( a1 ) ) / ( a0* ( a0 - a1 ) *a1 ) );
      }




      ///////////////////////////////////////////////////////////////////////////////
      // Implementation of inline functions and external operators                 //

      /**
      * Overloads 'operator<<' for displaying objects of class 'MeasureOfStraightLines'.
      * @param out the output stream where the object is written.
      * @param object the object of class 'MeasureOfStraightLines' to write.
      * @return the output stream after the writing.
      */
      std::ostream&
      ImaGene::operator<<( std::ostream & out,
  const MeasureOfStraightLines & object )
  {
    object.selfDisplay( out );
    return out;
  }

  //                                                                           //
  ///////////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////////
  // Internals - private :

  //                                                                           //
  ///////////////////////////////////////////////////////////////////////////////
