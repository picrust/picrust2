/** QR Decomposition.
<P>
   For an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
   orthogonal matrix Q and an n-by-n upper triangular matrix R so that
   A = Q*R.
<P>
   The QR decompostion always exists, even if the matrix does not have
   full rank, so the constructor will never fail.  The primary use of the
   QR decomposition is in the least squares solution of nonsquare systems
   of simultaneous linear equations.  This will fail if isFullRank()
   returns false.
*/

#ifndef _BOOST_UBLAS_QRDECOMPOSITION_
#define _BOOST_UBLAS_QRDECOMPOSITION_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/exception.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace boost { namespace numeric { namespace ublas {
    
class QRDecomposition {

    typedef vector<double> Vector;
    typedef matrix<double> Matrix;
        
/* ------------------------
   Class variables
 * ------------------------ */

   /** Array for internal storage of decomposition.
   @serial internal array storage.
   */
   Matrix QR;

   /** Row and column dimensions.
   @serial column dimension.
   @serial row dimension.
   */
   int m, n;

   /** Array for internal storage of diagonal of R.
   @serial diagonal of R.
   */
   Vector Rdiag;

public:
/* ------------------------
   Constructor
 * ------------------------ */

   /** QR Decomposition, computed by Householder reflections.
   @param A    Rectangular matrix
   @return     Structure to access R and the Householder vectors and compute Q.
   */

   QRDecomposition (const Matrix &A);

/* ------------------------
   Public Methods
 * ------------------------ */

   /** Is the matrix full rank?
   @return     true if R, and hence A, has full rank.
   */

   bool isFullRank () const;

   /** Return the Householder vectors
   @return     Lower trapezoidal matrix whose columns define the reflections
   */

   Matrix getH () const;

   /** Return the upper triangular factor
   @return     R
   */

   Matrix getR () const;

   /** Generate and return the (economy-sized) orthogonal factor
   @return     Q
   */

   Matrix getQ () const;

   /** Least squares solution of A*X = B
   @param B    A Matrix with as many rows as A and any number of columns.
   @return     X that minimizes the two norm of Q*R*X-B.
   @exception  IllegalArgumentException  Matrix row dimensions must agree.
   @exception  RuntimeException  Matrix is rank deficient.
   */

   Matrix solve (const Matrix &B);

   /** Matrix inverse
   @return     inverse(A).
   */

   Matrix inverse () {
      return solve(identity_matrix<double>(m,m));
   }
};

}}}
#endif
