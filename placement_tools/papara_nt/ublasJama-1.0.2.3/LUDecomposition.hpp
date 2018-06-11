   /** LU Decomposition.
   <P>
   For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n
   unit lower triangular matrix L, an n-by-n upper triangular matrix U,
   and a permutation vector piv of length m so that A(piv,:) = L*U.
   If m < n, then L is m-by-m and U is m-by-n.
   <P>
   The LU decompostion with pivoting always exists, even if the matrix is
   singular, so the constructor will never fail.  The primary use of the
   LU decomposition is in the solution of square systems of simultaneous
   linear equations.  This will fail if isNonsingular() returns false.
   */

#ifndef _BOOST_UBLAS_LUDECOMPOSITION_
#define _BOOST_UBLAS_LUDECOMPOSITION_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/exception.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace boost { namespace numeric { namespace ublas {
    
class LUDecomposition {

    typedef vector<double> Vector;
    typedef vector<std::size_t> PivotVector;
    typedef matrix<double> Matrix;

/* ------------------------
   Class variables
 * ------------------------ */

   /** Array for internal storage of decomposition.
   @serial internal array storage.
   */
   Matrix LU;

   /** Row and column dimensions, and pivot sign.
   @serial column dimension.
   @serial row dimension.
   @serial pivot sign.
   */
   int m, n, pivsign; 

   /** Internal storage of pivot vector.
   @serial pivot vector.
   */
   PivotVector piv;

public:
/* ------------------------
   Constructor
 * ------------------------ */

   /** LU Decomposition
   @param  A   Rectangular matrix
   @return     Structure to access L, U and piv.
   */

   LUDecomposition (const Matrix& A);

/* ------------------------
   Temporary, experimental code.
   ------------------------ *\

   \** LU Decomposition, computed by Gaussian elimination.
   <P>
   This constructor computes L and U with the "daxpy"-based elimination
   algorithm used in LINPACK and MATLAB.  In Java, we suspect the dot-product,
   Crout algorithm will be faster.  We have temporarily included this
   constructor until timing experiments confirm this suspicion.
   <P>
   @param  A             Rectangular matrix
   @param  linpackflag   Use Gaussian elimination.  Actual value ignored.
   @return               Structure to access L, U and piv.
   *\

   LUDecomposition (Matrix A, int linpackflag);

\* ------------------------
   End of temporary code.
 * ------------------------ */

/* ------------------------
   Public Methods
 * ------------------------ */

   /** Is the matrix nonsingular?
   @return     true if U, and hence A, is nonsingular.
   */

   bool isNonsingular () const {
      for (int j = 0; j < n; j++) {
         if (LU(j,j) == 0.)
            return false;
      }
      return true;
   }

   /** Return lower triangular factor
   @return     L
   */

   Matrix getL () const;

   /** Return upper triangular factor
   @return     U
   */

   Matrix getU () const;

   /** Return pivot permutation vector
   @return     piv
   */

   const PivotVector& getPivot () const {
      return piv;
   }

   /** Return pivot permutation vector as a one-dimensional double array
   @return     (double) piv
   */

   Vector getDoublePivot () const {
      Vector vals(m);
      for (int i = 0; i < m; i++) {
         vals(i) = (double) piv(i);
      }
      return vals;
   }

   /** Determinant
   @return     det(A)
   @exception  bad_size  Matrix must be square
   */

   double det () {
      BOOST_UBLAS_CHECK(m == n, bad_size("Matrix must be square."));
      double d = (double) pivsign;
      for (int j = 0; j < n; j++) {
         d *= LU(j,j);
      }
      return d;
   }

   /** Solve A*X = B
   @param  B   A Matrix with as many rows as A and any number of columns.
   @return     X so that L*U*X = B(piv,:)
   @exception  bad_size Matrix row dimensions must agree.
   @exception  singular  Matrix is singular.
   */

   Matrix solve (const Matrix& B) const;


   /** Matrix pseudoinverse
   @return     pseudoinverse(A).
   */

   Matrix pseudoinverse () {
       return solve(identity_matrix<double>(m,m));
   }
};

}}}
#endif
