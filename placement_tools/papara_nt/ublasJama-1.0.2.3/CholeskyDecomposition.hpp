   /** Cholesky Decomposition.
   <P>
   For a symmetric, positive definite matrix A, the Cholesky decomposition
   is an lower triangular matrix L so that A = L*L'.
   <P>
   If the matrix is not symmetric or positive definite, the constructor
   returns a partial decomposition and sets an internal flag that may
   be queried by the isSPD() method.
   */

#ifndef _BOOST_UBLAS_CHOLESKYDECOMPOSITION_
#define _BOOST_UBLAS_CHOLESKYDECOMPOSITION_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/exception.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace boost { namespace numeric { namespace ublas {
    
class CholeskyDecomposition {

    typedef vector<double> Vector;
    typedef matrix<double> Matrix;

/* ------------------------
   Class variables
 * ------------------------ */

   /** Array for internal storage of decomposition.
   @serial internal array storage.
   */
   Matrix L;

   /** Row and column dimension (square matrix).
   @serial matrix dimension.
   */
   int n;

   /** Symmetric and positive definite flag.
   @serial is symmetric and positive definite flag.
   */
   bool isspd;

public:
/* ------------------------
   Constructor
 * ------------------------ */

   /** Cholesky algorithm for symmetric and positive definite matrix.
   @param  A   Square, symmetric matrix.
   @return     Structure to access L and isspd flag.
   */

   CholeskyDecomposition (const Matrix& A);

/* ------------------------
   Temporary, experimental code.
 * ------------------------ *\

   \** Right Triangular Cholesky Decomposition.
   <P>
   For a symmetric, positive definite matrix A, the Right Cholesky
   decomposition is an upper triangular matrix R so that A = R'*R.
   This constructor computes R with the Fortran inspired column oriented
   algorithm used in LINPACK and MATLAB.  In Java, we suspect a row oriented,
   lower triangular decomposition is faster.  We have temporarily included
   this constructor here until timing experiments confirm this suspicion.
   *\

   \** Array for internal storage of right triangular decomposition. **\
   private transient double(,) R;

   \** Cholesky algorithm for symmetric and positive definite matrix.
   @param  A           Square, symmetric matrix.
   @param  rightflag   Actual value ignored.
   @return             Structure to access R and isspd flag.
   *\

   CholeskyDecomposition (Matrix Arg, int rightflag);

   \** Return upper triangular factor.
   @return     R
   *\

   public Matrix getR () {
      return new Matrix(R,n,n);
   }

\* ------------------------
   End of temporary code.
 * ------------------------ */

/* ------------------------
   Public Methods
 * ------------------------ */

   /** Is the matrix symmetric and positive definite?
   @return     true if A is symmetric and positive definite.
   */

   bool isSPD () const {
      return isspd;
   }

   /** Return triangular factor.
   @return     L
   */

   const Matrix& getL () const {
      return L;
   }

   /** Solve A*X = B
   @param  B   A Matrix with as many rows as A and any number of columns.
   @return     X so that L*L'*X = B
   @exception  bad_size  Matrix row dimensions must agree.
   @exception  singular  Matrix is not symmetric positive definite.
   */

   Matrix solve (const Matrix& B) const;
};

}}}
#endif
