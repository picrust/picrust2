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

#include <algorithm>
#include <cmath>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "LUDecomposition.hpp"

namespace boost { namespace numeric { namespace ublas {

/* ------------------------
   Constructor
 * ------------------------ */

   /** LU Decomposition
   @param  A   Rectangular matrix
   @return     Structure to access L, U and piv.
   */

LUDecomposition::LUDecomposition (const Matrix& A) {

   // Use a "left-looking", dot-product, Crout/Doolittle algorithm.

      assert( A.size1() < size_t(std::numeric_limits<int>::max()) && A.size2() < size_t(std::numeric_limits<int>::max()));

      LU = A;
      m = int(A.size1());
      n = int(A.size2()); // TODO: cast to int
      piv = PivotVector(m);
      for (int i = 0; i < m; i++) {
         piv(i) = i;
      }
      pivsign = 1;
      Vector LUcolj(m);

      // Outer loop.

      for (int j = 0; j < n; j++) {

         // Make a copy of the j-th column to localize references.

         for (int i = 0; i < m; i++) {
            LUcolj(i) = LU(i,j);
         }

         // Apply previous transformations.

         for (int i = 0; i < m; i++) {
             matrix_row<Matrix> LUrowi(LU,i);

            // Most of the time is spent in the following dot product.

            int kmax = std::min(i,j);
            double s = 0.0;
            for (int k = 0; k < kmax; k++) {
               s += LUrowi(k)*LUcolj(k);
            }

            LUrowi(j) = LUcolj(i) -= s;
         }
   
         // Find pivot and exchange if necessary.

         int p = j;
         for (int i = j+1; i < m; i++) {
            if (std::abs(LUcolj(i)) > std::abs(LUcolj(p))) {
               p = i;
            }
         }

         if (p != j) {
            for (int k = 0; k < n; k++) {
               double t = LU(p,k); LU(p,k) = LU(j,k); LU(j,k) = t;
            }
            size_t k = piv(p);  // TODO: SIM: fixed 32 bit sloppyness: this was an int...
			piv(p) = piv(j); 
			piv(j) = k;
            pivsign = -pivsign;
         }
         // Compute multipliers.
        
         if (j < m && LU(j,j) != 0.0) {
            for (int i = j+1; i < m; i++) {
               LU(i,j) /= LU(j,j);
            }
         }
      }
   }

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

   public LUDecomposition (Matrix A, int linpackflag) {
      // Initialize.
      LU = A.getArrayCopy();
      m = A.size1();
      n = A.size2();
      piv = new int(m);
      for (int i = 0; i < m; i++) {
         piv(i) = i;
      }
      pivsign = 1;
      // Main loop.
      for (int k = 0; k < n; k++) {
         // Find pivot.
         int p = k;
         for (int i = k+1; i < m; i++) {
            if (std::abs(LU(i,k)) > std::abs(LU(p,k))) {
               p = i;
            }
         }
         // Exchange if necessary.
         if (p != k) {
            for (int j = 0; j < n; j++) {
               double t = LU(p,j); LU(p,j) = LU(k,j); LU(k,j) = t;
            }
            int t = piv(p); piv(p) = piv(k); piv(k) = t;
            pivsign = -pivsign;
         }
         // Compute multipliers and eliminate k-th column.
         if (LU(k,k) != 0.0) {
            for (int i = k+1; i < m; i++) {
               LU(i,k) /= LU(k,k);
               for (int j = k+1; j < n; j++) {
                  LU(i,j) -= LU(i,k)*LU(k,j);
               }
            }
         }
      }
   }

\* ------------------------
   End of temporary code.
 * ------------------------ */

/* ------------------------
   Public Methods
 * ------------------------ */

   /** Return lower triangular factor
   @return     L
   */

LUDecomposition::Matrix LUDecomposition::getL () const {
      int d = std::min(m,n);
      Matrix L(m,d);
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < d; j++) {
            if (i > j) {
               L(i,j) = LU(i,j);
            } else if (i == j) {
               L(i,j) = 1.0;
            } else {
               L(i,j) = 0.0;
            }
         }
      }
      return L;
   }

   /** Return upper triangular factor
   @return     U
   */

LUDecomposition::Matrix LUDecomposition::getU () const {
      int d = std::min(m,n);
      Matrix U(d,n);
      for (int i = 0; i < d; i++) {
         for (int j = 0; j < n; j++) {
            if (i <= j) {
               U(i,j) = LU(i,j);
            } else {
               U(i,j) = 0.0;
            }
         }
      }
      return U;
   }


   /** Solve A*X = B
   @param  B   A Matrix with as many rows as A and any number of columns.
   @return     X so that L*U*X = B(piv,:)
   @exception  bad_size Matrix row dimensions must agree.
   @exception  singular  Matrix is singular.
   */

LUDecomposition::Matrix LUDecomposition::solve (const Matrix& B) const {
      BOOST_UBLAS_CHECK((int)B.size1() == m, bad_size("Matrix row dimensions must agree."));
      BOOST_UBLAS_CHECK(isNonsingular(), singular("Matrix is singular."));

      // Copy right hand side with pivoting
	  assert( B.size2() < size_t(std::numeric_limits<int>::max()));
      int nx = int(B.size2()); // TODO: cast to int
      Matrix X(m,nx);
      for (int i = 0; i < m; i++) {
          row(X,i) = row(B, piv(i));
      }

      // Solve L*Y = B(piv,:)
      for (int k = 0; k < n; k++) {
         for (int i = k+1; i < n; i++) {
            for (int j = 0; j < nx; j++) {
               X(i,j) -= X(k,j)*LU(i,k);
            }
         }
      }
      // Solve U*X = Y;
      for (int k = n-1; k >= 0; k--) {
         for (int j = 0; j < nx; j++) {
            X(k,j) /= LU(k,k);
         }
         for (int i = 0; i < k; i++) {
            for (int j = 0; j < nx; j++) {
               X(i,j) -= X(k,j)*LU(i,k);
            }
         }
      }
      return X;
   }

}}}
