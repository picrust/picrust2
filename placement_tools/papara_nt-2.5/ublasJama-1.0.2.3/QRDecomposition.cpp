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

#include <algorithm>
#include <cmath>
#include <boost/math/special_functions/hypot.hpp>
#include "QRDecomposition.hpp"

namespace boost { namespace numeric { namespace ublas {

   /** sqrt(a^2 + b^2) without under/overflow. **/
   /*
   static inline double hypot(double a, double b) {
      double r;
      if (std::abs(a) > std::abs(b)) {
         r = b/a;
         r = std::abs(a)*std::sqrt(1+r*r);
      } else if (b != 0) {
         r = a/b;
         r = std::abs(b)*std::sqrt(1+r*r);
      } else {
         r = 0.0;
      }
      return r;
   }
   */
    
    
QRDecomposition::QRDecomposition (const Matrix &A) {
      // Initialize.
      QR = A;
	  assert( A.size1() < size_t(std::numeric_limits<int>::max()) && A.size2() < size_t(std::numeric_limits<int>::max()));
      m = int(A.size1()); // TODO: cast to int
      n = int(A.size2());
      Rdiag = Vector(n);

      // Main loop.
      for (int k = 0; k < n; k++) {
         // Compute 2-norm of k-th column without under/overflow.
         double nrm = 0;
         for (int i = k; i < m; i++) {
             nrm = boost::math::hypot(nrm,QR(i,k));
         }

         if (nrm != 0.0) {
            // Form k-th Householder vector.
            if (QR(k,k) < 0) {
               nrm = -nrm;
            }
            for (int i = k; i < m; i++) {
               QR(i,k) /= nrm;
            }
            QR(k,k) += 1.0;

            // Apply transformation to remaining columns.
            for (int j = k+1; j < n; j++) {
               double s = 0.0; 
               for (int i = k; i < m; i++) {
                  s += QR(i,k)*QR(i,j);
               }
               s = -s/QR(k,k);
               for (int i = k; i < m; i++) {
                  QR(i,j) += s*QR(i,k);
               }
            }
         }
         Rdiag[k] = -nrm;
      }
   }
   /** Is the matrix full rank?
   @return     true if R, and hence A, has full rank.
   */

bool QRDecomposition::isFullRank () const {
      for (int j = 0; j < n; j++) {
         if (Rdiag(j) == 0.)
            return false;
      }
      return true;
   }

   /** Return the Householder vectors
   @return     Lower trapezoidal matrix whose columns define the reflections
   */

QRDecomposition::Matrix QRDecomposition::getH () const {
      Matrix H(m,n);
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            if (i >= j) {
               H(i,j) = QR(i,j);
            } else {
               H(i,j) = 0.0;
            }
         }
      }
      return H;
   }

   /** Return the upper triangular factor
   @return     R
   */

QRDecomposition::Matrix QRDecomposition::getR () const {
      Matrix R(n,n);
      for (int i = 0; i < n; i++) {
         for (int j = 0; j < n; j++) {
            if (i < j) {
               R(i,j) = QR(i,j);
            } else if (i == j) {
               R(i,j) = Rdiag(i);
            } else {
               R(i,j) = 0.0;
            }
         }
      }
      return R;
   }

   /** Generate and return the (economy-sized) orthogonal factor
   @return     Q
   */

QRDecomposition::Matrix QRDecomposition::getQ () const {
      Matrix Q(m,n);
      for (int k = n-1; k >= 0; k--) {
         for (int i = 0; i < m; i++) {
            Q(i,k) = 0.0;
         }
         Q(k,k) = 1.0;
         for (int j = k; j < n; j++) {
            if (QR(k,k) != 0) {
               double s = 0.0;
               for (int i = k; i < m; i++) {
                  s += QR(i,k)*Q(i,j);
               }
               s = -s/QR(k,k);
               for (int i = k; i < m; i++) {
                  Q(i,j) += s*QR(i,k);
               }
            }
         }
      }
      return Q;
   }
            
   /** Least squares solution of A*X = B
   @param B    A Matrix with as many rows as A and any number of columns.
   @return     X that minimizes the two norm of Q*R*X-B.
   @exception  IllegalArgumentException  Matrix row dimensions must agree.
   @exception  RuntimeException  Matrix is rank deficient.
   */

QRDecomposition::Matrix QRDecomposition::solve (const Matrix &B) {
      BOOST_UBLAS_CHECK ((int)B.size1() == m, bad_size ());
      BOOST_UBLAS_CHECK (isFullRank(), singular ());
      
      // Copy right hand side
	   assert( B.size2() < size_t(std::numeric_limits<int>::max()));

      int nx = int(B.size2()); // TODO: cast to int
      Matrix X(B);

      // Compute Y = transpose(Q)*B
      for (int k = 0; k < n; k++) {
         for (int j = 0; j < nx; j++) {
            double s = 0.0; 
            for (int i = k; i < m; i++) {
               s += QR(i,k)*X(i,j);
            }
            s = -s/QR(k,k);
            for (int i = k; i < m; i++) {
               X(i,j) += s*QR(i,k);
            }
         }
      }
      // Solve R*X = Y;
      for (int k = n-1; k >= 0; k--) {
         for (int j = 0; j < nx; j++) {
            X(k,j) /= Rdiag(k);
         }
         for (int i = 0; i < k; i++) {
            for (int j = 0; j < nx; j++) {
               X(i,j) -= X(k,j)*QR(i,k);
            }
         }
      }
      Matrix subX = subrange(X,0,n,0,nx);
      return subX;
   }

}}}
