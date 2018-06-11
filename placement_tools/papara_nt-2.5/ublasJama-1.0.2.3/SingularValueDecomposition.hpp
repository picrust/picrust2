   /** Singular Value Decomposition.
   <P>
   For an m-by-n matrix A, the singular value decomposition is
   an m-by-(m or n) orthogonal matrix U, an (m or n)-by-n diagonal matrix S, and
   an n-by-n orthogonal matrix V so that A = U*S*V'.
   <P>
   The singular values, sigma[k] = S[k][k], are ordered so that
   sigma[0] >= sigma[1] >= ... >= sigma[n-1].
   <P>
   The singular value decompostion always exists, so the constructor will
   never fail.  The matrix condition number and the effective numerical
   rank can be computed from this decomposition.
   */

// This version includes modifications and fixes by Andreas Kyrmegalos
// explanation: http://cio.nist.gov/esd/emaildir/lists/jama/msg01430.html
// final version: http://cio.nist.gov/esd/emaildir/lists/jama/msg01431.html

#ifndef _BOOST_UBLAS_SINGULARVALUEDECOMPOSITION_
#define _BOOST_UBLAS_SINGULARVALUEDECOMPOSITION_

#include <algorithm>
#include <cmath>
#include <limits>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/exception.hpp>
#include <boost/math/special_functions/hypot.hpp>

namespace boost { namespace numeric { namespace ublas {
            
template<class T>
class SingularValueDecomposition {

    typedef vector<T> vector_type;
    typedef matrix<T> matrix_type;

/* ------------------------
   Class variables
 * ------------------------ */

   /** Arrays for internal storage of U and V.
   @serial internal storage of U.
   @serial internal storage of V.
   */
   matrix_type U, V;

   /** Array for internal storage of singular values.
   @serial internal storage of singular values.
   */
   vector_type s;

   /** Row and column dimensions.
   @serial row dimension.
   @serial column dimension.
   @serial U column dimension.
   */
   int m, n, ncu;

   /** Column specification of matrix U
   @serial U column dimension toggle
   */
   
   bool thin;
   
   /** Construct the singular value decomposition
   @param A    Rectangular matrix
   @param thin  If true U is economy sized
   @param wantu If true generate the U matrix
   @param wantv If true generate the V matrix
   @return     Structure to access U, S and V.
   */
   void init (const matrix_type &Arg, bool thin, bool wantu, bool wantv);

   static matrix_vector_slice<matrix_type> subcolumn(matrix_type& M,size_t c,size_t start,size_t stop) {
      return matrix_vector_slice<matrix_type> (M, slice(start,1,stop-start), slice(c,0,stop-start));
   }
   static matrix_vector_slice<matrix_type> subrow(matrix_type& M,size_t r,size_t start,size_t stop) {
      return matrix_vector_slice<matrix_type> (M, slice(r,0,stop-start), slice(start,1,stop-start));
   }

public:
/* ------------------------
   Old Constructor
 * ------------------------ */
   /** Construct the singular value decomposition
   @param Arg  Rectangular matrix
   @return     Structure to access U, S and V.
   */
   
   SingularValueDecomposition (const matrix_type &Arg) {
      init(Arg,true,true,true);
   }
   
/* ------------------------
   Constructor
 * ------------------------ */

   /** Construct the singular value decomposition
   @param A    Rectangular matrix
   @param thin  If true U is economy sized
   @param wantu If true generate the U matrix
   @param wantv If true generate the V matrix
   @return     Structure to access U, S and V.
   */

   SingularValueDecomposition (const matrix_type &Arg, bool thin, bool wantu, bool wantv) {
      init(Arg,thin,wantu,wantv);
   }
    
/* ------------------------
   Public Methods
 * ------------------------ */

   /** Return the left singular vectors
   @return     U
   */

   const matrix_type& getU () const {
      return U;
   }

   /** Return the right singular vectors
   @return     V
   */

   const matrix_type& getV () const {
      return V;
   }

   /** Return the one-dimensional array of singular values
   @return     diagonal of S.
   */

   const vector_type& getSingularValues () const {
      return s;
   }

   /** Return the diagonal matrix of singular values
   @return     S
   */

   matrix_type getS () const {
      matrix_type S(m>=n?(thin?n:ncu):ncu,n,T/*zero*/());
      for (int i = std::min(m,n)-1; i >= 0; i--) {
         S(i,i) = s(i);
      }
      return S;
   }

   /** Return the diagonal matrix of the reciprocals of the singular values
   @return     S+
   */
   
   matrix_type getreciprocalS () const {
      matrix_type S(n,m>=n?(thin?n:ncu):ncu,T/*zero*/());
      for (int i = std::min(m,n)-1; i>=0; i--)
         S(i,i) = s(i)==T/*zero*/()?0.0:1.0/s(i);
      return S;
   }
   
   /** Return the rightmost column of V.
       Does not check to see whether or not the matrix actually was rank-deficient -
       the caller is assumed to have examined S and decided that to his or her satisfaction.
   @return     null vector
   */

   matrix_column<const matrix_type> getNullVector () const {
       return matrix_column<const matrix_type>(V, n-1);
   }

   /** Return the Moore-Penrose (generalized) inverse
    *  Slightly modified version of Kim van der Linde's code
   @param omit if true tolerance based omitting of negligible singular values
   @return     A+
   */
   
   matrix_type inverse(bool omit = true) const {
      matrix_type inverse(n,m);
      if(rank()> 0) {
         vector_type reciprocalS(s.size());
         if (omit) {
            T tol = std::max(m,n)*s(0)*std::numeric_limits<T>::epsilon();
            for (int i = s.size()-1;i>=0;i--)
               reciprocalS(i) = std::abs(s(i))<tol?0.0:1.0/s(i);
         }
         else
            for (int i=s.size()-1;i>=0;i--)
               reciprocalS(i) = s(i)==T/*zero*/()?0.0:1.0/s(i);
         int min = std::min(n, ncu);
         for (int i = n-1; i >= 0; i--)
            for (int j = m-1; j >= 0; j--)
               for (int k = min-1; k >= 0; k--)
                  inverse(i,j) += V(i,k) * reciprocalS(k) * U(j,k);
      } 
      return inverse;
   }

   /** Two norm
   @return     max(S)
   */

   T norm2 () const {
      return s(0);
   }

   /** Two norm condition number
   @return     max(S)/min(S)
   */

   T cond () const {
      return s(0)/s(std::min(m,n)-1);
   }

   /** Effective numerical matrix rank
   @return     Number of nonnegligible singular values.
   */

   inline unsigned int rank () const {
      T tol = std::max(m,n)*s[0]*std::numeric_limits<T>::epsilon();
      int r = 0;
      for (unsigned i = 0; i < s.size(); i++) {
         if (s(i) > tol) {
            r++;
         }
      }
      return r;
   }

};

template<class T>
void SingularValueDecomposition<T>::init (const matrix_type &Arg, bool thin, bool wantu, bool wantv) {

   // Derived from LINPACK code.
   // Initialize.
   matrix_type A = Arg;
   m = Arg.size1();
   n = Arg.size2();
   this->thin = thin;

   ncu = thin?std::min(m,n):m;
   s = vector_type(std::min(m+1,n));
   if (wantu) {
      U = matrix_type(m,ncu,T/*zero*/());
   }
   if (wantv) {
      V = matrix_type(n,n,T/*zero*/());
   }
   vector_type e(n);
   vector_type work(m);

   // Reduce A to bidiagonal form, storing the diagonal elements
   // in s and the super-diagonal elements in e.

   int nct = std::min(m-1,n);
   int nrt = std::max(0,std::min(n-2,m));
   int lu = std::max(nct,nrt);
   for (int k = 0; k < lu; k++) {
      if (k < nct) {

         // Compute the transformation for the k-th column and
         // place the k-th diagonal in s[k].
         // Compute 2-norm of k-th column without under/overflow.
         // s(k) = norm of elements k..m-1 of column k of A
         s(k) = norm_2(subcolumn(A,k,k,m));
         if (s(k) != T/*zero*/()) {
            if (A(k,k) < T/*zero*/()) {
               s(k) = -s(k);
            }
            // divide elements k..m-1 of column k of A by s(k)
            subcolumn(A,k,k,m) /= s(k);
            A(k,k) += 1.0;
         }
         s(k) = -s(k);
      }
      for (int j = k+1; j < n; j++) {
         if ((k < nct) && (s(k) != T/*zero*/()))  {

            // Apply the transformation.

            // t = dot-product of elements k..m-1 of columns k and j of A
            T t = inner_prod(subcolumn(A,k,k,m),subcolumn(A,j,k,m));
            t = -t/A(k,k);
            // elements k..m-1 of column j of A +=  t*(elements k..m-1 of column k of A)
            subcolumn(A,j,k,m) += t*subcolumn(A,k,k,m);
         }

         // Place the k-th row of A into e for the
         // subsequent calculation of the row transformation.

         e(j) = A(k,j);
      }
      if (wantu && (k < nct)) {

         // Place the transformation in U for subsequent back
         // multiplication.

         // elements k..m-1 of column k of U = elements k..m-1 of column k of A
         subcolumn(U,k,k,m) = subcolumn(A,k,k,m);
      }
      if (k < nrt) {

         // Compute the k-th row transformation and place the
         // k-th super-diagonal in e[k].
         // Compute 2-norm without under/overflow.
         // e(k) = norm of elements k+1..n-1 of e
         e(k) = norm_2(subrange(e,k+1,n));
         if (e(k) != T/*zero*/()) {
            if (e(k+1) < T/*zero*/()) {
               e(k) = -e(k);
            }
            // divide elements k+1..n-1 of e by e(k)
            subrange(e,k+1,n) /= e(k);
            e(k+1) += 1.0;
         }
         e(k) = -e(k);
         if ((k+1 < m) && (e(k) != T/*zero*/())) {

            // Apply the transformation.

            for (int i = k+1; i < m; i++) {
               work(i) = 0.0;
            }
            // elements k+1..m-1 of work = A.submatrix(k+1..m-1,k+1..n-1)*elements k+1..n-1 of e
            noalias(subrange(work, k+1, m)) = prod(subrange(A,k+1,m,k+1,n),subrange(e,k+1,n));
            for (int j = k+1; j < n; j++) {
               T t = -e(j)/e(k+1);
               // elements k+1..m-1 of column j of A += t*elements k+1..m-1 of work
               subcolumn(A,j,k+1,m) += t*subrange(work,k+1,m);
            }
         }
         if (wantv) {

            // Place the transformation in V for subsequent
            // back multiplication.

            // elements k+1..n-1 of column k of V = elements k+1..n-1 of e
            subcolumn(V,k,k+1,n) = subrange(e,k+1,n);
         }
      }
   }

   // Set up the final bidiagonal matrix or order p.

   int p = std::min(n,m+1);
   if (nct < n) {
      s(nct) = A(nct,nct);
   }
   if (m < p) {
      s(p-1) = 0.0;
   }
   if (nrt+1 < p) {
      e(nrt) = A(nrt,p-1);
   }
   e(p-1) = 0.0;

   // If required, generate U.

   if (wantu) {
      for (int j = nct; j < ncu; j++) {
         // set column j of U to zero
         std::fill(column(U,j).begin(),column(U,j).end(),T/*zero*/());
         U(j,j) = 1.0;
      }
      for (int k = nct-1; k >= 0; k--) {
         if (s(k) != T/*zero*/()) {
            for (int j = k+1; j < ncu; j++) {
               // t = dot-product of elements k..m-1 of columns k and j of U
               T t = inner_prod(subcolumn(U,k,k,m),subcolumn(U,j,k,m));
               t /= -U(k,k);
               // elements k..m-1 of column j of U +=  t*(elements k..m-1 of column k of U)
               subcolumn(U,j,k,m) += t*subcolumn(U,k,k,m);
            }
            // elements k..m-1 of column k of U *= -1.
            subcolumn(U,k,k,m) *= -1.0;
            U(k,k) += 1.0;
            if(k-1 > 0) {
               // set elements 0..k-2 of column k of U to zero.
               for (int i = 0; i < k-1; i++) { 
                  U(i,k) = T/*zero*/();
               }
            }
         } else {
            // set column k of U to zero
            std::fill(column(U,k).begin(),column(U,k).end(),T/*zero*/());
            U(k,k) = 1.0;
         }
      }
   }

   // If required, generate V.

   if (wantv) {
      for (int k = n-1; k >= 0; k--) {
         if ((k < nrt) && (e(k) != T/*zero*/())) {
            for (int j = k+1; j < n; j++) {
               // t = dot-product of elements k+1..n-1 of columns k and j of V
               T t = inner_prod(subcolumn(V,k,k+1,n),subcolumn(V,j,k+1,n));
               t /= -V(k+1,k);
               // elements k+1..n-1 of column j of V +=  t*(elements k+1..n-1 of column k of V)
               subcolumn(V,j,k+1,n) += t*subcolumn(V,k,k+1,n);
            }
         }
         // set column k of V to zero
         std::fill(column(V,k).begin(),column(V,k).end(),T/*zero*/());
         V(k,k) = 1.0;
      }
   }

   // Main iteration loop for the singular values.

   int pp = p-1;
   int iter = 0;
   T eps = std::numeric_limits<T>::epsilon();
   T tiny = std::numeric_limits<T>::min();
   while (p > 0) {
      int k,kase;

      // Here is where a test for too many iterations would go.

      // This section of the program inspects for
      // negligible elements in the s and e arrays.  On
      // completion the variables kase and k are set as follows.

      // kase = 1     if s(p) and e[k-1] are negligible and k<p
      // kase = 2     if s(k) is negligible and k<p
      // kase = 3     if e[k-1] is negligible, k<p, and
      //              s(k), ..., s(p) are not negligible (qr step).
      // kase = 4     if e(p-1) is negligible (convergence).

      for (k = p-2; k >= -1; k--) {
         if (k == -1) {
            break;
         }
         if (std::abs(e(k)) <=
             tiny + eps*(std::abs(s(k)) + std::abs(s(k+1)))) {
            e(k) = 0.0;
            break;
         }
      }
      if (k == p-2) {
         kase = 4;
      } else {
         int ks;
         for (ks = p-1; ks >= k; ks--) {
            if (ks == k) {
               break;
            }
            T t = (ks != p ? std::abs(e(ks)) : 0.) + 
                  (ks != k+1 ? std::abs(e(ks-1)) : 0.);
            if (std::abs(s(ks)) <= tiny + eps*t)  {
               s[ks] = 0.0;
               break;
            }
         }
         if (ks == k) {
            kase = 3;
         } else if (ks == p-1) {
            kase = 1;
         } else {
            kase = 2;
            k = ks;
         }
      }
      k++;

      // Perform the task indicated by kase.

      switch (kase) {

         // Deflate negligible s(p).

         case 1: {
            T f = e(p-2);
            e[p-2] = 0.0;
            for (int j = p-2; j >= k; j--) {
               T t = boost::math::hypot(s(j),f);
               T cs = s(j)/t;
               T sn = f/t;
               s(j) = t;
               if (j != k) {
                  f = -sn*e(j-1);
                  e(j-1) = cs*e(j-1);
               }
               if (wantv) {
                  for (int i = 0; i < n; i++) {
                     t = cs*V(i,j) + sn*V(i,p-1);
                     V(i,p-1) = -sn*V(i,j) + cs*V(i,p-1);
                     V(i,j) = t;
                  }
               }
            }
         }
            break;

            // Split at negligible s(k).

         case 2: {
            T f = e(k-1);
            e[k-1] = 0.0;
            for (int j = k; j < p; j++) {
               T t = boost::math::hypot(s(j),f);
               T cs = s(j)/t;
               T sn = f/t;
               s(j) = t;
               f = -sn*e(j);
               e(j) = cs*e(j);
               if (wantu) {
                  for (int i = 0; i < m; i++) {
                     t = cs*U(i,j) + sn*U(i,k-1);
                     U(i,k-1) = -sn*U(i,j) + cs*U(i,k-1);
                     U(i,j) = t;
                  }
               }
            }
         }
            break;

            // Perform one qr step.

         case 3: {

            // Calculate the shift.
   
            T scale = std::max(std::max(std::max(std::max(
                std::abs(s(p-1)),std::abs(s(p-2))),std::abs(e(p-2))), 
                                             std::abs(s(k))),std::abs(e(k)));
            T sp = s(p-1)/scale;
            T spm1 = s(p-2)/scale;
            T epm1 = e(p-2)/scale;
            T sk = s(k)/scale;
            T ek = e(k)/scale;
            T b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
            T c = (sp*epm1)*(sp*epm1);
            T shift = 0.0;
            if ((b != T/*zero*/()) || (c != T/*zero*/())) {
               shift = std::sqrt(b*b + c);
               if (b < T/*zero*/()) {
                  shift = -shift;
               }
               shift = c/(b + shift);
            }
            T f = (sk + sp)*(sk - sp) + shift;
            T g = sk*ek;
   
            // Chase zeros.
   
            for (int j = k; j < p-1; j++) {
               T t = boost::math::hypot(f,g);
               T cs = f/t;
               T sn = g/t;
               if (j != k) {
                  e(j-1) = t;
               }
               f = cs*s(j) + sn*e(j);
               e(j) = cs*e(j) - sn*s(j);
               g = sn*s(j+1);
               s(j+1) = cs*s(j+1);
               if (wantv) {
                  for (int i = 0; i < n; i++) {
                     t = cs*V(i,j) + sn*V(i,j+1);
                     V(i,j+1) = -sn*V(i,j) + cs*V(i,j+1);
                     V(i,j) = t;
                  }
               }
               t = boost::math::hypot(f,g);
               cs = f/t;
               sn = g/t;
               s(j) = t;
               f = cs*e(j) + sn*s(j+1);
               s(j+1) = -sn*e(j) + cs*s(j+1);
               g = sn*e(j+1);
               e(j+1) = cs*e(j+1);
               if (wantu && (j < m-1)) {
                  for (int i = 0; i < m; i++) {
                     t = cs*U(i,j) + sn*U(i,j+1);
                     U(i,j+1) = -sn*U(i,j) + cs*U(i,j+1);
                     U(i,j) = t;
                  }
               }
            }
            e(p-2) = f;
            iter++;
         }
            break;

            // Convergence.

         case 4: {
            // Make the singular values positive.
   
            if (s(k) <= T/*zero*/()) {
               s(k) = (s(k) < T/*zero*/() ? -s(k) : 0.0);
               if (wantv) {
                  // multiply column k of V by -1
                  column(V,k) *= -1.0;
               }
            }
   
            // Order the singular values.
   
            while (k < pp) {
               if (s(k) >= s(k+1)) {
                  break;
               }
               T t = s(k);
               s(k) = s(k+1);
               s(k+1) = t;
               if (wantv && (k < n-1)) {
                  // swap columns k and k+1 of V
                  column(V,k).swap(column(V,k+1));
               }
               if (wantu && (k < m-1)) {
                  // swap columns k and k+1 of U
                  column(U,k).swap(column(U,k+1));
               }
               k++;
            }
            iter = 0;
            p--;
         }
            break;
      }
   }
}

}}}

#endif
