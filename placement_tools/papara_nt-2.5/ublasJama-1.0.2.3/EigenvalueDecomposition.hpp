/** Eigenvalues and eigenvectors of a real matrix. 
<P>
    If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is
    diagonal and the eigenvector matrix V is orthogonal.
    I.e. A = V.times(D.times(V.transpose())) and 
    V.times(V.transpose()) equals the identity matrix.
<P>
    If A is not symmetric, then the eigenvalue matrix D is block diagonal
    with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
    lambda + i*mu, in 2-by-2 blocks, [lambda, mu; -mu, lambda].  The
    columns of V represent the eigenvectors in the sense that A*V = V*D,
    i.e. A.times(V) equals V.times(D).  The matrix V may be badly
    conditioned, or even singular, so the validity of the equation
    A = V*D*inverse(V) depends upon V.cond().
**/

#ifndef _BOOST_UBLAS_EIGENVALUEDECOMPOSITION_
#define _BOOST_UBLAS_EIGENVALUEDECOMPOSITION_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/exception.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace boost { namespace numeric { namespace ublas {
    
class EigenvalueDecomposition {

    typedef vector<double> Vector;
    typedef matrix<double> Matrix;

/* ------------------------
   Class variables
 * ------------------------ */

   /** Row and column dimension (square matrix).
   @serial matrix dimension.
   */
   int n;

   /** Symmetry flag.
   @serial internal symmetry flag.
   */
   bool issymmetric;

   /** Arrays for internal storage of eigenvalues.
   @serial internal storage of eigenvalues.
   */
   Vector d, e;

   /** Array for internal storage of eigenvectors.
   @serial internal storage of eigenvectors.
   */
   Matrix V;

   /** Array for internal storage of nonsymmetric Hessenberg form.
   @serial internal storage of nonsymmetric Hessenberg form.
   */
   Matrix H;

   /** Working storage for nonsymmetric algorithm.
   @serial working storage for nonsymmetric algorithm.
   */
   Vector ort;

/* ------------------------
   Private Methods
 * ------------------------ */

   // Symmetric Householder reduction to tridiagonal form.

   void tred2 ();

   // Symmetric tridiagonal QL algorithm.
   
   void tql2 ();

   // Nonsymmetric reduction to Hessenberg form.

   void orthes ();

   // Nonsymmetric reduction from Hessenberg to real Schur form.

   void hqr2 ();

public:
/* ------------------------
   Constructor
 * ------------------------ */

   /** Check for symmetry, then construct the eigenvalue decomposition
   @param A    Square matrix
   @return     Structure to access D and V.
   */

   EigenvalueDecomposition (const Matrix& A);

/* ------------------------
   Public Methods
 * ------------------------ */

   /** Is the matrix symmetric?
   @return     true if A is symmetric.
   */
   bool isSymmetric() const {
       return issymmetric;
   }
    
   /** Return the eigenvector matrix
   @return     V
   */

   const Matrix& getV () const {
      return V;
   }

   /** Return the real parts of the eigenvalues
   @return     real(diag(D))
   */

   const Vector& getRealEigenvalues () const {
      return d;
   }

   /** Return the imaginary parts of the eigenvalues
   @return     imag(diag(D))
   */

   const Vector& getImagEigenvalues () const {
      return e;
   }

   /** Return the block diagonal eigenvalue matrix
   @return     D
   */

   Matrix getD () const;
};

}}}
#endif
