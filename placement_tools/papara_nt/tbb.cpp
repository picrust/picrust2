/*
 * Copyright (C) 2009-2012 Simon A. Berger
 * 
 * This file is part of papara.
 * 
 *  papara is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  papara is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with papara.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <iostream>
#include <iterator>
// #include <tbb/task.h>

#include "ivymike/time.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <math.h>
#include "EigenvalueDecomposition.hpp"
#include "LUDecomposition.hpp"

// class fib_task : public tbb::task {
//     const long m_n;
//     long *const m_sum;
//     
// public:
//     fib_task( long n, long *sum ) : m_n(n), m_sum(sum) {}
// //     virtual ~fib_task() {
// //         std::cout << "term\n";
// //     }
//     
//     task * execute() {
// //         std::cout << "execute\n";
//         if( m_n < 2 ) {
//             *m_sum = m_n;
//         } else {
//             long x, y;
//             
//             fib_task &a = *new(allocate_child() ) fib_task(m_n-1, &x);
//             fib_task &b = *new(allocate_child() ) fib_task(m_n-2, &y);
//             
//             set_ref_count(3);
//             spawn(b);
//             spawn_and_wait_for_all(a);
//             
//             *m_sum = x + y;
//             
//         }
//         return 0;
//     }
//     
// };

// namespace {
// using namespace boost::numeric::ublas;
// 
// static void mytred2( matrix<double> &a, const int n, vector<double> &d, vector<double> &e)
// {
//     
//     
//     for ( int i = n; i > 1; i--)
//     {
//         int l = i - 1;
//         double h = 0.0;
//         double scale = 0.0;
// 
//         if (l > 1)
//         {
//             for (int k = 1; k <= l; k++)
//                 scale += fabs(a(k - 1, i - 1));
//             if (scale == 0.0)
//                 e[i - 1] = a(l - 1, i - 1);
//             else
//             {
//                 for ( int k = 1; k <= l; k++)
//                 {
//                     a(k - 1, i - 1) /= scale;
//                     h += a(k - 1, i - 1) * a(k - 1, i - 1);
//                 }
//                 double f = a(l - 1,i - 1);
//                 double g = ((f > 0) ? -sqrt(h) : sqrt(h)); /* diff */
//                 e[i - 1] = scale * g;
//                 h -= f * g;
//                 a(l - 1, i - 1) = f - g;
//                 f = 0.0;
//                 for ( int j = 1; j <= l; j++)
//                 {
//                     a(i - 1,j - 1) = a(j - 1,i - 1) / h;
//                     g = 0.0;
//                     for ( int k = 1; k <= j; k++)
//                         g += a(k - 1, j - 1) * a(k - 1, i - 1);
//                     for ( int k = j + 1; k <= l; k++)
//                         g += a(j - 1,k - 1) * a(k - 1,i - 1);
//                     e[j - 1] = g / h;
//                     f += e[j - 1] * a(j - 1, i - 1);
//                 }
//                 double hh = f / (h + h);
//                 for ( int j = 1; j <= l; j++)
//                 {
//                     f = a(j - 1, i - 1);
//                     g = e[j - 1] - hh * f;
//                     e[j - 1] = g;
//                     for ( int k = 1; k <= j; k++)
//                         a(k - 1, j - 1) -= (f * e[k - 1] + g * a(k - 1, i - 1));
//                 }
//             }
//         }
//         else
//             e[i - 1] = a(l - 1, i - 1);
//         d[i - 1] = h;
//     }
//     d[0] = 0.0;
//     e[0] = 0.0;
// 
//     for ( int i = 1; i <= n; i++)
//     {
//         int l = i - 1;
//         if (d[i - 1] != 0.0)
//         {
//             for ( int j = 1; j <= l; j++)
//             {
//                 double g = 0.0;
//                 for ( int k = 1; k <= l; k++)
//                     g += a(k - 1, i - 1) * a(j - 1, k - 1);
//                 for ( int k = 1; k <= l; k++)
//                     a(j - 1, k - 1) -= g * a(i - 1, k - 1);
//             }
//         }
//         d[i - 1] = a(i - 1, i - 1);
//        
//         
//         a(i - 1, i - 1) = 1.0;
//         for ( int j = 1; j <= l; j++)
//             a(i - 1, j - 1) = a(j - 1, i - 1) = 0.0;
//     }
// 
// 
// }
// 
// 
// // this seems to be the function tqli from NR chaper 11, but with rows/columns of z switched,
// // so that the output eigenvectors are in the rows of z
// static int mytqli( vector<double> &d, vector<double> &e, const int n, matrix<double> &z)
// {
//     int     m, l, iter, i, k;
//     double  s, r, p, g, f, dd, c, b;
// 
//     for (i = 2; i <= n; i++)
//         e[i - 2] = e[i - 1];
// 
//     e[n - 1] = 0.0;
// 
//     for (l = 1; l <= n; l++)
//     {
//         iter = 0;
//         do
//         {
//             for (m = l; m <= n - 1; m++)
//             {
//                 
//                 dd = fabs(d[m - 1]) + fabs(d[m]);
//                 //std::cout << "s: " << d[m-1] << "\n";
//                 
//                 if (fabs(e[m - 1]) + dd == dd)
//                     break;
//             }
// 
//             if (m != l)
//             {
//                 assert(iter < 30);
// 
//                 g = (d[l] - d[l - 1]) / (2.0 * e[l - 1]);
// 
//                 r = sqrt((g * g) + 1.0);
//                 g = d[m - 1] - d[l - 1] + e[l - 1] / (g + ((g < 0)?-fabs(r):fabs(r)));/*MYSIGN(r, g));*/
//                 s = c = 1.0;
//                 p = 0.0;
// 
//                 for (i = m - 1; i >= l; i--)
//                 {
//                     f = s * e[i - 1];
//                     b = c * e[i - 1];
//                     
//                     if (fabs(f) >= fabs(g))
//                     {
//                         c = g / f;
//                         r = sqrt((c * c) + 1.0);
//                         e[i] = f * r;
//                         c *= (s = 1.0 / r);
//                     }
//                     else
//                     {
//                         s = f / g;
//                         
// 
//                         r = sqrt((s * s) + 1.0);
//                         e[i] = g * r;
//                         s *= (c = 1.0 / r);
//                     }
//                     g = d[i] - p;
//                     r = (d[i - 1] - g) * s + 2.0 * c * b;
//                     p = s * r;
//                     d[i] = g + p;
//                     g = c * r - b;
//                     for (k = 1; k <= n; k++)
//                     {
//                         f = z(i,k-1);
//                         z(i,k-1) = s * z(i - 1, k - 1) + c * f;
//                         z(i - 1, k - 1) = c * z(i - 1, k - 1) - s * f;
//                     }
//                 }
// 
//                 d[l - 1] = d[l - 1] - p;
//                 e[l - 1] = g;
//                 e[m - 1] = 0.0;
//             }
//         }
//         while (m != l);
//     }
// 
// 
// 
//     return (1);
// }
// 
// 
// static void makeEigen( matrix<double> &_a, const int n, vector<double> &d, vector<double> &e)
// {
//     
//     mytred2(_a, n, d, e);
//     
//     
//     mytqli(d, e, n, _a);
// }
// }

template<typename T>
void bla( T x ) {
    std::cout << "bla1\n";
}

template<typename T>
void bla( T *x ) {
    std::cout << "bla2\n";
}


int main() {
    int x;
    int *y = &x;
    
    bla( x );
    bla( y );
    
//     {
//         long n = 10;
//         long sum = 0;
//         using namespace tbb;
//         fib_task &a = *new(task::allocate_root()) fib_task(n, &sum);
//         task::spawn_root_and_wait(a);
//         std::cout << sum << "\n";
//     }
    
    using namespace boost::numeric::ublas;
    using namespace boost::numeric;
    const size_t n = 2;
    
    vector<double> init_rates(n);
//     init_rates(0) = 0.5;
//     init_rates(1) = 0.5;
//     init_rates(2) = 0.5;
//     init_rates(3) = 1.0;
    
    vector<double> frequencies(n);
    
    
    frequencies(0) = 0.55;
    frequencies(1) = 0.45;
//     frequencies(2) = 0.25;
//     frequencies(3) = 0.25;
    
    //std::cout << frequencies << "\n";
//     throw std::runtime_error( "xxx" );
    
    matrix<double> rate_matrix(n,n,0.0);
    
    
    
    
    for ( size_t i=0; i < n; i++)
    {
        for ( size_t j=i+1;  j < n; j++)
        {
            //             double factor =  initialRates[m++];
//             double factor = 1.0; // initial rate = 1
#if 0
            rate_matrix(i,j) = rate_matrix(j,i) = factor * sqrt( frequencies(i) * frequencies(j));
            rate_matrix(i,i) -= factor * frequencies(j);
            rate_matrix(j,j) -= factor * frequencies(i);
#else
            rate_matrix(i,j) = frequencies(i);
            rate_matrix(j,i) = frequencies(j);
            
            rate_matrix(i,i) -= frequencies(i);
            rate_matrix(j,j) -= frequencies(j);
#endif
        }
    }
    
    ublas::EigenvalueDecomposition ed(rate_matrix);
    
    std::cout << "jama vecs: " <<ed.getV() << "\n";
    std::cout << "jama vals: " <<ed.getRealEigenvalues() << "\n";
    
    ublas::LUDecomposition lud(ed.getV());
    ublas::matrix<double> evecs_pinv = lud.pseudoinverse();
    std::cout << "jama pinv: " << evecs_pinv << "\n";
    
    ublas::vector<double>evals = ed.getRealEigenvalues();
    
    ublas::diagonal_matrix<double> diag_evals(n); 
    
    for( size_t i = 0 ; i < n; i++ ) {
        diag_evals(i,i) = exp( 0.1 * evals(i) );
    }
    
    ublas::matrix<double> t(n,n);
    ivy_mike::timer t1;
    
    t = prod( ed.getV(), diag_evals );
    t = prod( t, evecs_pinv );
    
    std::cout << "t: \n" << t1.elapsed() << "\n";
    
    std::cout << "res: " << t << "\n";
        
    ublas::matrix<double> pvm(2, 4);
    pvm(0,0) = 1;
    pvm(1,0) = 0;
    
    pvm(0,1) = 1;
    pvm(1,1) = 0;
    
    pvm(0,2) = 0;
    pvm(1,2) = 1;
    
    pvm(0,3) = 0;
    pvm(1,3) = 1;
    
    pvm = prod( t, pvm );
    
    std::cout << pvm << "\n";
    
    
//     std::cout << diag_evs(0,0) << " " << diag_evs(0,1) << "\n" << diag_evs(1,0) << " " << diag_evs(1,1) << "\n";
    
//     std::copy( m[0].begin1(), m[0].end1(), std::ostream_iterator<double>( std::cout, "\n" ));
//     std::copy( m[1].begin1(), m[1].end1(), std::ostream_iterator<double>( std::cout, "\n" ));
    
//     double a = m(0,0);
//     double b = m(0,1);
//     double c = m(1,0);
//     double d = m(1,1);
//     
//     double e1 = (a+d)/2 + sqrt(4 * b * c + (a-d) * (a - d))/2;
//     double e2 = (a+d)/2 - sqrt(4 * b * c + (a-d) * (a - d))/2;
//     
//     std::copy( m.begin2(), m.end2(), std::ostream_iterator<double>( std::cout, "\n" ));
//     
//     
//     std::cout << e1 << " " << e2 << "\n";
//     
}