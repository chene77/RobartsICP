//
// Elvis Chen
// chene@cs.queensu.ca
//
// Department of Computing and Information Science
// Queen's University, Kingston, Ontario, Canada
//
// Feb. 20, 2000
//

//
// Filename:  svd.h
//

// Initial implementation of SVD (singular value decomposition)
//

// given a matrix a[1..m][1..n], this routine computes its
// singular value decomposition, a = u*w*transpose(t).
// U is the same dimention as a, w is output vector of size [1..n],
// and the matrix v (not the transpose transpose(v)) is output as
// v[1..n][1..n]
//
// algorithms are taken from Numerical Recipes
//
//
// Usage:
//
// Matrix<double> A( 3, 3,
//                   " 7 2 3 "
//                   " 16 55 7 "
//                   " 7 8 9 " );
//
// Matrix<double> U, V;
// Vec<double> S;
//
// svdcmp(A, S, U, V);
//

#ifndef __SVD_H__
#define __SVD_H__

#include <cmath>
#include "matrix.h"

#ifndef SIGN
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#endif

template<class T>
T MAX( T a, T b ) 
{
  T arg1 = a, arg2 = b;
  
  return ( (arg1 > arg2) ? arg1 : arg2 );
}

template<class T>
T MIN( T a, T b ) 
{
  T arg1 = a, arg2 = b;
  
  return ( (arg1 < arg2) ? arg1 : arg2 );
}


template<class T>
T pythag( T a, T b ) 
{
  // computes (a^2 + b^2)^(1/2) without destructive underflow or verflow

  T absa = fabs(a);
  T absb = fabs(b);

  if (absa > absb)
    return ( absa * sqrt(1.0 + (absb/absa)*(absb/absa)) );
  else
    return ( absb == 0.0 ? 0.0 : absb*sqrt(1.0 + (absa/absb)*(absa/absb)) );
}


template< class T >
void svdcmp( const Matrix<T> &u, Vec<T> &w, Matrix<T> &a, Matrix<T> &v ) 
{

  // "u" is the given input matrix
  // "w" is a vector of size n that contains the singular values
  //     of the diagonal matrix
  // 
  // u = a*w*transpose(t);
  //
  Subscript m = u.num_rows();
  Subscript n = u.num_cols();
  Subscript i, j, jj, its, k, l, nm;
  int flag;
  T anorm, c, f, g, h, s, scale, x, y, z;
  
  if ( ( v.num_rows() != m ) || ( v.num_cols() != n ) ) v.newsize( m, n );
  if ( w.dim() != n ) w.newsize( n );
  
  Vec<T> rv1(n);
  
  a = u;
  
  g = scale = anorm = 0.0;

  //
  // Householder reduction to bidiagonal form.
  //
  
  for (i = 1; i <= n; i++) {
    l = i + 1;
    rv1(i) = scale * g;
    g = s = scale = 0.0;
    
    if (i <= m) {
      for (k = i; k <= m; k++) scale += fabs(a(k,i));
      
      if (scale) {
	for (k = i; k <= m; k++) {
	  a(k,i) /= scale;
	  s += a(k,i) * a(k,i);
	}
	
	f = a(i, i);
	g = -SIGN(sqrt(s), f);
	
	h = f * g - s;
	a(i,i) = f - g;
	for (j = l; j <= n; j++) {
	  for (s = 0.0, k = i; k <= m; k++) s += a(k,i) * a(k,j);
	  f = s/h;
	  for (k = i; k <= m; k++) a(k,j) += f * a(k,i);
	}
	for (k = i; k <= m; k++) a(k,i) *= scale;
      }
    }
    
    w(i) = scale * g;
    g = s = scale = 0.0;
    
    if ( (i <= m) && (i != n) ) {
      for (k = l; k <= n; k++) scale += fabs(a(i,k));
      if (scale) {
	for (k = l; k <= n; k++){
	  a(i,k) /= scale;
	  s += a(i,k) * a(i,k);
	}
	
	f = a(i,l);
	g = -SIGN(sqrt(s), f);
	h = f * g - s;
	a(i,l) = f - g;
	for (k = l; k <= n; k++) rv1(k) = a(i,k) / h;
	for (j = l; j <= m; j++) {
	  for (s = 0.0, k = l; k <= n; k++) s += a(j,k) * a(i,k);
	  for (k = l; k <= n; k++) a(j,k) += s * rv1(k);
	}
	for (k = l; k <= n; k++) a(i,k) *= scale;
      }
    }
    anorm = MAX(anorm, (fabs(w(i)) + fabs(rv1(i))));
  }
  

  //
  // Accumulation of right-hand transformations
  //
  for (i = n; i >= 1; i--) {
    if (i < n) {
      if (g) {
	for (j = l; j <= n; j++) // double division to avoid possible underflow
	  v(j,i) = (a(i,j) / a(i,l) ) / g;
	
	for (j = l; j <=n; j++) {
	  for (s = 0.0, k = l; k <= n; k++) s += a(i,k) * v(k,j);
	  for (k = l; k <= n; k++) v(k,j) += s*v(k,i);
	}
      }
      for (j = l; j <= n; j++) v(i,j) = v(j,i) = 0.0;
    }
    v(i,i) = 1.0;
    g = rv1(i);
    l = i;
  }
  
  for (i = MIN(m,n); i >= 1; i--) { 
    // Accumulation of left-hand transformations
    l = i + 1;
    g = w(i);
    for (j = l; j <= n; j++) a(i,j) = 0.0;
    if (g) {
      g = 1.0 / g;
      for (j = l; j <= n; j++) {
	for (s = 0.0, k = l; k <= m; k++) s += a(k,i) * a(k,j);
	f = (s / a(i,i)) * g;
	for (k = i; k <= m; k++) a(k,j) += f * a(k,i);
      }
      for (j = i; j <= m; j++) a(j,i) *= g;
    }
    else for (j = i; j <= m; j++) a(j,i) = 0.0;
    ++a(i,i);
  }
  
  //
  // Diagonalization of the bidiagonal form:
  // Loop over singular values, and over allowed iterations
  //
  for (k = n; k >= 1; k--) {
    for (its = 1; its <= 30; its++) {
      flag = 1;
      for (l = k; l >= 1; l--) {
	
	// test for splitting
	// note that rv1(1) is always zero

	nm = l - 1;
	if ( (T)(fabs(rv1(l))+anorm) == anorm ) {
	  flag = 0;
	  break;
	}
	if ( (T)(fabs(w(nm)) + anorm) == anorm) break;
      }
      
      if (flag) {
	c = 0.0;
	s = 1.0;
	for (i = l; i <= k; i++) {
	  f = s*rv1(i);
	  rv1(i) = c * rv1(i);
	  
	  if ( (T)(fabs(f) + anorm) == anorm ) break;
	  g = w(i);
	  h = pythag( f, g );
	  
	  w(i) = h;
	  h = 1.0 / h;
	  c = g * h;
	  s = -f * h;
	  
	  for (j = 1; j <= m; j++) {
	    y = a(j, nm);
	    z = a(j, i);
	    a(j, nm) = y * c + z * s;
	    a(j, i) = z * c - y * s;
	  }
	}
      }
      
      z = w(k);
      if ( l == k ) {
	// Convergence
	// singular value is made nonnegative
	if ( z < 0.0 ) {
	  w(k) = -z;
	  for (j = 1; j <= n; j++) v(j,k) = -v(j,k);
	}
	break;
      }
      
      assert ( its != 30 ); // check this
      
      // if (its == 30) break;

      x = w(l);
      nm = k - 1;
      y = w(nm);
      g = rv1(nm);
      h = rv1(k);
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = pythag( f, 1.0 );
      f = ((x - z) * (x + z) + h*((y / (f+SIGN(g,f))) - h)) / x;
      c = s = 1.0;

      // Next QR transformation
      for (j = l; j <= nm; j++) {
	i = j + 1;
	g = rv1(i);
	y = w(i);
	h = s * g;
	g = c * g;
	z = pythag(f, h);
	rv1(j) = z;
	c = f / z;
	s = h / z;
	f = x * c + g * s;
	g = g * c - x * s;
	h = y * s;
	y *= c;
	for (jj = 1; jj <= n; jj++) {
	  x = v(jj,j);
	  z = v(jj,i);
	  v(jj,j) = x * c + z * s;
	  v(jj,i) = z * c - x * s;
	}
	z = pythag( f, h );
	w(j) = z;  // Rotation can be arbitrary if z = 0
	if (z) {
	  z = 1.0 / z;
	  c = f * z;
	  s = h * z;
	}
	f = c * g + s * y;
	x = c * y - s * g;
	
	for (jj = 1; jj <= m; jj++) {
	  y = a(jj,j);
	  z = a(jj,i);
	  a(jj,j) = y * c + z * s;
	  a(jj,i) = z * c - y * s;
	}
      }
      rv1(l) = 0.0;
      rv1(k) = f;
      w(k) = x;
    }
  }
}

#endif // of __SVD_H__
