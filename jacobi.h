//
// Elvis Chen
// chene@cs.queensu.ca
// 
// Department of Computing and Information Science
// Queen's University, Kingston, Ontario, Canada
//
// August 05, 2000
//

//
// Filename:  jacobi.h
//
  
// Initial implementation of Jacobi Transformation
// that is used to find the eigen value/vector of
// a real symmetric matrix
//
// Algorithms are taken from Numerical Recipes
//

#ifndef __JACOBI_H__
#define __JACOBI_H__

#include <cmath>
#include <cstdlib>

template< class MAT >
void rot( MAT &a, const double s, const double tau, const int i,
	  const int j, const int k, const int l ) 
{
  double g, h;
  g = a[i][j];
  h = a[k][l];
  a[i][j] = g - s * ( h + g*tau );
  a[k][l] = h + s * ( g - h*tau );
}

template< class MAT >
void jacobi( MAT &a, MAT &d, MAT &v, int &nrot ) 
{
  //  Computes all eigen values and eigenvectors of a real symmetric
  // matrix a[0..n-1][0..n-1].  On output, elements of a above 
  // the diagonals are destroyed.  d[0][0..n-1] returns the
  // eigenvalues of a.  v[0..n-1][0..n-1] is a matrix whose columns contain, 
  // on output, the normalizes eivenvectors of a.  nrot returns
  // the number of Jacobi rotations that were required.  
    
  int i=0, j, ip, iq;
  double tresh, theta, tau, t, sm, s, h, g, c;
  
  int x = a.num_rows();
  int y = a.num_cols();
  
  if ( x == y ) {
    // make sure we have a square matrix first
    d.newsize( x, 1 );
    v.newsize( x, x );
    int n = x;
    MAT b( n, 1 ), z( n, 1 );
    
    for ( ip = 0; ip < n; ip ++ ) { // initialize to the identiay matrix
      for ( iq = 0; iq < n; iq ++ ) v[ip][iq] = 0.0;
      v[ip][ip] = 1.0;
    }
    
    for ( ip = 0; ip < n; ip++ ) { // initialize b and d to the
                                   // diagonal of a
      b[ip][0] = d[ip][0] = a[ip][ip];
      z[ip][0] = 0.0; // this vector will accumulate terms of the form
                      // ta_pq as in equation 11.1.14
    }
    nrot = 0;
    for ( i = 0; i <= 50; i++ ) {
      sm = 0.0;
      for ( ip = 0; ip < (n-1); ip++ ) {  // sum magnitude of off-diagonal 
	  for ( iq = (ip+1); iq < n; iq++ ) {  // elements
	    sm += fabs( a[ip][iq] );
	  }
      }
      
      if ( sm == 0.0 ) // the normal return, which relies on quadratic
        return; // convergence to machine underflow

      if ( i < 4 )
        tresh = 0.2 * sm / ( n*n ); // on the first three sweeps
      else
        tresh = 0.0; // .. thereafter.

      for ( ip = 0; ip < (n-1); ip ++ ) {
        for ( iq = (ip+1); iq < n; iq++ ) {
          g = 100.0 * fabs( a[ip][iq] );
	  // after 4 sweeps, skip the rotation if the off-diagonal
          // element is small

	  if ( (i > 4) &&
	       ( (fabs(d[ip][0])+g) == fabs(d[ip][0]) ) &&
	       ( (fabs(d[iq][0])+g) == fabs(d[iq][0]) ) )
            a[ip][iq] = 0.0;
	  else if ( fabs( a[ip][iq] ) > tresh ) {
            h = d[iq][0] - d[ip][0];
	    if ( (fabs(h)+g ) == fabs(h) )
              t = ( a[ip][iq] )/h;
	    else {
              theta = 0.5 * h / ( a[ip][iq] );
	      t = 1.0 / ( fabs(theta) + sqrt( 1.0+theta*theta ) );
	      if ( theta < 0.0 ) t = -t;
	    }
	    c = 1.0 / sqrt( 1.0 + t*t );
	    s = t * c;
	    tau = s/(1.0 + c);
	    h = t * a[ip][iq];
	    z[ip][0] -= h;
	    z[iq][0] += h;
	    d[ip][0] -= h;
	    d[iq][0] += h;
	    a[ip][iq] = 0.0;
	    for ( j = 0; j < ip; j++ ) 
		rot( a, s, tau, j, ip, j, iq );
	    for ( j = (ip+1); j < iq; j++ )
		rot( a, s, tau, ip, j, j, iq );
	    for ( j = (iq+1); j < n; j++ )
		rot( a, s, tau, ip, j, iq, j );
	    for ( j = 0; j < n; j++ ) 
		rot( v, s, tau, j, ip, j, iq );
	    ++nrot;
	  }
	}
      }
      for ( ip = 0; ip < n; ip++ ) {
	  b[ip][0] += z[ip][0];
	  d[ip][0] = b[ip][0];
	  z[ip][0] = 0.0;
      }
    }
    std::cerr << "Too many iterations in Jacobi" << std::endl;
    exit(1);
  }
}

#endif // of __JACOBI_H__
