/*=========================================================================

Program:   Robarts ICP
Module:    $RCSfile: mathUtils.cpp,v $
Creator:   Elvis C. S. Chen <chene@robarts.ca>
Language:  C++
Author:    $Author: Elvis Chen $
Date:      $Date: 2014/03/03 13:36:30 $
Version:   $Revision: 0.99 $

==========================================================================

This is an open-source copyright license based on BSD 2-Clause License,
see http://opensource.org/licenses/BSD-2-Clause

Copyright (c) 2013, Elvis Chia-Shin Chen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/


// C++ includes
#include <cmath>
#include <assert.h>

#include "mathUtils.h"
#include "jacobi.h"
#include "qsort.h"

// 
// determinant of a 3x3 matrix
//
double m3x3_det( Matrix<double> &m ) 
  {
  if ( m.num_cols() != 3 ||
    m.num_rows() != 3 ) 
    {
    std::cerr << "matrix is not 3x3" << std::endl;
    return( 0 );
    }
  else 
    {
    return( ( m[0][0] * m[1][1] * m[2][2] +
              m[0][1] * m[1][2] * m[2][0] +
              m[0][2] * m[1][0] * m[2][1] -
              m[0][0] * m[1][2] * m[2][1] -
              m[0][1] * m[1][0] * m[2][2] -
              m[0][2] * m[1][1] * m[2][0] ) );
    }
  }

// convert a quaternion to 3x3 rotation matrix
//
// a quaternion is defined as:
//
// q[0] = v1*sin(phi/2)
// q[1] = v2*sin(phi/2)
// q[2] = v3*sin(phi/2)
// q[3] =    cos(phi/2)
void q2m3x3( Vec<double> &qin, Matrix<double> &m )
{
  int nq = qin.dim();
  if ( nq != 4 )
    {
      std::cerr << "Expecting the quternion to be a 4x1 vector" << std::endl;
      exit(1);
    }


  // normalize the quaternion
  double l = sqrt( qin[0]*qin[0] +
		   qin[1]*qin[1] +
		   qin[2]*qin[2] +
		   qin[3]*qin[3] );

  Vec<double> q = 1.0/l * qin;

  m.newsize( 3, 3 );
  double xx = q[0]*q[0];
  double yy = q[1]*q[1];
  double zz = q[2]*q[2];

  double xy = q[0]*q[1];
  double xz = q[0]*q[2];

  double yz = q[1]*q[2];

  double wx = q[3]*q[0];
  double wy = q[3]*q[1];
  double wz = q[3]*q[2];

  
  m[0][0] = 1. - 2. * (yy + zz);
  m[0][1] =      2. * (xy - wz);
  m[0][2] =      2. * (xz + wy);

  m[1][0] =      2. * (xy + wz);
  m[1][1] = 1. - 2. * (xx + zz);
  m[1][2] =      2. * (yz - wx);

  m[2][0] =      2. * (xz - wy);
  m[2][1] =      2. * (yz + wx);
  m[2][2] = 1. - 2. * (xx + yy);
}


void estimateScalesFromPoints( Matrix<double> &p, Matrix<double> &m,
			       double &initScale,
			       Matrix<double> &initR )
{
  int nP = p.num_cols();
  int nM = m.num_cols();


  Matrix<double> l( 3, nP );

  // perform the initial rotation
  l = initR * p;


  Vec<double> cl( 3, 0.0 ), cr( 3, 0.0 );
  // find the centroid of both clouds
  for ( int i = 0; i < nP; i++ )
    {
    for ( int j = 0; j < 3; j++ ) 
      {
	cl[j] += l[j][i];
      }
    }

  for ( int i = 0; i < nM; i++ )
    {
      for ( int j = 0; j < 3; j++ )
	{
	  cr[j] += m[j][i];
	}
    }
  for ( int i = 0; i < 3; i++ )
    {
      cl[i] /= (double)nP;
      cr[i] /= (double)nM;
    }


  //
  // translate the input points by their centroids
  Matrix<double> lprime( 3, nP), rprime( 3, nM );
  for ( int i = 0; i < nP; i++ )
    {
      for ( int j = 0; j < 3; j++ )
	{
	  lprime[j][i] = l[j][i] - cl[j];
	}
    }
  for ( int i = 0; i < nM; i++ )
    {
      for ( int j = 0; j < 3; j++ )
	{
	  rprime[j][i] = m[j][i] - cr[j];
	}
    }


  // compute the covariance matrix
  Matrix<double> mu( 3, 3 ), lambda( 3, 3 );

  //
  // eigen velues are propotional to the number of the points
  //
  mu = 1.0/(double)nP * lprime * transpose( lprime );
  lambda = 1.0/(double)nM * rprime * transpose( rprime );


  Matrix<double> mu_eigvec, mu_eigval;
  Matrix<double> lambda_eigvec, lambda_eigval;
  int tempi;

  jacobi( mu, mu_eigval, mu_eigvec, tempi );
  jacobi( lambda, lambda_eigval, lambda_eigvec, tempi );

  double mu_min, mu_mid, mu_max;
  double lambda_min, lambda_mid, lambda_max;

  ordering3Numbers( mu_eigval[0][0], mu_eigval[1][0], mu_eigval[2][0],
		    mu_min, mu_mid, mu_max );
  ordering3Numbers( lambda_eigval[0][0], lambda_eigval[1][0], lambda_eigval[2][0],
		    lambda_min, lambda_mid, lambda_max );

  /*
  std::cerr << mu_min << " " << mu_mid << " " << mu_max << std::endl;
  std::cerr << lambda_min << " " << lambda_mid << " " << lambda_max << std::endl;
  */

  initScale = ( lambda_min/mu_min + lambda_mid/mu_mid + lambda_max/mu_max ) / 3.0;
}

//
// generate 40 rotations that includes the tetrahedral and
// the octahedral/hexahedral group as per the original ICP paper
// by Besl and McKay (page 247)
//
//
// q is a 4x40 quaterion matrix where each column of q is an unit quaternion
void FourtyRotations( Matrix<double> &q )
  {
  Matrix<double> tempq( 4, 54 );

  Vec<double> q0(2,"1 0");
  Vec<double> q1(3,"1 0 -1");
  Vec<double> q2(3,"1 0 -1");
  Vec<double> q3(3,"1 0 -1");

  int counter = 0;
  Matrix<double> temp(4,54);
  double ll = 1.0;
  for ( int i = 0; i < q0.dim(); i++ )
    for ( int j = 0; j < q1.dim(); j++ )
      for ( int k = 0; k < q2.dim(); k++ )
        for ( int l = 0; l < q3.dim(); l++ )
          {
          /*
          ll = 1.0/sqrt( q0[i]*q0[i] +
          q1[j]*q1[j] +
          q2[k]*q2[k] +
          q3[l]*q3[l] );
          ll=1.0;  // don't bother normalize it here
          */
          tempq[3][counter] = ll * q0[i]; // scalar component of the quaternion
          tempq[0][counter] = ll * q1[j];
          tempq[1][counter] = ll * q2[k];
          tempq[2][counter] = ll * q3[l];

          counter++;
          }


        // 41st quaternion is not valid
        // 42nd t0 54th are duplicates
        q.newsize(4,40);
        for ( int i = 0; i < 40; i++ )
          for ( int j = 0; j < 4; j++ )
            {
            q[j][i] = tempq[j][i];
            }
  }



//
// generate 60 rotations that includes the tetrahedral and
// the octahedral/hexahedral group as per the original ICP paper
// by Besl and McKay (page 247)
//
//
// q is a 4x60 quaterion matrix where each column of q is an unit quaternion
void SixtyRotations( Matrix<double> &q )
  {
  q.newsize( 4, 60 );

  double a = ( sqrt(5.0) - 1.0 ) / 4.0;
  double b = .5;
  double c = 1.0 / sqrt(2.0);
  double d = ( sqrt(5.0) + 1.0 ) / 4.0;

  q[1][ 0] = 1.0;  q[2][ 0] = 0.0;  q[3][ 0] = 0.0;  q[0][ 0] = 0.0;
  q[1][ 1] = 0.0;  q[2][ 1] = 1.0;  q[3][ 1] = 0.0;  q[0][ 1] = 0.0;
  q[1][ 2] = 0.0;  q[2][ 2] = 0.0;  q[3][ 2] = 1.0;  q[0][ 2] = 0.0;
  q[1][ 3] = 0.0;  q[2][ 3] = 0.0;  q[3][ 3] = 0.0;  q[0][ 3] = 1.0;

  q[1][ 4] = 0.0;  q[2][ 4] =   a;  q[3][ 4] =   b;  q[0][ 4] =   d;
  q[1][ 5] = 0.0;  q[2][ 5] =   a;  q[3][ 5] =   b;  q[0][ 5] =  -d;
  q[1][ 6] = 0.0;  q[2][ 6] =   a;  q[3][ 6] =  -b;  q[0][ 6] =   d;
  q[1][ 7] = 0.0;  q[2][ 7] =   a;  q[3][ 7] =  -b;  q[0][ 7] =  -d;

  q[1][ 8] = 0.0;  q[2][ 8] =   b;  q[3][ 8] =   d;  q[0][ 8] =   a;
  q[1][ 9] = 0.0;  q[2][ 9] =   b;  q[3][ 9] =   d;  q[0][ 9] =  -a;
  q[1][10] = 0.0;  q[2][10] =   b;  q[3][10] =  -d;  q[0][10] =   a;
  q[1][11] = 0.0;  q[2][11] =   b;  q[3][11] =  -d;  q[0][11] =  -a;

  q[1][12] = 0.0;  q[2][12] =   d;  q[3][12] =   a;  q[0][12] =   b;
  q[1][13] = 0.0;  q[2][13] =   d;  q[3][13] =   a;  q[0][13] =  -b;
  q[1][14] = 0.0;  q[2][14] =   d;  q[3][14] =  -a;  q[0][14] =   b;
  q[1][15] = 0.0;  q[2][15] =   d;  q[3][15] =  -a;  q[0][15] =  -b;

  q[1][16] =   a;  q[2][16] = 0.0;  q[3][16] =   d;  q[0][16] =   b;
  q[1][17] =   a;  q[2][17] = 0.0;  q[3][17] =   d;  q[0][17] =  -b;
  q[1][18] =   a;  q[2][18] = 0.0;  q[3][18] =  -d;  q[0][18] =   b;
  q[1][19] =   a;  q[2][19] = 0.0;  q[3][19] =  -d;  q[0][19] =  -b;

  q[1][20] =   b;  q[2][20] = 0.0;  q[3][20] =   a;  q[0][20] =   d;
  q[1][21] =   b;  q[2][21] = 0.0;  q[3][21] =   a;  q[0][21] =  -d;
  q[1][22] =   b;  q[2][22] = 0.0;  q[3][22] =  -a;  q[0][22] =   d;
  q[1][23] =   b;  q[2][23] = 0.0;  q[3][23] =  -a;  q[0][23] =  -d;

  q[1][24] =   d;  q[2][24] = 0.0;  q[3][24] =   b;  q[0][24] =   a;
  q[1][25] =   d;  q[2][25] = 0.0;  q[3][25] =   b;  q[0][25] =  -a;
  q[1][26] =   d;  q[2][26] = 0.0;  q[3][26] =  -b;  q[0][26] =   a;
  q[1][27] =   d;  q[2][27] = 0.0;  q[3][27] =  -b;  q[0][27] =  -a;

  q[1][28] =   a;  q[2][28] =   b;  q[3][28] = 0.0;  q[0][28] =   d;
  q[1][29] =   a;  q[2][29] =   b;  q[3][29] = 0.0;  q[0][29] =  -d;
  q[1][30] =   a;  q[2][30] =  -b;  q[3][30] = 0.0;  q[0][30] =   d;
  q[1][31] =   a;  q[2][31] =  -b;  q[3][31] = 0.0;  q[0][31] =  -d;

  q[1][32] =   b;  q[2][32] =   d;  q[3][32] = 0.0;  q[0][32] =   a;
  q[1][33] =   b;  q[2][33] =   d;  q[3][33] = 0.0;  q[0][33] =  -a;
  q[1][34] =   b;  q[2][34] =  -d;  q[3][34] = 0.0;  q[0][34] =   a;
  q[1][35] =   b;  q[2][35] =  -d;  q[3][35] = 0.0;  q[0][35] =  -a;

  q[1][36] =   d;  q[2][36] =   a;  q[3][36] = 0.0;  q[0][36] =   b;
  q[1][37] =   d;  q[2][37] =   a;  q[3][37] = 0.0;  q[0][37] =  -b;
  q[1][38] =   d;  q[2][38] =  -a;  q[3][38] = 0.0;  q[0][38] =   b;
  q[1][39] =   d;  q[2][39] =  -a;  q[3][39] = 0.0;  q[0][39] =  -b;

  q[1][40] =   a;  q[2][40] =   d;  q[3][40] =   b;  q[0][40] = 0.0;
  q[1][41] =   a;  q[2][41] =   d;  q[3][41] =  -b;  q[0][41] = 0.0;
  q[1][42] =   a;  q[2][42] =  -d;  q[3][42] =   b;  q[0][42] = 0.0;
  q[1][43] =   a;  q[2][43] =  -d;  q[3][43] =  -b;  q[0][43] = 0.0;

  q[1][44] =   b;  q[2][44] =   a;  q[3][44] =   d;  q[0][44] = 0.0;
  q[1][45] =   b;  q[2][45] =   a;  q[3][45] =  -d;  q[0][45] = 0.0;
  q[1][46] =   b;  q[2][46] =  -a;  q[3][46] =   d;  q[0][46] = 0.0;
  q[1][47] =   b;  q[2][47] =  -a;  q[3][47] =  -d;  q[0][47] = 0.0;

  q[1][48] =   d;  q[2][48] =   b;  q[3][48] =   a;  q[0][48] = 0.0;
  q[1][49] =   d;  q[2][49] =   b;  q[3][49] =  -a;  q[0][49] = 0.0;
  q[1][50] =   d;  q[2][50] =  -b;  q[3][50] =   a;  q[0][50] = 0.0;
  q[1][51] =   d;  q[2][51] =  -b;  q[3][51] =  -a;  q[0][51] = 0.0;

  q[1][52] =   b;  q[2][52] =   b;  q[3][52] =   b;  q[0][52] =   b;
  q[1][53] =   b;  q[2][53] =   b;  q[3][53] =   b;  q[0][53] =  -b;
  q[1][54] =   b;  q[2][54] =   b;  q[3][54] =  -b;  q[0][54] =   b;
  q[1][55] =   b;  q[2][55] =   b;  q[3][55] =  -b;  q[0][55] =  -b;

  q[1][56] =   b;  q[2][56] =  -b;  q[3][56] =   b;  q[0][56] =   b;
  q[1][57] =   b;  q[2][57] =  -b;  q[3][57] =   b;  q[0][57] =  -b;
  q[1][58] =   b;  q[2][58] =  -b;  q[3][58] =  -b;  q[0][58] =   b;
  q[1][59] =   b;  q[2][59] =  -b;  q[3][59] =  -b;  q[0][59] =  -b;

  }


void calFREMag( Matrix<double> &X,
                Matrix<double> &Y,
                Vec<double> &FREMag )
  {
  int Xx = X.num_rows();
  int Xy = X.num_cols();
  int Yx = Y.num_rows();
  int Yy = Y.num_cols();

  assert( Xx == Yx ); // make sure we are using column vectors and in 3D
  assert( Xx == 3 );
  assert( Xy == Yy ); // make sure we have homologous points

  if ( FREMag.size() != Xy )
    {
    FREMag.newsize( Xy );
    }

  double d;
  for ( int i = 0; i < Xy; i++ )
    {
    d = 0.0;
    for ( int j = 0; j < Xx; j++ )
      {
      d += ( X[j][i] - Y[j][i] ) * ( X[j][i] - Y[j][i] );
      } // j
    FREMag[i] = sqrt(d);
    }//i

  }
void closestPoint_with_EuclideanDistance( Matrix<double> &X,
                                          Matrix<double> &Y,
                                          Matrix<double> &out )
  {
  int Xx = X.num_rows();
  int Xy = X.num_cols();
  int Yx = Y.num_rows();
  int Yy = Y.num_cols();

  assert( Xx == Yx ); // make sure we are using column vectors and in 3D
  assert( Xx == 3 );

  if ( out.num_rows() != Xx ||
       out.num_cols() != Xy )
    { // resize output if necessary
    out.newsize( Xx, Xy );
    }

  double d, mind;
  int minidx;
  Vec<double> v( Xx );
  
  for ( int i = 0; i < Xy; i++ ) // loop through all X
    {
    d = 0.0;
    for ( int k = 0; k < Xx; k++ ) // compute the vector
      {
      v[k] = X[k][i] - Y[k][0]; // 1st point in Y
      
       d += ( v[k] * v[k] ); // squared Mahalanobis Distance
      }
    mind = sqrt( d );
    minidx = 0;

    // from 2nd point to the end
    for ( int j = 1; j < Yy; j++ )
      {
      d = 0.0;
      for ( int k = 0; k < Xx; k++ )
        {
        v[k] = X[k][i] - Y[k][j];

        d += ( v[k] * v[k] );
        }
      d = sqrt( d );
      if ( d < mind )
        {
        mind = d;
        minidx = j;
        }
      }

    for ( int k = 0; k < Xx; k++ ) // copy the closest point in Y to out
      {
      out[k][i] = Y[k][ minidx ];
      }
    }
  }
//
// find the closest points using Mahalanobis distance
//
// X, Y are point clouds using column vectors
// S is the covariance matrix, and in this particular case, diagonal
// out is the closest points of X in Y, hence out has the same dimension as X
void closestPoint_with_MahalanobisDistance( Matrix<double> &X,
                                            Matrix<double> &Y,
                                            Matrix<double> &S,
                                            Matrix<double> &out )
  {
  int Xx = X.num_rows();
  int Xy = X.num_cols();
  int Yx = Y.num_rows();
  int Yy = Y.num_cols();
  int Sx = S.num_rows();
  int Sy = S.num_cols();

  assert( Xx == Yx ); // make sure we are using column vectors and in 3D
  assert( Xx == 3 );

  assert( Sx == 3 && Sy == 3 );

  if ( out.num_rows() != Xx ||
       out.num_cols() != Xy )
    { // resize output if necessary
    out.newsize( Xx, Xy );
    }

  double d, mind;
  int minidx;
  Vec<double> v( Xx );
  
  for ( int i = 0; i < Xy; i++ ) // loop through all X
    {
    d = 0.0;
    for ( int k = 0; k < Xx; k++ ) // compute the vector
      {
      v[k] = X[k][i] - Y[k][0]; // 1st point in Y
      
       d += ( v[k] * v[k] ) / S[k][k]; // squared Mahalanobis Distance
      }
    mind = sqrt( d );
    minidx = 0;

    // from 2nd point to the end
    for ( int j = 1; j < Yy; j++ )
      {
      d = 0.0;
      for ( int k = 0; k < Xx; k++ )
        {
        v[k] = X[k][i] - Y[k][j];

        d += ( v[k] * v[k] ) / S[k][k];
        }
      d = sqrt( d );
      if ( d < mind )
        {
        mind = d;
        minidx = j;
        }
      }

    for ( int k = 0; k < Xx; k++ ) // copy the closest point in Y to out
      {
      out[k][i] = Y[k][ minidx ];
      }
    }
  }


//
// ordering of 3 numbers
//
// using the box trick
//
void ordering3Numbers( double a, double b, double c,
                       double &min, double &mid, double &max )
  {
  if ( a < b )
    {
    if ( a < c )
      {
      if ( b < c )
        {
        min = a;
        mid = b;
        max = c;
        }
      else
        {
        min = a;
        mid = c;
        max = c;
        } 
      }
    else
      {
      min = c;
      mid = a;
      max = b;
      }
    }
  else
    {
    if ( b < c )
      {
      if ( a < c )
        {
        min = b;
        mid = a; 
        max = c;
        }
      else
        {
        min = b;
        mid = c;
        max = a;
        }
      }
    else
      {
      min = c;
      mid = b;
      max = a;
      }
    }
  }

//
// find the mean value
//
double findMean( Vec<double> &m ) 
{  
  double v = 0.0;
  for ( int i = 0; i < m.dim(); i++ )
    {
      v += m[i];
    }
  
  return( v/(double)m.dim() );
}

//
// find the median value of a 1D array
//
double findMedian( Vec<double> &m ) 
{
  //
  // index-sort the input vector using qsort
  // 
  Vec<int> sorted_idx = qsort( m );

  int idx, idxx = -1;
  if ( m.dim() % 2 == 0 )
    {
      idx = m.dim() / 2 - 1;
      idxx = idx + 1;
    }
  else
    {
      idx = int( floor( m.dim()/2.0 ) );
    }
  
  if ( idxx == -1 )
    {
      // odd number of elements
      return ( m[ sorted_idx[ idx ] ] );
    }
  else 
    {
      // even number of elements
      return( ( m[ sorted_idx[ idx ] ] +
		m[ sorted_idx[ idxx ] ] ) / 2.0 );
    }
}
//
// find the MAD, median/mean absolute deviation
//
// y = mad(x)
//
// y = median(abs(x-median(x)))
// y = mean(abs(x-mean(x)))
//
double findMAD( Vec<double> &m, bool useMedian ) 
{  
  double mm;
  Vec<double> dem( m.dim() );
  if ( useMedian )
    mm = findMedian( m );
  else
    mm = findMean( m );
  
  for ( int i = 0; i < m.dim(); i++ )
    dem[i] = fabs( m[i] - mm );
  
  if ( useMedian )
    return( findMedian( dem ) );
  else
    return( findMean( dem ) );
}


double calcConfig( Matrix<double> &Xnew, Matrix<double> &Xold )
{
  int nx = Xnew.num_rows();
  int ny = Xnew.num_cols();
  int ox = Xold.num_rows();
  int oy = Xold.num_cols();
  
  if ( nx != 3 || ox != 3 ) 
    {
      std::cerr << "use column vectors" << std::endl;
      exit(1);
    }
  if ( ny != oy ) 
    {
      std::cerr << "number of points don't match " << std::endl;
      exit(1);
    }
  

  Vec<double> XoldMean(3,0.0);  
  for ( int i = 0; i < ny; i++ ) 
    for ( int j = 0; j < 3; j++ ) 
      XoldMean[j] += Xold[j][i];

  for ( int j = 0; j < 3; j++ ) XoldMean[j] /= ny;
  
  double S = 0.0, T = 0.0;
  double x, y, z;
  
  for ( int i = 0; i < ny; i++ ) 
    {
      x = Xnew[0][i] - Xold[0][i];
      y = Xnew[1][i] - Xold[1][i];
      z = Xnew[2][i] - Xold[2][i];
      
      S += ( x*x + y*y + z*z );
      
      x = Xold[0][i] - XoldMean[0];
      y = Xold[1][i] - XoldMean[1];
      z = Xold[2][i] - XoldMean[2];
   
      T += ( x*x + y*y + z*z );
    }
  
  return ( sqrt( S/T ) );
}
