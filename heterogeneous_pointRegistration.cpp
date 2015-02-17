/*=========================================================================

Program:   Robarts ICP
Module:    $RCSfile: heterogeneous_pointRegistration.cpp,v $
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

/*
 * FILE:  heterogeneous_pointRegistration.cpp
 *
 * Elvis C.S. Chen
 * chene@robarts.ca
 *
 * June 25, 2010
 *
 * Implementation to the iterative solution for the
 * heterogeneous-weighted Absolute Orientation Problem using the SVD method.
 *
 *
 */

// C++ includes
#include <cfloat>
#include <cmath>
#include <iostream>
#include <assert.h>

// matrix include
#include "matrix.h"
#include "svd.h"

// local include
#include "heterogeneous_pointRegistration.h"
#include "mathUtils.h"

//
// solves for the Absolute Orientation problem with heterogeneous
// weighting using the SVD solution
//

void heterogeneous_point_register( Matrix<double> &l, Matrix<double> &r,
                                   Matrix<double> &R, 
                                   Vec<double> &t,
                                   double &rms, double threshold,
                                   Vec<double> &residMag,
                                   Vec<double> &w )
{
  int lx = l.num_rows();
  int ly = l.num_cols();
  
  int rx = r.num_rows();
  int ry = r.num_cols();


  // double check the dimensions
  assert( lx == 3 );
  assert( rx == 3 );
  assert( ly == ry );
  assert( ly == w.dim() );
  
  //
  // find the weighted centroids of the input points
  //
  Vec<double> cl( 3, 0.0 ), cr( 3, 0.0 );
  double cw = 0.0;
  for ( int i = 0; i < ly; i++ ) 
    {
      cw += w[i];
      for ( int j = 0; j < 3; j++ ) 
	{
	  cl[j] += w[i] * l[j][i];
	  cr[j] += w[i] * r[j][i];
	}
    }

  for ( int i = 0; i < 3; i++ )
    {
      cl[i] /= cw;
      cr[i] /= cw;
    }

  //
  // translate the input points by their weighted centroids
  //
  Matrix<double> lprime( 3, ly ), rprime( 3, ry );
  for ( int i = 0; i < ly; i++ ) 
    {
      for ( int j = 0; j < 3; j++ ) 
	{
	  lprime[j][i] = l[j][i] - cl[j];
	  rprime[j][i] = r[j][i] - cr[j];
	}
    }

  //
  // compute the weighted cross correlation matrix
  //
  // Matrix<double> M = lprime * transpose( rprime );
  Matrix<double> M(3,3);
  for ( int i = 0; i < lx; i++ )
    {
      for ( int k = 0; k < rx; k++ ) 
	{
	  M[i][k] = 0.0;
	  for ( int j = 0; j < ly; j++ ) 
	    {
	      M[i][k] += ( w[j] * lprime[i][j] * rprime[k][j] );
	    }
	}
    }
  
  //
  // compute the SVD
  //
  Matrix<double> U, V;
  Vec<double> S;
  svdcmp( M, S, U, V );
  
  Matrix<double> dd = eye( 3, 1.0 ), VU = V * U;
  dd[2][2] = m3x3_det( VU );
  
  //
  // and the rotation is
  // 
  R = V * dd * transpose( U );
  
  //
  // the translation is:
  //
  t = cr - R * cl;

  
  
  Matrix<double> FREvect;
  FREvect = rprime - R * lprime;
  rms = 0.0;
  
  if ( residMag.size() != FREvect.num_cols() )
    { 
    residMag.newsize( FREvect.num_cols() );
    }

  for ( int i = 0; i < FREvect.num_cols(); i++ )
    {
    residMag[i] = sqrt( FREvect[0][i] * FREvect[0][i] +
                        FREvect[1][i] * FREvect[1][i] +
                        FREvect[2][i] * FREvect[2][i] );
    rms += residMag[i];
    } // i

  rms = ( rms/(double)FREvect.num_cols() );
}
