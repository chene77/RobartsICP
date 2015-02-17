/*=========================================================================

Program:   Robarts ICP
Module:    $RCSfile: tukey_pointRegistration.cpp,v $
Creator:   Elvis C. S. Chen <chene@robarts.ca>
Language:  C++
Author:    $Author: Elvis Chen $
Date:      $Date: 2014/03/04 10:45:30 $
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
#include <cfloat>
#include <cmath>
#include <iostream>
#include <vector>
#include <assert.h>

// matrix
#include "matrix.h"

// local includes
#include "tukey_pointRegistration.h"
#include "heterogeneous_pointRegistration.h"
#include "mathUtils.h"



//
//
void robust_heterogeneous_point_register( Matrix<double> &l, 
                                          Matrix<double> &r,
                                          Matrix<double> &R,
                                          Vec<double> &t,
                                          double &rms,
                                          Vec<double> &residMag,
                                          double sigma )
  {
  int nmax = 200;
  int lx = l.num_rows();
  int ly = l.num_cols();
  int rx = r.num_rows();
  int ry = r.num_cols();

  // double check the input dimensions
  assert( lx == 3 );
  assert( rx == 3 );
  assert( ly == ry );

  Vec<double> w( ly, 1.0 );
  double threshold = 1e-6;
  double config_change = FLT_MAX;

  // use the regular SVD solution to obtain an initial alignment
  heterogeneous_point_register( l, r, R, t, rms, threshold, residMag, w );

  
  //
  // transformed "l" points
  //
  Matrix<double> lold = R*l, lnew;
  for ( int i = 0; i < ly; i++ )
    for ( int j = 0; j < 3; j++ )
      lold[j][i] += t[j];
  
  double scale, median;
  Vec<double> x( ly ); // temp storage


  //
  // iteratively use Tukey biweight to re-weight the point-correspondance
  //
  // exit when no further improvement can be obtained
  //
  while ( config_change > threshold )
    {

    scale = findMAD( residMag, false );
    median = findMedian( residMag );

    if ( scale < sigma && sigma != FLT_MAX )
      scale = sigma;

    // assign the weights
    for ( int i = 0; i < ly; i++ )
      x[i] = residMag[i] - median;

    w = TukeyG( x, scale );

    heterogeneous_point_register( l, r, R, t, rms, threshold, residMag, w );


    lnew = R*l;
    for ( int i = 0; i < ly; i++ )
      for ( int j = 0; j < 3; j++ )
        lnew[j][i] += t[j];

    config_change = calcConfig( lnew, lold );
    lold = lnew;

    } // while
  }

//
// Tukey's biweight
//
Vec<double> TukeyRho( Vec<double> &z, double sigma ) 
  {
  Vec<double> p( z.dim() );

  double ss = sigma * sigma;
  double c = ss/6.0;
  double t;

  for ( int i = 0; i < z.dim(); i++ ) 
    {
    if ( fabs( z[i] ) <= sigma ) 
      {
      t = 1.0 - z[i] * z[i] / ss;
      p[i] = c * ( 1.0 - t * t * t );
      }
    else 
      {
      p[i] = c;
      }
    }

  return( p );
  }

Vec<double> TukeyPsi( Vec<double> &z, double sigma ) 
  {
  //
  // tukey PSI function is the tukey-G * z
  //
  Vec<double> w = TukeyG( z, sigma );

  for ( int i = 0; i < z.dim(); i++ ) 
    {
    if ( fabs( z[i] ) <= sigma ) 
      {
      w[i] = z[i] * w[i];
      }
    }

  return( w );
  }

Vec<double> TukeyG( Vec<double> &z, double sigma ) 
  {
  Vec<double> g( z.dim() );
  double t;
  double s = sigma * sigma;

  for ( int i = 0; i < z.dim(); i++ ) 
    {      
    if ( fabs( z[i] ) <= sigma ) 
      {
      t = z[i] * z[i] / s;
      g[i] = ( 1.0 - t ) * ( 1.0 - t );
      }
    else
      {
      g[i] = 0.0;
      }
    }

  return( g );
  }
