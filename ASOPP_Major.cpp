/*=========================================================================

Program:   Robarts Computer Vision
Module:    $RCSfile: ASOPP_Major.cpp,v $
Creator:   Elvis C. S. Chen <chene@robarts.ca>
Language:  C++
Author:    $Author: Elvis Chen $
Date:      $Date: 2013/09/09 17:07:30 $
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

// matrix include
#include "matrix.h"
#include "svd.h"

// local include
#include "ASOPP_Major.h"
#include "mathUtils.h"




void ASMajor_point_register( Matrix<double> &XX, Matrix<double> &YY,
                             Matrix<double> &Q, 
			     Matrix<double> &A,
			     Vec<double> &t,
			     double &FRE,
                 Vec<double> &FREMag )
{
  double threshold = 1e-9;
  ASMajor_point_register( XX, YY, Q, A, t, FRE, threshold, FREMag );
}


void ASMajor_point_register( Matrix<double> &XX, Matrix<double> &YY,
                             Matrix<double> &Q,
			     Matrix<double> &A,
			     Vec<double> &t,
			     double &FRE,
			     double threshold,
                 Vec<double> &FREMag )
{

  // double check the dimensions
  int XXx = XX.num_rows();
  int XXy = XX.num_cols();
  int YYx = YY.num_rows();
  int YYy = YY.num_cols();

  if ( XXx != 3 || YYx != 3 ) 
    { 
      std::cerr << "inputs must use column vectors" << std::endl;
      exit(1);
    }
  if ( XXy != YYy ) 
    { 
      std::cerr << "input point clouds must contain the same number of points" << std::endl;
      exit(1);
    }


  // for simplicity, we'll follow the notation of Dossee and Berge to use
  // row-vectors instead of column vectors.

  Matrix<double> X = transpose(XX);
  Matrix<double> Y = transpose(YY);


  // find the centroids
  Vec<double> Xbar( 3, 0.0 ), Ybar( 3, 0.0 );
  for ( int i = 0; i < X.num_rows(); i++ )
    {
      for ( int j = 0; j < X.num_cols(); j++ )
	{
	  Xbar[j] += X[i][j];
	  Ybar[j] += Y[i][j];
	}
    }
  for ( int i = 0; i < 3; i++ )
    {
      Xbar[i] /= (double)X.num_rows();
      Ybar[i] /= (double)Y.num_rows();
    }

  // translate the input by the centroids
  Matrix<double> Xtilde( X.num_rows(), X.num_cols() );
  Matrix<double> Ytilde( Y.num_rows(), Y.num_cols() );
  for ( int i = 0; i < X.num_rows(); i++ )
    {
      for ( int j = 0; j < X.num_cols(); j++ )
	{
	  Xtilde[i][j] = X[i][j] - Xbar[j];
	  Ytilde[i][j] = Y[i][j] - Ybar[j];
	}
    }


  // normlize Xtilde (page 115 of the paper)
  Vec<double> S(3, 0.0);
  for ( int i = 0; i < X.num_rows(); i++ )
    {
      for ( int j = 0; j < X.num_cols(); j++ )
	{
	  S[j] += ( Xtilde[i][j] * Xtilde[i][j] );
	}
    }
  for ( int i = 0; i < 3; i++ ) S[i] = sqrt( S[i] );

  for ( int i = 0; i < X.num_rows(); i++ )
    {
      for ( int j = 0; j < X.num_cols(); j++ )
        {
	  Xtilde[i][j] /= S[j];
        }
    }

  // compute the cross correlation matrix
  Matrix<double> B = transpose( Ytilde ) * Xtilde;

  // compute the intial rotation
  Matrix<double> U, V;
  svdcmp( B, S, U, V );

  Matrix<double> dd = eye( 3, 1.0 ), UV = U * V;
  dd[2][2] = m3x3_det( UV );
  Q = U * dd * transpose( V );
  // find the residual
  Matrix<double> FREvect = Xtilde * transpose( Q ) - Ytilde;

  FRE = 0.0;
  for ( int i = 0; i < FREvect.num_rows(); i++ )
    {
      FRE += ( FREvect[i][0] * FREvect[i][0] +
	       FREvect[i][1] * FREvect[i][1] +
	       FREvect[i][2] * FREvect[i][2] );

    }
  FRE = sqrt( FRE/(double)FREvect.num_rows() );

  double FRE_orig = 2.0 * ( FRE + threshold );

  Matrix<double> QB(3,3), I(3,3,0.0);
  while ( fabs( FRE_orig - FRE ) > threshold )
    {
      QB = transpose( Q ) * B;
      I[0][0] = QB[0][0];
      I[1][1] = QB[1][1];
      I[2][2] = QB[2][2];
      
      svdcmp( B*I, S, U, V );
      UV = U * V;
      dd[2][2] = m3x3_det( UV );
      
      Q = U * dd * transpose( V );
      
      // calculate the residual
      FREvect = Xtilde * transpose(Q) - Ytilde;
      FRE_orig = FRE;
      FRE = 0.0;
      for ( int i = 0; i < FREvect.num_rows(); i++ )
	{
	  FRE += ( FREvect[i][0] * FREvect[i][0] +
		   FREvect[i][1] * FREvect[i][1] +
		   FREvect[i][2] * FREvect[i][2] );
	}
      FRE = sqrt( FRE/(double)FREvect.num_rows() );

    }
    
  for ( int i = 0; i < Xtilde.num_rows(); i++ )
    {
      for ( int j = 0; j < Xtilde.num_cols(); j++ )
	{
	  Xtilde[i][j] = X[i][j] - Xbar[j];
	} 
    }
  B = transpose( Ytilde ) * Xtilde;


  // matlab:  A = diag( diag( B'*Q ) ./ diag( Xtilde' * Xtilde ) );
  //
  // reuse U and V
  U = transpose( B ) * Q;
  V = transpose( Xtilde ) * Xtilde;
  A.newsize( 3, 3, 0.0 );

  for ( int i = 0; i < 3; i++ )
    {
      A[i][i] = U[i][i] / V[i][i];
    }


  // now calculate the translation
  t.newsize( 3 );
  t[0] = t[1] = t[2] = 0.0;
  U = YY - Q * A * XX; // reuse U

  for ( int i = 0; i < U.num_cols(); i++ )
    {
      for ( int j = 0; j < 3; j++ )
	{
	  t[j] += U[j][i];
	}
    }
  for ( int i = 0; i < 3; i ++ )
    {
      t[i] /= (double)U.num_cols();
    }


  // calculate the FRE and the translation
  FREvect = Ytilde - Xtilde * A * transpose(Q);

  FRE = 0.0;

  if ( FREMag.size() != FREvect.num_rows() )
    {
    FREMag.newsize( FREvect.num_rows() );
    }

  for ( int i = 0; i < FREvect.num_rows(); i++ )
    {
    FREMag[i] = sqrt( FREvect[i][0] * FREvect[i][0] +
                      FREvect[i][1] * FREvect[i][1] +
                      FREvect[i][2] * FREvect[i][2] );
    FRE += FREMag[i];
    }
  FRE = ( FRE/(double)FREvect.num_rows() );

  t = Ybar - Q * A * Xbar;

}

