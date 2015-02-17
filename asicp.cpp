/*=========================================================================

Program:   Robarts ICP
Module:    $RCSfile: asicp.cpp,v $
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
#include <assert.h>
#include <cmath>


// ANN includes
#include <ANN/ANN.h>

// local includes
#include "matrix.h"
#include "asicp.h"
#include "mathUtils.h"
#include "ASOPP_Major.h"
#include "pointRegistration.h"


void asicp_md( Matrix<double> &points,             
              Matrix<double> &model,
              Matrix<double> &R,
              Matrix<double> &A,
              Vec<double> &t,
              double &FRE, double threshold,
              Vec<double> &FREmag )
  {
  assert( points.num_rows() == 3 );
  assert( model.num_rows() == 3 );


  int nPoints = points.num_cols();
  int nModel = model.num_cols();


  // build a KD-tree of the "model"
  ANNpointArray data_pts;   // data (model) points
  ANNpoint query_pt;        // query point
  ANNidxArray nn_idx;       // near neighbour indices
  ANNdistArray dists;       // near neighbour distances
  ANNkd_tree *model_tree=0;   // search structure (KD-tree)
  int ANNk = 1;             // find 1 near neighbor
  int ANNdim = 3;           // dimension
  double ANNeps = 0.0;      // exact near neighbor

  query_pt = annAllocPt( ANNdim );
  data_pts = annAllocPts( nModel, ANNdim );
  nn_idx = new ANNidx[ ANNk ];
  dists = new ANNdist[ ANNk ];

  int nIter = 0;
  Matrix<double> r( 3, nPoints ), l( 3, nPoints ), lprime( 3, nPoints );

  //
  // the ICP loop
  //
  Matrix<double> residual;
  Vec<double> residMag;
  double oldFRE = FRE;
  int maxIter = 2000;


  //
  // estimate the scales based on the eigenvalues of the covariance matrices
  //
  double initScale;
  lprime = R * points;
  for ( int i = 0; i < points.num_cols(); i++ )
    {
    for ( int j = 0; j < points.num_rows(); j++ )
      {
      lprime[j][i] += t[j];
      }
    }
  estimateScalesFromPoints( lprime, model, initScale, R ); // R is given
  A = eye( 3, initScale );

  std::cerr << "Initial Scales: " << A << std::endl;

  Matrix<double> RR(3,3), AA(3,3);
  Vec<double> tt(3);
  double residuals;
  bool changedScales;

  lprime = R * A * points;

  while ( nIter < maxIter )
    {

    std::cerr << "i: " << nIter << std::endl;

    // loop for a fixed number of times 'cus we may never converge

    // in this loop, "l" is the transformed points,
    // "r" is the nearest neighbour in the model
    l = R * A * points; // scaling followed by rotation
    for ( int i = 0; i < nPoints; i++ )
      {
      for ( int j = 0; j < 3; j++ )
        {
        l[j][i] += t[j]; // followed by translation
        } // j
      } // i

    // estimate the scales needed to be used in Mahalanobis distance
    //
    // what we do here is to use the correspondanding fiducials from
    // the previous iteration of ICP and register that to the current
    // iteration to get an idea of the scaling.
    //no
    ASMajor_point_register( lprime, l, RR, AA, tt, residuals, threshold, FREmag );
    // std::cerr << "i2: " << nIter << std::endl;
    lprime = l; // save l

    // find the closest point
    //
    // using the estimates scales from the previous iteration
    // for the calculation of Mahalanobis Distance

    // closestPoint_with_MahalanobisDistance( l, model, AA, r );
    if ( model_tree )
      {
      delete model_tree;
      model_tree = 0;
      }
    Matrix<double> sqrtAA;
    sqrtAA = AA;
    for ( int i = 0; i < 3; i++ )
      {
      sqrtAA[i][i] = sqrt( AA[i][i] );
      } // i

    //
    // copy all the model points into ANN structure
    //
    // NOTE: since we are interested in Mahalabonis Distance,
    // both the data/model points needs to be scaled
    for ( int i = 0; i < nModel; i++ )
      {
      for ( int j = 0; j < 3; j++ )
        {
        data_pts[i][j] = model[j][i] / sqrtAA[j][j];
        } // j
      } // i


    //
    // initialize the KD tree
    //
    model_tree = new ANNkd_tree( data_pts, nModel, ANNdim );

    // find the nearest neighbour of each point in 'l'
    for ( int i = 0; i < nPoints; i++ )
      {
      for ( int j = 0; j < 3; j++ )
        {
        query_pt[j] = l[j][i] / sqrtAA[j][j];
        } // j

      model_tree->annkSearch( query_pt, ANNk, nn_idx, dists, ANNeps );
      for ( int j = 0; j < 3; j++ )
        {
        r[j][i] = data_pts[ nn_idx[0] ][j] * sqrtAA[j][j];
        }
      } // i


    // once the correspondances are found, solve for the
    // anisotropic-scaled orthogonal procrustes analysis
    // using the solution by Dosse and Ten Berge
    ASMajor_point_register( points, r,
      R, A, t,
      FRE, threshold, FREmag );
    std::cerr << "Scales: " << A;


    // need to regularize scaling factors A here, and
    // recompute FRE if necessary.
    //
    //
    // In practice, scaling factors need to be bounded, otherwise
    // trivial (and wrong) solution such as scaling~=0 would occure
    // and the computed FRE would be minimal.

    changedScales = false;
    if ( A[0][0] <  .9 * initScale ) 
      {
      A[0][0] =  .9 * initScale;
      changedScales = true;
      }
    if ( A[0][0] > 1.1 * initScale )
      {
      A[0][0] = 1.1 * initScale;
      changedScales = true;
      }
    if ( A[1][1] <  .9 * initScale ) 
      {
      A[1][1] =  .9 * initScale;
      changedScales = true;
      }
    if ( A[1][1] > 1.1 * initScale ) 
      {
      A[1][1] = 1.1 * initScale;
      changedScales = true;
      }
    if ( A[2][2] <  .9 * initScale ) 
      {
      A[2][2] =  .9 * initScale;
      changedScales = true;
      }
    if ( A[2][2] > 1.1 * initScale ) 
      {
      A[2][2] = 1.1 * initScale;
      changedScales = true;
      }

    //
    // scales has been changed, so re-compute the rotation and FRE
    //
    if ( changedScales )
      {
      std::cerr << "Scales changed" << std::endl;
      l = A * points;

      point_register( l, r, R, t, FRE, threshold, FREmag );
      }

    if ( oldFRE == FRE )
      {
      nIter = maxIter + 1;
      }

    oldFRE = FRE;

    if ( FRE <= threshold )
      {
      nIter = maxIter + 1;
      }

    nIter++;

    } // while

  // kd tree cleanup
  delete [] nn_idx;
  delete [] dists;

  if ( model_tree )
    {
    delete model_tree;
    model_tree = 0;
    }
  annClose();
  }
