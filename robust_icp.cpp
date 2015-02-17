/*=========================================================================

Program:   Robarts ICP
Module:    $RCSfile: robust_icp.cpp,v $
Creator:   Elvis C. S. Chen <chene@robarts.ca>
Language:  C++
Author:    $Author: Elvis Chen $
Date:      $Date: 2014/03/04 12:51:30 $
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
* FILE:  robust_icp.cpp
*
* Elvis C.S. Chen
* chene@robarts.ca
*
* June 25, 2010
*
* rigid Iterative Closest Point algorithm.  The underlying point-to-point
* registration is based on rigid tukey-SVD solution.
*
*
*/
// C++ includes
#include <cfloat>
#include <cmath>
#include <iostream>
#include <vector>

// ANN includes
#include <ANN/ANN.h>

// local includes
#include "matrix.h"
#include "icp.h"
#include "tukey_pointRegistration.h"

//
// assumes R and t provided the initial guess
//
void robust_icp( Matrix<double> &points, Matrix<double> &model,
                Matrix<double> &R, Vec<double> &t,
                double &rms, double tau )
  {

  if ( points.num_rows() != 3 ||
    model.num_rows() != 3 ) 
    {
    std::cerr << "use column vectors" << std::endl
      << points.num_rows() << "x" << model.num_rows() << std::endl;

    exit(1);
    }

  int nPoints = points.num_cols();
  int nModel = model.num_cols();

  // build a KD-tree of the "model"
  ANNpointArray data_pts;   // data (model) points
  ANNpoint query_pt;        // query point
  ANNidxArray nn_idx;       // near neighbour indices
  ANNdistArray dists;       // near neighbour distances
  ANNkd_tree *model_tree;   // search structure (KD-tree)
  int ANNk = 1;             // find 1 near neighbor
  int ANNdim = 3;           // dimension
  double ANNeps = 0.0;      // exact near neighbor

  query_pt = annAllocPt( ANNdim );
  data_pts = annAllocPts( nModel, ANNdim );
  nn_idx = new ANNidx[ ANNk ];
  dists = new ANNdist[ ANNk ];

  //
  // copy all the model points into ANN structure
  //
  for ( int i = 0; i < nModel; i++ )
    {
    for ( int j = 0; j < 3; j++ ) 
      {
      data_pts[i][j] = model[j][i];
      }
    }

  //
  // initialize the KD tree
  //
  model_tree = new ANNkd_tree( data_pts, nModel, ANNdim );

  int nIter = 0;
  Matrix<double> r( 3, nPoints ), l( 3, nPoints );

  //
  // the ICP loop
  //
  Matrix<double> residual;
  Vec<double> residMag;
  double oldRMS = rms;
  int maxIter = 2000;


  //
  // double checks for R and t
  //
  if ( R.num_cols() != 3 || R.num_cols() != 3 )
    R = eye( 3, 1.0 );

  if ( t.dim() != 3 ) 
    {
    t.newsize( 3 );
    t[0] = t[1] = t[2] = 0.0;
    }


  while ( nIter < maxIter ) 
    {
    // loop for a fixed number of times 'cus we may never converge

    // in this loop, "l" is the transformed points, 
    // "r" is the nearest neighbour in the model

    // transform "points" into l
    l = R * points;
    for ( int i = 0; i < nPoints; i++ ) 
      {
      for ( int j = 0; j < 3; j++ ) 
        {
        l[j][i] += t[j];
        }
      }

    // find the nearest neighbour of each point in "l"
    for ( int i = 0; i < nPoints; i++ ) 
      {
      for ( int j = 0; j < 3; j++ ) 
        {
        query_pt[j] = l[j][i];
        }

      model_tree->annkSearch( query_pt, ANNk, nn_idx, dists, ANNeps );

      for ( int j = 0; j < 3; j++ )
        {
        r[j][i] = data_pts[ nn_idx[0] ][j];
        }
      }

    std::cerr << points.num_cols() << " " << points.num_rows() << std::endl
      << r.num_cols() << " " << r.num_rows() << std::endl;

    robust_heterogeneous_point_register( points, 
      r,
      R,
      t,
      rms,
      residMag,
      FLT_MAX );

    if ( oldRMS == rms ) 
      {
      nIter = maxIter + 1;
      }

    oldRMS = rms;

    if ( rms <= tau )
      {
      nIter = maxIter + 1;
      }

    nIter++;
    }

  
  // clean up
  delete model_tree;
  delete [] dists;
  delete [] nn_idx;

  annClose();
  //annDeallocPts( data_pts );
  //annDeallocPt( query_pt );
  }
