/*=========================================================================

Program:   Robarts ICP
Module:    $RCSfile: robust_icp_main.cpp,v $
Creator:   Elvis C. S. Chen <chene@robarts.ca>
Language:  C++
Author:    $Author: Elvis Chen $
Date:      $Date: 2014/03/04 12:50:30 $
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
#include <iostream>
#include <fstream>
#include <string>
#include <cfloat >

// local includes
#include "matrix.h"
#include "mathUtils.h"
#include "robust_icp.h"


int main( int argc, char *argv[] )
  {
  std::cerr << std::endl
    << "The Iterative Closest Point with Tukey-biweight reweighting (robust-ICP)" << std::endl << std::endl << std::endl;

  if ( argc < 3 )
    {
    std::cerr << "USAGE: " << argv[0]
              << " <points.txt> <model.txt> <initT -- optional> " << std::endl
              << std::endl << std::endl
              << "file format is given as: " << std::endl << std::endl
              << "n" << std::endl
              << "x1 y1 z1" << std::endl
              << "x1 y2 z2" << std::endl
              << "..." << std::endl
              << "xn yn zn" << std::endl;

    exit(1);
    }

  
  double m[16];
  Matrix<double> transform(4,4);;
  bool hasInitTransform = false;
  if ( argc > 2 )
    {
    std::cerr << "reading init" << std::endl;
    std::ifstream tfile;
    tfile.open( argv[3] );

    if ( !tfile.is_open() )
      {
      std::cerr << argv[3] << " not opened" << std::endl;
      exit(1);
      }
    for ( int i = 0; i < 16; i++ )
      {
      tfile >> m[i];
      }
    
    int count=0;
    for ( int i = 0; i < 4; i++ )
      for ( int j = 0; j < 4; j++ )
        transform[i][j] = m[count++];
    
    tfile.close();
    hasInitTransform = true;
    }
  
  // space allocation
  int nPoints, nModels;
  std::ifstream model_file;
  std::ifstream point_file;

  point_file.open( argv[1] );
  model_file.open( argv[2] );

  if ( !model_file.is_open() )
    {
      std::cerr << argv[2] << " not opened" << std::endl;
      exit(1);
    }

  if ( !point_file.is_open() )
    {
      std::cerr << argv[1] << " not opened" << std::endl;
      exit(1);
    }

  /*
   * Parse the input point clouds
   */
  model_file >> nModels;
  point_file >> nPoints;


  Matrix<double> models( 3, nModels ), points( 3, nPoints );

  for ( int i = 0; i < nModels; i++ )
    {
      model_file >> models[0][i]
                 >> models[1][i]
                 >> models[2][i];
    }
  for ( int i = 0; i < nPoints; i++ )
    {
      point_file >> points[0][i]
                 >> points[1][i]
                 >> points[2][i];

    }

  // close the files are they are not needed anymore
  model_file.close();
  point_file.close();

  std::cerr << "Finished reading files" << std::endl;
  
  // initial rotation
  double FRE = 0.0, tau = 1e-9;
  Matrix<double> R( 3, 3 );      // rotation
  Vec<double> t(3);              // translation

  Matrix<double> initialQuaternions; // rotation group
  FourtyRotations( initialQuaternions ); // fourty uniformally sampled rotations

  // loop through all the rotation group
  Vec<double> quat(4), minT(3);
  double minRMS = 0.0;
  Matrix<double> minR(3,3);
  Vec<double> FREMag, minFREMag;

  if ( hasInitTransform )
    { // run only 1 iteration of ICP based on the initial transform
    for ( int i = 0; i < 3; i++ )
      {
      for ( int j = 0; j < 3; j++ )
        {
        R[i][j] = transform[i][j];
        } // j
      t[i] = transform[i][3];
      } // i

    robust_icp( points, models, R, t, FRE, tau );

    minRMS = FRE;
    minR = R;
    minT = t;
    minFREMag = FREMag;
    }
  else
    {
    // go through the rotation group
    for ( int i = 0; i < initialQuaternions.num_cols(); i++ )
      {
      // go through all the rotations
      quat[0] = initialQuaternions[0][i];
      quat[1] = initialQuaternions[1][i];
      quat[2] = initialQuaternions[2][i];
      quat[3] = initialQuaternions[3][i];
      t[0] = t[1] = t[2] = 0.0;
      q2m3x3( quat, R ); // initial guess on rotation
      // translation does not matter much


      robust_icp( points, models, R, t, FRE, tau );

      std::cerr << i << " FRE: " << minRMS << std::endl
        << R << t << std::endl;
      if ( i == 0 )
        { 
        minRMS = FRE;
        minR = R;
        minT = t;
        minFREMag = FREMag;
        }

      if ( FRE < minRMS )
        {
        minRMS = FRE;
        minR = R;
        minT = t;
        minFREMag = FREMag;
        }
      }

    }
  
  std::cerr << "Final answer: " << minR << minT << std::endl;

  
  std::ofstream myfile;
  myfile.open( "robust_FRE.txt" );
  std::cerr << minFREMag.size() << std::endl;
  for ( int i = 0; i < minFREMag.size(); i++ )
    myfile << minFREMag[i] << " ";
  myfile << std::endl;
  myfile.close();
  return(0);


  return(0);
  }
