/*=========================================================================

Program:   Robarts ICP
Module:    $RCSfile: mathUtils.h,v $
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

#ifndef __MATHUTILS_H__
#define __MATHUTILS_H__

#include "matrix.h"


// 
// determinant of a 3x3 matrix
//
double m3x3_det( Matrix<double> &m );


// convert a quaternion to 3x3 rotation matrix
//
// a quaternion is defined as:
//
// q[0] = v1*sin(phi/2)
// q[1] = v2*sin(phi/2)
// q[2] = v3*sin(phi/2)
// q[3] =    cos(phi/2)
void q2m3x3( Vec<double> &qin, Matrix<double> &m );

//
// generate 40 rotations that includes the tetrahedral and
// the octahedral/hexahedral group as per the original ICP paper
// by Besl and McKay (page 247)
//
//
// q is a 4x40 quaterion matrix where each column of q is an unit quaternion
void FourtyRotations( Matrix<double> &q );
void SixtyRotations( Matrix<double> &q );

//
// find the closest points using Mahalanobis distance
//
// X, Y are point clouds using column vectors
// S is the covariance matrix, and in this particular case, diagonal
// out is the closest points of X in Y, hence out has the same dimension as X
void closestPoint_with_MahalanobisDistance( Matrix<double> &X,
                                            Matrix<double> &Y,
                                            Matrix<double> &S,
                                            Matrix<double> &out );
void closestPoint_with_EuclideanDistance( Matrix<double> &X,
                                          Matrix<double> &Y,
                                          Matrix<double> &out );

void calFREMag( Matrix<double> &X,
                Matrix<double> &Y,
                Vec<double> &FREMag );

//
// ordering of 3 numbers
//
// using the box trick
//
void ordering3Numbers( double a, double b, double c,
                       double &min, double &mid, double &max );

void estimateScalesFromPoints( Matrix<double> &p, Matrix<double> &m,
			       double &initScale,
			       Matrix<double> &initR );


// find the median value
double findMedian( Vec<double> &m );

// find the mean value
double findMean( Vec<double> &m );

// find the MAD, median/mean absolute deviation
double findMAD( Vec<double> &m, bool useMedian );

// small utility to calculate the exit condition
double calcConfig( Matrix<double> &Xnew, Matrix<double> &Xold );

#endif // of __MATHUTILS_H__
