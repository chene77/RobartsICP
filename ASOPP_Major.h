/*=========================================================================

Program:   Robarts Computer Vision
Module:    $RCSfile: ASOPP_Major.h,v $
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

/*
 * FILE: ASOPP_Major.h
 *
 * Elvis C.S. Chen
 * chene@robarts.ca
 *
 * January, 17, 2014
 *
 * This function implements the PRESCALING ANISOTROPIC Orthogonal Procrustes
 * Analysis algorithm based on the Majorization Principal algorithm presented
 * by Mohammed Bennani Dosse and Jos Ten Berge:
 *
 *
 *  @Article{springerlink:10.1007/s00357-010-9046-8,
 *    Author         = {Bennani Dosse, Mohammed and Ten Berge, Jos},
 *    Title          = {Anisotropic Orthogonal Procrustes Analysis},
 *    Journal        = {Journal of Classification},
 *    Volume         = {27},
 *    Pages          = {111-128},
 *    Note           = {10.1007/s00357-010-9046-8},
 *    affiliation    = {University of Rennes 2, Place du Recteur Henri Le
 *                     Moal, CS 24307 35043 Rennes Cedex, France},
 *    issn           = {0176-4268},
 *    issue          = {1},
 *    keyword        = {Computer Science},
 *    publisher      = {Springer New York},
 *    url            = {http://dx.doi.org/10.1007/s00357-010-9046-8},
 *    year           = 2010
 *  }
 *
 * We use FRE as the stopping criteria
 *
 *
 *
 *
 * The BR algorithm basically is a loop that solves for Rotation
 * and scaling independently in each iteration.  We use FRE as the
 * stopping criteria.
 *
 *
 * X and Y are array of 3D points where
 *
 *     size(X,2) = size(Y,2), and
 *     size(X,1) = size(Y,1) = 3
 *
 *
 * optional inputs:
 *     threshold -- the stopping criteria based on FRE
 *                  defaults to 1e-9
 *
 *
 * OUTPUTS:
 *     R, the orthonormal (rotation), 3x3
 *     A, the anisotropic scaling , 3x3
 *     t, the translation, 3x1
 *     FRE, Fiducial Registration Error, the residual
 *
 * We attempt to solve for the PRESCALING (scaling before rotation) case:
 *
 *     Y = R * A * X + t
 *
 * Based on my original Matlab code written on Oct. 22nd, 2010
 *
 */
#ifndef __ASOPP_MAJOR_H__
#define __ASOPP_MAJOR_H__

#include "matrix.h"

// with user-supplied threshold/termination value
void ASMajor_point_register( Matrix<double> &x, Matrix<double> &y,
			     Matrix<double> &Q,
			     Matrix<double> &A,
			     Vec<double> &t,
			     double &FRE, double threshold,
                 Vec<double> &FREMag );



// with the default threshold set to 1e-9
void ASMajor_point_register( Matrix<double> &x, Matrix<double> &y,
                             Matrix<double> &Q,
                             Matrix<double> &A,
                             Vec<double> &t,
                             double &FRE,
                             Vec<double> &FREMag );

#endif // of __ASOPP_MAJOR_H__
