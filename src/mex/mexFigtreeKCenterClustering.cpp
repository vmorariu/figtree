//-------------------------------------------------------------------
// File    : MexIfgtEvaluateDirectTree.cpp 
// Purpose : Interface between MATLAB and C++
// Author  : Vlad I. Morariu (morariu(at)umd(.)edu)
// Date    : 2007-06-25
//
// This file has been derived by Vlad I. Morariu from Vikas C.
// Raykar's original IFGT code (source file mexmain.cpp). 
//-------------------------------------------------------------------

//-------------------------------------------------------------------
// The code was written by Vlad Morariu and Vikas C. Raykar 
// and is copyrighted under the Lesser GPL: 
//
// Copyright (C) 2007 Vlad Morariu and Vikas C. Raykar
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 or later.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details. 
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, 
// MA 02111-1307, USA.  
//
// The author may be contacted via email at: 
// morariu(at)umd(.)edu, vikas(at)cs(.)umd(.)edu 
//-------------------------------------------------------------------


// preprocessor definitions used: FIGTREE_STATIC, USE_MATLAB_MEX
// The first makes sure that the compiled code does not require Ifgt.dll, but instead
// compiles the code directly in the code.  The second replaces all printf calls
// with mxPrintf and also causes memory operations such as new and delete to be managed
// by matlab via mxMalloc and mxFree.

#include "mex.h"
#include "figtree.h"
#include "string.h" // for memcpy

//The gateway function

void mexFunction(int nlhs,				// Number of left hand side (output) arguments
				 mxArray *plhs[],		// Array of left hand side arguments
				 int nrhs,              // Number of right hand side (input) arguments
				 const mxArray *prhs[])  // Array of right hand side arguments
{

  //check for proper number of arguments 
 
  if(nrhs != 4) mexErrMsgTxt("4 inputs required.");
  if(nlhs != 6) mexErrMsgTxt("6 outputs required.");

   //////////////////////////////////////////////////////////////
  // Input arguments
  //////////////////////////////////////////////////////////////
  
  // d
  int argu = 0;
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
    mexErrMsgTxt("Input 'd' must be a scalar.");
  int d = (int)
    mxGetScalar(prhs[argu]);
  if (d <= 0) 
    mexErrMsgTxt("Input 'd' must be a positive number.");

  // N
  argu++;
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
    mexErrMsgTxt("Input 'N' must be a scalar.");
  int N = (int) mxGetScalar(prhs[argu]);
  if (N <= 0) 
    mexErrMsgTxt("Input 'N' must be a positive number.");

  // x
  argu++;
  double *x = mxGetPr(prhs[argu]);
  int mrows = (int)mxGetM(prhs[argu]); //mrows
  int ncols = (int)mxGetN(prhs[argu]); //ncols
  if ( mrows != d && ncols != N )  
    mexErrMsgTxt("Input 'x' must be a d x N matrix.");
  
  // kMax
  argu++;
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
    mexErrMsgTxt("Input 'kMax' must be a scalar.");
  int kMax = (int) mxGetScalar(prhs[argu]);
  if (kMax <= 0) 
    mexErrMsgTxt("Input 'kMax' must be a positive number.");

  //////////////////////////////////////////////////////////////
  // Output arguments
  //////////////////////////////////////////////////////////////

  // K
  argu = 0;
  plhs[argu] = mxCreateDoubleMatrix(1, 1, mxREAL );
  double * K = (double*)mxGetPr(plhs[argu]);

  // rx
  argu++;
  plhs[argu] = mxCreateDoubleMatrix(1, 1, mxREAL);
  double * rx = (double*)mxGetPr(plhs[argu]);
  
  // clusterIndex
  argu++;
  plhs[argu] = mxCreateNumericMatrix(1, N, mxUINT32_CLASS, mxREAL);
  int * clusterIndex = (int*)mxGetPr(plhs[argu]);

  int kTemp; 
  int * numPointsTemp = (int*)mxMalloc( sizeof(int)*kMax );
  double * clusterCentersTemp = (double*)mxMalloc( sizeof(double)*kMax*d );
  double * clusterRadiiTemp = (double*)mxMalloc( sizeof(double)*kMax );

  // function call
  figtreeKCenterClustering( d, N, x, kMax, &kTemp, rx, clusterIndex, clusterCentersTemp, numPointsTemp, clusterRadiiTemp );
  *K = kTemp;

  // The sizes of these arrays depend on the number of final clusters, which can be less than kMax. 
  // clusterCenters  
  argu++;
  plhs[argu] = mxCreateDoubleMatrix(d, kTemp, mxREAL);
  double * clusterCenters = mxGetPr(plhs[argu]);
  memcpy( clusterCenters, clusterCentersTemp, sizeof(double)*kTemp*d );

  // numPoints
  argu++;
  plhs[argu] = mxCreateNumericMatrix(1, kTemp, mxUINT32_CLASS, mxREAL);
  int * numPoints = (int*)mxGetPr(plhs[argu]);
  memcpy( numPoints, numPointsTemp, sizeof(int)*kTemp );

  // clusterRadii
  argu++;
  plhs[argu] = mxCreateDoubleMatrix(1, kTemp, mxREAL);
  double * clusterRadii = mxGetPr(plhs[argu]);
  memcpy( clusterRadii, clusterRadiiTemp, sizeof(double)*kTemp );

  // free temporary arrays
  mxFree( clusterCentersTemp );
  mxFree( numPointsTemp );
  mxFree( clusterRadiiTemp );

  return;
  
}
