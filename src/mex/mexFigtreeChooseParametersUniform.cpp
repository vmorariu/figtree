//-------------------------------------------------------------------
// File    : MexFigtreeChooseParametersUniform.cpp 
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

// preprocessor definitions used: FIGTREE_STATIC, FIGTREE_USE_MATLAB_MEX
// The first makes sure that the compiled code does not require figtree.dll, but instead
// compiles the code directly in the code.  The second replaces all printf calls
// with mxPrintf and also causes memory operations such as new and delete to be managed
// by matlab via mxMalloc and mxFree.

#include "mex.h"
#include "figtree.h"


//The gateway function

void mexFunction(int nlhs,                // Number of left hand side (output) arguments
                 mxArray *plhs[],         // Array of left hand side arguments
                 int nrhs,                // Number of right hand side (input) arguments
                 const mxArray *prhs[])   // Array of right hand side arguments
{
  // check for proper number of arguments 
  if(nrhs != 5) 
    mexErrMsgTxt("5 input  arguments required.");
  if(nlhs != 3) 
    mexErrMsgTxt("3 output arguments required.");

  //////////////////////////////////////////////////////////////
  // Input arguments
  //////////////////////////////////////////////////////////////
  
  // d
  int argu = 0;
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu]) != 1 ) 
    mexErrMsgTxt("Input 'd' must be a scalar.");
  //  get the scalar input d
  int d = (int)mxGetScalar(prhs[argu]);
  if (d <= 0) 
    mexErrMsgTxt("Input 'd' must be a positive number.");

  // h
  argu++;
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu]) != 1 ) 
    mexErrMsgTxt("Input 'h' must be a scalar.");
  double h = (double)mxGetScalar(prhs[argu]);
  if ( h <= 0.0) 
    mexErrMsgTxt("Input 'h' must be a positive number.");

  // epsilon
  argu++;
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu]) != 1 ) 
    mexErrMsgTxt("Input 'epsilon' must be a scalar.");
  double epsilon = (double)mxGetScalar(prhs[argu]);
  if ( epsilon <= 0.0) 
    mexErrMsgTxt("Input 'epsilon' must be a positive number.");

  // kLimit
  argu++;
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu]) != 1 ) 
    mexErrMsgTxt("Input 'kLimit' must be a scalar.");
  int kLimit = (int)mxGetScalar(prhs[argu]);
  if (kLimit <= 0) 
    mexErrMsgTxt("Input 'kLimit' must be a positive number.");

  // maxRange
  argu++;
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu]) != 1 ) 
    mexErrMsgTxt("Input 'maxRange' must be a scalar.");
  double maxRange = (double)mxGetScalar(prhs[argu]);
  if ( maxRange <= 0.0) 
    mexErrMsgTxt("Input 'maxRange' must be a positive number.");

  //////////////////////////////////////////////////////////////
  // Output arguments
  //////////////////////////////////////////////////////////////

  // K
  argu = 0;
  plhs[argu] = mxCreateDoubleMatrix(1,1,mxREAL);
  double *K = mxGetPr(plhs[argu]);

  // pMax
  argu++;
  plhs[argu] = mxCreateDoubleMatrix(1,1,mxREAL);
  double *pMax = mxGetPr(plhs[argu]);

  // r
  argu++;
  plhs[argu] = mxCreateDoubleMatrix(1,1,mxREAL);
  double *r = mxGetPr(plhs[argu]);

  // Function call
  int kTmp, pMaxTmp;
  figtreeChooseParametersUniform( d, h, epsilon, kLimit, maxRange, &kTmp, &pMaxTmp, r );
  *K = kTmp;
  *pMax = pMaxTmp;

  return; 
}
