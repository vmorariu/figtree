//-------------------------------------------------------------------
// File    : MexFigtreeChooseEvaluateMethod.cpp 
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

//The gateway function

void mexFunction(int nlhs,                // Number of left hand side (output) arguments
                 mxArray *plhs[],         // Array of left hand side arguments
                 int nrhs,                // Number of right hand side (input) arguments
                 const mxArray *prhs[])   // Array of right hand side arguments
{
  //check for proper number of arguments 
 
  if(nrhs < 5) mexErrMsgTxt("at least 5 input  arguments required.");
  if(nlhs < 1) mexErrMsgTxt("at least 1 output argument  required.");

  //////////////////////////////////////////////////////////////
  // Input arguments
  //////////////////////////////////////////////////////////////
  
  int d, N, M, W;

  // x
  //  The 2D array is column-major: each column represents a point.
  int argu = 0;
  double *x = mxGetPr(prhs[argu]);
  int mrows = (int)  mxGetM(prhs[argu]); //mrows
  int ncols = (int) mxGetN(prhs[argu]); //ncols
  d = mrows;
  N = ncols;
  //if ( mrows != d && ncols != N)  mexErrMsgTxt("Input 'x' must be a d x N matrix");
  if( d <= 0 || N <= 0 ) 
    mexErrMsgTxt( "Input 'x' must be a d x N matrix, with d and N nonzero" );

  // h
  argu = 1;
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
    mexErrMsgTxt("Input 'h' must be a scalar.");
  double h = (double) mxGetScalar(prhs[argu]);
  if ( h <= 0.0) 
    mexErrMsgTxt("Input 'h' must be a positive number.");
 
  //  W
  argu = 2;
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
    mexErrMsgTxt("Input 'W' must be a scalar.");
  W = (int) mxGetScalar(prhs[argu]);
  if (W <= 0) 
    mexErrMsgTxt("Input 'W' must be a positive number.");

  // y
  argu = 3;
  double * y = mxGetPr(prhs[argu]);
  mrows = (int)mxGetM(prhs[argu]); //mrows
  ncols = (int)mxGetN(prhs[argu]); //ncols
  M = ncols;
  if( mrows != d ) 
    mexErrMsgTxt("Input 'y' must be a d x M matrix, with d the same as the d x N matrix 'x'.");
  if(  M <= 0 ) 
    mexErrMsgTxt( "Input 'y' must be a d x M matrix, with d and M nonzero" );

  // epsilon
  argu = 4;
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
    mexErrMsgTxt("Input 'epsilon' must be a scalar.");
  double epsilon = (double) mxGetScalar(prhs[argu]);
  if ( epsilon <= 0.0) 
    mexErrMsgTxt("Input 'epsilon' must be a positive number.");

  //////////////////////////////////////////////////////////////
  // Optional input arguments
  //////////////////////////////////////////////////////////////
  int ifgtParamMethod = FIGTREE_PARAM_NON_UNIFORM;
  int verbose = 0;

  argu++;
  if( nrhs >= argu+1 )
  {
    if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
      mexErrMsgTxt("Input 'ifgtParamMethod' must be an integer between 0 and 1" );
    ifgtParamMethod = (int) mxGetScalar(prhs[argu]);
    if( ifgtParamMethod < 0 || ifgtParamMethod >= FIGTREE_TRUNC_SIZE ) 
      mexErrMsgTxt("Input 'ifgtParamMethod' must be an integer between 0 and 1" );
  }

  argu++;
  if( nrhs >= argu+1 )
  {
    if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
      mexErrMsgTxt("Input 'verbose' must be an integer.");
    verbose = (int) mxGetScalar(prhs[argu]);
  }

  //////////////////////////////////////////////////////////////
  // Output arguments
  //////////////////////////////////////////////////////////////

  // g
  int bestEvalMethod = 0;
  double * flops = NULL;
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  if( nlhs > 1 )
  {
    plhs[1] = mxCreateDoubleMatrix(4,1,mxREAL);
    flops = mxGetPr(plhs[1]);
  }

  double * bestEvalMethodDouble = mxGetPr(plhs[0]);

  // Function call
  figtreeChooseEvaluationMethod( d, N, M, W, x, h, y, epsilon, ifgtParamMethod, verbose, &bestEvalMethod, flops );
  *bestEvalMethodDouble = (double)bestEvalMethod;

  return; 
}
