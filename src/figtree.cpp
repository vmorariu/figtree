// File: figtree.cpp
// Created:  11-03-06 by Vlad Morariu
//
// Modified:  6-22-07 by Vlad Morariu 
//   Initial changes from previous version of the IFGT code (written by Vikas C.
//   Raykar and Changjiang Yang) and FIGTree code (written by Vikas C. Raykar).  
//
//   Modifications include:
//   1) Code can compile into a dynamic library that provides C-style interface
//      without requiring Matlab.
//
//   2) Added an improved parameter selection method that removes assumption that 
//      sources are uniformly distributed (observed large speedup in cases where 
//      sources were not uniformly distributed, and often little slowdown from 
//      overhead when the sources were actually uniformly distributed).
//
//   3) Changed the IFGT code to take multiple sets of weights for same set of 
//      sources and targets instead of having to call IFGT code multiple times.  
//      By computing a set of coefficients for each weight set, much overhead is 
//      saved (eg. computing monomials, and so on), resulting in significant 
//      speedup.
//
//   4) Added function (figtree()) that performs all parameter selection/clustering 
//      using any choice of parameter selection and evaluation algorithms.
// 
//   5) Some bugs/problem cases were fixed (some bugs caused seg faults, others 
//      were certain problem cases that could result in bad parameter selection 
//      and, as a result, memory allocation errors or seg faults).
//
//   6) In the original implementation, most code resided in the constructor and 
//      Evaluate() functions of a class, and was actually called in sequential 
//      order as if it were a C function (thus not using any real advantages of 
//      C++ classes).  Thus, all code except for that of KCenterClustering, which
//      seems to fit better in a class, has been put in C-style functions inside 
//      of figtree.cpp.  The original location of the original source is indicated 
//      in figtree.cpp before each function.
//  
//   7) Stylistic changes (eg. variable naming conventions, function names, ...)
//    
// Modified:  9-23-07 by Vlad Morariu 
//   Change code to compile on linux and solaris.
//
// Modified: 10-03-07 by Vlad Morariu 
//   Remove requirement that data is in unit hypercube by adding 
//   maxRange parameter to figtreeChoose* functions.
//
// Modified: 01-22-08 by Vlad Morariu 
//   Rename library to FIGTree (and some
//   other function remanimg)
//
// Modified: 02-20-08 by Vlad Morariu
//   Added nchoosek_double function to use 'double' instead of 'int' to prevent
//   overflow issues.  The overflow would cause incorrect parameter estimation
//   which then resulted in out of memory errors.
//
// Modified: 02-21-08 by Vlad Morariu
//   Allow rx to be zero (each pt has a cluster center on it), and allow
//   figtreeChooseParametersNonUniform to choose a value of K that gives rx=0. 
//   In some cases in higher dimensions, it is significantly cheaper to have
//   a center at each pt (i.e. rx=0) than having even one cluster with nonzero
//   radius (since the radius might require excessively high pMax).
//   Also added FIGTREE_CHECK_POS_DOUBLE macro to allow rx to be zero when 
//   checking input parameters.
//

//------------------------------------------------------------------------------
// The code was written by Vlad Morariu, Vikas Raykar, and Changjiang Yang 
// and is copyrighted under the Lesser GPL: 
//
// Copyright (C) 2008 Vlad Morariu and Vikas Raykar and Changjiang Yang 
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
// morariu(at)umd(.)edu, vikas(at)umiacs(.)umd(.)edu, cyang(at)sarnoff(.)com
//------------------------------------------------------------------------------

#include "figtree.h"

#include <stddef.h>  // for definition of NULL
#include <math.h>    // for rounding (floor)
#include <stdio.h>   // for printf

#include "KCenterClustering.h" // provides class for KCenterClustering

#ifndef FIGTREE_NO_ANN
#include "ANN/ANN.h"           // ANN library used for kd-tree in FIGTree code
#endif

#ifndef INT_MAX
#include <limits.h>
#endif

#ifndef DBL_MAX
#include <float.h>
#endif 

// define the upper limit for the truncation number
#define P_UPPER_LIMIT 100

// define MAX and MIN if not yet defined
#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b) )
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b) )
#endif

//------------------------------------------------------------------------------
// change some functions (printf, new, delete) to use Matlab functions
//   instead if we are compiling library for use in matlab mex file.
//------------------------------------------------------------------------------
#ifdef FIGTREE_USE_MATLAB_MEX
#include "mex.h"

// use mexPrintf instead of printf for messages
#undef printf
#define printf mexPrintf

void * operator new (size_t size)
{
  void *p=mxMalloc(size);
  return p;
}

void operator delete (void *p)
{
  mxFree(p); 
}
#endif

//------------------------------------------------------------------------------
// define some macros to use for checking values of input args
// without the macros, all the if statements take up a lot of space
//------------------------------------------------------------------------------
#define FIGTREE_CHECK_POS_NONZERO_DOUBLE( VAR, FCN )                 \
  if( (VAR) <= 0.0 )                                                 \
  {                                                                  \
    printf( #FCN ": Input '" #VAR "' must be a positive number.\n"); \
    return -1;                                                       \
  }

#define FIGTREE_CHECK_POS_DOUBLE( VAR, FCN )                         \
  if( (VAR) < 0.0 )                                                 \
  {                                                                  \
    printf( #FCN ": Input '" #VAR "' must be a positive number.\n"); \
    return -1;                                                       \
  }

#define FIGTREE_CHECK_POS_NONZERO_INT( VAR, FCN )                    \
  if( (VAR) <= 0.0 )                                                 \
  {                                                                  \
    printf( #FCN ": Input '" #VAR "' must be a positive number.\n"); \
    return -1;                                                       \
  }

#define FIGTREE_CHECK_NONNULL_PTR( VAR, FCN )            \
  if( (VAR) == NULL )                                    \
  {                                                      \
  printf( #FCN ": Input pointer '" #VAR "' is NULL.\n"); \
    return -1;                                           \
  }

////////////////////////////////////////////////////////////////////////////////
// Helper functions (their prototpyes do not appear in the header file).
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
// Compute the combinatorial number nchoosek.
// Originally from ImprovedFastGaussTransform.cpp (IFGT source code)
//
// Modified by Vlad Morariu on 2007-06-19.
//------------------------------------------------------------------------------
int nchoosek(int n, int k)
{
  int n_k = n - k;
  
  if (k < n_k)
  {
    k = n_k;
    n_k = n - k;
  }

  int nchsk = 1; 
  for ( int i = 1; i <= n_k; i++)
  {
    nchsk *= (++k);
    nchsk /= i;
  }

  return nchsk;
}

//------------------------------------------------------------------------------
// Compute the combinatorial number nchoosek, using double precision.
// This prevents some overflow issues for large n.
//
// Created by Vlad Morariu on 2008-02-20.
//------------------------------------------------------------------------------
double nchoosek_double(int n, int k)
{
  int n_k = n - k;
  
  if (k < n_k)
  {
    k = n_k;
    n_k = n - k;
  }

  double nchsk = 1; 
  for ( int i = 1; i <= n_k; i++)
  {
    nchsk *= (++k);
    nchsk /= i;
  }

  return nchsk;
}

//------------------------------------------------------------------------------
// This function computes the constants  2^alpha/alpha!.
// Originally compute_constant_series from ImprovedFastGaussTransform.cpp (IFGT 
// source code).
//
// Modified by Vlad Morariu on 2007-06-19.
//------------------------------------------------------------------------------
void computeConstantSeries( int d, int pMaxTotal, int pMax, double * constantSeries )
{ 
  int *heads = new int[d+1];
  int *cinds = new int[pMaxTotal];
  
  for (int i = 0; i < d; i++)
    heads[i] = 0;
  heads[d] = INT_MAX;
  
  cinds[0] = 0;
  constantSeries[0] = 1.0;
  for (int k = 1, t = 1, tail = 1; k < pMax; k++, tail = t)
  {
    for (int i = 0; i < d; i++)
    {
      int head = heads[i];
      heads[i] = t;
      for ( int j = head; j < tail; j++, t++)
      {
        cinds[t] = (j < heads[i+1])? cinds[j] + 1 : 1;
        constantSeries[t] = 2.0 * constantSeries[j];
        constantSeries[t] /= (double) cinds[t];
      }
    }
  }
  
  delete [] cinds;
  delete [] heads; 
}

//------------------------------------------------------------------------------
// This function computes the monomials [(x_i-c_k)/h]^{alpha} and 
// norm([(x_i-c_k)/h])^2.
// Originally compute_source_center_monomials from 
// ImprovedFastGaussTransform.cpp (IFGT source code).
//
// Modified by Vlad Morariu on 2007-06-19.
//------------------------------------------------------------------------------
void computeSourceCenterMonomials( int d, double h, double * dx, 
                                   int p, double * sourceCenterMonomials )
{    
  int * heads = new int[d];

  for (int i = 0; i < d; i++)
  {
    dx[i]=dx[i]/h;
    heads[i] = 0;
  }
    
  sourceCenterMonomials[0] = 1.0;
  for (int k = 1, t = 1, tail = 1; k < p; k++, tail = t)
  {
    for (int i = 0; i < d; i++)
    {
      int head = heads[i];
      heads[i] = t;
      for ( int j = head; j < tail; j++, t++)
        sourceCenterMonomials[t] = dx[i] * sourceCenterMonomials[j];
    }            
  }          

  delete [] heads;
}

//------------------------------------------------------------------------------
// This function computes the monomials [(y_j-c_k)/h]^{alpha}
// Originally compute_target_center_monomials from 
// ImprovedFastGaussTransform.cpp (IFGT source code).
//
// Modified by Vlad Morariu on 2007-06-19.
//------------------------------------------------------------------------------
void computeTargetCenterMonomials( int d, double h, double * dy, 
                                   int pMax, double * targetCenterMonomials )
{    
  int *heads = new int[d];

  for (int i = 0; i < d; i++)
  {
    dy[i] = dy[i]/h;
    heads[i] = 0;
  }
    
  targetCenterMonomials[0] = 1.0;
  for (int k = 1, t = 1, tail = 1; k < pMax; k++, tail = t)
  {
    for (int i = 0; i < d; i++)
    {
      int head = heads[i];
      heads[i] = t;
      for ( int j = head; j < tail; j++, t++)
        targetCenterMonomials[t] = dy[i] * targetCenterMonomials[j];
    }            
  }          

  delete [] heads;
}

//------------------------------------------------------------------------------
// This function computes the coefficients C_k for all clusters.
// Originally compute_C from ImprovedFastGaussTransform.cpp (IFGT source code)
//
// Modified by Vlad Morariu on 2007-06-19.
//------------------------------------------------------------------------------
void computeC( int d, int N, int W, int K, int pMaxTotal, int pMax, 
               double h, int * clusterIndex, double * x, double * q,
               double * clusterCenter, double * C )
{
  double * sourceCenterMonomials = new double[pMaxTotal];
  double * constantSeries = new double[pMaxTotal];
  double hSquare = h*h;
  double * dx = new double[d];

  for (int i = 0; i < W*K*pMaxTotal; i++)
  {
    C[i] = 0.0;
  }

  for(int i = 0; i < N; i++)
  {
    int k = clusterIndex[i];
    int sourceBase = i*d;
    int centerBase = k*d;
    double sourceCenterDistanceSquare = 0.0;

    for (int j = 0; j < d; j++)
    {
      dx[j] = (x[sourceBase+j] - clusterCenter[centerBase+j]);
      sourceCenterDistanceSquare += (dx[j]*dx[j]);
    }
  
    computeSourceCenterMonomials( d, h, dx, pMax, sourceCenterMonomials );    
    
    for(int w = 0; w < W; w++ )
    {
      double f = q[N*w + i]*exp(-sourceCenterDistanceSquare/hSquare);
      for(int alpha = 0; alpha < pMaxTotal; alpha++)
      {
          C[(K*w + k)*pMaxTotal + alpha] += (f*sourceCenterMonomials[alpha]);
      }
    }  
  }

  computeConstantSeries( d, pMaxTotal, pMax, constantSeries );

  for(int w = 0; w < W; w++)
  {   
    for(int k = 0; k < K; k++)
    {
      for(int alpha = 0; alpha < pMaxTotal; alpha++)
      {
        C[(K*w + k)*pMaxTotal + alpha] *= constantSeries[alpha];
      } 
    }
  }

  delete [] sourceCenterMonomials;
  delete [] constantSeries;
  delete [] dx;
}


///////////////////////////////////////////////////////////////////////////////////
// Implementations of functions whose prototypes appear in figtree.h.
///////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
// This function provides a simple interface for evaluation of gauss transforms
// in several ways (direct, direct with approximate nearest-neighbors
// structure on sources, ifgt, and ifgt with approximate nearest-neighbors) and
// using different parameter selection methods (can assume uniform distribution
// or use the actual distribution in estimating parameters).
//
// Created by Vlad Morariu on 2006-11-03.
// Modified by Vlad Morariu on 2007-06-21.
//------------------------------------------------------------------------------
int figtree( int d, int N, int M, int W, double * x, double h, 
          double * q, double * y, double epsilon, double * g,
          int evalMethod, int paramMethod, 
          int verbose, int forceK )
{
  int ret = 0;

  // for FIGTREE_EVAL_DIRECT and FIGTREE_EVAL_DIRECT_TREE, we don't need to compute
  //   parameters, so we just run the fcns directly, once for each set of weights
  if( evalMethod == FIGTREE_EVAL_DIRECT )
  {
    for( int i = 0; i < W; i++ )
      figtreeEvaluateDirect( d, N, M, x, h, q+i*N, y, g+i*M );
  }

  if( evalMethod == FIGTREE_EVAL_DIRECT_TREE )
  {
    for( int i = 0; i < W; i++ )
      figtreeEvaluateDirectTree( d, N, M, x, h, q+i*N, y, epsilon, g+i*M );
  }

  // for FIGTREE_EVAL_IFGT and FIGTREE_EVAL_IFGT_TREE, we must first compute
  //   parameters
  if( evalMethod == FIGTREE_EVAL_IFGT || evalMethod == FIGTREE_EVAL_IFGT_TREE )
  {
    int kLimit = N, K, pMax, kMax;
    double r, maxRange=0;

    // calculate R, the diameter of the hypercube that encompasses sources and targets
    double * mins = new double[d];
    double * maxs = new double[d];
    figtreeCalcMinMax( d, N, x, mins, maxs, 0 );
    figtreeCalcMinMax( d, M, y, mins, maxs, 1 );
    figtreeCalcMaxRange( d, mins, maxs, &maxRange );
    delete [] mins;
    delete [] maxs;

    // choose parameters for IFGT
    if( paramMethod == FIGTREE_PARAM_NON_UNIFORM )
      ret = figtreeChooseParametersNonUniform( d, N, x, h, epsilon, kLimit, maxRange, &kMax, &pMax, &r );
    else
      ret = figtreeChooseParametersUniform( d, h, epsilon, kLimit, maxRange, &kMax, &pMax, &r );
    if( ret < 0 )
    {
      printf("figtree: figtreeChooseParameters%sUniform() failed.\n", 
             ((paramMethod == FIGTREE_PARAM_NON_UNIFORM) ? "Non" : ""));
      return ret;
    }

    verbose && printf("figtreeChooseParameters%sUniform() chose p=%i, k=%i.\n", 
                      ((paramMethod == FIGTREE_PARAM_NON_UNIFORM) ? "Non" : ""), pMax, kMax );

    // if forceK is set to nonzero value, then we set K to that value.
    if (forceK > 0)
      K = forceK;

    // do k-center clustering
    double rx;
    int    * clusterIndex = new int[N];
    int    * numPoints    = new int[kMax]; 
    double * clusterCenters = new double[d*kMax];
    double * clusterRadii = new double[kMax];

    ret = figtreeKCenterClustering( d, N, x, kMax, &K, &rx, clusterIndex, 
                                 clusterCenters, numPoints, clusterRadii );
    if( ret >= 0 )
    {
      // choose truncation number again now that clustering is done
      ret = figtreeChooseTruncationNumber( d, h, epsilon, rx, maxRange, &pMax );
      if( ret >= 0 )
      {
        // evaluate IFGT
        verbose && printf( "Eval IFGT(h= %f, pMax= %i, K= %i, r= %f, rx= %e, epsilon= %f)\n", 
                           h, pMax, K, r, rx, epsilon);
        if( evalMethod == FIGTREE_EVAL_IFGT_TREE  )
        {
          ret = figtreeEvaluateIfgtTree( d, N, M, W, x, h, q, y, pMax, K, clusterIndex, 
                                          clusterCenters, clusterRadii, r, epsilon, g );
        }
        else // evalMethod == FIGTREE_EVAL_IFGT
        {
          ret = figtreeEvaluateIfgt( d, N, M, W, x, h, q, y, pMax, K, clusterIndex, 
                                       clusterCenters, clusterRadii, r, epsilon, g );
        }

        if( ret < 0 )
        {
          printf("figtree: figtreeEvaluateIfgt%s() failed.\n", 
                 ((evalMethod == FIGTREE_EVAL_IFGT_TREE) ? "Tree" : ""));
        }
      }
      else // if figtreeChooseTruncationNumber fails...
      {
        printf("figtree: figtreeChooseTruncationNumber() failed.\n");
      }
    } // if figtreeKCenterClustering fails...
    else
    {
      printf("figtree: figtreeKCenterClustering() failed.\n");
    }

    delete [] clusterIndex;
    delete [] clusterCenters;
    delete [] numPoints;
    delete [] clusterRadii;
  }

  return ret;
}

//------------------------------------------------------------------------------
// Chooses minimum truncation number that satisfies desired error, given
//   the maximum radius of any cluster (rx).
//
// Originally constructor from 
// ImprovedFastGaussTransformChooseTruncationNumber.cpp (IFGT source code) by 
// Vikas C. Raykar.
//
// Modified by Vlad Morariu on 2007-06-20
// Modified by Vlad Morariu on 2007-10-03 - pass R (the max distance between
//     any source and target) as argument instead of assuming that data fits
//     in unit hypercube
//------------------------------------------------------------------------------
int figtreeChooseTruncationNumber( int d, double h, double epsilon, 
                                double rx, double maxRange, int * pMax )
{
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeChooseTruncationNumber );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( h, figtreeChooseTruncationNumber );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( epsilon, figtreeChooseTruncationNumber );
  FIGTREE_CHECK_POS_DOUBLE( rx, figtreeChooseTruncationNumber );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( maxRange, figtreeChooseTruncationNumber );
  FIGTREE_CHECK_NONNULL_PTR( pMax, figtreeChooseTruncationNumber );

  double R = maxRange*sqrt((double)d);
  double hSquare = h*h;
  double r = MIN(R, h*sqrt(log(1/epsilon)));
  double rxSquare = rx*rx;
  
  double error = 1;
  double temp = 1;
  int p = 0;
  while((error > epsilon) & (p <= P_UPPER_LIMIT)){
    p++;
    double b = MIN(((rx + sqrt((rxSquare) + (2*p*hSquare)))/2), rx + r);
    double c = rx - b;
    temp = temp*(((2*rx*b)/hSquare)/p);
    error = temp*(exp(-(c*c)/hSquare));      
  }  

  if( pMax != NULL )
    *pMax = p;    

  return 0;
}

//------------------------------------------------------------------------------
// Parameter selection for the Improved Fast Gauss Transform (IFGT).
//
// Implementation based on:
//
// Fast computation of sums of Gaussians in high dimensions. 
// Vikas C. Raykar, C. Yang, R. Duraiswami, and N. Gumerov,
// CS-TR-4767, Department of computer science,
// University of Maryland, Collegepark.
//
// Originally constructor from ImprovedFastGaussTransformChooseParameters.cpp 
//   by Vikas C. Raykar. (IFGT source code)
//
// Modified by Vlad Morariu on 2007-06-20
// Modified by Vlad Morariu on 2007-10-03 - pass R (the max distance between
//     any source and target) as argument instead of assuming that data fits
//     in unit hypercube
//------------------------------------------------------------------------------
int figtreeChooseParametersUniform( int d, double h, double epsilon, 
                                 int kLimit, double maxRange, int * K, int * pMax, double * r )
{
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeChooseParametersUniform );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( h, figtreeChooseParametersUniform );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( maxRange, figtreeChooseParametersUniform );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( epsilon, figtreeChooseParametersUniform );
  FIGTREE_CHECK_POS_NONZERO_INT( kLimit, figtreeChooseParametersUniform );

  double R = maxRange*sqrt((double)d);
  double hSquare = h*h;
  double complexityMin = DBL_MAX;

  // These variables will hold the values that will then be returned.
  // We use temporary variables in case caller does not care about a variable 
  // and passes a NULL pointer.
  int kTemp = 1;
  int pMaxTemp = P_UPPER_LIMIT + 1;
  double rTemp = MIN(R,h*sqrt(log(1/epsilon)));
  
  for(int i = 0; i < kLimit; i++)
  {
    double rx = maxRange*pow((double)i + 1,-1.0/(double)d);
    //double rx = pow((double)i + 1,-1.0/(double)d); // assumes unit hypercube
    double rxSquare = rx*rx;
    double n = MIN(i + 1, pow(rTemp/rx,(double)d));
    double error = epsilon + 1;
    double temp = 1;
    int p = 0;
    while((error > epsilon) & (p <= P_UPPER_LIMIT))
    {
      p++;
      double b = MIN(((rx + sqrt((rxSquare) + (2*p*hSquare)))/2), rx + rTemp);
      double c = rx - b;
      temp = temp*(((2*rx*b)/hSquare)/p);
      error = temp*(exp(-(c*c)/hSquare));      
    }  
    double complexity = (i + 1) + log((double)i + 1) + ((1 + n)*nchoosek_double(p - 1 + d, d));

     //printf("d=%d r=%f K=%d rx=%f n=%f p=%d terms=%d c=%f\n",d,r,i+1,rx,n,p,nchoosek(p-1+d,d),complexity);
    
    if (complexity < complexityMin )
    {
      complexityMin = complexity;
      kTemp = i + 1;
      pMaxTemp = p;
    }
  }

  // added this to catch case where desired error is never reached.
  // The best thing is to have as many clusters and terms in the taylor
  // series as possible (which will give lowest error)
  if( kTemp == 1 && pMaxTemp == P_UPPER_LIMIT + 1 )
  {
    kTemp = kLimit;
  }

  // set output variables to computed values
  if( K != NULL )
    *K = kTemp;
  if( pMax != NULL )
    *pMax = pMaxTemp;
  if( r != NULL )
    *r = rTemp;

  return 0;
}

//------------------------------------------------------------------------------
// Parameter selection scheme that does not assume uniform distribution.
// In cases where sources are not uniformly distribution, this can lead to 
// very large performance increases(observed as much as 10x speedup)
// because as the number of clusters increases, the max radius of any cluster
// decreases MUCH faster than it would if the sources were uniformly distributed.
// This function is based on ImprovedFastGaussTransformChooseParameters.cpp from
// the IFGT source code, by Vikas C. Raykar.  
//
// Initially created by Vlad Morariu on 2007-01-24.
// Modified by Vlad Morariu on 2007-06-20.
// Modified by Vlad Morariu on 2007-10-03 - pass R (the max distance between
//     any source and target) as argument instead of assuming that data fits
//     in unit hypercube
// Modified by Vlad Morariu on 2008-02-21 - allow rx to be zero.  In some cases,
//     if there isn't a center at each source pt to give rx=0, an excessively
//     large pMax is needed, and it is faster to just have a center at each pt.
//------------------------------------------------------------------------------
int figtreeChooseParametersNonUniform( int d, int N, double * x, 
                                    double h, double epsilon, int kLimit, double maxRange,
                                    int * K, int * pMax, double * r )
{
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeChooseParametersNonUniform );
  FIGTREE_CHECK_POS_NONZERO_INT( N, figtreeChooseParametersNonUniform );
  FIGTREE_CHECK_NONNULL_PTR( x, figtreeChooseParametersNonUniform );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( h, figtreeChooseParametersNonUniform );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( epsilon, figtreeChooseParametersNonUniform );
  FIGTREE_CHECK_POS_NONZERO_INT( kLimit, figtreeChooseParametersNonUniform );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( maxRange, figtreeChooseParametersNonUniform );

  // allocate temporary memory, and set some variables
  int * pClusterIndex = new int[N];
  KCenterClustering * kcc = new KCenterClustering( d, N, x, pClusterIndex, kLimit );

  double R = maxRange*sqrt((double)d);
  double hSquare = h*h;
  double complexityMin = DBL_MAX;

  double rTemp = MIN(R,h*sqrt(log(1/epsilon)));
  int kTemp = 1;
  int pMaxTemp = P_UPPER_LIMIT + 1;

  int numClusters;
  double rx;
  
  // Vlad 01/24/07 - add first cluster and get rx
  kcc->ClusterIncrement( &numClusters, &rx ); 

  // evaluate complexity for increasing values of K
  for(int i = 0; i < kLimit; i++)
  {
    //rx=pow((double)i+1,-1.0/(double)d); // Vlad 01/24/07 - not needed since we're using the real rx
    double rxSquare = rx*rx;
    double n = MIN(i + 1, pow(rTemp/rx,(double)d));
    double error = 1;
    double temp = 1;
    int p = 0;
    while((error > epsilon) & (p <= P_UPPER_LIMIT))
    {
      p++;
      double b = MIN(((rx + sqrt((rxSquare) + (2*p*hSquare)))/2), rx + rTemp);
      double c = rx - b;
      temp = temp*(((2*rx*b)/hSquare)/p);
      error = temp*(exp(-(c*c)/hSquare));      
    }  
    double complexity = (i + 1) + log((double)i + 1) + ((1 + n)*nchoosek_double(p - 1 + d, d));
    //printf("d=%d r=%f K=%d rx=%f n=%f p=%d terms=%d c=%f\n",d,r,i+1,rx,n,p,nchoosek(p-1+d,d),complexity);
    
    if ( (complexity < complexityMin) && (error <= epsilon))
    {
      complexityMin = complexity;
      kTemp = i + 1;
      pMaxTemp = p;    
    }
    
    // try to guess if we have gone past the minimum (the complexity function
    // zigzags as we try different number of clusters, but if it goes up enough,
    // we'll assume we've passed the global minimum).
    // Also stop if truncation number is only 1 or if the max number of unique 
    // clusters are reached (rx = 0).
    if( (p == 1) || (rx <= 0) || (complexity > 1.5*complexityMin ) )
      break;    

    // add another cluster center, and get new max cluster radius
    double rxOld = rx;
    kcc->ClusterIncrement( &numClusters, &rx );
  }  

  // copy results 
  if( K != NULL )
    *K = kTemp;
  if( pMax != NULL )
    *pMax = pMaxTemp;
  if( r != NULL )
    *r = rTemp;
  delete [] pClusterIndex;
  delete kcc;

  return 0;
}

//------------------------------------------------------------------------------
// This function is an interface to the C++ KCenterClustering class from the
// original IFGT library.
//
// Created by Vlad Morariu 2007-06-19.
//------------------------------------------------------------------------------
int figtreeKCenterClustering( int d, int N, double * x, int kMax, int * K,
                           double * rx, int * clusterIndex, double * clusterCenters, 
                           int * numPoints, double * clusterRadii )
{
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeKCenterClustering );
  FIGTREE_CHECK_POS_NONZERO_INT( N, figtreeKCenterClustering );
  FIGTREE_CHECK_NONNULL_PTR( x, figtreeKCenterClustering );
  FIGTREE_CHECK_POS_NONZERO_INT( kMax, figtreeKCenterClustering );
  FIGTREE_CHECK_NONNULL_PTR( K, figtreeKCenterClustering );
  FIGTREE_CHECK_NONNULL_PTR( rx, figtreeKCenterClustering );
  FIGTREE_CHECK_NONNULL_PTR( clusterIndex, figtreeKCenterClustering );
  FIGTREE_CHECK_NONNULL_PTR( clusterCenters, figtreeKCenterClustering );
  FIGTREE_CHECK_NONNULL_PTR( numPoints, figtreeKCenterClustering );
  FIGTREE_CHECK_NONNULL_PTR( clusterRadii, figtreeKCenterClustering );

  //k-center clustering
  KCenterClustering* pKCC = new KCenterClustering( d, N, x, clusterIndex, kMax );
  *K = pKCC->Cluster();
  if( rx != NULL )
    *rx = pKCC->MaxClusterRadius;
  pKCC->ComputeClusterCenters(*K, clusterCenters, numPoints, clusterRadii);

  delete pKCC;
  return 0;
}

//------------------------------------------------------------------------------
// Actual function to evaluate the exact Gauss Transform directly.
// Originally Evaluate() from GaussTransform.cpp, written by Vikas C. Raykar.
//
// Modified by Vlad Morariu on 2007-06-19.
//------------------------------------------------------------------------------
int figtreeEvaluateDirect( int d, int N, int M, double * x, double h, 
                        double * q, double * y, double * g )
{
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeEvaluateDirect );
  FIGTREE_CHECK_POS_NONZERO_INT( N, figtreeEvaluateDirect );
  FIGTREE_CHECK_POS_NONZERO_INT( M, figtreeEvaluateDirect );
  FIGTREE_CHECK_NONNULL_PTR( x, figtreeEvaluateDirect );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( h, figtreeEvaluateDirect );
  FIGTREE_CHECK_NONNULL_PTR( q, figtreeEvaluateDirect );
  FIGTREE_CHECK_NONNULL_PTR( y, figtreeEvaluateDirect );
  FIGTREE_CHECK_NONNULL_PTR( g, figtreeEvaluateDirect );

  // evaluate
  double hSquare = h*h;
  for(int j = 0; j < M; j++)
  {
    g[j] = 0.0;
    for(int i = 0; i < N; i++)
    {
      double norm = 0.0;
      for (int k = 0; k < d; k++)
      {
        double temp = x[(d*i) + k] - y[(d*j) + k];
        norm = norm + (temp*temp);
      }
      g[j] = g[j] + (q[i]*exp(-norm/hSquare));
    }
  }

  return 0;
}

//------------------------------------------------------------------------------
// This function approximates Gauss Transform.
// Originally constructor, Evaluate(), and destructor from 
// ImprovedFastGaussTransform.cpp (IFGT source code).
//
// Modified by Vlad Morariu on 2007-06-19.
//------------------------------------------------------------------------------
int figtreeEvaluateIfgt( int d, int N, int M, int W, double * x, 
                           double h, double * q, double * y, 
                           int pMax, int K, int * clusterIndex, 
                           double * clusterCenter, double * clusterRadii,
                           double r, double epsilon, double * g )
{
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeEvaluateIfgt );
  FIGTREE_CHECK_POS_NONZERO_INT( N, figtreeEvaluateIfgt );
  FIGTREE_CHECK_POS_NONZERO_INT( M, figtreeEvaluateIfgt );
  FIGTREE_CHECK_POS_NONZERO_INT( W, figtreeEvaluateIfgt );
  FIGTREE_CHECK_NONNULL_PTR( x, figtreeEvaluateIfgt );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( h, figtreeEvaluateIfgt );
  FIGTREE_CHECK_NONNULL_PTR( q, figtreeEvaluateIfgt );
  FIGTREE_CHECK_NONNULL_PTR( y, figtreeEvaluateIfgt );
  FIGTREE_CHECK_POS_NONZERO_INT( pMax, figtreeEvaluateIfgt );
  FIGTREE_CHECK_POS_NONZERO_INT( K, figtreeEvaluateIfgt );
  FIGTREE_CHECK_NONNULL_PTR( clusterIndex, figtreeEvaluateIfgt );
  FIGTREE_CHECK_NONNULL_PTR( clusterCenter, figtreeEvaluateIfgt );
  FIGTREE_CHECK_NONNULL_PTR( clusterRadii, figtreeEvaluateIfgt );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( r, figtreeEvaluateIfgt );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( epsilon, figtreeEvaluateIfgt );
  FIGTREE_CHECK_NONNULL_PTR( g, figtreeEvaluateIfgt );

  //Memory allocation
  int pMaxTotal = nchoosek(pMax - 1 + d, d);
  double hSquare=h*h;
  double * targetCenterMonomials = new double[pMaxTotal];
  double * dy = new double[d];
  double * C = new double[W*K*pMaxTotal];
  double * ry = new double[K];
  double * rySquare = new double[K];

  for(int i = 0; i < K; i++)
  {
    ry[i] = r + clusterRadii[i];
    rySquare[i] = ry[i]*ry[i];
  }   

  //////////////////////////////////////////////////////////////////////////////
  // Evaluate  
  //////////////////////////////////////////////////////////////////////////////
  computeC( d, N, W, K, pMaxTotal, pMax, h, clusterIndex, x, q, clusterCenter, C );  

  for(int j = 0; j < M; j++)
  {
    for( int w = 0; w < W; w++ )
    {
      g[M*w + j] = 0.0;
    }

    int targetBase = j*d;        
    for(int k = 0; k < K; k++)
    {
      int centerBase = k*d;
      double targetCenterDistanceSquare = 0.0;
      for(int i = 0; i < d; i++)
      {
        dy[i] = y[targetBase + i] - clusterCenter[centerBase + i];
        targetCenterDistanceSquare += dy[i]*dy[i];
        if(targetCenterDistanceSquare > rySquare[k]) break;
      }

      if(targetCenterDistanceSquare <= rySquare[k])
      {
        computeTargetCenterMonomials( d, h, dy, pMax, targetCenterMonomials );
        double f=exp(-targetCenterDistanceSquare/hSquare);
        for(int w = 0; w < W; w++ )
        {        
          for(int alpha = 0; alpha < pMaxTotal; alpha++)
          {
            g[M*w + j] += (C[(K*w + k)*pMaxTotal + alpha]*f*targetCenterMonomials[alpha]);
          }
        }                      
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Release memory
  //////////////////////////////////////////////////////////////////////////////
  delete [] rySquare;
  delete [] ry;
  delete [] C;
  delete [] dy;
  delete [] targetCenterMonomials;

  return 0;
}

//------------------------------------------------------------------------------
// This function approximates Gauss Transform using Approximate Nearest 
// Neighbors.
// Originally constructor, Evaluate(), and destructor from 
// ImprovedFastGaussTransform.cpp of FIGTree code, by Vikas C. Raykar.
//
// Modified by Vlad Morariu on 2007-06-19.
//------------------------------------------------------------------------------
int figtreeEvaluateIfgtTree( int d, int N, int M, int W, double * x, 
                              double h, double * q, double * y, 
                              int pMax, int K, int * clusterIndex, 
                              double * clusterCenter, double * clusterRadii,
                              double r, double epsilon, double * g )
{
#ifdef FIGTREE_NO_ANN
  printf("This code was not compiled with support for ANN.  Please recompile\n");
  printf("with 'FIGTREE_NO_ANN' not defined to enable ANN support.\n");
  return -1;
#else
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_POS_NONZERO_INT( N, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_POS_NONZERO_INT( M, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_POS_NONZERO_INT( W, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_NONNULL_PTR( x, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( h, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_NONNULL_PTR( q, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_NONNULL_PTR( y, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_POS_NONZERO_INT( pMax, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_POS_NONZERO_INT( K, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_NONNULL_PTR( clusterIndex, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_NONNULL_PTR( clusterCenter, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_NONNULL_PTR( clusterRadii, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( r, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( epsilon, figtreeEvaluateIfgtTree );
  FIGTREE_CHECK_NONNULL_PTR( g, figtreeEvaluateIfgtTree );

  //Memory allocation
  int pMaxTotal = nchoosek(pMax-1+d,d);
  double * targetCenterMonomials = new double[pMaxTotal];
  double * dy = new double[d];
  double * C = new double[W*K*pMaxTotal]; 
  double hSquare = h*h;

  //Find the maximum cluster radius
  double pcr_max = clusterRadii[0];
  for(int i = 0; i < K; i++)
  {
    if (clusterRadii[i] > pcr_max)
    {
      pcr_max = clusterRadii[i];
    }
  }
  double rSquare=(r+pcr_max)*(r+pcr_max);

  //Allocate storage using ANN procedures 
  ANNpointArray dataPts = annAllocPts(K,d);     // allocate data points
  ANNidxArray   nnIdx   = new ANNidx[K];        // allocate near neigh indices
  ANNdistArray  dists = new ANNdist[K];         // allocate near neighbor dists
  
  // Copy the cluster centers to the ANN data structure
  for (int k = 0; k < K; k++)
  {
    for ( int j = 0; j < d; j++ )
      dataPts[k][j]= clusterCenter[k*d + j];
  }

  // build search structure
  ANNkd_tree * kdTree = new ANNkd_tree(              
                                        dataPts,  // the data points
                                        K,        // number of points
                                        d,        // dimension of space
                                        1,
                                        ANN_KD_SUGGEST);

  ////////////////////////////////////////////////////////////////////
  // Evaluate
  ////////////////////////////////////////////////////////////////////
  computeC( d, N, W, K, pMaxTotal, pMax, h, clusterIndex, x, q, clusterCenter, C );  

  for(int j = 0; j < M; j++)
  {
    for( int w = 0; w < W; w++ )
    {
      g[M*w+j]=0.0;
    }

    int targetBase=j*d;        
    
    ANNpoint queryPt=&(y[targetBase]);

    int NN = kdTree->annkFRSearch(           // search
                                   queryPt,  // query point
                                   rSquare,  // squared radius
                                   0,        // number of near neighbors
                                   NULL,     // nearest neighbors (returned)
                                   NULL,     // distance (returned)
                                   0.0 );
    if (NN>0)
    {
      kdTree->annkFRSearch(           // search
                            queryPt,  // query point
                            rSquare,  // squared radius
                            NN,       // number of near neighbors
                            nnIdx,    // nearest neighbors (returned)
                            dists,    // distance (returned)
                            0.0);

      for(int l = 0; l < NN; l++)
      {
        int k = nnIdx[l];
        int centerBase = k*d;
        double  targetCenterDistanceSquare = dists[l];
        for(int i = 0; i < d; i++)
        {
          dy[i] = y[targetBase + i] - clusterCenter[centerBase + i];
        }
        computeTargetCenterMonomials( d, h, dy, pMax, targetCenterMonomials );
        double e = exp(-targetCenterDistanceSquare/hSquare);
        for(int w = 0; w < W; w++ )
        {        
          for(int alpha = 0; alpha < pMaxTotal; alpha++)
          {
            g[M*w + j] += (C[(K*w+k)*pMaxTotal + alpha]*e*targetCenterMonomials[alpha]);
          }    
        }
      }
    }
  }

  ////////////////////////////////////////////////////////////////////
  // Release Memory
  ////////////////////////////////////////////////////////////////////
  delete [] targetCenterMonomials;
  delete [] dy;
  delete [] C;

  annDeallocPts(dataPts);
  delete [] nnIdx;              
  delete [] dists;
  delete kdTree;
  annClose();        

  return 0;
#endif
}

//------------------------------------------------------------------------------
// Gauss Transform computed using the ANN library.
// Given a specified epsilon, the code computes the Gauss transform by summing 
// the sources only within a certain radius--whose contribution is at least 
// epsilon. The neighbors are found using the ANN library.
// http://www.cs.umd.edu/~mount/ANN/.
// Originally constructor, Evaluate(), and destructor from GaussTransformTree.cpp 
// of FIGTree code, by Vikas C. Raykar.
//
// Modified by Vlad Morariu on 2007-06-20
//------------------------------------------------------------------------------
int figtreeEvaluateDirectTree( int d, int N, int M, double * x, double h, 
                           double * q, double * y, double epsilon, double * g )
{
#ifdef FIGTREE_NO_ANN
  printf("This code was not compiled with support for ANN.  Please recompile\n");
  printf("with compiler flag 'FIGTREE_NO_ANN' not set to enable ANN support.\n");
  return -1;
#else
  // check input arguments 
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeEvaluateDirectTree );
  FIGTREE_CHECK_POS_NONZERO_INT( N, figtreeEvaluateDirectTree );
  FIGTREE_CHECK_POS_NONZERO_INT( M, figtreeEvaluateDirectTree );
  FIGTREE_CHECK_NONNULL_PTR( x, figtreeEvaluateDirectTree );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( h, figtreeEvaluateDirectTree );
  FIGTREE_CHECK_NONNULL_PTR( q, figtreeEvaluateDirectTree );
  FIGTREE_CHECK_NONNULL_PTR( y, figtreeEvaluateDirectTree );
  FIGTREE_CHECK_POS_NONZERO_DOUBLE( epsilon, figtreeEvaluateDirectTree );
  FIGTREE_CHECK_NONNULL_PTR( g, figtreeEvaluateDirectTree );

  double hSquare = h*h;
  double epsANN = 0.0;

  // Compute the cutoff radius
  double r = h*sqrt(log(1/epsilon));
  double rSquare=r*r;

  // Allocate storage using ANN procedures 
  ANNpointArray dataPts = annAllocPts(N,d);  // allocate data points
  ANNidxArray   nnIdx   = new ANNidx[N];     // allocate near neigh indices
  ANNdistArray  dists   = new ANNdist[N];    // allocate near neighbor dists
  
  // Copy the source points to the ANN data structure
  for (int i = 0; i < N; i++)
  {
    for ( int j = 0; j < d; j++ )
      dataPts[i][j]= x[i*d + j];
  }

  // build search structure
  ANNkd_tree * kdTree = new ANNkd_tree(              
                                        dataPts, // the data points
                                        N,       // number of points
                                        d,       // dimension of space
                                        1,
                                        ANN_KD_SUGGEST );

  ///////////////////////////////////////////////////////////////////////
  // Evaluate
  ///////////////////////////////////////////////////////////////////////
  for(int j = 0; j < M; j++)
  {
    g[j] = 0.0;  
    int targetBase = j*d;          
    ANNpoint queryPt = &(y[targetBase]);

    int NN = kdTree->annkFRSearch( // fixed radius search
                                   queryPt,  // query point
                                   rSquare,  // squared radius
                                   0,        // number of near neighbors
                                   NULL,     // nearest neighbors (returned)
                                   NULL,     // distance (returned)
                                   epsANN );
    if (NN > 0)
    {
      kdTree->annkFRSearch(           // fixed radius search
                            queryPt,  // query point
                            rSquare,  // squared radius
                            NN,       // number of near neighbors
                            nnIdx,    // nearest neighbors (returned)
                            dists,    // distance (returned)
                            epsANN );

      for(int l = 0; l < NN; l++)
      {
        int i = nnIdx[l];
        double sourceTargetDistanceSquare = dists[l];
        g[j] += (q[i]*exp(-sourceTargetDistanceSquare/hSquare));
      }
    } 
  }

  //////////////////////////////////////////////////////////////////////////////
  // Free memory
  //////////////////////////////////////////////////////////////////////////////
  annDeallocPts(dataPts);
  delete [] nnIdx;              
  delete [] dists;
  delete kdTree;
  annClose();   

  return 0;
#endif
}

int figtreeCalcMinMax( int d, int n, double * x, double * mins, double * maxs, int update )
{
  // check input arguments
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeCalcMinMax );
  FIGTREE_CHECK_POS_NONZERO_INT( n, figtreeCalcMinMax );
  FIGTREE_CHECK_NONNULL_PTR( x, figtreeCalcMinMax );
  FIGTREE_CHECK_NONNULL_PTR( mins, figtreeCalcMinMax );
  FIGTREE_CHECK_NONNULL_PTR( maxs, figtreeCalcMinMax );

  // use first sample values as current min and max if we're not updating 
  //   some previously computed min and max values.  
  if( update != 1 && n > 0 )
  {
    for( int i = 0; i < d; i++ )
    {
      mins[i] = x[i];
      maxs[i] = x[i];
    }
  }

  // go through each sample in x and update mins and maxs for each dimension
  for( int i = 0; i < n; i++ )
  {
    for( int j = 0; j < d; j++ )
    {
      mins[j] = MIN( mins[j], x[i*d+j] );
      maxs[j] = MAX( maxs[j], x[i*d+j] );
    }
  }

  return 0;
}

int figtreeCalcMaxRange( double d, double * mins, double * maxs, double * maxRange )
{
  FIGTREE_CHECK_POS_NONZERO_INT( d, figtreeCalcMaxRange );
  FIGTREE_CHECK_NONNULL_PTR( mins, figtreeCalcMaxRange );
  FIGTREE_CHECK_NONNULL_PTR( maxs, figtreeCalcMaxRange );
  FIGTREE_CHECK_NONNULL_PTR( maxRange, figtreeCalcMaxRange );

  double maxRangeTemp = maxs[0] - mins[0];
  for( int i = 0; i < d; i++ )
    maxRangeTemp = MAX( maxRangeTemp, maxs[i] - mins[i] );
  *maxRange = maxRangeTemp;
  return 0;
}
