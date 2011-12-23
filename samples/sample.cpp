//------------------------------------------------------------------------------
// The code was written by Vlad I. Morariu 
// and is copyrighted under the Lesser GPL: 
//
// Copyright (C) 2007 Vlad I. Morariu
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
// The author may be contacted via email at: morariu(at)cs(.)umd(.)edu 
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// File    : sample.cpp
// Purpose : Show example of how figtree() function can be used to evaluate
//           gauss transforms.  
// Author  : Vlad I. Morariu       morariu(at)cs(.)umd(.)edu
// Date    : 2007-06-25
// Modified: 2008-01-23 to add comments and make the sample clearer
// Modified: 2010-05-12 to add the automatic method selection function call
//             example.  
//------------------------------------------------------------------------------

#include <stdio.h>
#include <string.h> // for memset
#include <math.h>   // for abs
#include "figtree.h"

#ifdef WIN32
  // Link to figtree.dll if WIN32.  If not, then we assume that 
  // libfigtree.so or libfigtree.a are linked to from the Makefile.
  // The locations of figtree.dll and ANN.dll must either be in the PATH
  // environment variable, or in the same location from which sample.exe is
  // is executed.
  #pragma comment(lib,"../lib/figtree.lib")
#endif

// This file only shows examples of how the figtree() function can be used.
//
// See source code of figtree() in figtree.cpp for sample code of how 
// the functions figtreeEvaluate*(), figtreeChooseParameters*(), and figtreeKCenterClustering()
// can be used to evaluate gauss transforms.  Calling them directly instead of the 
// wrapper function, figtree(), might be useful if the same sources are always used so 
// the parameter selection and clustering steps do not have to be performed each time
// the gauss transform needs to be evaluated.

int main()
{
  // The dimensionality of each sample vector.
  int d = 7;

  // The number of targets (vectors at which gauss transform is evaluated).
  int M = 10;

  // The number of sources which will be used for the gauss transform. 
  int N = 20;

  // The bandwidth.  NOTE: this is not the same as standard deviation since 
  // the Gauss Transform sums terms exp( -||x_i - y_j||^2 / h^2 ) as opposed
  // to  exp( -||x_i - y_j||^2 / (2*sigma^2) ).  Thus, if sigma is known, 
  // bandwidth can be set to h = sqrt(2)*sigma.
  double h = .8;
  
  // Desired maximum absolute error after normalizing output by sum of weights.
  // If the weights, q_i (see below), add up to 1, then this is will be the 
  // maximum absolute error.
  // The smaller epsilon is, the more accurate the results will be, at the
  // expense of increased computational complexity.
  double epsilon = 1e-2;

  // The source array.  It is a contiguous array, where
  // ( x[i*d], x[i*d+1], ..., x[i*d+d-1] ) is the ith d-dimensional sample.
  // For example, below N = 20 and d = 7, so there are 20 rows, each
  // a 7-dimensional sample.
  double x[] = {0.7165, 0.5113, 0.7764, 0.4893, 0.1859, 0.7006, 0.9827,
                0.8066, 0.7036, 0.4850, 0.1146, 0.6649, 0.3654, 0.1400,
                0.5668, 0.8230, 0.6739, 0.9994, 0.9616, 0.0589, 0.3603,
                0.5485, 0.2618, 0.5973, 0.0493, 0.5711, 0.7009, 0.9623,
                0.7505, 0.7400, 0.4319, 0.6343, 0.8030, 0.0839, 0.9455,
                0.9159, 0.6020, 0.2536, 0.8735, 0.5134, 0.7327, 0.4222,
                0.1959, 0.1059, 0.3923, 0.1003, 0.6930, 0.2069, 0.6094,
                0.1059, 0.0396, 0.2093, 0.9693, 0.1059, 0.3029, 0.3069,
                0.9692, 0.6029, 0.2222, 0.2059, 0.3059, 0.6092, 0.2133,
                0.9614, 0.0721, 0.5534, 0.2920, 0.8580, 0.3358, 0.6802,
                0.2473, 0.3527, 0.1879, 0.4906, 0.4093, 0.4635, 0.6109,
                0.1865, 0.0395, 0.5921, 0.1853, 0.9963, 0.1953, 0.7659,
                0.0534, 0.3567, 0.4983, 0.4344, 0.5625, 0.6166, 0.1133,
                0.8983, 0.7546, 0.7911, 0.8150, 0.6700, 0.2009, 0.2731,
                0.6262, 0.5369, 0.0595, 0.0890, 0.2713, 0.4091, 0.4740,
                0.1332, 0.6926, 0.0009, 0.1532, 0.9632, 0.3521, 0.9692,
                0.9623, 0.3532, 0.7432, 0.0693, 0.2336, 0.6022, 0.2936,
                0.3921, 0.6023, 0.6323, 0.9353, 0.3963, 0.2835, 0.9868,
                0.2362, 0.6682, 0.2026, 0.0263, 0.1632, 0.9164, 0.1153,
                0.9090, 0.5962, 0.3290, 0.4782, 0.5972, 0.1614, 0.8295 };

  // The target array.  It is a contiguous array, where
  // ( y[j*d], y[j*d+1], ..., y[j*d+d-1] ) is the jth d-dimensional sample.
  // For example, below M = 10 and d = 7, so there are 10 rows, each
  // a 7-dimensional sample.
  double y[] = {0.9561, 0.5955, 0.0287, 0.8121, 0.6101, 0.7015, 0.0922,
                0.4249, 0.3756, 0.1662, 0.8332, 0.8386, 0.4516, 0.9566,
                0.1472, 0.8699, 0.7694, 0.4442, 0.6206, 0.9517, 0.6400,
                0.0712, 0.3143, 0.6084, 0.1750, 0.6210, 0.2460, 0.5874,
                0.5061, 0.4648, 0.5414, 0.9423, 0.3418, 0.4018, 0.3077,
                0.4116, 0.2859, 0.3941, 0.5030, 0.7220, 0.3062, 0.1122,
                0.4433, 0.4668, 0.0147, 0.6641, 0.7241, 0.2816, 0.2618,
                0.7085, 0.7839, 0.9862, 0.4733, 0.9028, 0.4511, 0.8045,
                0.8289, 0.1663, 0.3939, 0.5208, 0.7181, 0.5692, 0.4608,
                0.4453, 0.0877, 0.4435, 0.3663, 0.3025, 0.8518, 0.7595 };

  // The weight array.  The ith weight is associated with the ith source sample.
  // To evaluate the Gauss Transform with the same sources and targets, but 
  // different sets of weights, add another row of weights and set W = 2.  
  double q[] = {0.2280, 0.4496, 0.1722, 0.9688, 0.3557, 0.0490, 0.7553, 0.8948, 0.2861, 0.2512, 0.9327, 0.3353, 0.2530, 0.2532, 0.3352, 0.7235, 0.2506, 0.0235, 0.1062, 0.1061, 0.7234, 0.1532};

  // Number of weights.  For each set of weights a different Gauss Transform is computed, 
  // but by giving multiple sets of weights at once some overhead can be shared.
  int W = 1;  // in this case W = 1.

  // allocate array into which the result of the Gauss Transform will be stored for each
  // target sample.  The first M elements will correspond to the Gauss Transform computed
  // with the first set of weights, second M elements will correspond to the G.T. computed
  // with the second set of weights, etc.
  double * g_auto = new double[W*M];
  double * g_sf = new double[W*M];
  double * g_sf_tree = new double[W*M];
  double * g_ifgt_u = new double[W*M];
  double * g_ifgt_tree_u = new double[W*M];
  double * g_ifgt_nu = new double[W*M];
  double * g_ifgt_tree_nu = new double[W*M];
  
  // initialize all output arrays to zero
  memset( g_auto        , 0, sizeof(double)*W*M );
  memset( g_sf          , 0, sizeof(double)*W*M );
  memset( g_sf_tree     , 0, sizeof(double)*W*M );
  memset( g_ifgt_u      , 0, sizeof(double)*W*M );
  memset( g_ifgt_tree_u , 0, sizeof(double)*W*M );
  memset( g_ifgt_nu     , 0, sizeof(double)*W*M );
  memset( g_ifgt_tree_nu, 0, sizeof(double)*W*M );

  // 
  // RECOMMENDED way to call figtree().
  // 
  // Evaluates the Gauss transform using automatic method selection (the automatic
  // method selection function analyzes the inputs -- including source/target
  // points, weights, bandwidth, and error tolerance -- to automatically choose 
  // between FIGTREE_EVAL_DIRECT, FIGTREE_EVAL_DIRECT_TREE, FIGTREE_EVAL_IFGT,
  // FIGTREE_EVAL_IFGT_TREE.
  // This function call makes use of the default parameters for the eval method
  // and param method, and is equivalent to
  // figtree( d, N, M, W, x, h, q, y, epsilon, g_auto, FIGTREE_EVAL_AUTO, FIGTREE_PARAM_NON_UNIFORM ).
  figtree( d, N, M, W, x, h, q, y, epsilon, g_auto );

  //
  // MANUAL EVALUATION METHOD and PARAMETER METHOD selection.  If chosen 
  // incorrectly, this could cause run times to be several orders of 
  // magnitudes longer or to require several orders of magnitude more
  // memory (resulting in crashes).  The recommended way to call figtree
  // is using the automatic method selection, as shown above.
  //
  // evaluate gauss transform using direct (slow) method
  figtree( d, N, M, W, x, h, q, y, epsilon, g_sf, FIGTREE_EVAL_DIRECT );

  // evaluate gauss transform using direct method with approximate nearest neighbors
  figtree( d, N, M, W, x, h, q, y, epsilon, g_sf_tree, FIGTREE_EVAL_DIRECT_TREE );

  // evaluate gauss transform using FIGTREE (truncated series), estimating parameters with and without
  //   the assumption that sources are uniformly distributed
  figtree( d, N, M, W, x, h, q, y, epsilon, g_ifgt_u, FIGTREE_EVAL_IFGT, FIGTREE_PARAM_UNIFORM, 1 );
  figtree( d, N, M, W, x, h, q, y, epsilon, g_ifgt_nu, FIGTREE_EVAL_IFGT, FIGTREE_PARAM_NON_UNIFORM, 1 );

  // evaluate gauss transform using FIGTREE (truncated series), estimating parameters with and without
  //   the assumption that sources are uniformly distributed
  figtree( d, N, M, W, x, h, q, y, epsilon, g_ifgt_tree_u, FIGTREE_EVAL_IFGT_TREE, FIGTREE_PARAM_UNIFORM, 1 );
  figtree( d, N, M, W, x, h, q, y, epsilon, g_ifgt_tree_nu, FIGTREE_EVAL_IFGT_TREE, FIGTREE_PARAM_NON_UNIFORM, 1 );

  // compute absolute error of the Gauss Transform at each target and for all sets of weights.
  for( int i = 0; i < W*M; i++)
  {
    g_auto[i]         = fabs(g_auto[i]        -g_sf[i]);
    g_sf_tree[i]      = fabs(g_sf_tree[i]     -g_sf[i]);
    g_ifgt_u[i]       = fabs(g_ifgt_u[i]      -g_sf[i]);
    g_ifgt_nu[i]      = fabs(g_ifgt_nu[i]     -g_sf[i]);
    g_ifgt_tree_u[i]  = fabs(g_ifgt_tree_u[i] -g_sf[i]);
    g_ifgt_tree_nu[i] = fabs(g_ifgt_tree_nu[i]-g_sf[i]);
  }

  // print out results for all six ways to evaluate
  printf("Results:\n");
  printf("Direct result auto error    sf-tree error ifgt-u error  ifgt-nu error  ifgt-tree-u e ifgt-tree-nu error\n");
  for( int i = 0; i < W*M; i++ )
    printf("%13.4f %13.4e %13.4e %13.4e %13.4e %13.4e %13.4e\n", 
           g_sf[i], g_auto[i], g_sf_tree[i], g_ifgt_u[i], g_ifgt_nu[i], g_ifgt_tree_u[i], g_ifgt_tree_nu[i] );
  printf("\n");

  // deallocate memory
  delete [] g_auto;
  delete [] g_sf;
  delete [] g_sf_tree;
  delete [] g_ifgt_u;
  delete [] g_ifgt_nu;
  delete [] g_ifgt_tree_u;
  delete [] g_ifgt_tree_nu;
  return 0;
}
