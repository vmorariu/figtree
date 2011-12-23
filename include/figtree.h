// File: figtree.h
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

//------------------------------------------------------------------------------
// The code was written by Vlad Morariu, Vikas Raykar, and Changjiang Yang 
// and is copyrighted under the Lesser GPL: 
//
// Copyright (C) 2008 Vlad Morariu, Vikas Raykar, and Changjiang Yang 
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
#ifndef FIGTREE_H
#define FIGTREE_H

#ifdef WIN32
  //----------------------------------------------------------------------
  //  To compile the code into a windows DLL, you must define the 
  //  symbol FIGTREE_DLL_EXPORTS. 
  // 
  //  To compile the code statically into a windows executable 
  //    (i.e. not using a separate DLL) define FIGTREE_DLL_STATIC.
  //----------------------------------------------------------------------
  #ifdef FIGTREE_STATIC
    #define FIGTREE_DLL_API  // since FIGTREE_STATIC is defined, code is statically 
                          // linked and no exports or imports are needed
  #else
    #ifdef FIGTREE_DLL_EXPORTS
      #define FIGTREE_DLL_API __declspec(dllexport)
    #else
      #define FIGTREE_DLL_API __declspec(dllimport)
    #endif
  #endif
#else
  //----------------------------------------------------------------------
  // FIGTREE_DLL_API is ignored for all other systems
  //----------------------------------------------------------------------
  #define FIGTREE_DLL_API
#endif

// we want to export C functions if using DLL
#if defined(__cplusplus) && !defined(FIGTREE_STATIC)
extern "C" {
#endif

//------------------------------------------------------------------------------
// Some useful constants for the figtree() function call
#define FIGTREE_EVAL_DIRECT        0  // direct evaluation of gauss transform
#define FIGTREE_EVAL_IFGT          1  // truncated taylor series evaluation
#define FIGTREE_EVAL_DIRECT_TREE   2  // direct evaluation, with tree on sources
#define FIGTREE_EVAL_IFGT_TREE     3  // truncated taylor series, with tree on 
                                      //   cluster centers

#define FIGTREE_PARAM_UNIFORM      0  // estimate params assuming sources are 
                                      //   uniformly distributed
#define FIGTREE_PARAM_NON_UNIFORM  1  // estimate params by using actual source 
                                      //   distribution (runs k-center clustering 
                                      //   twice, but speedup during evaluation
                                      //   more than makes up for it in general)

// Note: All matrix pointers are assumed to point to a contiguous one 
//   dimensional array containing the entries of the matrx in row major format.
//   Thus, for an M x N matrix and a pointer ptr to its data, 
//   ptr[0] ... ptr[N-1] contains the first row of the matrix, and so on. 

// Evaluates gauss transform in one shot (chooses parameters and does clustering
// if necessary, and evaluates). This is the main function since it calls all 
// the functions provided below.  If more flexibility is required than this 
// function allows, use the functions provided below (eg. it is possible that 
// clustering needs to be done only once for a set of sources and then multiple 
// evaluations are done using same parameters and clustering).  
//
// Input
//    * d --> data dimensionality.
//    * N --> number of source points.
//    * M --> number of target points.
//    * W --> number of weights that will be used for each source point. 
//            This is useful if one needs multiple gauss transforms that have
//            the same sources, targets, and bandwidth, but different
//            weights/strengths (q). By computing coefficients for all W weight
//            sets at once, we avoid duplicating much of the overhead.  However,
//            more memory is needed to store a set of coefficients for each set
//            of weights.
//    * x --> N x d matrix of N source points in d dimensions.
//    * h --> the source scale or bandwidth.
//    * q --> W x N matrix of the source strengths.
//    * y --> M x d matrix of M target points in d dimensions.
//    * epsilon --> desired error
//    * evalMethod --> the evaluation method to use in evaluating gauss 
//            transform. Can be FIGTREE_EVAL_[DIRECT,IFGT,DIRECT_TREE,
//            IFGT_TREE], defined above. epsilon is needed for all but 
//            DIRECT method.  Parameter selection is done only in the IFGT
//            or IFGT_TREE case.  
//    * paramMethod --> the method to use for determining parameters.
//            Can be FIGTREE_PARAM_UNIFORM or FIGTREE_PARAM_NON_UNIFORM.  
//    * verbose --> if nonzero, prints parameters chosen for evaluation
//    * forceK --> if zero, the number of clusters(K) is determined by parameter
//            selection.  If nonzero, forceK overrides the K chosen by
//            parameter selection.
//
// Output
//    * g --> W x M vector of the Gauss Transform evaluated at each target point. 
//            The ith row is the result of the transform using the ith set of 
//            weights.
FIGTREE_DLL_API 
int figtree( int d, int N, int M, int W, double * x, double h, 
             double * q, double * y, double epsilon, double * g,
             int evalMethod = FIGTREE_EVAL_IFGT,
             int paramMethod = FIGTREE_PARAM_NON_UNIFORM, 
             int verbose = 0, int forceK = 0 );

// Given the maximum cluster radius, this function computes the maximum 
// truncation number that guarantees results within the desired error bound.
//Input
//    * d --> dimension of the points.
//    * h --> the source bandwidth.
//    * epsilon --> the desired error.
//    * rx --> maximum cluster radius
//    * maxRange --> max dimension range.  The range along a dimension is 
//          the difference between the max and min values that can ever
//          occur along that dimension.  The max dimension range is the 
//          maximum range among all dimensions.  For example, if all 
//          points lie in unit hypercube, then maxRange = 1.
//
//Output
//    * pMax --> maximum truncation number for the Taylor series.
FIGTREE_DLL_API 
int figtreeChooseTruncationNumber( int d, double h, double epsilon, 
                                   double rx, double maxRange, int * pMax );

// Chooses parameters for IFGT and FIGTree by assuming that sources are
//   uniformly distributed in a unit cube.
//Input
//    * d --> dimension of the points.
//    * h --> the source bandwidth.
//    * epsilon --> the desired error.
//    * kLimit --> upper limit on the number of clusters, kLimit.
//    * maxRange --> max dimension range.  The range along a dimension is 
//          the difference between the max and min values that can ever
//          occur along that dimension.  The max dimension range is the 
//          maximum range among all dimensions.  For example, if all 
//          points lie in unit hypercube, then maxRange = 1.
//
//Note : [ Use roughly kLimit=round(40*sqrt(d)/h) ]
//Output
//    * K --> number of clusters.
//    * pMax --> maximum truncation number for the Taylor series.
//    * r --> source cutoff radius.
FIGTREE_DLL_API 
int figtreeChooseParametersUniform( int d, double h, double epsilon, 
                                    int kLimit, double maxRange,
                                    int * K, int * pMax, double * r );

// Chooses parameters for IFGT and FIGTree without assumption that sources are
// uniformly distributed.  In this case, k-center clustering is done as part
// of the parameter selection so that the radius of each cluster is not 
// estimated but computed directly for the sources.  This results in an 
// additional k-center clustering operation, but the speedup when sources are
// not uniformly distributed can be large because
// optimal parameters are much more accurately estimated.  Even when sources are
// uniformly distributed, the slowdown is often small compared to 
// figtreeChooseParametersUniform.
// 
// Input
//    * d --> dimension of the points.
//    * N --> number of source points.
//    * x --> N x d matrix of N source points in d dimensions.
//    * h --> the source bandwidth.
//    * epsilon --> the desired error.
//    * kLimit --> upper limit on the number of clusters, K.
// Note : Use kLimit=N to allow for optimal param estimation.
//    * maxRange --> max dimension range.  The range along a dimension is 
//          the difference between the max and min values that can ever
//          occur along that dimension.  The max dimension range is the 
//          maximum range among all dimensions.  For example, if all 
//          points lie in unit hypercube, then maxRange = 1.
// Output
//    * K --> number of clusters.
//    * pMax --> maximum truncation number for the Taylor series.
//    * r --> source cutoff radius.
FIGTREE_DLL_API 
int figtreeChooseParametersNonUniform( int d, int N, double * x, 
                                       double h, double epsilon, int kLimit, double maxRange,
                                       int * K, int * pMax, double * r );

// Gonzalez's farthest-point clustering algorithm. O(N log K) version.
//
// Input
//    * d --> data dimensionality.
//    * N --> number of source points.
//    * x --> N x d matrix of N source points in d dimensions 
//       (in one contiguous array, row major format where each row is a point).
//    * kMax --> maximum number of clusters.
//
// Output
//    * K --> actual number of clusters (less than kMax if duplicate pts exist)
//    * rx --> maximum radius of the clusters (rx).
//    * clusterIndex --> vector of length N where the i th element is the 
//                cluster number to which the i th point belongs. 
//                ClusterIndex[i] varies between 0 to K-1.
//    * clusterCenters --> K x d matrix of K cluster centers 
//                (contiguous 1-d array, row major format).
//    * numPoints --> number of points in each cluster.
//    * clusterRadii --> radius of each cluster.
FIGTREE_DLL_API 
int figtreeKCenterClustering( int d, int N, double * x, int kMax, int * K,
                              double * rx, int * clusterIndex, double * clusterCenters, 
                              int * numPoints, double * clusterRadii );

// Computes exact gauss transform (within machine precision) using direct 
// evaluation. Provided for time/error comparison.
//
// Input
//    * d --> data dimensionality.
//    * N --> number of source points.
//    * M --> number of target points.
//    * x --> N x d matrix of N source points in d dimensions.
//    * h --> the source scale or bandwidth.
//    * q --> 1 x N or N x 1 vector of the source strengths.
//    * y --> M x d matrix of M target points in d dimensions.
//
// Output
//    * g --> 1 x M vector of the Gauss Transform evaluated at each target
//            point. 
FIGTREE_DLL_API 
int figtreeEvaluateDirect( int d, int N, int M, double * x, double h, 
                           double * q, double * y, double * g );

// Computes an approximation to Gauss Transform.  Implementation based on:
// Fast computation of sums of Gaussians in high dimensions. Vikas C. Raykar, 
//   C. Yang, R. Duraiswami, and N. Gumerov, CS-TR-4767, Department of computer
//   science, University of Maryland, College Park.
//
// Input
//    * d --> data dimensionality.
//    * N --> number of source points.
//    * M --> number of target points.
//    * W --> number of weights that will be used for each source point. 
//            This really does multiple transforms, with different weights each
//            time but with same sources and targets.  This saves a lot of time
//            since most of the work is not duplicated.  However, it requires 
//            more memory to store the coefficients for each set of weights.
//    * x --> N x d matrix of N source points in d dimensions.
//    * h --> the source scale or bandwidth.
//    * q --> W x N vector of the source strengths.
//    * y --> M x d matrix of M target points in d dimensions.
//    * pMax --> maximum truncation number for the Taylor series.
//    * K --> the number of clusters.
//    * clusterIndex --> N x 1 vector the i th element is the cluster number 
//            to which the i th point belongs. [ ClusterIndex[i] varies between
//            0 to K-1. ]
//    * clusterCenter --> K x d matrix of K cluster centers.
//    * clusterRadii  --> K x 1 matrix of the radius of each cluster.
//    * r --> cutoff radius
//    * epsilon --> desired error
//
//Output
//    * g --> W x M vector of the Gauss Transform evaluated at each target
//            point. Each row q is the result of the transform using the qth set
//            of weights.
FIGTREE_DLL_API 
int figtreeEvaluateIfgt( int d, int N, int M, int W, double * x, 
                         double h, double * q, double * y, 
                         int pMax, int K, int * clusterIndex, 
                         double * clusterCenter, double * clusterRadii,
                         double r, double epsilon, double * g );

// Computes an approximation to Gauss Transform using approximate 
// nearest-neighbors. Same as figtreeEvaluateIfgt() but uses Approximate 
// Nearest-Neighbors library to find source clusters which influence each target
// (part of FIGTree). Parameters for this function can be computed using 
// figtreeChooseParameters[Non]Uniform and figtreeKCenterClustering. 
//
// Input
//    * d --> data dimensionality.
//    * N --> number of source points.
//    * M --> number of target points.
//    * W --> number of weights that will be used for each source point.
//            This really does multiple transforms, with different weights each
//            time but with same sources and targets.  This saves a lot of time
//            since most of the work is not duplicated.  However, it requires
//            more memory to store the coefficients for each set of weights.
//    * x --> N x d matrix of N source points in d dimensions.
//    * h --> the source scale or bandwidth.
//    * q --> W x N vector of the source strengths.
//    * y --> M x d matrix of M target points in d dimensions.
//    * pMax --> maximum truncation number for the Taylor series.
//    * K --> the number of clusters.
//    * clusterIndex --> N x 1 vector the i th element is the cluster number 
//            to which the i th point belongs. [ ClusterIndex[i] varies between
//            0 to K-1. ]
//    * clusterCenter --> K x d matrix of K cluster centers.
//    * clusterRadii  --> K x 1 matrix of the radius of each cluster.
//    * r --> cutoff radius
//    * epsilon --> desired error
//
// Output
//    * g --> W x M vector of the Gauss Transform evaluated at each target
//            point.  Each row of g is the result of the transform using one set
//            of weights.
FIGTREE_DLL_API 
int figtreeEvaluateIfgtTree( int d, int N, int M, int W, double * x, 
                             double h, double * q, double * y, 
                             int pMax, int K, int * clusterIndex, 
                             double * clusterCenter, double * clusterRadii,
                             double r, double epsilon, double * g );

// Computes an approximation to Gauss Transform using approximate 
// nearest-neighbors.  Direct method (no taylor expansion is done), with tree 
// directly on samples (part of FIGTree). Requires Approximate 
// Nearest-Neighbor(ANN) library.
//
// Input
//    * d --> data dimensionality.
//    * N --> number of source points.
//    * M --> number of target points.
//    * x --> N x d matrix of N source points in d dimensions.
//    * h --> the source scale or bandwidth.
//    * q --> 1 x N vector of the source strengths.
//    * y --> M x d matrix of M target points in d dimensions.
//    * epsilon --> desired error
//
// Output
//    * g --> 1 x M vector of the Gauss Transform evaluated at each target
//            point.
FIGTREE_DLL_API 
int figtreeEvaluateDirectTree( int d, int N, int M, double * x, double h, 
                               double * q, double * y, double epsilon, double * g );

// Computes min and max values along each dimension.  Used to determine
//   the size of the hypercube (and the max distance R that any two pts 
//   can be from each other).
//
// Input
//    * d --> data dimensionality.
//    * n --> number of source points.
//    * x --> n x d matrix of n source points in d dimensions.
//    * mins --> d x 1 vector of minimum values; input values ignored if update == 0
//    * maxs --> d x 1 vector of maximum values; input values ignored if update == 0
//    * update --> if set to 1, then max[i] will contain 
//          max(max of values of all samples along dimension i, max[i] input value), and 
//          similarly for min[i].
//
// Output
//    * mins --> d x 1 vector of minimum values
//    * maxs --> d x 1 vector of maximum values
FIGTREE_DLL_API
int figtreeCalcMinMax( int d, int n, double * x, double * mins, double * maxs, int update=0 );

FIGTREE_DLL_API
int figtreeCalcMaxRange( double d, double * mins, double * maxs, double * R );

#if defined(__cplusplus) && !defined(FIGTREE_STATIC)
} // extern "C"
#endif

#endif //FIGTREE_H