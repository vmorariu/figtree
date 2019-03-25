# Author: Vlad Morariu
# Created: 2012-01-23
import numpy as np
cimport numpy as np


cdef extern from "figtree.h":
    cdef int c_figtree "figtree" (
       int d, int N, int M, int W, float * x, float h, float * q,
       float * y, float epsilon, float * g, int evalMethod,
       int ifgtParamMethod, int ifgtTruncMethod, int verbose)
    cdef int c_k_centers "figtreeKCenterClustering" (
       int d, int N, float * x, int kMax, int * K,
  float * rx, int * clusterIndex, float * clusterCenters, 
  int * numPoints, float * clusterRadii)


def figtree(np.ndarray[np.float64_t, ndim=2, mode='c'] X, float h,
            np.ndarray[np.float64_t, ndim=2, mode='c'] Q,
            np.ndarray[np.float64_t, ndim=2, mode='c'] Y,
            float epsilon, int evalMethod=4,
            int ifgtParamMethod=1, int ifgtTruncMethod=2, int verbose=0):
    """Wrapper for the figtree C function."""
    N, d = X.shape[0], X.shape[1]
    W = int(Q.size / N)
    M, dY = Y.shape[0], Y.shape[1]
    assert((d == dY) and (W*N == Q.size) and (Q.size == N or (
        (Q.ndim == 2) and (Q.shape[0] == W) and (Q.shape[1] == N))))
    assert(epsilon > 0)
    cdef np.ndarray G = np.zeros((W, M), dtype='float64')
    cdef int ret = c_figtree(d, N, M, W, <float*>X.data, h, <float*>Q.data,
                        <float*>Y.data, epsilon, <float*>G.data, evalMethod,
                        ifgtParamMethod, ifgtTruncMethod, verbose)
    if ret >= 0:
        return G
    return None


def k_centers(np.ndarray[np.float64_t, ndim=2, mode='c'] X, int K):
    """Wrapper for the K-center clustering function used by figtree.
       indexes, clusters, num_points, radii = k_centers(X, K)
       Input
          * X --> N x d matrix of N source points in d dimensions 
                 (in one contiguous array, row major format where each row is a point).
          * K --> maximum number of clusters.
       Output
          * indexes --> vector of length N where the i th element is the 
                        cluster number to which the i th point belongs. 
                        indexes[i] varies between 0 to K-1.
          * centers --> K x d matrix of K cluster centers 
                        (contiguous 1-d array, row major format).
          * num_points --> number of points in each cluster.
          * radii --> radius of each cluster.
    """    
    N, d = X.shape[0], X.shape[1]
    cdef int K_out = 0
    cdef float rx = 0
    cdef np.ndarray[np.int32_t, ndim=1] clusterIndex = np.zeros((N,), dtype='int32')
    cdef np.ndarray[np.float64_t, ndim=2] clusterCenters = np.zeros((K, d), dtype='float64')
    cdef np.ndarray[np.int32_t, ndim=1] numPoints = np.zeros((K,), dtype='int32')
    cdef np.ndarray[np.float64_t, ndim=1] clusterRadii = np.zeros((K,), dtype='float64')
    cdef int ret = c_k_centers(d, N, <float*>X.data, K, &K_out, &rx, <int*>clusterIndex.data,
                               <float*>clusterCenters.data, <int*>numPoints.data,
                               <float*>clusterRadii.data)
    if ret < 0:
        raise ValueError('k_centers returned failure code...')
    return (clusterIndex, clusterCenters[:K_out, :], numPoints[:K_out], clusterRadii[:K_out])
    
