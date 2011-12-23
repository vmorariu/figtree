function [g]=figtreeEvaluateIfgt(d, N, M, W, h, q, y, pMax, K, clusterIndex, clusterCenter, clusterRadii, r, epsilon)
%  Computes Gauss Transform by evaluating truncated Taylor series (IFGT).
%
%  Computes and approximation $$\hat{G}(y_j)$$ to $$G(y_j)=\sum_{i=1}^{N} q_i
%  e^{\|x_i-y_j\|^2/h^2},\:\:j=1...M$$  such that
%  $$|\hat{G}(y_j)-G(y_j)| \leq Q \epsilon$$ , where
%  $$Q=\sum_{i=1}^{N}q_i$$.
%
%  Implementation based on:
%  Fast computation of sums of Gaussians in high dimensions. Vikas C. Raykar, 
%  C. Yang, R. Duraiswami, and N. Gumerov, CS-TR-4767, Department of computer
%  science, University of Maryland, College Park.%
% 
%% Input
%
%    * d --> data dimensionality.
%    * N --> number of source points.
%    * M --> number of target points.
%    * W --> number of weights that will be used for each source point. 
%            This really does multiple transforms, with different weights each
%            time but with same sources and targets.  This saves a lot of time
%            since most of the work is not duplicated.  However, it requires 
%            more memory to store the coefficients for each set of weights.
%    * x --> d x N matrix of N source points in d dimensions.
%    * h --> the source scale or bandwidth.
%    * q --> N x W vector of the source strengths.
%    * y --> d x M matrix of M target points in d dimensions.
%    * pMax --> maximum truncation number for the Taylor series.
%    * K --> the number of clusters.
%    * clusterIndex --> 1 x N vector the i th element is the cluster number 
%            to which the i th point belongs. [ ClusterIndex[i] varies between
%            0 to K-1. ]
%    * clusterCenter --> d x K matrix of K cluster centers.
%    * clusterRadii  --> 1 x K matrix of the radius of each cluster.
%    * r --> cutoff radius
%    * epsilon --> desired error
%
%% Ouput
%
%    * g --> W x M vector of the Gauss Transform evaluated at each target
%            point. Each row q is the result of the transform using the qth set
%            of weights.
%
%% Signature
%
% Author: Vlad I Morariu 
%         (original implementation by Vikas C. Raykar, vikas@cs.umd.edu)
% E-Mail: morariu@cs.umd.edu
% Date:  2007-06-26
%
