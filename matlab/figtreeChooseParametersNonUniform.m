function [K, pMax, r] = figtreeChooseParametersNonUniform(d, h, epsilon, Klimit, maxRange)
%
%  Choose the parameters for the Improved Fast Gauss Transform.
%
%  Chooses parameters for IFGT and FIGTree without assumption that sources are
%  uniformly distributed.  In this case, k-center clustering is done as part
%  of the parameter selection so that the radius of each cluster is not 
%  estimated but computed directly for the sources.  This results in an 
%  additional k-center clustering operation, but the speedup when sources are
%  not uniformly distributed can be large (as much as 10 times faster) because
%  optimal parameters are much more accurately estimated.  Even when sources are
%  uniformly distributed, the slowdown is often small compared to 
%  figtreeChooseParametersUniform (however, there are cases where this is slower).
%
%
%% Input
%
%    * d --> dimension of the points.
%    * N --> number of source points.
%    * x --> d x N matrix of N source points in d dimensions.
%    * h --> the source bandwidth.
%    * epsilon --> the desired error.
%    * kLimit --> upper limit on the number of clusters, K.
% Note : Use kLimit=N to allow for optimal param estimation.
%    * maxRange --> max dimension range.  The range along a dimension is 
%          the difference between the max and min values that can ever
%          occur along that dimension.  The max dimension range is the 
%          maximum range among all dimensions.  For example, if all 
%          points lie in unit hypercube, then maxRange = 1.%
%% Ouput
%
%    * K -->number of clusters.
%    * pMax --> maximum truncation number for the Taylor series.
%    * r --> source cutoff radius.
%
%% Signature
%
% Author: Vlad I Morariu 
%         (original implementation by Vikas C. Raykar, vikas@cs.umd.edu)
% E-Mail: morariu@cs.umd.edu
% Date:  2007-06-26
%
