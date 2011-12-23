function [K, pMax, r] = figtreeChooseParametersUniform(d, h, epsilon, Klimit)
%
%  Choose the parameters for the Improved Fast Gauss Transform.
%
%  Chooses parameters for IFGT and FIGTree by assuming that sources are
%  uniformly distributed in a unit cube.
%
%  Implementation based on:
%  
%  Fast computation of sums of Gaussians in high dimensions. 
%  Vikas C. Raykar, C. Yang, R. Duraiswami, and N. Gumerov,
%  CS-TR-4767, Department of computer science,
%  University of Maryland, Collegepark.
%
%% Input
%
%    * d --> dimension of the points.
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
%    * K --> number of clusters.
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
