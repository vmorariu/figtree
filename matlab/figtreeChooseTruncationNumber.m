function [pMax] = figtreeChooseTruncationNumber(d, h, epsilon, rx)
%
%  Choose the parameters for the Improved Fast Gauss Transform.
%
%  Given the maximum cluster radius, this function computes the maximum 
%  truncation number that guarantees results within the desired error bound.
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
%    * rx --> maximum cluster radius
%    * maxRange --> max dimension range.  The range along a dimension is 
%          the difference between the max and min values that can ever
%          occur along that dimension.  The max dimension range is the 
%          maximum range among all dimensions.  For example, if all 
%          points lie in unit hypercube, then maxRange = 1.%
%
%% Ouput
%
%    * pMax --> maximum truncation number for the Taylor series.%
%
%% Signature
%
% Author: Vlad I Morariu 
%         (original implementation by Vikas C. Raykar, vikas@cs.umd.edu)
% E-Mail: morariu@cs.umd.edu
% Date:  2007-06-26
%
