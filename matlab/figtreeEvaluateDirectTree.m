function [ g ] = figtreeEvaluateDirectTree(d, N, M, x, h, q, y, epsilon)
%
% Gauss Transform with fast nearest neighbor searching.
%
% Given a specified epsilon, the code computes the Gauss transform by 
% summing the sources only within a certain radius--whos contribution 
% is alteast epsilon.
%
% The neighbors are found using the ANN library
% http://www.cs.umd.edu/~mount/ANN/.
%
% Computes and approximation $$\hat{G}(y_j)$$ to $$G(y_j)=\sum_{i=1}^{N} q_i
% e^{\|x_i-y_j\|^2/h^2},\:\:j=1...M$$  such that
% $$|\hat{G}(y_j)-G(y_j)| \leq Q \epsilon$$ , where
% $$Q=\sum_{i=1}^{N}q_i$$.
%
% This method is extremely fast for small bandwidths and different source
% and target distributions.
%
%% Input
%
% * d --> data dimensionality.
% * N --> number of source points.
% * M --> number of target points.
% * x --> d x N matrix of N source points in d dimensions.
% * h --> source bandwidth or scale.
% * q --> N x 1 or 1 x N vector of the source strengths.
% * y --> d x M matrix of M target points in d dimensions.
% * epsilon --> desired error.
%
%% Ouput
%
% * g --> M x 1 vector of the Gauss Transform evaluated at  each target point.
%         Note that unlike GaussTransform from the original IFGT library, 
%         the output is a column vector.
%
%% Signature
%
% Author: Vlad I Morariu 
%         (original implementation by Vikas C. Raykar, vikas@cs.umd.edu)
% E-Mail: morariu@cs.umd.edu
% Date:  2007-06-26
%