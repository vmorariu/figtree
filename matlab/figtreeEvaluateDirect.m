function [ g ] = figtreeEvaluateDirect(d, N, M, x, h, q, y)
%
%     Direct computation of the Gauss Transform.
%
%     $$G(y_j)=\sum_{i=1}^{N} q_i e^{\|x_i-y_j\|^2/h^2},\:\:j=1...M$$
%
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