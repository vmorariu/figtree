function [G] = figtree(d, N, M, W, x, h, q, y, epsilon)
%
%     Fast computation of the Gauss Transform.
%
%     Computes and approximation $$\hat{G}(y_j)$$ to $$G(y_j)=\sum_{i=1}^{N} q_i
%     e^{\|x_i-y_j\|^2/h^2},\:\:j=1...M$$  such that
%     $$|\hat{G}(y_j)-G(y_j)| \leq Q \epsilon$$ , where
%     $$Q=\sum_{i=1}^{N}q_i$$.
%
%     C++ Implementation.
%
%
%% Input
%
%    * d --> data dimensionality.
%    * N --> number of source points.
%    * M --> number of target points.
%    * W --> number of weights that will be used for each source point. 
%            This is useful if one needs multiple gauss transforms that have
%            the same sources, targets, and bandwidth, but different
%            weights/strengths (q). By computing coefficients for all W weight
%            sets at once, we avoid duplicating much of the overhead.  However,
%            more memory is needed to store a set of coefficients for each set
%            of weights.
%    * x --> d x N matrix of N source points in d dimensions.
%    * h --> the source scale or bandwidth.
%    * q --> N x W matrix of the source strengths (one set of weights is a column vector).
%    * y --> d x M matrix of M target points in d dimensions.
%    * epsilon --> desired error
%    * evalMethod --> (optional, default = 1) the evaluation method to use in evaluating gauss 
%            transform. Can be 0, 1, 2, or 3 for DIRECT, TRUNCATED, DIRECT_ANN,
%            TRUNCATED_ANN evaluation, respectively. epsilon is needed for all but 
%            DIRECT method.  Parameter selection is done only in the TRUNCATED
%            or TRUNCATED_ANN case.  
%    * paramMethod --> (optional, default = 1) the method to use for determining parameters.
%            Can be 0 (IFGT_PARAM_UNIFORM) or 1 (IFGT_PARAM_NON_UNIFORM).  
%    * verbose --> (optional, default = 0) if nonzero, prints parameters chosen for evaluation
%    * forceK --> (optional, default = 0) if zero, the number of clusters(K) is determined by parameter
%            selection.  If nonzero, forceK overrides the K chosen by
%            parameter selection.
%
%% Ouput
%
%    * g --> W x M vector of the Gauss Transform evaluated at each target point. 
%            The ith column is the result of the transform using the ith set of 
%            weights (ith column of q).
%
%% Signature
%
% Author: Vlad I Morariu 
%         (original implementation by Vikas C. Raykar, vikas@cs.umd.edu)
% E-Mail: morariu@cs.umd.edu
% Date:  2007-06-26
%
%% See also
%
%  IfgtChooseParametersUniform, IfgtChooseParametersNonUniform,
%  IfgtChooseTruncationNumber, IfgtEvaluateDirect, IfgtEvaluateDirectAnn,
%  IfgtEvaluateTruncated, IfgtEvaluateTruncatedAnn, IfgtKCenterClustering
%