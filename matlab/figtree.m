function [G] = figtree( x, h, q, y, epsilon, evalMethod, ifgtParamMethod, ifgtTruncMethod, verbose )
%
%     Fast computation of the Gauss Transform.
%
%     Computes and approximation $$\hat{G}(y_j)$$ to $$G(y_j)=\sum_{i=1}^{N} q_i
%     e^{\|x_i-y_j\|^2/h^2},\:\:j=1...M$$  such that
%     $$|\hat{G}(y_j)-G(y_j)| \leq Q \epsilon$$ , where
%     $$Q=\sum_{i=1}^{N}q_i$$.
%
%     C++ Implementation of algorithms described in the following publication:
%
%     Vlad I. Morariu, Balaji Vasan Srinivasan, Vikas C. Raykar, Ramani Duraiswami, and Larry S. Davis. 
%     Automatic online tuning for fast Gaussian summation. Advances in Neural Information Processing 
%     Systems (NIPS), 2008.
%
%% Input
%
%
%    * x --> d x N matrix of N source points in d dimensions.
%    * h --> the source scale or bandwidth.
%    * q --> N x W matrix of the source strengths (one set of weights is a column vector).
%            Here W is the number of weights that will be used for each source point. 
%            This is useful if one needs multiple gauss transforms that have
%            the same sources, targets, and bandwidth, but different
%            weights/strengths (q). By computing coefficients for all W weight
%            sets at once, we avoid duplicating much of the overhead.  However,
%            more memory is needed to store a set of coefficients for each set
%            of weights.
%    * y --> d x M matrix of M target points in d dimensions.
%    * epsilon --> desired error
%    * evalMethod --> (optional, default = 4) the evaluation method to use in evaluating gauss 
%            transform. Can be the following:
%                 0 - DIRECT
%                 1 - IFGT
%                 2 - DIRECT_TREE
%                 3 - IFGT_TREE
%                 4 - AUTO (automatically chooses between 0-3 based on data and parameters)
%            epsilon is needed for all but 'direct' method.  Parameter selection is done only 
%            in the ifgt-based methods. Setting to 4 (AUTO) allows figtree to estimate which 
%            method is likely to be faster based on source and target distributions (NOTE: the 
%            best method to run changes with bandwidth, epsilon and source/target distributions).
%    * ifgtParamMethod --> (optional, default = 1) the method to use for determining parameters.
%            Can be 0 (UNIFORM) or 1 (NON_UNIFORM).  The first assumes that points are
%            uniformly distributed, and the second adapts to whatever distribution is given when
%            choosing IFGT parameters K (number of clusters) and p_max (max truncation number).
%            NOTE: NON_UNIFORM is recommended if data is not uniformly distributed.
%    * ifgtTruncMethod -> (optional, default = 2) the method to use for computing truncation numbers.
%                 0 - MAX, use max truncation number for all points
%                 1 - POINT, point-wise truncation
%                 2 - CLUSTER, cluster-wise truncation
%    * verbose --> (optional, default = 0) if nonzero, prints parameters chosen for evaluation
%
%% Ouput
%
%    * g --> M x W vector of the Gauss Transform evaluated at each target point. 
%            The ith column is the result of the transform using the ith set of 
%            weights (ith column of q).
%
%% Signature
%
% Author: Vlad I Morariu 
% E-Mail: morariu@cs.umd.edu
% Date:  2008-12-05
%