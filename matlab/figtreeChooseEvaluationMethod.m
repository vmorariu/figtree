function [bestMethod, estCost] = figtreeChooseEvaluationMethod( x, h, W, y, epsilon, ifgtParamMethod, verbose )
%
%  This function chooses the best evaluation method given the source and target points as well as the parameters
%  such as bandwidth and epsilon.  This function analyzes the source and target distributions and estimates
%  the cost of each method.  In the interest of saving time, some estimates are not done if it is obvious
%  that the corresponding methods are slower, in which case an entry of -1 will be returned for est_cost.
%
%
%% Input
%
%    * x --> d x N matrix of N source points in d dimensions.
%    * h --> the source scale or bandwidth.
%    * W --> number of weights that will be used for each source point. 
%            This is useful if one needs multiple gauss transforms that have
%            the same sources, targets, and bandwidth, but different
%            weights/strengths (q). By computing coefficients for all W weight
%            sets at once, we avoid duplicating much of the overhead.  However,
%            more memory is needed to store a set of coefficients for each set
%            of weights.
%    * y --> d x M matrix of M target points in d dimensions.
%    * epsilon --> desired error
%    * ifgtParamMethod --> (optional, default = 1) the method to use for determining parameters.
%            Can be 0 (UNIFORM) or 1 (NON_UNIFORM).  The first assumes that points are
%            uniformly distributed, and the second adapts to whatever distribution is given when
%            choosing IFGT parameters K (number of clusters) and p_max (max truncation number).
%            NOTE: NON_UNIFORM is recommended if data is not uniformly distributed.
%    * verbose --> (optional, default = 0) if nonzero, prints parameters chosen for evaluation
%
%% Ouput
%
%    * bestMethod --> the evaluation method to use in evaluating gauss 
%            transform. Can be 0, 1, 2, or 3 for DIRECT, IFGT, DIRECT_TREE,
%            IFGT_TREE evaluation, respectively.
%    * estCost --> the estimated cost of each method.  If in the interest of time the cost 
%            was not estimate for some method that was determined to be too slow, a -1 is returned.
%
%% Signature
%
% Author: Vlad I Morariu 
% E-Mail: morariu@cs.umd.edu
% Date:  2008-12-05
%
%% See also
%
%  figtree
%