% Script the demonstrate the use of FIGTree.  Note, only 
% (FIGTREE-DIRECTORY)/matlab/figtree.dll is needed for the figtree() function to 
% function.
% 
% Created by Vlad Morariu on 2008-12-03.
clear all;
close all;
clear functions;
clc;

% the MEX binaries for the FIGTree library reside in ../matlab
% we assume that the script is run from the figtree/samples/ directory,
% otherwise the correct path for the mex files must be set.
addpath('../matlab');

disp('---------------------------------------------');
disp(sprintf('Example to demonstrate the use of FIGTree'));
disp('---------------------------------------------');

% the distribution
disp(sprintf('Sources uniformly distributed\nTargets uniformly distributed'));

% the number of sources
N = 10000;
disp(sprintf('Number of source points N=%d',N));

% the number of targets
M = 10000;
disp(sprintf('Number of target points M=%d',M));

% the desired error
epsilon = 1e-2;     
disp(sprintf('Target error epsilon=%e',epsilon));

% the data dimensionality
d = 3;
disp(sprintf('Dimensionality d=%d',d));

% the bandwidth
h = .02;
disp(sprintf('Bandwidth h=%f',h));

% the source points
% d x N matrix of N source points in d dimensions.
X = rand(d,N);

% the target points
% d x M matrix of M source points in d dimensions.
Y = rand(d,M);

% the source weights
% 1 x N row vector
q = rand(1,N);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate gauss transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call figtree, using default parameters.  The default is to choose the
% fastest method automatically by predicting estimated evaluation cost of
% each method (evalMethod=4), to choose ifgt parameters (if an ifgt method
% is chosen) without assuming a uniform distribution (ifgtParamMethod=1),
% and to use cluster-wise truncations for ifgt (ifgtTruncMethod=2).  These
% default parameters allow for the fastest evaluation and best parameter
% tuning to the actual distributions.
%
% If we wanted to, we could include these default parameters as follows:
%
% evalMethod = 4; ifgtParamMethod = 1; ifgtTruncMethod = 2; verbose=0;
% [ g_auto ] = figtree( X, h, q, Y, epsilon, evalMethod, ifgtParamMethod, ifgtTruncMethod, verbose );
%
% However, the function call below is equivalent to the two lines shown
% above, and is simpler to write out.
fprintf('Running figtree (automatic method + param selection) ... ');
tic;
[ g_auto ] = figtree( X, h, q, Y, epsilon );
t_auto = toc;
fprintf('done.\n');

%
% To evaluate error and speedup, we will also call figtree with
% evalMethod=0 so that direct evaluation is performed.
%
fprintf('Running exact method, for error and time comparison ... ');
tic;
evalMethod=0;
[ g_exact ] = figtree( X, h, q, Y, epsilon, evalMethod );
t_exact = toc;
fprintf('done.\n');

disp('---------------------------------------------');
fprintf('exact  : %3.2e seconds\n', t_exact );
fprintf('auto   : %3.2e seconds, speedup = %3.2f, error = %3.2e\n', ... 
    t_auto, t_exact/t_auto, max(abs(g_exact-g_auto))/sum(q) );