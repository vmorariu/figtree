% Script the demonstrate the use of FIGTree.  Note, only 
% (FIGTREE-DIRECTORY)/matlab/figtree.dll is needed for the figtree() function to 
% function.  All other mex files are provided in case users require
% advanced control over parameter selection and clustering processes.  See
% sample2.m for example of how to use the additional mex files provided.
clear all;
close all;
clear functions;
clc;

% the MEX binaries for the FIGTree library reside in ../matlab
% we assume that the script is run from the figtree/samples/ directory,
% otherwise the correct path for the mex files must be set.
addpath('../matlab');

disp('---------------------------------------------');
disp(sprintf('Examples to demonstrate the use of FIGTree'));
disp('---------------------------------------------');

% the number of sources
N = 2000;
disp(sprintf('Number of source points N=%d',N));

% the number of targets
M = 2000;
disp(sprintf('Number of target points M=%d',M));

% the desired error
epsilon = 1e-3;     
disp(sprintf('Target error epsilon=%e',epsilon));

% verbosity
verbose = 0;

ds = 1:6;
hs = [.01 .02 .04 .08 .16 .32 .64];
ntrials = 1;

timings_direct = zeros( ntrials, length(hs), length(ds) );
timings_direct_tree = zeros( ntrials, length(hs), length(ds) );
timings_ifgt = zeros( ntrials, length(hs), length(ds) );
timings_ifgt_nu = zeros( ntrials, length(hs), length(ds) );
timings_ifgt_tree = zeros( ntrials, length(hs), length(ds) );
timings_ifgt_tree_nu = zeros( ntrials, length(hs), length(ds) );

for k = 1:length(ds)
    d = ds(k);
    for j = 1:length(hs)
        h = hs(j);
        for i = 1:ntrials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Generate sources and targets
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('\n\n');
            disp('----------------------------------------------------------------------');
            disp(sprintf('Sources drawn from multiple gaussians\nTargets uniformly distributed'));

            % the data dimensionality
            disp(sprintf('Dimensionality d=%d',d));

            % the bandwidth
            %h = rand(1);
            disp(sprintf('Bandwidth h=%f',h));
            disp('----------------------------------------------------------------------');

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

            % direct method (exact up to machine precision)
            tic; 
            [ g_direct           ] = figtree( d, N, M, 1, X, h, q, Y, epsilon, 0, 0, verbose ); 
            t_direct = toc;
            fprintf('direct       : %3.2e seconds\n', t_direct );

            % direct evaluation with approximate nearest neighbor structure 
            % (parameter selection not relevant, so it can be anything) - 
            % implemented as in the FIGTree paper (Raykar et al, submitted to 
            % Joural of Machine Learning).
            tic;
            [ g_direct_tree       ] = figtree( d, N, M, 1, X, h, q, Y, epsilon, 2, 0, verbose );
            t_direct_tree = toc;
            fprintf('direct-tree  : %3.2e seconds, speedup = %3.2f\n', t_direct_tree, t_direct/t_direct_tree );

            % truncated series, uniform parameter selection - implemented exactly
            % as described in IFGT paper (Raykar et al, 2005).
            tic; 
            [ g_ifgt_u      ] = figtree( d, N, M, 1, X, h, q, Y, epsilon, 1, 0, verbose ); 
            t_ifgt_u = toc;
            fprintf('ifgt-u       : %3.2e seconds, speedup = %3.2f\n', t_ifgt_u, t_direct/t_ifgt_u );

            % truncated series, with nonuniform parameter selection - the
            % default method when calling function without last two arguments
            % that select the method: figtree( d, N, M, 1, X, h, q, Y, epsilon );
            % works well in general
            tic; 
            [ g_ifgt_nu     ] = figtree( d, N, M, 1, X, h, q, Y, epsilon, 1, 1, verbose );
            t_ifgt_nu = toc;
            fprintf('ifgt-nu      : %3.2e seconds, speedup = %3.2f\n', t_ifgt_nu, t_direct/t_ifgt_nu );

            % truncated series with approximate nearest neighbor structure,
            % with uniform parameter selection - implemented as in the FIGTree
            % paper (Raykar et al, submitted to Joural of Machine Learning).
            tic;
            [ g_ifgt_tree_u  ] = figtree( d, N, M, 1, X, h, q, Y, epsilon, 3, 0, verbose );
            t_ifgt_tree_u = toc;
            fprintf('ifgt-tree-u  : %3.2e seconds, speedup = %3.2f\n', t_ifgt_tree_u, t_direct/t_ifgt_tree_u );

            % truncated series, with non-uniform parameter selection.  Is faster 
            % than ifgt-nu if targets are not identically distributed to 
            % sources.
            tic;
            [ g_ifgt_tree_nu ] = figtree( d, N, M, 1, X, h, q, Y, epsilon, 3, 1, verbose );
            t_ifgt_tree_nu = toc;
            fprintf('ifgt-tree-nu : %3.2e seconds, speedup = %3.2f\n', t_ifgt_tree_nu, t_direct/t_ifgt_tree_nu );

            % show errors
            fprintf('\nMaximum absolute error from direct (exact) method:\n');
            fprintf('direct-tree  : %3.2e\n', max(abs(g_direct-g_direct_tree))/sum(q));
            fprintf('ifgt-u       : %3.2e\n', max(abs(g_direct-g_ifgt_u))/sum(q));
            fprintf('ifgt-nu      : %3.2e\n', max(abs(g_direct-g_ifgt_nu))/sum(q));
            fprintf('ifgt-tree-u  : %3.2e\n', max(abs(g_direct-g_ifgt_tree_u))/sum(q));
            fprintf('ifgt-tree-nu : %3.2e\n', max(abs(g_direct-g_ifgt_tree_nu))/sum(q));
            
            timings_direct(i,j,k) = t_direct;
            timings_direct_tree(i,j,k) = t_direct_tree;
            timings_ifgt(i,j,k) = t_ifgt_u;
            timings_ifgt_nu(i,j,k) = t_ifgt_nu;
            timings_ifgt_tree(i,j,k) = t_ifgt_tree_u;
            timings_ifgt_tree_nu(i,j,k) = t_ifgt_tree_nu;
        end;
    end;
    
    figure;
    loglog( hs, sum(timings_direct(:,:,k), 1), 'r.-', ...
            hs, sum(timings_direct_tree(:,:,k), 1), 'g.-', ...
            hs, sum(timings_ifgt(:,:,k), 1), 'b.-', ...
            hs, sum(timings_ifgt_nu(:,:,k), 1), 'c.-', ...
            hs, sum(timings_ifgt_tree(:,:,k), 1), 'm.-', ...
            hs, sum(timings_ifgt_tree_nu(:,:,k), 1), 'k.-' ...
        );
    legend('direct', 'direct-tree', 'ifgt', 'ifgt-nu', 'ifgt-tree', 'ifgt-tree-nu');
    xlabel('bandwidth(h)');
    ylabel('CPU Time');
    title(sprintf('CPU Time vs h (d = %i, N = %i, M = %i)', d, N, M));
    %pause;
end;