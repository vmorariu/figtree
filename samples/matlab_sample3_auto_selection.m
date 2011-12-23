% Script the demonstrate the use of FIGTree.  Note, only 
% (FIGTREE-DIRECTORY)/matlab/figtree.dll is needed for the figtree() function to 
% function.
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
N = 4000;
disp(sprintf('Number of source points N=%d',N));

% the number of targets
M = 2000;
disp(sprintf('Number of target points M=%d',M));

% the desired error
epsilon = 1e-2;     
disp(sprintf('Target error epsilon=%e',epsilon));

% verbosity
verbose = 0;

ds = 1:4;
hs = [.01 .02 .04 .08 .16 .32 .64];
ntrials = 1;

timings_direct      = zeros( ntrials, length(hs), length(ds) );
timings_direct_tree = zeros( ntrials, length(hs), length(ds) );
timings_ifgt        = zeros( ntrials, length(hs), length(ds) );
timings_ifgt_tree   = zeros( ntrials, length(hs), length(ds) );
timings_auto_method = zeros( ntrials, length(hs), length(ds) );

flops_direct      = zeros( ntrials, length(hs), length(ds) );
flops_direct_tree = zeros( ntrials, length(hs), length(ds) );
flops_ifgt        = zeros( ntrials, length(hs), length(ds) );
flops_ifgt_tree   = zeros( ntrials, length(hs), length(ds) );

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
            % the arguments 'ifgtParamMethod' and 'ifgtTruncMethod' are
            % irrelevant since we are not evaluating using ifgt.
            evalMethod = 0; ifgtParamMethod = 1; ifgtTruncMethod = 2;          
            tic; 
            [ g_direct           ] = figtree( X, h, q, Y, epsilon, evalMethod, ifgtParamMethod, ifgtTruncMethod, verbose ); 
            t_direct = toc;
            fprintf('direct      : %3.2e seconds\n', t_direct );

            % direct evaluation with tree data structure
            % the arguments 'ifgtParamMethod' and 'ifgtTruncMethod' are
            % irrelevant since we are not evaluating using ifgt.
            evalMethod = 2; ifgtParamMethod = 1; ifgtTruncMethod = 2;  
            tic;
            [ g_direct_tree       ] = figtree( X, h, q, Y, epsilon, 2, 0, 2, verbose );
            t_direct_tree = toc;
            fprintf('direct+tree : %3.2e seconds, speedup = %3.2f\n', t_direct_tree, t_direct/t_direct_tree );

            % ifgt evaluation, with non-uniform parameter selection (does
            % not assume uniform distribution), with cluster-wise
            % truncations
            evalMethod = 1; ifgtParamMethod = 1; ifgtTruncMethod = 2;
            tic; 
            [ g_ifgt     ] = figtree( X, h, q, Y, epsilon, evalMethod, ifgtParamMethod, ifgtTruncMethod, verbose );
            t_ifgt = toc;
            fprintf('ifgt        : %3.2e seconds, speedup = %3.2f\n', t_ifgt, t_direct/t_ifgt );

            % ifgt evaluation, with non-uniform parameter selection (does
            % not assume uniform distribution), with cluster-wise
            % truncations
            evalMethod = 3; ifgtParamMethod = 1; ifgtTruncMethod = 2;
            tic;
            [ g_ifgt_tree ] = figtree( X, h, q, Y, epsilon, evalMethod, ifgtParamMethod, ifgtTruncMethod, verbose );
            t_ifgt_tree = toc;
            fprintf('ifgt-tree   : %3.2e seconds, speedup = %3.2f\n', t_ifgt_tree, t_direct/t_ifgt_tree );

            % we run 'figtreeChooseEvaluationMethod' just so we can obtain
            % the estimated flops and to time method selection
            tic;
            [ best_method, flops ] = figtreeChooseEvaluationMethod( X, h, 1, Y, epsilon, 1, verbose );
            t_auto_method = toc;
           
            flops_direct(i,j,k)      = flops(1);
            flops_ifgt(i,j,k)        = flops(2);
            flops_direct_tree(i,j,k) = flops(3);
            flops_ifgt_tree(i,j,k)   = flops(4);            
            
            % ifgt evaluation, with non-uniform parameter selection (does
            % not assume uniform distribution), with cluster-wise
            % truncations
            evalMethod = 4; ifgtParamMethod = 1; ifgtTruncMethod = 2;          
            tic;
            [ g_auto ] = figtree( X, h, q, Y, epsilon, evalMethod, ifgtParamMethod, ifgtTruncMethod, verbose );
            t_auto = toc;
            fprintf('auto        : %3.2e seconds, speedup = %3.2f  (%3.2e seconds to choose method)\n', t_auto_method+t_auto, t_direct/(t_auto_method+t_auto), t_auto_method );
            
                       
            % show errors
            fprintf('\nMaximum absolute error from direct (exact) method:\n');
            fprintf('direct+tree  : %3.2e\n', max(abs(g_direct-g_direct_tree))/sum(q));
            fprintf('ifgt         : %3.2e\n', max(abs(g_direct-g_ifgt))/sum(q));
            fprintf('ifgt+tree    : %3.2e\n', max(abs(g_direct-g_ifgt_tree))/sum(q));
            fprintf('auto(%i)      : %3.2e\n', best_method, max(abs(g_direct-g_auto))/sum(q));
            
            timings_direct(i,j,k)      = t_direct;
            timings_direct_tree(i,j,k) = t_direct_tree;
            timings_ifgt(i,j,k)        = t_ifgt;
            timings_ifgt_tree(i,j,k)   = t_ifgt_tree;
            timings_auto_method(i,j,k) = t_auto_method + t_auto;
        end;
    end;
    
    figure;
    loglog( hs, sum(timings_direct(:,:,k), 1), 'r.-', ...
            hs, sum(timings_direct_tree(:,:,k), 1), 'g.-', ...
            hs, sum(timings_ifgt(:,:,k), 1), 'c.-', ...
            hs, sum(timings_ifgt_tree(:,:,k), 1), 'k.-', ...
            hs, sum(timings_auto_method(:,:,k), 1), 'y.-' ...
        );
    legend('direct', 'direct+tree', 'ifgt', 'ifgt+tree', 'auto');
    xlabel('bandwidth(h)');
    ylabel('CPU Time');
    title(sprintf('CPU Time vs h (d = %i, N = %i, M = %i)', d, N, M));

    figure;
    loglog( hs, sum(flops_direct(:,:,k), 1), 'r.-', ...
            hs, sum(flops_direct_tree(:,:,k), 1), 'g.-', ...
            hs, sum(flops_ifgt(:,:,k), 1), 'c.-', ...
            hs, sum(flops_ifgt_tree(:,:,k), 1), 'k.-' ...
        );
    legend('direct', 'direct+tree', 'ifgt', 'ifgt+tree');
    xlabel('bandwidth(h)');
    ylabel('Estimated Flops');
    title(sprintf('Estimated Flops vs h (d = %i, N = %i, M = %i)', d, N, M));

    %pause;
end;